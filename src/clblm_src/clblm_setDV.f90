!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_DV
   USE Module_ConstParam ,ONLY: r8=>kind_r8, FILLINT, FILLREAL

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: setDV, &
             CLBLM_SpectGrid, &
             CLBLM_SpectGrid_init, &
             CLBLM_SpectGrid_final


   TYPE :: CLBLM_SpectGrid
      integer           :: nLay   = 0
      real(r8)          :: V1     =FILLREAL !beginning wavenumber for the calculation
      real(r8)          :: V2     =FILLREAL !ending wavenumber for the calculation
      real ,allocatable :: DVnormal(:)      ![nLay];(cm-1);DV at coarse resolution for normal lines
      real ,allocatable :: DVnarrow(:)      ![nLay];(cm-1);DV at fine resolution for narrow lines
      real              :: DVOut  =FILLREAL !(cm-1);Output DV for OD from ODLAY; if<=-999 the output DV is either DVnormal or DVnarrow depends on the existence of narrow lines.
   END TYPE


   INTERFACE setDV
      module procedure setDV_struct
      module procedure setDV_array
   END INTERFACE


CONTAINS !==================== MODULE CONTAINS =========================


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_SpectGrid_init(this, nLay)
   !-------------------------------------------------------------------
      type(CLBLM_SpectGrid),    intent(inout) :: this
      integer,                  intent(in)    :: nLay

      if (allocated(this%DVnormal)) deallocate(this%DVnormal)
      allocate(this%DVnormal(nLay)) ;this%DVnormal(:) = 0.

      if (allocated(this%DVnarrow)) deallocate(this%DVnarrow)
      allocate(this%DVnarrow(nLay)) ;this%DVnarrow(:) = 0.

      this%nLay = nLay

   END SUBROUTINE

   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_SpectGrid_final(this)
   !-------------------------------------------------------------------
      type(CLBLM_SpectGrid),  intent(inout) :: this

      if (allocated(this%DVnormal)) deallocate(this%DVnormal)
      if (allocated(this%DVnarrow)) deallocate(this%DVnarrow)
      this%nLay = 0

   END SUBROUTINE


!-----------------------------------------------------------------------
! A wrapper subroutine to call setDV_struct to output DV in array instead
! of structure.
!-----------------------------------------------------------------------
   SUBROUTINE setDV_array( path, dvCtrl, gridType, V1,V2, DVarray, DVout )
!-----------------------------------------------------------------------
      USE Module_ConstParam    ,ONLY: r8=>kind_r8
      USE Module_AtmPath       ,ONLY: CLBLM_Path
      USE Module_Config        ,ONLY: CLBLM_DV_Ctrl
      IMPLICIT NONE

      type(CLBLM_Path)    ,intent(in)    :: path
      type(CLBLM_DV_Ctrl) ,intent(in)    :: dvCtrl
      integer             ,intent(in)    :: gridType
      real(r8)            ,intent(in)    :: V1,V2
      real                ,intent(out)   :: DVarray(:,:)
      real      ,optional ,intent(in)    :: DVout  !input user provided DV and output uniform output DV

      integer                :: nLay
      type(CLBLM_SpectGrid)  :: spectGrid


      call setDV_struct( spectGrid, path, dvCtrl, gridType, V1,V2,DVout )

      nLay = path%nLay
      DVarray(1:nLay,1) = spectGrid%DVnormal(1:nLay)
      DVarray(1:nLay,2) = spectGrid%DVnarrow(1:nLay)

      !if (present(DVout)) DVout = spectGrid%DVout

   END SUBROUTINE


!-----------------------------------------------------------------------
!
! Input/output combinations:
!         Input                Output
!   ----------------   -----------------------
!   gridType  DVout     DVnarror+DVnormal DVOUT
!     -1        Y            exact        -999.9
!     -1        N            exact        -999.9
!      0        Y         fixRat,DVSET    -999.9
!      0        N            fixRat       -999.9
!      1        Y         exact,DVSET     DVO
!      1        N            exact        min(DVnarrow)
!   DVSET happened when difference between DV and DVO is less than 20%DV
!
!-----------------------------------------------------------------------
   SUBROUTINE setDV_struct( spectGrid, path, dvCtrl, gridType, V1,V2, DVout )
!-----------------------------------------------------------------------
      USE Module_ConstParam    ,ONLY: r8=>kind_r8, FILLREAL
      USE Module_AtmPath       ,ONLY: CLBLM_Path
      USE Module_Config        ,ONLY: CLBLM_DV_Ctrl
      IMPLICIT NONE

      type(CLBLM_SpectGrid) ,intent(out) :: spectGrid
      type(CLBLM_Path)      ,intent(in)  :: path
      type(CLBLM_DV_Ctrl)   ,intent(in)  :: dvCtrl
      integer               ,intent(in)  :: gridType
      real(r8)              ,intent(in)  :: V1,V2
      real        ,optional ,intent(in)  :: DVout  !input user provided DV


      integer  :: nLay
      real(r8) :: Vbar
      real     :: DVC(path%nLay) !exact DVs
      real     :: DVL(path%nLay)
      real     :: DVC_fine(path%nLay)
      real     :: DVL_fine(path%nLay)
      real     :: DVO

      nLay = path%nLay
      DVO = FILLREAL
      if (present(DVout)) DVO=DVout


      !--- Calculate exact DV (DVC) for each layer

      ! AERI SGP case was used to test the spectral range.
      ! We determined the spectral range is not limited in CLBLM.
      ! The maximum spectral range for a same DV is 2020 cm-1.
      !IF ((V2-V1).GT.2020.) THEN
        ! STOP '--- setDV(): V2-V1 .GT. 2020'
      !ENDIF

      Vbar = 0.5*(V1 + V2)
      IF (V2.LT.V1) Vbar = V1

      !---Calculate DV for normal lines
      call calculateDV( DVC, Vbar, path, dvCtrl%ALFAL0, dvCtrl%SAMPLE )

      !---Calculate finer DV for narrow lines
      call calculateDV( DVC_fine, Vbar, path, dvCtrl%ALFAL0_fine, dvCtrl%SAMPLE )



      !--- If needed, adjust DVs
      !
      DVL = DVC
      DVL_fine = DVC_fine

      if (gridType >=0) then !not asking for exact DV

         if ( gridType ==0 ) then  !adjusted DV
            CALL fixDVRatio(DVL,      DVC,      nLay)
            CALL fixDVRatio(DVL_fine, DVC_fine, nLay)
         endif

         if ( DVO >0. ) then !User provided DV and it is not very much away from calculated DV, adjusted calculated DV according to userDV
            if ( abs( 1.-DVO/DVL(     nLay) ) <0.2 ) CALL DVSET(  DVL,      DVL,      DVO, nLay )
            if ( abs( 1.-DVO/DVL_fine(nLay) ) <0.2 ) CALL DVSET(  DVL_fine, DVL_fine, DVO, nLay )
         endif

      endif



      !--- Load the spectGrid object.
      !
      call CLBLM_SpectGrid_init( spectGrid, nLay )

      !Load V1,V2 and nLay
      spectGrid%nLay = nLay
      spectGrid%V1 = V1
      spectGrid%V2 = V2

      !Load DVnormal and DVnarrow
      spectGrid%DVnormal(1:nLay) = DVL(1:nLay)
      spectGrid%DVnarrow(1:nLay) = DVL_fine(1:nLay)

      !Load DV
      ! If uniform grid is requested, the output DV is a same value for
      ! all layers. The value could be a user input DV or the finest DV
      ! for all layers. If gridType/=1, the output DV could be either
      ! DVnormal or DVnarrow depending on the existence of narrow lines
      ! within the spectral range.
      !
      if (gridType == 1) then !uniform output DV
         if (DVO >0.) then
            spectGrid%DVout = DVO !uniform DV, all interpolate to DVout
         else
            spectGrid%DVout = minval( spectGrid%DVnarrow(1:nLay) ) !uniform DV, all interpolate to the finest DV
         endif
      else
         spectGrid%DVout = FILLREAL !non-uniform DV, output DV could be either DVnormal or DVnarrow
      endif

      !spectGrid%dvCtrl = dvCtrl !Save the control information used for this spectGrid.

   END SUBROUTINE


!-----------------------------------------------------------------------
! Calculate exact DV for each layer
!-----------------------------------------------------------------------
   SUBROUTINE calculateDV( DVC, Vbar, path, ALFAL0,SAMPLE  )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8, molIndex, &
                                       P0=>Press1013, T0=>Temp296, &
                                       AVMASS=>AVMWT
      USE Module_AtmPath        ,ONLY: CLBLM_Path

      IMPLICIT NONE

      real,             intent(out) :: DVC(:) !(path%nLay) !out
      real(r8),         intent(in)  :: Vbar
      type(CLBLM_Path), intent(in)  :: path
      real,             intent(in)  :: ALFAL0
      real,             intent(in)  :: SAMPLE

      integer :: nLay, nMol, il
      real    :: SUMWK, WTOTL, FRH2O, ALFCOR, H2OSLF, ALBAR,ADBAR,AVBAR
      real    :: DV
      !integer :: IHIRAC


      nLay = path%nLay
      nMol = path%nMol

      !ALFAL0 = dvCtrl%ALFAL0
      !SAMPLE = dvCtrl%SAMPLE !sampling rate
      !IHIRAC = odCtrl%IHIRAC !line shape

      do il = 1,nLay

!wtot         SUMWK = sum( path%W(1:nMol,il) )
!wtot         WTOTL = SUMWK + path%Wbroad(il)
         WTOTL = path%Wtot(il)
         FRH2O = path%W( molIndex('H2O', path%molID), il ) / WTOTL

         ALFCOR = ( path%Pave(il)/P0) * SQRT(T0/path%Tave(il))
         H2OSLF = (1.-FRH2O+5.*FRH2O) ! CORRECT FOR WATER SELF BROADENING

         ALBAR = ALFAL0*ALFCOR*H2OSLF
         !3.58115E-07 = SQRT( 2.* LOG(2.)*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) )
         ADBAR = 3.58115E-07*Vbar*SQRT(path%Tave(il)/AVMASS)
         AVBAR = 0.5*(ALBAR+SQRT(ALBAR*ALBAR+4.*ADBAR*ADBAR))

         DV = AVBAR/SAMPLE
         !IF (IHIRAC.EQ.2) DV = ALBAR/SAMPLE
         !IF (IHIRAC.EQ.3) DV = ADBAR/SAMPLE

         DVC(il) = DV
      enddo

   END SUBROUTINE

!-----------------------------------------------------------------------
! This subroutine adjusts DV for each layer to maintain a fixed DV ratio
! between adjacent layers.
!-----------------------------------------------------------------------
   SUBROUTINE fixDVRatio( DVout, DVin, nLay )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real,    intent(out) :: DVout(:)  !(nLay)
      real,    intent(in)  :: DVin(:)  !(nLay)
      integer, intent(in)  :: nLay

      integer, PARAMETER :: I_2 = 2
      real,    PARAMETER :: TYPMAX = 2.5

      integer :: il
      real    :: DV, OLDDV
      real    :: TYPE
      integer :: ITYPE
      real    :: SCAL
      integer :: ISCAL,IDV

      do il =1,nLay

         DV = DVin(il)
         TYPE = 0.
         ITYPE = 99

         ! DV is assumed to be less than 1
         ! Set DV to 3 significant figures
         IF (il.EQ.1) THEN

            ISCAL = LOG10(DV)-3.
            SCAL = 10.**ISCAL
            IDV = (DV/SCAL)+0.5

            ! Set IDV to be even
            IF (MOD(IDV,I_2).GT.0) IDV = IDV+1
            DV = SCAL* REAL(IDV)

         ELSE

            TYPE = OLDDV/DV

            IF (TYPE.GT.TYPMAX) THEN
               WRITE (*,962) TYPMAX
962            FORMAT(20X,'  DV RATIO  .GT. ',F10.2)
               STOP
            ELSEIF (TYPE.GE.1.2) THEN

               ! TYPE is between 1.2 and TYPMAX
               DV = OLDDV
               ITYPE = 1./(TYPE-1.)+0.5
               IF (ITYPE.EQ.3) ITYPE = 2
               DV = OLDDV* REAL(ITYPE)/ REAL(ITYPE+1)

            ELSEIF (TYPE.GE.0.8) THEN

               ! TYPE is between 0.8 and 1.2 (set to 1.0)
               DV = OLDDV
               ITYPE = 0

            ELSE

               ! TYPE is less than 0.8
               DV = OLDDV
               ITYPE = 0

            ENDIF

         ENDIF !IF (layer.EQ.1) THEN

         OLDDV = DV
         DVout(il) = DV
      enddo !do il =1,nLay

   END SUBROUTINE


!-----------------------------------------------------------------------
! Set the final layer DV to the value of userDV and scale all the DVs based on
! ratio=userDV/DV(nLay). userDV must be within +/-20% of the calculated DV.
!-----------------------------------------------------------------------
   SUBROUTINE DVSET(  DVout, DVin, userDV, nLay )
!-----------------------------------------------------------------------
      real    ,intent(out) :: DVout(:)
      real    ,intent(in)  :: DVin(:)
      real    ,intent(in)  :: userDV
      integer ,intent(in)  :: nLay

      character(*), parameter :: routineName='DVSET'
      real    :: RATIO
      integer :: L


      RATIO = userDV/DVin(nLay)

      !--- Test to be sure that userDV is not more than
      !     20% different than the monochromatic DV.
      IF (RATIO.GT.1.2.OR.RATIO.LT.0.8 .or. userDV<0.) THEN
         !WRITE (IPR,967) RATIO,DVSET,DV
         STOP '--- '//routineName//'(): userDV must be >0. and within +/-20% of calculated final layer DV.'
      ENDIF

      !--- Scale DVs
      do L = 1, nLay
         DVout(L) = DVin(L)*RATIO
      enddo

   END SUBROUTINE

END MODULE
