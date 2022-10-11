!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_Spectrum
   USE Module_ConstParam, ONLY: r8=>kind_r8

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_Spectrum, &
             CLBLM_Spectrum_init, &
             CLBLM_Spectrum_final, &
             spectraHaveSameGrid, &
             copySpectrum, &
             moveSpectrum, &
             trimSpectrum, &
             resizeSpectrum, &
             interpSpectrum, &
             CLBLM_SpectrumPointer, &
             interpToFinestGrid, &
             XINT, &
             pre_xint


   !---------------------------------------------
   ! * A Spectrum:
   !                ........................
   !                |                      |
   ! frequency:     V1                     V2
   !                |                      |
   ! array index: indV1                  indV2
   !
   !---------------------------------------------
   TYPE :: CLBLM_Spectrum
      real, allocatable :: spect(:)
      real              :: DV        =0.
      real(r8)          :: V1        =0.
      real(r8)          :: V2        =0.
      integer           :: NLIM      =0  !num. of grid point from V1 to V2 inclusively
      integer           :: indV1     =0  !index of V1
      integer           :: indV2     =0  !index of V2
   END TYPE



   !INTERFACE operator(+)
   !   module procedure vector_add
   !END INTERFACE
   !
   !INTERFACE operator(*)
   !   module procedure vector_mult
   !END INTERFACE


   TYPE CLBLM_SpectrumPointer
      type(CLBLM_Spectrum), POINTER :: ptr=>null()
   END TYPE

   INTERFACE interpSpectrum
      module procedure interpSpectrum_array2array
      module procedure interpSpectrum_array2spect
      module procedure interpSpectrum_spect2spect
      module procedure interpSpectrum_spect
   END INTERFACE

   INTERFACE interpToFinestGrid
      module procedure interpToFinestGrid_argument
      module procedure interpToFinestGrid_pointerArray
      module procedure interpToFinestGrid_array
   END INTERFACE


CONTAINS !======================= MODULE CONTAINS ======================

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE CLBLM_Spectrum_init( this, V1,DV,NLIM, indV1 )
   !--------------------------------------------------------------------
      type(CLBLM_Spectrum), intent(inout) :: this
      real(r8),             intent(in)    :: V1
      real,                 intent(in)    :: DV
      integer,              intent(in)    :: NLIM
      integer,    OPTIONAL, intent(in)    :: indV1

      integer :: iV1

      iV1 = 1
      if (present(indV1)) iV1 = indV1

      this%DV    = DV
      this%V1    = V1
      this%V2    = V1 + DV*(NLIM-1)
      this%NLIM  = NLIM                   !num. of grid point from V1 to V2 inclusively
      this%indV1 = iV1                    !index of V1
      this%indV2 = iV1 + NLIM - 1         !index of V2

      if (allocated(this%spect))  deallocate(this%spect)
      allocate( this%spect( this%indV1 : this%indV2 ) )
!yma 170320: array initialization is very time consuming, comment this out.
      !this%spect(:) = 0.

   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE CLBLM_Spectrum_final( this )
   !--------------------------------------------------------------------
      type(CLBLM_Spectrum), intent(inout) :: this

      if (allocated(this%spect))  deallocate(this%spect)

      this%DV    =0.
      this%V1    =0.
      this%V2    =0.
      this%NLIM  =0  !num. of grid point from V1 to V2 inclusively
      this%indV1 =0  !index of V1
      this%indV2 =0  !index of V2

   END SUBROUTINE

   !!--------------------------------------------------------------------
   !!--------------------------------------------------------------------
   !FUNCTION vector_add( sp1, sp2 ) result(sumSp)
   !!--------------------------------------------------------------------
   !   type(CLBLM_Spectrum), intent(in) :: sp1,sp2
   !   type(CLBLM_Spectrum)             :: sumSp
   !
   !   integer :: ip
   !
   !
   !   if (.not.spectraHaveSameGrid(sp1,sp2)) &
   !      STOP '--- vector_add(): Right-hand-side spectra not in same grid.'
   !
   !   if (allocated(sumSp%spect)) then
   !   !--- sumSp already exist.
   !
   !      if (.not.spectraHaveSameGrid(sp1,sumSp)) &
   !         STOP '--- vector_add(): Exist left-hand-side spectrum not in same grid as Right-hand-side.'
   !
   !      do ip = sumSp%indV1 , sumSp%indV2
   !         sumSp%spect(ip) = sp1%spect(ip) + sp2%spect(ip)
   !      enddo
   !
   !   else
   !   !--- sumSp not yet allocated
   !
   !      call CLBLM_Spectrum_init( sumSp, sp1%V1, sp1%DV, sp1%NLIM, sp1%indV1 )
   !
   !      do ip = sumSp%indV1 , sumSp%indV2
   !         sumSp%spect(ip) = sp1%spect(ip) + sp2%spect(ip)
   !      enddo
   !
   !   endif
   !
   !END FUNCTION


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   logical FUNCTION spectraHaveSameGrid( sp1, sp2 )
   !--------------------------------------------------------------------
      type(CLBLM_Spectrum), intent(in) :: sp1,sp2

      spectraHaveSameGrid = &
          (allocated(sp1%spect)   .and. &
           allocated(sp2%spect)   .and. &
           abs(sp1%V1 - sp2%V1) <epsilon(sp1%V1) .and. &
           abs(sp1%V2 - sp2%V2) <epsilon(sp1%V2) .and. &
           abs(sp1%DV - sp2%DV) <epsilon(sp1%DV) .and. &
           sp1%NLIM  == sp2%NLIM  .and. &
           sp1%indV1 == sp2%indV1 .and. &
           sp1%indV2 == sp2%indV2 )

   END FUNCTION


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE copySpectrum( srcSpect, dstSpect )
   !--------------------------------------------------------------------
      type(CLBLM_Spectrum), intent(in)  :: srcSpect !source spectrum
      type(CLBLM_Spectrum), intent(out) :: dstSpect !destination spectrum

      dstSpect = srcSpect

      !call CLBLM_Spectrum_init( dstSpect, &
      !                          srcSpect%V1, &
      !                          srcSpect%DV, &
      !                          srcSpect%NLIM, &
      !                          indV1=srcSpect%indV1 )
      !
      ! dstSpect%spect( dstSpect%indV1 : dstSpect%indV2 ) = &
      ! srcSpect%spect( srcSpect%indV1 : srcSpect%indV2 )

   END SUBROUTINE


   !--------------------------------------------------------------------
   ! * 'moveSpectrum' copy the content of 'srcSpect' to 'desSpect' and
   !   finalize the 'srcSpect" (delete the content of the 'srcSpect').
   ! * 'moveSpectrum' may be faster than 'copySpectrum', if the source
   !   spectrum is not need to be kept.
   !--------------------------------------------------------------------
   SUBROUTINE moveSpectrum( srcSpect, dstSpect )
   !--------------------------------------------------------------------
      type(CLBLM_Spectrum), intent(inout) :: srcSpect !source spectrum
      type(CLBLM_Spectrum), intent(out)   :: dstSpect !destination spectrum

      dstSpect%DV         = srcSpect%DV
      dstSpect%V1         = srcSpect%V1
      dstSpect%V2         = srcSpect%V2
      dstSpect%NLIM       = srcSpect%NLIM        !num. of grid point from V1 to V2 inclusively
      dstSpect%indV1      = srcSpect%indV1       !index of V1
      dstSpect%indV2      = srcSpect%indV2       !index of V2

      call move_alloc( srcSpect%spect, dstSpect%spect )

      call CLBLM_Spectrum_final(srcSpect)

   END SUBROUTINE

   !--------------------------------------------------------------------
   ! When reading in the spectra data from LBLRTM file. The initial number
   ! of points my not be the exact NLIM. To keep the consistency inside
   ! a CLBLM_Spectrum object, this subroutine trim the empty points at the
   ! end of the spectrum. This subroutine adjusts the indexes only, no
   ! data copy occurred.
   !--------------------------------------------------------------------
   SUBROUTINE trimSpectrum( this, actNLIM )
   !--------------------------------------------------------------------
      type(CLBLM_Spectrum) ,intent(inout) :: this
      integer              ,intent(in)    :: actNLIM

      integer :: dNLIM

      !--- Already in exact size, do nothing.
      if (this%NLIM == actNLIM)  RETURN

      !--- Can only trim short
      if (this%NLIM < actNLIM) STOP '--- trimSpectrum(): actNLIM must be less than the current spectral length.'

      !---
      dNLIM = this%NLIM - actNLIM

      !---
      !this%DV         = this%DV
      !this%V1         = this%V1
      this%V2         = this%V2 - dNLIM*this%DV
      this%NLIM       = actNLIM                              !num. of grid point from V1 to V2 inclusively
      !this%indV1      = this%indV1                          !index of V1
      this%indV2      = this%indV1 + this%NLIM - 1           !index of V2

   END SUBROUTINE

   !--------------------------------------------------------------------
   ! Resize the spectrum and add the input dataArr to the end of the spectrum.
   !--------------------------------------------------------------------
   SUBROUTINE resizeSpectrum( spect, dataArr,V1,DV,NP )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      type(CLBLM_Spectrum) ,intent(inout) :: spect
      real                 ,intent(in)    :: dataArr(:)
      real(r8)             ,intent(in)    :: V1
      real                 ,intent(in)    :: DV
      integer              ,intent(in)    :: NP

      character(*) ,parameter :: routineName = 'resizeSpectrum'
      real         ,parameter :: TOL = 1.e-8
      integer           :: NLIM
      real ,allocatable :: tempArr(:)

      !--- If spect is empty, fill the spect and return
      if ( spect%NLIM==0 ) then
         call CLBLM_Spectrum_init( spect, V1,DV,NP )
         spect%spect(1:NP) = dataArr(1:NP)
         RETURN
      endif


      if ( abs(DV-spect%DV) >TOL ) STOP '--- '//routineName//'(): The two DVs are different, can not resize spectrum.'
      if ( abs(V1-spect%V2-DV) >TOL ) STOP '--- '//routineName//'(): The gap between V2 of the old spectrum and V1 of the new data larger than DV, can not resize spectrum.'

      NLIM = spect%NLIM
      allocate(tempArr( NLIM + NP ))

      tempArr(1:NLIM) = spect%spect(:)
      tempArr(NLIM+1:NLIM+NP) = dataArr(1:NP)

      call CLBLM_Spectrum_init( spect, spect%V1,spect%DV, NLIM+NP )
      call move_alloc( tempArr, spect%spect )

   END SUBROUTINE

   !--------------------------------------------------------------------
   ! * Apply to aligned (with same V1) and not-aligned cases
   ! * Apply to low to high resolution and high to low resolution.
   ! * The output grid resolution is specified in "DVout"
   ! * It is assumed that the input spectral range is wider than output
   !   spectral range.
   ! * Input V1 and out V1 need not to the same, but input V1 has to be
   !   less than the output V1
   !
   !   V1in
   !    |      V1out
   !    |       |
   !    |       |
   !    |       .   .   .   .   .   .   .     spectOut%spect
   !    .     .     .     . |   .     .     . spInExt
   !               -1     0     1     2
   !
   !    if DVin > DVout, the spectOut%spect(2) is interpolated as:
   !       .  .  .  .  .  intervalLen_out = 4
   !   .   .  |.   .   .  intervalLen_in  = 3
   !  -1   0   1   2
   !
   !    if DVin < DVout, the spectOut%spect(2) is interpolated as:
   !       .   .   .   .  intervalLen_out = 3
   !       .  .| .  .  .  intervalLen_in  = 4
   !      -1  0  1  2
   !
   !--------------------------------------------------------------------
   SUBROUTINE interpSpectrum_array2array( SPin,DVin,V1in,NLIMin, &
                                          SPout,DVout,V1out,NLIMout )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real     ,intent(in)  :: spIn(1:)
      real     ,intent(in)  :: DVin
      real(r8) ,intent(in)  :: V1in
      integer  ,intent(in)  :: NLIMin
      real     ,intent(out) :: spOut(1:)
      real     ,intent(in)  :: DVout
      real(r8) ,intent(in)  :: V1out
      integer  ,intent(in)  :: NLIMout

      character(*) ,parameter :: routineName='interpSpectrum'

      real     ,PARAMETER :: TOL = 1.e-8
      integer  ,PARAMETER :: NUMCOF = 201
      real          ,SAVE :: A1(0:202),A2(0:202),A3(0:202),A4(0:202), XNUMCF
      logical       ,SAVE :: firsttime=.TRUE.
      real                :: B1(10),B2(10),B3(10),B4(10)
      real                :: c1(0:100),c2(0:100),c3(0:100),c4(0:100)

      integer :: I,J,J1,J2,iIn,iOut,JP,IP,NL,JJ
      logical :: spectraHaveSameGrid
      real(r8):: V2in, V2out, Vo, Vi
      integer :: I1out, I2out
      real    :: dvRatio
      logical :: sameV1
      integer :: intervalLen_in
      integer :: intervalLen_out
      integer :: ITYPE !if partially aligned case, ITYPE gives the number points need to be interpolated in between two aligned points.
      real    :: AP, P, PP
      real    :: A1N,A2N,A3N,A4N
      real    :: FJ1DIF, FJJ

      integer ,PARAMETER :: nExtPoints=4
      real  ,allocatable :: spInExt(:)



      !--- Return if the two spectra have same grid
      spectraHaveSameGrid = &
          (abs(V1in - V1out) <epsilon(V1in) .and. &
           abs(DVin - DVout) <epsilon(DVin) .and. &
           NLIMin  == NLIMout )
      if (spectraHaveSameGrid) RETURN


      ! SET UP FOUR POINT INTERPOLATION COEFICIENTS FOR P FOR 201
      ! POINTS BETWEEN 0 AND 1.0, with an extra point at each end
      IF (firsttime) THEN

         XNUMCF = REAL(NUMCOF)
         DO IP = 0, NUMCOF+1
            P = ( REAL(IP)-1.0)/(XNUMCF-1.0)
            PP = P**2
            A1(IP) = -P/2.0*(1-P)**2
            A2(IP) = 1.0-PP*(3.0-2.0*P)+PP/2.0*(1.0-P)
            A3(IP) = PP*(3.0-2.0*P)+P/2.0*(1.0-P)**2
            A4(IP) = -PP/2*(1.0-P)
         ENDDO

         firsttime = .FALSE.
      ENDIF

      V2in  = V1in + real(NLIMin-1)*DVin
      V2out = V1out + real(NLIMout-1)*DVout
      ! 10/3/22 Betsy Berry commented out the two following lines (and uncommened the one above)
      ! Reason Igor discovered: It was causing an error during interpolation to the finest grid.
      ! Since the last point was omitted the last point mrgRadDn%spect was set to zero and thus a delta function was created.
      ! Consequently, it results in a sinc-like pattern in the vicinity of the last point.
      !V2out = V1out + real(NLIMout-2)*DVout
      !I2out = NLIMout-1
      I1out = 1
      I2out = NLIMout


      !--- Check to make sure the input range covers the output range at starting point
      if ( V1out < V1in ) then
         print *,'V1in, V1out:', V1in, V1out
         STOP '--- '//routineName//'(): Output V1 can not be less than input V1.'
      endif

      !--- No more than 1 DV distance extrapolation will occur.
      if ( V2in + DVin < V2out ) then
         print *,'V2in, V2out:', V2in, V2out

         ! STOP '--- '//routineName//'(): Output V2 can not exceed input V2 by more than one DVin'
      endif

      !--- Check if the two spectra has the same starting point
      ! Fixed ratio interpolation only allowed when input and output spectra
      ! have the same starting point.
      sameV1 = abs(V1in - V1out) <epsilon(V1in)


      !---
      dvRatio = DVOUT/DVIN
      IF (ABS(DVOUT-DVIN).LT.TOL) dvRatio = 1.

      if ( sameV1 .and. abs(dvRatio - 1.) <epsilon(dvRatio) ) then  !same resolution

         intervalLen_in  = 1
         intervalLen_out = 1
         ITYPE = intervalLen_out -1 !=0

      elseif (sameV1 .and. any( abs(dvRatio-[ 2./1., 3./2., 4./3., 5./4., 6./5. ]) <TOL ) ) then !high to low resolution

         intervalLen_in  = int(DVOUT/(DVOUT-DVIN)+0.5)
         intervalLen_out = intervalLen_in - 1 ! =int(DVIN/(DVOUT-DVIN)+0.5)
         ITYPE = intervalLen_out -1 !ITYPE={0,1,2,3}; =int(DVIN/(DVOUT-DVIN)+0.5) - 1

      elseif (sameV1 .and. any( abs(dvRatio-[ 1./2., 2./3., 3./4., 4./5., 5./6. ]) <TOL ) ) then !low to high resolution

         intervalLen_in  = int(DVOUT/(DVIN-DVOUT)+0.5)
         intervalLen_out = intervalLen_in + 1 ! =int(DVIN/(DVIN-DVOUT)+0.5)
         ITYPE = intervalLen_out -1 !ITYPE={1,2,3,4}; =int(DVOUT/(DVIN-DVOUT)+0.5)

      else
         intervalLen_in  = 0
         intervalLen_out = 0
         ITYPE = intervalLen_out -1 !=-1
      endif



      !---
      ! * Extend the SPin to include a header and a tail of size of 4 points
      !   to facilitate the 4-points Lagrangian interpolation.
      allocate( spInExt( 1-nExtPoints : NLIMin+nExtPoints ) )
      spInExt(1:NLIMin) = spIn( 1:NLIMin )

      !--- Add points to the head and tail of the input spectrum for interpolation
      !spInExt(-1) = spInExt(1)
      !spInExt( 0) = spInExt(1)
      !spInExt(NLIMin+1) = spInExt(NLIMin)
      !spInExt(NLIMin+2) = spInExt(NLIMin)
      spInExt( 0) = 2.0*spInExt(1) - spInExt(2)
      spInExt(-1) = spInExt(0)
      spInExt(NLIMin+1) = 2.0*spInExt(NLIMin) - spInExt(NLIMin-1)
      spInExt(NLIMin+2) = spInExt(NLIMin+1)


      !--- Perform interpolation
      !
      if (ITYPE == 0) then !IYTPE=0 when dvRatio==1 or dvRatio==2/1

         iIn = 1 - nint(dvRatio)
         do I = I1out, I2out
            iIn = iIn + nint(dvRatio)
            spOut(I) = spInExt(iIn)
         enddo

      elseif (ITYPE >0) then !partially aligned, ratio is n/(n+1) or (n+1)/n

         AP = 1.0/(real(ITYPE)+1.0) != 1./intervalLen_out
         DO JP = 1, ITYPE
            IP = JP+1
            P = AP*real(JP)
            if ( DVin > DVout ) P = 1.0-P
            B1(IP) = -P*(P-1.0)*(P-2.0)/6.0
            B2(IP) = (P**2-1.0)*(P-2.0)*0.5
            B3(IP) = -P*(P+1.0)*(P-2.0)*0.5
            B4(IP) = P*(P**2-1.0)/6.0
         END DO

         do NL = 1,intervalLen_out

            iIn = NL - intervalLen_in
            if ( NL > 1 .and. DVin > DVout ) iIn = iIn-1

            if (NL ==1) then !EXACT FREQUENCY - NO INTERPOLATION

               do I = I1out, I2out, intervalLen_out
                  iIn = iIn + intervalLen_in
                  spOut(I) = spInExt(iIn)
               enddo

            else !NOT EXACT FREQUENCY - INTERPOLATE RESULT

               A1N = B1(NL)
               A2N = B2(NL)
               A3N = B3(NL)
               A4N = B4(NL)
               do I = I1out+NL-1, I2out, intervalLen_out

                  iIn = iIn + intervalLen_in
                  spOut(I) = A1N * spInExt(iIn-1) + &
                             A2N * spInExt(iIn)   + &
                             A3N * spInExt(iIn+1) + &
                             A4N * spInExt(iIn+2)

               enddo

            endif !if (NL ==1)

         enddo !do NL = 1,intervalLen_out

      else  !(ITYPE <0), no alignment assumed.

         !
         !                      J1
         !      V1out           Vo
         !      .   .   .   .   .   .
         !   .       .       .       .       .
         !   V1in            Vi
         !                   iIn
         !
         ! P IS INCREMENTED BY ADDING DVout/DVin BUT WILL BE REINITIALIZED
         ! HERE FOR EACH PANEL TO AVOID THE ACCUMULATION OF TRUNCATION ERRORS

         Vo = V1out
         DO WHILE ( Vo < V2out+TOL )

            iIn = (Vo-V1in)/DVin+1.00001  !--iIn INDEXES THE LARGEST Vi .LE. Vo
            Vi = V1in + REAL(iIn-1)*DVin
            P = (Vo-Vi)/DVin

            J1 = INT((Vo-V1out)/DVout+1.001) + I1out-1
            J2 = MIN( I2out, J1+2400-1 ) + I1out-1
            DO J = J1, J2
               IP = P*XNUMCF+1.00001
               spOut(J) = A1(IP)*spInExt(iIn-1) + &
                          A2(IP)*spInExt(iIn)   + &
                          A3(IP)*spInExt(iIn+1) + &
                          A4(IP)*spInExt(iIn+2)
               ! INCREMENT P AND iIn
               P = P + dvRatio
               IF (P.GE.1.0) THEN
                  iIn = iIn + P
                  P = P- REAL(INT(P))
               ENDIF
            ENDDO !DO J = J1, J2

            ! Vo IS THE FREQUENCY OF THE NEXT OUTPUT POINT
            Vo = V1out + J2*DVout

         ENDDO !do


         !DO JP = 0, 100
         !   P = 0.01*real(JP)
         !   c1(JP) = -P*(P-1.0)*(P-2.0)/6.0
         !   c2(JP) = (P**2-1.0)*(P-2.0)*0.5
         !   c3(JP) = -P*(P+1.0)*(P-2.0)*0.5
         !   c4(JP) = P*(P**2-1.0)/6.0
         !END DO
         !
         !!RATDV  = DVOUT/DVIN
         !FJ1DIF = (V1out-V1in)/DVIN+1.+2. !FJJ is offset by 2. for rounding purposes
         !DO iOut = I1out, I2out
         !   FJJ = FJ1DIF + dvRatio*REAL(iOut-I1out)  !FJJ = FJ1DIF + dvRatio*REAL(II-1)
         !   JJ = INT(FJJ)-2
         !   JP = (FJJ- REAL(JJ))*100.-199.5
         !   iIn = JJ + 1 - 1
         !   spOut(iOut) = c1(JP)*spInExt(iIn-1) + &
         !                 c2(JP)*spInExt(iIn)   + &
         !                 c3(JP)*spInExt(iIn+1) + &
         !                 c4(JP)*spInExt(iIn+2)
         !END DO

      endif !if (ITYPE == 0) then

      !where (spOut <0.) spOut = 0. ! IF (ODLAYI.LT.0.) ODLAYI = 0.
   END SUBROUTINE



   !--------------------------------------------------------------------
   ! Interpolate a array to finer grid and output to a CLBLM_Spectrum object.
   !--------------------------------------------------------------------
   SUBROUTINE interpSpectrum_array2spect( SPin,DVin,V1in,NLIMin, spectOut )
   !--------------------------------------------------------------------
      IMPLICIT NONE

      real                 ,intent(in)    :: SPin(:)
      real                 ,intent(in)    :: DVin
      real(r8)             ,intent(in)    :: V1in
      integer              ,intent(in)    :: NLIMin
      type(CLBLM_Spectrum) ,intent(inout) :: spectOut


      call interpSpectrum_array2array( SPin, DVin, V1in, NLIMin, &
                                       spectOut%spect( spectOut%indV1:spectOut%indV2 ), &
                                       spectOut%DV, &
                                       spectOut%V1, &
                                       spectOut%NLIM )

   END SUBROUTINE


   !--------------------------------------------------------------------
   ! Interpolate spectIn to spectOut. Both arguments are CLBLM_Spectum objects.
   ! The output grid is defined by spectOut. spectOut will be filled with the
   ! interpolated values from spectIn.
   !--------------------------------------------------------------------
   SUBROUTINE interpSpectrum_spect2spect( spectIn, spectOut )
   !--------------------------------------------------------------------
      IMPLICIT NONE
      type(CLBLM_Spectrum), intent(in)    :: spectIn
      type(CLBLM_Spectrum), intent(inout) :: spectOut


      if ( spectraHaveSameGrid(spectIn, spectOut) ) RETURN

      call interpSpectrum_array2array( spectIn%spect( spectIn%indV1:spectIn%indV2 ), &
                                       spectIn%DV, &
                                       spectIn%V1, &
                                       spectIn%NLIM, &
                                       spectOut%spect( spectOut%indV1:spectOut%indV2 ), &
                                       spectOut%DV, &
                                       spectOut%V1, &
                                       spectOut%NLIM )

   END SUBROUTINE

   !--------------------------------------------------------------------
   ! Interpolate spect to a new grid given by DV,V1,NLIM and indV1. The
   ! content of spect will replaced with the newly interpolated values.
   !--------------------------------------------------------------------
   SUBROUTINE interpSpectrum_spect( spect, DV,V1,NLIM,indV1 )
   !--------------------------------------------------------------------
      IMPLICIT NONE
      type(CLBLM_Spectrum) ,intent(inout) :: spect
      real                 ,intent(in)    :: DV
      real(r8)             ,intent(in)    :: V1
      integer              ,intent(in)    :: NLIM
      integer    ,optional ,intent(in)    :: indV1

      integer :: ind1
      type(CLBLM_Spectrum) :: tmpSp


      ind1 = 1
      if (present(indV1)) ind1=indV1

      call CLBLM_Spectrum_init( tmpSp, V1, DV, NLIM, ind1 )

      if ( spectraHaveSameGrid(spect,tmpSp) ) RETURN

      call interpSpectrum_array2array( spect%spect( spect%indV1:spect%indV2 ), &
                                       spect%DV, &
                                       spect%V1, &
                                       spect%NLIM, &
                                       tmpSp%spect( tmpSp%indV1:tmpSp%indV2 ), &
                                       tmpSp%DV, &
                                       tmpSp%V1, &
                                       tmpSp%NLIM )

      call moveSpectrum( tmpSp, spect )

   END SUBROUTINE

!-----------------------------------------------------------------------
! This is a wrapper of the interpToFinestGrid_pointerArray(). It make a
! pointer array pointing to the input spectra and call interpToFinestGrid_pointerArray()
! to do interpolation. All input arguments are optional, so it is like
! a subroutine that can take variable number of inputs.
! The maximum number of inputs allowed is 10 for now.
!-----------------------------------------------------------------------
   SUBROUTINE interpToFinestGrid_argument( sp1,sp2,sp3,sp4,sp5, &
                                           sp6,sp7,sp8,sp9,sp10 )
!-----------------------------------------------------------------------
      type(CLBLM_Spectrum) ,OPTIONAL,TARGET ,intent(inout) :: sp1,sp2,sp3,sp4,sp5, &
                                                              sp6,sp7,sp8,sp9,sp10

      integer ,parameter :: MXARG = 10

      integer :: i,j,k,nSp
      type(CLBLM_SpectrumPointer) ,allocatable :: spArray(:)

      !--- Count the number of input spectrum
      nSp = 0
      if (present(sp1)) nSp = nSp + 1
      if (present(sp2)) nSp = nSp + 1
      if (present(sp3)) nSp = nSp + 1
      if (present(sp4)) nSp = nSp + 1
      if (present(sp5)) nSp = nSp + 1
      if (present(sp6)) nSp = nSp + 1
      if (present(sp7)) nSp = nSp + 1
      if (present(sp8)) nSp = nSp + 1
      if (present(sp9)) nSp = nSp + 1
      if (present(sp10)) nSp = nSp + 1

      !--- Allocate pointer array
      allocate( spArray(nSp) )
      do i=1,nSp; NULLIFY(spArray(i)%ptr); enddo

      !--- Pointing to input spectra
      i=0
      if (present(sp1)) then; i=i+1; spArray(i)%ptr => sp1; endif
      if (present(sp2)) then; i=i+1; spArray(i)%ptr => sp2; endif
      if (present(sp3)) then; i=i+1; spArray(i)%ptr => sp3; endif
      if (present(sp4)) then; i=i+1; spArray(i)%ptr => sp4; endif
      if (present(sp5)) then; i=i+1; spArray(i)%ptr => sp5; endif
      if (present(sp6)) then; i=i+1; spArray(i)%ptr => sp6; endif
      if (present(sp7)) then; i=i+1; spArray(i)%ptr => sp7; endif
      if (present(sp8)) then; i=i+1; spArray(i)%ptr => sp8; endif
      if (present(sp9)) then; i=i+1; spArray(i)%ptr => sp9; endif
      if (present(sp10)) then; i=i+1; spArray(i)%ptr => sp10; endif

      call interpToFinestGrid( spArray )

   END SUBROUTINE

!-----------------------------------------------------------------------
! Given an array of spectrum pointers, find the one with minimum DV and interpolate
! all the rest to the minimum-DV grid.
!-----------------------------------------------------------------------
   SUBROUTINE interpToFinestGrid_pointerArray( spArray )
!-----------------------------------------------------------------------
      type(CLBLM_SpectrumPointer) ,intent(inout) :: spArray(:) !an array of pointers to CLBLM_Spectrum objects.

      character(*), parameter :: routineName='interpToFinestGrid'

      integer              :: numSp, isp, loc(1), locMinDV
      real, allocatable    :: spDV(:)
      real                 :: DV_min
      real(r8)             :: V1_min
      integer              :: NLIM_min
      type(CLBLM_Spectrum) :: newSpect


      !--- number of input spectra
      numSp = size(spArray)

      !--- Check existence of each spectrum
      do isp = 1,numSp
         if (.not.associated(spArray(isp)%ptr)) CYCLE !if pointing to NULL, cycle.
         if (.not.allocated( spArray(isp)%ptr%spect ))  &
            STOP '--- '//routineName//'(): Spectrum not allocated.'
      enddo

      !--- Find the minimum DV and its location in spArray
      allocate( spDV(numSp) )
      spDV(:) = huge(1.0)
      do isp=1,numSp
         if (.not.associated(spArray(isp)%ptr)) CYCLE
         spDV(isp) = spArray(isp)%ptr%DV
      enddo
      loc = minloc(spDV)
      locMinDV = loc(1)  !location of the spectrum with minimum DV
      DV_min   = spArray(locMinDV)%ptr%DV   !minimum DV
      !print *, 'DV_min in Spectrum.f90', DV_min
      V1_min   = spArray(locMinDV)%ptr%V1
      NLIM_min = spArray(locMinDV)%ptr%NLIM


      !--- Interpolate to the finest DV
      do isp = 1,numSp
         if (.not.associated(spArray(isp)%ptr)) CYCLE

         !--- Check if the spectra are in same grid
         if (isp == locMinDV) CYCLE
         if (spectraHaveSameGrid( spArray(isp)%ptr, &
                                  spArray(locMinDV)%ptr )) CYCLE

         !--- Starting V1 must be same
         if ( abs(spArray(isp)%ptr%V1 - V1_min) >epsilon(V1_min) ) &
            STOP '--- '//routineName//'(): Starting wavenumbers need to be same.'


         !--- Interpolate to finer grid
         call CLBLM_Spectrum_init( newSpect, spArray(isp)%ptr%V1, &
                                             DV_min, &
                                             NLIM_min, &
                                             spArray(isp)%ptr%indV1 )
         call interpSpectrum( spArray(isp)%ptr, newSpect )
         call moveSpectrum( newSpect, spArray(isp)%ptr )

      enddo

      if (allocated(spDV)) deallocate(spDV)
   END SUBROUTINE




!-----------------------------------------------------------------------
! Given an array of spectra, find the one with minimum DV and interpolate
! all the rest to the minimum-DV grid.
!-----------------------------------------------------------------------
   SUBROUTINE interpToFinestGrid_array( spArray )
!-----------------------------------------------------------------------
      type(CLBLM_Spectrum) :: spArray(:)

      character(*), parameter :: routineName='interpToFinestGrid'

      integer              :: numSp, isp, loc(1), locMinDV
      real, allocatable    :: spDV(:)
      real                 :: DV_min
      real(r8)             :: V1_min
      integer              :: NLIM_min
      type(CLBLM_Spectrum) :: newSpect


      !--- number of input spectra
      numSp = size(spArray)

      !--- Check existence of each spectrum
      do isp = 1,numSp
         if (.not.allocated( spArray(isp)%spect ))  &
            STOP '--- '//routineName//'(): Spectrum not allocated.'
      enddo

      !--- Find the minimum DV and its location in spArray
      allocate( spDV(numSp) )
      do isp=1,numSp
         spDV(isp) = spArray(isp)%DV
      enddo
      loc = minloc(spDV)
      locMinDV = loc(1)  !location of the spectrum with minimum DV
      DV_min   = spArray(locMinDV)%DV   !minimum DV
      V1_min   = spArray(locMinDV)%V1
      NLIM_min = spArray(locMinDV)%NLIM


      !--- Interpolate to the finest DV
      do isp = 1,numSp

         !--- Check if the spectra are in same grid
         if (isp == locMinDV) CYCLE
         if (spectraHaveSameGrid( spArray(isp), &
                                  spArray(locMinDV) ) ) CYCLE

         !--- Starting V1 must be same
         if ( abs(spArray(isp)%V1 - V1_min) >0. ) &
            STOP '--- '//routineName//'(): Starting wavenumbers need to be same.'


         !--- Interpolate to finer grid
         call CLBLM_Spectrum_init( newSpect, spArray(isp)%V1, &
                                             DV_min, &
                                             NLIM_min, &
                                             spArray(isp)%indV1 )
         call interpSpectrum( spArray(isp), newSpect )
         call moveSpectrum( newSpect, spArray(isp) )

      enddo

      if (allocated(spDV)) deallocate(spDV)
   END SUBROUTINE

!-----------------------------------------------------------------------
!     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED
!     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE
!     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN
!     INCREMENTS OF DVR3
!-----------------------------------------------------------------------
      SUBROUTINE XINT( V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3 )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8, ONEPL,ONEMI
      IMPLICIT NONE !REAL*8           (V)

      real(r8) ,intent(in)    :: V1A, V2A
      real     ,intent(in)    :: DVA
      real     ,intent(in)    :: A(:)
      real     ,intent(in)    :: AFACT
      real(r8) ,intent(in)    :: VFT
      real     ,intent(in)    :: DVR3
      real     ,intent(inout) :: R3(:)
      integer  ,intent(in)    :: N1R3,N2R3

      INTEGER  :: I,     IHI, ILO,    J
      REAL     :: B,     B1,  B2,     C
      REAL     :: CONTI, P,   RECDVA
      REAL(r8) :: VI,    VJ


      RECDVA = 1./DVA
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI
      ILO = MAX(ILO,N1R3)
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI
      IHI = MIN(IHI,N2R3)

      DO 10 I = ILO, IHI
         VI = VFT+DVR3* REAL(I-1)
         J = (VI-V1A)*RECDVA+ONEPL
         VJ = V1A+DVA* REAL(J-1)
         P = RECDVA*(VI-VJ)
         C = (3.-2.*P)*P*P
         B = 0.5*P*(1.-P)
         B1 = B*(1.-P)
         B2 = B*P
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2
         R3(I) = R3(I)+CONTI*AFACT
   10 END DO

      RETURN
      END  SUBROUTINE

!     --------------------------------------------------------------
      SUBROUTINE pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      USE Module_ConstParam ,ONLY: r8=>kind_r8, ONEPL,ONEMI
      IMPLICIT NONE

      real(r8) ,intent(in)    :: v1ss, v2ss, v1abs
      real     ,intent(in)    :: dvabs
      integer  ,intent(in)    :: nptabs
      integer  ,intent(inout) :: ist, last

      REAL     :: v1abs_loc
      INTEGER  :: nbnd_v1c, nbnd_v2c

!   Set up needed variables for call to XINT
!   Output variables
!     v1abs_loc - wavenumber of first value to be processed in XINT
!     nptabs_loc - number of values to be processed in XINT
!     ist - index of first value to be processed in XINT

      nbnd_v1c =  2 +  (v1ss-v1abs)/dvabs + 1.e-5
      ist = max(1,nbnd_v1c)
      v1abs_loc = v1abs + dvabs * float(ist-1)

      nbnd_v2c = 1 + (v2ss-v1abs)/dvabs + 1.e-5
      last = min(nptabs,nbnd_v2c)

     RETURN
     END SUBROUTINE

!     --------------------------------------------------------------

END MODULE
