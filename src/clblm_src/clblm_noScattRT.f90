!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!


MODULE Module_noScattRT

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mode1_lookDn_mergeDn, &
             mode2_lookDn_mergeUp, &
             mode3_lookUp_mergeUP

CONTAINS !======================  MODULE CONTAINS ======================



!-----------------------------------------------------------------------
!
!!* Three types of path
!                     TOA   (Sun)
!                     ---    ---
!                     /      /
!                    /      /
! obsLevel          /   solPath
!   ---            /      /
!     \       downPath   /
!      \         /      /
!      path     /      /
!        \     /      /
!         \   /      /
!    ______\_/______/________
!
!
!  -------------------------- obsLevel = N
!         layer N-1
!  -------------------------- level N-1; level# start from 1
!
!-----------------------------------------------------------------------
   SUBROUTINE  mode1_lookDn_mergeDn( outCtrl, rtCtrl, fluxCtrl, input_scene,&
                                     paths, surf, solRadTOA, &
                                     odCtrl, spectGrid, dvCtrl, flxttdOut, &
                                     SfcRadd, SfcUFlx, totalRadUp, TxArray, &
                                     ODarray, NLTE, postCtrl, doPreBox)
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, FILLINT,FILLREAL
      USE Module_ConstParam  ,ONLY: EPS, PI, GWGO1, GWGD1, GWGD2, GWGT1, GWGT2, GWGT3
      USE Module_Scene        ,ONLY: CLBLM_Scene
      USE Module_Config      ,ONLY: CLBLM_Output_Ctrl, &
                                    CLBLM_OD_Ctrl, &
                                    CLBLM_DV_Ctrl, &
                                    CLBLM_RT_Ctrl, &
                                    CLBLM_FLUX_Ctrl, &
                                    CLBLM_Post_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid, &
                                    moveSpectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Path, &
                                    CLBLM_PathCompound, &
                                    CLBLM_Layer, &
                                    loadLayerFromPath
      USE Module_Scene       ,ONLY: CLBLM_Surface, &
                                    surfThmEmis
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      USE Module_ScanFilter  ,ONLY: clblm_SCANFN
      IMPLICIT NONE


      type(CLBLM_Output_Ctrl)    ,intent(in)    :: outCtrl
      type(CLBLM_Scene)          ,intent(in)    :: input_scene
      type(CLBLM_RT_Ctrl)        ,intent(in)    :: rtCtrl
      type(CLBLM_FLUX_Ctrl)       ,intent(in)    :: fluxCtrl           !input flux options
      type(CLBLM_PathCompound)   ,intent(in)    :: paths          !input paths
      type(CLBLM_Surface)        ,intent(in)    :: surf           !input surface properties
      type(CLBLM_Spectrum)       ,intent(inout) :: solRadTOA      !input solar radiance at TOA (=S0*umu0/pi). It may be interpolated to finer grid.
      type(CLBLM_OD_Ctrl)        ,intent(in)    :: odCtrl         !input ODLAY options
      type(CLBLM_SpectGrid)      ,intent(in)    :: spectGrid      !input spectral grid info.
      type(CLBLM_DV_Ctrl)        ,intent(in)    :: dvCtrl
      type(CLBLM_Spectrum)       ,intent(out)   :: totalRadUp     !output upwelling radiance at observer level
      type(CLBLM_Spectrum)       ,intent(out)   :: TxArray(:)     !output Tx profile
      type(CLBLM_Spectrum)       ,intent(out)   :: ODarray(:)     !output layer ODs
      logical                    ,intent(in)    :: NLTE
      type(CLBLM_Post_Ctrl)      ,intent(in)    :: postCtrl
      logical                    ,intent(out)   :: doPreBox
      real, dimension(:,:),allocatable ,intent(out)   :: flxttdOut  !output downwelling flux to write out in clblm_driver
      real, dimension(:,:),allocatable ,intent(out)   :: SfcRadd    !output surface downwelling radiance for input into upwelling flux calc
      real, dimension(:),allocatable ,intent(out)     :: SfcUFlx    !output downwelling surface flux to write out in clblm_driver

      character(*) ,parameter :: routineName = 'mode1_lookDn_mergeDn'

      integer              :: i,j,k,m, il, lyrNo, iv, initNLIM, TxLay
      real(r8)             :: V1,V2, V1B,V2B, V1out,V2out
      real                 :: DV, DVB,DVout
      integer              :: lyrLo, lyrHi, nLayers !Number of layers between observer and target
      integer              :: nLyr2          !Number of layers with IPATH=2
      type(CLBLM_Path)     :: vertpath       !input vertical path
      type(CLBLM_Path)     :: uppath         !input viewing path
      type(CLBLM_Path)     :: downPath       !input downwelling path
      type(CLBLM_Path)     :: solPath        !input solar path
      integer              :: obsLevel
      logical              :: belowObs
      type(CLBLM_Layer)    :: refLayer, upLayer, downLayer, solLayer
      logical              :: ThermalOn, SolarOn, ODonly, Flux_flag
      integer              :: linInTau
      logical              :: scalThermalPath  !if=.true. uppath and downpath ODs are scaled from reference path OD
      logical              :: scalSolarPath    !if=.true. solar path OD is scaled from viewing path OD
      integer              :: refPath
      real                 :: scalUpFac, scalDownFac, scalSolarFac

      type(CLBLM_Spectrum) :: mrgTxA       !merged total Tx to observer
      type(CLBLM_Spectrum) :: mrgRadUp     !merged total upward radiance
      type(CLBLM_Spectrum) :: mrgRadDn     !merged total downward radiance
      type(CLBLM_Spectrum) :: mergeRadDn1     !merged total downward radiance at an angle
      type(CLBLM_Spectrum) :: mergeRadDn2     !merged total downward radiance at an angle
      type(CLBLM_Spectrum) :: mergeRadDn3     !merged total downward radiance at an angle
      type(CLBLM_Spectrum) :: mergeRadUp1     !merged total upward radiance at an angle
      type(CLBLM_Spectrum) :: mergeRadUp2     !merged total upward radiance at an angle
      type(CLBLM_Spectrum) :: mergeRadUp3     !merged total upward radiance at an angle
      type(CLBLM_Spectrum) :: mergeTxA1     !merged TxA at an angle
      type(CLBLM_Spectrum) :: mergeTxA2     !merged TxA at an angle
      type(CLBLM_Spectrum) :: mergeTxA3     !merged TxA at an angle
      type(CLBLM_Spectrum) :: solMrgTx     !merged total solar beam transmittance
      !type(CLBLM_Spectrum) :: surfRad      !surface emission
      !type(CLBLM_Spectrum) :: surfRefl     !surface reflection
      type(CLBLM_Spectrum) :: thmRadUp
      type(CLBLM_Spectrum) :: lyrOD
      integer              :: boxcarSize
      ! Added following for flux calc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type(CLBLM_Spectrum) :: surfRad, surfEmis
      integer                      :: out, nout, nlev, ndl, ilev, l, nang, outinrat, icount, iout, istart, kk, nlim, ii,iAng
      real                         :: dv_flux, factor, scalODfac, fsum
      real, dimension(:,:), allocatable :: radd, sradd, flxttd, radu, txa, sfc_sradd
      real, dimension(:), allocatable :: bound, secants, dfluxdv

      vertPath        = paths%vert
      uppath          = paths%view
      downPath        = paths%down
      solPath         = paths%sun
      obsLevel        = paths%obsLev
      scalThermalPath = paths%AMScalingThermal
      scalSolarPath   = paths%AMScalingSolar
      refPath         = paths%refPath

      ! Inputs for flux calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      flux_flag = fluxCtrl%flux_flag
      nang = fluxCtrl%nang
      DV_flux = fluxCtrl%DV_flux
      scalODfac=1.0

      ! Assign secants of the quadrature angles based on cosine values from Radsum (Clough et al. 1992 Table 2)
      allocate (secants(nang))
      !     For one angle (cosine = 0.66666666667)
      if (nang .EQ. 1) then
         secants(1)=1./0.66666666667
      endif
      !     For two angles (cosines are 0.84494897 and 0.35505103)
      if (nang .EQ. 2) then
         secants(1)=1.0/0.84494897
         secants(2)=1.0/0.35505103
      endif
      !     For three angles  (0.91141204,0.59053314,0.21234054)
      if (nang .EQ. 3) then
         secants(1)=1.0/0.91141204
         secants(2)=1.0/0.59053314
         secants(3)=1.0/0.21234054
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ThermalOn = rtCtrl%ThermalOn
      SolarOn   = rtCtrl%SolarOn
      linInTau  = rtctrl%linInTau

      if (outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD==0) &
         STOP '--- '//routineName//'(): No output requested.'

      if (outCtrl%Rad/=0 .and. .not.ThermalOn .and. .not.SolarOn) &
         STOP '--- '//routineName//'(): Request radiance output but no sources.'


      ODonly = (outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD/=0)
      !--- Find layout of IPATH
      !call findHiLoLayers( path%IPATH, lyrLo, lyrHi, nLayers, nLyr2 )
      lyrLo = 1
      lyrHi = max( uppath%nLay, downPath%nLay, solPath%nLay )

      !--- Initialization
      !
      V1 = spectGrid%V1
      V2 = spectGrid%V2
      initNLIM = 2 !ceiling( (spectGrid%V2 - spectGrid%V1) / spectGrid%DV(lyrHi) + 1. )

      !Check consistency of input for flux calc.  NOUT is number of output groups.
      OUT = (V2 - V1)/(DV_flux)
      NOUT = INT (OUT + EPS)
      IF (ABS(FLOAT(NOUT)-OUT) .GT. EPS) THEN
         STOP 'V1, V2, (OUT DV)/(IN DV)  ARE INCONSISTENT'
      ENDIF

      ! Assign some variables for flux calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Srad stores the radiances summed over a user-supplied spectral bin for flux calculations, for each angle
      allocate (sradd(nout+1,nang))
      allocate (sfc_sradd(nout+1,nang))
      allocate (bound(nout+1))
      bound(:)=0.
      ! Flxtt store the fluxes at each level, on the user-supplied spectral resolution
      allocate (flxttd((lyrhi-lyrlo)+1,nout+1))
      FLXTTD(:,:) = 0.0
      allocate (flxttdout((lyrhi-lyrlo)+1,nout+1))
      FLXTTDout(:,:) = 0.0
      allocate (SfcUFlx(nout+1))
      SfcUFlx(:) = 0.0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ThermalOn .and. SolarOn) then

         call CLBLM_Spectrum_init( mrgTxA,   V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( mrgRadUp, V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( mrgRadDn, V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( solMrgTx, V1,(V2-V1),initNLIM )
         mrgTxA%spect(:) =1.
         mrgRadUp%spect(:) =0.
         mrgRadDn%spect(:) =0.
         solMrgTx%spect(:) =1.

      elseif (ThermalOn) then

         call CLBLM_Spectrum_init( mrgTxA,   V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( mrgRadUp, V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( mrgRadDn, V1,(V2-V1),initNLIM )
         mrgTxA%spect(:) =1.
         mrgRadUp%spect(:) =0.
         mrgRadDn%spect(:) =0.

      elseif (SolarOn) then

         call CLBLM_Spectrum_init( mrgTxA,   V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( mrgRadUp, V1,(V2-V1),initNLIM )
         call CLBLM_Spectrum_init( solMrgTx, V1,(V2-V1),initNLIM )
         mrgTxA%spect(:) =1.
         mrgRadUp%spect(:) =0.
         solMrgTx%spect(:) =1.

      else !TxOnly

         if (.not.ODonly) then
            call CLBLM_Spectrum_init( mrgTxA,   V1,(V2-V1),initNLIM )
            mrgTxA%spect(:) =1.
         endif
      endif

      !--- Downward layer merging
      !
      do il = lyrHi, lyrLo, -1

         lyrNo = il
         belowObs = (lyrNo < obsLevel)

         !--- Load layer temperature, pressure and concentration ...
         ! * lyrNo must be the actual subscript of the path layer array.
         if (scalThermalPath .or. scalSolarPath) then
            if     (refPath==1) then; call loadLayerFromPath( refLayer, upPath, lyrNo );
            elseif (refPath==0) then; call loadLayerFromPath( refLayer, vertPath, lyrNo );
            endif
         endif

         if (                .not.scalThermalPath) call loadLayerFromPath( upLayer,   uppath, lyrNo )
         if (ThermalOn .and. .not.scalThermalPath) call loadLayerFromPath( downLayer, downPath, lyrNo )
         if (SolarOn   .and. .not.scalSolarPath)   call loadLayerFromPath( solLayer,  solPath, lyrNo )

         !--- Calculate scaling factors
         if (scalThermalPath) scalUpFac    = upPath%Wtot(lyrNo)   / refLayer%Wtot
         if (scalThermalPath) scalDownFac  = downPath%Wtot(lyrNo) / refLayer%Wtot
         if (scalSolarPath)   scalSolarFac = solPath%Wtot(lyrNo)  / refLayer%Wtot
         ! Setup flux calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! reset theses flux counters for each level
         ISTART = 1
         kk = 0
         ICOUNT = 0
         !        Compute width of output groups and wavenumbers of boundaries of flux output groups.
         !         SRADU (the radiance summed for flux) gets reset at each level
         DO K = 1, NOUT
            BOUND(K) = V1 + dv_flux * FLOAT(K-1)
            DO iAng = 1, NANG
               SRADD(K,iAng) = 0.0
               SFC_SRADD(K,iAng)= 0.0
            ENDDO
         ENDDO
         BOUND(NOUT+1) = V2

         IF (flux_flag .eqv. .true.) then
            IOUT = 1
            !loop over the angles

            DO iAng = 1, NANG
               !scalODfac will be used in the subroutine layerMerge_mode to scale the OD by the secant of the quadrature angle
               scalODfac=secants(iAng)

               IF (il .LT. lyrHi) then
                  ! set the incoming radiance to the radiance at the given quadrature angle
                  mrgRadDn%spect(:)=RADD(:,iAng)
                  mrgRadUp%spect(:)=RADU(:,iAng)
                  mrgTxA%spect(:)=TxA(:,iAng)
               ENDIF
               call layerMerge_mode1( refLayer, upLayer, downLayer, solLayer, &
                                odCtrl, spectGrid, dvCtrl, &
                                ODonly, ThermalOn, SolarOn, &
                                scalThermalPath, scalSolarPath, &
                                scalUpFac, scalDownFac, scalSolarFac, &
                                linInTau, NLTE, belowObs, &
                                mrgTxA, mrgRadUp, mrgRadDn, solMrgTx, lyrOD, scalODfac )
               ! Assign some variables for flux calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF (il .EQ. lyrHi .and. iAng .EQ. 1) then
                  ! factor is used to multiply radiance in W/cm2 ster to obtain fluxes in W/m2.
                  FACTOR = mrgRadDn%DV * 1.E04 * 2. * PI
                  outinrat=dv_flux/mrgRadDn%dv
                  ! Rad stores radiances for each angle, same spectral resolution as the radiances
                  allocate (radd(mrgRadDn%nlim,nang))
                  radd(:,:)=0.0
                  allocate (radu(mrgRadUp%nlim,nang))
                  radu(:,:)=0.0
                  allocate (txa(mrgTxA%nlim,nang))
                  txa(:,:)=0.0
                  ! Store downwelling surface flux on same grid radiances come in on
                  allocate (dfluxdv(mrgRadDn%nlim))
                  dfluxdv(:)=0.0
                  allocate (SfcRadd(mrgRadDn%nlim,nang))
                  SfcRadd(:,:)=0.0
               ENDIF ! ends allocating variables once the spectrum has been interpolated

               ! Keep track of the radiances and transmittance at this angle
               RADD(:,iAng)=mrgRadDn%spect(:)
               RADU(:,iAng)=mrgRadUp%spect(:)
               TxA(:,iAng)=mrgTxA%spect(:)

               ! keep track of the mrgRadDn structure, not just radiances,
               ! for the surface term from subroutine addThermalBoundary
               if (iAng .EQ. 1) then
                  mergeRadDn1=mrgRadDn
                  mergeRadUp1=mrgRadUp
                  mergeTxA1=mrgTxA
               endif
               if (iAng .EQ. 2) then
                  mergeRadDn2=mrgRadDn
                  mergeRadUp2=mrgRadUp
                  mergeTxA2=mrgTxA
               endif
               if (iAng .EQ. 3) then
                  mergeRadDn3=mrgRadDn
                  mergeRadUp3=mrgRadUp
                  mergeTxA3=mrgTxA
               endif

               !Before adding the radiances from the lowest downwelling layer they must be multiplied
               !by the surface reflectance then passed on to the first upwelling calculation;
               !once they have been passed they can then be summed and multiplied, assuming Lambertian reflection.
               !--- Add surface surface emission and reflected solar beam radiation
               if (il .EQ. lyrLo) then
                  ! SfcRadd stores surface radiances for each angle, same spectral resolution as the radiances
                  if (iAng .EQ. 1) then
                     call addThermalBoundary( thmRadUp, surf, mergeTxA1, mergeRadUp1, mergeRadDn1,fluxCtrl)
                     ! Keep track of the surface radiance at this angle
                     SfcRadd(:,iAng)=thmRadUp%spect(:)
                  endif
                  if (iAng .EQ. 2) then
                     call addThermalBoundary( thmRadUp, surf, mergeTxA2, mergeRadUp2,mergeRadDn2,fluxCtrl)
                     ! Keep track of the surface radiance at this angle
                     SfcRadd(:,iAng)=thmRadUp%spect(:)
                  endif
                  if (iAng .EQ. 3) then
                     call addThermalBoundary( thmRadUp, surf, mergeTxA3, mergeRadUp3, mergeRadDn3,fluxCtrl)
                     ! Keep track of the surface radiance at this angle
                     SfcRadd(:,iAng)=thmRadUp%spect(:)
                  endif
               endif ! ends surface radiance calculation

            ENDDO !ends loop over angles for radiance calculation at the given level

            !     Keep a running total of radiances in each desired output group.
            DO K = ISTART, mrgRadDn%nlim
               !  kk = kk + 1
               DO iAng = 1, NANG
                  SRADD(IOUT,iAng) = SRADD(IOUT,iAng) + RADD(K,iAng)
               ENDDO !ends loop over the quadrature angles
               ICOUNT = ICOUNT + 1
               IF (ICOUNT .GE. OUTINRAT) THEN
                  !           Current output group is complete.
                  ICOUNT = 0
                  IOUT = IOUT + 1
                  !    kk = 0
               ENDIF
               ! Save surface downwelling fluxes on the grid the radiances come in on.
               IF (il .eq. lyrlo) then
                  IF (NANG .EQ. 1) THEN
                     dfluxdv(k) = factor * GWGO1 * RADD(k,1)
                  ELSEIF (NANG .EQ. 2) THEN
                     dfluxdv(k) = factor * (GWGD1 * RADD(k,1) + GWGD2 * RADD(k,2))
                  ELSEIF (NANG .EQ. 3) THEN
                     dfluxdv(k) = factor * (GWGT1 * RADD(k,1) + GWGT2 * RADD(k,2)+ GWGT3 * RADD(k,3))
                  ENDIF
                  ! if (kk .eq. OUTINRAT) kk = 0
               ENDIF
            ENDDO ! ends K loop over the spectral range for adding the radiances

            ! All needed radiances have been summed for this level.  Time to
            !     calculate fluxes.
            ! QUADRATURE METHOD: First-order
            !     NDL = NLEV - ILEV
            DO L = 1, NOUT
               IF (NANG .EQ. 1) THEN
                  FLXTTD(IL,L) = GWGO1 * SRADD(L,1) * FACTOR
               ELSEIF (NANG .EQ. 2) THEN
                  FLXTTD(IL,L) = (GWGD1 * SRADD(L,1) + GWGD2 * SRADD(L,2)) * FACTOR
               ELSEIF (NANG .EQ. 3) THEN
                  FLXTTD(IL,L) = (GWGT1 * SRADD(L,1) + GWGT2 * SRADD(L,2)+ GWGT3 * SRADD(L,3)) * FACTOR
               ELSE
               STOP ' ERROR IN NANG '
               ENDIF
            ENDDO ! ends flux calculation loop over output groups

         ENDIF ! ends flux calculation for a given level

         scalODfac=1.0
         call layerMerge_mode1( refLayer, upLayer, downLayer, solLayer, &
                                odCtrl, spectGrid, dvCtrl, &
                                ODonly, ThermalOn, SolarOn, &
                                scalThermalPath, scalSolarPath, &
                                scalUpFac, scalDownFac, scalSolarFac, &
                                linInTau, NLTE, belowObs, &
                                mrgTxA, mrgRadUp, mrgRadDn, solMrgTx, lyrOD, scalODfac )

         !--- Calculate the spectral grid used for preboxing.
         !  Assuming the highest layer has finest DV.
         !
         if (il==lyrHi) then
            doPreBox = postCtrl%boxcarHW >0. .and. &
                       postCtrl%boxcarHW > mrgTxA%DV .and. &
                       postCtrl%functID /=FILLINT .and. &
                       (outCtrl%Rad <0 .or. outCtrl%Tx <0)


            if (doPreBox) then
               V1B = mrgTxA%V1
               V2B = mrgTxA%V2
               DVB = mrgTxA%DV
               boxcarSize = int( 2.0 * postCtrl%boxcarHW / DVB )
               V1out = V1B + (boxcarSize-1)*DVB/2.
               V2out = V2B - (boxcarSize-1)*DVB/2.
               DVout = (boxcarSize-1)*DVB
            endif
         endif

         !--- Save Tx
         TxLay = 0
         if (abs(outCtrl%Tx)==2 .and. belowObs) then !Save cumulative Tx
            TxArray(obsLevel-il) = mrgTxA
            TxLay = obsLevel-il
         elseif (abs(outCtrl%Tx)==1 .and. il==lyrLo) then !Save total Tx
            TxArray(1) = mrgTxA
            TxLay = 1
         endif

         !--- If needed do prebox smoothing
         if (doPreBox .and. TxLay >0 ) then
            call clblm_SCANFN( TxArray(TxLay), V1out,V2out,DVout, &
                               functID=1, functHWHM=postCtrl%boxcarHW )
         endif


         !--- Save OD
         if (outCtrl%OD/=0 .and. belowObs) then
            call moveSpectrum( lyrOD, ODarray(obsLevel-il) )
         endif
      enddo ! ends loop over layers

      if (flux_flag .eqv. .true.) then
         !     All incoming data have been processed.
         !     For each output group, calculate the surface upwelling flux
         !     by summing over the interval the sum of 1) the Planck function
         !     at the surface temperature times the emissivity and 2) the
         !     downwelling surface flux times the reflectivity.
         call surfThmEmis( surf, V1,mrgRadDn%dv,mrgRadDn%nlim, &
         surfRad, surfEmis )
         ICOUNT = 0
         ISTART=1
         IOUT=1
         DO K = ISTART, mrgRadDn%nlim-1
         !For Lambertian option, initialize the upwelling radiances at each of the three angles
         ! with the surface upwelling FLUX divided by PI, which would convert it to radiance
         ! then run the upwelling calculation as before.
              if (input_scene%sfc%ThmReflMode == 0) then !lambertian reflection
                 SfcRadd(k,:)=surfRad%spect(k) + ((1.-surfEmis%spect(k)) * (dfluxdv(k)/(PI*mrgRadDn%DV*1.E04)))
              endif
            SfcUFlx(IOUT) = SfcUFlx(IOUT) + surfRad%spect(k) * mrgRadDn%DV * 1.E04 * PI + (1.-surfEmis%spect(k)) * dfluxdv(k)
            ICOUNT = ICOUNT + 1
            IF (ICOUNT .GE. OUTINRAT) THEN
               !           Current output group is complete.
               ICOUNT = 0
               IOUT = IOUT + 1
            endif
         enddo
      endif
      !---Output Flux
      flxttdOut = flxttd
      !--- Add surface surface emission and reflected solar beam radiation
      !
      if (outCtrl%Rad/=0) then

         call addThermalBoundary( thmRadUp, surf, mrgTxA, mrgRadUp, mrgRadDn,fluxCtrl )

         if (SolarOn) then
            call addSolarBoundary( totalRadUp, surf, mrgTxA, thmRadUp, solMrgTx, solRadTOA )
         else
            totalRadUp = thmRadUp
         endif

         !--- If needed do prebox smoothing
         if (doPreBox) then
            call clblm_SCANFN( totalRadUp, V1out,V2out,DVout, &
                               functID=1, functHWHM=postCtrl%boxcarHW )
         endif

      endif

   END SUBROUTINE mode1_lookDn_mergeDn


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine layerMerge_mode1( refLayer, upLayer, downLayer, solLayer, &
                                odCtrl, spectGrid, dvCtrl, &
                                ODonly, ThermalOn, SolarOn, &
                                scalThermalPath, scalSolarPath, &
                                scalUpFac, scalDownFac, scalSolarFac, &
                                linInTau, NLTE, belowObs, &
                                mrgTxA, mrgRadUp, mrgRadDn, solMrgTx, lyrOD, scalODfac )
   !--------------------------------------------------------------------
      USE Module_Config      ,ONLY: CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid, &
                                    moveSpectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Layer
      USE Module_ODLAY       ,ONLY: ODLAY
      USE Module_EMLAY       ,ONLY: layerEmis
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      IMPLICIT NONE

      type(CLBLM_Layer)         ,intent(in)    :: refLayer
      type(CLBLM_Layer)         ,intent(in)    :: upLayer
      type(CLBLM_Layer)         ,intent(in)    :: downLayer
      type(CLBLM_Layer)         ,intent(in)    :: solLayer
      type(CLBLM_OD_Ctrl)       ,intent(in)    :: odCtrl        !input ODLAY options
      type(CLBLM_SpectGrid)     ,intent(in)    :: spectGrid     !input spectral grid info.
      type(CLBLM_DV_Ctrl)       ,intent(in)    :: dvCtrl
      logical                   ,intent(in)    :: ODonly        !Flag to indicate OD only calculation
      logical                   ,intent(in)    :: ThermalOn
      logical                   ,intent(in)    :: SolarOn
      logical                   ,intent(in)    :: scalThermalPath, scalSolarPath
      real                      ,intent(in)    :: scalUpFac, scalDownFac, scalSolarFac
      integer                   ,intent(in)    :: linInTau      ! =0 linear-in-tau not used; =1 standard; =2 LBLRTM version
      logical                   ,intent(in)    :: NLTE
      logical                   ,intent(in)    :: belowObs
      type(CLBLM_Spectrum)      ,intent(inout) :: mrgTxA        !merged total Tx to observer
      type(CLBLM_Spectrum)      ,intent(inout) :: mrgRadUp      !merged upwelling radiance
      type(CLBLM_Spectrum)      ,intent(inout) :: mrgRadDn      !merged downwelling raidance
      type(CLBLM_Spectrum)      ,intent(inout) :: solMrgTx      !merged solar beam transmittance
      type(CLBLM_Spectrum)      ,intent(out)   :: lyrOD         !output layer ODs
      ! Added scaling factor for quadrature angles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real                      ,intent(in)    :: scalODfac

      integer              :: i,j,k, iv
      real                 :: Ttop, Tbot, Tave
      type(CLBLM_Spectrum) :: refPathOD, refPathNLTEFac !layer OD along reference path
      type(CLBLM_Spectrum) :: lyrODA          !layer OD, up path
      type(CLBLM_Spectrum) :: lyrODB          !layer OD, down path
      type(CLBLM_Spectrum) :: lyrTxA          !layer Tx corresponding to layer OD.
      type(CLBLM_Spectrum) :: lyrTxB          !layer Tx, if Lambertian surface lyrTxB is the Tx at diffusivity angel, otherwise it equals to lyrTxA
      type(CLBLM_Spectrum) :: lyrEmUp         !layer upward emission
      type(CLBLM_Spectrum) :: lyrEmDn         !layer downward emission, if Lambertian surface, it is emission at diffusivity angle
      type(CLBLM_Spectrum) :: solLyrOD        !layer OD along solar path
      type(CLBLM_Spectrum) :: solLyrTx        !layer Tx along solar path
      type(CLBLM_Spectrum) :: nlteEmisFacUp   !NLTE parameters for emission calculation
      type(CLBLM_Spectrum) :: nlteEmisFacDn   !NLTE parameters for emission calculation
      type(CLBLM_Spectrum) :: dummySpect      !dummy argument, not used



      !--- Calculate reference path OD, if necessary
      if ( (scalThermalPath .and. ThermalOn) .or. &
           (scalSolarPath .and. SolarOn) .or. &
           (scalThermalPath .and. belowObs) .or. &
           (scalThermalPath .and. ODonly) ) then

        call ODLAY( refLayer, odCtrl, spectGrid, dvCtrl, &
                    refPathOD, NLTE, refPathNLTEFac )
      endif


      !--- Calculate up path OD
      if ( belowObs .or. ODonly ) then

         if ( scalThermalPath ) then
            call CLBLM_Spectrum_init( lyrODA, refPathOD%V1, &
                                              refPathOD%DV, &
                                              refPathOD%NLIM )
            do iv = lyrODA%indV1, lyrODA%indV2
               lyrODA%spect(iv) = refPathOD%spect(iv) * scalUpFac
            enddo

            ! For flux calculations, scale the optical depth by secant of the quadrature angle
            ! Otherwise scalODfac=1
            do iv = lyrODA%indV1, lyrODA%indV2
               lyrODA%spect(iv) = lyrODA%spect(iv) * scalODfac
            enddo

            if (NLTE) then
               call CLBLM_Spectrum_init( nlteEmisFacUp, refPathOD%V1, &
                                                        refPathOD%DV, &
                                                        refPathOD%NLIM )
               do iv = nlteEmisFacUp%indV1, nlteEmisFacUp%indV2
                  nlteEmisFacUp%spect(iv) = refPathNLTEFac%spect(iv) * scalUpFac
               enddo
            endif
         else
            call ODLAY( upLayer, odCtrl, spectGrid, dvCtrl, &
                        lyrODA, NLTE,nlteEmisFacUp )
         endif
      endif !( belowObs .or. ODonly )


      !--- OD for down welling path
      if (ThermalOn) then

         if (scalThermalPath) then !scale down path OD from up path OD

            call CLBLM_Spectrum_init( lyrODB, refPathOD%V1, &
                                              refPathOD%DV, &
                                              refPathOD%NLIM )
            do iv = lyrODB%indV1, lyrODB%indV2
               lyrODB%spect(iv) = refPathOD%spect(iv) * scalDownFac
            enddo

            ! For flux calculations, scale the optical depth by secant of the quadrature angle "scalODfac"
            ! Otherwise scalODfac=1
            do iv = lyrODB%indV1, lyrODB%indV2
               lyrODB%spect(iv) = lyrODB%spect(iv) * scalODfac
            enddo

            if (NLTE) then
               call CLBLM_Spectrum_init( nlteEmisFacDn, refPathOD%V1, &
                                                        refPathOD%DV, &
                                                        refPathOD%NLIM )
               do iv = nlteEmisFacDn%indV1, nlteEmisFacDn%indV2
                  nlteEmisFacDn%spect(iv) = refPathNLTEFac%spect(iv) * scalDownFac
               enddo
            endif
         else
            call ODLAY( downLayer, odCtrl, spectGrid, dvCtrl, &
                        lyrODB, NLTE,nlteEmisFacDn )
         endif
      endif !if (ThermalOn) then



      !--- Calculate transmittance for upwelling path
      if (belowObs .and. .not.ODonly) then

         call CLBLM_Spectrum_init( lyrTxA, lyrODA%V1, &
                                           lyrODA%DV, &
                                           lyrODA%NLIM )
         do iv = lyrTxA%indV1, lyrTxA%indV2
            lyrTxA%spect(iv) = exp( -lyrODA%spect(iv) )
         enddo
      endif


      if (ThermalOn) then
      ! * Calculate Tx along downwelling path.
      ! * Calculate layer downwelling thermal radiance
      ! * Merge downwelling radiance
      ! * Calculate layer upwelling thermal radiance
      ! * Merge upwelling radiance
      ! * If belowObs, only downwelling radiance is needed.

         !--- Transmittance for thermal down path
         call CLBLM_Spectrum_init( lyrTxB, lyrODB%V1, &
                                           lyrODB%DV, &
                                           lyrODB%NLIM )
         do iv = lyrTxB%indV1, lyrTxB%indV2
            lyrTxB%spect(iv) = exp( -lyrODB%spect(iv) )
         enddo


         !--- Downwelling thermal emission
         if (scalThermalPath) then
            Ttop = refLayer%Ttop
            Tbot = refLayer%Tbot
            Tave = refLayer%Tave
         else
            Ttop = downLayer%Ttop
            Tbot = downLayer%Tbot
            Tave = downLayer%Tave
         endif

         call layerEmis( lyrEmDn, Tbot, Ttop, Tave,&
                                  lyrTxB, lyrODB, linInTau, &
                                  NLTE, nlteEmisFacDn )

         !--- Merge downwelling radiance
         !
         !call interpSpectrum( lyrEmDn, mrgRadDn%DV, mrgRadDn%V1, mrgRadDn%NLIM )
         !call interpSpectrum( lyrTxB,  mrgRadDn%DV, mrgRadDn%V1, mrgRadDn%NLIM )
         call interpToFinestGrid( lyrEmDn, lyrTxB, mrgRadDn )

         do iv = mrgRadDn%indV1, mrgRadDn%indV2
            mrgRadDn%spect(iv) = lyrEmDn%spect(iv) + &
                                 mrgRadDn%spect(iv) * lyrTxB%spect(iv)
         enddo


         if (belowObs) then
            !--- Upwelling thermal emission
            if (linInTau/=0) then
               if (scalThermalPath) then
                  Ttop = refLayer%Ttop
                  Tbot = refLayer%Tbot
                  Tave = refLayer%Tave
               else
                  Ttop = upLayer%Ttop
                  Tbot = upLayer%Tbot
                  Tave = upLayer%Tave
               endif
               call layerEmis( lyrEmUp, Ttop, Tbot, Tave,&
                                        lyrTxA, lyrODA, linInTau, &
                                        NLTE, nlteEmisFacUp )
            else
               lyrEmUp = lyrEmDn !what about NLTE factor?
            endif


            !--- Merge upwelling radiance
            !
            !call interpSpectrum( lyrEmUp, mrgTxA%DV, mrgTxA%V1, mrgTxA%NLIM )
            call interpToFinestGrid( lyrEmUp, mrgTxA, mrgRadUp )

            do iv = mrgRadUp%indV1, mrgRadUp%indV2
               mrgRadUp%spect(iv) = mrgRadUp%spect(iv) + &
                                    mrgTxA%spect(iv) * lyrEmUp%spect(iv)
            enddo
         endif

      endif !if (ThermalOn) then


      if (belowObs .and. .not.ODonly) then
      !--- Merge transmittance for upward radiance

         !call interpSpectrum( lyrTxA, mrgTxA%DV, mrgTxA%V1, mrgTxA%NLIM )
         call interpToFinestGrid( lyrTxA, mrgTxA )

         do iv = mrgTxA%indV1, mrgTxA%indV2
            mrgTxA%spect(iv) = mrgTxA%spect(iv) * lyrTxA%spect(iv)
         enddo
      endif


      if (SolarOn) then

         !--- Transmittacne for solar path
         if ( scalSolarPath ) then  !solar path OD scaled from up path OD
            call CLBLM_Spectrum_init( solLyrTx, refPathOD%V1, &
                                                refPathOD%DV, &
                                                refPathOD%NLIM )
            do iv = solLyrTx%indV1, solLyrTx%indV2
               solLyrTx%spect(iv) = exp( -refPathOD%spect(iv) * scalSolarFac )
            enddo

         else !not scaling approximation

            call ODLAY( solLayer, odCtrl, spectGrid, dvCtrl, &
                        solLyrOD, NLTE,dummySpect )

            call CLBLM_Spectrum_init( solLyrTx, solLyrOD%V1, &
                                                solLyrOD%DV, &
                                                solLyrOD%NLIM )
            do iv = solLyrTx%indV1, solLyrTx%indV2
               solLyrTx%spect(iv) = exp( -solLyrOD%spect(iv) )
            enddo

         endif

         !--- Merge solar beam transmittance
         !
         !call interpSpectrum(solLyrTx, solMrgTx%DV,solMrgTx%V1,solMrgTx%NLIM)
         call interpToFinestGrid( solLyrTx, solMrgTx )

         do iv = solMrgTx%indV1, solMrgTx%indV2
            solMrgTx%spect(iv) = solMrgTx%spect(iv) * solLyrTx%spect(iv)
         enddo

      endif !if (SolarOn) then

      !--- Output layer OD
      call moveSpectrum( lyrODA, lyrOD )

   end subroutine



!-----------------------------------------------------------------------
! * Calculate total downwelling radiance and cumulative transmittance by
!   merging upward.
! * Cumulative transmittances are output in uniform DV.
!
!!* Viewing path
!
!                     TOA/(Sun)
!                     ---
!                     /
!                    /
!                   /
!                  /
!                path
!                /
!               /
!              /
!             /
!        ____/___
!        observer
!
!  -------------------------- N+1
!         layer N
!  -------------------------- obsLevel = N, N start from 1
!
!-----------------------------------------------------------------------
      SUBROUTINE mode3_lookUp_mergeUp( outCtrl, rtCtrl, &
                                       path, obsLevel, solRadTOA, &
                                       odCtrl, spectGrid, dvCtrl, &
                                       totalRadDn, TxArray, ODarray, NLTE, postCtrl, doPreBox )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, FILLINT
      USE Module_Config      ,ONLY: CLBLM_Output_Ctrl, &
                                    CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl, &
                                    CLBLM_RT_Ctrl, &
                                    CLBLM_Post_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    CLBLM_SpectrumPointer, &
                                    interpSpectrum, &
                                    interpToFinestGrid, &
                                    moveSpectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Path, &
                                    CLBLM_Layer, &
                                    loadLayerFromPath
      USE Module_ODLAY       ,ONLY: ODLAY
      USE Module_EMLAY       ,ONLY: layerEmis
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      USE Module_ScanFilter  ,ONLY: clblm_SCANFN
      IMPLICIT NONE


      type(CLBLM_Output_Ctrl)     ,intent(in)    :: outCtrl
      type(CLBLM_RT_Ctrl)         ,intent(in)    :: rtCtrl
      type(CLBLM_Path)            ,intent(in)    :: path          !input viewing path
      integer                     ,intent(in)    :: obsLevel
      type(CLBLM_Spectrum)        ,intent(inout) :: solRadTOA     !optional input solar radiance at TOA (=S0*umu0/pi)
      type(CLBLM_OD_Ctrl)         ,intent(in)    :: odCtrl        !input ODLAY options
      type(CLBLM_SpectGrid)       ,intent(in)    :: spectGrid     !input spectral grid info.
      type(CLBLM_DV_Ctrl)         ,intent(in)    :: dvCtrl
      type(CLBLM_Spectrum)        ,intent(out)   :: totalRadDn    !output upwelling radiance at observer level
      type(CLBLM_Spectrum)        ,intent(out)   :: TxArray(:)    !output Tx profile
      type(CLBLM_Spectrum)        ,intent(out)   :: ODarray(:)    !output layer OD
      logical                     ,intent(in)    :: NLTE
      type(CLBLM_Post_Ctrl)       ,intent(in)    :: postCtrl
      logical                     ,intent(out)   :: doPreBox

      character(*) ,parameter :: routineName = 'mode3_lookUp_mergeUp'

      integer              :: i,j,k, il, lyrNo, iv, initNLIM, l1,l2, Lmin
      real(r8)             :: V1,V2,V1B,V2B,V1out,V2out
      real                 :: DVB,DVout
      integer              :: lyrLo, lyrHi, nLayers !Number of layers between observer and target
      integer              :: nLyr2   !Number of layers with IPATH=2
      type(CLBLM_Layer)    :: aLayer
      logical              :: ThermalOn, SolarOn, ODonly
      integer              :: linInTau
      real                 :: DVmin
      type(CLBLM_Spectrum) :: lyrOD            !layer OD
      type(CLBLM_Spectrum) :: lyrTx            !layer Tx corresponding to layer OD.
      type(CLBLM_Spectrum) :: lyrEmDn          !layer downward emission, if Lambertian surface, it is emission at diffusivity angle
      type(CLBLM_Spectrum) :: mrgTx            !merged total Tx to observer
      type(CLBLM_Spectrum) :: mrgRadDn         !merged total downward radiance
      type(CLBLM_Spectrum) :: nlteEmisFac      !NLTE parameters for emission calculation
      integer              :: boxcarSize


      ThermalOn = rtCtrl%ThermalOn
      SolarOn   = rtCtrl%SolarOn
      linInTau  = rtCtrl%linInTau

      if (outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD==0) &
         STOP '--- '//routineName//'(): No output requested.'

      if (outCtrl%Rad/=0 .and. .not.ThermalOn .and. .not.SolarOn) &
         STOP '--- '//routineName//'(): Request radiance output but no sources specified.'


      ODonly = ( outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD/=0 )


      !--- Find layout of IPATH
      !call findHiLoLayers( path%IPATH, lyrLo, lyrHi, nLayers, nLyr2 )
      lyrLo = obsLevel !Since uplooking path contains layers from obsLevel to TOA only, obsLevel is actually =1
      lyrHi = path%nLay


      !--- Initializatoin
      !
      V1 = spectGrid%V1
      V2 = spectGrid%V2
      initNLIM = 2 !ceiling( (spectGrid%V2 - spectGrid%V1) / spectGrid%DV(lyrLo) + 1. )

      if (.not.ODonly) then
         call CLBLM_Spectrum_init( mrgTx, V1,(V2-V1),initNLIM )
         mrgTx%spect(:) = 1.0
      endif

      if (ThermalOn) then
         call CLBLM_Spectrum_init( mrgRadDn, V1,(V2-V1),initNLIM )
         mrgRadDn%spect(:) = 0.0
      endif



      !--- Merging upward
      !
      DO il = lyrLo,lyrHi

         lyrNo = il !In order to get the right value of DV, lyrNo must be the actual subscript of the path layer array.

         !--- Load layer temperature, pressure and concentration ...
         call loadLayerFromPath( aLayer, path, lyrNo )

         !--- Calculate OD
         call ODLAY( aLayer, odCtrl, spectGrid, dvCtrl, &
                     lyrOD, NLTE,nlteEmisFac )


         !--- Transmittance for down path
         if (.not.ODonly) then
            call CLBLM_Spectrum_init( lyrTx, lyrOD%V1, &
                                             lyrOD%DV, &
                                             lyrOD%NLIM )
            do iv = lyrTx%indV1, lyrTx%indV2
               lyrTx%spect(iv) = exp( -lyrOD%spect(iv) )
            enddo
         endif


         if (ThermalOn) then

            !--- Downwelling thermal emission
            call layerEmis( lyrEmDn, aLayer%Tbot, aLayer%Ttop, aLayer%Tave,&
                                     lyrTx, lyrOD, linInTau, &
                                     NLTE, nlteEmisFac )

            !--- Merge downwelling thermal
            !
            !call interpSpectrum( mrgTx,    lyrEmDn%DV,lyrEmDn%V1,lyrEmDn%NLIM)
            !call interpSpectrum( mrgRadDn, lyrEmDn%DV,lyrEmDn%V1,lyrEmDn%NLIM)
            call interpToFinestGrid( mrgTx, mrgRadDn, lyrEmDn )

            do iv = mrgRadDn%indV1, mrgRadDn%indV2
               mrgRadDn%spect(iv) = mrgRadDn%spect(iv) + &
                                    lyrEmDn%spect(iv) * mrgTx%spect(iv)
            enddo

         endif


         !--- Merge transmittance
         !
         if (.not.ODonly) then
            !call interpSpectrum( mrgTx, lyrTx%DV,lyrTx%V1,lyrTx%NLIM)
            call interpToFinestGrid( mrgTx, lyrTx )

            do iv = mrgTx%indV1, mrgTx%indV2
               mrgTx%spect(iv) = mrgTx%spect(iv) * lyrTx%spect(iv)
            enddo
         endif


         !--- Save Tx
         if ( abs(outCtrl%Tx) ==2 .and. il>obsLevel ) then !Save cumulative Tx
            TxArray(il-obsLevel) = mrgTx
         elseif ( abs(outCtrl%Tx)==1 .and. il==lyrHi ) then !Save total Tx
            TxArray(1) = mrgTx
         endif

         !--- Save OD
         if ( outCtrl%OD/=0 ) then
            call moveSpectrum( lyrOD, ODarray(il-obsLevel+1) )
         endif

      ENDDO !DO il = 1,size(mrgSeq)



      !--- Final merge
      !
      if (ThermalOn .and. SolarOn) then
         totalRadDn = mrgRadDn

         !call interpSpectrum( solRadTOA, mrgTx%DV, mrgTx%V1, mrgTx%NLIM )
         call interpToFinestGrid( solRadTOA, mrgTx, totalRadDn)

         do iv = totalRadDn%indV1, totalRadDn%indV2
            totalRadDn%spect(iv) = totalRadDn%spect(iv) + &
                                   mrgTx%spect(iv) * solRadTOA%spect(iv)
         enddo

      elseif (ThermalOn) then

         totalRadDn = mrgRadDn

      elseif (SolarOn) then

         call CLBLM_Spectrum_init( totalRadDn, mrgTx%V1, mrgTx%DV, mrgTx%NLIM )

         !call interpSpectrum( solRadTOA, mrgTx%DV, mrgTx%V1, mrgTx%NLIM )
         call interpToFinestGrid( solRadTOA, mrgTx )

         do iv = totalRadDn%indV1, totalRadDn%indV2
            totalRadDn%spect(iv) = mrgTx%spect(iv) * solRadTOA%spect(iv)
         enddo
      else
      endif



      !--- Regrid Tx to uniform grid and do prebox smoothing if needed
      !
      if (outCtrl%Tx /=0) then
         l1 = lbound(TxArray,1)
         l2 = ubound(TxArray,1)
         Lmin  = minloc( [(TxArray(i)%DV, i=l1,l2)], 1 )
         V1B = TxArray(Lmin)%V1
         V2B = TxArray(Lmin)%V2
         DVB = TxArray(Lmin)%DV
         DVmin = DVB
      elseif (outCtrl%Rad /=0) then
         V1B = totalRadDn%V1
         V2B = totalRadDn%V2
         DVB = totalRadDn%DV
         DVmin = DVB
      endif

      doPreBox = .FALSE.
      if ( outCtrl%Tx <0 .or. outCtrl%Rad <0 ) then
         doPreBox = postCtrl%boxcarHW > 0. .and. &
                    postCtrl%boxcarHW > DVmin .and. &
                    postCtrl%functID /=FILLINT
      endif

      if (doPreBox) then
         boxcarSize = int( 2.0 * postCtrl%boxcarHW / DVmin )
         V1out = V1B + (boxcarSize-1)*DVmin/2.
         V2out = V2B - (boxcarSize-1)*DVmin/2.
         DVout = (boxcarSize-1)*DVmin
      endif


      !--- cumTx will be output in uniform DV
      if ( abs(outCtrl%Tx)==2 ) then

         l1 = lbound(TxArray,1)
         l2 = ubound(TxArray,1)

         !DVmin = minval( [(TxArray(i)%DV, i=l1,l2)] )
         do il = l1,l2
            call interpSpectrum( TxArray(il), DVmin, &
                                              TxArray( l2 )%V1, &
                                              TxArray( l2 )%NLIM )

            if (doPreBox) then
               call clblm_SCANFN( TxArray(il), V1out,V2out,DVout, &
                                  functID=1, functHWHM=postCtrl%boxcarHW )
            endif

         enddo
      endif

      !--- Prebox smoothing total radiance
      if (outCtrl%Rad/=0 .and. doPreBox) then
         call clblm_SCANFN( totalRadDn, V1out,V2out,DVout, &
                            functID=1, functHWHM=postCtrl%boxcarHW )
      endif

   END SUBROUTINE mode3_lookUp_mergeUp


!-----------------------------------------------------------------------
!
!!* Three types of paths
!                     TOA   (Sun)
!                     ---    ---
!                     /      /
!                    /      /
! obsLevel          /   solPath
!   ---            /      /
!     \       downPath   /
!      \         /      /
!      path     /      /
!        \     /      /
!         \   /      /
!    ______\_/______/________
!
!
!
!  -------------------------- obsLevel = N
!         layer N-1
!  -------------------------- level N-1; level# start from 1
!
!-----------------------------------------------------------------------t
   SUBROUTINE mode2_lookDn_mergeUp( outCtrl, rtCtrl, fluxCtrl, &
                                    paths, surf, solRadTOA, &
                                    odCtrl, spectGrid, dvCtrl, flxttuOut, &
                                    SfcRad, totalRadUp, ODarray, NLTE, postCtrl, doPreBox )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, FILLINT
      USE Module_ConstParam  ,ONLY: EPS, PI, GWGO1, GWGD1, GWGD2, GWGT1, GWGT2, GWGT3
      USE Module_Config      ,ONLY: CLBLM_Output_Ctrl, &
                                    CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl, &
                                    CLBLM_RT_Ctrl, &
                                    CLBLM_FLUX_Ctrl, &
                                    CLBLM_Post_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid, &
                                    spectraHaveSameGrid, &
                                    moveSpectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Path, &
                                    CLBLM_PathCompound, &
                                    CLBLM_Layer, &
                                    loadLayerFromPath
      USE Module_Scene       ,ONLY: CLBLM_Surface, &
                                    surfThmEmis,&
                                    getSurfRefl
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      USE Module_ScanFilter  ,ONLY: clblm_SCANFN
      IMPLICIT NONE


      type(CLBLM_Output_Ctrl)     ,intent(in)    :: outCtrl
      type(CLBLM_RT_Ctrl)         ,intent(in)    :: rtCtrl
      type(CLBLM_FLUX_Ctrl)       ,intent(in)    :: fluxCtrl           !input flux options
      type(CLBLM_PathCompound)    ,intent(in)    :: paths
      type(CLBLM_Surface)         ,intent(in)    :: surf               !input surface properties
      type(CLBLM_Spectrum)        ,intent(inout) :: solRadTOA          !optional input solar radiance at TOA (=S0*umu0/pi)
      type(CLBLM_OD_Ctrl)         ,intent(in)    :: odCtrl             !input ODLAY options
      type(CLBLM_SpectGrid)       ,intent(in)    :: spectGrid          !input spectral grid info.
      type(CLBLM_DV_Ctrl)         ,intent(in)    :: dvCtrl
      type(CLBLM_Spectrum)        ,intent(out)   :: totalRadUp         !output upwelling radiance at observer level
      !type(CLBLM_Spectrum)        ,intent(out)   :: totalTx            !output total Tx
      type(CLBLM_Spectrum)        ,intent(out)   :: ODarray(:)         !output layer ODs
      logical                     ,intent(in)    :: NLTE
      type(CLBLM_Post_Ctrl)       ,intent(in)    :: postCtrl
      logical                     ,intent(out)   :: doPreBox
      real, dimension(:,:),allocatable ,intent(out)   :: flxttuOut   !output upwelling flux to write out in clblm_driver
      real, dimension(:,:)        ,intent(in)   :: SfcRad    !input surface downwelling radiance into upwelling flux calc


      character(*) ,parameter :: routineName = 'mode2_lookDn_mergeUp'


      integer                      :: i,j,k, lyrNo, il, iv, initNLIM
      real(r8)                     :: V1,V2, V1B,V2B, V1out,V2out
      real                         :: DV,DVB,DVout
      integer                      :: lyrLo, lyrHi, nLayers !Number of layers between observer and target
      integer                      :: nLyr2   !Number of layers with IPATH=2
      type(CLBLM_Path)             :: vertpath           !input vertical path
      type(CLBLM_Path)             :: uppath             !input viewing path
      type(CLBLM_Path)             :: downPath           !input downwelling path
      type(CLBLM_Path)             :: solPath            !optional input solar path
      integer                      :: obsLevel           !optional input of observer level that is lower than TOA
      logical                      :: belowObs
      type(CLBLM_Layer)            :: refLayer,upLayer,downLayer,solLayer
      logical                      :: ThermalOn, SolarOn, ODonly, flux_flag
      integer                      :: linInTau
      logical                      :: scalThermalPath    !if=.T. thermal path OD scaled from reference path OD
      logical                      :: scalSolarPath      !if=.T. solar path OD scaled from up path OD
      integer                      :: refPath
      real                         :: scalUpFac, scalDownFac, scalSolarFac
      type(CLBLM_Spectrum)         :: mrgRadUp     !merged total upward radiance
      type(CLBLM_Spectrum)         :: mrgThmRefl   !
      type(CLBLM_Spectrum)         :: solMrgRefl   !merged total solar beam transmittance
      type(CLBLM_SpectrumPointer)  :: spArray(10)
      type(CLBLM_Spectrum)         :: surfRad, surfEmis
      type(CLBLM_Spectrum)         :: thmRefl, solRefl
      type(CLBLM_Spectrum)         :: lyrOD
      integer                      :: boxcarSize
      ! Added following for flux calc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer                      :: out, nout, nlev, ndl, ilev, l, nang, outinrat, icount, iout, istart, kk, nlim, ii,iAng
      real                         :: dv_flux, factor, scalODfac
      real, dimension(:,:), allocatable :: radu, sradu, flxttu, thmref_ang
      real, dimension(:), allocatable :: bound, secants
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      vertpath        = paths%vert
      uppath          = paths%view
      downPath        = paths%down
      solPath         = paths%sun
      obsLevel        = paths%obsLev
      scalThermalPath = paths%AMScalingThermal
      scalSolarPath   = paths%AMScalingSolar
      refPath         = paths%refPath

      ! Inputs for flux calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      flux_flag = fluxCtrl%flux_flag
      nang = fluxCtrl%nang
      DV_flux = fluxCtrl%DV_flux

      scalODfac=1.0

      ! Assign secants of the quadrature angles based on cosine values from Radsum (Clough et al. 1992 Table 2)
      allocate (secants(nang))
      !     For one angle (cosine = 0.66666666667)
      if (nang .EQ. 1) then
         secants(1)=1./0.66666666667
      endif
      !     For two angles (cosines are 0.84494897 and 0.35505103)
      if (nang .EQ. 2) then
         secants(1)=1.0/0.84494897
         secants(2)=1.0/0.35505103
      endif
      !     For three angles  (0.91141204,0.59053314,0.21234054)
      if (nang .EQ. 3) then
         secants(1)=1.0/0.91141204
         secants(2)=1.0/0.59053314
         secants(3)=1.0/0.21234054
      endif

      ThermalOn = rtCtrl%ThermalOn
      SolarOn   = rtCtrl%SolarOn
      linInTau  = rtCtrl%linInTau

      if (outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD==0) &
         STOP '--- '//routineName//'(): No output requested.'

      if (outCtrl%Rad/=0 .and. .not.ThermalOn .and. .not.SolarOn) &
         STOP '--- '//routineName//'(): Request radiance output but no sources.'


      ODonly = (outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD/=0)


      !--- Find layout of IPATH
      !call findHiLoLayers( path%IPATH, lyrLo, lyrHi, nLayers, nLyr2 )
      lyrLo = 1
      lyrHi = max( uppath%nLay, downPath%nLay, solPath%nLay )

      !--- Initialzation

      V1 = spectGrid%V1
      V2 = spectGrid%V2
      DV = spectGrid%DVnormal(lyrLo)

      if (flux_flag .eqv. .true.) then
         DV = spectGrid%DVout
      endif

      initNLIM = ceiling( (V2-V1)/DV + 1. )

      if ( ThermalOn ) then
         call surfThmEmis( surf, V1,DV,initNLIM, &
                           surfRad, surfEmis )
         thmRefl = surfEmis !initialize the thmRefl structure
         thmRefl%spect(:) = 1.0-surfEmis%spect(:)
         ! For flux calculations, set thmRefl to zero, as it has already been accounted for in the surface radiances from the downwelling branch
         if (flux_flag .eqv. .true.) then
            thmRefl%spect(:) = 0.0
         endif
      endif

      if (SolarOn) then
         call getSurfRefl( surf, V1,DV,initNLIM, &
                           solRefl )
      endif

      if (ThermalOn .and. SolarOn) then
         mrgRadUp = surfRad
         mrgThmRefl = thmRefl
         solMrgRefl = solRefl
         elseif (ThermalOn) then
            mrgRadUp = surfRad
            mrgThmRefl = thmRefl
         elseif (SolarOn) then
            call CLBLM_Spectrum_init( mrgRadUp, V1,DV,initNLIM )
            mrgRadUp%spect(:) = 0.0
            solMrgRefl = solRefl
         else
      endif

      !Check consistency of input for flux calc.  NOUT is number of output groups.
      if (flux_flag .eqv. .true.) then
         OUT = (V2 - V1)/(DV_flux)
         NOUT = INT (OUT + EPS)
         IF (ABS(FLOAT(NOUT)-OUT) .GT. EPS) THEN
            STOP 'V1, V2, (OUT DV)/(IN DV)  ARE INCONSISTENT'
         endif
	   ENDIF


      ! Sradu stores the upwelling radiances summed over a user-supplied spectral bin for flux calculations, for each angle
      allocate (sradu(nout+1,nang))
      allocate (bound(nout+1))
      bound(:)=0.
      ! Flxttu store the upwelling fluxes at each level, on the user-supplied spectral resolution
      allocate (flxttu((lyrhi-lyrlo)+1,nout+1))
      FLXTTU(:,:) = 0.0
      allocate (flxttuout((lyrhi-lyrlo)+1,nout+1))
      FLXTTUout(:,:) = 0.0

      !--- Merging layers
      ! * Calculate total upwelling thermal radiance
      ! * Calculate merged solar reflectance
      do il = lyrLo,lyrHi

         lyrNo = il

         belowObs = (lyrNo < obsLevel)

         !--- Load reference layer
         if (scalThermalPath .or. scalSolarPath) then
            if     (refPath==1) then; call loadLayerFromPath( refLayer, upPath, lyrNo )
            elseif (refPath==0) then; call loadLayerFromPath( refLayer, vertPath, lyrNo )
            endif
         endif

         !--- Load layer temperature, pressure and concentration ...
         if (                .not.scalThermalPath) call loadLayerFromPath( upLayer,   upPath, lyrNo )
         if (ThermalOn .AND. .not.scalThermalPath) call loadLayerFromPath( downLayer, downPath, lyrNo )
         if (SolarOn   .AND. .not.scalSolarPath)   call loadLayerFromPath( solLayer,  solPath, lyrNo )

         !--- Scaling factors
         if (scalThermalPath) scalUpFac    = upPath%Wtot(lyrNo)   / refLayer%Wtot
         if (scalThermalPath) scalDownFac  = downPath%Wtot(lyrNo) / refLayer%Wtot
         if (scalSolarPath)   scalSolarFac = solPath%Wtot(lyrNo)  / refLayer%Wtot

         ! Setup flux calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! reset theses flux counters for each level
         ISTART = 1
         kk = 0
         ICOUNT = 0
         !        Compute width of output groups and wavenumbers of boundaries of flux output groups.
         !         SRADU (the radiance summed for flux) gets reset at each level
         if (flux_flag .eqv. .true.) then
            DO K = 1, NOUT
               BOUND(K) = V1 + dv_flux * FLOAT(K-1)
               DO iAng = 1, NANG
                  SRADU(K,iAng) = 0.0
               enddo
            enddo
         BOUND(NOUT+1) = V2
         endif

         if (flux_flag .eqv. .true.) then
            IOUT = 1
            !loop over the angles
            do iAng = 1, NANG
               !scalODfac will be used in the subroutine layerMerge_mode to scale the OD by the secant of the quadrature angle
               scalODfac=secants(iAng)

               ! Initialize mrgRadup with downwelling surface radiances
               if (il .EQ. lyrlo) then
                  mrgRadUp%spect(:)=SfcRad(:,iAng)
               endif
               if (il .GT. 1) then
                  ! set the incoming radiance to the radiance at the given quadrature angle
                  mrgRadUp%spect(:)=RADU(:,iAng)
                  mrgThmRefl%spect(:)=thmref_ang(:,iAng)
               endif
               call layerMerge_mode2( refLayer, upLayer, downLayer, solLayer, &
                                odCtrl, spectGrid, dvCtrl, &
                                ODonly, ThermalOn, SolarOn, &
                                scalThermalPath, scalSolarPath, &
                                scalUpFac, scalDownFac, scalSolarFac, &
                                linInTau, NLTE, belowObs, &
                                mrgRadUp, mrgThmRefl, solMrgRefl, lyrOD, scalODfac )

               ! Assign some variables for flux calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (il .EQ. lyrLo .and. iAng .EQ. 1) then
                  ! factor is used to multiply radiance in W/cm2 ster to obtain fluxes in W/m2.
                  FACTOR = mrgRadUp%DV * 1.E04 * 2. * PI
                  outinrat=dv_flux/mrgRadUp%dv
                  ! Rad stores radiances for each angle, same spectral resolution as the radiances
                  allocate (radu(mrgRadUp%nlim,nang))
                  radu(:,:)=0.0
                  allocate (thmref_ang(mrgRadUp%nlim,nang))
                  thmref_ang(:,:)=0.0
               endif ! ends allocating variables once the spectrum has been interpolated

               ! Keep track of the radiance at this angle
               RADU(:,iAng)=mrgRadUp%spect(:)
               thmref_ang(:,iAng)=mrgThmRefl%spect(:)

            ENDDO !ends loop over angles for radiance calulation at the given level

            !     Keep a running total of radiances in each desired output group.
            DO K = ISTART, mrgRadUp%nlim
               !  kk = kk + 1
               DO iAng = 1, NANG
                  SRADU(IOUT,iAng) = SRADU(IOUT,iAng) + Radu(K,iAng)
               enddo !ends loop over the quadrature angles
               ICOUNT = ICOUNT + 1
               IF (ICOUNT .GE. OUTINRAT) THEN
                  !           Current output group is complete.
                  ICOUNT = 0
                  IOUT = IOUT + 1
                  !    kk = 0
               ENDIF
            enddo ! ends K loop over the spectral range for adding the radiances

            ! All needed radiances have been summed for this level.  Time to
            !     calculate fluxes.
            DO L = 1, NOUT
               IF (NANG .EQ. 1) THEN
                  FLXTTU(IL,L) = GWGO1 * SRADU(L,1) * FACTOR
               ELSEIF (NANG .EQ. 2) THEN
                  FLXTTU(IL,L) = (GWGD1 * SRADU(L,1) + GWGD2 * SRADU(L,2)) * FACTOR
               ELSEIF (NANG .EQ. 3) THEN
                  FLXTTU(IL,L) = (GWGT1 * SRADU(L,1)+ GWGT2 * SRADU(L,2) + GWGT3 *SRADU(L,3)) * FACTOR
               ELSE
                  STOP ' ERROR IN NANG '
               ENDIF
            enddo ! ends flux calculation loop over output groups
         endif ! ends if statement for flux calculation at this level

         scalODfac=1.0
         call layerMerge_mode2( refLayer, upLayer, downLayer, solLayer, &
                                odCtrl, spectGrid, dvCtrl, &
                                ODonly, ThermalOn, SolarOn, &
                                scalThermalPath, scalSolarPath, &
                                scalUpFac, scalDownFac, scalSolarFac, &
                                linInTau, NLTE, belowObs, &
                                mrgRadUp, mrgThmRefl, solMrgRefl, lyrOD, scalODfac )

         !---Output OD
         if (outCtrl%OD/=0 .and. belowObs) then
            call moveSpectrum( lyrOD, ODarray(il) )
         endif

      enddo ! ends loop over layers

      !---Output Flux
      flxttuOut = flxttu

      if (SolarOn) then
      ! Add attenuated solar beam into total upwelling radiance

         call interpToFinestGrid( mrgRadUp, solMrgRefl, solRadTOA )

         do iv = mrgRadUp%indV1, mrgRadUp%indV2
            mrgRadUp%spect(iv) = mrgRadUp%spect(iv) + &
                                 solMrgRefl%spect(iv) * solRadTOA%spect(iv)
         enddo

      endif


      !---Output Rad
      if (outCtrl%Rad/=0) then

         totalRadUp = mrgRadUp

         doPreBox = postCtrl%boxcarHW > 0. .and. &
                    postCtrl%boxcarHW > totalRadUp%DV .and. &
                    postCtrl%functID /=FILLINT .and. &
                    outCtrl%Rad <0

         if (doPreBox) then
            V1B = totalRadUp%V1
            V2B = totalRadUp%V2
            DVB = totalRadUp%DV
            boxcarSize = int( 2.0 * postCtrl%boxcarHW / DVB )
            V1out = V1B + (boxcarSize-1)*DVB/2.
            V2out = V2B - (boxcarSize-1)*DVB/2.
            DVout = (boxcarSize-1)*DVB
            call clblm_SCANFN( totalRadUp, V1out,V2out,DVout, &
                               functID=1, functHWHM=postCtrl%boxcarHW )
         endif
      endif

   END SUBROUTINE mode2_lookDn_mergeUp


  !-----------------------------------------------------------------------
   SUBROUTINE writeSpectFile_2D_array( dataOut,  fileName, varName, dimName, fileHdr)
      !-----------------------------------------------------------------------
            USE netcdf
            USE Module_Spectrum  ,ONLY: CLBLM_Spectrum
            USE Module_FileIO    ,ONLY: createNetCDFFile, &
                                        checkNetCDFcall, &
                                        CLBLM_SpectFileHeader
            IMPLICIT NONE

            character(*) ,parameter :: routineName='writeSpectFile_2D_array'

            real                           ,intent(in) :: dataOut(:,:) ![nPoints, nSpectr]
            character(*)                   ,intent(in) :: fileName
            character(*)                   ,intent(in) :: varName
            character(*)                   ,intent(in) :: dimName(:)
            !type(CLBLM_SpectFileHeader)    ,intent(in) :: fileHdr
            character(*)                   ,intent(in) :: fileHdr

            integer(4) :: ncid, dim1ID,dim2ID, dim1,dim2, varID


            !--- Create output file
            if (.not. createNetCDFFile(trim(fileName), ncid)) then
               STOP '--- '//routineName//'(): cannot create netcdf file, this should not happen'
            endif
            call checkNetCDFcall( nf90_enddef(ncid) )

            !--- Write file header and radiance data
            !call writeSpectFileHeader(ncid, fileHdr)

            dim1 = size(dataOut,1)
            dim2 = size(dataOut,2)
            call checkNetCDFcall(nf90_def_dim( ncid, dimName(1), dim1, dim1ID))
            call checkNetCDFcall(nf90_def_dim( ncid, dimName(2), dim2, dim2ID))
            call checkNetCDFcall(nf90_def_var( ncid, varName, NF90_DOUBLE, [dim1ID,dim2ID], varID))

            call checkNetCDFcall(nf90_enddef( ncid))

            ! write the data
            call checkNEtCDFcall(nf90_put_var( ncid, varID, dataOut(:,:) ))

            call checkNetCDFcall(nf90_close( ncid))

         END SUBROUTINE



      !----------------------------------------------------------------------- !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine layerMerge_mode2( refLayer, upLayer, downLayer, solLayer, &
                                odCtrl, spectGrid, dvCtrl, &
                                ODonly, ThermalOn, SolarOn, &
                                scalThermalPath, scalSolarPath, &
                                scalUpFac, scalDownFac, scalSolarFac, &
                                linInTau, NLTE, belowObs, &
                                mrgRadUp, mrgThmRefl, solMrgRefl, lyrOD, scalODfac )
   !--------------------------------------------------------------------
      USE Module_Config      ,ONLY: CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    spectraHaveSameGrid, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid, &
                                    moveSpectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Layer
      USE Module_ODLAY       ,ONLY: ODLAY
      USE Module_EMLAY       ,ONLY: layerEmis
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      IMPLICIT NONE

      type(CLBLM_Layer)         ,intent(in)    :: refLayer
      type(CLBLM_Layer)         ,intent(in)    :: upLayer
      type(CLBLM_Layer)         ,intent(in)    :: downLayer
      type(CLBLM_Layer)         ,intent(in)    :: solLayer
      type(CLBLM_OD_Ctrl)       ,intent(in)    :: odCtrl       !input ODLAY options
      type(CLBLM_SpectGrid)     ,intent(in)    :: spectGrid     !input spectral grid info.
      type(CLBLM_DV_Ctrl)       ,intent(in)    :: dvCtrl
      logical                   ,intent(in)    :: ODonly
      logical                   ,intent(in)    :: ThermalOn
      logical                   ,intent(in)    :: SolarOn
      logical                   ,intent(in)    :: scalThermalPath, scalSolarPath
      real                      ,intent(in)    :: scalUpFac, scalDownFac, scalSolarFac
      integer                   ,intent(in)    :: linInTau      ! =0 linear-in-tau not used; =1 standard; =2 LBLRTM version
      logical                   ,intent(in)    :: NLTE
      logical                   ,intent(in)    :: belowObs
      type(CLBLM_Spectrum)      ,intent(inout) :: mrgRadUp      !merged total upwelling radiance
      type(CLBLM_Spectrum)      ,intent(inout) :: mrgThmRefl    !
      type(CLBLM_Spectrum)      ,intent(inout) :: solMrgRefl    !merged total solar beam reflectance
      type(CLBLM_Spectrum)      ,intent(out)   :: lyrOD
      ! Added scaling factor for quadrature angles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real                      ,intent(in)    :: scalODfac

      integer              :: i,j,k,iv
      real                 :: Ttop, Tbot, Tave
      type(CLBLM_Spectrum) :: refPathOD, refPathNLTEFac  !layer OD along reference path
      type(CLBLM_Spectrum) :: lyrODA          !layer OD, along up path
      type(CLBLM_Spectrum) :: lyrODB          !layer OD, along donw path
      type(CLBLM_Spectrum) :: lyrTxA          !layer Tx corresponding to layer OD.
      type(CLBLM_Spectrum) :: lyrTxB          !layer Tx, if Lambertian surface lyrTxB is the Tx at diffusivity angel, otherwise it equals to lyrTxA
      type(CLBLM_Spectrum) :: lyrEmUp         !layer upward emission
      type(CLBLM_Spectrum) :: lyrEmDn         !layer downward emission, if Lambertian surface, it is emission at diffusivity angle
      type(CLBLM_Spectrum) :: solLyrOD        !layer OD along solar path
      type(CLBLM_Spectrum) :: solLyrTx        !layer Tx along solar path
      type(CLBLM_Spectrum) :: nlteEmisFacUp   !NLTE parameters for emission calculation
      type(CLBLM_Spectrum) :: nlteEmisFacDn   !NLTE parameters for emission calculation
      type(CLBLM_Spectrum) :: dummySpect      !dummy argument, not used


      !--- Reference path OD
      if ( (scalThermalPath .and. ThermalOn) .or. &
           (scalSolarPath .and. SolarOn) .or. &
           (scalThermalPath .and. belowObs) .or. &
           (scalThermalPath .and. ODonly) ) then

         call ODLAY( refLayer, odCtrl, spectGrid, dvCtrl, &
                     refPathOD, NLTE, refPathNLTEFac )
      endif



      !--- Calculate OD for upwelling path
      if ( belowObs .or. ODonly ) then

         if (scalThermalPath) then
            call CLBLM_Spectrum_init( lyrODA, refPathOD%V1, &
                                              refPathOD%DV, &
                                              refPathOD%NLIM )
            do iv = lyrODA%indV1, lyrODA%indV2
               lyrODA%spect(iv) = refPathOD%spect(iv) * scalUpFac
            enddo

            ! For flux calculations, scale the optical depth by secant of the quadrature angle
            ! Otherwise scalODfac=1
            do iv = lyrODA%indV1, lyrODA%indV2
               lyrODA%spect(iv) = lyrODA%spect(iv) * scalODfac
            enddo

            if (NLTE) then
               call CLBLM_Spectrum_init( nlteEmisFacUp, refPathOD%V1, &
                                                        refPathOD%DV, &
                                                        refPathOD%NLIM )
               do iv = nlteEmisFacUp%indV1, nlteEmisFacUp%indV2
                  nlteEmisFacUp%spect(iv) = refPathNLTEFac%spect(iv) * scalUpFac
               enddo
            endif
         else

            call ODLAY( upLayer, odCtrl, spectGrid, dvCtrl, &
                        lyrODA, NLTE,nlteEmisFacUp )
         endif
      endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !--- OD for downwelling path
      if (ThermalOn) then ! .and. reflecting surface?

         if (scalThermalPath) then

            call CLBLM_Spectrum_init( lyrODB, refPathOD%V1, &
                                              refPathOD%DV, &
                                              refPathOD%NLIM )
            do iv = lyrODB%indV1, lyrODB%indV2
               lyrODB%spect(iv) = refPathOD%spect(iv) * scalDownFac
            enddo
            ! For flux calculations, scale the optical depth by secant of the quadrature angle
            ! Otherwise scalODfac=1
            do iv = lyrODB%indV1, lyrODB%indV2
               lyrODB%spect(iv) = lyrODB%spect(iv) * scalODfac
            enddo

            if (NLTE) then
               call CLBLM_Spectrum_init( nlteEmisFacDn, refPathOD%V1, &
                                                        refPathOD%DV, &
                                                        refPathOD%NLIM )
               do iv = nlteEmisFacDn%indV1, nlteEmisFacDn%indV2
                  nlteEmisFacDn%spect(iv) = refPathNLTEFac%spect(iv) * scalDownFac
               enddo
            endif

         else

            call ODLAY( downLayer, odCtrl, spectGrid, dvCtrl, &
                        lyrODB, NLTE,nlteEmisFacDn )
         endif
      endif !if (ThermalOn) then




      !--- Calculate transmittance for upwelling path
      if (belowObs .and. .not.ODonly) then

         call CLBLM_Spectrum_init( lyrTxA, lyrODA%V1, &
                                           lyrODA%DV, &
                                           lyrODA%NLIM )
         do iv = lyrTxA%indV1, lyrTxA%indV2
            lyrTxA%spect(iv) = exp( -lyrODA%spect(iv) )
         enddo
      endif


      if (ThermalOn)  then
      ! * Calculate Tx along downwelling path.
      ! * Calculate layer downwelling thermal radiance
      ! * Calculate layer upwelling thermal radiance
      ! * Add the layer emission to the previous upwelling radiance
      ! * Update the thermal reflectance.
      ! * If belowObs, only downwelling radiance is needed.


         !--- Transmittance for thermal down path
         call CLBLM_Spectrum_init( lyrTxB, lyrODB%V1, &
                                           lyrODB%DV, &
                                           lyrODB%NLIM )
         do iv = lyrTxB%indV1, lyrTxB%indV2
            lyrTxB%spect(iv) = exp( -lyrODB%spect(iv) )
         enddo


         !--- Downwelling thermal emission
         if (scalThermalPath) then
            Ttop = refLayer%Ttop
            Tbot = refLayer%Tbot
            Tave = refLayer%Tave
         else
            Ttop = downLayer%Ttop
            Tbot = downLayer%Tbot
            Tave = downLayer%Tave
         endif
         call layerEmis( lyrEmDn, Tbot, Ttop, Tave,&
                                  lyrTxB, lyrODB, linInTau, &
                                  NLTE, nlteEmisFacDn )


         !--- Upwelling thermal emission
         if (belowObs) then

            if (linInTau /=0) then !use linear-in-tau approximation
               if (scalThermalPath) then
                  Ttop = refLayer%Ttop
                  Tbot = refLayer%Tbot
                  Tave = refLayer%Tave
               else
                  Ttop = upLayer%Ttop
                  Tbot = upLayer%Tbot
                  Tave = upLayer%Tave
               endif
               call layerEmis( lyrEmUp, Ttop, Tbot, Tave,&
                                        lyrTxA, lyrODA, linInTau, &
                                        NLTE, nlteEmisFacUp )
            else !what about NLTE case?
               lyrEmUp = lyrEmDn
            endif
         endif


         !--- Merge upwelling radiance and update thermal reflectance
         !
         if (belowObs) then
            !call interpSpectrum( mrgRadUp,   lyrODA%DV,lyrODA%V1,lyrODA%NLIM )
            !call interpSpectrum( mrgThmRefl, lyrODA%DV,lyrODA%V1,lyrODA%NLIM )

            call interpToFinestGrid( mrgRadUp, mrgThmRefl, lyrTxA, &
                                     lyrTxB, lyrEmUp, lyrEmDn )

            !--- Merge upwelling radiance

            do iv = mrgRadUp%indV1, mrgRadUp%indV2
               mrgRadUp%spect(iv) = lyrEmUp%spect(iv) + &
                                      lyrTxA%spect(iv) * &
                                      ( mrgRadUp%spect(iv) + &
                                        lyrEmDn%spect(iv) * mrgThmRefl%spect(iv) )
            enddo

            !--- Update lower layer reflectance
            do iv = mrgThmRefl%indV1, mrgThmRefl%indV2
               mrgThmRefl%spect(iv) = lyrTxA%spect(iv) * &
                                      lyrTxB%spect(iv) * &
                                      mrgThmRefl%spect(iv)
            enddo

         else !layer is above observer

            call interpToFinestGrid( mrgRadUp, mrgThmRefl, lyrTxB, lyrEmDn)

             !--- Merge downwelling radiance
             do iv = mrgRadUp%indV1, mrgRadUp%indV2
                mrgRadUp%spect(iv) = mrgRadUp%spect(iv) + &
                                     lyrEmDn%spect(iv) * mrgThmRefl%spect(iv)
             enddo

             !--- Update thermal reflectance
             do iv = mrgThmRefl%indV1, mrgThmRefl%indV2
                mrgThmRefl%spect(iv) = lyrTxB%spect(iv) * mrgThmRefl%spect(iv)
             enddo

         endif !if (belowObs) then

      endif !if (ThermalOn)



      if (SolarOn) then
      ! * Update lower layer solar reflectance

         !--- Transmittacne for solar path
         if ( scalSolarPath ) then

            call CLBLM_Spectrum_init( solLyrTx, refPathOD%V1, &
                                                refPathOD%DV, &
                                                refPathOD%NLIM )
            do iv = solLyrTx%indV1, solLyrTx%indV2
               solLyrTx%spect(iv) = exp( -refPathOD%spect(iv) * scalSolarFac )
            enddo

         else !not airmass scaling approximation

            call ODLAY( solLayer, odCtrl, spectGrid, dvCtrl, &
                        solLyrOD, NLTE,dummySpect )

            call CLBLM_Spectrum_init( solLyrTx, solLyrOD%V1, &
                                                solLyrOD%DV, &
                                                solLyrOD%NLIM )
            do iv = solLyrTx%indV1, solLyrTx%indV2
               solLyrTx%spect(iv) = exp( -solLyrOD%spect(iv) )
            enddo
         endif


         !--- Update solar reflectance
         !
         if (belowObs) then

            call interpToFinestGrid( solMrgRefl, solLyrTx, lyrTxA)

            do iv = solMrgRefl%indV1, solMrgRefl%indV2
               solMrgRefl%spect(iv) = lyrTxA%spect(iv) * &
                                      solLyrTx%spect(iv) * &
                                      solMrgRefl%spect(iv)
            enddo
         else

            call interpToFinestGrid( solMrgRefl, solLyrTx )

            do iv = solMrgRefl%indV1, solMrgRefl%indV2
               solMrgRefl%spect(iv) = solLyrTx%spect(iv) * &
                                      solMrgRefl%spect(iv)
            enddo
         endif

      endif !if (SolarOn) then

      !---Output layer OD
      call moveSpectrum( lyrODA, lyrOD )

   end subroutine !layerMerge_mode2()


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE findHiLoLayers( IPATH, lyrLo, lyrHi, nLayers, nLyr2 )
!-----------------------------------------------------------------------
      integer    ,intent(in)  :: IPATH(:)
      integer    ,intent(out) :: lyrLo
      integer    ,intent(out) :: lyrHi
      integer    ,intent(out) :: nLayers !Number of layers between observer and target
      integer    ,intent(out) :: nLyr2   !Number of layers with IPATH=2

      integer :: i,j,k,il
      integer :: pathLength, hi2


      !---  Find the lowest and highest layers.
      !     CLBLM layer ordering is from bottom up
      !
      pathLength = size(IPATH) !total number of layers along the path

      lyrLo=0
      do i=1,pathLength
         if ( IPATH(i)==1 .or. IPATH(i)==2 .or. IPATH(i)==3 ) then
            lyrLo = i
            exit
         endif
      enddo

      lyrHi=0
      do i=pathLength,1,-1
         if ( IPATH(i)==1 .or. IPATH(i)==2 .or. IPATH(i)==3 ) then
            lyrHi = i
            exit
         endif
      enddo

      if ( lyrLo == 0 .or.&
           lyrHi == 0 .or.&
           lyrHi < lyrLo ) then
         STOP '--- findHiLoLayers(): Invalid IPATH values.'
      endif

      !--- Number of layers between observer and target.
      nLayers = lyrHi-lyrLo+1


      !--- Find the lowest and highest IPATH=2 layers, if any.
      ! * CLBLM layer ordering is from bottom up
      !   Possible cases are like [22222333] or [22222111], it is
      !   impossible to have a case like [11122222] or [33322222]
      !
      nLyr2 = 0
      if ( any(IPATH(lyrLo:lyrHi)==2) ) then

         hi2 = 0
         do i=lyrHi,lyrLo,-1
            if (IPATH(i)==2) then
               hi2=i
               exit
            endif
         enddo

         !lo2 = lyrLo
         nLyr2 = hi2-lyrLo+1
      endif

   END SUBROUTINE


!-----------------------------------------------------------------------
! Add boundary contributions (thermal emission and thermal reflection)
! to the merged radiances
!-----------------------------------------------------------------------
   SUBROUTINE addThermalBoundary( totRadUp, surf, atmTxA, atmRadUp, atmRadDn, &
                                  fluxCtrl, calcThmRefl, mrgThmRefl, atmTxB)
!-----------------------------------------------------------------------
      USE Module_Scene       ,ONLY: CLBLM_Surface, &
                                    surfThmEmis
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    spectraHaveSameGrid, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid
      USE Module_Config      ,ONLY: CLBLM_FLUX_Ctrl

      type(CLBLM_Spectrum)            ,intent(out)   :: totRadUp
      type(CLBLM_Surface)             ,intent(in)    :: surf
      type(CLBLM_Spectrum)            ,intent(inout) :: atmTxA
      type(CLBLM_Spectrum)  ,optional ,intent(inout) :: atmRadUp
      type(CLBLM_Spectrum)  ,optional ,intent(inout) :: atmRadDn
      logical               ,optional ,intent(in)    :: calcThmRefl  !if =.true. calculate the reflectance of the merged layer including surface reflection.
      type(CLBLM_Spectrum)  ,optional ,intent(out)   :: mrgThmRefl
      type(CLBLM_Spectrum)  ,optional ,intent(inout) :: atmTxB
      type(CLBLM_FLUX_Ctrl)           ,intent(in)    :: fluxCtrl

      character(*) ,parameter :: routineName='addThermalBoundary'

      integer                      :: i,j,k, iv, loc(1)
      logical                      :: mergeRefl, isReflSurf
      type(CLBLM_Spectrum)         :: surfRad, surfEmis, surfRefl
      type(CLBLM_SpectrumPointer)  :: spArray(10)
      logical                      :: flux_flag

      flux_flag = fluxCtrl%flux_flag
      mergeRefl = .FALSE.
      if (present(calcThmRefl)) mergeRefl = calcThmRefl

      isReflSurf = ( allocated(surf%surfEm) .and. any(surf%surfEm <1.))

      !---
      !
      if ( isReflSurf )  then ! surface emission + reflection

         if (.not.present(atmRadDn)) &
            STOP '--- '//routineName//'(): Downwelling radiance must be present for reflecting surface.'

         !--- Find the finest DV and interpolate the others to that grid
         call interpToFinestGrid( atmTxA, atmRadUp, atmRadDn )

         !--- Get surface emission and thermal reflectance
         call surfThmEmis( surf, &
                           atmRadUp%V1,&
                           atmRadUp%DV,&
                           atmRadUp%NLIM, &
                           surfRad, surfEmis )
         call CLBLM_Spectrum_init( surfRefl, &
                                   surfEmis%V1, &
                                   surfEmis%DV, &
                                   surfEmis%NLIM )
         do iv = 1,surfEmis%NLIM
            surfRefl%spect(iv) = 1.0 - surfEmis%spect(iv)
         enddo


         !--- Merge the surface and atmospheric layer
         call CLBLM_Spectrum_init( totRadUp, &
                                   atmRadUp%V1, &
                                   atmRadUp%DV, &
                                   atmRadUp%NLIM )

         IF (flux_flag .eqv. .false.) then
            do iv = 1, totRadUp%NLIM
               totRadUp%spect(iv) = atmRadUp%spect(iv) + &
                                 atmTxA%spect(iv) * &
                                 ( surfRad%spect(iv) + surfRefl%spect(iv)*atmRadDn%spect(iv) )
            enddo
         endif

         IF (flux_flag .eqv. .true.) then
            do iv = 1, totRadUp%NLIM
            totRadUp%spect(iv) = surfRad%spect(iv) + surfRefl%spect(iv)*atmRadDn%spect(iv)
            enddo
         endif

         !--- Combine the surface reflectance and layer transmittances.
         if ( mergeRefl ) then

            if (.not.present(atmTxB)) &
               STOP '--- '//routineName//'(): To merge reflectance (surface+layer), downward Tx must be present.'

            !--- Interpolation.
            call interpToFinestGrid( atmTxA, atmTxB, surfRefl )

            !--- Merge reflection
            call CLBLM_Spectrum_init( mrgThmRefl, &
                                      atmTxA%V1, &
                                      atmTxA%DV, &
                                      atmTxA%NLIM )
            do iv = 1, mrgThmRefl%NLIM
               mrgThmRefl%spect(iv) = atmTxB%spect(iv) * &
                                      atmTxA%spect(iv) * &
                                      surfRefl%spect(iv)
            enddo

         endif ! if ( mergeRefl==.TRUE. ) then


      elseif ( .NOT.isReflSurf .AND. surf%Tskin >0. ) then !--- emission only

         !--- Interpolate to finer grid
         call interpToFinestGrid( atmRadUp, atmTxA  )

         !--- Get surface emission
         call surfThmEmis( surf, &
                           atmRadUp%V1,&
                           atmRadUp%DV,&
                           atmRadUp%NLIM, &
                           surfRad )

         !--- Combine the surface and the atmospheric layer
         call CLBLM_Spectrum_init( totRadUp, &
                                   atmRadUp%V1, &
                                   atmRadUp%DV, &
                                   atmRadUp%NLIM )
         do iv = 1,totRadUp%NLIM
            totRadUp%spect(iv) = atmRadUp%spect(iv) + &
                                 atmTxA%spect(iv)*surfRad%spect(iv)
         enddo

      endif !if ( isReflSurf )

   END SUBROUTINE


!-----------------------------------------------------------------------
! Add surface reflected solar beam radiation to merged upwelling radiances.
! Update combined (surface+layer) solar reflectance, if request.
!-----------------------------------------------------------------------
   SUBROUTINE addSolarBoundary( totRadUp, surf, atmTxA, atmRadUp, solAtmTx, solRadTOA,&
                                calcSolRefl, solMrgRefl )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: PI
      USE Module_Scene       ,ONLY: CLBLM_Surface,&
                                    getSurfRefl
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    spectraHaveSameGrid, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid

      type(CLBLM_Spectrum)           ,intent(out)   :: totRadUp
      type(CLBLM_Surface)            ,intent(in)    :: surf
      type(CLBLM_Spectrum)           ,intent(inout) :: atmTxA
      type(CLBLM_Spectrum)           ,intent(inout) :: atmRadUp
      type(CLBLM_Spectrum)           ,intent(inout) :: solAtmTx
      type(CLBLM_Spectrum)           ,intent(inout) :: solRadTOA  !Solar radiance at TOA. = S0*umu0/pi
      logical              ,optional ,intent(in)    :: calcSolRefl
      type(CLBLM_Spectrum) ,optional ,intent(out)   :: solMrgRefl


      integer              :: i,j,k, iv, loc(1)
      logical              :: mergeRefl
      type(CLBLM_Spectrum) :: surfRefl
      type(CLBLM_SpectrumPointer) :: spArray(10)


      mergeRefl = .FALSE.
      if (present(calcSolRefl)) mergeRefl = calcSolRefl

      !--- Interpolation
      call interpToFinestGrid( atmTxA, atmRadUp, solAtmTx, solRadTOA )

      !--- Read in the solar reflection function
      call getSurfRefl( surf, &
                        atmRadUp%V1,&
                        atmRadUp%DV,&
                        atmRadUp%NLIM, &
                        surfRefl )

      !--- Add reflected solar radiation
      call CLBLM_Spectrum_init( totRadUp, &
                                atmTxA%V1, &
                                atmTxA%DV, &
                                atmTxA%NLIM )

      do iv = 1, totRadUp%NLIM
         totRadUp%spect(iv) = atmRadUp%spect(iv) + &
                              atmTxA%spect(iv) * &
                              surfRefl%spect(iv) * &
                              solAtmTx%spect(iv) * solRadTOA%spect(iv)
      enddo


      !--- Combine the surface reflectance and layer transmittance
      if ( mergeRefl ) then

         call CLBLM_Spectrum_init( solMrgRefl, atmTxA%V1, &
                                               atmTxA%DV, &
                                               atmTxA%NLIM )
         do iv = solMrgRefl%indV1 , solMrgRefl%indV2
            solMrgRefl%spect(iv) = atmTxA%spect(iv) * &
                                   solAtmTx%spect(iv) * &
                                   surfRefl%spect(iv)
         enddo
      endif

   END SUBROUTINE

!-----------------------------------------------------------------------
! Combine the surface reflection and atmospheric transmittances to form a
! merged solar beam reflectance for combined atmos. layer and surface.
!-----------------------------------------------------------------------
   SUBROUTINE mergeSolBoundaryRefl( solMrgRefl, surf, atmTxA, solAtmTx )
!-----------------------------------------------------------------------
      USE Module_Scene       ,ONLY: CLBLM_Surface,&
                                    getSurfRefl
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init, &
                                    spectraHaveSameGrid, &
                                    CLBLM_SpectrumPointer, &
                                    interpToFinestGrid

      type(CLBLM_Spectrum)         ,intent(out)   :: solMrgRefl
      type(CLBLM_Surface)          ,intent(in)    :: surf
      type(CLBLM_Spectrum)         ,intent(inout) :: atmTxA
      type(CLBLM_Spectrum)         ,intent(inout) :: solAtmTx
      type(CLBLM_SpectrumPointer) :: spArray(10)


      integer              :: i,j,k, iv, loc(1)
      type(CLBLM_Spectrum) :: surfRefl


      !--- Interpolation
      call interpToFinestGrid( solAtmTx,atmTxA )

      !--- Read in the solar reflection function
      call getSurfRefl( surf, &
                        atmTxA%V1,&
                        atmTxA%DV,&
                        atmTxA%NLIM, &
                        surfRefl )

      !--- Combine the surface reflectance and layer transmittances
      call CLBLM_Spectrum_init( solMrgRefl, atmTxA%V1, &
                                            atmTxA%DV, &
                                            atmTxA%NLIM )
      do iv = solMrgRefl%indV1 , solMrgRefl%indV2
         solMrgRefl%spect(iv) = atmTxA%spect(iv) * &
                                solAtmTx%spect(iv) * &
                                surfRefl%spect(iv)
      enddo

   END SUBROUTINE


END MODULE
