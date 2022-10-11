
!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_noScattJac

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: noScattJacob


CONTAINS !======================  MODULE CONTAINS ======================


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine noScattJacob( paths, surf, spectGrid, odCtrl, dvCtrl, postCtrl,&
                            SolarOn, solRadTOA, linInTau, ajMolNames, &
                            doJacMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, doRad, &
                            JacMol, JacTemp, JacTskin, JacEmis, JacRsfc, JacRad, doPreBox )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, molIndex
      USE Module_Utility     ,ONLY: getLun
      USE Module_Config      ,ONLY: CLBLM_DV_Ctrl, &
                                    CLBLM_OD_Ctrl, &
                                    CLBLM_Post_Ctrl
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum,&
                                    CLBLM_Spectrum_init, &
                                    interpSpectrum, &
                                    moveSpectrum
      USE Module_AtmPath     ,ONLY: CLBLM_PathCompound, &
                                    CLBLM_Path, &
                                    CLBLM_Layer
      USE Module_Scene       ,ONLY: CLBLM_Scene, &
                                    CLBLM_Surface, &
                                    surfThmEmis, &
                                    getSurfRefl
      IMPLICIT NONE

      type(CLBLM_PathCompound)  ,intent(in)  :: paths
      type(CLBLM_Surface)       ,intent(in)  :: surf
      type(CLBLM_SpectGrid)     ,intent(in)  :: spectGrid     !input spectral grid info.
      type(CLBLM_OD_Ctrl)       ,intent(in)  :: odCtrl        !input ODLAY options
      type(CLBLM_DV_Ctrl)       ,intent(in)  :: dvCtrl
      type(CLBLM_Post_Ctrl)     ,intent(in)  :: postCtrl
      logical                   ,intent(in)  :: SolarOn
      type(CLBLM_Spectrum)      ,intent(in)  :: solRadTOA
      integer                   ,intent(in)  :: linInTau      ! =0 linear-in-tau not used; =2 standard; =1 LBLRTM version
      character(*)              ,intent(in)  :: ajMolNames(:)
      logical                   ,intent(in)  :: doJacMol, doJacTemp, doJacTskin, &
                                                doJacEmis, doJacRsfc, doRad
      type(CLBLM_Spectrum)      ,intent(out) :: JacMol(:,:)  ![nLev,nMol], A.J. w.r.t molecular contentration
      type(CLBLM_Spectrum)      ,intent(out) :: JacTemp(:)   ![nLev],      A.J. w.r.t air temperature
      type(CLBLM_Spectrum)      ,intent(out) :: JacTskin     !             A.J. w.r.t surface temperature
      type(CLBLM_Spectrum)      ,intent(out) :: JacEmis      !             A.J. w.r.t surface emissivity
      type(CLBLM_Spectrum)      ,intent(out) :: JacRsfc      !             A.J. w.r.t surface reflectance
      type(CLBLM_Spectrum)      ,intent(out) :: JacRad       !             total radiance
      logical                   ,intent(out) :: doPreBox     ! output doPreBox flag for deconvolution use.

      !--- Local variables
      !
      character(*) ,parameter :: routinename='noScattJacob'

      integer                           :: i, iv, il, LL, im, nLay, nlv, N1,N2
      type(CLBLM_Path)                  :: vertPath
      type(CLBLM_Path)                  :: viewPath
      type(CLBLM_Path)                  :: downPath
      type(CLBLM_Path)                  :: solPath
      integer                           :: obsLevel
      logical                           :: scalThermalPath !if=.T. down path OD obtained from scaling the up path OD
      logical                           :: scalSolarPath !if=.T. solar path OD obtained from scaling the up path OD
      integer                           :: refPath       !if=0, user vertical as reference path, if=1, use uppath as reference
      logical                           :: UpLook  !if=.T. up looking sensor
      integer              ,allocatable :: ajMolInd(:)
      integer                           :: nAJMol
      integer                           :: blockSize, NLIM, boxcarSize
      real(r8)                          :: V1B, V2B, vn
      real                              :: DVB
      type(CLBLM_SpectGrid)             :: subGrid
      type(CLBLM_DV_Ctrl)               :: subDvCtrl
      type(CLBLM_Spectrum) ,allocatable :: pathOD(:), pathOD_dn(:), pathOD_sun(:)
      type(CLBLM_Spectrum) ,allocatable :: pathNLTE(:), pathNLTE_dn(:)
      real                 ,allocatable :: viewOD(:), downOD(:), solOD(:)
      real                 ,allocatable :: amfView(:), amfDown(:), amfSun(:)
      real                              :: Tsfc, emis, Rsfc_sun, solRad
      type(CLBLM_Spectrum)              :: surfEmis, surfRefl
      type(CLBLM_Spectrum)              :: subSolRadTOA
      real                 ,allocatable :: dRad_dTlev(:,:), dRad_dTlev_v(:)
      real                 ,allocatable :: dRad_dTau(:,:),  dRad_dTau_v(:)
      real                 ,allocatable :: Rad(:)
      real                 ,allocatable :: dRad_dTsfc(:)
      real                 ,allocatable :: dRad_dEmis(:)
      real                 ,allocatable :: dRad_dRsfc(:)
      real                              :: rad_v
      real                              :: dRad_dTsfc_v
      real                              :: dRad_dEmis_v
      real                              :: dRad_dRsfc_v
      real                 ,allocatable :: dT_dTup(:)
      real                 ,allocatable :: dT_dTlo(:)
      real                 ,allocatable :: dW_dQup(:,:)
      real                 ,allocatable :: dW_dQlo(:,:)
      real                 ,allocatable :: QLev(:,:)


      !---
      scalThermalPath = paths%AMscalingThermal !if=.T. uppath and downpath OD obtained from scaling the reference path OD
      scalSolarPath   = paths%AMscalingSolar   !if=.T. solar path OD obtained from scaling the reference path OD
      refPath         = paths%refPath          !if=0, user vertical as reference path, if=1, use uppath as reference
      vertPath = paths%vert
      viewPath = paths%view
      downPath = paths%down
      solPath  = paths%sun
      obsLevel = paths%obsLev

      !--- Get indexes of A.J. molecules in path%molID array.
      if (doJacMol) then
         nAJMol = size(ajMolNames)
         allocate( ajMolInd( nAJMol))
         do im = 1,nAJMol
            ajMolInd(im) = molIndex( ajMolNames(im), viewPath%molID ) !location of a.j. molecule in active molecular array
         enddo
         if (any(ajMolInd<0)) STOP '--- '//routineName//'(): Invalid molecular name for Jacobian calculation.'
      else
         nAJMol = 0
      endif

      UpLook = viewPath%geom%ANGLE <=90.
      nLay = viewPath%nLay
      blockSize = 100000!30000


      !--- Allocate memory for block arrays
      !
      allocate( pathOD(      nLay))
      allocate( pathOD_dn(   nLay))
      allocate( pathOD_sun(  nLay))
      allocate( pathNLTE(    nLay))
      allocate( pathNLTE_dn( nLay))
      allocate( viewOD(  nLay))
      allocate( downOD(  nLay))
      allocate( solOD(   nLay))
      allocate( amfView( nLay))
      allocate( amfDown( nLay))
      allocate( amfSun(  nLay))
      allocate( dRad_dTlev( blockSize, nLay+1), dRad_dTlev_v(nLay+1))
      allocate( dRad_dTau(  blockSize, nLay  ), dRad_dTau_v( nLay  ))
      allocate( Rad( blockSize))
      allocate( dRad_dTsfc( blockSize))
      allocate( dRad_dEmis( blockSize))
      allocate( dRad_dRsfc( blockSize))
      allocate( dT_dTup( nLay))
      allocate( dT_dTlo( nLay))
      allocate( dW_dQup( nLay, nAJMol))
      allocate( dW_dQlo( nLay, nAJMol))
      allocate( QLev(  nLay+1, nAJMol))

      !--- Calculate dT_dTup, dT_dTlo, dW_dQup, dW_dQlo
      do im =1,nAJMol
         QLev( 1:nLay+1,im ) = viewPath%Q( ajMolInd(im), 1:nLay+1 )
      enddo
!viewPath or vertPath?
      CALL level2layer( viewPath%pRT, viewPath%T, viewPath%Tave, QLev, &
                        nAJMol, viewPath%nLay, &
                        dT_dTup, dT_dTlo, dW_dQup, dW_dQlo )


      !--- Load atmospheric and surface parameters for input to layerDerivative()
      !
      N1 = 1
      N2 = nLay
                                   amfView(:) = viewPath%Wtot(:)/vertPath%Wtot(:)
      if (.not.UpLook)             amfDown(:) = downPath%Wtot(:)/vertPath%Wtot(:)
      if (.not.UpLook.and.SolarOn) amfSun( :) = solPath%Wtot( :)/vertPath%Wtot(:)


      !---
      if (spectGrid%DVout<0) STOP '--- '//routineName//'(): Spectral grid is not uniform.'
      DVB = spectGrid%DVout
      V1B = spectGrid%V1
      V2B = V1B - DVB

      doPreBox =  postCtrl%boxcarHW >0. .and. postCtrl%boxcarHW >DVB
      boxcarSize = int( 2.0 * postCtrl%boxcarHW / DVB )
      !if (mod(boxcarSize,2)/=0) boxcarSize=boxcarSize-1  ! boxcarSize to even number
      if (boxcarSize <=1) doPreBox = .FALSE.             ! If boxcar width <= DVB, don't do prebox smoothing


      !--- Loop over sub-intervals to calculate analytic Jacobians
      !
      DO WHILE (V2B < spectGrid%V2)

         V1B = V2B + DVB
         V2B = V1B + (blockSize-2)*DVB !ODLAY will set the size to be ceiling((V2B-V1B)/DVB)+1.
         if (V2B > spectGrid%V2) then
            V2B = spectGrid%V2
         endif

         subGrid = spectGrid
         subGrid%V1 = V1B
         subGrid%V2 = V2B

         subDvCtrl = dvCtrl
         !subDvCtrl%gridType = 1 !uniform DV interp to finest DV for the path


         !--- Calculate path ODs
         ! * for      UpLook: pathOD and pathNLTE contains layer ODs from obsLevel to TOA only
         !   for .not.UpLook: pathOD and pathNLTE contains layer ODs from obsLevel-1 to surface only
         ! * pathOD_dn and pathNLTE_dn only available for .not.UpLook caes
         ! * pathOD_sun only available for .not.UpLook and SolarOn case
         CALL pathOD_up_dn_sun ( vertPath, viewPath, downPath, solPath, obsLevel,UpLook, &
                                 SolarOn, scalThermalPath, scalSolarPath, refPath, &
                                 odCtrl, subGrid, subDvCtrl, &
                                 pathOD, pathNLTE, &
                                 pathOD_dn, pathNLTE_dn, &
                                 pathOD_sun)


         !--- Check if the returned V2 differ from initial V2
         if (UpLook) then
            LL=nLay
         else
            LL=1
         endif
         if ( V2B /= pathOD(LL)%V2 ) V2B = pathOD(LL)%V2
         NLIM = pathOD(LL)%NLIM

         !--- Recalculate the size if prebox averaging is requested
         if ( doPreBox ) then
            if (NLIM>=boxcarSize) then
               NLIM = NLIM - mod(NLIM,boxcarSize)
               V2B = V1B + (blockSize-1)*DVB
            else
               NLIM = mod(NLIM,boxcarSize)
               V2B = V1B + (NLIM-1)*DVB
            endif
         endif



         !--- If NLIM differ from initial value, resize the arrays
         if (NLIM /= blockSize) then
            deallocate(dRad_dTlev) ;allocate( dRad_dTlev( NLIM, nLay+1))
            deallocate(dRad_dTau)  ;allocate( dRad_dTau(  NLIM, nLay))
            deallocate(Rad)        ;allocate( Rad(        NLIM))
            deallocate(dRad_dTsfc) ;allocate( dRad_dTsfc( NLIM))
            deallocate(dRad_dEmis) ;allocate( dRad_dEmis( NLIM))
            deallocate(dRad_dRsfc) ;allocate( dRad_dRsfc( NLIM))
         endif


         !--- Interpolate surface emissivity, reflectance and solar irradiance to subGrid
         !
         if (.not.UpLook)               call surfThmEmis( surf, V1B,DVB,NLIM, sfcEmis=surfEmis )
         if (.not.UpLook .and. SolarOn) call getSurfRefl( surf, V1B,DVB,NLIM, surfRefl )

         if (SolarOn) then
            call CLBLM_Spectrum_init( subSolRadTOA, V1B,DVB,NLIM )
            call interpSpectrum( solRadTOA, subSolRadTOA )
         endif


         !--- Call layer derivative subroutine
         !
         DO iv = 1,NLIM

            vn = V1B + (iv - 1) * DVB
            if (UpLook) then
               viewOD( obsLevel:nLay ) = [(pathOD(i)%spect(iv), i=obsLevel,nLay)] !viewOD( N1:obsLevOSS-1) = [(pathOD(i)%spect(iv), i=nLay,obsLevel,-1)]
            else
               viewOD( 1:obsLevel-1 ) = [(pathOD(i)%spect(iv), i=1,obsLevel-1)]  !viewOD( obsLevOSS:N2) = [(pathOD(i)%spect(iv), i=obsLevel-1,1,-1)] !pathOD(obsLevel-1:1:-1)%spect(iv)
               downOD( 1:nLay ) = [(pathOD_dn(i)%spect(iv), i=1,nLay)]  !downOD( N1:N2)        = [(pathOD_dn(i)%spect(iv), i=nLay,1,-1)] !pathOD_dn(   nLay:1:-1)%spect(iv)
               if (SolarOn) then
                  solOD( 1:nLay ) = [(pathOD_sun(i)%spect(iv), i=1,nLay)]  !solOD(N1:N2) = [(pathOD_sun(i)%spect(iv), i=nLay,1,-1)] !pathOD_sun(nLay:1:-1)%spect(iv)
                  Rsfc_sun = surfRefl%spect(iv)
                  solRad = subSolRadTOA%spect(iv)
               endif
               emis = surfEmis%spect(iv)
            endif


            CALL Jacobians( viewOD, downOD, solOD, amfView, amfDown, amfSun, &                         !input
                            viewPath%Tave, viewPath%T, surf%Tskin, emis, SolarOn, solRad, Rsfc_sun, &  !input
                            vn, N1, N2, obsLevel, UpLook, linInTau, &                                  !input
                            nAJMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, doRad, &              !input
                            dT_dTup, dT_dTlo, &                                                        !input
                            dRad_dTlev_v, dRad_dTau_v, &                                               !out
                            rad_v, dRad_dTsfc_v, dRad_dEmis_v, dRad_dRsfc_v )                          !out

            if (doJacTemp .or. nAJMol>0) dRad_dTau( iv, N1:N2  ) = dRad_dTau_v( N1:N2)
            if (doJacTemp)               dRad_dTlev(iv, N1:N2+1) = dRad_dTlev_v(N1:N2+1)
            if (doRad)                   Rad(       iv         ) = rad_v
            if (doJacTskin)              dRad_dTsfc(iv         ) = dRad_dTsfc_v
            if (doJacEmis)               dRad_dEmis(iv         ) = dRad_dEmis_v
            if (doJacRsfc)               dRad_dRsfc(iv         ) = dRad_dRsfc_v

         ENDDO !do iv = 1,numFreq, loop over frequency points

         !--- Load the output arrays
         call output_SurfJac_Rad( V1B,V2B,DVB,NLIM, doPreBox, boxcarSize, postCtrl%boxcarHW, & !this line used for output results
                                  dRad_dTsfc, dRad_dEmis, dRad_dRsfc, Rad, &
                                  doJacTskin, doJacEmis,  doJacRsfc,  doRad, &
                                  JacTskin,   JacEmis,    JacRsfc,    JacRad )


         !--- Calculate and write derivatives w.r.t level temperature
!viewPath or vertPath?
         if (doJacTemp) then
            CALL drv_LevelTemp( viewPath, pathOD, obsLevel, UpLook, &
                                odCtrl, subGrid, subDvCtrl, &
                                dRad_dTlev, dRad_dTau, dT_dTup, dT_dTlo, &
                                linInTau, &
                                V1B,V2B,DVB,NLIM, doPreBox, boxcarSize, postCtrl%boxcarHW, & !this line used for output results
                                JacTemp )
         endif

         !--- Calculate and write derivatives w.r.t molecular concentration
         if (nAJMol >0) then
!(1)viewPath or vertPath?  (2)if vertpath, uplooking case need from obslev to TOA only, be careful the lyrNo
            CALL drv_LevelConc( ajMolInd, viewPath, obsLevel, UpLook, &
                                odCtrl, subGrid, subDvCtrl, &
                                dRad_dTau, dW_dQup, dW_dQlo, &
                                V1B,V2B,DVB,NLIM, doPreBox, boxcarSize, postCtrl%boxcarHW, & !this line used for output results
                                JacMol )
         endif


!         !--- Update the number of points done.
!         numPts = numPts + NLIMwrite

      ENDDO !do while (V2B <= spectGrid%V2)


   END SUBROUTINE


!----------------------------------------------------------------------------
! PURPOSE: Compute radiances (in mw/m2/str/cm-1) and derivatives of radiances
!          with respect to layer OD, layer temperature and surface parameters
!----------------------------------------------------------------------------
   SUBROUTINE Jacobians( viewOD, downOD, solOD, amfView, amfDown, amfSun, &
                         Tlay, Tlev, Tsfc, emis, SolarOn, solRadTOA, Rsfc_sun, &
                         vn, N1, N2, obsLev, UpLook, linInTau, &
                         nAJMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, doRad, &
                         dT_dTup, dT_dTlo, &
                         dRad_dTlev, dRad_dTau, &
                         rad, dRad_dTsfc, dRad_dEmis, dRad_dRsfc )
!----------------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_EMLAY       ,ONLY: planckFnDrv, rlin
      IMPLICIT NONE

      !---Input variables
      real      ,intent(in)  :: viewOD(:)       ![:],       OD along upwelling path
      real      ,intent(in)  :: downOD(:)       ![N1:N2],   OD along downwelling path
      real      ,intent(in)  :: solOD(:)        ![N1:N2],   OD along solar path
      real      ,intent(in)  :: amfView(:)      ![N1:N2],   air mass factors
      real      ,intent(in)  :: amfDown(:)      ![N1:N2],   air mass factors
      real      ,intent(in)  :: amfSun(:)       ![N1:N2],   air mass factors
      real      ,intent(in)  :: Tlay(:)         ![N1:N2],   layer averade temperature
      real      ,intent(in)  :: Tlev(:)         ![N1:N2+1], level temperature
      real      ,intent(in)  :: Tsfc            !           surface temperature
      real      ,intent(in)  :: emis            !           surface emissivity
      logical   ,intent(in)  :: SolarOn         !           if true SolarOn contribution is computed
      real      ,intent(in)  :: solRadTOA       !           solar radiance at TOA
      real      ,intent(in)  :: Rsfc_sun        !           surface reflectivity
      real(r8)  ,intent(in)  :: vn              !           wavenumber
      integer   ,intent(in)  :: N1,N2           !           first layer to include, last layer to include
      integer   ,intent(in)  :: obsLev
      logical   ,intent(in)  :: UpLook
      integer   ,intent(in)  :: linInTau        !           if ture, use linear-in-tau approximation
      integer   ,intent(in)  :: nAJMol
      logical   ,intent(in)  :: doJacTemp, doJacTskin, doJacEmis, doJacRsfc, doRad
      real      ,intent(in)  :: dT_dTup(:)      ![N1:N2]
      real      ,intent(in)  :: dT_dTlo(:)      ![N1:N2]
      real      ,intent(out) :: dRad_dTlev(:)   ![N1:N2+1]
      real      ,intent(out) :: dRad_dTau(:)    ![N1:N2]
      real      ,intent(out) :: rad             !           computed upwelling radiance
      real      ,intent(out) :: dRad_dTsfc      !           radiance derivatives with respect to surface temperature
      real      ,intent(out) :: dRad_dEmis      !           radiance derivatives with respect to surface emissivity
      real      ,intent(out) :: dRad_dRsfc      !           radiance derivatives with respect to surface solar reflectance

      !---Local variables:
      integer ,parameter :: MxLev = 200
      INTEGER   :: l,lv,L1,L2,obsLay
      REAL      :: rsfc, radsun
      real      :: cumODA, cumODB, cumTxA(MxLev), cumTxB(MxLev), cumTxB2(MxLev)
      REAL      :: totODUp, totTxUp, totODDn, totODSun, totTxSun2
      REAL      :: dbavgdb, dbavgdbdod, Blt
      real      :: dTran
      real      :: dTsfc_dTup
      real      :: dTsfc_dTlo
      real      :: Bsfc, dBsfc, sfcEmit
      !REAL      :: Bbar(MxLev), dBbar(MxLev)
      !REAL      :: Blev(MxLev), dBlev(MxLev)
      REAL      :: Bbar(N1:N2), dBbar(N1:N2)
      REAL      :: Blev(N1:N2+1), dBlev(N1:N2+1)
      real      :: dRad_dT_up(MxLev)
      real      :: dRad_dT_dn(MxLev)
      real      :: dRad_dTbar(MxLev)

      !--- Compute Planck function and its derivative wrt temperature
      ! * for LinInTau Blev and dBlev should be determined
      !
      if (linInTau/=0) then
         call planckFnDrv(vn, Tlev, Blev, dBlev)
         dTsfc_dTup=0.
         dTsfc_dTlo=1.
      endif
      call planckFnDrv(vn, Tlay, Bbar, dBbar )
      call planckFnDrv(vn, Tsfc, Bsfc, dBsfc)

      !--- Initialize radiance and derivative arrays
      rad                  = 0.
      radsun               = 0.
      dRad_dTau(N1:N2)     = 0.
      dRad_dTbar(N1:N2)    = 0.
      dRad_dT_up(N1:N2+1)  = 0.
      dRad_dT_dn(N1:N2+1)  = 0.
      dRad_dTLev(N1:N2+1)  = 0.
      dRad_dEmis           = 0.
      dRad_dTsfc           = 0.
      dRad_dRsfc           = 0.

      obsLay = obsLev


      IF ( UpLook ) THEN !Up looking sensor

         !        |
         !      | |
         !    | | | ...
         ! ------------- obsLev
         !
         cumTxB(N1:obsLev) = 1.0 !cumTxB(obsLev:N2+1)  = 1.0
         cumODB            = 0.0
         obsLay            = max(obsLay,N1) !min(obsLev-1,N2)
         DO l=obsLay,N2 !N2prim,1,-1
            cumODB      = cumODB + viewOD(l) !tauTot(l)*sec(l)
            cumTxB(l+1) = EXP(-cumODB)
         END DO

         if ( linInTau/=0 ) then

            DO l=N2,obsLay,-1  !1,N2prim
               dTran = cumTxB(l)-cumTxB(l+1) !cumTxB(l+1)-cumTxB(l)
               !odsec = tauTot(l)*sec(l)

               call rlin( viewOD(l), dbavgdb, dbavgdbdod )

               if ( linInTau == 2 ) then !Standard linear-in-tau

                  Blt             = Blev(l)*(1-dbavgdb) + Blev(l+1)*dbavgdb           !Blev(l+1)*(1-dbavgdb) + Blev(l)*dbavgdb
                  dRad_dTau(l)    = (cumTxB(l+1)*Blt - rad + dTran * dbavgdbdod * &   !(cumTxB(l)*Blt - rad + dTran * dbavgdbdod * &
                                    (Blev(l+1)-Blev(l)))* amfView(l) !sec(l)          !(Blev(l)-Blev(l+1)))* amfView(l) !sec(l)
                  if (doJacTemp) then
                     dRad_dT_dn(l+1) = dRad_dT_dn(l+1) + dTran*dBlev(l+1)*dbavgdb     !dRad_dT_dn(l)   = dRad_dT_dn(l) + dTran*dBlev(l)*dbavgdb
                     dRad_dT_dn(l)   =                   dTran*dBlev(l)*(1.0-dbavgdb) !dRad_dT_dn(l+1) =                 dTran*dBlev(l+1)*(1.0-dbavgdb)
                     dRad_dTbar(l)   = 0.0                                            !dRad_dTbar(l)   = 0.0
                  endif
               elseif ( linInTau == 1 ) then !LBLRTM linear-in-tau

                  Blt             = 2.0*Bbar(l)*dbavgdb + Blev(l)*(1.0-2.0*dbavgdb)    !2.0*Bbar(l)*dbavgdb + Blev(l+1)*(1.0-2.0*dbavgdb)
                  dRad_dTau(l)    = (cumTxB(l+1)*Blt - rad + 2.0*dTran*dbavgdbdod * &  !(cumTxB(l)*Blt - rad + 2.0*dTran*dbavgdbdod * &
                                    (Bbar(l)-Blev(l))) * amfView(l) !sec(l)          !(Bbar(l)-Blev(l+1))) * amfView(l) !sec(l)
                  if (doJacTemp) then
                     dRad_dT_dn(l) = dTran*dBlev(l)*(1.0-2.0*dbavgdb)  !dRad_dT_dn(l+1) = dTran*dBlev(l+1)*(1.0-2.0*dbavgdb)
                     dRad_dTbar(l) = dTran*2.0*dBbar(l)*dbavgdb        !dRad_dTbar(l)   = dTran*2.0*dBbar(l)*dbavgdb
                  endif
               end if

               rad = rad + dTran*Blt
            END DO

         else !Not linear-in-tau

            DO l=N2,obsLay,-1 !1,N2prim
               dTran         = cumTxB(l)-cumTxB(l+1)                            !cumTxB(l+1)-cumTxB(l)
               dRad_dTau(l)  = (cumTxB(l+1)*Bbar(l) -rad) * amfView(l) !sec(l)  !(cumTxB(l)*Bbar(l) -rad) * amfView(l) !sec(l)
               rad           = dTran*Bbar(l) + rad

               if (doJacTemp) dRad_dTbar(l) = dTran*dBbar(l)
            ENDDO
         endif

         !DO l=1,N2prim
         !   dRad_dTbar(l) = dRad_dTbar(l) + dRad_dTau(l)*dTau_dT !dTaudTmp(l)
         !end do

      ELSE  !Down looking sensor

         ! 0- Transmittance profile along upwelling path from TOA down to surface
         ! ------------- TOA
         !    | | | ...
         !      | |
         !        |
         !
         cumODA  = 0.                                      !cumulative OD, cumulated down ward from TOA to surface.
         cumTxA(obsLev:N2+1) = 1.  !cumTxA(1:obsLev) = 1.  !cumulative Tx, cumulated down ward from TOA to surface.
         DO l=obsLay-1,N1,-1 !obsLev,N2
            cumODA = cumODA + viewOD(l) !lyrOD(l)*sec
            cumTxA(l) = EXP(-cumODA)    !cumTxA(l+1)  = EXP(-cumODA)
         END DO
         totODUp = cumODA
         totTxUp = cumTxA(N1) !cumTxA(N2+1)


         ! 1- Downwelling thermal radiance calculation:
         !------------------------------------------------------------------
         !
         IF(totTxUp.GT.1.e-06)THEN

            !        |
            !      | |
            !    | | | ...
            ! ------------- Surface
            !
            cumTxB2(N1) = totTxUp  !cumTxB2(N2+1) = totTxUp !two-path transmittance, equals totTxUp times cumulative Tx, cumulateed from surface upward.
            cumODB = 0.                                     !cumulative OD, cumulated upward from surface to level 1.
            DO l=N1,N2 !N2,N1,-1
               cumODB = cumODB + downOD(l)
               cumTxB2(l+1) = EXP( -(totODUp+cumODB) ) !cumTxB2(l) = EXP( -(totODUp+cumODB) ) !EXP( -(totODUp + cumODB * secdif) )
            END DO

            if (linInTau/=0) then

               DO l=N2,N1,-1 !N1,N2

                  dTran = cumTxB2(l) - cumTxB2(l+1) !cumTxB2(l+1) - cumTxB2(l)
                  !odsec = secdif * lyrOD(l)

                  call rlin( downOD(l), dbavgdb, dbavgdbdod )

                  if ( linInTau==2 ) then !Standard linear-in-tau

                     Blt            = Blev(l) + (Blev(l+1) - Blev(l)) * dbavgdb           !Blev(l+1) + (Blev(l) - Blev(l+1)) * dbavgdb
                     dRad_dTau(l)   = (cumTxB2(l+1) * Blt - rad + dTran * dbavgdbdod * &  !(cumTxB2(l) * Blt - rad + dTran * dbavgdbdod * &
                                      (Blev(l+1)-Blev(l)) ) * amfDown(l)                  !(Blev(l)-Blev(l+1)) ) * amfDown(l)
                     if (doJacTemp) then
                        dRad_dT_dn(l+1) = dRad_dT_dn(l+1) + dTran * dBlev(l+1) * dbavgdb    !dRad_dT_dn(l)   = dRad_dT_dn(l) + dTran * dBlev(l) * dbavgdb
                        dRad_dT_dn(l)   =                   dTran * dBlev(l) * (1.-dbavgdb) !dRad_dT_dn(l+1) =                 dTran * dBlev(l+1) * (1.-dbavgdb)
                        dRad_dTbar(l)   = 0.0
                     endif
                  elseif ( lininTau==1 ) then !LBLRTM linear-in-tau

                     Blt             = 2.0*Bbar(l)*dbavgdb + blev(l)*(1.0-2.0*dbavgdb)  !2.0*Bbar(l)*dbavgdb + blev(l+1)*(1.0-2.0*dbavgdb)
                     dRad_dTau(l)    = (cumTxB2(l+1)*Blt - rad + 2.0*dTran*dbavgdbdod * &   !(cumTxB2(l)*Blt - rad + 2.0*dTran*dbavgdbdod * &
                                       (Bbar(l)-Blev(l))) * amfDown(l) !secRefl(l)      !(Bbar(l)-Blev(l+1))) * amfDown(l) !secRefl(l)
                     if (doJacTemp) then
                        dRad_dT_dn(l) = dTran*dblev(l)*(1.0-2.0*dbavgdb)  !dRad_dT_dn(l+1) = dTran*dblev(l+1)*(1.0-2.0*dbavgdb)
                        dRad_dTbar(l) = dTran*2.0*dBbar(l)*dbavgdb
                     endif
                  endif

                  rad = rad + dTran * Blt
               ENDDO

               !dRad_dT_dn(N2)   = dRad_dT_dn(N2) + dRad_dT_dn(N2+1) * dTsfc_dTup
               !dRad_dT_dn(N2+1) = dRad_dT_dn(N2+1) * dTsfc_dTlo

            else

               DO l=N2,N1,-1 !N1,N2
                  dTran        = cumTxB2(l)-cumTxB2(l+1)                     !cumTxB2(l+1)-cumTxB2(l)
                  dRad_dTau(l) = (cumTxB2(l+1) * Bbar(l) - rad) * amfDown(l) !(cumTxB2(l) * Bbar(l) - rad) * amfDown(l)
                  rad          = rad + dTran * Bbar(l)

                  if (doJacTemp) dRad_dTbar(l) = dTran * dBbar(l)
               ENDDO

            endif !if (linInTau/=0) then

         END IF !IF(totTxUp.GT.1.e-06)THEN


         ! 2- dRad_dTau * (1-emis)
         !------------------------------------------------------------------
         !
         rsfc = (1.-emis)
         DO l=N1,N2
                           dRad_dTau(l)  = dRad_dTau(l)*rsfc
            if (doJacTemp) dRad_dTbar(l) = dRad_dTbar(l)*rsfc
         END DO


         ! 3- Add Surface terms and Derivatives wrt emissivity and sfc skin temperature:
         !------------------------------------------------------------------
         !
         if (doJacEmis)  dRad_dEmis = totTxUp * Bsfc - rad
         if (doJacTskin) dRad_dTsfc = emis * totTxUp * dBsfc
         sfcEmit = emis * totTxUp * Bsfc
         rad     = rad*rsfc + sfcEmit


         ! 4- Add solar component
         !------------------------------------------------------------------
         !
         IF(SolarOn)THEN

            totODSun = 0.   !total OD along solar path
            do l=N1,N2
               totODSun = totODSun + solOD(l)
            enddo
            totTxSun2 = EXP( -(totODUp+totODSun) ) !EXP(-(totODUp + cumODB*sec0))

            dRad_dRsfc = totTxSun2 * solRadTOA
            radsun     = rsfc_sun * dRad_dRsfc
            do l=N1,N2
               dRad_dTau(l) = dRad_dTau(l) - radsun * amfSun(l) !-radsun * (sec+sec0); amfView term will be handled below through rad in upwelling part.
            enddo

            rad = rad + radsun
         END IF


         ! 5- Upwelling thermal radiance calculation
         !------------------------------------------------------------------
         !
         if (linInTau/=0) then

            DO l=N1,obsLay-1 !N2,obsLev,-1

               dTran = cumTxA(l+1) - cumTxA(l) !cumTxA(l) - cumTxA(l+1)
               !odsec = sec * lyrOD(l)

               call rlin( viewOD(l), dbavgdb, dbavgdbdod )

               if ( linInTau==2 ) then !Standard linear-in-tau

                  Blt             = Blev(l+1) * (1-dbavgdb) + Blev(l) * dbavgdb              !Blev(l) * (1-dbavgdb) + Blev(l+1) * dbavgdb
                  dRad_dTau(l)    = dRad_dTau(l) + ( cumTxA(l)*Blt - rad + &                 !dRad_dTau(l) + ( cumTxA(l+1)*Blt - rad + &
                                    dTran * dbavgdbdod * (Blev(l)-Blev(l+1)) ) * amfView(l)  !dTran * dbavgdbdod * (Blev(l+1)-Blev(l)) ) * amfView(l)  !dTran * dbavgdbdod * (Blev(l+1)-Blev(l)) ) * sec + dRad_dTau_sun

                  if (doJacTemp) then
                     dRad_dT_up(l+1) =                 dTran * dBlev(l+1) * (1-dbavgdb)  !dRad_dT_up(l)   =                   dTran * dBlev(l) * (1-dbavgdb)
                     dRad_dT_up(l)   = dRad_dT_up(l) + dTran * dBlev(l) * dbavgdb        !dRad_dT_up(l+1) = dRad_dT_up(l+1) + dTran * dBlev(l+1) * dbavgdb
                  endif
               elseif ( linInTau==1 ) then !LBLRTM linear-in-tau

                  Blt           = 2.0*Bbar(l)*dbavgdb + Blev(l+1)*(1.0-2.0*dbavgdb)            !2.0*Bbar(l)*dbavgdb + Blev(l)*(1.0-2.0*dbavgdb)
                  dRad_dTau(l)  = dRad_dTau(l) + ( cumTxA(l)*Blt - rad + &                     !dRad_dTau(l) + ( cumTxA(l+1)*Blt - rad + &
                                  dTran*dbavgdbdod*2.0 * (Bbar(l) - Blev(l+1)) ) * amfView(l)  !dTran*dbavgdbdod*2.0 * (Bbar(l) - Blev(l)) ) * amfView(l) !sec(l)

                  if (doJacTemp) then
                     dRad_dT_up(l+1) = dTran*dBlev(l+1)*(1.0-2.0*dbavgdb) !dRad_dT_up(l) = dTran*dBlev(l)*(1.0-2.0*dbavgdb)
                     dRad_dTbar(l)   = dRad_dTbar(l) + dTran*2.0*dBbar(l)*dbavgdb
                  endif
               endif

               rad = rad + dTran * Blt
            ENDDO

            !dRad_dT_up(N2)   = dRad_dT_up(N2) + dRad_dT_up(N2+1) * dTsfc_dTup
            !dRad_dT_up(N2+1) = dRad_dT_up(N2+1) * dTsfc_dTlo

         else

            DO l=N1,obsLay-1 !N2,obsLev,-1
               dTran         = cumTxA(l+1) - cumTxA(l) !cumTxA(l) - cumTxA(l+1)
               dRad_dTau(l)  = dRad_dTau(l) + (cumTxA(l) * Bbar(l)-rad) * amfView(l) !dRad_dTau(l) + (cumTxA(l+1) * Bbar(l)-rad) * amfView(l) !dRad_dTau(l)*rsfc + (cumTxA(l+1) * Bbar(l)-rad) * sec + dRad_dTau_sun
               rad           = rad + dTran*Bbar(l)
               !dRad_dT(l)    = dRad_dT(l)*rsfc + dTran*dBbar(l) + dRad_dTau(l)*dTau_dT(l)

               if (doJacTemp) dRad_dTbar(l) = dRad_dTbar(l) + dTran*dBbar(l)
            ENDDO

         endif !if (linInTau/=0) then

         !DO l=1,numLay
         !   dRad_dTbar(l) = dRad_dTbar(l) + dRad_dTau(l)*dTau_dT !dTaudTmp(l)
         !END DO

      ENDIF !IF ( UpLook ) THEN


      ! dRad_dEmission * dEmission*dTlev
      !-----------------------------------------------------------------
      if (doJacTemp) then

         if (linInTau/=0) then

            if (UpLook) then
               do lv = obsLev,N2+1 !N1,obsLev
                  dRad_dTlev(lv) = dRad_dT_dn(lv)
               enddo
            else
               do lv = N1,N2+1 !N1,N2+1  !obsLev,N2+1
                  dRad_dTlev(lv) = dRad_dT_up(lv) + dRad_dT_dn(lv) * rsfc
               enddo
            endif
         endif


         if (linInTau==0 .or. linInTau==1) then

            if (UpLook) then
               L1 = N2      !L1 = N1
               L2 = obsLay  !L2 = obsLev-1
            else
               L1 = N2      !L1 = N1  !obsLev
               L2 = N1      !L2 = N2
            endif

            do l = L1,L2,-1 !L1,L2
               dRad_dTlev(l+1) = dRad_dTlev(l+1) + dRad_dTbar(l)*dT_dTup(l)   !dRad_dTlev(l)   = dRad_dTlev(l)   + dRad_dTbar(l)*dT_dTup(l)
               dRad_dTlev(l)   = dRad_dTlev(l)   + dRad_dTbar(l)*dT_dTlo(l)   !dRad_dTlev(l+1) = dRad_dTlev(l+1) + dRad_dTbar(l)*dT_dTlo(l)
            enddo
         endif
      endif

   END SUBROUTINE



!----------------------------------------------------------------------------
! Calculate optical depth for up, down and solar path in a sub-interval
!----------------------------------------------------------------------------
   SUBROUTINE pathOD_up_dn_sun ( &
                  vertPath, viewPath, downPath, solPath, obsLevel, UpLook, &
                  SolarOn, scalThermalPath, scalSolarPath, refPath, &
                  odCtrl, subGrid, subDvCtrl, &
                  pathOD, pathNLTE, pathOD_dn, pathNLTE_dn, pathOD_sun)
!----------------------------------------------------------------------------
      USE Module_Config      ,ONLY: CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init
      USE Module_AtmPath     ,ONLY: CLBLM_Path, &
                                    CLBLM_Layer, &
                                    loadLayerFromPath
      USE Module_ODLAY       ,ONLY: ODLAY
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      IMPLICIT NONE

      type(CLBLM_Path)          ,intent(in)  :: vertPath
      type(CLBLM_Path)          ,intent(in)  :: viewPath
      type(CLBLM_Path)          ,intent(in)  :: downPath
      type(CLBLM_Path)          ,intent(in)  :: solPath
      integer                   ,intent(in)  :: obsLevel
      logical                   ,intent(in)  :: UpLook
      logical                   ,intent(in)  :: SolarOn
      logical                   ,intent(in)  :: scalThermalPath, scalSolarPath
      integer                   ,intent(in)  :: refPath      !if=0, ref. to vertical, if=1, ref. to up path
      type(CLBLM_OD_Ctrl)       ,intent(in)  :: odCtrl       !input ODLAY options
      type(CLBLM_SpectGrid)     ,intent(in)  :: subGrid      !input spectral grid info.
      type(CLBLM_DV_Ctrl)       ,intent(in)  :: subDvCtrl
      type(CLBLM_Spectrum)      ,intent(out) :: pathOD(:)
      type(CLBLM_Spectrum)      ,intent(out) :: pathOD_dn(:)
      type(CLBLM_Spectrum)      ,intent(out) :: pathOD_sun(:)
      type(CLBLM_Spectrum)      ,intent(out) :: pathNLTE(:)
      type(CLBLM_Spectrum)      ,intent(out) :: pathNLTE_dn(:)

      !--- Local variables
      integer                           :: iv, il, L1,L2, nLay
      type(CLBLM_Layer)                 :: aLayer
      type(CLBLM_Spectrum)              :: dummySpect
      type(CLBLM_Spectrum) ,allocatable :: refPathOD(:), refPathNLTEFac(:)
      real                 ,allocatable :: scalDownFac(:), &
                                           scalUpFac(:), &
                                           scalSolarFac(:)



      !--- If scaling approximation used, calculate OD along the reference path
      if ((scalThermalPath .or. scalSolarPath) .and. .not.upLook) then
         L1 = 1
         if (refPath==1) then !use uppath as reference
            L2 = viewPath%nLay
         else
            L2 = vertPath%nLay
         endif

         if (allocated(refPathOD)) deallocate(refPathOD); allocate(refPathOD(L1:L2))
         !if (NLTE)then
         !   if allocated(refPathNLTEFac) deallocate(refPathNLTEFac);  allocate(refPathNLTEFac(L1:L2))
         !endif

         if (allocated(scalDownFac))  deallocate(scalDownFac);  allocate(scalDownFac( L1:L2))
         if (allocated(scalUpFac))    deallocate(scalUpFac);    allocate(scalUpFac(   L1:L2))
         if (allocated(scalSolarFac)) deallocate(scalSolarFac); allocate(scalSolarFac(L1:L2))


         do il = L1, L2
            !--- Load layer temperature, pressure and concentration ...
            if     (refPath==1) then; call loadLayerFromPath( aLayer, viewPath, il )
            elseif (refPath==0) then; call loadLayerFromPath( aLayer, vertPath, il )
            endif

            !--- OD along reference path for this sub-interval
            call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, &
                        refPathOD(il), NLTE=.false. )

            scalUpFac(   il) = viewPath%Wtot(il) / aLayer%Wtot
            scalDownFac( il) = downPath%Wtot(il) / aLayer%Wtot
            scalSolarFac(il) = solPath%Wtot(il)  / aLayer%Wtot
         enddo
      endif



      !--- OD along viewing path
      !
      if (UpLook) then
         L1 = obsLevel !Since uplooking viewPath contains layers from obsLevel to TOA only, obsLevel =1 for uplooking viewPath
         L2 = viewPath%nLay
      else
         L1 = 1
         L2 = obsLevel-1
      endif

      do il = L1, L2

         if ( .not.scalThermalPath .or. upLook ) then

            !--- Load layer temperature, pressure and concentration ...
            call loadLayerFromPath( aLayer, viewPath, il )
            !--- OD along upwelling path for this sub-interval
            call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, &
                        pathOD(il), NLTE=.false. )

         else
            call CLBLM_Spectrum_init( pathOD(il), refPathOD(il)%V1, &
                                                  refPathOD(il)%DV, &
                                                  refPathOD(il)%NLIM )
            if ( refPath==1) then
               pathOD(il) = refPathOD(il)
            else
               do iv = 1,pathOD(il)%NLIM
                  pathOD(il)%spect(iv) = refPathOD(il)%spect(iv) * scalUpFac(il)
               enddo
            endif
         endif
      enddo



      !--- OD along downwelling path
      if (.not.UpLook) then
         do il = 1,downPath%nLay

            if ( .not.scalThermalPath ) then

               call loadLayerFromPath( aLayer, downPath, il )

               call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, &
                           pathOD_dn(il), NLTE=.false. )
            else

               call CLBLM_Spectrum_init( pathOD_dn(il), refPathOD(il)%V1, &
                                                        refPathOD(il)%DV, &
                                                        refPathOD(il)%NLIM )
               do iv = 1,pathOD_dn(il)%NLIM
                  pathOD_dn(il)%spect(iv) = refPathOD(il)%spect(iv) * scalDownFac(il)
               enddo

               !/// pathNLTE_dn?
            endif
         enddo
      endif


      !--- OD along solar path
      if (.not.UpLook .and. SolarOn) then
         do il = 1,solPath%nLay

            if ( .not.scalSolarPath ) then

               call loadLayerFromPath( aLayer, solPath, il )

               call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, &
                           pathOD_sun(il), NLTE=.false. )
            else
               call CLBLM_Spectrum_init( pathOD_sun(il), refPathOD(il)%V1, &
                                                         refPathOD(il)%DV, &
                                                         refPathOD(il)%NLIM )
               do iv = 1,pathOD_sun(il)%NLIM
                  pathOD_sun(il)%spect(iv) = refPathOD(il)%spect(iv) * scalSolarFac(il)
               enddo
            endif
         enddo
      endif


      !--- If .not.LooiUp, pathOD in between obsLevel and TOA are not
      ! needed anymore, release the memory
      !///

   END SUBROUTINE


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
   SUBROUTINE drv_LevelTemp( path, pathOD, obsLevel, UpLook, &
                             odCtrl, subGrid, subDvCtrl, &
                             dRad_dTlev, dRad_dTau, dT_dTup, dT_dTlo, &
                             linInTau, &
                             V1B,V2B,DVB,NLIM, doPreBox, boxcarSize, boxcarHW, &
                             JacTemp )
!----------------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Config      ,ONLY: CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Path, &
                                    CLBLM_Layer, &
                                    loadLayerFromPath
      USE Module_ODLAY       ,ONLY: ODLAY
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      IMPLICIT NONE

      type(CLBLM_Path)       ,intent(in)    :: path
      type(CLBLM_Spectrum)   ,intent(in)    :: pathOD(:)
      integer                ,intent(in)    :: obsLevel
      logical                ,intent(in)    :: UpLook
      !logical                ,intent(in)    :: scalThmUpPath
      !real                   ,intent(in)    :: amfView
      type(CLBLM_OD_Ctrl)    ,intent(in)    :: odCtrl           !input ODLAY options
      type(CLBLM_SpectGrid)  ,intent(in)    :: subGrid          !input spectral grid info.
      type(CLBLM_DV_Ctrl)    ,intent(in)    :: subDvCtrl
      real                   ,intent(in)    :: dRad_dTlev(:,:)  ![NLIM,nLev]
      real                   ,intent(in)    :: dRad_dTau(:,:)   ![NLIM,nLay]
      real                   ,intent(in)    :: dT_dTup(:)       ![nLay] layer temperature derivatives with respect to upper level temperature
      real                   ,intent(in)    :: dT_dTlo(:)       ![nLay] layer temperature derivatives with respect to lower level temperature
      integer                ,intent(in)    :: linInTau         ! =0 linear-in-tau not used; =1 standard; =2 LBLRTM version
      real(r8)               ,intent(inout) :: V1B, V2B
      real                   ,intent(inout) :: DVB
      integer                ,intent(inout) :: NLIM
      logical                ,intent(in)    :: doPreBox
      integer                ,intent(in)    :: boxcarSize
      real                   ,intent(in)    :: boxcarHW
      type(CLBLM_Spectrum)   ,intent(inout) :: JacTemp(:) ![nLev]

      !--- Local variables
      integer                                  :: iv, il, L1,L2,lyrCnt
      type(CLBLM_Layer)                        :: aLayer
      real                                     :: TaveSave, WtotSave, fac
      type(CLBLM_Spectrum)                     :: lyrOD_m, lyrOD_p
      type(CLBLM_Spectrum)                     :: nlteEmisFac
      real               ,TARGET  ,allocatable :: tempLevDrv_A(:)
      real               ,TARGET  ,allocatable :: tempLevDrv_B(:)
      real ,dimension(:) ,POINTER              :: upLevDrv=>null()
      real ,dimension(:) ,POINTER              :: loLevDrv=>null()
      real                                     :: deltaT, deltaW, molAmt
      real                        ,allocatable :: dTau_dT(:) !i-th layer optical thickness derivatives with respect of temperature
      real                        ,allocatable :: dTau_dW(:) !i-th layer optical thickness derivatives with respect of Wtot
      real                                     :: dRaddT



      allocate( tempLevDrv_A(NLIM) )
      allocate( tempLevDrv_B(NLIM) )
      allocate( dTau_dT(NLIM) )
      allocate( dTau_dW(NLIM) )

      tempLevDrv_A(:) = 0.
      tempLevDrv_B(:) = 0.

      !--- Determine the number of layers to be calculated
      if (UpLook) then
         L1 = obsLevel !Since uplooking viewPath contains layers from obsLevel to TOA only, obsLevel =1 for uplooking viewPath
         L2 = path%nLay
      else
         L1 = 1
         L2 = path%nLay !obsLevel-1
      endif


      do il = L1,L2
         lyrCnt = il-L1+1

         !--- Switch pointers
         if (mod(lyrCnt,2)/=0) then
            upLevDrv => tempLevDrv_A
            loLevDrv => tempLevDrv_B
         else
            upLevDrv => tempLevDrv_B
            loLevDrv => tempLevDrv_A
         endif

         !--- Load layer temperature, pressure and concentration ...
         call loadLayerFromPath( aLayer, path, il )
         TaveSave = aLayer%Tave

         !!--- Calculate Tau for Tave-0.5K
         !deltaT = -0.5
         !aLayer%Tave =  aLayer%Tave + deltaT
         !call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, lyrOD_m, .false.,nlteEmisFac )
         !aLayer%Tave = TaveSave

         !--- Calculate Tau for original Tave
         ! Unless downlooking case and il>=obslevel, layer OD has been calcualted before.
         if (pathOD(il)%NLIM==0) then
            call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, lyrOD_m, .false.,nlteEmisFac )
         else
            lyrOD_m = pathOD(il)
         endif

         !--- Calculate Tau for Tave+1K
         deltaT = 1.0
         aLayer%Tave =  aLayer%Tave + deltaT
         call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, lyrOD_p, .false.,nlteEmisFac )
         aLayer%Tave = TaveSave

         !--- Calculate dTau_dT using finite difference
         do iv = 1,NLIM
            dTau_dT(iv) = lyrOD_p%spect(iv) - lyrOD_m%spect(iv)
         enddo


         !--- Calculate dTau_dW using finite difference
         !
         WtotSave =aLayer%Wtot
         deltaW = 0.01 * aLayer%Wtot
         aLayer%Wtot = aLayer%Wtot + deltaW
         call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, lyrOD_p, .false.,nlteEmisFac )
         aLayer%Wtot = WtotSave

         do iv = 1,NLIM
            dTau_dW(iv) = (lyrOD_p%spect(iv) - lyrOD_m%spect(iv)) / deltaW
         enddo


         !--- Adjust dTau_dT for change in Wtot
         fac = aLayer%Wtot/aLayer%Tave
         do iv = 1,NLIM
            dTau_dT(iv) = dTau_dT(iv) - dTau_dW(iv)*fac
         enddo


         !---Rad derivatives w.r.t level temperatures
         do iv = 1,NLIM
            dRaddT = dRad_dTau(iv,il) * dTau_dT(iv)
            upLevDrv(iv) =                dRaddT * dT_dTup(il)
            loLevDrv(iv) = loLevDrv(iv) + dRaddT * dT_dTlo(il) +  dRad_dTlev(iv,il)
         enddo

         if (il==L2) then
            do iv = 1,NLIM
               upLevDrv(iv) = upLevDrv(iv) + dRad_dTlev(iv,il+1)
            enddo
         endif


         !--- Load in output array
         !
         call output_JacArray( V1B,V2B,DVB,NLIM, &
                               doPreBox, boxcarSize, boxcarHW, &
                               loLevDrv(1:NLIM), JacTemp( lyrCnt ) )

         if (il==L2) then
            call output_JacArray( V1B,V2B,DVB,NLIM, &
                                  doPreBox, boxcarSize, boxcarHW, &
                                  upLevDrv(1:NLIM), JacTemp( lyrCnt+1 ) )
         endif

      enddo !do il = 1, nLay

   END SUBROUTINE

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
   SUBROUTINE drv_LevelConc( ajMolInd, path, obsLevel, UpLook, &
                             odCtrl, subGrid, subDvCtrl, &
                             dRad_dTau, dW_dQup, dW_dQlo, &
                             V1B,V2B,DVB,NLIM, doPreBox, boxcarSize, boxcarHW, &
                             JacMol )
!----------------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, molIndex
      USE Module_Config      ,ONLY: CLBLM_OD_Ctrl,&
                                    CLBLM_DV_Ctrl
      USE Module_LineData    ,ONLY: CLBLM_LineData
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum
      USE Module_AtmPath     ,ONLY: CLBLM_Path, &
                                    CLBLM_Layer, &
                                    loadLayerFromPath
      USE Module_ODLAY       ,ONLY: ODLAY
      USE Module_DV          ,ONLY: CLBLM_SpectGrid
      IMPLICIT NONE

      integer               ,intent(in)    :: ajMolInd(:)
      type(CLBLM_Path)      ,intent(in)    :: path
      integer               ,intent(in)    :: obsLevel
      logical               ,intent(in)    :: UpLook
      type(CLBLM_OD_Ctrl)   ,intent(in)    :: odCtrl         !input ODLAY options
      type(CLBLM_SpectGrid) ,intent(in)    :: subGrid        !input spectral grid info.
      type(CLBLM_DV_Ctrl)   ,intent(in)    :: subDvCtrl
      real                  ,intent(in)    :: dRad_dTau(:,:) ![NLIM,nLay]
      real                  ,intent(in)    :: dW_dQup(:,:)   ![nLay,nMol] layer molecular amount  derivatives with respect to upper level molecular mixing ratio
      real                  ,intent(in)    :: dW_dQlo(:,:)   ![nLay,nMol] layer molecular amount  derivatives with respect to lower level molecular mixing ratio
      real(r8)              ,intent(inout) :: V1B, V2B
      real                  ,intent(inout) :: DVB
      integer               ,intent(inout) :: NLIM
      logical               ,intent(in)    :: doPreBox
      integer               ,intent(in)    :: boxcarSize
      real                  ,intent(in)    :: boxcarHW
      type(CLBLM_Spectrum)  ,intent(inout) :: JacMol(:,:) ![nLev,nMol]

      !--- Local variables
      integer                                  :: iv, il, L1,L2,lyrCnt, im, imol
      type(CLBLM_Layer)                        :: aLayer
      type(CLBLM_Spectrum)                     :: lyrOD
      type(CLBLM_Spectrum)                     :: nlteEmisFac
      real               ,TARGET  ,allocatable :: tempLevDrv_A(:)
      real               ,TARGET  ,allocatable :: tempLevDrv_B(:)
      real ,dimension(:) ,POINTER              :: upLevDrv=>null()
      real ,dimension(:) ,POINTER              :: loLevDrv=>null()
      real                                     :: deltaT, molAmt
      real                        ,allocatable :: dTau_dW(:) !i-th layer optical thickness derivatives with respect of jth molecular amount
      real                        ,allocatable :: dRad_dW(:)
      character                                :: molName*20



      allocate( tempLevDrv_A(NLIM) )
      allocate( tempLevDrv_B(NLIM) )
      allocate( dTau_dW(NLIM) )
      allocate( dRad_dW(NLIM) )

      !--- Determine layers to be calcualted
      if (UpLook) then
         L1 = obsLevel !Since uplooking viewPath contains layers from obsLevel to TOA only, obsLevel=1 for uplooking viewPath
         L2 = path%nLay
      else
         L1 = 1
         L2 = path%nLay !obsLevel-1
      endif


      do im = 1, size(ajMolInd)

         imol = ajMolInd(im) !location of a.j. molecule in active molecular array
         molName = path%molID( imol )
         tempLevDrv_A(:) = 0.
         tempLevDrv_B(:) = 0.

         do il = L1,L2
            lyrCnt = il-L1+1

            !--- Switch pointers
            if (mod(lyrCnt,2)/=0) then
               upLevDrv => tempLevDrv_A
               loLevDrv => tempLevDrv_B
            else
               upLevDrv => tempLevDrv_B
               loLevDrv => tempLevDrv_A
            endif

            !--- Load layer temperature, pressure and concentration ...
            call loadLayerFromPath( aLayer, path, il )

            !--- set all layer amount other then the target one to zero
            molAmt = aLayer%W( imol )
            aLayer%W(1:aLayer%nMol) = 0.0
            aLayer%W(imol) = molAmt

            !--- Calculate Tau for perturbed path
            call ODLAY( aLayer, odCtrl, subGrid, subDvCtrl, &
                        lyrOD, .false.,nlteEmisFac, molName )

            !--- Calcualte absorption coefficient
            !  If W is in log scale, dTau/d(ln(W))= OD
            do iv = 1,NLIM
               dTau_dw(iv) = lyrOD%spect(iv) !/ molAmt
            enddo


            !--- From dRad_dTau to dRad_dW
            do iv = 1,NLIM
               dRad_dW(iv) = dRad_dTau(iv,il) * dTau_dW(iv) !W is in log scale.
            enddo

            !--- Derivative w.r.t level molecular concentration
            ! In dW_dQ, both W and Q are in log scale.
            do iv = 1,NLIM
               upLevDrv(iv) =                dRad_dW(iv) * dW_dQup(il,im)
               loLevDrv(iv) = loLevDrv(iv) + dRad_dW(iv) * dW_dQlo(il,im)
            enddo


            !--- Load into output array
            !
            call output_JacArray( V1B,V2B,DVB,NLIM, &
                                  doPreBox, boxcarSize, boxcarHW, &
                                  loLevDrv(1:NLIM), JacMol( lyrCnt, im ) )

            if (il==L2) then
               call output_JacArray( V1B,V2B,DVB,NLIM, &
                                     doPreBox, boxcarSize, boxcarHW, &
                                     upLevDrv(1:NLIM), JacMol( lyrCnt+1, im ) )
            endif


         enddo !do il = 1, nLay
      enddo !do im = 1, nAJMol

   END SUBROUTINE

!-----------------------------------------------------------------------
! subroutine to convert layer derivatives to level derivatives
!-----------------------------------------------------------------------
   SUBROUTINE level2layer( Plev,TLev,Tave,QLev, nMol,nLay, &
                           dT_dTup, dT_dTlo, dlnW_dlnQup, dlnW_dlnQlo )
!-----------------------------------------------------------------------
      real    ,intent(in)  :: PLev(:)           ![nLev]
      real    ,intent(in)  :: TLev(:)           ![nLev]
      real    ,intent(in)  :: Tave(:)           ![nLay]
      real    ,intent(in)  :: QLev(:,:)         ![nLev,nMol]
      real    ,intent(out) :: dT_dTup(:)        ![nLay]
      real    ,intent(out) :: dT_dTlo(:)        ![nLay]
      real    ,intent(out) :: dlnW_dlnQup(:,:)  ![nLay,nMol]  d(ln(W))/d(ln(Qup))
      real    ,intent(out) :: dlnW_dlnQlo(:,:)  ![nLay,nMol]  d(ln(W))/d(ln(Qlo))
      integer ,intent(in)  :: nMol, nLay


      integer           :: l,k
      real              :: Pup,Plo,Tup,Tlo,Tbar
      real, allocatable :: denUp(:),denLo(:)
      real              :: rhoU, rhoL, alpha, alphaT
      real              :: ratU, ratL

      allocate(denUp(nMol))
      allocate(denLo(nMol))


      DO l=1,nLay

         Pup  = Plev(l+1)
         Plo  = Plev(l)
         Tup  = Tlev(l+1)
         Tlo  = Tlev(l)
         Tbar = Tave(l)
         denUp(1:nMol) = QLev( l+1, 1:nMol )
         denLo(1:nMol) = QLev( l,   1:nMol )


         rhoU   = Pup/(Tup*1.3806503E-19)
         rhoL   = Plo/(Tlo*1.3806503E-19)
         alpha  = rhoU/rhoL
         alphaT = -(Tup-Tlo) / alog(alpha)

         !--- temperature
         dT_dTlo(l) = ( (Tbar-alphaT) / Tlo ) &
                     *( rhoL / (rhoL-rhoU) )  &
                     +( 1.0 - alphaT/Tlo ) / alog(alpha)

         dT_dTup(l) = ( (Tbar-alphaT) / Tup ) &
                     *( -rhoU / (rhoL-rhoU) ) &
                     -( 1.0 - alphaT/Tup ) / alog(alpha)

         !--- molecules
         do  k=1,nmol

            if ( denLo(k).ne.0.0 ) then

               ratU = denUp(k)/rhoU
               ratL = denLo(k)/rhoL

               dlnW_dlnQlo(l,k) = (ratL / (ratL-alpha*ratU)) &
                             +1.0 / alog(alpha*ratU/ratL)

               dlnW_dlnQup(l,k) = ( (-alpha*ratU) / (ratL-alpha*ratU) ) &
                             -1.0 / alog(alpha*ratU/ratL)
            else
               dlnW_dlnQlo(l,k)=0.0
               dlnW_dlnQup(l,k)=0.0
            endif

         enddo

      ENDDO !DO l=1,nLay

   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE output_SurfJac_Rad( V1B,V2B,DVB,NLIM, doPreBox, boxcarSize, boxcarHW, &
                                  dRad_dTsfc, dRad_dEmis, dRad_dRsfc, Rad, &
                                  doJacTskin, doJacEmis,  doJacRsfc,  doRad, &
                                  JacTskin,   JacEmis,    JacRsfc,    JacRad )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum
      IMPLICIT NONE

      real(r8)             ,intent(in)    :: V1B, V2B
      real                 ,intent(in)    :: DVB
      integer              ,intent(in)    :: NLIM
      logical              ,intent(in)    :: doPreBox
      integer              ,intent(in)    :: boxcarSize
      real                 ,intent(in)    :: boxcarHW
      real                 ,intent(in)    :: dRad_dTsfc(:)
      real                 ,intent(in)    :: dRad_dEmis(:)
      real                 ,intent(in)    :: dRad_dRsfc(:)
      real                 ,intent(in)    :: Rad(:)
      logical              ,intent(in)    :: doJacTskin, doJacEmis, doJacRsfc, doRad
      type(CLBLM_Spectrum) ,intent(inout) :: JacTskin  !A.J. w.r.t surface temperature
      type(CLBLM_Spectrum) ,intent(inout) :: JacEmis   !A.J. w.r.t surface emissivity
      type(CLBLM_Spectrum) ,intent(inout) :: JacRsfc   !A.J. w.r.t surface reflectance
      type(CLBLM_Spectrum) ,intent(inout) :: JacRad



      if (doJacTskin) then
         call output_JacArray( V1B,V2B,DVB,NLIM, &
                               doPreBox, boxcarSize, boxcarHW, &
                               dRad_dTsfc, JacTskin )
      endif


      if (doJacEmis) then
         call output_JacArray( V1B,V2B,DVB,NLIM, &
                               doPreBox, boxcarSize, boxcarHW, &
                               dRad_dEmis, JacEmis )
      endif


      if (doJacRsfc)  then
         call output_JacArray( V1B,V2B,DVB,NLIM, &
                               doPreBox, boxcarSize, boxcarHW, &
                               dRad_dRsfc, JacRsfc )
      endif


      if (doRad)  then
         call output_JacArray( V1B,V2B,DVB,NLIM, &
                               doPreBox, boxcarSize, boxcarHW, &
                               Rad, JacRad )
      endif

   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE output_JacArray( V1B,V2B,DVB,NLIM, &
                               doPreBox, boxcarSize, boxcarHW, &
                               JacArray, JacOut )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      USE Module_ScanFilter ,ONLY: clblm_SCANFN
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum, &
                                   resizeSpectrum
      IMPLICIT NONE

      real(r8)             ,intent(in)    :: V1B, V2B
      real                 ,intent(in)    :: DVB
      integer              ,intent(in)    :: NLIM
      logical              ,intent(in)    :: doPreBox
      integer              ,intent(in)    :: boxcarSize
      real                 ,intent(in)    :: boxcarHW
      real                 ,intent(in)    :: JacArray(:)
      type(CLBLM_Spectrum) ,intent(inout) :: JacOut

      real(r8)          :: V1out, V2out
      real              :: DVout
      integer           :: NPout
      real ,allocatable :: scanData(:)


      if ( doPreBox ) then

         V1out = V1B + (boxcarSize-1)*DVB/2.
         V2out = V2B - (boxcarSize-1)*DVB/2.
         DVout = (boxcarSize-1)*DVB

         call clblm_SCANFN( JacArray, V1B, V2B, DVB, NLIM, &
                            scanData, V1out,V2out,DVout,NPout, 1, boxcarHW )

         call resizeSpectrum( JacOut, scanData,V1out,DVout,NPout)
         !V1B = V1out
         !V2B = V2out
         !DVB = DVout
         !NLIM = NPout

         !if (NLIM >=boxcarSize) then
         !
         !   !call clblm_FFTSCN( JacArray, V1B, V2B, DVB, NLIM, 0., & !padding=0.
         !   !                   spV1want,spV2want,spDVwant, 1, boxcarHW, & !functID =1, preboxcar averaging
         !   !                   -1., .false., [0.,0.,0.] ) !The last three arguments are no use.
         !   call clblm_SCANFN( JacArray, V1B, V2B, DVB, NLIM, &
         !                      scanData, V1out,V2out,DVout,NPout, 1, boxcarHW )
         !
         !   call resizeSpectrum( JacOut, scanData,V1out,DVout,NPout)
         !   !V1B = V1out
         !   !V2B = V2out
         !   !DVB = DVout
         !   !NLIM = NPout
         !else
         !
         !   allocate(scanData(1))
         !   scanData(1) = sum( JacArray(1:NLIM) )/NLIM
         !
         !   V1B = V1B + DVout !2.*boxcarHW
         !   V2B = V1B
         !   NLIM = 1
         !   call resizeSpectrum( JacOut, scanData,V1B,DVB,NLIM)
         !
         !endif
      else
         call resizeSpectrum( JacOut, JacArray,V1B,DVB,NLIM)
      endif

   END SUBROUTINE

END MODULE
