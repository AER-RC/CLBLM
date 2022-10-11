!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_Drivers
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: RT_basic, &
             RT_Jac, &
             convolPreviousResults


   INTERFACE writeSpectFile
      module procedure writeSpectFile_1D_spect
      module procedure writeSpectFile_1D_array
      module procedure writeSpectFile_2D_spect
      module procedure writeSpectFile_2D_array
   END INTERFACE

CONTAINS !======================  MODULE CONTAINS ======================



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE RT_basic( Rad, Tx, OD, input_scene, &
                        path_flags, DV_flags, OD_flags, RT_flags, Flux_flags, &
                        VLo, VHi, DV, flxuOut, flxdOut, &
                        RadOut,nSampRad, TxOut,nSampTx, &
                        postProc_flags, RTgrid,nRTgrid, &
                        NLTE, grid_type, fileHdr )
!-----------------------------------------------------------------------
   USE Module_ConstParam   ,ONLY: r8=>kind_r8
   USE Module_ConstParam  ,ONLY: EPS, HEATFC
   USE Module_AtmPath      ,ONLY: CLBLM_Path, &
                                  CLBLM_PathCompound, &
                                  buildRTgrid, &
                                  calculatePath, &
                                  calcCLBLMPaths
   USE Module_Scene        ,ONLY: CLBLM_Scene, &
                                  CLBLM_Surface
   USE Module_Config       ,ONLY: CLBLM_Output_Ctrl, &
                                  CLBLM_Path_Ctrl, &
                                  CLBLM_OD_Ctrl, &
                                  CLBLM_DV_Ctrl, &
                                  CLBLM_RT_Ctrl, &
                                  CLBLM_Flux_Ctrl, &
                                  CLBLM_Post_Ctrl, &
                                  ioFiles, noPrnt
   USE Module_Spectrum     ,ONLY: CLBLM_Spectrum
   USE Module_DV           ,ONLY: setDV, &
                                  CLBLM_SpectGrid
   USE Module_noScattRT    ,ONLY: mode1_lookDn_mergeDn, &
                                  mode2_lookDn_mergeUp, &
                                  mode3_lookUp_mergeUP
   USE Module_PostProc     ,ONLY: postProc_scan, &
                                  expandBound, &
                                  CLBLM_FilterFunct, &
                                  readFilterFunct
   USE Module_Solar        ,ONLY: solarRadiance
   USE Module_FileIO       ,ONLY: CLBLM_SpectFileHeader
   IMPLICIT NONE


   integer                              ,intent(in)    :: Rad
   integer                              ,intent(in)    :: Tx
   integer                              ,intent(in)    :: OD
   type(CLBLM_Scene)                    ,intent(in)    :: input_scene
   type(CLBLM_Path_Ctrl)                ,intent(in)    :: path_flags
   type(CLBLM_DV_Ctrl)                  ,intent(in)    :: DV_flags
   type(CLBLM_OD_Ctrl)                  ,intent(in)    :: OD_flags
   type(CLBLM_RT_Ctrl)                  ,intent(in)    :: RT_flags
   type(CLBLM_Flux_Ctrl)                ,intent(in)    :: Flux_flags
   real(r8)                             ,intent(in)    :: VLo
   real(r8)                             ,intent(in)    :: VHi
   real                                 ,intent(in)    :: DV
   real          ,allocatable ,optional ,intent(out)   :: RadOut(:)
   real          ,allocatable ,optional ,intent(out)   :: TxOut(:,:)
   real          ,allocatable ,optional ,intent(out)   :: FlxuOut(:,:)
   real          ,allocatable ,optional ,intent(out)   :: FlxdOut(:,:)
   integer                    ,optional ,intent(out)   :: nSampTx, nSampRad
   type(CLBLM_Post_Ctrl)      ,optional ,intent(in)    :: postProc_flags
   real                       ,optional ,intent(in)    :: RTgrid(:)
   integer                    ,optional ,intent(in)    :: nRTgrid
   logical                    ,optional ,intent(in)    :: NLTE
   integer                    ,optional ,intent(in)    :: grid_type
   type(CLBLM_SpectFileHeader),optional ,intent(inout) :: fileHdr

!todo: output OD arrays. OD sizes are different for different layers.

   !---Local variables
   !
   character(*) ,parameter :: routineName = 'RT_basic'
   real         ,parameter :: TOL = 5.e-4

   integer                           :: il, is, np, nl, lyrNo, obsLevel
   logical                           :: NLTE_flag
   real(r8)                          :: V1,V2
   real                 ,allocatable :: zRT(:)
   integer                           :: zRTsize
   integer                           :: gridType
   integer                           :: functID
   real                              :: bound
   type(CLBLM_Path)                  :: path
   type(CLBLM_PathCompound)          :: paths
   type(CLBLM_Surface)               :: surf
   type(CLBLM_SpectGrid)             :: spGrid, spectGrid
   type(CLBLM_Spectrum)              :: solRadTOA
   type(CLBLM_Output_Ctrl)           :: outCtrl
   logical                           :: upLook,dnLook
   logical                           :: ODonly
   type(CLBLM_Spectrum)              :: totRad
   type(CLBLM_Spectrum) ,allocatable :: TxArray(:)
   type(CLBLM_Spectrum) ,allocatable :: ODarray(:)
   real                              :: padding
   type(CLBLM_Post_Ctrl)             :: postCtrl
   logical                           :: didPreBox
   integer                           :: nchan
   type(CLBLM_FilterFunct)           :: filtFunct
   real                 ,allocatable :: filtOut_rad(:)
   real                 ,allocatable :: filtOut_tx(:,:)
   real                 ,allocatable :: flxu(:,:) !upwelling flux
   real                 ,allocatable :: flxd(:,:) !downwelling flux
   real                 ,allocatable :: netflx(:,:) !net flux flux
   real                 ,allocatable :: htr(:,:) !heating rates
   real                 ,allocatable :: prethk(:) !pressure thickness
   real                 ,allocatable :: preslv(:) !pressure level
   real                 ,allocatable :: SfcRad(:,:) !downwelling surface radiances
   real                 ,allocatable :: SfcFlx(:) !Surface flux
   logical                           :: flux_flag
   character                         :: dim1Name*20, dim2Name*20, var1Name*20, var2Name*20, var3Name*20, var4Name*20, var5Name*20, var6Name*20, var7Name*20
   integer                           :: out, nout, k
   real, dimension(:), allocatable :: flux_bound
   integer, dimension(:), allocatable :: level
   real                              :: dv_flux


   flux_flag = flux_flags%flux_flag
   print *,'driver flux flag', flux_flag
   ODonly = Rad==0.and.Tx==0.and.OD/=0

   NLTE_flag = .FALSE.
   if (present(NLTE)) NLTE_flag = NLTE

   if ( ( Tx<0 .or. Rad<0 ) .and. .not.(present(postProc_flags)) ) then
      STOP '--- '//routineName//'(): Post process requested, but postProc_flags not present.'
   endif

   !--- Prepare path, downPath and solPath
   if (ODonly) then

      call buildRTgrid( input_scene%prfl, input_scene%geom, path_flags, zRT,zRTsize, RTgrid,nRTgrid )

      call calculatePath( path, input_scene%prfl, input_scene%geom, path_flags, zRT,zRTsize )

      obsLevel = -1
      do il =1,path%nRTlev
         if ( abs( path%zRT(il) - input_scene%geom%obs%obsAlt ) <TOL ) then
            obsLevel = il
            EXIT
         endif
      enddo

      paths%view = path
      paths%obsLev = obsLevel

   else
      CALL calcCLBLMPaths( input_scene, path_flags, paths, RTgrid,nRTgrid,flux_Flags )
   endif


   !--- Scale the path amount. This need the user to provide a subroutine "scalePathAmount".
   !call scalePathAmount_interface( path )

   !--- If entry for TBOUND < 0, use TZ(O) as boundary temperature. (From LBLATM.f90/FPACK())
   !IF (TBOUND.LT.0.) TBOUND = TZ(0)
   surf = input_scene%sfc
   IF (surf%Tskin.LT.0.) surf%Tskin = paths%view%T(1)


   !--- Set up spectral sampling interval(DV)
   ! Expand the spectral range if convolution is requested
   !
   V1 = VLo
   V2 = VHi

   if ( present(postProc_flags) .and. .not.ODonly ) then
      functID = postProc_flags%functID
      if ( functID >=1) then !For all non filtering options, expand the spectral limits
         bound = postProc_flags%HWHM * expandBound( functID )

         !--- If pre-boxcar smoothing is requested, make expansion for boxcar
         if (postProc_flags%boxcarHW >0.) then
            bound = bound + postProc_flags%boxcarHW
         endif

         V1 = VLo - bound
         V2 = VHi + bound
      endif

      !TODO expansion for filer functions
   endif


   !CALL setDV( spGrid1, paths%vert, V1,V2, DV_flags )
   !CALL setDV( spGrid2, paths%view,   V1,V2, DV_flags )
   !CALL setDV( spGrid3, paths%donw, V1,V2, DV_flags )
   !CALL setDV( spGrid4, paths%sun,  V1,V2, DV_flags )
   !minDVLoc = minloc( [spGrid1%DVnorm(1), spGrid2%DVnorm(1), spGrid3%DVnorm(1), spGrid4%DVnorm(1)] )
   !if     (minDVLoc ==1) then; spectGrid=spGrid1;
   !elseif (minDVLoc ==2) then; spectGrid=spGrid2;
   !elseif (minDVLoc ==2) then; spectGrid=spGrid3;
   !elseif (minDVLoc ==3) then; spectGrid=spGrid4; endif
   if (ODonly) then

      if (.not.present(grid_type)) STOP '--- '//routineName//'(): grid_type must be present for OD-only mode.'

      gridType = grid_type
      if (DV>0.) then
         CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2,DV )
      else
         CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2 )
      endif

   elseif (Rad<0 .or. Tx<0) then !results will be post processed later. use fixed-ratio DV grid

      gridType = 0 !Use fixed ratio DVs
      CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2 )

   else
      if (DV>0.) then !Uer user provided DV
         gridType = 1 !Uniform DV of user input value
 print *, 'uniform dv of user input value'
         CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2, DV )
 print *, 'DV', DV
      else
         gridType = 0 !Fixed ratio DV
         CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2 )
      endif
   endif



   !--- If solar on, read TOA solar radaiance data
   if ( RT_flags%SolarOn .and. .not.ODonly ) then
      call solarradiance( solRadTOA, V1,V2, &
                          RT_flags%JulDay, RT_flags%solarSource, RT_flags%solarVarOption, &
                          RT_flags%solarConst, RT_flags%SOLCYCFRAC, RT_flags%faculaVar, RT_flags%spotVar )
   endif


   !--- Call no-scattering RT solver
   !
   upLook = paths%view%geom%ANGLE <=90.
   dnLook = .not.upLook      !including both down-to-the-earth path and limb path.
   !limb   = paths%view%geom%Hend > paths%view%geom%HMIN

   !--- Allocate TxArray
   if (allocated(TxArray)) deallocate(TxArray)
   if ( abs(Tx) ==1 ) then
      allocate( TxArray( 1 ) );  !total Tx
   elseif (abs(Tx) ==2 ) then
      if (    upLook) then;  nl=paths%view%nRTLev - paths%obsLev
      elseif (dnLook) then;  nl=paths%obsLev-1
      endif
      allocate( TxArray( nl ) );  !Tx profile
   else
      allocate( TxArray(0)) !dummy array
   endif

   !--- Allocated ODarray
   if (allocated(ODarray)) deallocate(ODarray)
   if (OD/=0) then
      if (ODonly) then
         allocate( ODarray(paths%view%nLay) ) !all layers
      else
         if (    upLook) then;  nl=paths%view%nLay - paths%obsLev +1
         elseif (dnLook) then;  nl=paths%obsLev-1
         endif
         allocate( ODArray( nl ) ); !layers in between observer and surface/TOA
      endif
   else
      allocate( ODarray(0)) !dummy array
   endif


   !--- Call RT subroutines
   !
   outCtrl%Rad = Rad
   outCtrl%Tx  = Tx
   outCtrl%OD  = OD
 print *, 'Rad', Rad
   didPreBox = .FALSE.
   if (present(postProc_flags)) postCtrl = postProc_flags

   ! If running a radiative flux calculation, we will run the RT subroutine for the
   ! downwelling fluxes (Mode 1) and then upwelling fluxes (Mode 2)
   if (flux_flag .eqv. .true.) then
      print *, 'flux case'
         CALL mode1_lookDn_mergeDn( outCtrl, RT_flags, Flux_flags, input_scene,& !RT_flags, &
                                    paths, surf, solRadTOA, &
                                    OD_flags, spectGrid, DV_flags, &
                                    flxd, SfcRad, SfcFlx, totRad, TxArray, ODarray, &
                                    NLTE_flag, postCtrl, didPreBox )

         CALL mode2_lookDn_mergeUp( outCtrl, RT_flags, Flux_flags, & !RT_flags, &
                                    paths, surf, solRadTOA, &
                                    OD_flags, spectGrid, DV_flags, &
                                    flxu, SfcRad, totRad, ODarray, &
                                    NLTE_flag, postCtrl, didPreBox )
   elseIF  ( upLook ) then !Up looking cases    !if ( all(path%IPATH) ==3 )
  print *, 'up looking case'

      CALL mode3_lookUp_mergeUp( outCtrl, RT_flags, & !RT_flags, &
                                 paths%view, paths%obsLev, solRadTOA, &
                                 OD_flags, spectGrid, DV_flags, &
                                 totRad, TxArray, ODarray, &
                                 NLTE_flag, postCtrl, didPreBox)

   ELSEIF ( dnLook ) then

      if (Tx==0) then !Down-looking, request radiance. Mode2 doesn't have transmittance output.
  print *, 'down looking case, merge up'
         CALL mode2_lookDn_mergeUp( outCtrl, RT_flags, Flux_flags, & !RT_flags, &
                                    paths, surf, solRadTOA, &
                                    OD_flags, spectGrid, DV_flags, &
                                    flxu, SfcRad, totRad, ODarray, &
                                    NLTE_flag, postCtrl, didPreBox )

      else !Down-looking, request radiance and total transmittance or transmittance profile
  print *, 'down looking case, merge down'
         CALL mode1_lookDn_mergeDn( outCtrl, RT_flags, Flux_flags, input_scene, & !RT_flags, &
                                    paths, surf, solRadTOA, &
                                    OD_flags, spectGrid, DV_flags, flxd, &
                                    SfcRad, SfcFlx, totRad, TxArray, ODarray, &
                                    NLTE_flag, postCtrl, didPreBox )

      endif

   ENDIF


   !--- Post processing
   !

   if ( Tx<0 .or. Rad<0 ) then

      if (postProc_flags%functID /=0) then ! Scan the results

         if (Tx <0) then !Scan Tx
            padding = 1.
            do is = 1,size(TxArray)
               call postProcRTdata( TxArray(is), VLo,VHi,DV, &
                                    postProc_flags, padding, didPreBox )
            enddo
         endif

         if (Rad <0) then !Scan radiance
            padding = 0.
            call postProcRTdata( totRad, VLo,VHi,DV, &
                                 postProc_flags, padding, didPreBox )
         endif

      else !filter the resutls

         !Get filter function and allocate output array if filtering is requested.
         call readFilterFunct( trim(ioFiles%filterFunctPath)//&
                               trim(ioFiles%filterFunctFile), &
                               filtFunct )

         !--- Post process transmittance
         if (Tx <0) then !Scan Tx
            allocate( filtOut_tx( filtFunct%nChan, size(TxArray)  ) )
            padding = 1.
            do is = 1,size(TxArray)
               call postProcRTdata( TxArray(is), VLo,VHi,DV, &
                                    postProc_flags, padding, didPreBox, &
                                    filtFunct, filtOut_tx(:,is) )
            enddo
         endif

         if (Rad <0) then !Scan radiance
            allocate( filtOut_rad( filtFunct%nChan ) )
            padding = 0.
            call postProcRTdata( totRad, VLo,VHi,DV, &
                                 postProc_flags, padding, didPreBox, &
                                 filtFunct, filtOut_rad )
         endif

      endif !(postProc_flags%functID/=0)

   endif !( Tx<0 .or. Rad<0 )


   !--- Copy result to output array
   if (present(TxOut) .and. Tx/=0) then
      if (Tx>0 .or. postProc_flags%functID /=0) then !scanning
         np = TxArray(1)%NLIM
         nl = size(TxArray)
         allocate( TxOut(np,nl) )
         do is = 1,nl
            TxOut(:,is) = TxArray(is)%spect(:)
         enddo
      else !filtering
         allocate( TxOut(size(filtOut_tx,1),size(filtOut_tx,2)) )
         TxOut = filtOut_tx
      endif
   endif
   if (present(nSampTx))  nSampTx = size(Txout,1)

   !print *,'postproc flag functID', postProc_flags%functID
   !print *,'read filter function path', ioFiles%filterFunctPath
   !print *,'read filter function file', ioFiles%filterFunctFile
   !--- Copy radiance to output array
   if (present(RadOut) .and. Rad/=0) then
      if (Rad>0 .or. postProc_flags%functID /=0) then !scanning
         allocate( RadOut(totRad%NLIM) )
         RadOut(:) = totRad%spect(:)
      else !filtering
         allocate(RadOut(size(filtOut_rad)))
         RadOut = filtOut_rad
      endif
   endif
   if (present(nSampRad)) nSampRad = size(RadOut)

   !--- Copy flux to output array !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (present(flxuOut) .and. flux_flag .eqv. .true.) then
      dv_flux = flux_flags%dv_flux
      print *, 'P', input_scene%prfl%P
      print *, 'nlev', input_scene%prfl%nlev

      !print *,'sfc flux', SfcFlx
      OUT = (V2 - V1)/(DV_flux)
      NOUT = INT (OUT + EPS)
      allocate (flxuOut(input_scene%prfl%nlev,nout+1))
      allocate (flxdOut(input_scene%prfl%nlev,nout+1))
      allocate (netflx(input_scene%prfl%nlev,nout+1))
      allocate (htr(input_scene%prfl%nlev,nout+1))
      allocate (prethk(input_scene%prfl%nlev))
      allocate (preslv(input_scene%prfl%nlev))
      allocate (level(input_scene%prfl%nlev))
      flxuOut(1,:)=SfcFlx(:)
      flxdOut(input_scene%prfl%nlev,:)=0.
      level(1)=0
      do is=2,input_scene%prfl%nlev
         flxuOut(is,:)=flxu(is-1,:)
         level(is)=is-1
      enddo
      do is=1,input_scene%prfl%nlev-1
         flxdOut(is,:)=flxd(is,:)
      enddo
      do is=1,input_scene%prfl%nlev
         preslv(is)=input_scene%prfl%P(is)
      enddo
      print *,'level', level
      !flxdOut = flxd
      ! print *,'flxdOut', flxdOut
      OUT = (V2 - V1)/(DV_flux)
      NOUT = INT (OUT + EPS)
      allocate (flux_bound(nout))
      flux_bound(:)=0.
   ! For output, calculate wavenumber, net flux and heating rate
      DO K = 1, NOUT-1
         flux_BOUND(K) = V1 + dv_flux * FLOAT(K-1)
      enddo
      flux_BOUND(NOUT) = V2

!Compute net fluxes and heating rates, then output fluxes and
!heating rates from top of atmosphere down for each level.

      DO K = 1, NOUT
         DO  il = input_scene%prfl%nlev,1,-1
            netflx(il,K) = flxuOut(il,K) - flxdOut(il,K)
            IF (il.EQ.input_scene%prfl%nlev) THEN
               HTR(il,K) = 0.
!               PRESLV(il) = 0.
               else
               HTR(il,K) = NETFLX(il,K) - NETFLX(il+1,K)
               PRETHK(il) = preslv(il) - preslv(il+1)
               HTR(il,K) = HEATFC * HTR(il,K) / PRETHK(il)
            ENDIF
         enddo
      enddo

         print *,'write fluxes'
         dim1Name='levels'
         dim2Name='fluxSpectralBins'
         var1Name='upwellingFluxes'
         var2Name='SpectralBinBoundary'
         var3Name='downwellingFluxes'
         var4Name='Pressure'
         var5Name='Level'
         var6Name='NetFlux'
         var7Name='HeatingRate'
         fileHdr%V1               = V1
         fileHdr%V2               = V2
         fileHdr%DV               = DV_flux
      call  writeSpectFile_2D_flux(flxuOut,flux_bound,flxdOut,input_scene%prfl%P,level,netflx,htr, [var1Name,var2Name,var3Name,var4Name,var5Name,var6Name,var7Name], [dim1name,dim2name], fileHdr)
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !--- Write out results
   !
   !fileHdr%sceneFile         = trim(sceneFile)         !filled already
   !fileHdr%fileID            = UID                     !filled already
   !fileHdr%sceneNo           = sceneNo                 !filled already
   !fileHdr%productName                                 !to be filled later
   !fileHdr%procuctSpectType                            !to be filled later
   !fileHdr%V1                                          !to be filled later
   !fileHdr%V2                                          !to be filled later
   !fileHdr%DV                                          !to be filled later
   !fileHdr%nSamp                                       !to be filled later
   !fileHdr%filterFunctFile   = filterFunctFile         !filled already
   !fileHdr%convolParam       = postCtrl                !filled already
   !fileHdr%obsAlt                                      !filled already
   !fileHdr%viewAng                                     !filled already

   if (.not.present(fileHdr)) then
      print*, '--- '//routineName//'(): fileHdr structure is not present, no results will be written.'
   endif

   if (present(fileHdr)) then
      if (flux_flag .eqv. .false.) then
         if (Rad<0 .and. postProc_flags%functID==0) then !Filtered result
            call writeRadiance( fileHdr, Rad, filtOut=filtOut_rad )
         elseif (Rad/=0) then !mono or scanned result
 print *, 'write Radiance mono or scanned'
 print *, 'Rad', Rad
            call writeRadiance( fileHdr, Rad, totRad )
         endif
      endif

      if (Tx<0 .and. postProc_flags%functID==0) then !Filtered result
         call writeTransmittance( fileHdr, Tx, filtOut=filtOut_tx )
      elseif (Tx/=0) then !mono or scanned result
         call writeTransmittance( fileHdr, Tx, TxArray )
      endif


      if (OD /=0) then
         do il = 1,size(ODarray)

            if (upLook) then
               lyrNo = paths%obsLev + il-1
            elseif (dnLook) then
               if (Tx==0) then !merge-up mode called
                  lyrNo = il
               else !merge down mode callded
                  if (ODonly) then
                     lyrNo = paths%view%nLay - il+1
                  else
                     lyrNo = paths%obsLev-il
                  endif
               endif
            endif

            call writeOD( fileHdr, ODarray(il), lyrNo )
         enddo
      endif

   endif ! present(fileHdr)

END SUBROUTINE



!-----------------------------------------------------------------------
! Rad:
!  = 0, No radiance calculated
!  = 1, Returns monochromatic radiances
!  =-1, Returns convolved radiances
! Jac:
!  = 1, Returns monochromatic Jacobians
!  =-1, Returns convolved Jacobians
!
!-----------------------------------------------------------------------
   SUBROUTINE RT_Jac( Rad, Jac, Jacobians_list, input_scene, &
                      path_flags, DV_flags, OD_flags, RT_flags, Flux_flags,&
                      VLO, VHi, DV, &
                      JacOut_mol, JacOut_temp, JacOut_Tskin, JacOut_emis, JacOut_Rsfc, RadOut, &
                      nSampJac, postProc_flags, RTgrid,nRTgrid, fileHdr )
!-----------------------------------------------------------------------
   USE Module_ConstParam   ,ONLY: r8=>kind_r8
   USE Module_Utility      ,ONLY: upper
   USE Module_AtmPath      ,ONLY: CLBLM_Path, &
                                  CLBLM_PathCompound, &
                                  calcCLBLMPaths
   USE Module_Scene        ,ONLY: CLBLM_Scene, &
                                  CLBLM_Surface
   USE Module_Config       ,ONLY: CLBLM_Path_Ctrl, &
                                  CLBLM_OD_Ctrl, &
                                  CLBLM_DV_Ctrl, &
                                  CLBLM_RT_Ctrl, &
                                  CLBLM_Flux_Ctrl, &
                                  CLBLM_Post_Ctrl, &
                                  ioFiles, noPrnt
   USE Module_Spectrum     ,ONLY: CLBLM_Spectrum
   USE Module_DV           ,ONLY: setDV, &
                                  CLBLM_SpectGrid
   USE Module_noScattJac   ,ONLY: noScattJacob
   USE Module_PostProc     ,ONLY: postProc_scan, &
                                  expandBound, &
                                  CLBLM_FilterFunct, &
                                  readFilterFunct
   USE MODULE_Solar        ,ONLY: solarradiance
   USE Module_FileIO       ,ONLY: CLBLM_SpectFileHeader
   IMPLICIT NONE


   integer                              ,intent(in)    :: Rad
   integer                              ,intent(in)    :: Jac
   character(*)                         ,intent(in)    :: Jacobians_list(:)
   type(CLBLM_Scene)                    ,intent(in)    :: input_scene
   type(CLBLM_Path_Ctrl)                ,intent(in)    :: path_flags
   type(CLBLM_DV_Ctrl)                  ,intent(inout) :: DV_flags !gridType may be changed
   type(CLBLM_OD_Ctrl)                  ,intent(in)    :: OD_flags
   type(CLBLM_RT_Ctrl)                  ,intent(inout) :: RT_flags
   real(r8)                             ,intent(in)    :: VLO
   real(r8)                             ,intent(in)    :: VHi
   real                                 ,intent(in)    :: DV
   real          ,allocatable ,optional ,intent(out)   :: JacOut_mol(:,:,:)  ![NLIM,nLay,nMol], A.J. w.r.t molecular contentration
   real          ,allocatable ,optional ,intent(out)   :: JacOut_temp(:,:)   ![NLIM,nLay],      A.J. w.r.t air temperature
   real          ,allocatable ,optional ,intent(out)   :: JacOut_Tskin(:)    ![NLIM],           A.J. w.r.t surface temperature
   real          ,allocatable ,optional ,intent(out)   :: JacOut_emis(:)     ![NLIM],           A.J. w.r.t surface emissivity
   real          ,allocatable ,optional ,intent(out)   :: JacOut_Rsfc(:)     ![NLIM],           A.J. w.r.t surface reflectance
   real          ,allocatable ,optional ,intent(out)   :: radOut(:)
   integer                    ,optional ,intent(out)   :: nSampJac
   type(CLBLM_Post_Ctrl)      ,optional ,intent(in)    :: postProc_flags
   real                       ,optional ,intent(in)    :: RTgrid(:)
   integer                    ,optional ,intent(in)    :: nRTgrid
   type(CLBLM_SpectFileHeader),optional ,intent(inout) :: fileHdr
   type(CLBLM_Flux_Ctrl)                ,intent(in)    :: Flux_flags



   !---Local variables
   !
   character(*) ,parameter :: routineName = 'RT_Jac'
   real         ,parameter :: TOL = 5.e-4

   type(CLBLM_Spectrum) ,allocatable :: JacMol(:,:)  ![nLay,nMol], A.J. w.r.t molecular contentration
   type(CLBLM_Spectrum) ,allocatable :: JacTemp(:)   ![nLay],      A.J. w.r.t air temperature
   type(CLBLM_Spectrum)              :: JacTskin     !             A.J. w.r.t surface temperature
   type(CLBLM_Spectrum)              :: JacEmis      !             A.J. w.r.t surface emissivity
   type(CLBLM_Spectrum)              :: JacRsfc      !             A.J. w.r.t surface reflectance
   type(CLBLM_Spectrum)              :: JacRad
   type(CLBLM_Post_Ctrl)             :: postCtrl
   type(CLBLM_FilterFunct)           :: filtFunct
   real                 ,allocatable :: filtOut_mol(:,:,:)
   real                 ,allocatable :: filtOut_temp(:,:)
   real                 ,allocatable :: filtOut_Tskin(:)
   real                 ,allocatable :: filtOut_emis(:)
   real                 ,allocatable :: filtOut_Rsfc(:)
   real                 ,allocatable :: filtOut_rad(:)

   integer                    :: il, im, np,nl
   integer                    :: gridType
   real(r8)                   :: V1,V2
   integer                    :: functID
   real                       :: bound
   type(CLBLM_PathCompound)   :: paths
   type(CLBLM_Surface)        :: surf
   type(CLBLM_SpectGrid)      :: spectGrid
   type(CLBLM_Spectrum)       :: solRadTOA
   integer                    :: nAJMol
   character(20) ,allocatable :: ajMolNames(:), tempMolNames(:)
   logical                    :: doJacTskin, doJacEmis, doJacRsfc, &
                                 doJacTemp, doJacMol, doRad
   real                       :: padding
   logical                    :: didPreBox
   logical                    :: flux_flag




   !--- Prepare path, downPath and solPath
   CALL calcCLBLMPaths( input_scene, path_flags, paths, RTgrid,nRTgrid,flux_Flags )


   !--- Scale the path amount. This need the user to provide a subroutine "scalePathAmount".
   !call scalePathAmount_interface( path )

   !--- If entry for TBOUND < 0, use TZ(1) as boundary temperature. (From LBLATM.f90/FPACK())
   !IF (TBOUND.LT.0.) TBOUND = TZ(1)
   surf = input_scene%sfc
   IF (surf%Tskin.LT.0.) surf%Tskin = paths%view%T(1)


   !--- Set up spectral sampling interval(DV)
   !
   V1 = VLo
   V2 = VHi
   !--- Expand the spectral range if convolution is requested
   if (present(postProc_flags)) then

      postCtrl = postProc_flags

      functID = postCtrl%functID
      if ( functID >=1 ) then !For all non filtering options, expand the spectral limits
         bound = postCtrl%HWHM * expandBound( functID )

         !--- If pre-boxcar smoothing is requested, make expansion for boxcar
         if (postCtrl%boxcarHW >0.) then
            bound = bound + postCtrl%boxcarHW
         endif

         V1 = VLo - bound
         V2 = VHi + bound
      endif

      !TODO expansion for filer function
   endif


   !--- Set up spectral sampling interval(DV)
   !
   !CALL setDV( spGrid1, paths%vert, V1,V2, DV_flags )
   !CALL setDV( spGrid2, paths%view, V1,V2, DV_flags )
   !CALL setDV( spGrid3, paths%donw, V1,V2, DV_flags )
   !CALL setDV( spGrid4, paths%sun,  V1,V2, DV_flags )
   !minDVLoc = minloc( [spGrid1%DVnorm(1), spGrid2%DVnorm(1), spGrid3%DVnorm(1), spGrid4%DVnorm(1)] )
   !if     (minDVLoc ==1) then; spectGrid=spGrid1;
   !elseif (minDVLoc ==2) then; spectGrid=spGrid2;
   !elseif (minDVLoc ==2) then; spectGrid=spGrid3;
   !elseif (minDVLoc ==3) then; spectGrid=spGrid4; endif
   !
   gridType = 1 !User uniform DV grid
   if ( DV>0 .and. .not.(Jac<0 .or. Rad<0) ) then  !uniform grid DV=userDV when post processing not requested.
      CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2, DV )
   else !uniform grid DV=finest DV
      CALL setDV( spectGrid, paths%view, DV_flags, gridType, V1,V2 )
   endif


   !--- If solar on, read TOA solar radaiance data
   if ( RT_flags%SolarOn ) then
      call solarradiance( solRadTOA, V1,V2, &
                          RT_flags%JulDay, RT_flags%solarSource, RT_flags%solarVarOption, &
                          RT_flags%solarConst, RT_flags%SOLCYCFRAC, RT_flags%faculaVar, RT_flags%spotVar )
   endif


   !--- Set A.J. indicators and extract the molecular names from JacList array, and
   !  get indexes of A.J. molecules in path%molID array.
   !
   doJacTskin = .FALSE.
   doJacEmis  = .FALSE.
   doJacRsfc  = .FALSE.
   doRad      = .FALSE.
   doJacTemp  = .FALSE.
   doJacMol   = .FALSE.
   nAJMol     = 0

   if (Rad/=0)  doRad = .TRUE.

   if (Jac/=0) then
      allocate(tempMolNames(size(Jacobians_list)))
      do im = 1,size(Jacobians_list)
         select case ( upper(trim(adjustl( Jacobians_list(im) ))) )
         case ("T")     ;doJacTemp  = .TRUE.
         case ('TSKIN') ;doJacTskin = .TRUE.;
         case ('EMIS')  ;doJacEmis  = .TRUE.;
         case ('RSFC')  ;if (RT_flags%SolarOn) doJacRsfc = .TRUE.; !must be SolarOn
         case default
            nAJMol = nAJMol+1
            tempMolNames(nAJMol) = Jacobians_list(im)
         end select
      enddo

      if (nAJMol>0) then
         allocate( ajMolNames( nAJMol))
         do im = 1,nAJMol
            ajMolNames(im) = tempMolNames(im)
         enddo
         doJacMol = .TRUE.
      endif
      deallocate(tempMolNames)
   endif

   if (.not.any([doRad,doJacTskin,doJacEmis,doJacRsfc,doJacTemp,doJacMol]) )  RETURN !nothing to do, return.


   !--- Call Jacobians solver
   !
   allocate( JacMol(  paths%view%nLay+1,nAJMol ))
   allocate( JacTemp( paths%view%nLay+1 ))

   CALL noScattJacob( paths, surf, spectGrid, OD_flags, DV_flags, postCtrl, &
                      RT_flags%Solaron, solRadTOA,  RT_flags%linInTau, ajMolNames, &
                      doJacMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, doRad, &
                      JacMol, JacTemp, JacTskin, JacEmis, JacRsfc, JacRad, didPreBox )



   !--- Post process
   !
   !Get filter function if filtering is requested.
   if ((Jac<0 .or. Rad<0) .and. postCtrl%functID==0) then
      call readFilterFunct( trim(ioFiles%filterFunctPath)//&
                            trim(ioFiles%filterFunctFile), &
                            filtFunct )
   endif

   if (Jac<0) then
      if (postCtrl%functID /=0) then !Scanning

         padding = 0.

         if (doJacMol) then
            do im =1,nAJMol
            do il =1,paths%view%nLay+1
               call postProcRTdata( JacMol(il,im), VLo,VHi,DV, postCtrl, padding, didPreBox )
            enddo
            enddo
         endif

         if (doJacTemp) then
            do il =1,paths%view%nLay+1
               call postProcRTdata( JacTemp(il), VLo,VHi,DV, postCtrl, padding, didPreBox )
            enddo
         endif

         if (doJacTskin) call postProcRTdata( JacTskin, VLo,VHi,DV, postCtrl, padding, didPreBox )

         if (doJacEmis)  call postProcRTdata( JacEmis, VLo,VHi,DV, postCtrl, padding, didPreBox )

         if (doJacRsfc)  call postProcRTdata( JacRsfc, VLo,VHi,DV, postCtrl, padding, didPreBox )

      else  !Filtering

         padding = 0.

         if (doJacMol) then
            allocate( filtOut_mol( filtFunct%nChan, paths%view%nLay+1, nAJMol ) )
            do im =1,nAJMol
            do il =1,paths%view%nLay+1
               call postProcRTdata( JacMol(il,im), VLo,VHi,DV, &
                                    postCtrl, padding, didPreBox, &
                                    filtFunct, filtOut_mol(:,il,im) )
            enddo
            enddo
         endif

         if (doJacTemp) then
            allocate( filtOut_temp( filtFunct%nChan, paths%view%nLay+1 ) )
            do il =1,paths%view%nLay+1
               call postProcRTdata( JacTemp(il), VLo,VHi,DV, &
                                    postCtrl, padding, didPreBox, &
                                    filtFunct, filtOut_temp(:,il) )
            enddo
         endif


         if (doJacTskin) then
            allocate( filtOut_Tskin( filtFunct%nChan ) )
            call postProcRTdata( JacTskin, VLo,VHi,DV, &
                                 postCtrl, padding, didPreBox, &
                                 filtFunct, filtOut_Tskin )
         endif

         if (doJacEmis) then
            allocate( filtOut_emis( filtFunct%nChan ) )
            call postProcRTdata( JacEmis, VLo,VHi,DV, &
                                 postCtrl, padding, didPreBox,&
                                 filtFunct, filtOut_emis )
         endif

         if (doJacRsfc) then
            allocate( filtOut_Rsfc( filtFunct%nChan ) )
            call postProcRTdata( JacRsfc, VLo,VHi,DV, &
                                 postCtrl, padding, didPreBox, &
                                 filtFunct, filtOut_Rsfc )
         endif

      endif !postCtrl%functID /=0
   endif !if (Jac<0) then


   if (Rad<0 .and. doRad)  then
      padding = 0.

      if (postCtrl%functID/=0) then
         call postProcRTdata( JacRad, VLo,VHi,DV, postCtrl, padding, didPreBox )
      else
         allocate( filtOut_rad( filtFunct%nChan ) )
         call postProcRTdata( JacRad, VLo,VHi,DV, &
                              postCtrl, padding, didPreBox, &
                              filtFunct, filtOut_rad )
      endif
   endif



   !--- Copy result to output arrays
   !
   if (present(JacOut_mol  ) .and. doJacMol  ) then
      if (Jac>0 .or. postCtrl%functID /=0) then
         np = JacMol(1,1)%NLIM
         nl = size(JacMol,1)
         allocate( JacOut_mol( np,nl,nAJMol ) )
         do im =1,nAJMol
         do il =1,nl
            JacOut_mol(1:np,il,im) = JacMol(il,im)%spect(1:np)
         enddo
         enddo
      else
         allocate( JacOut_mol( size(filtOut_mol,1), size(filtOut_mol,2), size(filtOut_mol,3)) )
         JacOut_mol = filtOut_mol
      endif
   endif

   if (present(JacOut_temp ) .and. doJacTemp ) then
      if (Jac>0 .or. postCtrl%functID /=0) then
         np = JacTemp(1)%NLIM
         nl = size(JacTemp)
         allocate( JacOut_temp( np,nl ) )
         do il = 1,nl
            JacOut_temp(1:np,il) = JacTemp(il)%spect(1:np)
         enddo
      else
         allocate( JacOut_temp( size(filtOut_temp,1), size(filtOut_temp,2) ) )
         JacOut_temp = filtOut_temp
      endif
   endif

   if (present(JacOut_Tskin) .and. doJacTskin) then
      if (Jac>0 .or. postCtrl%functID /=0) then
         np = JacTskin%NLIM
         allocate( JacOut_Tskin( np ) )
         JacOut_Tskin(1:np) = JacTskin%spect(1:np)
      else
         allocate( JacOut_Tskin( size(filtOut_Tskin)) )
         JacOut_Tskin = filtOut_Tskin
      endif
   endif

   if (present(JacOut_emis ) .and. doJacEmis ) then
      if (Jac>0 .or. postCtrl%functID /=0) then
         np = JacEmis%NLIM
         allocate( JacOut_emis( np ) )
         JacOut_emis(1:np) = JacEmis%spect(1:np)
      else
         allocate( JacOut_emis( size(filtOut_emis)) )
         JacOut_emis = filtOut_emis
      endif
   endif

   if (present(JacOut_Rsfc ) .and. doJacRsfc ) then
      if (Jac>0 .or. postCtrl%functID /=0) then
         np = JacRsfc%NLIM
         allocate( JacOut_Rsfc( np ) )
         JacOut_Rsfc(1:np) = JacRsfc%spect(1:np)
      else
         allocate( JacOut_Rsfc( size(filtOut_Rsfc)))
         JacOut_Rsfc = filtOut_Rsfc
      endif
   endif

   if (present(RadOut) .and. doRad ) then
      if (Jac>0 .or. postCtrl%functID /=0) then
         np = JacRad%NLIM
         allocate( RadOut( np ) )
         RadOut(1:np) = JacRad%spect(1:np)
      else
         allocate( RadOut( size(filtOut_rad)) )
         RadOut = filtOut_rad
      endif
   endif

   if (present(nSampJac)) then
      if (Jac>0 .or. postCtrl%functID/=0) then
         nSampJac=np
      else
         nSampJac=filtFunct%nChan
      endif
   endif



   !--- Write out results
   !
   !fileHdr%sceneFile         = trim(sceneFile)         !filled already
   !fileHdr%fileID            = UID                     !filled already
   !fileHdr%sceneNo           = sceneNo                 !filled already
   !fileHdr%productName                                 !to be filled later
   !fileHdr%procuctSpectType                            !to be filled later
   !fileHdr%V1                                          !to be filled later
   !fileHdr%V2                                          !to be filled later
   !fileHdr%DV                                          !to be filled later
   !fileHdr%nSamp                                       !to be filled later
   !fileHdr%filterFunctFile   = filterFunctFile         !filled already
   !fileHdr%convolParam       = postCtrl                !filled already
   !fileHdr%obsAlt                                      !filled already
   !fileHdr%viewAng                                     !filled already

   if (.not.present(fileHdr)) then
      print*, '--- '//routineName//'(): fileHdr structure is not present, no results will be written.'
   endif

   if (present(fileHdr)) then

      if (Jac<0 .and. postProc_flags%functID==0) then !Filtered result
         call writeJacobians( fileHdr, JacMol, JacTemp, JacTskin, JacEmis, JacRsfc, &
                              doJacMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, &
                              Jac, Jacobians_list, input_scene%prfl%molID, &
                              filtOut_mol, filtOut_temp, filtOut_Tskin, filtOut_emis, filtOut_Rsfc )
      else
         call writeJacobians( fileHdr, JacMol, JacTemp, JacTskin, JacEmis, JacRsfc, &
                              doJacMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, &
                              Jac, Jacobians_list, input_scene%prfl%molID  )
      endif


      if (doRad) then
         if (Rad<0 .and. postProc_flags%functID==0) then !Filtered result
            call writeRadiance( fileHdr, Rad, filtOut=filtOut_rad, extStr='JacRad' )
         else
            call writeRadiance( fileHdr, Rad, JacRad, extStr='JacRad' )
         endif
      endif

   endif

   END SUBROUTINE



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE convolPreviousResults( inputFile_Rad, &
                                     inputFile_TxTot, &
                                     inputFile_TxPrfl, &
                                     inputFile_Jac, &
                                     postCtrl, VLo,VHi,DV )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, FILLREAL, FILLINT
      USE Module_Config      ,ONLY: CLBLM_Post_Ctrl, &
                                    ioFiles
      USE Module_Utility     ,ONLY: splitFilePath, &
                                    find_replace_string, &
                                    upper
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum
      USE Module_FileIO      ,ONLY: CLBLM_SpectFileHeader
      USE Module_PostProc    ,ONLY: postProc_scan, &
                                    postProc_filter, &
                                    CLBLM_FilterFunct,&
                                    readFilterFunct
      IMPLICIT NONE

      character(*)          ,intent(in) :: inputFile_Rad, &
                                           inputFile_TxTot, &
                                           inputFile_TxPrfl, &
                                           inputFile_Jac
      type(CLBLM_Post_Ctrl) ,intent(in) :: postCtrl
      real(r8)              ,intent(in) :: VLo,VHi
      real                  ,intent(in) :: DV

      character(*) ,parameter :: routineName = 'convolPreviousResults'

      integer                           :: i1,i2, il, ivar
      character                         :: jacStr*10, varName*20, dim1Name*20, dim2Name*20
      character                         :: dirName*200, baseName*128, inFile*256, outFile*256
      real(r8)                          :: V1O,V2O
      real                              :: DVO
      real                              :: padding
      type(CLBLM_Spectrum)              :: spect
      type(CLBLM_Spectrum) ,allocatable :: spectArr(:)
      type(CLBLM_SpectFileHeader)       :: fileHdr
      type(CLBLM_FilterFunct)           :: filtFunct
      real                 ,allocatable :: filtOut_1D(:)
      real                 ,allocatable :: filtOut_2D(:,:)


      !fileHdr%sceneFile         = trim(sceneFile)         !filled already
      !fileHdr%fileID            = UID                     !filled already
      !fileHdr%sceneNo           = sceneNo                 !filled already
      !fileHdr%productName                                 !need to be filled here
      !fileHdr%procuctSpectType                            !need to be filled here
      !fileHdr%V1                                          !need to be filled here
      !fileHdr%V2                                          !need to be filled here
      !fileHdr%DV                                          !need to be filled here
      !fileHdr%nSamp                                       !need to be filled here
      !fileHdr%filterFunctFile   = filterFunctFile         !need to be filled here
      !fileHdr%convolParam       = postCtrl                !need to be filled here
      !fileHdr%obsAlt                                      !filled already
      !fileHdr%viewAng                                     !filled already

      DO ivar = 1,3

         if (ivar==1) inFile = inputFile_Rad
         if (ivar==2) inFile = inputFile_TxTot
         if (ivar==3) inFile = inputFile_TxPrfl
         if (inFile=='') CYCLE

         if (ivar==1) call readSpectFile_real_1D( inFile, 'radiance',      'numPoints', spect, fileHdr )
         if (ivar==2) call readSpectFile_real_1D( inFile, 'transmittance', 'numPoints', spect, fileHdr )
         if (ivar==3) call readSpectFile_real_1D( inFile, 'transmittance', 'numPoints', spect, fileHdr )

         if (ivar==1) padding = 0.
         if (ivar==2) padding = 1.
         if (ivar==3) padding = 1.

         !--- Process spetral data
         !
         if (postCtrl%functID == 0) then !Filter using instrument SRF

            !Get filter function if filtering is requested.
            call readFilterFunct( trim(ioFiles%filterFunctPath)//&
                                  trim(ioFiles%filterFunctFile), &
                                  filtFunct )

            if (allocated(filtOut_1D)) deallocate(filtOut_1D)
            allocate( filtOut_1D( filtFunct%nChan ) )
            call postProc_filter( spect, filtFunct, filtOut_1D)

            fileHdr%productSpectType = 'Filtered'
            fileHdr%V1               = FILLREAL
            fileHdr%V2               = FILLREAL
            fileHdr%DV               = FILLREAL
            fileHdr%nSamp            = FILLINT
            fileHdr%filterFunctFile  = ioFiles%filterFunctFile
            fileHdr%convolParam      = postCtrl

            outFile = inFile
            if ( .not.find_replace_string( outFile, '_mono', '_filtered' )) STOP '--- '//routineName//'(): Not a monochromatic spectral file.'

            call splitFilePath( outFile, dirName, basename)

            if (ivar==1 .and. ioFiles%outPath_Rad    /='') outFile = trim(ioFiles%outPath_Rad   )//trim(baseName)
            if (ivar==2 .and. ioFiles%outPath_TxTot  /='') outFile = trim(ioFiles%outPath_TxTot )//trim(baseName)
            if (ivar==3 .and. ioFiles%outPath_TxPrfl /='') outFile = trim(ioFiles%outPath_TxPrfl)//trim(baseName)

            if (ivar==1) call writeSpectFile( filtOut_1D, outfile, 'radiance',      'numPoints', fileHdr)
            if (ivar==2) call writeSpectFile( filtOut_1D, outfile, 'transmittance', 'numPoints', fileHdr)
            if (ivar==3) call writeSpectFile( filtOut_1D, outfile, 'transmittance', 'numPoints', fileHdr)

         else !Convolve with spectrally invariant scaning function

            V1O = VLo
            V2O = VHi
            DVO = DV

            call postProc_scan( spect, V1O,V2O,DVO, postCtrl, padding )

            fileHdr%productSpectType = 'Convolved'
            fileHdr%V1               = spect%V1
            fileHdr%V2               = spect%V2
            fileHdr%DV               = spect%DV
            fileHdr%nSamp            = spect%NLIM
            !fileHdr%filterFunctFile   = filterFunctFile
            fileHdr%convolParam       = postCtrl

            outFile = inFile
            if ( .not.find_replace_string( outFile, '_mono', '_convolved' )) STOP '--- '//routineName//'(): Not a monochromatic spectral file.'

            call splitFilePath( outFile, dirName, basename)

            if (ivar==1 .and. ioFiles%outPath_Rad    /='') outFile = trim(ioFiles%outPath_Rad   )//trim(baseName)
            if (ivar==2 .and. ioFiles%outPath_TxTot  /='') outFile = trim(ioFiles%outPath_TxTot )//trim(baseName)
            if (ivar==3 .and. ioFiles%outPath_TxPrfl /='') outFile = trim(ioFiles%outPath_TxPrfl)//trim(baseName)

            if (ivar==1) call writeSpectFile( spect, outfile, 'radiance',      'numPoints', fileHdr)
            if (ivar==2) call writeSpectFile( spect, outfile, 'transmittance', 'numPoints', fileHdr)
            if (ivar==3) call writeSpectFile( spect, outfile, 'transmittance', 'numPoints', fileHdr)
         endif

      ENDDO !DO ivar = 1,3



      IF (inputFile_Jac /='' ) THEN

         inFile = inputFile_Jac
         call splitFilePath( inFile, dirName, baseName)

         !---Find out the Jacobians string
         i1 = index( upper(baseName), 'DRAD-D', back=.TRUE. )
         i2 = index( upper(baseName), '_MONO', back=.TRUE. )
         if (i1<0) STOP '--- '//routineName//'(): Not a valide Jacobians file.'
         if (i2<0) STOP '--- '//routineName//'(): Not a monochromatic Jacobians file.'
         jacStr = baseName(i1+6:i2-1)

         !--- Read data
         !
         dim1Name = 'numPoints'
         dim2Name = 'numLayers'

         select case (upper(jacStr))
         case ('TSKIN')
            varName = 'dRad/dTskin'
            allocate(spectArr(1))
            call readSpectFile_real_1D( inFile, varName, dim1Name, spectArr(1), fileHdr )
         case ('EMIS')
            varName = 'dRad/dEmis'
            allocate(spectArr(1))
            call readSpectFile_real_1D( inFile, varName, dim1Name, spectArr(1), fileHdr )
         case ('RSFC')
            varName = 'dRad/dRsfc'
            allocate(spectArr(1))
            call readSpectFile_real_1D( inFile, varName, dim1Name, spectArr(1), fileHdr )
         case ('T')
            varName = 'dRad/dT'
            allocate(spectArr(1))
            call readSpectFile_real_1D( inFile, varName, dim1Name, spectArr(1), fileHdr )
         case default ! molecular Jac
            varName = 'dRad/d'//upper(trim(jacStr))
            call readSpectFile_real_2D( inFile, varName, [dim1Name,dim2Name], spectArr, fileHdr )
         end select


         !--- Process and output
         !
         if ( postCtrl%functID ==0) then !Filter with instrument SRF

            !Get filter function if filtering is requested.
            call readFilterFunct( trim(ioFiles%filterFunctPath)//&
                                  trim(ioFiles%filterFunctFile), &
                                  filtFunct )

            allocate( filtOut_2D( filtFunct%nChan, size(spectArr) ) )
            call postProc_filter( spectArr, size(spectArr), filtFunct, filtOut_2D)

            fileHdr%productSpectType = 'Filtered'
            fileHdr%V1               = FILLREAL
            fileHdr%V2               = FILLREAL
            fileHdr%DV               = FILLREAL
            fileHdr%nSamp            = FILLINT
            fileHdr%filterFunctFile  = ioFiles%filterFunctFile
            fileHdr%convolParam      = postCtrl

            outFile = inFile
            if ( .not.find_replace_string( outFile, '_mono', '_filtered' )) STOP '--- '//routineName//'(): Not a monochromatic spectral file.'

            if (ioFiles%outPath_Jac /='') then
               call splitFilePath( outFile, dirName, basename)
               outFile = trim(ioFiles%outPath_Jac)//trim(baseName)
            endif

            call writeSpectFile( filtOut_2D, outfile, varName, [dim1Name,dim2Name], fileHdr)

         else !Convolvle with spectrally invariant sanning function

            do il=1,size(spectArr)
               V1O = VLo
               V2O = VHi
               DVO = DV
               padding = 0.
               call postProc_scan( spectArr(il), V1O,V2O,DVO, postCtrl, padding )

               fileHdr%productSpectType = 'Convolved'
               fileHdr%V1               = spectArr(il)%V1
               fileHdr%V2               = spectArr(il)%V2
               fileHdr%DV               = spectArr(il)%DV
               fileHdr%nSamp            = spectArr(il)%NLIM
               !fileHdr%filterFunctFile   = filterFunctFile
               fileHdr%convolParam       = postCtrl
            enddo

            outFile = inFile
            if ( .not.find_replace_string( outFile, '_mono', '_convolved' )) STOP '--- '//routineName//'(): Not a monochromatic spectral file.'

            if (ioFiles%outPath_Jac /='') then
               call splitFilePath( outFile, dirName, basename)
               outFile = trim(ioFiles%outPath_Jac)//trim(baseName)
            endif


            if (size(spectArr)==1) then

               call writeSpectFile( spectArr(1), outfile, varName, dim1Name, fileHdr)

            elseif (size(spectArr)>1) then

               call writeSpectFile( spectArr, outfile, varName, [dim1Name,dim2Name], fileHdr)
            endif

         endif

      ENDIF

   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProcRTdata( spect, V1out,V2out,DVout, &
                              postProc_flags, padding, didPreBox, &
                              filtFunct, filtOut )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum
      USE Module_Config     ,ONLY: CLBLM_Post_Ctrl
      USE Module_PostProc   ,ONLY: postProc_scan,&
                                   postProc_filter, &
                                   CLBLM_FilterFunct
      IMPLICIT NONE

      type(CLBLM_Spectrum)              ,intent(inout) :: spect
      real(r8)                          ,intent(in)    :: V1out,V2out
      real                              ,intent(in)    :: DVout
      type(CLBLM_Post_Ctrl)             ,intent(in)    :: postProc_flags
      real                              ,intent(in)    :: padding
      logical                           ,intent(in)    :: didPreBox !flag to indicate if pre-smoothing was done
      type(CLBLM_FilterFunct) ,optional ,intent(in)    :: filtFunct
      real                    ,optional ,intent(out)   :: filtOut(:)

      character(*), parameter :: routineName='postProcRTdata'

      real(r8)              :: V1O, V2O, DVO
      type(CLBLM_Post_Ctrl) :: postCtrl



      if (postProc_flags%functID /=0) then ! Normal scan

         ! If didProBox, boxcar smoothing was done before in noScattJacob(),
         ! set boxcarHW=-1 and call post processing do normal convolution and
         ! then reset boxcarHW and do deconvolution
         postCtrl = postProc_flags
         if (didPreBox) postCtrl%boxcarHW = -1. !turn off boxcar smoothing and deconvolution in clblm_FFTSCN

         V1O = V1out
         V2O = V2out
         DVO = DVout
         call postProc_scan( spect, V1O,V2O,DVO, postCtrl, padding )

      else !if functID==0, no scan applied

         V1O = spect%V1
         V2O = spect%V2
         DVO = spect%DV
      endif

      !--- Deconvolution
      if ( didPreBox ) then !pre-smoothing was done, do deconvolution.
         postCtrl = postProc_flags

         postCtrl%FFT = .TRUE.
         postCtrl%functID = -1 !Call clblm_FFTSCN to do deconvolution
         DVO = -1. !No interpolation in clblm_FFTSCN after deconvolution

         call postProc_scan( spect, V1O,V2O,DVO, postCtrl, padding )
      endif

      !--- Filtering
      if ( postProc_flags%functID ==0 ) then
         if (.not.present(filtFunct) .or. .not.(present(filtOut))) then
            STOP '--- '//routineName//'(): Filtering requested, but filtFunct or filtOut not present.'
         endif

         call postProc_filter( spect, filtFunct, filtOut )
      endif

   END SUBROUTINE




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeRadiance( fileHdr, radFlag, &
                             radOut, filtOut, &
                             extStr, outFilePath, fullFileName )
!-----------------------------------------------------------------------
      USE netcdf
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum
      USE Module_Config     ,ONLY: ioFiles
      USE Module_FileIO     ,ONLY: createNetCDFFile, &
                                   checkNetCDFcall, &
                                   CLBLM_SpectFileHeader
      USE Module_ConstParam ,ONLY: FILLREAL, FILLINT
      IMPLICIT NONE

      type(CLBLM_SpectFileHeader)    ,intent(inout) :: fileHdr
      integer                        ,intent(in)    :: radFlag
      type(CLBLM_Spectrum) ,OPTIONAL ,intent(in)    :: radOut
      real                 ,OPTIONAL ,intent(in)    :: filtOut(:)
      character(*)         ,OPTIONAL ,intent(in)    :: extStr
      character(*)         ,OPTIONAL ,intent(in)    :: outFilePath
      character(*)         ,OPTIONAL ,intent(in)    :: fullFileName


      character(20)  :: altStr, angstr, altNang, sNNN, extraStr
      character(256) :: outfile, outPath
      integer        :: ncid_out, spectDimId, radoutId



      !fileHdr%sceneFile         = trim(sceneFile)         !filled already
      !fileHdr%fileID            = UID                     !filled already
      !fileHdr%sceneNo           = sceneNo                 !filled already
      !fileHdr%productName                                 !need to be filled here
      !fileHdr%procuctSpectType                            !need to be filled here
      !fileHdr%V1                                          !need to be filled here
      !fileHdr%V2                                          !need to be filled here
      !fileHdr%DV                                          !need to be filled here
      !fileHdr%nSamp                                       !need to be filled here
      !fileHdr%filterFunctFile   = filterFunctFile         !filled already
      !fileHdr%convolParam       = postCtrl                !filled already
      !fileHdr%obsAlt                                      !filled already
      !fileHdr%viewAng                                     !filled already

      altStr=''
      angStr=''
      altNang=''
      if( fileHdr%obsAlt  >=0.) write(altStr, '("_a",I3.3,"km-")') int(fileHdr%obsAlt)
      if( fileHdr%viewAng >=0.) write(angStr, '("o",I3.3,"deg")') int(fileHdr%viewAng)
      altNang = trim(altStr)//trim(angStr)
      write(sNNN, '("s",I3.3)') fileHdr%sceneNo

      extraStr=''
      if (present(extStr)) extraStr='_'//trim(extStr)

      outPath = ioFiles%outPath_Rad
      if (outPath=='') outPath = ioFiles%clblmOutPath
      if (present(outFilePath)) outPath = outFilePath


      if (radFlag ==1) then

         if (present(fullFileName)) then
            outFile = fullFileName
         else
            outFile = trim(outPath)//trim(iofiles%rootName_Rad)//'_mono'//&
                      trim(altNang)//'_'//trim(sNNN)//trim(extraStr)//'.nc'
         endif

         fileHdr%productName = 'Radiance'
         fileHdr%productSpectType = 'Monochromatic'

      elseif (radFlag ==-1) then

         if (fileHdr%convolParam%functID /=0) then

            if (present(fullFileName)) then
               outFile = fullFileName
            else
               outFile = trim(outPath)//trim(iofiles%rootName_Rad)//'_convolved'//&
                         trim(altNang)//'_'//trim(sNNN)//trim(extraStr)//'.nc'
            endif

            fileHdr%productName = 'Radiance'
            fileHdr%productSpectType = 'Convolved'

         else

            if (present(fullFileName)) then
               outFile = fullFileName
            else
               outFile = trim(outPath)//trim(iofiles%rootName_Rad)//'_filtered'//&
                         trim(altNang)//'_'//trim(sNNN)//trim(extraStr)//'.nc'
            endif

            fileHdr%productName = 'Radiance'
            fileHdr%productSpectType = 'Filtered'
         endif

      endif


      if (radFlag>0 .or. fileHdr%convolParam%functID /=0) then

         fileHdr%V1    = radOut%V1
         fileHdr%V2    = radOut%V2
         fileHdr%DV    = radOut%DV
         fileHdr%nSamp = radOut%NLIM

         call writeSpectFile( radout, outfile, 'radiance', 'numPoints', fileHdr)

      else

         fileHdr%V1    = FILLREAL
         fileHdr%V2    = FILLREAL
         fileHdr%DV    = FILLREAL
         fileHdr%nSamp = FILLINT

         call writeSpectFile( filtOut, outfile, 'radiance', 'numPoints', fileHdr)

      endif

   END SUBROUTINE




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeTransmittance( fileHdr, TxFlag, &
                                  TxArray, filtOut,&
                                  outFilePath, fullFileName )
!-----------------------------------------------------------------------
      USE netcdf
      USE Module_ConstParam ,ONLY: FILLREAL, FILLINT
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum
      USE Module_Config     ,ONLY: ioFiles
      USE Module_FileIO     ,ONLY: createNetCDFFile, &
                                   checkNetCDFcall, &
                                   CLBLM_SpectFileHeader
      IMPLICIT NONE

      type(CLBLM_SpectFileHeader)    ,intent(inout) :: fileHdr
      integer                        ,intent(in)    :: TxFlag
      type(CLBLM_Spectrum) ,optional ,intent(in)    :: TxArray(:)
      real                 ,optional ,intent(in)    :: filtOut(:,:)
      character(*)         ,optional ,intent(in)    :: outFilePath
      character(*)         ,optional ,intent(in)    :: fullFileName


      character(80)  :: altStr, angstr, altNang, sNNN, LLL
      character(256) :: outfile,outPath
      integer        :: ncid_out, npDimId, txoutid, nLay, il


      !fileHdr%sceneFile         = trim(sceneFile)         !filled already
      !fileHdr%fileID            = UID                     !filled already
      !fileHdr%sceneNo           = sceneNo                 !filled already
      !fileHdr%productName                                 !need to be filled here
      !fileHdr%procuctSpectType                            !need to be filled here
      !fileHdr%V1                                          !need to be filled here
      !fileHdr%V2                                          !need to be filled here
      !fileHdr%DV                                          !need to be filled here
      !fileHdr%nSamp                                       !need to be filled here
      !fileHdr%filterFunctFile   = filterFunctFile         !filled already
      !fileHdr%convolParam       = postCtrl                !filled already
      !fileHdr%obsAlt                                      !filled already
      !fileHdr%viewAng                                     !filled already


      altStr=''
      angStr=''
      altNang=''
      if( fileHdr%obsAlt  >=0.) write(altStr, '("_a",I3.3,"km-")') int(fileHdr%obsAlt)
      if( fileHdr%viewAng >=0.) write(angStr, '("o",I3.3,"deg")') int(fileHdr%viewAng)
      altNang = trim(altStr)//trim(angStr)
      write(sNNN, '("s",I3.3)') fileHdr%sceneNo


      !  0 = none
      !  1 = mono, total
      ! -1 = scanned, total
      !  2 = mono, profile
      ! -2 = scanned, profile

      IF (abs(TxFlag)==1)  then

         outPath = ioFiles%outPath_TxTot
         if (outPath=='') outPath = ioFiles%clblmOutPath
         if (present(outFilePath)) outPath = outFilePath

         if ( TxFlag ==1) then

            if (present(fullFileName)) then
               outFile = fullFileName
            else
               outFile = trim(outPath)//trim(iofiles%rootName_txtot)//'_mono'//&
                         trim(altNang)//'_'//trim(sNNN)//'.nc'
            endif

            fileHdr%productName = 'Total_transmittance'
            fileHdr%productSpectType = 'Monochromatic'

            fileHdr%V1    = TxArray(1)%V1
            fileHdr%V2    = TxArray(1)%V2
            fileHdr%DV    = TxArray(1)%DV
            fileHdr%nSamp = TxArray(1)%NLIM

         elseif ( TxFlag ==-1) then

            if (fileHdr%convolParam%functID /=0) then

               if (present(fullFileName)) then
                  outFile = fullFileName
               else
                  outFile = trim(outPath)//trim(iofiles%rootName_txtot)//'_convolved'//&
                            trim(altNang)//'_'//trim(sNNN)//'.nc'
               endif

               fileHdr%productName = 'Total_transmittance'
               fileHdr%productSpectType = 'Convolved'

               fileHdr%V1    = TxArray(1)%V1
               fileHdr%V2    = TxArray(1)%V2
               fileHdr%DV    = TxArray(1)%DV
               fileHdr%nSamp = TxArray(1)%NLIM

            else

               if (present(fullFileName)) then
                  outFile = fullFileName
               else
                  outFile = trim(outPath)//trim(iofiles%rootName_txtot)//'_filtered'//&
                            trim(altNang)//'_'//trim(sNNN)//'.nc'
               endif

               fileHdr%productName = 'Total_transmittance'
               fileHdr%productSpectType = 'Filtered'

               fileHdr%V1    = FILLREAL
               fileHdr%V2    = FILLREAL
               fileHdr%DV    = FILLREAL
               fileHdr%nSamp = FILLINT
            endif

         endif

         if (TxFlag>0 .or. fileHdr%convolParam%functID /=0) then
            call writeSpectFile( TxArray(1), outfile, 'transmittance', 'numPoints', fileHdr)
         else
            call writeSpectFile( filtOut(:,1), outfile, 'transmittance', 'numPoints', fileHdr)
         endif

      ELSEIF (abs(TxFlag)==2) then

         outPath = ioFiles%outPath_TxPrfl
         if (outPath=='') outPath = ioFiles%clblmOutPath
         if (present(outFilePath)) outPath = outFilePath

         if (TxFlag>0 .or. fileHdr%convolParam%functID /=0) then
            nLay = size(TxArray)
         else
            nLay = size(filtOut,2)
         endif

         do il = 1, nLay

            write(LLL, '(I3.3)') il

            if ( TxFlag ==2) then

               outFile = trim(outPath)//trim(iofiles%rootName_txprfl)//'_mono'//&
                         trim(altNang)//'_'//trim(sNNN)//'_'//trim(LLL)//'.nc'

               fileHdr%productName = 'Transmittance_profile'
               fileHdr%productSpectType = 'Monochromatic'

               fileHdr%V1    = TxArray(1)%V1
               fileHdr%V2    = TxArray(1)%V2
               fileHdr%DV    = TxArray(1)%DV
               fileHdr%nSamp = TxArray(1)%NLIM

            elseif ( TxFlag ==-2) then

               if (fileHdr%convolParam%functID /=0) then

                  outFile = trim(outPath)//trim(iofiles%rootName_txprfl)//'_convolved'//&
                            trim(altNang)//'_'//trim(sNNN)//'_'//trim(LLL)//'.nc'

                  fileHdr%productName = 'Transmittance_profile'
                  fileHdr%productSpectType = 'Convolved'

                  fileHdr%V1    = TxArray(1)%V1
                  fileHdr%V2    = TxArray(1)%V2
                  fileHdr%DV    = TxArray(1)%DV
                  fileHdr%nSamp = TxArray(1)%NLIM

               else

                  outFile = trim(outPath)//trim(iofiles%rootName_txprfl)//'_filtered'//&
                            trim(altNang)//'_'//trim(sNNN)//'_'//trim(LLL)//'.nc'

                  fileHdr%productName = 'Transmittance_profile'
                  fileHdr%productSpectType = 'Filtered'

                  fileHdr%V1    = FILLREAL
                  fileHdr%V2    = FILLREAL
                  fileHdr%DV    = FILLREAL
                  fileHdr%nSamp = FILLINT
               endif

            endif

            if (TxFlag>0 .or. fileHdr%convolParam%functID /=0) then
               call writeSpectFile( TxArray(il), outfile, 'transmittance', 'numPoints', fileHdr)
            else
               call writeSpectFile( filtOut(:,il), outfile, 'transmittance', 'numPoints', fileHdr)
            endif

         enddo !dol il = 1,nLay

      ENDIF !IF (abs(TxFlag)==1)  then

   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeJacobians( fileHdr, JacMol, JacTemp, JacTskin, JacEmis, JacRsfc, &
                              doJacMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc, &
                              JacFlag, JacList, molList, &
                              filtOut_mol,filtOut_temp,filtOut_Tskin,filtOut_emis,filtOut_Rsfc, &
                              outFilePath)
!-----------------------------------------------------------------------
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum
      USE Module_Config     ,ONLY: ioFiles
      USE Module_FileIO     ,ONLY: CLBLM_SpectFileHeader
      USE Module_Utility    ,ONLY: upper
      USE Module_ConstParam ,ONLY: molIndex, FILLREAL,FILLINT
      USE Module_PostProc   ,ONLY: CLBLM_FilterFunct
      IMPLICIT NONE

      type(CLBLM_SpectFileHeader)       ,intent(inout) :: fileHdr
      type(CLBLM_Spectrum)              ,intent(in)    :: JacMol(:,:)
      type(CLBLM_Spectrum)              ,intent(in)    :: JacTemp(:)
      type(CLBLM_Spectrum)              ,intent(in)    :: JacTskin
      type(CLBLM_Spectrum)              ,intent(in)    :: JacEmis
      type(CLBLM_Spectrum)              ,intent(in)    :: JacRsfc
      logical                           ,intent(in)    :: doJacMol, doJacTemp, doJacTskin, doJacEmis, doJacRsfc
      integer                           ,intent(in)    :: JacFlag
      character(*)                      ,intent(in)    :: JacList(:)
      character(*)                      ,intent(in)    :: molList(:)
      real                    ,optional ,intent(in)    :: filtOut_mol(:,:,:)
      real                    ,optional ,intent(in)    :: filtOut_temp(:,:)
      real                    ,optional ,intent(in)    :: filtOut_Tskin(:)
      real                    ,optional ,intent(in)    :: filtOut_emis(:)
      real                    ,optional ,intent(in)    :: filtOut_Rsfc(:)
      character(*)            ,optional ,intent(in)    :: outFilePath

      character(*) ,parameter :: routineName='writeJacobians'

      integer              :: iJ, il, iMol
      integer              :: nPts,nLay,nMol
      character(10)        :: altStr, angstr, altNang, sNNN
      character(20)        :: JacStr,JacNameStr,varDefName, dim1Nm, dim2Nm
      character(256)       :: outfile, outPath
      integer ,allocatable :: varDimArray(:)
      real    ,allocatable :: tempArr(:,:)
      logical              :: isAirTempJac, isTskinJac, isEmisJac, isRsfcJac, isMolJac


      !fileHdr%sceneFile         = trim(sceneFile)         !filled already
      !fileHdr%fileID            = UID                     !filled already
      !fileHdr%sceneNo           = sceneNo                 !filled already
      !fileHdr%productName                                 !need to be filled here
      !fileHdr%procuctSpectType                            !need to be filled here
      !fileHdr%V1                                          !need to be filled here
      !fileHdr%V2                                          !need to be filled here
      !fileHdr%DV                                          !need to be filled here
      !fileHdr%nSamp                                       !need to be filled here
      !fileHdr%filterFunctFile   = filterFunctFile         !filled already
      !fileHdr%convolParam       = postCtrl                !filled already
      !fileHdr%obsAlt                                      !filled already
      !fileHdr%viewAng                                     !filled already


      altStr=''
      angStr=''
      altNang=''
      if( fileHdr%obsAlt  >=0.) write(altStr, '("_a",I3.3,"km-")') int(fileHdr%obsAlt)
      if( fileHdr%viewAng >=0.) write(angStr, '("o",I3.3,"deg")') int(fileHdr%viewAng)
      altNang = trim(altStr)//trim(angStr)
      write(sNNN, '("s",I3.3)') fileHdr%sceneNo

      !--- Loop over Jacobians, one file per Jacobian
      !
      DO iJ = 1,size(JacList)

         isAirTempJac = .FALSE.
         isTskinJac = .FALSE.
         isEmisJac = .FALSE.
         isRsfcJac = .FALSE.
         isMolJac = .FALSE.

         JacStr = upper(trim(adjustl(JacList(iJ))))
         select case (trim(JacStr))
         case('T')
            if (doJacTemp) then
               nPts = JacTemp(1)%NLIM
               nLay = size( JacTemp )
               nMol = 0
               JacNameStr = 'T'                !"drad-dX" in file name
               varDefName = 'dRad_dT'          !variable name in NetCDF
               if(allocated(varDimArray)) deallocate(varDimArray)
               allocate(varDimArray(2))
               isAirTempJac = .TRUE.
               fileHdr%V1    = JacTemp(1)%V1
               fileHdr%V2    = JacTemp(1)%V2
               fileHdr%DV    = JacTemp(1)%DV
               fileHdr%nSamp = JacTemp(1)%NLIM
            endif
         case('TSKIN')
            if (doJacTskin) then
               nPts = JacTskin%NLIM
               nLay = 0
               nMol = 0
               JacNameStr = 'Tskin'                !"drad-dX" in file name
               varDefName = 'dRad_dTskin'          !variable name in NetCDF
               if(allocated(varDimArray)) deallocate(varDimArray)
               allocate(varDimArray(1))
               isTskinJac = .TRUE.
               fileHdr%V1    = JacTskin%V1
               fileHdr%V2    = JacTskin%V2
               fileHdr%DV    = JacTskin%DV
               fileHdr%nSamp = JacTskin%NLIM
            endif
         case('EMIS')
            if (doJacEmis) then
               nPts = JacEmis%NLIM
               nLay = 0
               nMol = 0
               JacNameStr = 'emis'                !"drad-dX" in file name
               varDefName = 'dRad_dEmis'          !variable name in NetCDF
               if(allocated(varDimArray)) deallocate(varDimArray)
               allocate(varDimArray(1))
               isEmisJac = .TRUE.
               fileHdr%V1    = JacEmis%V1
               fileHdr%V2    = JacEmis%V2
               fileHdr%DV    = JacEmis%DV
               fileHdr%nSamp = JacEmis%NLIM
            endif
         case('RSFC')
            if (doJacRsfc) then
               nPts = JacRsfc%NLIM
               nLay = 0
               nMol = 0
               JacNameStr = 'Rsfc'                !"drad-dX" in file name
               varDefName = 'dRad_dRsfc'          !variable name in NetCDF
               if(allocated(varDimArray)) deallocate(varDimArray)
               allocate(varDimArray(1))
               isRsfcJac = .TRUE.
               fileHdr%V1    = JacRsfc%V1
               fileHdr%V2    = JacRsfc%V2
               fileHdr%DV    = JacRsfc%DV
               fileHdr%nSamp = JacRsfc%NLIM
            endif
         case default
            if ( doJacMol ) then
               iMol = molIndex(JacStr,molList)
               if (iMol>0) then !is one of the molecules in the scene file.
                  nPts = JacMol(1,iMol)%NLIM
                  nLay = size( JacMol, 1 )
                  nMol = 0 !Only one molecule is written each time, no need dimension for molecules.
                  JacNameStr = JacStr                        !"drad-dX" in file name
                  varDefName = 'dRad_d'//upper(trim(JacStr)) !variable name in NetCDF
                  if(allocated(varDimArray)) deallocate(varDimArray)
                  allocate(varDimArray(2))
                  isMolJac = .TRUE.
                  fileHdr%V1    = JacMol(1,iMOl)%V1
                  fileHdr%V2    = JacMol(1,iMOl)%V2
                  fileHdr%DV    = JacMol(1,iMOl)%DV
                  fileHdr%nSamp = JacMol(1,iMOl)%NLIM
               endif
            else
               CYCLE !Nother to write, try next one
            endif
         end select

         outPath = ioFiles%outPath_Jac
         if (outPath=='') outPath = ioFiles%clblmOutPath
         if (present(outFilePath)) outPath = outFilePath

         if (JacFlag ==1) then

            outFile = trim(outPath)//trim(iofiles%rootName_jac)//trim(JacNameStr)//&
                      '_mono'//trim(altNang)//'_'//trim(sNNN)//'.nc'

            fileHdr%productName = 'Concentration_Jacobians'
            fileHdr%productSpectType = 'Monochromatic'

         elseif (JacFlag ==-1) then

            if (fileHdr%convolParam%functID /=0) then

               outFile = trim(outPath)//trim(iofiles%rootName_jac)//trim(JacNameStr)//&
                         '_convolved'//trim(altNang)//'_'//trim(sNNN)//'.nc'

               fileHdr%productName = 'Concentraton_Jacobians'
               fileHdr%productSpectType = 'Convolved'

            else

               outFile = trim(outPath)//trim(iofiles%rootName_jac)//trim(JacNameStr)//&
                         '_filtered'//trim(altNang)//'_'//trim(sNNN)//'.nc'

               fileHdr%productName = 'Concentraton_Jacobians'
               fileHdr%productSpectType = 'Filtered'

               fileHdr%V1    = FILLREAL
               fileHdr%V2    = FILLREAL
               fileHdr%DV    = FILLREAL
               fileHdr%nSamp = FILLINT

            endif

         endif


         !--- Write out data
         !
         dim1Nm = 'numPoints'
         dim2Nm = 'numLayers'

         if (JacFlag ==-1 .and. fileHdr%convolParam%functID ==0) then !Filtered results

            if (    isMolJac)     then; call writeSpectFile( filtOut_mol(:,:,iMol), trim(outfile), trim(varDefName), [dim1Nm,dim2Nm], fileHdr)
            elseif (isAirTempJac) then; call writeSpectFile( filtOut_temp,          trim(outfile), trim(varDefName), [dim1Nm,dim2Nm], fileHdr)
            elseif (isTskinJac)   then; call writeSpectFile( filtOut_Tskin,         trim(outfile), trim(varDefName), dim1Nm, fileHdr)
            elseif (isEmisJac)    then; call writeSpectFile( filtOut_emis,          trim(outfile), trim(varDefName), dim1Nm, fileHdr)
            elseif (isRsfcJac)    then; call writeSpectFile( filtOut_Rsfc,          trim(outfile), trim(varDefName), dim1Nm, fileHdr)
            endif

         else !mono or scanned results

            if (isMolJac) then

               if (allocated(tempArr)) deallocate(tempArr)
               allocate( tempArr(nPts,nLay) )
               do il=1,nLay
                  tempArr(1:nPts,il) = JacMol(il,iMol)%spect(1:nPts)
               enddo

               call writeSpectFile( tempArr, trim(outfile), trim(varDefName), [dim1Nm,dim2Nm], fileHdr)

            elseif (isAirTempJac) then
               call writeSpectFile( JacTemp, trim(outfile), trim(varDefName), [dim1Nm,dim2Nm], fileHdr)
            elseif (isTskinJac) then
               call writeSpectFile( JacTskin, trim(outfile), trim(varDefName), dim1Nm, fileHdr)
            elseif (isEmisJac) then
               call writeSpectFile( JacEmis, trim(outfile), trim(varDefName), dim1Nm, fileHdr)
            elseif (isRsfcJac) then
               call writeSpectFile( JacRsfc, trim(outfile), trim(varDefName), dim1Nm, fileHdr)
            endif

         endif !(JacFlag ==-1 .and. fileHdr%convolParam%functID ==0)

      ENDDO !DO iJ = 1,size(JacList)

   END SUBROUTINE



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeOD( fileHdr, monoOD, layNo, nlteEmisFac, outFilePath, fullFileName )
!-----------------------------------------------------------------------
      USE netcdf
      !USE netcdfbase       ,only: createNetCDFFile, checkNetCDFcall
      !USE netcdfwriter     ,only: writeSpectFileHeader
      USE Module_Config    ,ONLY: ioFiles
      USE Module_Spectrum  ,ONLY: CLBLM_Spectrum
      USE Module_FileIO ,ONLY: createNetCDFFile, &
                                  checkNetCDFcall, &
                                  CLBLM_SpectFileHeader
      IMPLICIT NONE

      type(CLBLM_SpectFileHeader)    ,intent(inout) :: fileHdr
      type(CLBLM_Spectrum)           ,intent(in)    :: monoOD
      integer                        ,intent(in)    :: layNo
      type(CLBLM_Spectrum) ,optional ,intent(in)    :: nlteEmisFac
      character(*)         ,optional ,intent(in)    :: outFilePath
      character(*)         ,optional ,intent(in)    :: fullFileName


      character(20)  :: altStr, angstr, altNang, sNNN, LLL
      character(256) :: outfile, outPath
      integer(4)     :: ncid_out, spectDimId, odId, nlteId



      !fileHdr%sceneFile         = trim(sceneFile)         !filled already
      !fileHdr%fileID            = UID                     !filled already
      !fileHdr%sceneNo           = sceneNo                 !filled already
      !fileHdr%productName                                 !need to be filled here
      !fileHdr%procuctSpectType                            !need to be filled here
      !fileHdr%V1                                          !need to be filled here
      !fileHdr%V2                                          !need to be filled here
      !fileHdr%DV                                          !need to be filled here
      !fileHdr%nSamp                                       !need to be filled here
      !fileHdr%filterFunctFile   = filterFunctFile         !filled already
      !fileHdr%convolParam       = postCtrl                !filled already
      !fileHdr%obsAlt                                      !filled already
      !fileHdr%viewAng                                     !filled already

      altStr=''
      angStr=''
      altNang=''
      if( fileHdr%obsAlt  >=0.) write(altStr, '("_a",I3.3,"km-")') int(fileHdr%obsAlt)
      if( fileHdr%viewAng >=0.) write(angStr, '("o",I3.3,"deg")') int(fileHdr%viewAng)
      altNang = trim(altStr)//trim(angStr)
      write(sNNN, '("s",I3.3)') fileHdr%sceneNo
      write(LLL, '(I3.3)') layNo


      !--- Output file name
      outPath = ioFiles%outPath_OD
      if (outPath=='') outPath = ioFiles%clblmOutPath
      if (present(outFilePath)) outPath = outFilePath

      if (present(fullFileName)) then
         outFile = fullFileName
      else
         outFile = trim(outPath)//trim(iofiles%rootName_OD)//&
                   trim(altNang)//'_'//trim(sNNN)//'_'//trim(LLL)//'.nc'
      endif

      fileHdr%productName = 'Layer_OD'
      fileHdr%productSpectType = 'Monochromatic'

      fileHdr%V1     = monoOD%V1
      fileHdr%V2     = monoOD%V2
      fileHdr%DV     = monoOD%DV
      fileHdr%nSamp  = monoOD%NLIM


      !--- Create output file
      if (.not. createNetCDFFile(trim(outfile), ncid_out)) then
          STOP "cannot create netcdf file for scene, this should not happen"
      endif
      call checkNetCDFcall( nf90_enddef(ncid_out) )


      !--- Write file header and radiance data
      call writeSpectFileHeader(ncID_out, fileHdr)

      call checkNetCDFcall(nf90_def_dim(ncid_out, "numPoints", int(size(monoOD%spect),4), spectDimId))
      call checkNetCDFcall(nf90_def_var(ncid_out, "layerOD", NF90_DOUBLE, (/ spectDimId /), odId))
      if (present(nlteEmisFac)) then
         call checkNetCDFcall(nf90_def_var(ncid_out, "nlteEmisFac", NF90_DOUBLE, (/ spectDimId /), nlteId))
      endif


      ! write the data
      call checkNetCDFcall(nf90_enddef(ncid_out))

      call checkNEtCDFcall(nf90_put_var(ncid_out, odId, monoOD%spect))
      if (present(nlteEmisFac)) then
         call checkNEtCDFcall(nf90_put_var(ncid_out, nlteId, nlteEmisFac%spect))
      endif

      call checkNetCDFcall(nf90_close(ncid_out))

   END SUBROUTINE




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeSpectFile_1D_spect( spect, fileName, varName, dimName, fileHdr)
!-----------------------------------------------------------------------
      USE Module_Spectrum  ,ONLY: CLBLM_Spectrum
      USE Module_FileIO    ,ONLY: CLBLM_SpectFileHeader
      IMPLICIT NONE

      type(CLBLM_Spectrum)           ,intent(in) :: spect
      character(*)                   ,intent(in) :: fileName
      character(*)                   ,intent(in) :: varName
      character(*)                   ,intent(in) :: dimName
      type(CLBLM_SpectFileHeader)    ,intent(in) :: fileHdr


      call writeSpectFile_1D_array( spect%spect( 1:spect%NLIM ), &
                                    fileName, varName, dimName, fileHdr)

   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeSpectFile_1D_array( dataOut, fileName, varName, dimName, fileHdr)
!-----------------------------------------------------------------------
      USE netcdf
      USE Module_Spectrum  ,ONLY: CLBLM_Spectrum
      USE Module_FileIO    ,ONLY: createNetCDFFile, &
                                  checkNetCDFcall, &
                                  CLBLM_SpectFileHeader
      IMPLICIT NONE

      character(*) ,parameter :: routineName='writeSpectFile_1D_array'

      real                           ,intent(in) :: dataOut(:) ![nPoints]
      character(*)                   ,intent(in) :: fileName
      character(*)                   ,intent(in) :: varName
      character(*)                   ,intent(in) :: dimName
      type(CLBLM_SpectFileHeader)    ,intent(in) :: fileHdr

      integer(4) :: ncid, dimID, varID


      !--- Create output file
      if (.not. createNetCDFFile(trim(fileName), ncid)) then
         STOP '--- '//routineName//'(): cannot create netcdf file, this should not happen'
      endif
      call checkNetCDFcall( nf90_enddef(ncid) )

      !--- Write file header and radiance data
      call writeSpectFileHeader(ncid, fileHdr)

      call checkNetCDFcall(nf90_def_dim( ncid, dimName, int(size(dataOut),4), dimID))
      call checkNetCDFcall(nf90_def_var( ncid, varName, NF90_DOUBLE, [dimID], varID))

      ! write the data
      call checkNetCDFcall(nf90_enddef( ncid))
      call checkNEtCDFcall(nf90_put_var( ncid, varID, dataOut ))
      call checkNetCDFcall(nf90_close( ncid))

   END SUBROUTINE



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeSpectFile_2D_spect( spect, fileName, varName, dimName, fileHdr)
!-----------------------------------------------------------------------
      USE Module_Spectrum  ,ONLY: CLBLM_Spectrum
      USE Module_FileIO    ,ONLY: CLBLM_SpectFileHeader
      IMPLICIT NONE

      type(CLBLM_Spectrum)           ,intent(in) :: spect(:)
      character(*)                   ,intent(in) :: fileName
      character(*)                   ,intent(in) :: varName
      character(*)                   ,intent(in) :: dimName(:)
      type(CLBLM_SpectFileHeader)    ,intent(in) :: fileHdr

      integer :: dim1,dim2, il
      real ,allocatable :: tempArr(:,:)


      dim1 = spect(1)%NLIM
      dim2 = size(spect)
      allocate(tempArr( dim1,dim2 ))
      do il=1,dim2
         tempArr(1:dim1,il) = spect(il)%spect(1:dim1)
      enddo

      call writeSpectFile_2D_array( tempArr, fileName, varName, dimName, fileHdr)

   END SUBROUTINE

!-----------------------------------------------------------------------
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
      type(CLBLM_SpectFileHeader)    ,intent(in) :: fileHdr

      integer(4) :: ncid, dim1ID,dim2ID, dim1,dim2, varID


      !--- Create output file
      if (.not. createNetCDFFile(trim(fileName), ncid)) then
         STOP '--- '//routineName//'(): cannot create netcdf file, this should not happen'
      endif
      call checkNetCDFcall( nf90_enddef(ncid) )

      !--- Write file header and radiance data
      call writeSpectFileHeader(ncid, fileHdr)

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



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE writeSpectFile_2D_flux( dataOut,  flux_bound, data2Out,data3Out,data4Out,data5Out,data6Out, varName, dimName, fileHdr)
      !-----------------------------------------------------------------------
            USE netcdf
            USE Module_Config    ,ONLY: ioFiles
            USE Module_Spectrum  ,ONLY: CLBLM_Spectrum
            USE Module_FileIO    ,ONLY: createNetCDFFile, &
                                        checkNetCDFcall, &
                                        CLBLM_SpectFileHeader
            IMPLICIT NONE

            character(*) ,parameter :: routineName='writeSpectFile_2D_array'

            real                           ,intent(in) :: dataOut(:,:) ![nlevels, bins] upwelling fluxes
            real                           ,intent(in) :: data2Out(:,:) ![nlevels, bins] downwelling fluxes
            real                           ,intent(in) :: data3Out(:) ![nlevels] pressure
            integer                        ,intent(in) :: data4Out(:) ![nlevels] levels
            real                           ,intent(in) :: data5Out(:,:) ![nlevels,bins] net flux
            real                           ,intent(in) :: data6Out(:,:) ![nlevels,bins] heating rate
            real                           ,intent(in) :: flux_bound(:) ![bins]
            character(*)                   ,intent(in) :: varName(:)
            character(*)                   ,intent(in) :: dimName(:)
            type(CLBLM_SpectFileHeader)    ,intent(in) :: fileHdr

            integer(4) :: ncid, dim1ID,dim2ID, dim1,dim2, varID, varID2, varID3, varID4, varID5, varID6, varID7
            character(256) :: outfile, outPath

            outPath = ioFiles%outPath_Rad
            outFile = trim(outPath)//trim(iofiles%rootName_Rad)//'_FLUX_OUTPUT'//'.nc'

            !--- Create output file
            if (.not. createNetCDFFile(trim(outFile), ncid)) then
               STOP '--- '//routineName//'(): cannot create netcdf file, this should not happen'
            endif
            call checkNetCDFcall( nf90_enddef(ncid) )

            !--- Write file header and radiance data
            call writeSpectFileHeader(ncid, fileHdr)

            dim1 = size(dataOut,1)
            dim2 = size(dataOut,2)
            print *,'dim1', dim1
            print *,'dim2', dim2
            call checkNetCDFcall(nf90_def_dim( ncid, dimName(1), dim1, dim1ID))
            call checkNetCDFcall(nf90_def_dim( ncid, dimName(2), dim2, dim2ID))
            call checkNetCDFcall(nf90_def_var( ncid, varName(1), NF90_DOUBLE, [dim1ID,dim2ID], varID))
            call checkNetCDFcall(nf90_def_var( ncid, varName(2), NF90_DOUBLE, [dim2ID], varID2))
            call checkNetCDFcall(nf90_def_var( ncid, varName(3), NF90_DOUBLE, [dim1ID,dim2ID], varID3))
            call checkNetCDFcall(nf90_def_var( ncid, varName(4), NF90_DOUBLE, [dim1ID], varID4))
            call checkNetCDFcall(nf90_def_var( ncid, varName(5), NF90_DOUBLE, [dim1ID], varID5))
            call checkNetCDFcall(nf90_def_var( ncid, varName(6), NF90_DOUBLE, [dim1ID,dim2ID], varID6))
            call checkNetCDFcall(nf90_def_var( ncid, varName(7), NF90_DOUBLE, [dim1ID,dim2ID], varID7))

            call checkNetCDFcall(nf90_enddef( ncid))

            ! write the data
            call checkNEtCDFcall(nf90_put_var( ncid, varID, dataOut(:,:) ))
            call checkNEtCDFcall(nf90_put_var( ncid, varID2, flux_bound(:) ))
            call checkNEtCDFcall(nf90_put_var( ncid, varID3, data2Out(:,:) ))
            call checkNEtCDFcall(nf90_put_var( ncid, varID4, data3Out(:) ))
            call checkNEtCDFcall(nf90_put_var( ncid, varID5, data4Out(:) ))
            call checkNEtCDFcall(nf90_put_var( ncid, varID6, data5Out(:,:) ))
            call checkNEtCDFcall(nf90_put_var( ncid, varID7, data6Out(:,:) ))

            call checkNetCDFcall(nf90_close( ncid))

         END SUBROUTINE



!-----------------------------------------------------------------------
   SUBROUTINE readSpectFile_real_1D( fileName, varName, dimName, spect, fileHdr )
!-----------------------------------------------------------------------
      USE netcdf
      USE Module_FileIO     ,ONLY: openNetCDFFile, &
                                   check, dimByName, &
                                   CLBLM_SpectFileHeader,&
                                   readArray_real_1D
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum, &
                                   CLBLM_Spectrum_init
      IMPLICIT NONE

      character(*)                ,intent(in)  :: fileName
      character(*)                ,intent(in)  :: varName
      character(*)                ,intent(in)  :: dimName
      type(CLBLM_Spectrum)        ,intent(out) :: spect
      type(CLBLM_SpectFileHeader) ,intent(out) :: fileHdr

      integer(4) :: ncid, dim1


      call openNetCDFFile( trim(fileName), ncID )

      call readSpectFileHeader( ncID, fileHdr )

      dim1 = dimByName( ncid, dimName)

      call CLBLM_Spectrum_init( spect, fileHdr%V1,fileHdr%DV,fileHdr%nSamp )

      call readArray_real_1D( ncid, trim(varName), spect%spect )

   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE readSpectFile_real_2D( fileName, varName, dimName, spect, fileHdr )
!-----------------------------------------------------------------------
      USE netcdf
      USE Module_FileIO     ,ONLY: openNetCDFFile, &
                                   check, dimByName, &
                                   CLBLM_SpectFileHeader, &
                                   readArray_real_2D
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum, &
                                   CLBLM_Spectrum_init
      IMPLICIT NONE

      character(*)                      ,intent(in)  :: fileName
      character(*)                      ,intent(in)  :: varName
      character(*)                      ,intent(in)  :: dimName(2)
      type(CLBLM_Spectrum) ,allocatable ,intent(out) :: spect(:)
      type(CLBLM_SpectFileHeader)       ,intent(out) :: fileHdr

      integer(4) :: ncid, ndim, dim1,dim2, il
      real, allocatable :: tempArr(:,:)

      !ndim = size(dimName)

      call openNetCDFFile( trim(fileName), ncID )

      call readSpectFileHeader( ncID, fileHdr )

      dim1 = dimByName( ncid, trim(dimName(1)))
      dim2 = dimByName( ncid, trim(dimName(2)))

      allocate(tempArr( dim1,dim2 ))
      call readArray_real_2D( ncid, trim(varName), tempArr )


      !--- Load the output spectrum array
      !
      if (allocated(spect)) deallocate(spect)
      allocate(spect(dim2))
      do il =1,dim2
         call CLBLM_Spectrum_init( spect(il), fileHdr%V1,fileHdr%DV,fileHdr%nSamp )
         spect(il)%spect(1:dim1) = tempArr(1:dim1,il)
      enddo

   END SUBROUTINE





   ! -------------------------------------------------------------------------
   ! Writes a variety of parameters as global attributes to a NetCDF output file.
   !--------------------------------------------------------------------------
   SUBROUTINE writeSpectFileHeader( ncid, fileHdr )
   !--------------------------------------------------------------------------
      USE NETCDF
      USE Module_FileIO    ,ONLY: CLBLM_SpectFileHeader, check
      USE Module_Utility   ,ONLY: upper
      IMPLICIT NONE

      integer*4                   ,intent(in) :: ncid
      type(CLBLM_SpectFileHeader) ,intent(in) :: fileHdr

      integer :: fftFlag, deconvFlag

      !fileHdr%sceneFile
      !fileHdr%fileID
      !fileHdr%sceneNo
      !fileHdr%productName
      !fileHdr%procuctSpectType
      !fileHdr%V1
      !fileHdr%V2
      !fileHdr%DV
      !fileHdr%nSamp
      !fileHdr%filterFunctFile
      !fileHdr%convolParam
      !fileHdr%obsAlt
      !fileHdr%viewAng

      ! Put file pointer back into define mode
      call check( nf90_redef(ncid))

      call check( nf90_put_att( ncid, NF90_GLOBAL, "title",            "CLBLM Output File"))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "sceneFileName",    fileHdr%sceneFile))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "UID",              fileHdr%fileID))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "sceneNumber",      fileHdr%sceneNo))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "productName",      fileHdr%productName))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "spectralDataType", fileHdr%productSpectType))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "V1",               fileHdr%V1))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "V2",               fileHdr%V2))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "DV",               fileHdr%DV))
      call check( nf90_put_att( ncid, NF90_GLOBAL, "NSAMP",            fileHdr%nSamp))

      if ( upper(trim(fileHdr%productSpectType)) .eq. "CONVOLVED") then
         call check( nf90_put_att( ncid, NF90_GLOBAL, "functionID",    fileHdr%convolParam%functID))
         if ( fileHdr%convolParam%functID ==0) then
            call check( nf90_put_att( ncid, NF90_GLOBAL, "convolutionFunction", trim(fileHdr%filterFunctFile)))
         else
            if (fileHdr%convolParam%FFT) then
               fftFlag = 1
            else
               fftFlag = 0
            endif

            if (fileHdr%convolParam%deconvPreScan) then
               deconvFlag = 1
            else
               deconvFlag = 0
            endif

            call check( nf90_put_att( ncid, NF90_GLOBAL, "FFT_flag",   fftFlag))
            call check( nf90_put_att( ncid, NF90_GLOBAL, "functHWHM",  fileHdr%convolParam%HWHM))
            call check( nf90_put_att( ncid, NF90_GLOBAL, "functParam", fileHdr%convolParam%functParams(1:3) ))
            call check( nf90_put_att( ncid, NF90_GLOBAL, "boxcarHW",   fileHdr%convolParam%boxcarHW))
            call check( nf90_put_att( ncid, NF90_GLOBAL, "deconvFlag", deconvFlag))
         endif
      endif

      return
    END SUBROUTINE

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
   SUBROUTINE readSpectFileHeader( ncid, fileHdr )
!--------------------------------------------------------------------------
      USE NETCDF
      USE Module_FileIO    ,ONLY: CLBLM_SpectFileHeader, check
      USE Module_Utility   ,ONLY: upper
      IMPLICIT NONE

      integer*4                   ,intent(in)  :: ncid
      type(CLBLM_SpectFileHeader) ,intent(out) :: fileHdr

      character(128) :: title
      integer        :: fftFlag, deconvFlag


      !fileHdr%sceneFile
      !fileHdr%fileID
      !fileHdr%sceneNo
      !fileHdr%productName
      !fileHdr%procuctSpectType
      !fileHdr%V1
      !fileHdr%V2
      !fileHdr%DV
      !fileHdr%nSamp
      !fileHdr%filterFunctFile
      !fileHdr%convolParam
      !fileHdr%obsAlt
      !fileHdr%viewAng

      call check( nf90_get_att( ncid, NF90_GLOBAL, "title",            title))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "sceneFileName",    fileHdr%sceneFile))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "UID",              fileHdr%fileID))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "sceneNumber",      fileHdr%sceneNo))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "productName",      fileHdr%productName))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "spectralDataType", fileHdr%productSpectType))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "V1",               fileHdr%V1))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "V2",               fileHdr%V2))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "DV",               fileHdr%DV))
      call check( nf90_get_att( ncid, NF90_GLOBAL, "NSAMP",            fileHdr%nSamp))

      if ( upper(trim(fileHdr%productSpectType)) .eq. "CONVOLVED") then
         call check( nf90_get_att( ncid, NF90_GLOBAL, "functionID",    fileHdr%convolParam%functID))
         if ( fileHdr%convolParam%functID ==0) then
            call check( nf90_get_att( ncid, NF90_GLOBAL, "convolutionFunction", fileHdr%filterFunctFile))
         else
            call check( nf90_get_att( ncid, NF90_GLOBAL, "FFT_flag",   fftFlag))
            call check( nf90_get_att( ncid, NF90_GLOBAL, "functHWHM",  fileHdr%convolParam%HWHM))
            call check( nf90_get_att( ncid, NF90_GLOBAL, "functParam", fileHdr%convolParam%functParams(1:3) ))
            call check( nf90_get_att( ncid, NF90_GLOBAL, "boxcarHW",   fileHdr%convolParam%boxcarHW))
            call check( nf90_get_att( ncid, NF90_GLOBAL, "deconvFlag", deconvFlag))

            if (fftFlag>0) then
               fileHdr%convolParam%FFT = .TRUE.
            else
               fileHdr%convolParam%FFT = .FALSE.
            endif

            if (deconvFlag>0) then
               fileHdr%convolParam%deconvPreScan = .TRUE.
            else
               fileHdr%convolParam%deconvPreScan = .FALSE.
            endif
         endif
      endif

      return
   END SUBROUTINE

END MODULE
