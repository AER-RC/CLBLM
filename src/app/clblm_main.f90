
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Program clblm_main
   !-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: FILLREAL
      USE Module_Scene        ,ONLY: CLBLM_Scene, &
                                     readSceneFileHeader, &
                                     inputSceneData,&
                                     check_sceneData
      USE Module_Config       ,ONLY: CLBLM_Output_Ctrl, &
                                     CLBLM_OutputSpectGrid, &
                                     CLBLM_Path_Ctrl, &
                                     CLBLM_OD_Ctrl, &
                                     CLBLM_DV_Ctrl, &
                                     CLBLM_RT_Ctrl, &
                                     CLBLM_Flux_Ctrl, &
                                     CLBLM_Post_Ctrl, &
                                     SceneFileOverRide, &
                                     openLogFile, &
                                     NLTE, &
                                     ioFiles, noPrnt, inputUserDirectives
      USE Module_Spectrum     ,ONLY: CLBLM_Spectrum, &
                                     CLBLM_Spectrum_init
      USE Module_Drivers      ,ONLY: RT_basic, RT_Jac, &
                                     convolPreviousResults
      USE Module_FileIO       ,ONLY: CLBLM_SpectFileHeader, &
                                     openNetCDFFile
      IMPLICIT NONE

      integer                           :: isn,iz,ia
      logical                           :: convolOnly
      character(256)                    :: configFile
      integer                           :: nSelectScenes, nObsAlt, nObsAng
      integer(4)                        :: ncID
      integer                           :: nscenes, UID
      integer              ,allocatable :: sceneNums(:)
      integer                           :: sceneNo
      type(CLBLM_Scene)                 :: scene
      type(SceneFileOverRide)           :: sceneOverride
      real                              :: ovrdObsAlt, ovrdViewAng
      real                 ,allocatable :: pRT(:)
      integer                           :: pRTsize
      real                 ,allocatable :: zRT(:)
      integer                           :: zRTsize
      character(20)        ,allocatable :: molID(:)
      integer                           :: numMol
      type(CLBLM_Path_Ctrl)             :: pathCtrl
      type(CLBLM_DV_Ctrl)               :: dvCtrl
      type(CLBLM_OD_Ctrl)               :: odCtrl
      type(CLBLM_RT_Ctrl)               :: rtCtrl
      type(CLBLM_Flux_Ctrl)             :: fluxCtrl
      type(CLBLM_Output_Ctrl)           :: outctrl
      type(CLBLM_OutputSpectGrid)       :: outGrid
      type(CLBLM_Post_Ctrl)             :: postCtrl
      integer                           :: gridType
      real                 ,allocatable :: TxOut(:,:)
      real                 ,allocatable :: RadOut(:)
      real                 ,allocatable :: flxuOut(:,:)
      real                 ,allocatable :: flxdOut(:,:)
      real                 ,allocatable :: JacOut_mol(:,:,:)  ![NLIM,nLay,nMol], A.J. w.r.t molecular contentration
      real                 ,allocatable :: JacOut_temp(:,:)   ![NLIM,nLay],      A.J. w.r.t air temperature
      real                 ,allocatable :: JacOut_Tskin(:)    ![NLIM],           A.J. w.r.t surface temperature
      real                 ,allocatable :: JacOut_emis(:)     ![NLIM],           A.J. w.r.t surface emissivity
      real                 ,allocatable :: JacOut_Rsfc(:)     ![NLIM],           A.J. w.r.t surface reflectance
      real                 ,allocatable :: ODout(:,:)
      real                 ,allocatable :: nlteEmisFacOut(:,:)
      integer                           :: nSampRad, nSampTx, nSampOD, nSampJac
      type(CLBLM_SpectFileHeader)       :: fileHdr
      logical                           :: isExist



      !--- Get user's directive file name from command line arguments
      ! If not specified, use the default config file 'clblm_config.json'
      CALL get_command_argument( number=1, value=configFile )
      if (configFile=='') then
         configFile='clblm_config.json'
         inquire(file=configFile, exist=isExist)
         if (.not. isExist) then
           print *,'not found default config file: ', trim(configFile)
           call exit(1)
         endif
      else
         inquire(file=configFile, exist=isExist)
         if (.not. isExist) then
           print *,'not found user config file: ', trim(configFile)
           call exit(1)
         endif
      endif

      CALL openLogFile()

      CALL inputUserDirectives( trim(configFile), &
                               outCtrl, outGrid, pathCtrl, dvCtrl, &
                               odCtrl, rtCtrl, postCtrl, sceneOverride, FluxCtrl )

      IF ( outCtrl%OD/=0 .or. &
           outCtrl%Tx/=0 .or. &
           outCtrl%Rad/=0 .or. &
           outCtrl%Jac/=0) THEN

         nObsAlt = 0
         nObsAng = 0
         if (allocated( sceneOverride%obsAlt  )) nObsAlt = size(sceneOverride%obsAlt)
         if (allocated( sceneOverride%viewAng )) nObsAng = size(sceneOverride%viewAng)


         !---Read scene data and call RT
         call openNetCDFFile( trim(ioFiles%sceneFile), ncID )
         call readSceneFileHeader( ncID, UID, nScenes,sceneNums, &
                                   pRT,zRT,pRTsize,zRTsize, molID, numMol)


         do isn = 1,nScenes

            sceneNo = sceneNums(isn)

            !--- Check if this scene is a selected one
            if ( allocated( sceneOverride % selectedScene_ID )) then
               if (.not. any( sceneOverride % selectedScene_ID - sceneNo ==0)) CYCLE
            elseif ( sceneOverride % numSelectedScenes >0) then
               if (sceneNo > sceneOverride % numSelectedScenes) CYCLE
            endif

            call inputSceneData( ncID, sceneNo, scene)

            !--- If user input alternative viewing geometry, use these geometry instead.
            ! Otherwise use viewing geometry come with the scene data (index==0).
            !
            do iz = 0, nObsAlt

               !--- If there are alternative viewing geometry exist, skip the
               ! viewing geomtry from scene file.
               if (iz==0.and.nObsAlt>0) CYCLE
               if (iz/=0) then
                  ovrdObsAlt = sceneOverride%obsAlt(iz)
               else
                  ovrdObsAlt = FILLREAL
               endif


               do ia = 0, nObsAng

                  !--- If there are alternative viewing geometry exist, skip the
                  ! viewing geomtry from scene file.
                  if (ia==0.and.nObsAng>0) CYCLE
                  if (ia/=0) then
                     ovrdViewAng = sceneOverride%viewAng(ia)
                  else
                     ovrdViewAng = FILLREAL
                  endif


                  !--- Check scene data
                  call check_sceneData( scene, &
                                        molID, numMol, &
                                        pRT, pRTsize, zRT,zRTsize, pathCtrl%RTgridFlag, &
                                        ovrdObsAlt, ovrdViewAng, outCtrl )


                  !Fill out file header structure used for output files
                  fileHdr%sceneFile         = trim(ioFiles%sceneFile)
                  fileHdr%fileID            = UID
                  fileHdr%sceneNo           = sceneNo
                  !fileHdr%productName                        !to be filled later
                  !fileHdr%procuctSpectType                   !to be filled later
                  !fileHdr%V1                                 !to be filled later
                  !fileHdr%V2                                 !to be filled later
                  !fileHdr%DV                                 !to be filled later
                  !fileHdr%nSamp                              !to be filled later !number of spectral points
                  fileHdr%filterFunctFile   = trim(ioFiles%filterFunctFile)
                  fileHdr%convolParam       = postCtrl
                  fileHdr%obsAlt            = scene%geom%obs%obsAlt
                  fileHdr%viewAng           = scene%geom%obs%viewAng


                  if (outCtrl%Jac/=0) then

                     CALL RT_Jac( outCtrl%Rad, outCtrl%Jac, & !Jacobians_list, scene, &
                                  outctrl%jaclist, scene, &
                                  pathCtrl, dvCtrl, odCtrl, rtCtrl, fluxCtrl,&
                                  outGrid%V1, outGrid%V2, outGrid%DV, &
                                  JacOut_mol, JacOut_temp, JacOut_Tskin, &
                                  JacOut_emis, JacOut_Rsfc, RadOut, nSampJac, &
                                  postCtrl, zRT,zRTsize, fileHdr=fileHdr )

                  elseif ( outCtrl%Rad/=0 .or. outCtrl%Tx/=0 .or. outCtrl%OD/=0 ) then

                     CALL RT_basic( outCtrl%Rad, outCtrl%Tx, outCtrl%OD, scene, &
                                    pathCtrl, dvCtrl, odCtrl, rtCtrl, fluxCtrl, &
                                    outGrid%V1, outGrid%V2, outGrid%DV, flxuOut, flxdOut, &
                                    RadOut,nSampRad, TxOut,nSampTx, postCtrl, &
                                    zRT,zRTsize, NLTE, grid_type=outGrid%gridType, fileHdr=fileHdr )

                  endif

               enddo !do ia = 0, nObsAng
            enddo !do iz = 0, nObsAlt
         enddo !do isn = 1,nScenes
      ENDIF


      !Do convolution only processing
      IF ( ioFiles%inFile_Rad    /='' .or. &
           ioFiles%inFile_TxTot  /='' .or. &
           ioFiles%inFile_TxPrfl /='' .or. &
           ioFiles%inFile_Jac    /='' ) THEN

         CALL convolPreviousResults( ioFiles%inFile_Rad, &
                                     ioFiles%inFile_TxTot, &
                                     ioFiles%inFile_TxPrfl, &
                                     ioFiles%inFile_Jac, &
                                     postCtrl, outGrid%V1, outGrid%V2, outGrid%DV )
      ENDIF

   END Program
