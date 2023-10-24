!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE Module_Config
   USE Module_ConstParam ,ONLY: r8=>kind_r8, FILLREAL,FILLINT

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_Output_Ctrl, &
             CLBLM_OutputSpectGrid, &
             CLBLM_Path_Ctrl, &
             CLBLM_DV_Ctrl, &
             CLBLM_OD_Ctrl, &
             CLBLM_RT_Ctrl, &
             CLBLM_FLUX_Ctrl, &
             CLBLM_Post_Ctrl, &
             CLBLM_IO_Files, &
             SceneFileOverRide, &
             inputUserDirectives, &
             openLogFile

   PUBLIC :: & !--- Single-variable control options and global variables
             NLTE, &
             ioFiles, IPR, noPrnt



   !---------------------------CLBLM_Output_Ctrl--------------------------
   TYPE :: CLBLM_Output_Ctrl
      integer                    :: OD  = 0     != 0: No OD output
                                                != 1: output layer OD
      integer                    :: Tx  = 0     != 0: No transmittance calculated
                                                != 1: Returns monochromatic total Tx (together with radiances, if Rad≠0)
                                                !=-1: Returns scanned total Tx
                                                != 2: Computes profile of monochromatic cumulative transmittance between observer and all levels up to TOA or down to surface depending on viewing direction
                                                !=-2: Returns scanned Tx profile
      integer                    :: Rad = 0     != 0: No radiance calculated
                                                != 1: monochromatic radiances
                                                !=-1: scanned radiances
      integer                    :: Jac = 0     != 0: No Jacobians calculated
                                                != 1: Returns monochromatic Jacobians
                                                !=-1: Returns convolved Jacobians
      character(20) ,allocatable :: JacList(:)  !like ["T","O3","CH4"], is a list of derivatives desired
   END TYPE


   !---------------------------CLBLM_OutputSpectGrid--------------------------
   TYPE CLBLM_OutputSpectGrid
       real    :: V1       = FILLREAL         ! from
       real    :: V2       = FILLREAL         ! to
       real    :: DV       = FILLREAL         ! resolution
       integer :: gridType = 0                ! if =-1, no adjustment; =0 using fixed DV ratios to adjust the calculated DV; if=1, uniform output DV.
   END TYPE


   !-------------------------CLBLM_Path_Ctrl-----------------------------
   TYPE :: CLBLM_Path_Ctrl
      real     :: V_refrac         = FILLREAL    !Frequency for refractive geometry calculations (default = (V1+V2)/2 )
      integer  :: zeroAbsAmt       = 0           !if =1 zeroes absorber amount which are <%0.1
      real     :: maxAlphaRat      = 1.5         !     Maximum Voigt width ratio across a layer
      real     :: maxTmpDiff(2)    = [5.,8.]     !(K)  Maximum layer temperature difference at TmpDiffAlt1 and 2TmpDiffAlt2
      real     :: maxTmpRefAlt(2)  = [0.,100.]   !(Km) Altitude of maxTmpDiff1 and maxTmpDiff2
      logical  :: AMscalingThermal = .TRUE.
      logical  :: AMScalingSolar   = .FALSE.
      integer  :: refPath          = 1           !if=0 vertical; if=1 up path
      integer  :: RTgridFlag       = 2           !if=1 RT grid from scene file header; if=2 auto layering; if=3 RT grid same as the input profile grid
   END TYPE


   !------------------ CLBLM_DV_Ctrl---------------------
   !
   TYPE :: CLBLM_DV_Ctrl
      real   :: SAMPLE      = 4.           !number of smaples per mean half-width
      real   :: ALFAL0      = 0.04         !average collision broadened halfwidth
      real   :: ALFAL0_fine = 0.005        !ALFAL0_fine used for narrow line DV for double line convolution,
   END TYPE


   !--------- CLBLM_OD_Ctrl-------------
   !
   TYPE :: CLBLM_OD_Ctrl
      !--- Line-by-line
      logical :: lineOff       = .FALSE.         !Turns calculation of lines contribution off
      logical :: speciesBroad  = .FALSE.         !use species by species broadening parameters.
      logical :: lineRejec     = .TRUE.
      real    :: DPTMIN        = 0.0000          !Minimum mol. optical depth. if <0 will be set to 0.0000
      real    :: DPTFAC        = 0.000           !Factor applying to F4 OD to determine line rejection; if <0 will be set to 0.001

      !--- Continuum
      logical :: contnmOff = .FALSE.             !Turns continuum contribution off. If=.TRUE. turn on continua and use individual scaling factors
      real    :: XSELF =1.
      real    :: XFRGN =1.
      real    :: XCO2C =1.
      real    :: XO3CN =1.
      real    :: XO2CN =1.
      real    :: XN2CN =1.
      real    :: XRAYL =1.

      !--- X-sections
      logical :: xsPressConv =.TRUE.             !=.TRUE., cross-section convolved with pressure
      logical :: XSectOff   =.FALSE.             !Turns calculation of cross-sections contribution off

      !--- Additional options
      !integer :: ILBLF4 =1                       !flag for F4; !=0, LBLF4 not activated, =1 line-by-line bound is 25cm-1; =2 25cm-1 bound for layers with pressures>0.5mb and 5cm-1 for other layers.
      integer :: JRAD   =1                       ! JRAD=  1  RADIATION TERM INCLUDED IN LINE STRENGTHS
                                                 ! JRAD=  0  RADIATION TERM PUT IN BY PANEL
                                                 ! JRAD= -1  NO RADIATION TERM IN ABSORPTION COEFFICIENTS
   END TYPE



   !------------ CLBLM_RT_Ctrl------------------
   !
   TYPE :: CLBLM_RT_Ctrl

      logical :: ThermalOn = .TRUE.          ! If =.F., No thermal emission; =1 Include thermal emission
      integer :: linInTau = 1                ! If =0, isothermal layer; =1 LBLRTM linear-in-tau approximation; =2, standard linear-in-tau formular.


      logical :: SolarOn = .FALSE.           ! If =.F., No solar source; =.T. Include solar source
                                             !
      integer :: solarSource = 1             ! If=1, use Kurucz solar irradiance data; one component only; Solar constant equals 1368.22 Wm-2, unless scaled by "solarConst".
                                             ! If=2, use NRLSSI2 solar irradiance data; 1 or 3 components, depending on solarVarOption.
                                             !       Solar constant for 1 component is 1360.85 Wm-2, unless scaled by "solarConst".
                                             !       ( The NRLSSI2 1-component solar constant is for the spectral range 100-50000 cm-1 with quiet sun,
                                             !         facular and sunspot contributions fixed to the mean of solar Cycles 13-24
                                             !         and averaged over the mean solar cycle )
                                             !
      integer :: solarVarOption = 1          ! =1 No variability.
                                             ! =2 Facular brightening and sunspot blocking amplitudes
                                             !    are by default determined by the fraction of the
                                             !    way into the solar cycle (see SOLCYCFRAC) or can
                                             !    be scaled (see faculaVar and spotVar).
                                             !    This option needs the NRLSSI2 3-component data to be selected.
                                             ! =3 Facular brightening and sunspot blocking amplitudes
                                             !    are determined by the Mg and SB indeces (see faculaVar and spotVar)
                                             !    This option needs the NRLSSI2 3-component data to be selected.
                                             !
      real    :: solarConst = 0.             ! Specifies solar constant
                                             !    =0.0 no scaling of internal solar irradiance
                                             !    >0.0 and 1-component solar data, total solar irradiance is scaled to solarConst
                                             !    >0.0 and 3-component solar data and solarVarOption = 2, integral of total solar irradiance
                                             !         averaged over solar cycle is scaled to solarConst
                                             !
      real    :: solCycFrac = FILLREAL       ! (0-1) Solar cycle fraction (solarVarOption=2 only).
                                             !       fraction of the way through the mean 11-year
                                             !       cycle with 0 and 1 defined as the minimum phase of the solar cycle
                                             !
      real    :: faculaVar = FILLREAL        ! Solar variability scaling factors or indices (solarVarOption=2,3 only)
      real    :: spotVar   = FILLREAL        !    if solarVarOption = 2 =>
                                             !       faculaVar: Facular (Mg) index amplitude scale factor
                                             !       spotVar:   Sunspot (SB) index amplitude scale factor
                                             !    if solarVarOption = 3 =>
                                             !       faculaVar: Facular (Mg) index as defined in the NRLSSI2 model;
                                             !                  used for modeling specific solar activity
                                             !       spotVar:   Sunspot (SB) index as defined in the NRLSSI2 model;
                                             !                  used for modeling specific solar activity
      integer :: JulDay = FILLINT            ! Julian day, range from 0~366; if = 0, no scaling of solar source function by Julian day


      !--- For RT with scattering
      !integer :: multiScatt     = 0
      !integer :: numStreams
      !integer :: numPhaseMoments

   END TYPE



   !------------ CLBLM_Post_Ctrl--------------
   ! "functID" is a integer from 1 to 14, for built-in scanning functions .
   ! If functID=0: filter function is read in from file.
   ! 1          Boxcar
   ! 2          Triangle
   ! 3          Gauss
   ! 4          Sinc
   ! 5          Sinc2
   ! 6          Beer
   ! 7          Hamming
   ! 8          Hanning
   ! 9          Norton-Beer
   ! 10,11,12,  Brault
   ! 13         Kaiser-Bessel
   ! 14         Kiruna
   !
   TYPE :: CLBLM_Post_Ctrl
      integer    :: functID           =FILLINT  !used by FFT,scan.filt,  Select scanning/apodization function (selection depends on FFT) or filtering.
      logical    :: FFT               =.FALSE.  !used by FFT             If .TRUE. performs convolution using Fast Fourier Transform,  otherwise subroutine performs convolution in wavenumber space.
      real       :: HWHM              =FILLREAL !used by FFT,scan.       Half Width Half Maximum of scanning function
      real       :: functParams(3)    =FILLREAL !used by FFT.            Only required if FFT= .TRUE. and FuncID =11, 12 or 13
      real       :: boxcarHW          =FILLREAL !used by FFT.            If<0, no boxcar averaging performed; if >0 real number indicating the width of the boxcar function in cm-1 when pre-boxcaring
      logical    :: deconvPreScan     =.FALSE.  !used by FFT.            IF=T, deconvolve the scanned spectrum with the boxcar
   END TYPE

!------------ CLBLM_FLUX_Ctrl------------------
   !
   TYPE :: CLBLM_FLUX_Ctrl

      logical :: FLUX_FLAG      = .FALSE.          ! If =.F., No flux calc; =1 Include flux calculation
      real    :: DV_FLUX        = FILLREAL         ! resolution for flux calculations
      integer :: nang           = 1                ! =1 perform flux calculations for one quadrature angle
                                              ! =2 perform flux calculations for quadrature angles
                                              ! =3 perform flux calculations for quadrature angles
   END TYPE

   !---------------------------CLBLM_IO_Files--------------------------------
   !
   TYPE CLBLM_IO_Files
      !pdc = '/' !path delimiter chracter


      !--- Static input data
      character(256) :: lineFile            = 'clblm_data/spectroscopy/TAPE3'                       ! Line data parameters
      character(200) :: xsFilePath          = 'clblm_data/spectroscopy/xs/'                         ! Cross-section data. !path string plus delimiting character, '/' for Linux
      character(256) :: nlteStatPopFile     = 'clblm_data/TAPE4'                                    ! NLTE vibrational temperature profiles
      character(200) :: solarPath           = 'clblm_data/solar_irradiance/'
      character(56)  :: solarFile_Kurucz    = 'SOLAR.RAD.nc'
      character(56)  :: solarFile_1comp_NRL = 'SOLAR.RAD.nc'
      character(56)  :: solarFile_3comp_NRL = 'SOLAR.RAD.nc'

      !---Scratch files
      character(200) :: scratchPath     = 'scratch/'
      character(256) :: lineF4File      = 'scratch/TAPE9'                  ! Scratch file for line F4 convolution
      character(256) :: narrowLineFile  = 'scratch/narrowTAPE3.dat'        ! Scratch file for narrow line convolution

      !--- Input data
      character(256) :: sceneFile        = 'user_archive/scene_files/scenes.nc'   ! Scene file
      character(200) :: cloudAerosolPath = 'user_archive/cloud_aerosol/'
      character(200) :: surfacePath      = 'user_archive/surface/'
      character(200) :: filterFunctPath  = 'user_archive/sensor/'
      character(256) :: filterFunctFile  = ''                               !user must input a filter file name


      !---
      character(200) :: clblmInPath   = 'clblm_out/'
      character(200) :: inFile_OD     = ''
      character(200) :: inFile_Rad    = ''
      character(200) :: inFile_TxTot  = ''
      character(200) :: inFile_TxPrfl = ''
      character(200) :: inFile_Jac    = ''


      !--- Output data
      character(200) :: clblmOutPath    = 'clblm_out/'
      character(200) :: outPath_OD      = ''
      character(200) :: outPath_Rad     = ''
      character(200) :: outPath_TxTot   = ''
      character(200) :: outPath_TxPrfl  = ''
      character(200) :: outPath_Jac     = ''
      character(20)  :: rootName_OD     = 'od'        !e.g. od_sNNN_LLL.nc
      character(20)  :: rootName_Rad    = 'rad'       !e.g. rad_mono_o10km-a40deg_sNNN.nc
      character(20)  :: rootName_TxTot  = 'tx-total'  !e.g. tx-total_mono_sNNN.nc
      character(20)  :: rootName_TxPrfl = 'tx-prfl'   !e.g. tx-prfl_mono_o10km-a40deg_sNNN_LLL.nc
      character(20)  :: rootName_Jac    = 'drad-d'    !e.g. drad-dT_mono_sNNN_LLL.nc

      character(256) :: logFile         = 'clblm_out/TAPE6'
   END TYPE



   !---------------------- CLBLM_SceneFileOverRide -----------------------------
   TYPE SceneFileOverRide
      !character(256)        :: sceneFile = ''
      integer               :: numSelectedScenes = FILLINT
      integer ,allocatable  :: selectedScene_ID(:)
      real    ,allocatable  :: obsAlt(:)
      real    ,allocatable  :: viewAng(:)
   END TYPE


   !----------------------Global variables -----------------------------
   !
   ! "ioFiles" is a globle object to store all the files.
   ! "IPR" and "noPrnt" are used so often, so make it a global variable.
   TYPE(CLBLM_IO_Files) ,SAVE :: ioFiles
   integer              ,SAVE :: IPR
   integer              ,SAVE :: noPrnt = 0

   !--- Single variable options
   logical :: NLTE = .FALSE.


CONTAINS !=================== MODULE CONTAINS ==========================


 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   SUBROUTINE openLogFile( logFileName )
 !----------------------------------------------------------------------
      USE Module_Utility ,ONLY: getLun
      IMPLICIT NONE

      character(*) ,intent(in),OPTIONAL :: logFileName
      logical :: op
      character(256) :: logFile

      if (present(logFileName)) then
         logFile = logFileName
      else
         logFile = ioFiles%logFile
      endif

      inquire( FILE=trim(logFile), opened=op )
      if (.not.op) then
         IPR = getLun() !initialize the global variable IPR
         OPEN (IPR,FILE=trim(logFile),STATUS='UNKNOWN')
      endif
   END SUBROUTINE

 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   SUBROUTINE inputUserDirectives( cfgFile, outCtrl, outGrid, pathCtrl,dvCtrl,odCtrl, &
                                   rtCtrl,postCtrl, sceneOverride,fluxCtrl )
 !----------------------------------------------------------------------
      USE Module_Utility ,ONLY: getLun
      USE jsonconfig
      use json_data_types
      IMPLICIT NONE

      character(*)                  ,intent(in)  :: cfgFile
      type(CLBLM_Output_Ctrl)       ,intent(out) :: outCtrl
      type(CLBLM_OutputSpectGrid)   ,intent(out) :: outGrid
      type(CLBLM_Path_Ctrl)         ,intent(out) :: pathCtrl
      type(CLBLM_DV_Ctrl)           ,intent(out) :: dvCtrl
      type(CLBLM_OD_Ctrl)           ,intent(out) :: odCtrl
      type(CLBLM_RT_Ctrl)           ,intent(out) :: rtCtrl
      type(CLBLM_Post_Ctrl)         ,intent(out) :: postCtrl
      type(SceneFileOverRide)       ,intent(out) :: sceneOverride
      type(CLBLM_FLUX_Ctrl)         ,intent(out) :: fluxCtrl


      integer                       :: lun
      logical                       :: res
      character(len=1024)           :: errStr
      type(json_configuration_type) :: config


      !type json_configuration_Type
      !    logical                               :: NLTE
      !    type(clblm_in_type)                   :: clblm_in
      !    type(clblm_out_type)                  :: clblm_out
      !    type(rt_flags_type)                   :: rt_flags
      !    type(path_type)                       :: path
      !    type(od_flags_Type)                   :: od_flags
      !    type(spectral_convolution_flags_type) :: spectral_convolution_flags
      !    type(output_spectral_grid_type)       :: output_grid
      !    type(target_viewing_type)             :: target_viewing
      !    type(scene_selection_type)            :: scene_selection
      !    type(geometry_Type)                   :: geometry
      !    type(solar_irradiance_type)           :: solar_irradiance
      !end type


      lun = getLun()
      open(lun, file=cfgFile, status='old')

      !--- Read in Jason file contents
      res = getConfigFileStruct( cfgFile, config, errStr )
      if (.not. res) then
          write(IPR, *)'Unable to find and/or load JSON configuration file.  Error = ', &
              errStr(1:len(trim(errStr))), '.  Exiting...'
          print*, 'Unable to read JSON file. '//trim(errStr)
          STOP
      endif


      !--- sceneSelectin and geometry
      call covertSceneOverride( config%scene_selection, config%geometry, sceneOverride )

      !--- clblm_out
      call convertOutCtrl( config%clblm_out, outctrl)

      !--- output_spectral_grid
      call convertOutputGrid( config%output_grid, outGrid)

      !--- path
      call convertPath( config%path, pathctrl)

      !--- od_flags
      call convertOdFlags( config%od_flags, odctrl)

      !--- clblm_in
      call convertInCtrl( config%clblm_in )

      !--- rt_flags
      call convertRtFlags( config%rt_flags, config%solar_variability, rtctrl)

      !--- flux_flags
      call convertFluxFlags( config%flux_flags, fluxctrl)

      !--- convolution_flags
      call convertSpectralFlags( config%spectral_convolution_flags, postctrl)

      !--- Assign the singular key word, NLTE
      if (    config%nlte_flag%NLTE==0) then; NLTE=.FALSE.;
      elseif (config%nlte_flag%NLTE==1) then; NLTE=.TRUE.;
      else
         NLTE = .FALSE.
      endif



      !--- Modify default values
      !
      ! Set default value for V_refrac
      if ( pathCtrl%V_refrac <0) then !User not provided the V_refrac, set to default value.
         pathCtrl%V_refrac = 0.5*( outGrid%V1 + outGrid%V2 )
      endif

      ! IF OD-only or Tx-only mode turn off thermal and solar sources.
      if ( outCtrl%Rad==0 .and. outCtrl%Jac==0) then
         rtCtrl%ThermalOn = .FALSE.
         rtCtrl%SolarOn = .FALSE.
      endif

      ! If solarVarOption==3, solarConst is ignored
      if (rtCtrl%solarVarOption ==3) rtCtrl%solarConst = 0.

      ! If "clblm_in" group exist, user is requesting to do convolution for
      ! previous results, turn off OD/Tx/Rad/Jac calculations in outCtrl.
      ! outCtrl%JacList kept unchanged because convolution-only mode may need it.
      if ( ioFiles%inFile_Rad    /='' .or. &
           ioFiles%inFile_TxTot  /='' .or. &
           ioFiles%inFile_TxPrfl /='' .or. &
           ioFiles%inFile_Jac    /='' ) then  !Do convolution only processing
         outCtrl%OD  = 0
         outCtrl%Tx  = 0
         outCtrl%Rad = 0
         outCtrl%Jac = 0
      endif


      !--- Check the user directives
      CALL check_userDirectives( outCtrl, outGrid, pathCtrl, dvCtrl, odCtrl, &
                                 rtCtrl, postCtrl, sceneOverride, NLTE, fluxCtrl )

   END SUBROUTINE



 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   subroutine covertSceneOverride( scene_selection, geometry, sceneOverride )
 !----------------------------------------------------------------------
      use json_data_types
      implicit none

      type(scene_selection_type) ,intent(in)  :: scene_selection
      type(geometry_type)        ,intent(in)  :: geometry
      type(SceneFileOverRide)    ,intent(out) :: sceneOverride


      ! Set scene file name
      if ( trim( scene_selection%scene_file ) /='' ) then !user provided scene file name
        ioFiles%sceneFile = trim( scene_selection%scene_file )
      endif
      !sceneOverride%sceneFile = trim(scene_selection%scene_file)

      sceneOverride % numSelectedScenes = scene_selection%nscenes

      if ( allocated(sceneOverride % selectedScene_ID)) deallocate(sceneOverride % selectedScene_ID)
      if ( allocated( scene_selection%scene_id ) ) then
         allocate( sceneOverride % selectedScene_ID( size(scene_selection%scene_id) ) )
         sceneOverride % selectedScene_ID(:) = scene_selection%scene_id(:)
      endif

      if ( allocated(sceneOverride % obsAlt)) deallocate(sceneOverride % obsAlt)
      if ( allocated(geometry%obs_altitudes)) then
         allocate( sceneOverride % obsAlt( size(geometry%obs_altitudes) ) )
         sceneOverride % obsAlt(:) = geometry%obs_altitudes(:)
      endif

      if ( allocated(sceneOverride % viewAng)) deallocate(sceneOverride % viewAng)
      if ( allocated(geometry%view_angles)) then
         allocate( sceneOverride % viewAng( size(geometry%view_angles) ) )
         sceneOverride % viewAng(:) = geometry%view_angles(:)
      endif

   end subroutine


 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   subroutine convertInCtrl(input)
 !----------------------------------------------------------------------
      use json_data_types
      implicit none

       type(clblm_in_type), intent(in) :: input
       character(len=MAX_STRING_LENGTH) :: basename

       if (len(trim(input%od_in)) .gt. 0) then
           !call splitFilePath(input%od_in, iofiles%outPath_OD, basename)
           iofiles%inFile_OD = trim(input%od_in)
       !else
       !    iofiles%inFile_OD = trim(iofiles%clblmInPath)
       endif

       if (len(trim(input%total_tx_in)) .gt. 0) then
           !call splitFilePath(input%od_in, iofiles%outPath_OD, basename)
           iofiles%inFile_TxTot = trim(input%total_tx_in)
       !else
       !    iofiles%inFile_TxTot = trim(iofiles%clblmInPath)
       endif

       if (len(trim(input%tx_profile_in)) .gt. 0) then
           !call splitFilePath(input%od_in, iofiles%outPath_OD, basename)
           iofiles%inFile_TxPrfl = trim(input%tx_profile_in)
       !else
       !    iofiles%inFile_TxPrfl = trim(iofiles%clblmInPath)
       endif

       if (len(trim(input%rad_in)) .gt. 0) then
           !call splitFilePath(input%od_in, iofiles%outPath_OD, basename)
           iofiles%inFile_Rad = trim(input%rad_in)
       !else
       !    iofiles%inFile_Rad = trim(iofiles%clblmInPath)
       endif

       if (len(trim(input%jacobians_in)) .gt. 0) then
           !call splitFilePath(input%od_in, iofiles%outPath_OD, basename)
           iofiles%inFile_Jac = trim(input%jacobians_in)
       !else
       !    iofiles%inFile_Jac = trim(iofiles%clblmInPath)
       endif
   end subroutine


 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
    SUBROUTINE convertOutCtrl(clblm_out, outCtrl)
 !----------------------------------------------------------------------
      use json_data_types
      implicit none

      type(clblm_out_type),    intent(in)  :: clblm_out
      type(CLBLM_Output_Ctrl), intent(out) :: outCtrl

      character(len=MAX_STRING_LENGTH) :: string
      integer :: i, ilen, kret


      outCtrl%od = 0
      if (trim(clblm_out%od%solver_type) .eq. "mono") then
        outCtrl%od = 1
        if (len(trim(clblm_out%od%path)) .gt. 0) &
              iofiles%outPath_OD = trim(clblm_out%od%path)
      endif


      ! TX
      outCtrl%tx = 0
      if ((trim(clblm_out%total_tx%solver_type) .eq. "mono") .or. &
          (trim(clblm_out%total_tx%solver_Type) .eq. "convolved")) then

          if (trim(clblm_out%total_tx%solver_Type) .eq. "mono") then
              outCtrl%tx = 1
          else
              outCtrl%tx = -1
          endif
          if (len(trim(clblm_out%total_tx%path)) .gt. 0) &
              iofiles%outPath_TxTot = trim(clblm_out%total_tx%path)
      endif

      if ((trim(clblm_out%tx_profile%solver_type) .eq. "mono") .or. &
          (trim(clblm_out%tx_profile%solver_Type) .eq. "convolved")) then

          if (trim(clblm_out%tx_profile%solver_Type) .eq. "mono") then
              outCtrl%tx = 2
          else
              outCtrl%tx = -2
          endif
          if (len(trim(clblm_out%tx_profile%path)) .gt. 0) &
              iofiles%outPath_TxPrfl = trim(clblm_out%tx_profile%path)
      endif

      ! Rad
      outCtrl%rad = 0
      if ((trim(clblm_out%rad%solver_type) .eq. "mono") .or. &
          (trim(clblm_out%rad%solver_Type) .eq. "convolved")) then

          if (trim(clblm_out%rad%solver_Type) .eq. "mono") then
              outCtrl%rad = 1
          else
              outCtrl%rad = -1
          endif
          if (len(trim(clblm_out%rad%path)) .gt. 0) &
              iofiles%outPath_Rad = trim(clblm_out%rad%path)
      endif

      ! JAcobians
      outCtrl%jac = 0
      if ((trim(clblm_out%jacobians%solver_type) .eq. "mono") .or. &
          (trim(clblm_out%jacobians%solver_Type) .eq. "convolved")) then

          if (trim(clblm_out%jacobians%solver_Type) .eq. "mono") then
              outCtrl%jac = 1
          else
              outCtrl%jac = -1
          endif
          if (len(trim(clblm_out%jacobians%path)) .gt. 0) &
              iofiles%outPath_Jac = trim(clblm_out%jacobians%path)
      endif

      if (allocated(clblm_out%jacobian_list)) then
          if (allocated(outCtrl%jaclist)) deallocate(outCtrl%jaclist)
          allocate(outCtrl%jaclist(size(clblm_out%jacobian_list)), stat=kret)
          if (kret .ne. 0) then
              write(ipr, *) "unable to allocate clblm_output_ctrl%jaclist"
              stop
          endif
          do i = 1, size(clblm_out%jacobian_list)
            outCtrl%jaclist(i) = '' !initialize the string
              ilen = len(trim(clblm_out%jacobian_list(i)))
              if (ilen .gt. len(outCtrl%jaclist(i))) ilen = len(outCtrl%jaclist(i))
              outCtrl%jaclist(i)(1:ilen) = clblm_out%jacobian_list(i)(1:ilen)
          enddo
      endif

    END SUBROUTINE




 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   subroutine convertRtFlags(rt_flags, solar_variability, rtCtrl)
 !----------------------------------------------------------------------
   USE Module_Utility, ONLY: upper
   use json_data_types
   implicit none

   type(rt_flags_type)          ,intent(in)  :: rt_flags
   type(solar_variability_type) ,intent(in)  :: solar_variability
   type(clblm_rt_ctrl)          ,intent(out) :: rtCtrl

   if (    rt_flags%thermal_source==0) then; rtCtrl%thermalon=.FALSE.;
   elseif (rt_flags%thermal_source==1) then; rtCtrl%thermalon=.TRUE.;
   endif

   if (rt_flags%linear_in_tau >=0) rtCtrl%linintau  = rt_flags%linear_in_tau

   if( upper(trim(rt_flags%solar_source)) == "KURUCZ") then
      rtCtrl%solarOn     = .TRUE.
      rtCtrl%solarSource = 1
   elseif( upper(trim(rt_flags%solar_source)) == "NRL") then
      rtCtrl%solarOn     = .TRUE.
      rtCtrl%solarSource = 2
   else !if( upper(trim(rt_flags%solar_source)) == "") then
      rtCtrl%solarOn = .FALSE.
   endif

   if (rt_flags%julday >=0)               rtCtrl%JulDay         = rt_flags%julday
   if (rt_flags%solar_cnst >=0.)          rtCtrl%SolarConst     = rt_flags%solar_cnst
   if (solar_variability%option >=1)      rtCtrl%solarVarOption = solar_variability%option
   if (solar_variability%cycle_frac >=0.) rtCtrl%solCycFrac     = solar_variability%cycle_frac

   rtCtrl%faculaVar = solar_variability%facula_var
   rtCtrl%spotVar   = solar_variability%spot_var

   end subroutine

 !----------------------------------------------------------------------

 !----------------------------------------------------------------------
   subroutine convertSpectralFlags(spectral_convolution_flags, postCtrl)
 !----------------------------------------------------------------------
      use json_data_types
      implicit none

      type(spectral_convolution_flags_type), intent(in)  :: spectral_convolution_flags
      type(clblm_post_ctrl),                 intent(out) :: postCtrl


      ! If we have a filter filename, copy that in
      if ( trim( spectral_convolution_flags%filter_file) /='' ) then
         iofiles%filterFunctFile = trim( spectral_convolution_flags%filter_file)
      endif

      postCtrl%functid     = spectral_convolution_flags%function_id
      postCtrl%hwhm        = spectral_convolution_flags%hwhm
      postCtrl%boxcarHW    = spectral_convolution_flags%averaging_width

      if (    spectral_convolution_flags%fft==0) then; postCtrl%fft=.FALSE.
      elseif (spectral_convolution_flags%fft==1) then; postCtrl%fft=.TRUE.
      endif

      if (postCtrl%fft .and. (postCtrl%functid==12 .or. postCtrl%functid==13)) then
        postCtrl%functparams(1) = spectral_convolution_flags%function_params(1)
      elseif (postCtrl%fft .and. postCtrl%functid==14) then
        postCtrl%functparams(1:3) = spectral_convolution_flags%function_params(2:4)
      endif

   end subroutine

 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   subroutine convertPath(path_flags, pathCtrl)
 !----------------------------------------------------------------------
    use json_data_types
    implicit none

    type(path_type), intent(in) :: path_flags
    type(clblm_path_ctrl), intent(out) :: pathCtrl

    pathCtrl%v_refrac = path_flags%v_refrac

    if (    path_flags%airmass_scaling(1)==0) then;  pathCtrl%AMscalingThermal = .FALSE.;
    elseif (path_flags%airmass_scaling(1)==1) then;  pathCtrl%AMscalingThermal = .TRUE.;
    endif

    if (    path_flags%airmass_scaling(2)==0) then; pathCtrl%AMScalingSolar = .FALSE.;
    elseif (path_flags%airmass_scaling(2)==1) then; pathCtrl%AMScalingSolar = .TRUE.;
    endif

    if (path_flags%reference_path >=0) pathCtrl%refPath    = path_flags%reference_path !if=0, vertical; if=1, up path
    if (path_flags%RT_grid >=1)        pathCtrl%RTgridFlag = path_flags%RT_grid

    !pathCtrl%zeroAbsAmt       =   !if =1 zeroes absorber amount which are <%0.1
    !pathCtrl%maxAlphaRat      =   !     Maximum Voigt width ratio across a layer
    !pathCtrl%maxTmpDiff1      =   !(K)  Maximum layer temperature difference at TmpDiffAlt1
    !pathCtrl%maxTmpDiff2      =   !(K)  Maximum layer temperature difference at TmpDiffAlt2
    !pathCtrl%TmpDiffAlt1      =   !(Km) Altitude of maxTmpDiff1
    !pathCtrl%TmpDiffAlt2      =   !(Km) Altitude of maxTmpDiff2

   end subroutine

 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   subroutine convertOdFlags(od_flags, odCtrl)
 !----------------------------------------------------------------------
    use json_data_types
    implicit none

    type(od_flags_Type), intent(in)  :: od_flags
    type(clblm_od_ctrl), intent(out) :: odCtrl

    integer :: i


    if (    od_flags%lines_contribution==0) then; odCtrl%lineoff = .TRUE.;
    elseif (od_flags%lines_contribution==1) then; odCtrl%lineoff = .FALSE.;
    endif

    !if (od_flags%dptmin>0.) odCtrl%dptmin = od_flags%dptmin
    !if (od_flags%dptfac>0.) odCtrl%dptfac = od_flags%dptfac
    if (od_flags%line_rejection_params(1)>=0.) odCtrl%dptmin = od_flags%line_rejection_params(1)
    if (od_flags%line_rejection_params(2)>=0.) odCtrl%dptfac = od_flags%line_rejection_params(2)

    if (    od_flags%collision_partners_broadening==0) then; odCtrl%speciesbroad=.FALSE.;
    elseif (od_flags%collision_partners_broadening==1) then; odCtrl%speciesbroad=.TRUE.;
    endif

    if (    od_flags%line_rejection==0) then; odCtrl%linerejec=.FALSE.;
    elseif (od_flags%line_rejection==1) then; odCtrl%linerejec=.TRUE.;
    endif

    if (    od_flags%continuum_contribution==0) then; odCtrl%contnmOff=.TRUE.;
    elseif (od_flags%continuum_contribution==1) then; odCtrl%contnmOff=.FALSE.;
    endif

    if (    od_flags%p_convolution==0) then; odCtrl%xsPressConv=.FALSE.;
    elseif (od_flags%p_convolution==1) then; odCtrl%xsPressConv=.TRUE.;
    endif

    !--- Continuum scaling factors. Must be in the order of
    !   [XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL]
    if (od_flags%continuum_scaling(1)>=0.) odCtrl%XSELF = od_flags%continuum_scaling(1)
    if (od_flags%continuum_scaling(2)>=0.) odCtrl%XFRGN = od_flags%continuum_scaling(2)
    if (od_flags%continuum_scaling(3)>=0.) odCtrl%XCO2C = od_flags%continuum_scaling(3)
    if (od_flags%continuum_scaling(4)>=0.) odCtrl%XO3CN = od_flags%continuum_scaling(4)
    if (od_flags%continuum_scaling(5)>=0.) odCtrl%XO2CN = od_flags%continuum_scaling(5)
    if (od_flags%continuum_scaling(6)>=0.) odCtrl%XN2CN = od_flags%continuum_scaling(6)
    if (od_flags%continuum_scaling(7)>=0.) odCtrl%XRAYL = od_flags%continuum_scaling(7)

   end subroutine

   !----------------------------------------------------------------------
   subroutine convertFluxFlags(flux_flags, fluxCtrl)
   !----------------------------------------------------------------------
        use json_data_types
        implicit none

        type(flux_flags_type)            ,intent(in)  :: flux_flags
        type(clblm_flux_ctrl)            ,intent(out) :: fluxCtrl

        if (    flux_flags%flux_flag==0) then; fluxCtrl%flux_flag=.FALSE.;
        elseif (flux_flags%flux_flag==1) then; fluxCtrl%flux_flag=.TRUE.;
        endif

        fluxCtrl%dv_flux = flux_flags%dv_flux
        fluxCtrl%nang = flux_flags%nang

        end subroutine
 !----------------------------------------------------------------------

 !----------------------------------------------------------------------
   subroutine convertOutputGrid(source, outGrid)
 !----------------------------------------------------------------------
    use json_data_types
    implicit none

    type(output_spectral_grid_type),intent(in)  :: source
    type(CLBLM_OutputSpectGrid)  ,intent(out) :: outGrid

    integer :: gTyp

    if (trim(source%grid_type) .eq. "uniformDV") then;
        gTyp = 1
    else if ((trim(source%grid_type) .eq. "adjustedDV") .or. &
             (trim(source%grid_type) .eq. '')) then
        gTyp = 0
    else
        gTyp = -1
    endif

    outGrid%v1 = source%v1
    outGrid%v2 = source%v2
    outGrid%dv = source%dv
    outGrid%gridType = gTyp

   end subroutine





!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE check_userDirectives( outCtrl, outGrid, pathCtrl, dvCtrl, odCtrl, &
                                    rtCtrl, postCtrl, sceneOverride, NLTE, fluxCtrl )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY:FILLINT
      use json_data_types
      IMPLICIT NONE

      type(CLBLM_Output_Ctrl)      ,intent(in) :: outCtrl
      type(CLBLM_OutputSpectGrid)  ,intent(in) :: outGrid
      type(CLBLM_Path_Ctrl)        ,intent(in) :: pathCtrl
      type(CLBLM_DV_Ctrl)          ,intent(in) :: dvCtrl
      type(CLBLM_OD_Ctrl)          ,intent(in) :: odCtrl
      type(CLBLM_RT_Ctrl)          ,intent(in) :: rtCtrl
      type(CLBLM_FLUX_Ctrl)        ,intent(in) :: fluxCtrl
      type(CLBLM_Post_Ctrl)        ,intent(in) :: postCtrl
      type(SceneFileOverRide)      ,intent(in) :: sceneOverride
      logical                      ,intent(in) :: NLTE


      character(*) ,parameter :: routineName='check_userDirectives'
      logical :: fileExist



      !--- Check outCtrl ---------------------------------------------
      !
      if ( .NOT.( outCtrl%OD ==0 .or. &    != 0: No OD output
                  outCtrl%OD ==1 ) ) then  != 1: output layer OD;
         STOP '--- '//routineName//'(): Invalid value for OD output flag.'//' Program stopped.'
      endif
      if ( .NOT.( abs(outCtrl%Tx) <=2) ) then
         STOP '--- '//routineName//'(): Invalid value for Tx output flag.'//' Program stopped.'
      endif
      if ( .NOT.( abs(outCtrl%Rad) <=1) ) then
         STOP '--- '//routineName//'(): Invalid value for Rad output flag.'//' Program stopped.'
      endif
      if ( .NOT.( abs(outCtrl%Jac) <=1) ) then
         STOP '--- '//routineName//'(): Invalid value for Jac output flag.'//' Program stopped.'
      endif

      !if ( (outCtrl%Jac==0 .and. outCtrl%OD/=0 .and. outCtrl%Tx/=0 .and. outCtrl%Rad/=0) ) then
      !   STOP '--- '//routineName//'(): No RT calculation selected.'//' Program stopped.'
      !endif

      if ( outCtrl%Jac/=0 .and. (outCtrl%OD/=0 .or. outCtrl%Tx/=0) ) then
         STOP '--- '//routineName//'(): Jacobian calculation will not output OD or Tx.'//' Program stopped.'
      endif

      if ( outCtrl%Jac/=0 .and. .not.allocated(outCtrl%JacList) .and. outCtrl%Rad==0 ) then
         STOP '--- '//routineName//'(): Jacobians requested but Jacobians list is missing.'//' Program stopped.'
      endif



      !--- Check outGrid ---------------------------------------------
      !
      if ( outGrid%V1 <0. .or. outGrid%V2 <0. ) then
         STOP '--- '//routineName//'(): Invalid value for V1 and V2.'//' Program stopped.'
      endif
      if ( outGrid%V2 < outGrid%V1 ) then
         STOP '--- '//routineName//'(): Spectral limits V2 < V1.'//' Program stopped.'
      endif
      if ( .NOT.( outGrid%gridType ==-1 .or. &    !if =-1, exact DV, no adjustment;
                  outGrid%gridType == 0 .or. &    !if =0 using fixed DV ratios to adjust the calculated DV;
                  outGrid%gridType == 1 ) ) then  !if=1, uniform output DV.
         STOP '--- '//routineName//'(): Invalid value for gridType flag.'//' Program stopped.'
      endif



      !--- Check atmCtrl------------------------------------------------
      ! V_refrac, zeroAbsAmt
      !
      if ( pathCtrl%V_refrac<0. ) then
         STOP '--- '//routineName//'(): Invalid value for V_refrac.'//' Program stopped.'
      endif
      if ( .NOT.( pathCtrl%zeroAbsAmt==0 .or. &  !donot zero absorber amount
                  pathCtrl%zeroAbsAmt==1) ) then !zeroes absorber amount which are <%0.1
         STOP '--- '//routineName//'(): Invalid value for zeroAbsAmt.'//' Program stopped.'
      endif
      if ( .NOT.( pathCtrl%refPath==0 .or. &  !vertical
                  pathCtrl%refPath==1) ) then !upPath
         STOP '--- '//routineName//'(): Invalid value for refPath.'//' Program stopped.'
      endif
      if ( .NOT.( pathCtrl%RTgridFlag==1 .or. &  !use pRT from scene file header
                  pathCtrl%RTgridFlag==2 .or. &  !auto layering
                  pathCtrl%RTgridFlag==3 ) )then  !user profile z grid
         STOP '--- '//routineName//'(): Invalid value for RTgridFlag.'//' Program stopped.'
      endif



      !--- Check odCtrl -----------------------------------------------
      !
      if ( .not.odCtrl%contnmOff .and. &
           any( [odCtrl%XSELF,&
                 odCtrl%XFRGN,&
                 odCtrl%XCO2C,&
                 odCtrl%XO3CN,&
                 odCtrl%XO2CN,&
                 odCtrl%XN2CN,&
                 odCtrl%XRAYL] <0. ) ) then
         STOP '--- '//routineName//': Invalid value for one or more continuum absorption multiplicative factors.'//' Program stopped.'
      endif

      !if ( odCtrl%contnmOff .and. &
      !     any( abs([odCtrl%XSELF,&
      !               odCtrl%XFRGN,&
      !               odCtrl%XCO2C,&
      !               odCtrl%XO3CN,&
      !               odCtrl%XO2CN,&
      !               odCtrl%XN2CN,&
      !               odCtrl%XRAYL]) >0. ) ) then
      !   STOP '--- '//routineName//'(): Invalid value for one or more continuum absorption multiplicative factors.'//' Program stopped.'
      !endif



      !--- Check dvCtrl ------------------------------------------------
      !
      if ( dvCtrl%ALFAL0 <0. ) then
         STOP '--- '//routineName//'(): Invalid value for ALFAL0.'//' Program stopped.'
      endif
      if ( dvCtrl%ALFAL0_fine <0. ) then
         STOP '--- '//routineName//'(): Invalid value for ALFAL0_fine.'//' Program stopped.'
      endif



      !--- Check rtCtrl ------------------------------------------------
      !
      if ( .NOT.( rtCtrl%linInTau ==0 .or. &    ! =0, isothermal layer;
                  rtCtrl%linInTau ==1 .or. &    ! =1 LBLRTM linear-in-tau approximation;
                  rtCtrl%linInTau ==2 ) ) then  ! =2, standard linear-in-tau formular.
         STOP '--- '//routineName//'(): Invalid value for linear-in-tau flag.'//' Program stopped.'
      endif
      if ( .NOT.( rtCtrl%solarSource == 1 .or. &
                  rtCtrl%solarSource == 2 ) ) then
         STOP '--- '//routineName//'(): Invalid value for solar source flag.'//' Program stopped.'
      endif
      if ( .NOT.( rtCtrl%solarVarOption == 1 .or. &
                  rtCtrl%solarVarOption == 2 .or. &
                  rtCtrl%solarVarOption == 3 ) ) then
         STOP '--- '//routineName//'(): Invalid value for solar variability option flag.'//' Program stopped.'
      endif
      if ( (rtCtrl%solarVarOption==2 .or. rtCtrl%solarVarOption==3) .and. rtCtrl%solarSource /=2 ) then
         STOP '--- '//routineName//'(): Temporal variability needs 3-component solar irradiance data.'//' Program stopped.'
      endif
      !if ( .NOT.( rtCtrl%solCycFrac>=0. .and. rtCtrl%solCycFrac<=1. ) ) then
      !   STOP '--- '//routineName//'(): Invalid value for solar cycle fraction.'//' Program stopped.'
      !endif
      if ( rtCtrl%solCycFrac >=0. .and. (rtCtrl%solarVarOption/=2 .or. rtCtrl%solarSource/=2) ) then
         STOP '--- '//routineName//'(): Solar cycle fraction applies only if solar source ="NRL" and option =2.'//' Program stopped.'
      endif

      !--- Check fluxCtrl ---------------------------------------------
      !
      if (fluxCtrl%flux_flag .eqv. .True.) then
         if ( .NOT.( fluxCtrl%nang == 1 .or. &    !if =1 one quadrature angle
                  fluxCtrl%nang == 2 .or. &    !if =2 two quadrature angles
                  fluxCtrl%nang == 3 ) ) then  !if =3 three quadrature angles
         STOP '--- '//routineName//'(): Invalid value for number of quadrature angles.'//' Program stopped.'
         endif
      endif


      !--- Check postCtrl ----------------------------------------------
      !
      if ( (postCtrl%functID <0 .or. postCtrl%functID>14) .and. postCtrl%functID/=FILLINT ) then
         STOP '--- '//routineName//'(): Function ID value has to be >=0 and <=14.'//' Program stopped.'
      endif

      if ( .not.postCtrl%FFT .and. postCtrl%functID >6 ) then
         STOP '--- '//routineName//'(): Invalid function ID for non-FFT scan.'//' Program stopped.'
      endif

      if ( ( postCtrl%functID ==12 ) .AND. &
           .NOT.( postCtrl%functParams(1) >0. .and. &
                  postCtrl%functParams(1) <1. ) ) then
         STOP '--- '//routineName//'(): Invalid parameter value for function ID 12, need to be in [0,1].'//' Program stopped.'
      endif

      if ( postCtrl%functID == 13 .AND. &
           .NOT.( postCtrl%functParams(1) >2. .and. &
                  postCtrl%functParams(1) <4. ) ) then
         STOP '--- '//routineName//'(): Invalid parameter value for function ID 13, need to be in [2,4].'//' Program stopped.'
      endif

      !if ( postCtrl%functID ==14 .AND. &
      !     all( abs(postCtrl%functParams(:)) <tiny(0.) ) ) then
      !   STOP '--- '//routineName//'(): fftscan_JFN=13 but all fftscan_PARMS(:) = 0.'//' Program stopped.'
      !endif



      !--- Check i/o files ---------------------------------------------
      !
      inquire( FILE=ioFiles%sceneFile, EXIST=fileExist )
      if (.NOT.fileExist) STOP '--- '//routineName//'(): Scene file is not exist.'//' Program stopped.'

      inquire( FILE=ioFiles%lineFile, EXIST=fileExist )
      if (.NOT.fileExist) STOP '--- '//routineName//'(): Line data file not exist.'//' Program stopped.'

      if ( NLTE ) then
         inquire( FILE=ioFiles%nlteStatPopFile, EXIST=fileExist )
         if ( .NOT. fileExist ) then
            STOP '--- '//routineName//'(): State populations as a function of altitude required whne NLTE==1.'//' Program stopped.'
         endif
      endif

      if ( .not.odCtrl%XsectOff ) then
         inquire( FILE=trim(ioFiles%xsFilePath)//'FSCDXS', EXIST=fileExist )
         if ( .NOT. fileExist ) then
            STOP '--- '//routineName//'(): .not.XSecOff, "FSCDXS" file must be existed under xsFilePath.'//' Program stopped.'
         endif
      endif

      if ( postCtrl%functID==0 ) then !Instrument response function filtering
         inquire( FILE=trim(ioFiles%filterFunctPath)//trim(ioFiles%filterFunctFile), EXIST=fileExist )
         if (.not. fileExist) then
            STOP '--- '//routineName//'(): Filter function file not exist.'//' Program stopped.'
         endif
      endif



      !--- Checks involve multiple groups
      !
      if ( (outCtrl%OD==0 .and. outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%Jac==0) .and. &
           (ioFiles%inFile_Rad=='' .and. &
            ioFiles%inFile_TxTot=='' .and. &
            ioFiles%inFile_TxPrfl=='' .and. &
            ioFiles%inFile_Jac=='' ) ) then
         STOP '--- '//routineName//'(): No calculation requested, nothing to do.'//' Program stopped.'
      endif

      if ( (outCtrl%Rad<0 .or. outCtrl%Tx<0 .or. outCtrl%Jac<0) .and. postCtrl%functID<0 ) then
         STOP '--- '//routineName//'(): Request post processing but no convolution function provided.'//' Program stopped.'
      endif

      if ( (outCtrl%Rad==0 .and. outCtrl%Jac==0) .and. & !OD-only or Tx-only mode
           (rtCtrl%ThermalOn .or. rtCtrl%SolarOn) ) then
         STOP '--- '//routineName//'(): OD-only or Tx-only mode needs the thermal source and solar source to be OFF.'
      endif

      if ( postCtrl%functID>=0 .and. (outGrid%DV<0.)) then
         STOP '--- '//routineName//'(): Post convolution requested but no output DV specified.'
      endif

      if ( (outCtrl%Rad/=0 .or. outCtrl%Jac/=0) .and. ( .not.rtCtrl%ThermalOn .and. .not.rtCtrl%SolarOn) ) then
         STOP '--- '//routineName//'(): Rad or Jac requested but there are no emitting sources.'//' Program stopped.'
      endif

      if ( rtCtrl%SolarOn ) then
         if ( rtCtrl%solarSource ==1 ) then
            inquire( FILE=trim(ioFiles%solarPath)//trim(ioFiles%solarFile_Kurucz), EXIST=fileExist )
            if (.not. fileExist) then
               STOP '--- '//routineName//'(): Solar irradiance file (Kurucz data) is missing.'//' Program stopped.'
            endif
         endif

         if ( rtCtrl%solarSource ==2 .and. rtCtrl%solarVarOption==1 ) then
            inquire( FILE=trim(ioFiles%solarPath)//trim(ioFiles%solarFile_1comp_NRL), EXIST=fileExist )
            if (.not. fileExist) then
               STOP '--- '//routineName//'(): Solar irradiance file (1 component NRL data) is missing.'//' Program stopped.'
            endif
         endif

         if ( rtCtrl%solarSource ==2 .and. rtCtrl%solarVarOption==2 .or. rtCtrl%solarVarOption==3 ) then
            inquire( FILE=trim(ioFiles%solarPath)//trim(ioFiles%solarFile_3comp_NRL), EXIST=fileExist )
            if (.not. fileExist) then
               STOP '--- '//routineName//'(): Solar irradiance file (3 component NRL data) is missing.'//' Program stopped.'
            endif
         endif
      endif


   END SUBROUTINE


END MODULE
