module jsonconfig
use json_data_types
use stringmethods
use jsonmodule, ONLY : jsonlist_t, getJsonValue, printjsonlist, clearJsonList, value_t
use jsonfile,   only : readFile
!use printjson,  only : printJsonConfig, MAX_STRING_LENGTH

! For the not-so-fortran-savy reader, this set of interfaces
! can be thought of as like C/C++ function overloading.
! You call `getConfigFileStruct()`, and the compiler picks the
! right procedure to use based upon the arguments given.
interface getConfigFileStruct
    module procedure getConfigFileStruct_ClblmIn
    module procedure getConfigFileStruct_ClblmOut
    module procedure getConfigFileStruct_RtFlags
    module procedure getConfigFileStruct_FluxFlags
    module procedure getConfigFileStruct_SolarVariability
    module procedure getConfigFileStruct_Path
    module procedure getConfigFileStruct_OdFlags
    module procedure getConfigFileStruct_SpectralConvolutionFlags
    module procedure getConfigFileStruct_OutputGrid
    module procedure getConfigFileStruct_TargetViewing
    module procedure getConfigFileStruct_SceneSelection
    module procedure getConfigFileStruct_Geometry
    module procedure getConfigFileStruct_NLTE
    module procedure getConfigFileStruct_SolarSelection
    module procedure getConfigFileStruct_All
end interface

interface assignArray
    module procedure assignArray_Integer, assignArray_Float, & !assignArray_Double, &
            assignArray_Logical, assignArray_String
end interface

interface assignScalar
    module procedure assignInteger, assignLogical, assignFloat, &
            assignString
end interface

! TODO's
interface copyStruct
    module procedure copySceneSelection
    module procedure copyGeometry
    ! TODO
end interface


! TODO's
interface validateJson
!     !module procedure validate_DataPath
!     module procedure validate_ClblmIn
    module procedure validate_ClblmOut
!     module procedure validate_RtFlags
!     !module procedure validate_Mode
!     module procedure validate_Path
!     module procedure validate_OdFlags
!     !module procedure validate_VerticalGrid
!     module procedure validate_SpectralConvolutionFlags
!     module procedure validate_OutputGrid
!     module procedure validate_TargetViewing
!     module procedure validate_SceneSelection
!     module procedure validate_Geometry
!     !module procedure validate_All
end interface

contains

subroutine copyGeometry(source, destination)
    type(geometry_type), intent(in) :: source
    type(geometry_type), intent(out) :: destination

    if (allocated(destination%obs_altitudes)) deallocate(destination%obs_altitudes)
    if (allocated(destination%view_angles)) deallocate(destination%view_angles)

    if (allocated(source%obs_altitudes)) then
        allocate(destination%obs_altitudes(size(source%obs_altitudes)), stat=iret)
        if (iret .ne. 0) stop "cannot allocate obs altitudes; should not happen"
        destination%obs_altitudes(:) = source%obs_altitudes(:)
    endif
    if (allocated(source%view_angles)) then
        allocate(destination%view_angles(size(source%view_angles)), stat=iret)
        if (iret .ne. 0) stop "cannot allocate view angles; should not happen"
        destination%view_angles(:) = source%view_angles(:)
    endif

end subroutine

subroutine copySceneSelection(source, destination)
    implicit none

    type(scene_selection_type), intent(in) :: source
    type(scene_selection_type), intent(out) :: destination
    integer :: kret

    !destination%scene_file = repeat(' ', len(destination%scene_file))//trim(&
    !        source%scene_file)
    destination%scene_file = trim(source%scene_file)
    destination%nscenes = source%nscenes

    if (allocated(source%scene_id)) then
        if (.not. allocated(destination%scene_id)) then
            allocate(destination%scene_id(size(source%scene_id)), stat=kret)
            if (kret .ne. 0) stop "Unable to allocate SCENE ID's, this should not happen"
        else
            if (size(destination%scene_id) .ne. size(source%scene_id)) then
                deallocate(destination%scene_id)
                allocate(destination%scene_id(size(source%scene_id)), stat=kret)
                if (kret .ne. 0) stop "Unable to allocate SCENE ID's, this should not happen"
            endif
        endif
        destination%scene_id(:) = source%scene_id(:)
    else
        if (allocated(destination%scene_id)) deallocate(destination%scene_id)
    endif
    return
end subroutine



! TODO
logical function validate_ClblmIn(struct, errStr)
    implicit none
    character(len=*), intent(out) :: errstr
    type(clblm_in_type), intent(in) :: struct
    
    errStr = ''
    validate_ClblmIn = .false.
end function

! TODO
logical function validate_ClblmOut(struct, errStr)
    implicit none
    character(len=*), intent(out) :: errstr
    type(clblm_out_type), intent(in) :: struct

    errStr = ''
    validate_ClblmOut = .true.
end function

! Function overrides for scalar data assignments
subroutine assignInteger(jsonlist, keys, destination, errStr)
    implicit none
    integer :: source
    integer, intent(out) :: destination
    type(jsonlist_t), intent(in) :: jsonlist
    character(len=*), intent(out) :: errstr
    character(len=maxNameLength) :: keys(:)

    errStr = repeat(' ', len(errStr))
    if (getJsonValue(jsonlist, source, keys)) then
        destination = source
    else
        errStr = "Unable to parse or find '"//trim(keys(ubound(keys,1)))//"' key"
    endif
end subroutine

!------>>> ytma    
subroutine assignLogical(jsonlist, keys, logicDestination, intDestination, errStr)
    implicit none
    type(jsonlist_t) , intent(in)  :: jsonlist
    character(len=maxNameLength)   :: keys(:)
    logical                        :: logicDestination
    integer          , intent(out) :: intDestination
    character(len=*) , intent(out) :: errstr
    
    logical :: logicSource
    integer :: intSource

    errStr = repeat(' ', len(errStr))
    if (getJsonValue(jsonlist, logicSource,intSource, keys)) then
        intDestination = intSource
    else
        errStr = "Unable to parse or find '"//trim(keys( ubound(keys,1) ))//"' key"
    endif
end subroutine
!------>>><<<    

subroutine assignFloat(jsonlist, keys, destination, errStr)
    implicit none
    real :: source
    real, intent(out) :: destination
    type(jsonlist_t), intent(in) :: jsonlist
    character(len=*), intent(out) :: errstr
    character(len=maxNameLength) :: keys(:)

    errStr = repeat(' ', len(errStr))
    if (getJsonValue(jsonlist, source, keys)) then
        destination = source
    else
        errStr = "Unable to parse or find '"//trim(keys(ubound(keys,1)))//"' key"
    endif
end subroutine
subroutine assignString(jsonlist, keys, destination, errStr)
    implicit none
    character(len=MAX_STRING_LENGTH) :: source
    character(len=*), intent(out) :: destination
    type(jsonlist_t), intent(in) :: jsonlist
    character(len=*), intent(out) :: errstr
    character(len=maxNameLength) :: keys(:)

    errStr = repeat(' ', len(errStr))
    destination = repeat(' ', len(destination))
    source = repeat(' ', len(source))
    if (getJsonValue(jsonlist, source, keys)) then
        destination = trim(source)
    else
        errStr = "Unable to parse or find '"//trim(keys(ubound(keys,1)))//"' key"
    endif
end subroutine

! Overrides for Array assignment
subroutine assignArray_String(destination, source)
    character(len=*), intent(in) :: source(:)
    character(len=*), intent(out) :: destination(:)
    do i = 1, size(destination)
        destination(i) = repeat(' ', len(destination(i)))
        if (i .gt. size(source)) exit
        ilength = len(trim(source(i)))
        destination(i)(1:ilength) = source(i)(1:ilength)
    enddo
end subroutine
subroutine assignArray_Integer(destination, source)
    integer, intent(out) :: destination(:)
    integer, intent(in) :: source(:)

    do i = 1, size(destination)
        if (i .gt. size(source)) exit
        destination(i) = source(i)
    enddo
end subroutine
subroutine assignArray_Float(destination, source)
    real, intent(out) :: destination(:)
    real, intent(in) :: source(:)

    do i = 1, size(destination)
        if (i .gt. size(source)) exit
        destination(i) = source(i)
    enddo
end subroutine
subroutine assignArray_Double(destination, source)
    real*8, intent(out) :: destination(:)
    real*8, intent(in) :: source(:)

    do i = 1, size(destination)
        if (i .gt. size(source)) exit
        destination(i) = source(i)
    enddo
end subroutine

!------>>> yma
subroutine assignArray_Logical(destination, logicSource,intSource)
    integer, intent(out) :: destination(:)
    logical, intent(in)  :: logicSource(:)
    integer, intent(in)  :: intSource(:)

    do i = 1, size(destination)
        if (i .gt. size(intSource)) exit
        destination(i) = intSource(i)
    enddo
end subroutine
!------>>><<<

logical function getConfigFileStruct_Geometry(jsonlist, config, errStr)
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(geometry_Type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    character(len=MAX_STRING_LENGTH) :: charVal
    integer :: kret

    real, allocatable :: tempArr(:)

    charVal = repeat(' ', len(charVal))
    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_Geometry = .false.
    keys(1) = GEOMETRY_TRIGGER

    keys(2) = 'obs-altitudes'
    if (getJsonValue(jsonlist, tempArr, keys)) then
        allocate(config%obs_altitudes(size(tempArr)), stat=kret)
        if (kret .ne. 0) then
            errStr = 'Unable to allocate scene_selection%scene_id'
            return
        endif
        call assignArray(config%obs_altitudes, tempArr)
    endif
    if (allocated(tempArr)) deallocate(tempArr)

    keys(2) = 'view-angles'
    if (getJsonValue(jsonlist, tempArr, keys)) then
        allocate(config%view_angles(size(tempArr)), stat=kret)
        if (kret .ne. 0) then
            errStr = 'Unable to allocate scene_selection%scene_id'
            return
        endif
        call assignArray(config%view_angles, tempArr)
    endif
    if (allocated(tempArr)) deallocate(tempArr)

    getConfigFileStruct_Geometry = .true.
    return
end function

logical function getConfigFileStruct_SceneSelection(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(scene_selection_type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    character(len=MAX_STRING_LENGTH) :: charVal
    integer :: kret
    integer, allocatable :: tempArr(:)

    charVal = repeat(' ', len(charVal))
    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_SceneSelection = .false.
    keys(1) = returnString(SCENE_SELECTION_TRIGGER)

    keys(2) = returnString('scene-file')
    call assignScalar(jsonlist, keys, config%scene_file, errstr)

    keys(2) = returnString('nscenes')
    call assignScalar(jsonlist, keys, config%nscenes, errStr)

    keys(2) = returnString('scene-ID')
    if (getJsonValue(jsonlist, tempArr, keys)) then
        allocate(config%scene_id(size(tempArr)), stat=kret)
        if (kret .ne. 0) then
            errStr = 'Unable to allocate scene_selection%scene_id'
            return
        endif
        call assignArray(config%scene_id, tempArr)
        deallocate(tempArr)
    endif

    getConfigFileStruct_SceneSelection = .true.
    return
end function

logical function getConfigFileStruct_TargetViewing(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(target_viewing_type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_TargetViewing = .false.
    keys(1) = returnString(TARGET_VIEWING_TRIGGER)

    keys(2) = returnString('target-alt')
    call assignScalar(jsonlist, keys, config%target_alt, errStr)

    keys(2) = returnString('obs-alt')
    call assignScalar(jsonlist, keys, config%obs_alt, errStr)

    keys(2) = returnString('view-angle')
    call assignScalar(jsonlist, keys, config%view_angle, errStr)

    keys(2) = returnString('horizontal')
    call assignScalar(jsonlist, keys, config%horizontal, errStr)

    getConfigFileStruct_TargetViewing = .true.
    return
end function

logical function getConfigFileStruct_OutputGrid(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(output_spectral_grid_type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    character(len=MAX_STRING_LENGTH) :: charVal

    charVal = repeat(' ', len(charVal))
    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_OutputGrid = .false.
    keys(1) = returnString(OUTPUT_GRID_TRIGGER)

    keys(2) = returnString('from')
    call assignScalar(jsonlist, keys, config%v1, errStr)

    keys(2) = returnString('to')
    call assignScalar(jsonlist, keys, config%v2, errStr)

    keys(2) = returnString('DV')
    call assignScalar(jsonlist, keys, config%dv, errStr)

    keys(2) = returnString('grid-type')
    call assignScalar(jsonlist, keys, config%grid_type, errStr)

    getConfigFileStruct_OutputGrid = .true.
    return
end function

logical function getConfigFileStruct_SpectralConvolutionFlags(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(spectral_convolution_flags_type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    character(len=MAX_STRING_LENGTH) :: charVal
    real, allocatable :: tempFloatArray(:)
    logical :: logicDummy
    
    charVal = repeat(' ', len(charVal))
    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_SpectralConvolutionFlags = .false.
    keys(1) = returnString(SPECTRAL_CONVOLUTION_FLAGS_TRIGGER)

    keys(2) = returnString('FFT')
    call assignScalar(jsonlist, keys, logicDummy,config%fft, errStr)

    keys(2) = returnString('function ID')
    call assignScalar(jsonlist, keys, config%function_id, errStr)

    keys(2) = returnString('function-params')
    if (getJsonValue(jsonlist, tempFloatArray, keys)) then
        call assignArray(config%function_params, tempFloatArray)
        deallocate(tempFloatArray)
    else
        errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif

    keys(2) = returnString('HWHM')
    call assignScalar(jsonlist, keys, config%hwhm, errStr)

    keys(2) = returnString('averaging-width')
    call assignScalar(jsonlist, keys, config%averaging_width, errStr)

    keys(2) = returnString('filter-file')
    call assignScalar(jsonlist, keys, config%filter_file, errStr)

    getConfigFileStruct_SpectralConvolutionFlags = .true.
    return
end function

logical function getConfigFileStruct_OdFlags(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(od_Flags_type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    real, allocatable :: tempFloatArray(:)
    logical :: logicDummy
    integer :: intDummy

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_OdFlags = .false.
    keys(1) = returnString(OD_FLAGS_TRIGGER)

    keys(2) = returnString('lines-contribution')
    call assignScalar(jsonlist, keys, logicDummy, config%lines_contribution, errStr)

    keys(2) = returnString('continuum-contribution')
    call assignScalar(jsonlist, keys, logicDummy, config%continuum_contribution, errStr)

    keys(2) = returnString('collision-partners-broadening')
    call assignScalar(jsonlist, keys, logicDummy, config%collision_partners_broadening, errStr)

    keys(2) = returnString('speed-dependent-voigt')
    call assignScalar(jsonlist, keys, config%sdep_voigt_flag, errStr)
!------>
! IP
! activate reading of line-rejection status
    keys(2) = returnString('line-rejection')
    call assignScalar(jsonlist, keys, logicDummy,config%line_rejection, errStr)

    keys(2) = returnString('line-rejection-params')
    if (getJsonValue(jsonlist, tempFloatArray, keys)) then
        call assignArray(config%line_rejection_params, tempFloatArray)
        deallocate(tempFloatArray)
    else
        errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif


    keys(2) = returnString('x-sections-p-convolution')
    call assignScalar(jsonlist, keys, logicDummy,config%p_convolution, errStr)

    keys(2) = returnString('continuum-scaling')
    if (getJsonValue(jsonlist, tempFloatArray, keys)) then
        call assignArray(config%continuum_scaling, tempFloatArray)
        deallocate(tempFloatArray)
    else
        errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif
    getConfigFileStruct_OdFlags = .true.
end function

logical function getConfigFileStruct_Path(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(path_type), intent(inout) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    logical, allocatable :: logicTempArray(:)
    integer, allocatable :: intTempArray(:)

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_Path = .false.
    keys(1) = returnString(PATH_TRIGGER)

    keys(2) = returnString('v-refrac')
    call assignScalar(jsonlist, keys, config%v_refrac, errStr)

    keys(2) = returnString('airmass-scaling')
!------>>> yma
    if (getJsonValue(jsonlist, logicTempArray,intTempArray, keys)) then
        call assignArray(config%airmass_scaling, logicTempArray,intTempArray)
        deallocate(logicTempArray)
        deallocate(intTempArray)
    else
        errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif
    
    keys(2) = returnString('reference-path')
    call assignScalar(jsonlist, keys, config%reference_path, errStr)
    
    keys(2) = returnString('RT-grid')
    call assignScalar(jsonlist, keys, config%RT_grid, errStr)
!------>>><<<    
    getConfigFileStruct_Path = .true.
    return
end function

logical function getConfigFileStruct_RtFlags(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(rt_flags_type), intent(out) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    logical :: logicDummy

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_RtFlags = .false.
    keys(1) = returnString(RT_FLAGS_TRIGGER)

    keys(2) = returnString('thermal-source')
    call assignScalar(jsonlist, keys, logicDummy,config%thermal_source, errStr)

    keys(2) = returnString('linear-in-tau')
    call assignScalar(jsonlist, keys, config%linear_in_tau, errStr)

    keys(2) = returnString('solar-source')
    call assignScalar(jsonlist, keys, config%solar_source, errStr)

    keys(2) = returnString('solar-cnst')
    call assignScalar(jsonlist, keys, config%solar_cnst, errStr)

    keys(2) = returnString('julday')
    call assignScalar(jsonlist, keys, config%julday, errStr)

    getConfigFileStruct_RtFlags = .true.
    return
end function

logical function getConfigFileStruct_FluxFlags(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(flux_flags_type), intent(out) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    logical :: logicDummy

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_FluxFlags = .false.
    keys(1) = returnString(FLUX_FLAGS_TRIGGER)

    keys(2) = returnString('flux_flag')
    call assignScalar(jsonlist, keys, logicDummy,config%flux_flag, errStr)

    keys(2) = returnString('dv_flux')
    call assignScalar(jsonlist, keys, config%dv_flux, errStr)

    keys(2) = returnString('nang')
    call assignScalar(jsonlist, keys, config%nang, errStr)

    getConfigFileStruct_FluxFlags = .true.
    return
end function


logical function getConfigFileStruct_SolarVariability(jsonlist, config, errStr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t)             ,intent(in)    :: jsonlist
    type(solar_variability_type) ,intent(out)   :: config
    character(len=*)             ,intent(inout) :: errstr
    
    character(len=maxNameLength) :: keys(2)
    logical :: logicDummy

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_SolarVariability = .false.
    keys(1) = returnString(solar_variability_TRIGGER)

    keys(2) = returnString('option')
    call assignScalar(jsonlist, keys, config%option, errStr)

    keys(2) = returnString('cycle-frac')
    call assignScalar(jsonlist, keys, config%cycle_frac, errStr)

    keys(2) = returnString('facula-var')
    call assignScalar(jsonlist, keys, config%facula_var, errStr)

    keys(2) = returnString('spot-var')
    call assignScalar(jsonlist, keys, config%spot_var, errStr)

    getConfigFileStruct_SolarVariability = .true.
    return
end function


logical function getConfigFileStruct_ClblmOut(jsonlist, config, errstr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t)      , intent(in)    :: jsonlist
    type(clblm_out_type)  , intent(inout) :: config
    character(len=*)      , intent(inout) :: errstr
    
    character(len=maxNameLength)              :: keys(2)
    character(len=maxNameLength), allocatable :: tempCharray(:)
    character(len=MAX_STRING_LENGTH)          :: charVal
    character(len=9)                          :: solver_types(2)
    integer, allocatable                      :: tempArr(:)
    integer                                   :: kret, i

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_ClblmOut = .false.
    charVal = repeat(' ', len(charval))
    keys(1) =  returnString(CLBLM_OUT_TRIGGER)

    solver_types(1) =  returnString("convolved")
    solver_types(2) =  returnString("mono     ")

    
    ! This can apparently be a string or an array.
    !keys(2) = 'od'
    !if (getJsonValue(jsonlist, tempArr, keys)) then
    !    allocate(config%list_of_layer_numbers(size(tempArr)), stat=kret)
    !    if (kret .ne. 0) then
    !        errStr = 'Unable to allocate clblm_out_type%list_of_layer_numbers'
    !        return
    !    endif
    !    call assignArray(config%list_of_layer_numbers, tempArr)
    !    deallocate(tempArr)
    !    config%od = "(see layer numbers)"
    !else if (getJsonValue(jsonlist, charVal, keys)) then
    !    if (trim(charVal) .eq. 'all') then
    !        config%od = "all"
    !    else
    !        errStr = 'Invalid value for clblm_out_type%od; must be "all" or a list of layer numbers'
    !        return
    !    endif
    !else
    !    config%od = "none"
    !endif
    charVal = repeat(' ', len(charval))
    keys(2) =  returnString('od')
    if (getJsonValue(jsonlist, charval, keys)) then
      config%od%solver_type =  returnString('mono')
      config%od%path = trim(charval)
    else
      errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif

    
    do i = 1, 2
        charVal = repeat(' ', len(charval))
        keys(2) = trim(solver_types(i))// returnString(' rad')
        if (getJsonValue(jsonlist, charval, keys)) then
            config%rad%solver_type = trim(solver_types(i))
            config%rad%path = trim(charval)
        else
            errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
        endif
        charVal = repeat(' ', len(charval))
        keys(2) = trim(solver_types(i))// returnString(' total-tx')
        if (getJsonValue(jsonlist, charval, keys)) then
            config%total_tx%solver_type = trim(solver_types(i))
            config%total_tx%path = trim(charval)
        else
            errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
        endif
        charVal = repeat(' ', len(charval))
        keys(2) = trim(solver_types(i))// returnString(' tx-profile')
        if (getJsonValue(jsonlist, charval, keys)) then
            config%tx_profile%solver_type = trim(solver_types(i))
            config%tx_profile%path = trim(charval)
        else
            errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
        endif
        charVal = repeat(' ', len(charval))
        keys(2) = trim(solver_types(i))// returnString(' jacobians')
        if (getJsonValue(jsonlist, charval, keys)) then
            config%jacobians%solver_type = trim(solver_types(i))
            config%jacobians%path = trim(charval)
        else
            errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
        endif
    end do
    
    charVal = repeat(' ', len(charval))
    keys(2) =  returnString('jacobian-list')
    if (getJsonValue(jsonlist, tempCharray, keys)) then
        allocate(config%jacobian_list(size(tempCharray)), stat=kret)
        if (kret .ne. 0) then
            errStr = 'Unable to allocate jacobian-list'
            return
        endif
        call assignArray(config%jacobian_list, tempCharray)
        deallocate(tempCharray)
    else
        errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif

    getConfigFileStruct_ClblmOut = .true.
    return
end function

logical function getConfigFileStruct_ClblmIn(jsonlist, config, errstr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t), intent(in) :: jsonlist
    type(clblm_in_type), intent(out) :: config
    character(len=*), intent(inout) :: errstr
    character(len=maxNameLength) :: keys(2)
    character(len=MAX_STRING_LENGTH) :: charVal
    character(len=maxNameLength), allocatable :: tempCharray(:)
    integer :: kret

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_ClblmIn = .false.

    charVal = repeat(' ', len(charval))
    keys(1) = returnString(CLBLM_IN_TRIGGER)
    keys(2) = returnString('od')
    call assignScalar(jsonlist, keys, config%od_in, errstr)

    keys(2) = returnString('rad')
    call assignScalar(jsonlist, keys, config%rad_in, errstr)

    keys(2) = returnString('total-tx')
    call assignScalar(jsonlist, keys, config%total_tx_in, errstr)

    keys(2) = returnString('tx-profile')
    call assignScalar(jsonlist, keys, config%tx_profile_in, errstr)

    keys(2) = returnString('jacobians')
    call assignScalar(jsonlist, keys, config%jacobians_in, errstr)

    keys(2) = returnString('jacobian-list')
    if (getJsonValue(jsonlist, tempCharray, keys)) then
        allocate(config%jacobian_list(size(tempCharray)), stat=kret)
        if (kret .ne. 0) then
            errStr = 'Unable to allocate jacobian-list'
            return
        endif
        call assignArray(config%jacobian_list, tempCharray)
        deallocate(tempCharray)
    else
        errStr = "Unable to parse or find '"//trim(keys(2))//"' key"
    endif

    getConfigFileStruct_ClblmIn = .true.
    return
end function


logical function getConfigFileStruct_NLTE(jsonlist, config, errstr)
    use JsonModule, only : returnString
    implicit none
    type(jsonlist_t)     ,intent(in)    :: jsonlist
    type(nlte_flag_type) ,intent(out)   :: config
    character(len=*)     ,intent(inout) :: errstr
    
    character(len=maxNameLength) :: keys(1)
    logical :: logicDummy

    errStr = repeat(' ', len(errStr))
    getConfigFileStruct_NLTE = .false.
    
    keys(1) = returnString('nlte')
    call assignScalar(jsonlist, keys, logicDummy, config%nlte, errStr)

    getConfigFileStruct_NLTE = .true.
    return
end function


logical function validate_RtFlags(config, errstr)
    type(rt_flags_Type), intent(in) :: config
    character(len=*), intent(out) :: errstr
    validate_RtFlags = .false.
    if (config%linear_in_tau .lt. 0 .or. config%linear_in_tau .gt. 2) then
        errstr = "linear_in_tau should be 0, 1 or 2.  Your value is "//trim(&
                toString(config%linear_in_tau))//"."
        return
    endif
    validate_RtFlags = .true.
end function

logical function validate_FluxFlags(config, errstr)
    type(flux_flags_Type), intent(in) :: config
    character(len=*), intent(out) :: errstr
    validate_FluxFlags = .false.
    if (config%nang .lt. 1 .or. config%nang .gt. 3) then
        errstr = "nang should be 1, 2 or 3.  Your value is "//trim(&
                toString(config%nang))//"."
        return
    endif
    validate_FluxFlags = .true.
end function

! TODO: fill in debugging
logical function getConfigFileStruct_All(configFilePath, config, errstr)
    use validateJSONmodule, only: isJSONvalid
    implicit none
    character(len=*)              , intent(in)    :: configFilePath
    type(json_configuration_Type) , intent(inout) :: config
    character(len=1024)           , intent(inout) :: errstr

    type(jsonlist_t) :: jsonlist

    
    getConfigFileStruct_All = .false.
    errStr = repeat(' ', len(errStr))

    if ( readFile(trim(configFilePath), jsonlist)) then

        call printJsonList(jsonlist)
        
        if (.not. isJSONvalid(jsonList)) then
          print *, 'ERR: processing JSON input file', trim(configFilePath)
          call exit(1)
        endif

        if (.not. getConfigFileStruct(jsonlist, config%clblm_out, errStr)) then
            errStr = "CLBLM Error: Could not validateJson for CLBLM_OUT "//&
                        "group.  Error='"//trim(errstr)//"'."
             print *, 'getConfigFileStruct config%clblm_out failed'           
             goto 33
        endif

        ! Not required
        if (.not. getConfigFileStruct(jsonlist, config%clblm_in, errStr)) then
             print *, 'getConfigFileStruct config%clblm_in failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%rt_flags, errStr)) then
             print *, 'getConfigFileStruct config%rt_flags failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%flux_flags, errStr)) then
            print *, 'getConfigFileStruct config%flux_flags failed'           
           !goto 33
       endif
        if (.not. getConfigFileStruct(jsonlist, config%solar_variability, errStr)) then
             print *, 'getConfigFileStruct config%solar_variability failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%path, errStr)) then
             print *, 'getConfigFileStruct config%path failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%od_flags, errStr)) then
             print *, 'getConfigFileStruct config%od_flags failed'           
            !goto 33
        endif

        ! If at least one product is convolved, this is REQUIRED
        if (.not. getConfigFileStruct(jsonlist, config%spectral_convolution_flags, errStr)) then
            ! TODO, fill in logic
            !if ((trim(config%clblm_out%rad) == "convolved") .or. &
            !        (trim(config%clblm_out%total_tx) == "convolved") .or. &
            !        (trim(config%clblm_out%tx_profile) == "convolved") .or. &
            !        (trim(config%clblm_out%jacobians) == "convolved")) then
            !    errStr = "CLBLM Error: 'spectral_convolution_flags' JSON block is missing or invalid; "//&
            !        "REQUIRED with convolved products!  Error - '"//errStr(1:len(trim(errStr)))//"'."
            !    goto 33
            !endif
             print *, 'getConfigFileStruct config%spectral_convolution_flags failed'           
        endif

        ! If we have a v1 and v2, set the default v-refrac if it wasn't read in
        if (.not. getConfigFileStruct(jsonlist, config%output_grid, errStr)) then
             print *, 'getConfigFileStruct config%output_grid failed'           
            goto 33
        else
            if (config%path%V_refrac .le. -999) then
                config%path%v_refrac = 0.5*(config%output_grid%v1 + config%output_grid%v2)
            endif
        endif

        if (.not. getConfigFileStruct(jsonlist, config%target_viewing, errStr)) then
             print *, 'getConfigFileStruct config%target_viewing failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%scene_selection, errStr)) then
             print *, 'getConfigFileStruct config%scene_selection failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%geometry, errstr)) then
             print *, 'getConfigFileStruct config%geometry failed'           
            !goto 33
        endif
        if (.not. getConfigFileStruct(jsonlist, config%nlte_flag, errstr)) then
             print *, 'getConfigFileStruct config%nlte_flag failed'           
            !goto 33
        endif
        
    else
        errStr = 'CLBLM Error: Unable to parse JSON config file "'//trim(configFilePath)//'".'
        goto 33
    endif
    getConfigFileStruct_All = .true.    ! if you reached here, you parsed correctly.
 33 continue
    call clearJsonList(jsonlist)
    return
end function

! Load a JSON config file for constructing a NetCDF file
! made of selected wavenumber subsets of one or more other NetCDF files.
!
! Schema:
! {
!   "inputs": [{
!     "path": "input_file1.nc",
!     "start-wavenumber": 240.0,
!     "end-wavenumber": 1000.3
!   },{
!   ...
!   },{
!     "path": "input_fileN.nc",
!     "start-wavenumber": 2030.31,
!     "end-wavenumber": 406000.4
!   }],
!   "path": "output_file.nc"
! }
!
! Returns a .TRUE. for a successful load. The `config` object
! will have two members: 'output_file` (a string), and `input_selections`
! (an array of `solar_select_type`), which is guaranteed to be at
! at least of length one. The owner is responsible for deallocating
! this array when done.
!
! Each of these elements have three members: `input_file` (a string),
! `start_wavenumber` (a float), and `end_wavenumber` (also a float).
! `start_wavenumber` is guaranteed to be less-than-or-equal-to
! `end_wavenumber`, and it represents the inclusive bounds of the
! desired subset from the associated input file. While ill-advised,
! it is possible to specify overlapping bounds between two or more
! input files.
!
! It is perfectly valid to specify two or more subsets from the same file.
! Use an entry for each subset, each with the same "path" value.
!
! If the given config file fails to load, `errstr` will contain the
! error message.
!
logical function getConfigFileStruct_SolarSelection(configFilePath, config, errstr)
    use validateJSONmodule, only: isSolarSelectionValid
    use JsonModule, only : returnString
    implicit none
    character(len=*)               , intent(in)    :: configFilePath
    type(solar_configuration_Type) , intent(inout) :: config
    character(len=1024)            , intent(inout) :: errstr

    character(len=maxNameLength) :: keys(1)
    type(value_t) :: tempJson
    integer :: kret
    type(jsonlist_t) :: jsonlist
    integer :: i, j

    getConfigFileStruct_SolarSelection = .FALSE.
    errStr = repeat(' ', len(errStr))

    if ( readFile(trim(configFilePath), jsonlist)) then

        call printJsonList(jsonlist)

        if (.not. isSolarSelectionValid(jsonList)) then
            errStr = 'Invalid solar selection configuration: '//trim(configFilePath)
            return
        endif

        keys(1) = returnString("inputs")
        if (getJsonValue(jsonList, tempJson, keys)) then
            allocate(config%input_selections(size(tempJson%arrDat)), stat=kret)
            if (kret .ne. 0) then
                errStr = 'Unable to allocate solarbuild%input_selections'
                return
            endif

            do i = 1, size(config%input_selections)
                keys(1) = returnString("path")
                call assignScalar(tempJson%arrDat(i)%jsnLst, keys, config%input_selections(i)%input_file, errStr)
                keys(1) = returnString("start-wavenumber")
                call assignScalar(tempJson%arrDat(i)%jsnLst, keys, config%input_selections(i)%start_wavenumber, errStr)
                keys(1) = returnString("end-wavenumber")
                call assignScalar(tempJson%arrDat(i)%jsnLst, keys, config%input_selections(i)%end_wavenumber, errStr)
            enddo
        else
            errStr = "Unable to parse or find '"//trim(keys(1))//"' key"
        endif

        keys(1) = returnString("path")
        call assignScalar(jsonlist, keys, config%output_file, errStr)
    else
        errStr = 'SolarSelection Error: Unable to parse JSON config file "'//trim(configFilePath)//'".'
    endif
    getConfigFileStruct_SolarSelection = ( trim(errStr) .eq. "")
    call clearJsonList(jsonlist)
    return
end function getConfigFileStruct_SolarSelection


end module
