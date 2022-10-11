! DataTypes.f90
! Author: Christopher Brodowski, AER. Inc (cbrodows@aer.com)

module json_data_types
use jsonmodule, ONLY: maxNameLength

integer, parameter :: MAX_STRING_LENGTH      = 256
integer, parameter :: NUM_CONTINUUM_ELEMENTS = 7
real,    parameter :: FILLER_FLOAT           = -999.9
integer, parameter :: FILLER_INT             = -999

!character(len=*), parameter :: SCENE_FILE_DEFAULT = "scenes.nc"


type scene_selection_type
    !sequence
    character(len=MAX_STRING_LENGTH) :: scene_file = repeat(' ', MAX_STRING_LENGTH) ! default is scenes.nc
    integer                          :: nscenes    = FILLER_INT
    integer, allocatable             :: scene_id(:)
end type


type target_viewing_type
    !sequence
    real :: target_alt = FILLER_FLOAT
    real :: obs_alt    = FILLER_FLOAT
    real :: view_angle = FILLER_FLOAT
    real :: horizontal = FILLER_FLOAT
end type


type geometry_type
    !sequence
    real, allocatable :: obs_altitudes(:)
    real, allocatable :: view_angles(:)
end type


type clblm_out_solver_type
    !sequence
    character(len=maxNameLength) :: solver_type = repeat(' ', maxNameLength)    ! mono or convolved
    character(len=maxNameLength) :: path        = repeat(' ', maxNameLength)    ! filename or path
end type
type clblm_out_type
    !sequence

    !! "All", or a list of layer numbers
    !character(len=MAX_STRING_LENGTH)          :: od = repeat(' ', MAX_STRING_LENGTH)
    !! TODO: This will get populated depending on the value of od; if "all", we
    !!       need to add all the layers.  If an array, we populate the values in the
    !!       array with everything.
    !integer,                      allocatable :: list_of_layer_numbers(:)
    type(clblm_out_solver_type)               :: od
    type(clblm_out_solver_type)               :: rad
    type(clblm_out_solver_type)               :: total_tx
    type(clblm_out_solver_type)               :: tx_profile
    type(clblm_out_solver_type)               :: jacobians
    character(len=maxNameLength), allocatable :: jacobian_list(:)    ! TODO
end type


type clblm_in_type
    !sequence
    character(len=MAX_STRING_LENGTH)          :: od_in          = repeat(' ', MAX_STRING_LENGTH)
    character(len=MAX_STRING_LENGTH)          :: rad_in         = repeat(' ', MAX_STRING_LENGTH)
    character(len=MAX_STRING_LENGTH)          :: total_tx_in    = repeat(' ', MAX_STRING_LENGTH)
    character(len=MAX_STRING_LENGTH)          :: tx_profile_in  = repeat(' ', MAX_STRING_LENGTH)
    character(len=MAX_STRING_LENGTH)          :: jacobians_in   = repeat(' ', MAX_STRING_LENGTH)
    character(len=maxNameLength), allocatable :: jacobian_list(:)    ! TODO
end type


type rt_flags_type
    integer                          :: thermal_source = -1
    integer                          :: linear_in_tau  = FILLER_INT
    character(len=MAX_STRING_LENGTH) :: solar_source   = repeat(' ', MAX_STRING_LENGTH)
    real                             :: solar_cnst     = FILLER_FLOAT
    integer                          :: julDay         = FILLER_INT
end type


!------>>> yma
type solar_variability_type
   integer :: option     = FILLER_INT
   real    :: cycle_frac = FILLER_FLOAT
   real    :: facula_var = FILLER_FLOAT
   real    :: spot_var   = FILLER_FLOAT
end type
!------>>><<<


type spectral_convolution_flags_type
    !sequence
    integer                          :: fft                = -1
    integer                          :: function_id        = FILLER_INT     ! if 0, filter function read from file
    character(len=MAX_STRING_LENGTH) :: filter_file        = repeat(' ', MAX_STRING_LENGTH) ! required if id = 0
    real                             :: function_params(4) = FILLER_FLOAT
    real                             :: HWHM               = FILLER_FLOAT
    real                             :: averaging_width    = FILLER_FLOAT
end type


type output_spectral_grid_type
    !sequence
    real                         :: V1        = FILLER_FLOAT    ! from
    real                         :: V2        = FILLER_FLOAT   ! to
    real                         :: dv        = FILLER_FLOAT   ! resolution
    character(len=maxNameLength) :: grid_type = repeat(' ', maxNameLength) ! adjustedDV/exactDV/uniformDV
end type


!------>>> yma
type path_type
    !sequence
    real    :: V_refrac           = FILLER_FLOAT   ! defaults to 0.5*(v1+v2); the -999.9 triggers it to set the default
    integer :: airmass_scaling(2) = -1  !logical :: airmass_scaling(2) = .FALSE.
    integer :: reference_path     = FILLER_INT
    integer :: RT_grid            = FILLER_INT
end type
!------>>><<<


type od_flags_type
    !sequence
    integer :: lines_contribution            = -1
    integer :: continuum_contribution        = -1
    integer :: collision_partners_broadening = -1
    integer :: line_rejection                = -1
    real    :: line_rejection_params(2)      = FILLER_FLOAT
    !real    :: dptmin                        = FILLER_FLOAT
    !real    :: DPTFAC                        = FILLER_FLOAT
    integer :: p_convolution                 = -1
    ! Tables say this should be a 7-element array; JSON files have this as a dictionary. FIXME
    ! Keeping as an array.  It's easier that way.
    real :: continuum_scaling(NUM_CONTINUUM_ELEMENTS) = FILLER_FLOAT
end type


type nlte_flag_type
    integer :: nlte = -1
end type

type flux_flags_type
    integer                          :: flux_flag      = -1
    integer                          :: nang           = FILLER_INT
    real                             :: dv_flux        = FILLER_FLOAT   ! spectral sum interval for flux calculation
end type

type solar_select_type
    character(len=MAX_STRING_LENGTH) :: input_file         = repeat(' ', MAX_STRING_LENGTH)
    real                             :: start_wavenumber   = FILLER_FLOAT
    real                             :: end_wavenumber     = FILLER_FLOAT
end type

type solar_configuration_Type
    type(solar_select_type), allocatable  :: input_selections(:)
    character(len=MAX_STRING_LENGTH)      :: output_file   = repeat(' ', MAX_STRING_LENGTH)
end type

type json_configuration_Type    
    type(nlte_flag_type)                  :: nlte_flag     !A single key to turn on/off Non-LTE mode.
    type(clblm_in_type)                   :: clblm_in
    type(clblm_out_type)                  :: clblm_out
    type(rt_flags_type)                   :: rt_flags
    type(path_type)                       :: path
    type(od_flags_Type)                   :: od_flags
    type(spectral_convolution_flags_type) :: spectral_convolution_flags
    type(output_spectral_grid_type)       :: output_grid
    type(target_viewing_type)             :: target_viewing
    type(scene_selection_type)            :: scene_selection
    type(geometry_Type)                   :: geometry
    type(solar_variability_type)          :: solar_variability
    type(flux_flags_type)                 :: flux_flags
end type

character(len=*), parameter :: NLTE_TRIGGER                       = "nlte"
character(len=*), parameter :: CLBLM_IN_TRIGGER                   = "clblm-in"
character(len=*), parameter :: CLBLM_OUT_TRIGGER                  = "clblm-out"
character(len=*), parameter :: SPECTRAL_CONVOLUTION_FLAGS_TRIGGER = "spectral-convolution-flags"
character(len=*), parameter :: OUTPUT_GRID_TRIGGER                = "output-spectral-grid"
character(len=*), parameter :: TARGET_VIEWING_TRIGGER             = "target-viewing"
character(len=*), parameter :: SCENE_SELECTION_TRIGGER            = "scenes"
character(len=*), parameter :: PATH_TRIGGER                       = "path-calculation-ctrl"
character(len=*), parameter :: RT_FLAGS_TRIGGER                   = "rt-flags"
character(len=*), parameter :: OD_FLAGS_TRIGGER                   = "od-flags"
character(len=*), parameter :: GEOMETRY_TRIGGER                   = "geometry"
character(len=*), parameter :: solar_variability_trigger          = "solar-irradiance"
character(len=*), parameter :: FLUX_FLAGS_trigger                 = "flux-flags"
! This is for a separate config file, so it doesn't go into the set of main triggers
character(len=*), parameter :: SOLARSELECTION_IN_TRIGGER          = "inputs"

character(len=256) :: VALID_MAIN_TRIGGERS(13)
data VALID_MAIN_TRIGGERS /NLTE_TRIGGER, CLBLM_IN_TRIGGER, CLBLM_OUT_TRIGGER, SPECTRAL_CONVOLUTION_FLAGS_TRIGGER, &
        OUTPUT_GRID_TRIGGER, TARGET_VIEWING_TRIGGER, SCENE_SELECTION_TRIGGER, & !           = "scenes"
        PATH_TRIGGER, RT_FLAGS_TRIGGER, OD_FLAGS_TRIGGER, GEOMETRY_TRIGGER, & !                  = "geometry"
        solar_variability_trigger, & !         = "solar-irradiance"
        FLUX_FLAGS_trigger / !         = "flux-flags"
save VALID_MAIN_TRIGGERS

end module
