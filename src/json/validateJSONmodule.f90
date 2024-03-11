!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: yhe@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2016-, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

module validateJSONmodule
   use JsonModule, only: value_t, typeJson, typeString, typeLogic, typeArray, typeNull, typeFloat, typeInt, &
                  maxStringLength, compareStrings

   use json_data_types, only : CLBLM_IN_TRIGGER, SPECTRAL_CONVOLUTION_FLAGS_TRIGGER, OUTPUT_GRID_TRIGGER,&
                  RT_FLAGS_TRIGGER,  solar_variability_trigger, PATH_TRIGGER, OD_FLAGS_TRIGGER, &
                  SCENE_SELECTION_TRIGGER, GEOMETRY_TRIGGER, CLBLM_OUT_TRIGGER, NLTE_TRIGGER, &
                  TARGET_VIEWING_TRIGGER, FLUX_FLAGS_TRIGGER, &
                  SOLARSELECTION_IN_TRIGGER
   implicit none
   public :: isJSONvalid
   public :: isSOLARSELECTIONValid
   private
   ! form json_data_types
   ! character(len=*), parameter :: CLBLM_IN_TRIGGER                   = "clblm-in"
   ! character(len=*), parameter :: SPECTRAL_CONVOLUTION_FLAGS_TRIGGER = "spectral-convolution-flags"
   ! character(len=*), parameter :: OUTPUT_GRID_TRIGGER                = "output-spectral-grid"
   ! character(len=*), parameter :: RT_FLAGS_TRIGGER                   = "rt-flags"

   ! character(len=*), parameter :: solar_variability_trigger          = "solar-irradiance"
   ! character(len=*), parameter :: PATH_TRIGGER                       = "path-calculation-ctrl"
   ! character(len=*), parameter :: OD_FLAGS_TRIGGER                   = "od-flags"

   ! character(len=*), parameter :: SCENE_SELECTION_TRIGGER            = "scenes"
   ! character(len=*), parameter :: GEOMETRY_TRIGGER                   = "geometry"
   ! character(len=*), parameter :: CLBLM_OUT_TRIGGER                  = "clblm-out"
   ! character(len=*), parameter :: NLTE_TRIGGER                       = "nlte"

   ! character(len=*), parameter :: TARGET_VIEWING_TRIGGER             = "target-viewing"

   logical, parameter :: dbg = .false.
   logical, parameter :: caseSensitiveValidation=.true.

contains    
logical function isJSONvalid(this)
    use JsonModule, only: JsonList_t
   type(JsonList_t), intent(in) :: this
   integer   jj
   character(len=maxStringLength) :: bufStr
   isJSONvalid = .TRUE.
   if (dbg) print *, '>>>>>>>>>>>>>>>> JSON LIST'
   do jj = 1, size(this%jsnObj)
      if (dbg) print *, 'validating: ', trim(this%jsnObj(jj)%name)
      bufStr = trim(this%jsnObj(jj)%name)
      if (compareStrings(bufStr, CLBLM_IN_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_CLBLM_IN(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, SPECTRAL_CONVOLUTION_FLAGS_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_SPECTRAL_CONVOLUTION_FLAGS(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, OUTPUT_GRID_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_OUTPUT_GRID_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, RT_FLAGS_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_OUTPUT_RT_FLAGS_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, Flux_FLAGS_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_OUTPUT_Flux_FLAGS_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, solar_variability_trigger)) then
         isJSONvalid = isJSONvalid .and. validate_solar_variability_trigger(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, PATH_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_PATH_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, OD_FLAGS_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_OD_FLAGS_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, SCENE_SELECTION_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_SCENE_SELECTION_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, GEOMETRY_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_GEOMETRY_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, CLBLM_OUT_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_CLBLM_OUT_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      if (compareStrings(bufStr, NLTE_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_NLTE_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif

      if (compareStrings(bufStr, TARGET_VIEWING_TRIGGER)) then
         isJSONvalid = isJSONvalid .and. validate_TARGET_VIEWING_TRIGGER(this%jsnObj(jj)%value)         
         cycle
      endif
      print *,'ERR: validateJSONmodule::isJSONvalid - root key unknown : ',  trim(this%jsnObj(jj)%name)
      isJSONvalid = .false.
   end do
   if (dbg) print *, '<<<<<<<<<<<<<<<< JSON LIST'
end function isJSONvalid

! Checks the given JsonList_t object to see if it is
! a valid Solar Selection configuration. Detailed notes
! on the schema can be found in the documentation for
! getConfigFileStruct_SOLARSELECTION() in ./JsonConfig.f90
!
logical function isSOLARSELECTIONValid(this)
    use JsonModule, only: JsonList_t
   type(JsonList_t), intent(in) :: this
   integer   jj
   logical hasInputs, hasOutput
   character(len=maxStringLength) :: bufStr
   isSOLARSELECTIONValid = .TRUE.
   hasInputs = .FALSE.
   hasOutput = .FALSE.
   if (dbg) print *, '>>>>>>>>>>>>>>>> JSON LIST'
   do jj = 1, size(this%jsnObj)
       if (dbg) print *, 'validating: ', trim(this%jsnObj(jj)%name)
       bufStr = trim(this%jsnObj(jj)%name)
       if (compareStrings(bufStr, SOLARSELECTION_IN_TRIGGER)) then
           hasInputs = validate_SOLARSELECTION_IN_TRIGGER(this%jsnObj(jj)%value)
           cycle
       endif
       if (compareStrings(bufStr, "path")) then
           hasOutput = .TRUE.
           cycle
       endif
       print *,'ERR: validateJSONmodule::isSOLARSELECTIONValid - root key unknown : '//trim(this%jsnObj(jj)%name)
       isSOLARSELECTIONValid = .FALSE.
   end do
   if (.not. hasInputs ) then
       print *, 'ERR: validateJSONmodule::isSOLARSELECTIONValid - missing or invalid "'//SOLARSELECTION_IN_TRIGGER//'" section'
       isSOLARSELECTIONValid = .FALSE.
   endif
   if (.not. hasOutput ) then
       print *, 'ERR: validateJSONmodule::isSOLARSELECTIONValid - missing "path" value for output filename'
       isSOLARSELECTIONValid = .FALSE.
   endif

   if (dbg) print *, '<<<<<<<<<<<<<<<< JSON LIST'
end function isSOLARSELECTIONValid

!
! "SOLARSELECTION-in"
! Primary keys:  list of "path", "start-wavenumber", "end-wavenumber" objects
logical function validate_SOLARSELECTION_IN_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj, ii
   logical hasPath, hasStart, hasEnd
   real(8) :: v1, v2
   character(len=maxStringLength) :: bufStr

   if (dbg) print *, 'validating: validate_SOLARSELECTION_IN_TRIGGER'
   validate_SOLARSELECTION_IN_TRIGGER = .true.
   if ( dat%elementType /= typeArray ) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 SOLARSELECTION_IN_TRIGGER ,' key. must be Array. use [...] to enclose'
         validate_SOLARSELECTION_IN_TRIGGER = .FALSE.
         return
   endif
   do jj = 1, size(dat%arrDat)
      if (dbg) print *, 'validating input index: ', jj
      if ( dat%arrDat(jj)%elementType /= typeJson ) then
          print *, 'ERR: validateJSONmodule - wrong type of ', &
                  SOLARSELECTION_IN_TRIGGER, '[', jj, ']. Must be an Object. Use {...} to enclose'
          validate_SOLARSELECTION_IN_TRIGGER = .FALSE.
          return
      endif

      ! TODO: check that start < end, and maybe even no overlapping selections?
      hasPath = .FALSE.
      hasStart = .FALSE.
      hasEnd = .FALSE.
      do ii = 1, size(dat%arrDat(jj)%jsnLst%jsnObj)
          ! Make sure there are "path", "start-wavenumber", and "end-wavenumber"
          ! and only these.
          bufStr = trim( adjustl(dat%arrDat(jj)%jsnLst%jsnObj(ii)%name) )
          if (compareStrings(bufStr, 'path')) then
              hasPath = .TRUE.
              cycle
          endif
          if (compareStrings(bufStr, 'start-wavenumber')) then
              hasStart = .TRUE.
              v1 = dat%arrDat(jj)%jsnLst%jsnObj(ii)%value%fltDat
              cycle
          endif
          if (compareStrings(bufStr, 'end-wavenumber')) then
              hasEnd = .TRUE.
              v2 = dat%arrDat(jj)%jsnLst%jsnObj(ii)%value%fltDat
              cycle
          endif
          print *,'ERR: validateJSONmodule - unknown key : ',  &
                 bufStr , ' for root key ', SOLARSELECTION_IN_TRIGGER, '[', ii, ']'
          validate_SOLARSELECTION_IN_TRIGGER = .FALSE.
      enddo
      if ( .not. (hasPath .and. hasStart .and. hasEnd) ) then
          print *, 'ERR: validateJSONmodule - missing key for root key ', &
              SOLARSELECTION_IN_TRIGGER
          validate_SOLARSELECTION_IN_TRIGGER = .FALSE.
      else if ( v1 .gt. v2 ) then
          print *, 'ERR: validateJSONmodule - start-wavenumber must be less-than-or-equal-to end-wavenumber'
          validate_SOLARSELECTION_IN_TRIGGER = .FALSE.
      endif
   enddo
   if (size(dat%arrDat) .eq. 0) then
       print *, 'ERR: validateJSONmodule - Need at least one input entry'
       validate_SOLARSELECTION_IN_TRIGGER = .FALSE.
   endif
   if (dbg) print *, 'validating: validate_SOLARSELECTION_IN_TRIGGER DONE'
end function validate_SOLARSELECTION_IN_TRIGGER


!
! "clblm-out"
! Primary keys:  "rad", "total-tx", " tx-profile", "jacobians"
! Comments: If preceded by ‘mono’: monochromatic data requested
!           If preceded by ‘convolved’: spectral convolution is applied
! Secondary keys: "jacobian-list"
!           Only used if "jacobians" key is included

logical function validate_CLBLM_OUT_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj, idx, sz
   character(len=maxStringLength) :: bufStr

   if (dbg) print *, 'validating: validate_CLBLM_OUT_TRIGGER'
   validate_CLBLM_OUT_TRIGGER = .true.
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 CLBLM_OUT_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_CLBLM_OUT_TRIGGER = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'od'                  )) cycle
      if (compareStrings(bufStr, 'mono rad'            )) cycle
      if (compareStrings(bufStr, 'mono total-tx'       )) cycle
      if (compareStrings(bufStr, 'mono tx-profile'     )) cycle
      if (compareStrings(bufStr, 'mono jacobians'      )) cycle
      if (compareStrings(bufStr, 'convolved rad'       )) cycle
      if (compareStrings(bufStr, 'convolved total-tx'  )) cycle
      if (compareStrings(bufStr, 'convolved tx-profile')) cycle
      if (compareStrings(bufStr, 'convolved jacobians' )) cycle
      select case( trim( adjustl(dat%jsnLst%jsnObj(jj)%name)  ) ) 
      case('od')
         continue
      case('mono rad')
         continue
      case('mono total-tx')
         continue
      case('mono tx-profile')
         continue
      case('mono jacobians')
         continue
      case('convolved rad')
         continue
      case('convolved total-tx')
         continue
      case('convolved tx-profile')
         continue
      case('convolved jacobians')
         continue

      case('jacobian-list')
         continue
      case default
         print *,'ERR: validateJSONmodule - key unknown : ',  &
                 trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', CLBLM_OUT_TRIGGER
         validate_CLBLM_OUT_TRIGGER = .false.
      end select   
   enddo   
   if (dbg) print *, 'validating: validate_CLBLM_OUT_TRIGGER DONE'
end function validate_CLBLM_OUT_TRIGGER

!=================================================================================================
! "output-spectral-grid"
! "from", "to", "DV"
! Not used if convolving with user-supplied instrument function
! "grid-type"
! (OD-only mode)
logical function validate_OUTPUT_GRID_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_OUTPUT_GRID_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_OUTPUT_GRID_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 OUTPUT_GRID_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_OUTPUT_GRID_TRIGGER = .false.
         return
   endif

   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'from'     )) cycle
      if (compareStrings(bufStr, 'to'       )) cycle
      if (compareStrings(bufStr, 'DV'       )) cycle
      if (compareStrings(bufStr, 'grid-type')) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', OUTPUT_GRID_TRIGGER
      validate_OUTPUT_GRID_TRIGGER = .false.
   enddo   
   if (dbg) print *, 'validating: validate_OUTPUT_GRID_TRIGGER DONE'
end function validate_OUTPUT_GRID_TRIGGER
!=================================================================================================
! "clblm-in"
! "rad", "total-tx", " tx-profile", "jacobians"
! N/A
!=================================================================================================
logical function validate_CLBLM_IN(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr

   validate_CLBLM_IN = .true.
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 CLBLM_IN_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_CLBLM_IN = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'od'           )) cycle
      if (compareStrings(bufStr, 'rad'          )) cycle
      if (compareStrings(bufStr, 'total-tx'     )) cycle
      if (compareStrings(bufStr, 'tx-profile'   )) cycle
      if (compareStrings(bufStr, 'jacobians'    )) cycle
      if (compareStrings(bufStr, 'jacobian-list')) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', CLBLM_IN_TRIGGER
      validate_CLBLM_IN = .false.
   enddo   
end function validate_CLBLM_IN

!=================================================================================================
! "geometry"
! "obs-altitudes", "view-angles"
! N/A
!=================================================================================================
logical function validate_GEOMETRY_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_GEOMETRY_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_GEOMETRY_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 OUTPUT_GRID_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_GEOMETRY_TRIGGER = .false.
         return
   endif

   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'obs-altitudes')) cycle
      if (compareStrings(bufStr, 'view-angles'  )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', GEOMETRY_TRIGGER
      validate_GEOMETRY_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_GEOMETRY_TRIGGER DONE'
end function validate_GEOMETRY_TRIGGER

!=================================================================================================
!"scenes"
! "scene-file", "nscenes", "scene-ID"
! N/A
!=================================================================================================
logical function validate_SCENE_SELECTION_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_SCENE_SELECTION_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_SCENE_SELECTION_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 SCENE_SELECTION_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_SCENE_SELECTION_TRIGGER = .false.
         return
   endif

   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'scene-file')) cycle
      if (compareStrings(bufStr, 'nscenes'   )) cycle
      if (compareStrings(bufStr, 'scene-ID'  )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', SCENE_SELECTION_TRIGGER
      validate_SCENE_SELECTION_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_SCENE_SELECTION_TRIGGER DONE'
end function validate_SCENE_SELECTION_TRIGGER

!=================================================================================================
! "od-flags"
! "lines-contribution", "continuum-contribution", "line-rejection"
! "dptmin", dptfac"
! ignored if lines contribution or line rejection =. false.
! "p-convolution", "collision-partners-broadening"
! ignored if lines contribution = .false.
!=================================================================================================
logical function validate_OD_FLAGS_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_OD_FLAGS_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_OD_FLAGS_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 OD_FLAGS_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_OD_FLAGS_TRIGGER = .false.
         return
   endif

   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'lines-contribution'           )) cycle
      if (compareStrings(bufStr, 'continuum-contribution'       )) cycle
      if (compareStrings(bufStr, 'collision-partners-broadening')) cycle
      if (compareStrings(bufStr, 'line-rejection'               )) cycle
      if (compareStrings(bufStr, 'line-rejection-params'        )) cycle
      if (compareStrings(bufStr, 'speed-dependent-voigt'        )) cycle
      if (compareStrings(bufStr, 'x-sections-p-convolution'     )) cycle
      if (compareStrings(bufStr, 'continuum-scaling'            )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', OD_FLAGS_TRIGGER
      validate_OD_FLAGS_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_OD_FLAGS_TRIGGER DONE'
end function validate_OD_FLAGS_TRIGGER

!=================================================================================================
! "spectral-convolution-flags" (only used if one or more output RT product is convolved
! "FFT", "function ID", "function-params", "HWHM", "averaging-width"
! "boxcar-width"
! (if included, activates pre-averaging of monochromatic RT product)
! "filter-file"
! Must be included if convolving with user-supplied instrument function
!=================================================================================================
logical function validate_SPECTRAL_CONVOLUTION_FLAGS(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_SPECTRAL_CONVOLUTION_FLAGS = .true.
   if (dbg) print *, 'validating: validate_SPECTRAL_CONVOLUTION_FLAGS'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 SPECTRAL_CONVOLUTION_FLAGS_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_SPECTRAL_CONVOLUTION_FLAGS = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'FFT'            )) cycle
      if (compareStrings(bufStr, 'function ID'    )) cycle
      if (compareStrings(bufStr, 'function-params')) cycle
      if (compareStrings(bufStr, 'HWHM'           )) cycle
      if (compareStrings(bufStr, 'averaging-width')) cycle
      if (compareStrings(bufStr, 'filter-file'    )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', SPECTRAL_CONVOLUTION_FLAGS_TRIGGER
      validate_SPECTRAL_CONVOLUTION_FLAGS = .false.
   enddo      
   if (dbg) print *, 'validating: validate_SPECTRAL_CONVOLUTION_FLAGS DONE'   
end function validate_SPECTRAL_CONVOLUTION_FLAGS

!=================================================================================================
! "rt-flags"
! "thermal-source", "linear-in-tau", "solar-source", "solar-cnst", "julday" 
! "linear-in-tau"
!=================================================================================================
logical function validate_OUTPUT_RT_FLAGS_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_OUTPUT_RT_FLAGS_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_OUTPUT_RT_FLAGS_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 RT_FLAGS_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_OUTPUT_RT_FLAGS_TRIGGER = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'thermal-source')) cycle
      if (compareStrings(bufStr, 'linear-in-tau' )) cycle
      if (compareStrings(bufStr, 'solar-source'  )) cycle
      if (compareStrings(bufStr, 'solar-cnst'    )) cycle
      if (compareStrings(bufStr, 'julday    '    )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', RT_FLAGS_TRIGGER
      validate_OUTPUT_RT_FLAGS_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_OUTPUT_RT_FLAGS_TRIGGER DONE'   
end function validate_OUTPUT_RT_FLAGS_TRIGGER

!=================================================================================================
! "flux-flags"
! "flux_flag", "dv_flux", "nang"
!=================================================================================================
logical function validate_OUTPUT_Flux_FLAGS_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_OUTPUT_Flux_FLAGS_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_OUTPUT_Flux_FLAGS_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 Flux_FLAGS_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_OUTPUT_Flux_FLAGS_TRIGGER = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'flux_flag')) cycle
      if (compareStrings(bufStr, 'dv_flux' )) cycle
      if (compareStrings(bufStr, 'nang'  )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', Flux_FLAGS_TRIGGER
      validate_OUTPUT_Flux_FLAGS_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_OUTPUT_Flux_FLAGS_TRIGGER DONE'   
end function validate_OUTPUT_Flux_FLAGS_TRIGGER

!=================================================================================================
! "solar-irradiance"
! "option"
! Ignored if use solar irradiance data from Kuruscz
! "cycle-frac", "facula-var", "spot-var"
! Ignored if option = 1
!=================================================================================================
logical function validate_solar_variability_trigger(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_solar_variability_trigger = .true.
   if (dbg) print *, 'validating: validate_solar_variability_trigger'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 solar_variability_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_solar_variability_trigger = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'option'    )) cycle
      if (compareStrings(bufStr, 'cycle-frac')) cycle
      if (compareStrings(bufStr, 'facula-var')) cycle
      if (compareStrings(bufStr, 'spot-var'  )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', solar_variability_TRIGGER
      validate_solar_variability_trigger = .false.
   enddo      
   if (dbg) print *, 'validating: validate_solar_variability_trigger DONE'   
end function validate_solar_variability_trigger

!=================================================================================================
! "path-calculation-ctrl"
! "RT-grid", "airmass-scaling", "reference-path", "v-refrac"
! N/A
!=================================================================================================
logical function validate_PATH_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_PATH_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_PATH_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 PATH_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_PATH_TRIGGER = .false.
         return
   endif
   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'v-refrac'       )) cycle
      if (compareStrings(bufStr, 'airmass-scaling')) cycle
      if (compareStrings(bufStr, 'reference-path' )) cycle
      if (compareStrings(bufStr, 'RT-grid'        )) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', PATH_TRIGGER
      validate_PATH_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_PATH_TRIGGER DONE'     
end function validate_PATH_TRIGGER

!=================================================================================================
! "nlte" (single key)
! N/A
! N/A
!=================================================================================================
logical function validate_NLTE_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   validate_NLTE_TRIGGER = .true.
end function validate_NLTE_TRIGGER

!=================================================================================================
! Nothing in the carrent document
    ! 'target-alt'
    ! 'obs-alt'
    ! 'view-angle'
    ! 'horizontal'
!=================================================================================================
logical function validate_TARGET_VIEWING_TRIGGER(dat)
   type(value_t), intent(in)  :: dat
   integer jj
   character(len=maxStringLength) :: bufStr
   validate_TARGET_VIEWING_TRIGGER = .true.
   if (dbg) print *, 'validating: validate_TARGET_VIEWING_TRIGGER'
   if ( dat%elementType /= typeJson) then
         print *,'ERR: validateJSONmodule - wrong type of ', &
                 TARGET_VIEWING_TRIGGER ,' key. must be JSON. use {....} to enclose'
         validate_TARGET_VIEWING_TRIGGER = .false.
         return
   endif

   do jj = 1, size(dat%jsnLst%jsnObj)
      if (dbg) print *, 'validating: ', trim(dat%jsnLst%jsnObj(jj)%name)
      ! split the word
      bufStr = trim( adjustl(dat%jsnLst%jsnObj(jj)%name) )
      if (compareStrings(bufStr, 'target-alt')) cycle
      if (compareStrings(bufStr, 'obs-alt'   )) cycle
      if (compareStrings(bufStr, 'view-angle')) cycle
      if (compareStrings(bufStr, 'horizontal')) cycle
      print *,'ERR: validateJSONmodule - key unknown : ',  &
              trim(dat%jsnLst%jsnObj(jj)%name) , ' for root key ', TARGET_VIEWING_TRIGGER
      validate_TARGET_VIEWING_TRIGGER = .false.
   enddo      
   if (dbg) print *, 'validating: validate_TARGET_VIEWING_TRIGGER DONE'   
end function validate_TARGET_VIEWING_TRIGGER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function checkAllElements(dat, elType)
   type(value_t), intent(in)  :: dat
   integer      , intent(in)  :: elType
   integer                    :: jj
   checkAllElements = .true.
   do jj=1,size(dat%arrDat)
      if (dat%arrDat(jj)%elementType /= elType) checkAllElements = .false.
   end do
end function checkAllElements
end module validateJSONmodule
