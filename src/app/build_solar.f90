      program solar_selection
      ! Subsets one or more model netcdf files into a single netcdf
      ! readable by clblm.
      ! The schema of the output netcdf file will match the schema
      ! of the input netcdf file(s).
      !
      ! This program takes one command-line argument, which is the name
      ! of the JSON configuration file. The config has the following
      ! schema:
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
      ! More detail about this config file can be found in the
      ! documentation for getConfigFileStruct_SolarSelection()
      ! in src/json/JsonConfig.f90.
      !
      ! Note that the wavenumber selection range is inclusive.
      ! The new file will have a wavenumber domain starting at the
      ! smallest "start-wavenumber" and ending at the largest
      ! "end-wavenumber". Any gaps in the selections will result in NANs
      ! for that portion of the output.
      use netcdf
      use jsonconfig
      use json_data_types
      ! Needed for NaNs
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Defined by command-line argument
      character (len=1024) :: configname

      character(len=8) :: ncomp_str
      integer :: idx, ncomp, vgnu_ncomp, sz, nout, out_idx
      integer :: ncid_in, ncid_out, start_varid, end_varid, &
                 delta_varid, sol_varid, sol_in_varid, &
                 dimids(2), comp_dimid, pt_dimid
      ! Indices for start/end wavenumbers in an input file
      integer, allocatable :: v1_idx(:), npts(:)

      real(real64) :: vgnu_min, vgnu_max, vgnu_delta, &
                      delta, v1, v2
      real(real32),allocatable :: sol(:,:)
      real(real32) :: NaNf

      character(len=1024) :: errstr
      type(solar_configuration_Type) :: config
      type(solar_select_type) :: selection

      NaNf = ieee_value(NaNf, ieee_quiet_nan)

      if (command_argument_count() /= 1) then
          write(*,*) "1 argument required: config"
          stop
      end if

      call get_command_argument(1, configname)

      if ( .not. getConfigFileStruct( &
                    trim(configname), config, errstr) ) then
          write(*,*) "ERROR: Invalid solar selection config file: ", &
                     trim(configname), ".\nReason: ", trim(errstr)
          stop
      endif

      ! Establish output domain from selections (we are guaranteed at
      ! least one input selection after validation).
      sz = size(config%input_selections)
      vgnu_min = config%input_selections(1)%start_wavenumber
      vgnu_max = config%input_selections(1)%end_wavenumber
      do idx=2,sz
          if (config%input_selections(idx)%start_wavenumber < vgnu_min) then
              vgnu_min = config%input_selections(idx)%start_wavenumber
          endif
          if (config%input_selections(idx)%end_wavenumber > vgnu_max) then
              vgnu_max = config%input_selections(idx)%end_wavenumber
          endif
      enddo

      allocate(v1_idx(sz))
      allocate(npts(sz))
      vgnu_delta = 0.0
      vgnu_ncomp = 0
      do idx=1,sz
          selection = config%input_selections(idx)
          call check( nf90_open(trim(selection%input_file), &
                                mode=NF90_NOWRITE, ncid=ncid_in) )
          call check( nf90_inq_varid(ncid_in, "StartWavenumber", &
                                     start_varid) )
          call check( nf90_inq_varid(ncid_in, "EndWavenumber", &
                                     end_varid) )
          call check( nf90_inq_varid(ncid_in, "SpectralIncrement", &
                                     delta_varid) )
          call check( nf90_inq_dimid(ncid_in, "numComponents", &
                                     comp_dimid) )
          call check( nf90_get_var(ncid_in, start_varid, v1) )
          call check( nf90_get_var(ncid_in, end_varid, v2) )
          call check( nf90_get_var(ncid_in, delta_varid, delta) )
          call check( nf90_inquire_dimension(ncid_in, comp_dimid, &
                                             len=ncomp) )
          if ( v2 .lt. v1 ) then
              write (*,*) "invalid input file.", selection%input_file, &
                    ". End wavenumber is less than the start wavenumber"
              call check( nf90_close(ncid_in) )
              stop
          endif
          if ( delta .le. 0.0 ) then
              write (*,*) "invalid input file.", selection%input_file, &
                    ". Spectral Incrment is less than or equal to zero."
              call check( nf90_close(ncid_in) )
              stop
          endif

          if ( (vgnu_delta /= delta) .and. (vgnu_delta /= 0.0) ) then
              write (*,*) "input spectral increments do not match"
              call check( nf90_close(ncid_in) )
              stop
          endif
          vgnu_delta = delta
          if ( (vgnu_ncomp /= ncomp) .and. (vgnu_ncomp /= 0) ) then
              write (*,*) "input component counts do not match"
              call check( nf90_close(ncid_in) )
              stop
          endif
          vgnu_ncomp = ncomp
          if ( v1 .gt. selection%start_wavenumber ) then
              write (*,*) "input start wavenumber (", v1, &
                        ") is greater than the selection start", &
                        " wavenumber (", selection%start_wavenumber, ")"
              call check( nf90_close(ncid_in) )
              stop
          endif
          if ( v2 .lt. selection%end_wavenumber ) then
              write (*,*) "input end wavenumber (", v2, &
                        ") is less than the selection end", &
                        " wavenumber (", selection%end_wavenumber, ")"
              call check( nf90_close(ncid_in) )
              stop
          endif

          v1_idx(idx) = nint((selection%start_wavenumber - v1) &
                             / delta) + 1
          npts(idx) = nint((selection%end_wavenumber &
                            - selection%start_wavenumber) / delta) + 1
          call check( nf90_close(ncid_in) )
      enddo

      nout = nint((vgnu_max - vgnu_min) / delta) + 1
      ! Convert integer to string
      write(ncomp_str,'(I0)') vgnu_ncomp

      ! Start prepping the creation of the output netcdf file
      call check( nf90_create(trim(config%output_file), &
                              OR(NF90_CLOBBER, NF90_NETCDF4), &
                              ncid_out) )
      call check( nf90_def_dim(ncid_out, "numComponents", vgnu_ncomp, &
                               comp_dimid) )
      call check( nf90_def_dim(ncid_out, "numPoints", nout, pt_dimid) )
      dimids = (/ pt_dimid, comp_dimid /)

      call check( nf90_def_var(ncid_out, "StartWavenumber", NF90_DOUBLE, &
                  varid=start_varid) )
      call check( nf90_def_var(ncid_out, "EndWavenumber", NF90_DOUBLE, &
                  varid=end_varid) )
      call check( nf90_def_var(ncid_out, "SpectralIncrement", NF90_FLOAT, &
                  varid=delta_varid) )
      call check( nf90_def_var(ncid_out, "SolarData", NF90_FLOAT, &
                  dimids, sol_varid) )
      call check( nf90_put_att(ncid_out, sol_varid, &
                               '_FillValue', NaNf) )
      call check( nf90_put_att(ncid_out, NF90_GLOBAL, "SolarDataInfor", &
           trim(ncomp_str)// &
           "-component solar irradiance data from NRLSSI2 model.") )

      call check( nf90_enddef(ncid_out) )

      call check( nf90_put_var(ncid_out, start_varid, vgnu_min) )
      call check( nf90_put_var(ncid_out, end_varid, vgnu_max) )
      call check( nf90_put_var(ncid_out, delta_varid, vgnu_delta) )

      ! For each selection, write out the data.
      do idx=1,sz
          selection = config%input_selections(idx)
          call check( nf90_open(trim(selection%input_file), &
                                mode=NF90_NOWRITE, ncid=ncid_in) )
          call check( nf90_inq_varid(ncid_in, "SolarData", &
                                     sol_in_varid) )

          allocate(sol(npts(idx), vgnu_ncomp))
          call check( nf90_get_var(ncid_in, sol_in_varid, sol, &
                                   (/ v1_idx(idx), 1 /), &
                                   (/ npts(idx), vgnu_ncomp /)) )
          out_idx = nint( (selection%start_wavenumber - vgnu_min) &
                          / vgnu_delta ) + 1
          call check( nf90_put_var(ncid_out, sol_varid, sol, &
                                   (/ out_idx, 1 /), &
                                   (/ npts(idx), vgnu_ncomp /)) )

          deallocate(sol)
          call check( nf90_close(ncid_in) )
      enddo
      deallocate(v1_idx)
      deallocate(npts)
      call check( nf90_close(ncid_out) )

      contains

      subroutine check(status)
       integer, intent (in) :: status

       if (status /= nf90_noerr) then 
         write(*,*) trim(nf90_strerror(status))
         stop "Stopped"
       end if
      end subroutine check
      end program
