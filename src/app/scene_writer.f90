!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Program scene_writer
!-----------------------------------------------------------------------
   USE NetCDF
   USE Module_Scn_TAPE5     ,ONLY: CLBLM_TAPE5, &
                                   readTAPE5, &
                                   separateXSectDataPresent, &
                                   isotplPresent, &
                                   userRTgridPresent
   USE Module_Scn_Geometry  ,ONLY: CLBLM_ScnGeom, &
                                   calcH1H2FromPress, &
                                   readSceneGeometry_tape5
   USE Module_Scn_Profile   ,ONLY: CLBLM_Profile, &
                                   CLBLM_LayerGrid, &
                                   readProfile_tape5, &
                                   read_XSect_Profile_tape5, &
                                   readRTgrid_tape5, &
                                   calcAltGridFromPress
   USE Module_Scn_Surface   ,ONLY: CLBLM_Surface, &
                                   readSurface_tape5
   USE Module_Scn_Config    ,ONLY: openLogFile, IPR
   USE Module_Scn_netCDF    ,ONLY: createSceneFile, &
                                   writeSceneHeader, &
                                   writeScene
   !USE SceneTools
   IMPLICIT NONE

   character(*), parameter :: routineName='scene_writer'

   character(256)         :: tape5file, nctcdfSceneFile
   character              :: YorN*1
   character              :: scnName*80, aNum*5
   integer                :: nArg, it5, st, numScn, scnNo, nMol
   integer(4)             :: ncID
   logical                :: xsIsPresent, fileExists
   type(CLBLM_TAPE5)      :: tape5
   type(CLBLM_ScnGeom)    :: scnGeom
   type(CLBLM_LayerGrid)  :: layGrid
   type(CLBLM_Profile)    :: prfl
   type(CLBLM_Profile)    :: xPrfl
   type(CLBLM_Surface)    :: surf
   character(len(prfl%molID)) ,allocatable :: molID(:)



   !--- Open the logfile
   CALL openLogfile( 'scene_writer.log' )


   !--- Get scence file name from command line arguments
   !
   nArg = command_argument_count()
   if (nArg==0) STOP '--- Usage: executable-name  TAPE5-file-name1 TAPE5-file-name2 ...  netCDF-file-name'

   CALL get_command_argument( number=nArg, value=nctcdfSceneFile, status=st)
   if (st/=0) STOP '--- Usage: executable-name  TAPE5-file-name1 TAPE5-file-name2 ...  netCDF-file-name'

   inquire(file=nctcdfSceneFile, exist=fileExists)
   if (fileExists) then
      write(*,'(A)',advance='no') 'Output NetCDF file exists, overwrite? (Y/N): '
      read(*,'(A)') YorN
      if (.not.(YorN=='Y' .or. YorN=='y')) STOP
   endif



   !--- Create NetCDF file for output
   CALL createSceneFile( nctcdfSceneFile, ncID, createUID=.TRUE. )


   !--- Loop over the input TAPE5's. Write them all in a single NetCDF file.
   DO it5 = 1,nArg-1

      CALL get_command_argument( number=it5, value=tape5file, status=st)
      if (st/=0) then
         STOP '--- Usage: executable-name  TAPE5-file-name TAPE5-file-name2 ...   netCDF-file-name'
      endif

      !--- Open the TAPE5 file and read the scene data into an internal structure.
      CALL readTAPE5( tape5, tape5file )

      !---First, get scene geometry information.
      !   "scnGeom%latitude", "scnGeom%earthRadius" and "scnGeom%HSPACE" are
      !   needed to process profile.  Once the profile processing is done,
      !   the completed profile will passed to geometry processing subroutine to
      !   do unit conversion for geometry variables.
      CALL readSceneGeometry_tape5( scnGeom, tape5 )

      !--- Input the profile from TAPE5 file.
      CALL readProfile_tape5( prfl, tape5 )

      ! If x-section data is input separately from the line molecules,
      ! read them in and merged the data into the "prfl" object.
      ! The separately input x-section profile may be on different vertical grid,
      ! interpolation may be needed to merge it into the line molecular profile.
      xsIsPresent = separateXSectDataPresent( tape5 )
      if ( xsIsPresent ) then
         CALL read_XSect_Profile_tape5( xPrfl, tape5 )
      endif


      !--- preprocssing profile
      CALL processProfile( prfl, scnGeom, xsIsPresent, xPrfl )


      !--- Read in user provided RT grid.
      ! Atmospheric layering is based on RT grid levels.
      if ( userRTgridPresent(tape5) ) then
         CALL readRTgrid_tape5( layGrid, tape5 )

         !--- If RT grid was given in pressure, interpolate PBND to ZBND
         CALL calcAltGridFromPress( layGrid, &              !in/out
                                    prfl, &                 !in
                                    scnGeom%latitude, &     !in
                                    scnGeom%earthRadius )   !in
      endif


      !--- If the input H1 and H2 in scene geometry structure are pressure values, covert them to Km.
      CALL calcH1H2FromPress( scnGeom, &  !in/out
                              prfl )       !in

      !--- Get surface information
      CALL readSurface_tape5( surf, tape5 )



      !--- Check to make sure all scene has same molecues
      if (it5==1) then
         if (allocated(molID)) deallocate(molID)
         allocate( molID(prfl%nMol) )
         nMol = prfl%nMol
         molID(:) = prfl%molID(:)
      else
         if ( prfl%nMol /= nMOl) then
            STOP '---'//routineName//'(): All scenes in a scene file must have same molecules.'
         endif
         !///to do: check the moleculse are the same
      endif


      !--- Write output file header
      numScn = nArg-1
      !if ( userRTgridPresent(tape5) .and. layGrid%P_unit/=0 ) then
      if ( userRTgridPresent(tape5) ) then
         CALL writeSceneHeader(ncID, numScn, molID, layGrid )
      else
         CALL writeSceneHeader(ncID, numScn, molID )
      endif

      !---Write scecne data
      scnNo = it5
      write(aNum, "(I5)") scnNo
      scnName = 'scene_'//trim(adjustl(aNum))
      CALL writeScene( ncID, prfl, surf, scnGeom, scnName, scnNo)

   ENDDO


CONTAINS !------------------ Internal subroutine -----------------------

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   subroutine processProfile( prfl, scnGeom, xsIsPresent, xPrfl )
   !-----------------------------------------------------------------------
      USE Module_Scn_Geometry  ,ONLY: CLBLM_ScnGeom
      USE Module_Scn_Profile   ,ONLY: CLBLM_Profile, &
                                      mergeStdAtmos, &
                                      CMPALT_profile, &
                                      setTopOfAtmos, &
                                      interpAltTempOnPressGrid_profile, &
                                      mergeStdAtmos_xs, &
                                      convertMolUnit, &
                                      interpMolProfile, &
                                      scaleMolProfile, &
                                      calcAltGridFromPress
      USE Module_UnitConv      ,ONLY: press2mb, &
                                      temp2Kelvin
      !USE SceneTools
      IMPLICIT NONE

      type(CLBLM_Profile) ,intent(inout) :: prfl
      type(CLBLM_ScnGeom) ,intent(in)    :: scnGeom
      logical             ,intent(in)    :: xsIsPresent
      type(CLBLM_Profile) ,intent(inout) :: xPrfl

      integer :: il, nLnMol
      real    :: newP, newT


      !--- Unit conversion for P and T
      !
      if ( prfl%P_unit /=10 ) then !if not in mb
         do il = 1,prfl%nLev

            CALL press2mb( newP, prfl%P(il), prfl%P_unit )
            prfl%P(il) = newP
         enddo
      endif

      if ( prfl%T_unit /=10 ) then !if not in Kelvin
         do il = 1,prfl%nLev

            CALL temp2Kelvin( newT, prfl%T(il), prfl%T_unit )
            prfl%T(il) = newT
         enddo
      endif


      !---
      ! If needed, fill up missing pressure, temperature or
      ! constituent concentration using one of the 6 AFGL
      ! standard atmospheres
      !
      nLnMol = prfl%nMol

      if ( any( [ prfl%Z_unit, prfl%P_unit, prfl%T_unit ] <=6 ) .or. &
           any( prfl%molUnit(1:nLnMol) <=6 ) ) then

         CALL mergeStdAtmos( prfl )
      endif


      !--- CONVERSION OF GENERIC UNITS TO DENSITIES FOR LBLRTM RUNS
      if ( any( prfl%molUnit(1:nLnMol) /=11 ) ) then !not in number density
         CALL convertMolUnit( prfl )
      endif


      !---
      ! If altitude not provided, calculate altitude using the hydrostatic equation
      ! If IMMAX<0 Z(1) must be present.
      if (prfl%Z_unit <=0) then
         CALL CMPALT_profile( prfl, scnGeom%latitude, scnGeom%earthRadius )
      endif


      !--- Find the top level (or the maximum number of profile levels) based on HSPACE.
      !    The top level is the highest level that has &
      !    Z(toaLev)<HSPACE+0.001. prfl%toaLev <= prfl%nLev
      CALL setTopOfAtmos( prfl, scnGeom%HSPACE )


      !---
      ! If x-section data is input separately from the line molecules,
      ! read them in and merged the data into the "prfl" object.
      ! The separately input x-section profile may be on different vertical grid,
      ! interpolation may be needed to merge it into the line molecular profile.
      if ( xsIsPresent ) then

         !--- If needed, calculate altitude on the pressure grid from prfl
         !    using hydrostatic equation and by interpolation
         !
         if (xPrfl%Z_unit <=0) then  !yma: what is Z_unit>0? no interpolation?
            CALL interpAltTempOnPressGrid_profile( xPrfl, &              !in/out
                                                   prfl, &               !in
                                                   scnGeom%latitude, &   !in
                                                   scnGeom%earthRadius ) !in
         endif

         !--- Fill in default x-section profiles
         CALL mergeStdAtmos_xs( xPrfl )


         !---
         ! Interpolate the x-section profiles onto the altitude grid from "prfl" object.
         ! Also it converts the volume mixing ratio to number density.
         CALL interpMolProfile( xPrfl, prfl )

         !--- Transfer the x-section profiles from xPrfl to prfl structure.
         CALL groupXsMolAndLnMol( prfl, xPrfl )

      endif !x-section profile processing


      !--- scaling profile if necessary
      !   'l',or 'L' or
      !    number '1'     scaling factor used directly to scale profile
      !   'c' or 'C'      column amount in mol/cm^2 units to which the profile is to be scaled
      !   'd' or 'D'      column amount in Dobson units to which the profile is to be scaled
      !   'm' or 'M'      volume mixing ratio (ppv) wrt dry air for the total column to which the profile will be scaled!
      !   'p' or 'P'      value of Precipitable Water Vapor (cm) to which the profile will be scaled (water vapor only)
      !   CALL scaleMolProfile( prfl, ['ccl4'], [2.187421511738269E+015/2], ['c'] )
      !CALL scaleMolProfile( prfl, ['ccl4'], [0.5], ['l'] )

   end subroutine


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine groupXsMolAndLnMol( prfl, xPrfl )
   !--------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: xsMolNum, XsName_LnName_areEqual
      USE MOdule_Scn_Profile ,ONLY: CLBLM_Profile, CLBLM_Profile_resizeNumMol

      type(CLBLM_Profile) ,intent(inout) :: prfl
      type(CLBLM_Profile) ,intent(in)    :: xPrfl

      integer :: nLev, nLnMol, nXsMol, nMol, ix,il,nx,nm
      integer :: xsNum( xPrfl%nMol )

      nLev   = prfl%nLev
      nLnMol = prfl%nMol
      nXsMol = xPrfl%nMol


      !--- Check if the xs molecules already exist in the line molecules
      nx = nXsMol
      do ix = 1,nXsMol
         xsNum(ix) = xsMolNum( xPrfl%molID(ix) )
         do il = 1,nLnMol
            if (XsName_LnName_areEqual( xPrfl%molID(ix), prfl%molID(il) )) then
               xsNum(ix) = 0
               nx = nx-1
               exit
            endif
         enddo
      enddo


      !---Resize the molecular arrays in prfl structrue
      if (nx>0) then

         nMol = nLnMol+nx
         CALL CLBLM_Profile_resizeNumMol( prfl, nMol )

         !--- Move the xs arrays to prfl structure
         nm = nLnMol
         do ix = 1,nXsMol
            if (xsNum(ix) /=0) then
               nm = nm+1
               prfl%molID(nm) = xPrfl%molID(ix)
               prfl%molUnit(nm) = xPrfl%molUnit(ix)
               prfl%Q(nm,1:nLev) = xPrfl%Q(ix,1:nLev)
            endif
         enddo
      endif

   end subroutine

END PROGRAM

