!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_Scn_netCDF
   IMPLICIT NONE   
   PRIVATE 
   PUBLIC :: createSceneFile, &
             writeSceneHeader, &
             writeScene
   
   
CONTAINS ! --------------------- Module procedures ---------------------

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine createSceneFile(fname, ncID, title, createUID, enddef)
   !--------------------------------------------------------------------
      USE NetCDF
   
      character(len=*) ,intent(in)           :: fname
      integer(4)       ,intent(out)          :: ncid
      character(len=*) ,intent(in), optional :: title
      logical          ,intent(in), optional :: createUID
      logical          ,intent(in), optional :: enddef
   
      integer :: datetime(8), UID
      logical :: enddefFlag
      
      enddefFlag = .TRUE.
      if (present(enddef)) enddefFlag = enddef
      
      call check( nf90_create(fname, NF90_NETCDF4, ncid) )
      
      if (present(title)) &
            call check( nf90_put_att(ncid, NF90_GLOBAL, 'title',title) )
            
      if (present(createUID)) then
         if (createUID) then
            call date_and_time( values=datetime )
            UID = 1e9+datetime(5)*1e7+datetime(6)*1e5+datetime(7)*1e3+datetime(8) ![1hhmmssuuu]
            call check( nf90_put_att(ncid, NF90_GLOBAL, 'UID',UID) )
         endif
      endif

      if (enddefFlag) call check( nf90_enddef(ncid) )

   end subroutine


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine writeSceneHeader(ncID, nScn, molID, layGrid )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_Scn_Profile   ,ONLY: CLBLM_LayerGrid
      
      integer(4)                     ,intent(in) :: ncid
      integer                        ,intent(in) :: nScn
      character(*)                   ,intent(in) :: molID(:)
      type(CLBLM_LayerGrid),optional ,intent(in) :: layGrid
      
      
      call check( nf90_redef(ncid))

      !--- Number of scenes
      call check( nf90_put_att(ncid, NF90_GLOBAL, "numScenes", nScn) )      
      
      !--- RT grid
      if (present(layGrid)) then
         call check( nf90_put_att( ncID, NF90_GLOBAL, 'num_RT_grid_lev', size(layGrid%Z)) )
         if (layGrid%P_unit == 10) then 
            call check( nf90_put_att( ncID, NF90_GLOBAL, 'RT_grid_lev_in_mb',layGrid%P) )      
         else
            call check( nf90_put_att( ncID, NF90_GLOBAL, 'RT_grid_lev_in_km', layGrid%Z) )
         endif
      endif
            
      !--- Molecular ID's
      call check( nf90_put_att( ncID, NF90_GLOBAL, 'numMol',size(molID) ))
      call writeHeaderAtt_molID(ncID, 'molecules', molID)
      
      
      call check( nf90_enddef(ncid) )
      
   end subroutine


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine writeHeaderAtt_molID(ncID, attName, molID)
   !--------------------------------------------------------------------
      USE NetCDF
      
      integer(4)    ,intent(in) :: ncid
      character(*)  ,intent(in) :: attName
      character(*)  ,intent(in) :: molID(:)

      integer :: im
      character :: molIDStr*1000, str*20
            
      molIDStr=''
      do im=1,size(molID)
         write(str,'(A)'), trim(adjustl(molID(im)))
         molIDStr = trim(molIDStr)//', '//str
      enddo
      molIDStr(1:1)='' !remove the first ","
      call check( nf90_put_att( ncid, NF90_GLOBAL, attName, trim(adjustl(molIDStr))) )      
      
   end subroutine


   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE writeScene( ncID, prfl, surf, scnGeom, scnName, scnNo )
   !--------------------------------------------------------------------     
      USE NetCDF
      USE Module_Scn_Profile   ,ONLY: CLBLM_Profile, &
                                      CLBLM_LayerGrid
      USE Module_Scn_Surface   ,ONLY: CLBLM_Surface
      USE Module_Scn_Geometry  ,ONLY: CLBLM_ScnGeom

      integer(4)           ,intent(in) :: ncid
      type(CLBLM_Profile)  ,intent(in) :: prfl
      type(CLBLM_Surface)  ,intent(in) :: surf
      type(CLBLM_ScnGeom)  ,intent(in) :: scnGeom      
      character(*)         ,intent(in) :: scnName
      integer              ,intent(in) :: scnNo
  
      character(*) ,parameter :: routineName='writeScene'
  
      integer    :: i,j,k
      integer(4) :: grpid
      integer(4) :: molDimId, levDimId, EmRfDimID, dimIds(2)
      integer(4) :: zId, pId, tId, qId, firstZId
      integer(4) :: obsAltId, obsPressId, obsZenId, tgtAltId, lenId, spaceId, earthRId, obsLatId
      integer(4) :: sfcTempId, sfcModeID
      integer    :: numEmisReflPts, numEmisPts, numReflPts
      integer(4) :: sfcSpectrGridID, sfcEmisID, sfcReflID 
      integer    :: sfcPropInputMode
      real ,allocatable :: emisNodeFreq(:), reflNodeFreq(:), sfcEmis(:), sfcRefl(:)

      
      !--- Surface thermal reflection mode
      if (surf%surf_refl=='l') then
         sfcPropInputMode = 2
      else
         sfcPropInputMode = 1
      endif

      ! Surface emittance and reflectivity in CLBLM are either input as look-up table 
      ! or as a spectrally independent constant, no quadratic formula used anymore. Only 
      ! the first coefficient from the TAPE5 are output to NetCDF scene file.      
      !
      if (surf%emisCoef(1)<0) then 
         numEmisPts = surf%emisTbl%NLIM
         allocate( emisNodeFreq( numEmisPts ) )
         allocate( sfcEmis( numEmisPts ) )
         do i=1,numEmisPts
            emisNodeFreq(i) = surf%emisTbl%V1 + (i-1)*surf%emisTbl%DV
         enddo
         sfcEmis(1:numEmisPts) = surf%emisTbl%ZEMIS( 1:numEmisPts )
      else
         numEmisPts = 1
         allocate( sfcEmis( numEmisPts ) )
         sfcEmis(1) = surf%emisCoef(1)
      endif
      
      if (surf%reflCoef(1)<0) then
         numReflPts = surf%reflTbl%NLIM         
         allocate( reflNodeFreq( numReflPts ) )
         allocate( sfcRefl( numReflPts ))
         do i=1,numReflPts
            reflNodeFreq(i) = surf%reflTbl%V1 + (i-1)*surf%reflTbl%DV
         enddo
         sfcRefl(1:numReflPts) = surf%reflTbl%ZRFLT( 1:numReflPts )
      else
         numReflPts = 1
         allocate( sfcRefl( numReflPts ))
         sfcRefl(1) = surf%reflCoef(1)
      endif                        
      
      if (numEmisPts>1 .or. numReflPts>1) then
         if ( numEmisPts /= numReflPts .or. any( emisNodeFreq /= reflNodeFreq) ) then
            STOP '---'//routineName//'(): In CLBLM, surface emissivity and reflectivity should always on same spectral grid.'
         endif
      endif
      numEmisReflPts = numEmisPts
      
      
      !--- Re-open the define mode
      call check( nf90_redef(ncid))
   
      !--- Create a group and write group attributes
      call check( nf90_def_grp(ncid, scnName, grpID ) )
      call check( nf90_put_att(grpID, NF90_GLOBAL, 'sceneNumber', scnNo) )       
      call writeHeaderAtt_molID(grpID, 'molecules', prfl%molID)
      
      
      !--- Define dimensions in the scene
      call check( nf90_def_dim( grpID, "numPrflLev"     ,int(prfl%nLev,4)      ,levDimId ) )
      call check( nf90_def_dim( grpId, "numMol"         ,int(prfl%nMol,4)      ,molDimId ) )
      call check( nf90_def_dim( grpId, "numEmisReflPts" ,int(numEmisReflPts,4) ,EmRfDimID ) )


      !--- Define variables: Prfofiles
      !
      call check(nf90_def_var(grpId, "altitude", NF90_DOUBLE, levDimId, zId ))
      call check( nf90_put_att(grpId, zId, "long_name", "Level altitude") )
      call check( nf90_put_att(grpId, zId, "units", "km") )
      !call check( nf90_put_att(grpId, zId, "valid_range",(/-1.0, 120.0/) ) )

      call check(nf90_def_var(grpId, "sfcAltitude", NF90_DOUBLE, firstZId ))
      call check( nf90_put_att(grpId, firstZId, "long_name", "Altitude for the first pressure level.") )
      call check( nf90_put_att(grpId, firstZId, "units", "km") )
      
      call check(nf90_def_var(grpId, "pressure", NF90_DOUBLE, levDimId, pId ))
      call check( nf90_put_att(grpId, pId, "long_name", "Level pressure") )
      call check( nf90_put_att(grpId, pId, "units", "mb") )
      call check( nf90_put_att(grpId, pId, "valid_range",(/0.0, 1200.0/) ) )

      call check(nf90_def_var(grpId, "temperature", NF90_DOUBLE, levDimId, tId ))
      call check( nf90_put_att(grpId, tId, "long_name", "Level temperature") )
      call check( nf90_put_att(grpId, tId, "units", "K") )
      !call check( nf90_put_att(grpId, tId, "valid_range",(/150.0, 400.0/) ) )
            
      dimids(1) = levDimId
      dimids(2) = molDimId
      call check(nf90_def_var(grpId, "molDensities", NF90_DOUBLE, dimids(1:2), qId ))
      call check( nf90_put_att(grpId, qId, "long_name", "Molecular concentration on levels") )
      call check( nf90_put_att(grpId, qId, "units", "[cm^{-3}]") )
      call check( nf90_put_att(grpId, qId, "valid_range",(/0., 1e26/) ) )
      
      
      !--- Define variables: Geometry
      !
      call check(nf90_def_var(grpId, "obsAltitude", NF90_DOUBLE, obsAltId ))
      call check( nf90_put_att(grpId, obsAltId, "long_name", "Observer level altitude") )
      call check( nf90_put_att(grpId, obsAltId, "units", "km") )
      !call check( nf90_put_att(grpId, obsAltId, "valid_range", (/0., /) ) )

      !call check(nf90_def_var(grpId, "obsPressure", NF90_DOUBLE, obsPressId ))
      !call check( nf90_put_att(grpId, obsPressId, "long_name", "Observer level pressure") )
      !call check( nf90_put_att(grpId, obsPressId, "units", "mb") )
      
      call check(nf90_def_var(grpId, "viewAngle", NF90_DOUBLE, obsZenID ))
      call check( nf90_put_att(grpId, obsZenID, "long_name", "Viewing zenith angle at observer altitude") )
      call check( nf90_put_att(grpId, obsZenID, "units", "deg") )
      call check( nf90_put_att(grpId, obsZenID, "valid_range", (/0., 180.0/) ) )

      !call check(nf90_def_var(grpId, "targetAlt", NF90_DOUBLE, tgtAltId))
      !call check( nf90_put_att(grpId, tgtAltId, "long_name", "Target altitude") )
      !call check( nf90_put_att(grpId, tgtAltId, "units", "km") )
      !!call check( nf90_put_att(grpId, tgtAltId, "valid_range", (/0., 120.0/) ) )

      !call check( nf90_def_var( grpID, "pathLengthType", nf90_int, lenId ) )
      !call check( nf90_put_att( grpID, lenId, 'long_name', 'Indicator of path length through a tangent height. if = 0, short path, if = 1, long path through a tangent height.') )
      !call check( nf90_put_att( grpID, lenId, "units", "none") )
      !call check( nf90_put_att( grpID, lenId, "valid_range",(/0, 1/) ) )
      !call check( nf90_put_att( grpId, lenId, "default", 0 ) )
      !call check( nf90_put_att( grpID, lenId, 'comment', 'pathLengthType is needed when H1>H2 and obsAngle>90.') )
      
      !call check(nf90_def_var(grpId, "spaceAlt", NF90_DOUBLE, spaceID))
      !call check( nf90_put_att(grpId, spaceID, "long_name", "Altitude definition for space") )
      !call check( nf90_put_att(grpId, spaceID, "units", "km") )
      !call check( nf90_put_att(grpId, spaceID, "valid_range", (/0., 120.0/) ) )
      !call check( nf90_put_att(grpId, spaceID, "default", 100.0 ) )

      call check(nf90_def_var(grpId, "earthRadius", NF90_DOUBLE, earthRId ))
      call check( nf90_put_att(grpId, earthRId, "long_name", "Earth radius") )
      call check( nf90_put_att(grpId, earthRId, "units", "km") )
      !call check( nf90_put_att(grpId, earthRId, "valid_range", (/0., 120.0/) ) )
      call check( nf90_put_att(grpId, earthRId, "default", 6371.23 ) )

      call check(nf90_def_var(grpId, "latitude", NF90_DOUBLE, obsLatID ))
      call check( nf90_put_att(grpId, obsLatID, "long_name", "Observer Latitude") )
      call check( nf90_put_att(grpId, obsLatID, "units", "deg") )
      call check( nf90_put_att(grpId, obsLatID, "valid_range", (/-90.0, 90.0/) ) )
      call check( nf90_put_att(grpId, obsLatID, "default", 45.0 ) )

      
      !--- Define variables: Surface
      !
      call check(nf90_def_var(grpId, "sfcSkinTemp", NF90_DOUBLE, sfcTempID ))
      call check( nf90_put_att(grpId, sfcTempID, "long_name", "Surface temperature") )
      call check( nf90_put_att(grpId, sfcTempID, "units", "K") )
      call check( nf90_put_att(grpId, sfcTempID, "valid_range", (/200.0, 400.0/) ) )
      
      call check(nf90_def_var(grpId, "sfcPropInputMode", NF90_INT, sfcModeID ))
      call check( nf90_put_att(grpId, sfcModeID, "long_name", "Surface property input mode" ))
      call check( nf90_put_att(grpId, sfcModeID, "units", "none") )
      call check( nf90_put_att(grpId, sfcModeID, "valid_range", (/1, 3/) ) )
      call check( nf90_put_att(grpId, sfcModeID, "default", 1 ) )
      call check( nf90_put_att(grpId, sfcModeID, "comment", &
         "1: thermal reflectivity = (1.- emissivity) with surface assumed specular"//&
         "2: thermal reflectivity = (1.- emissivity) with surface assumed Lambertian"//&
         "In modes 1 and 2, reflectivity of surface-incident solar beam in the direction of the observer is specified in separate sfcRefl array."))
      
      if (numEmisReflPts>1) then
         call check(nf90_def_var(grpId, "sfcPropSpectrGrid", NF90_DOUBLE, EmRfDimID, sfcSpectrGridID ))
         call check( nf90_put_att(grpId, sfcSpectrGridID, "long_name", "Surface optical property spectral hinge points") )
         call check( nf90_put_att(grpId, sfcSpectrGridID, "units", "cm^{-1}") )
         call check( nf90_put_att(grpId, sfcSpectrGridID, "valid_range", (/0.0, 50000.0/) ) )
      endif
      
      call check(nf90_def_var(grpId, "sfcEmis", NF90_DOUBLE, EmRfDimID, sfcEmisID ))
      call check( nf90_put_att(grpId, sfcEmisID, "long_name", "Surface emissivity at hinge points on sfcPropSpectrGrid if numEmisReflPts>1."//&
                                                                 "Or spectrally independent surface emissivty if numEmisReflPts==1.") )
      call check( nf90_put_att(grpId, sfcEmisID, "valid_range", (/0., 1.0/) ) )         
      
      call check(nf90_def_var(grpId, "sfcRefl", NF90_DOUBLE, EmRfDimID, sfcReflID ))
      call check( nf90_put_att(grpId, sfcReflID, "long_name", "Surface reflectivity at hinge points on sfcPropSpectrGrid if numEmisReflPts>1."//&
                                                                      "Or spectrally independent surface reflectivity if numEmisReflPts==1.") )
      call check( nf90_put_att(grpId, sfcReflID, "valid_range", (/0., 1.0/) ) )
      
   
      !--- End of define mode
      call check( nf90_enddef(ncid) )

      
      !--- Write data into variables
      !      
      call check( nf90_put_var(grpId, zId,  prfl%Z) )
      call check( nf90_put_var(grpId, pId,  prfl%P) )
      call check( nf90_put_var(grpId, tId,  prfl%T) )
      call check( nf90_put_var(grpId, qId,  transpose(prfl%Q) ) )![nMol,nLay] -> [nLay,nMol]
      call check( nf90_put_var(grpId, firstZId,  prfl%Z(1)) )

      call check( nf90_put_var(grpId, obsAltId,   scnGeom%Hobs) )
      !call check( nf90_put_var(grpId, obsPressId, -999) )
      call check( nf90_put_var(grpId, obsZenId,   scnGeom%obsAng) )
      !call check( nf90_put_var(grpId, tgtAltId,   scnGeom%Hend) )
      !call check( nf90_put_var(grpId, lenId,      scnGeom%LEN) )
      !call check( nf90_put_var(grpId, spaceID,    scnGeom%HSPACE) )
      call check( nf90_put_var(grpId, earthRId,   scnGeom%earthRadius) )
      call check( nf90_put_var(grpId, obsLatID,   scnGeom%latitude) )

      call check( nf90_put_var(grpId, sfcTempID, surf%Ts) )
      call check( nf90_put_var(grpId, sfcModeID, sfcPropInputMode) )

      if (numEmisReflPts>1) call check( nf90_put_var(grpId, sfcSpectrGridID, emisNodeFreq) )
      call check( nf90_put_var(grpId, sfcEmisID, sfcEmis ))
      call check( nf90_put_var(grpId, sfcReflId, sfcRefl ))
         
   END SUBROUTINE   

   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine closeFile(ncID)
   !--------------------------------------------------------------------
      USE NetCDF
      IMPLICIT NONE      
      integer(4), intent(in)  :: ncid
      call check( nf90_close(ncid) )
   end subroutine closeFile

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine check(status)
   !--------------------------------------------------------------------
      USE NetCDF
      IMPLICIT NONE      
      integer(4), intent (in) :: status

      if(status /= nf90_noerr) then
        print *, ' netcdfSceneFile error:', trim(nf90_strerror(status))
        call exit(1)
      end if
   end subroutine check

END MODULE
