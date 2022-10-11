!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_FileIO
   USE Module_ConstParam  ,ONLY: r8=>kind_r8
   USE Module_Config      ,ONLY: CLBLM_Post_Ctrl

   IMPLICIT NONE   
   PRIVATE
   PUBLIC :: createNetCDF4File, &
             openNetCDFFile, &
             closeNetCDFFile, &
             attIsPresent, &
             varIsPresent, &
             dimIsPresent, &
             dimByName, &
             readArray_char, &
             readScalar_real, &
             readScalar_int, &
             readArray_real_1D, &
             readArray_real_2D, &
             check
   PUBLIC :: checkNetcdfCall, &
             createNetCDFFile
   PUBLIC :: CLBLM_SpectFileHeader
   PUBLIC :: NWDL, endfil_4, SKIPFL


   
   TYPE :: CLBLM_SpectFileHeader
        character(256)        :: sceneFile
        integer               :: fileID       
        integer               :: sceneNo         !scene #
        character(32)         :: productName
        character(32)         :: productSpectType
        real(r8)              :: V1
        real(r8)              :: V2
        real                  :: DV
        integer               :: nSamp    !number of spectral points
        character(256)        :: filterFunctFile
        type(CLBLM_Post_Ctrl) :: convolParam
        real                  :: obsAlt
        real                  :: viewAng
   END TYPE
   

CONTAINS !=================== MODULE CONTAINS ==========================

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   subroutine createNetCDF4File(fname, ncID, title, createUID, enddef)
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
            call check( nf90_put_att(ncid, NF90_GLOBAL, 'fileID',UID) )
         endif
      endif

      if (enddefFlag) call check( nf90_enddef(ncid) )

   end subroutine createNetCDF4File

   
   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE openNetCDFFile( fname, ncID )
   !--------------------------------------------------------------------
      USE NetCDF
      character(len=*), intent(in)  :: fname
      integer(4),       intent(out) :: ncid

      call check( nf90_open(fname, nf90_nowrite, ncid) )
   END SUBROUTINE

   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   FUNCTION attIsPresent( grpId, attName )
   !--------------------------------------------------------------------
      USE NetCDF
      logical                       :: attIsPresent
      integer(4),        intent(in) :: grpId
      character(len=*),  intent(in) :: attName
      
      integer(4) :: status, attnum
   
      status =  nf90_inquire_attribute( grpId, NF90_GLOBAL, attName, attnum=attnum)
      if(status /= nf90_noerr) then
         attIsPresent = .FALSE.
      else
         attIsPresent = .TRUE.
      end if
   END FUNCTION

   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   FUNCTION varIsPresent( grpId, varName, varID )
   !--------------------------------------------------------------------
      USE NetCDF
      logical                       :: varIsPresent
      integer(4),        intent(in) :: grpId
      character(len=*),  intent(in) :: varName
      integer  ,optional,intent(out):: varID
      
      integer(4) :: status, vID
   
      status =  nf90_inq_varid( grpId, varName, varid=vID)
      if(status /= nf90_noerr) then
         varIsPresent = .FALSE.
      else
         varIsPresent = .TRUE.
         if (present(varID)) varID = vID
      end if
   END FUNCTION


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   FUNCTION dimIsPresent( grpId, dimName, dimID )
   !--------------------------------------------------------------------
      USE NetCDF
      logical                         :: dimIsPresent
      integer(4),          intent(in) :: grpId
      character(len=*),    intent(in) :: dimName
      integer(4) ,optional,intent(out):: dimID
      
      integer(4) :: status,dID
   
      status =  nf90_inq_dimid( grpId, dimName, dimid=dID)
      if(status /= nf90_noerr) then
         dimIsPresent = .FALSE.
      else
         dimIsPresent = .TRUE.
         if (present(dimID)) dimID=dID
      end if
   END FUNCTION


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   FUNCTION dimByName( grpId, dimName )
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4)                    :: dimByName
      integer(4),        intent(in) :: grpId
      character(len=*),  intent(in) :: dimName
      
      integer(4) :: dimID
      
      if ( dimIsPresent( grpID, dimName, dimID ) ) then
         call check( nf90_inquire_dimension( grpID, dimID, len=dimByName ) )         
      else
         dimByName = -1
      endif

   END FUNCTION
   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE  readArray_char( grpId, varName, val )
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4)       ,intent(in)    :: grpId
      character(len=*) ,intent(in)    :: varName
      character(len=*) ,intent(inout) :: val(:)
      
      integer :: status,  jj
      integer(4) :: varID
      
      call check( nf90_inq_varid( grpId, varName, varid=varId ) )         
      !call check( nf90_inquire_variable( grpId, varid, dimids=dimids ) )
      !!call check( nf90_inquire_dimension( grpId, dimids(1), len=lenMolName ) )
      !call check( nf90_inquire_dimension( grpId, dimids(2), len=nMol ) )                 
      !allocate( val(nMol) )         
      call check( nf90_get_var( grpId, varId, val ) )
      
      !do jj=1, tScene%nMol
      !   tscene%prfl%molID(jj) = trim( adjustl( val(jj) ) )
      !end do
      
   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE  readScalar_real( grpid, varName, val )
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4),           intent(in)   :: grpid
      character(len=*),     intent(in)   :: varName
      real,                 intent(inout):: val
               
      integer(4) :: varId

      call check( nf90_inq_varid( grpId, varName, varid=varId ) )
      call check( nf90_get_var( grpId, varId, val ) )
      
   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE  readScalar_int( grpid, varName, val )
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4),           intent(in)   :: grpid
      character(len=*),     intent(in)   :: varName
      integer,              intent(inout):: val

      integer(4) :: varId

      call check( nf90_inq_varid( grpId, varName, varid=varId ) )
      call check( nf90_get_var( grpId, varId, val ) )
      
   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE  readArray_real_1D( grpid, varName, val)
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4),        intent(in)   :: grpid
      character(len=*),  intent(in)   :: varName
      real,              intent(inout):: val(:)
      
      integer(4)  :: status, varId

      call check( nf90_inq_varid( grpId, varName, varid=varId ) )
      call check( nf90_get_var( grpId, varId, val ) )
      
   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE  readArray_real_2D( grpid, varName, val)
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4),       intent(in)   :: grpid
      character(len=*), intent(in)   :: varName
      real,             intent(inout):: val(:,:)
      
      integer(4) :: varId

      call check( nf90_inq_varid( grpId, varName, varid=varId ) )
      call check( nf90_get_var( grpId, varId, val ) )
      
   END SUBROUTINE
         
   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE closeNetCDFFile(ncID)
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4), intent(in)         :: ncid
        call check( nf90_close(ncid) )
   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE check(status)
   !--------------------------------------------------------------------
      USE NetCDF
      integer(4), intent (in) :: status

      if(status /= nf90_noerr) then
         print *, ' netcdfSceneFile error:', trim(nf90_strerror(status))
         call exit(1)
      end if
   END SUBROUTINE check
   


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   logical function createNetCDFFile(sFileName, iNcid)
   !--------------------------------------------------------------------
      USE NetCDF
      implicit none
      
      character(len=*) ,intent(in)  :: sFileName
      integer(4)       ,intent(out) :: iNcid

      createNetCDFFile = .false.
      call checkNetcdfCall(nf90_create(sFileName(1:len(trim(sFileName))), NF90_64BIT_OFFSET, iNcid), &
              sinStringCommand="nf90_create("//sFileName(1:len(trim(sFileName)))//")")
      createNetCDFFile = .true.
      return
    end function


   ! -----------------------------------------------------------------------
   ! Checks the return status of a NetCDF call.  If the call failed, logs the
   ! error message and (by default) kills the program outright.
   !
   ! @param {integer} iStat - The status code returned from a NetCDF call.
   ! @param {string} [sInStringCommand="UNKNOWN"] - The name of the NetCDF command being executed.
   ! @param {integer} [iLogUnit=6] - The file unit in which to log the error message.
   ! @param {logical} [bInToStop=.true.] - Whether or not to kill the program if error occurs.
   ! -----------------------------------------------------------------------
   subroutine checkNetcdfCall(iStat, sInStringCommand, iLogUnit, bInToStop)
   ! -----------------------------------------------------------------------
      USE NETCDF
      
      integer(4)       ,intent(in)           :: iStat
      character(len=*) ,intent(in) ,optional :: sInStringCommand
      integer          ,intent(in) ,optional :: iLogUnit
      logical          ,intent(in) ,optional :: bInToStop

      integer, parameter :: MAX_NC_BUFFER_LENGTH = 1024      
      integer :: iUnit
      logical :: toStop
      character(len=MAX_NC_BUFFER_LENGTH) :: sStringCommand

      sStringCommand = repeat(' ', len(sStringCommand))
      if (present(sInStringCommand)) then
          sStringCommand = sInStringCommand(1:len(trim(sInStringCommand)))
      else
          sStringCommand = 'UNKNOWN'
      endif

      toStop = .true.
      if (present(bInToStop)) toStop = bInToStop

      ! Warning: This does not check if the unit is open, conflicting,
      !          or anything like that.  That is the user's responsibility!
      iUnit = 6
      if (present(iLogUnit)) iUnit = iLogUnit

      if (iStat /= nf90_noerr) then
          write(iUnit, *) 'Error: NetCDF operation failed.  Error message: "', &
                  trim(nf90_strerror(iStat)), '" when attempting to run the NetCDF command "', &
                  sStringCommand(1:len(trim(sStringCommand))), '"".  Exiting with STOP 2...'
          if (toStop) STOP 2
      endif
    end subroutine







   



   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
      integer FUNCTION NWDL (IWD,ILAST) 
   !--------------------------------------------------------------------
      integer ,intent(in) :: IWD(*) !Con't be IWD(:), has to be the assumed-size array to allow the LBLRTM trick to count the size of a common block.
      integer ,intent(in) :: ILAST
      
      integer :: I

      DO 10 I = 1, 9000 
         IF (IWD(I).EQ.ILAST) THEN 
            NWDL = I-1 
            RETURN 
         ENDIF 
   10 END DO 

      STOP ' NWDL - IWD,ILAST ' 
      END FUNCTION

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
      SUBROUTINE ENDFIL (IFILE) 
   !--------------------------------------------------------------------
      integer , intent(in) :: IFILE
      
      !DIMENSION IDUM(6) 
      !DATA IDUM / 6*-99 / 
      integer ,PARAMETER :: IDUM(6)=-99

      CALL BUFOUT (IFILE,IDUM(1),6) 

      RETURN 
      END SUBROUTINE

      
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
      subroutine endfil_4 (ifile) 
   !--------------------------------------------------------------------
      integer*4 ,intent(in) :: ifile 
      
      !integer*4 idum(6) 
      !data idum / 6*-99 / 
      integer*4, PARAMETER :: idum(6)=-99
      integer :: i

      write(ifile) (idum(i),i=1,6) 

      return 
      END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
      SUBROUTINE SKIPFL (NUMFL,IFILE,IEOF) 
   !--------------------------------------------------------------------
      integer ,intent(in)  :: NUMFL,IFILE
      integer ,intent(out) :: IEOF

      integer :: ISKIP      
      real :: DUM(1) !DIMENSION DUM(1) 

      IF (NUMFL.LE.0) RETURN 
      ISKIP = 0 
   10 CALL BUFIN (IFILE,IEOF,DUM(1),1) 
      IF (IEOF.EQ.1) GO TO 10 
      ISKIP = ISKIP+1 
      IF (ISKIP.LT.NUMFL) GO TO 10 

      RETURN 
      END SUBROUTINE
   
   
   
END MODULE 
