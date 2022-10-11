!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE Module_Scn_Config
!-----------------------------------------------------------------------
   IMPLICIT NONE   
   PRIVATE             
   PUBLIC :: IPR, openLogFile !openLogFile open logfile and set unit number to IPR
   
   !--- Global variables
   ! "IPR" is used so often, so make it a global variable.
   integer  ,SAVE :: IPR
   !character(256) :: logFile   = 'sceneToolLog.txt'
   !character(256) :: sceneFile = 'clblm_scene.nc'           
   
CONTAINS !=================== MODULE CONTAINS ==========================

 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
   SUBROUTINE openLogFile( logFile )
 !----------------------------------------------------------------------
      USE Module_Utility ,ONLY: getLun
      IMPLICIT NONE

      character(*) ,intent(in) :: logFile
      
      logical :: op
   
      inquire( FILE=logFile, opened=op )
      if (.not.op) then 
         IPR = getLun() !initialize the global variable IPR
         OPEN (IPR,FILE=trim(logFile),STATUS='UNKNOWN') 
      endif
      
   END SUBROUTINE
      
END MODULE
