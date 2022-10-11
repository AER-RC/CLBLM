!     path:      %P%                                                    
!     revision:  $Revision: 10548 $                                     
!     created:   $Date: 2011-10-07 13:05:24 -0400 (Fri, 07 Oct 2011) $  
!     presently: %H%  %T%                                               
!_______________________________________________________________________
!yma                                                                        
!yma      SUBROUTINE LBLDAT(HDATE) 
!yma!                                                                       
!yma      character*8 hdate 
!yma      integer*4 iday,imon,iyear 
!yma!                                                                       
!yma      CHARACTER GDATE*8 
!yma                                                                        
!yma      COMMON /CVRUTL/ HNAMUTL,HVRUTL 
!yma!                                                                       
!yma      CHARACTER*18 HNAMUTL,HVRUTL 
!yma!                                                                       
!yma!     ASSIGN NAME and CVS VERSION NUMBER TO MODULE                      
!yma!                                                                       
!yma      HNAMUTL= 'util_linux_intel.f:' 
!yma      HVRUTL = '$Revision: 10548 $' 
!yma                                                                        
!yma!                                                                       
!yma      CALL IDATE(iday,imon,iyear) 
!yma!                                                                       
!yma      iyear=iyear-(100*(iyear/100)) 
!yma!                                                                       
!yma      WRITE (GDATE,900) iyear,imon,iday 
!yma!                                                                       
!yma      READ (GDATE,901) HDATE 
!yma!                                                                       
!yma!     CHANGE THESE TO WORD SIZE AND FORMAT OF ROUTINE DATE              
!yma!                                                                       
!yma      RETURN 
!yma!                                                                       
!yma  900 FORMAT (   I2.2,2('/',I2.2)) 
!yma  901 FORMAT (A8) 
!yma!                                                                       
!yma      END                                           
!yma      SUBROUTINE FTIME (HTIME) 
!yma!                                                                       
!yma!                                                                       
!yma      CHARACTER*8 gtime,htime 
!yma!                                                                       
!yma      INTEGER*4 IARRAY(3) 
!yma!                                                                       
!yma      CALL ITIME (IARRAY) 
!yma!                                                                       
!yma      WRITE (GTIME,900) IARRAY 
!yma!                                                                       
!yma      READ (GTIME,901) HTIME 
!yma!                                                                       
!yma!     CHANGE THESE TO WORD SIZE AND FORMAT OF ROUTINE GTIME             
!yma!                                                                       
!yma      RETURN 
!yma!                                                                       
!yma!>900 FORMAT (1X,I2,2(':',I2.2))                                        
!yma  900 FORMAT (   I2,2(':',I2.2)) 
!yma!>901 FORMAT (1X,A8)                                                    
!yma  901 FORMAT (A8) 
!yma!                                                                       
!yma      END                                           
      SUBROUTINE CPUTIM (TIME) 
!                                                                       
      COMMON /TIMIN/ A1 
!                                                                       
      REAL*4 ETIME,TARRAY(2) 
!                                                                       
!     THIS SUBROUTINE OBTAINS CPU TIME                                  
!                                                                       
      IF (A1.LE.0.) THEN 
         A1 = ETIME(TARRAY) 
         TIME = a1 
      ELSE 
         TIME = ETIME(TARRAY) 
      ENDIF 
!                                                                       
      RETURN 
!                                                                       
      END                                           
      BLOCK DATA BTIM 
!                                                                       
      COMMON /TIMIN/ A1 
!                                                                       
      DATA A1 / 0. / 
!                                                                       
      END                                           
                                                                        
      SUBROUTINE BUFIN (IFILE,IEOF,IARRAY,IWORDS) 
!                                                                       
!     THIS SUBROUTINE BUFFERS IN (READS) IWORDS INTO  IARRAY STARTING   
!     AT LOCATION IARRAY                                                
!                                                                       
!     IFILE IS THE FILE DESIGNATION                                     
!                                                                       
      DATA i_4 / 4 / 
!                                                                       
      DIMENSION IARRAY(IWORDS) 
!                                                                       
!                                                                       
      IEOF = 1 
!                                                                       
!#    BUFFER IN (IFILE,1) (IARRAY(ILO),IARRAY(IHI))                     
!#    IF (UNIT(IFILE).EQ.0.) GO TO 10                                   
!                                                                       
      READ (IFILE,END=10) IARRAY 
      ITEST = MIN(IWORDS,i_4) 
      IF (IARRAY(ITEST).EQ.-99) IEOF = -99 
!                                                                       
      RETURN 
!                                                                       
   10 IEOF = 0 
!                                                                       
      RETURN 
!                                                                       
      END                                           
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                         
                                                                        
!         note the name change                                          
                                                                        
      SUBROUTINE BUFIN_sgl (IFILE,IEOF,IARRAY,IWORDS) 
!                                                                       
!     THIS SUBROUTINE BUFFERS IN (READS) IWORDS INTO  IARRAY STARTING   
!     AT LOCATION IARRAY                                                
!                                                                       
!     IFILE IS THE FILE DESIGNATION                                     
!                                                                       
      implicit integer*4 (i-n) 
      implicit real*4    (a-h,o-z) 
                                                                        
      DATA i_4 / 4 / 
                                                                        
      DIMENSION IARRAY(IWORDS) 
!                                                                       
      IEOF = 1 
!                                                                       
!#    BUFFER IN (IFILE,1) (IARRAY(ILO),IARRAY(IHI))                     
!#    IF (UNIT(IFILE).EQ.0.) GO TO 10                                   
!                                                                       
      READ (IFILE,END=10) IARRAY 
      ITEST = MIN(IWORDS,i_4) 
      IF (IARRAY(ITEST).EQ.-99) IEOF = -99 
!                                                                       
      RETURN 
!                                                                       
   10 IEOF = 0 
!                                                                       
      RETURN 
!                                                                       
      END                                           
!_______________________________________________________________________
                                                                        
      SUBROUTINE BUFOUT (IFILE,IARRAY,IWORDS) 
!                                                                       
!     THIS SUBROUTINE BUFFERS OUT (WRITES) IWORDS FROM IARRAY STARTING  
!     AT LOCATION IARRAY                                                
!                                                                       
!     IFILE IS THE FILE DESIGNATION                                     
!                                                                       
      DIMENSION IARRAY(IWORDS) 
!                                                                       
!#    BUFFER OUT (IFILE,1) (IARRAY(ILO),IARRAY(IHI))                    
!#    IF (UNIT(IFILE).EQ.0.) STOP ' ERROR IN BUFOUT '                   
!                                                                       
      WRITE (IFILE) IARRAY 
!                                                                       
      RETURN 
!                                                                       
      END                                           
!_______________________________________________________________________
                                                                        
      SUBROUTINE BUFOUT_sgl (IFILE,IARRAY,IWORDS) 
!                                                                       
!     THIS SUBROUTINE BUFFERS OUT (WRITES) IWORDS FROM IARRAY STARTING  
!     AT LOCATION IARRAY                                                
!                                                                       
!     IFILE IS THE FILE DESIGNATION                                     
!                                                                       
!                                                                       
      implicit integer*4 (i-n) 
      implicit real*4    (a-h,o-z) 
!                                                                       
      DIMENSION IARRAY(IWORDS) 
!                                                                       
!#    BUFFER OUT (IFILE,1) (IARRAY(ILO),IARRAY(IHI))                    
!#    IF (UNIT(IFILE).EQ.0.) STOP ' ERROR IN BUFOUT '                   
!                                                                       
      WRITE (IFILE) IARRAY 
!                                                                       
      RETURN 
!                                                                       
      END                                           
