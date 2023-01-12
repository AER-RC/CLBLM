!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
Module Module_Utility       

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: FTIME, &
             getlun, &
             EXPINT, &
             EXPINT_real8, &
             upper, lower, &
             find_replace_string, &
             splitFilePath
   

   !yma: This doesn't work since Fortran can not distinguish real and real8 when compiling with '-r8' option.
   !INTERFACE EXPINT      
   !   module procedure EXPINT_real
   !   module procedure EXPINT_real8
   !END INTERFACE
   
   
   CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
   CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

   
CONTAINS !===================== MODULE CONTAIS =========================

      
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE FTIME (HTIME) 
   !--------------------------------------------------------------------
      CHARACTER(8) ,intent(out) :: htime 

      CHARACTER(8) :: gtime
      INTEGER*4 IARRAY(3) 

      CALL ITIME (IARRAY) 
      WRITE (GTIME,900) IARRAY 
      
      !CHANGE THESE TO WORD SIZE AND FORMAT OF ROUTINE GTIME             
      READ (GTIME,901) HTIME 
      RETURN 

!>900 FORMAT (1X,I2,2(':',I2.2))                                        
  900 FORMAT (   I2,2(':',I2.2)) 
!>901 FORMAT (1X,A8)                                                    
  901 FORMAT (A8) 

   END SUBROUTINE
      
   !--------------------------------------------------------------------
   ! Returen a unused file lun.
   !--------------------------------------------------------------------
   integer FUNCTION getlun()
   !--------------------------------------------------------------------
      !--- Save the saveLun to deal with the case when getting LUNs for a group of 
      ! files together and then opening them latterly. Otherwise it'll return a same
      ! LUN number for all files in the group.
      integer ,parameter :: minLun=21, maxLun = 99
      integer      ,SAVE :: saveLun=minLun-1
      integer            :: loopCount
      logical            :: op
      
      loopCount=0
      do
         saveLun = saveLun + 1
         if (saveLun>maxLun) saveLun = minLun
         
         INQUIRE( saveLun, OPENED = op )         
         IF ( .not.op ) then
            getlun = saveLun
            EXIT
         endif

         loopCount = loopCount + 1
         if ( loopCount > (maxLun-minLun+1) ) then
            STOP '--- getlun(): LUN numbers (11~99) used up.'
         endif         
      enddo
      
   END FUNCTION

   !--------------------------------------------------------------------
   !     THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY      
   !     THE FACTOR A. This is a default real precision version of EXPINT.                                                      
   !--------------------------------------------------------------------
   PURE SUBROUTINE EXPINT (X,X1,X2,A) 
   !--------------------------------------------------------------------
      real ,intent(in)  :: X1,X2,A
      real ,intent(out) :: X
      
      IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10 
      X = X1*(X2/X1)**A 
   
      RETURN 
   
   10 X = X1+(X2-X1)*A 
   
      RETURN 
   
   END SUBROUTINE

   !--------------------------------------------------------------------
   !     THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY      
   !     THE FACTOR A, This is a double precision version or EXPINT.
   !--------------------------------------------------------------------
   PURE SUBROUTINE EXPINT_real8 (X,X1,X2,A) 
   !--------------------------------------------------------------------
      integer, parameter :: r8=selected_real_kind(14) !yma: Better put this into a separate kind definition module.

      real(r8) ,intent(in)  :: X1,X2
      real     ,intent(in)  :: A
      real(r8) ,intent(out) :: X
      
      IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10 
      X = X1*(X2/X1)**A 

      RETURN 

   10 X = X1+(X2-X1)*A 

      RETURN 

   END SUBROUTINE

   


   !-------------------------------------------------------------------
   !Both of these were obtained from:
   !
   !   Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
   !   1995 Springer-Verlag, New York.
   !
   !--Paul van Delst     
   !-------------------------------------------------------------------
   PURE FUNCTION upper( Input_String ) RESULT ( Output_String )
   !-------------------------------------------------------------------
     ! -- Argument and result
     CHARACTER( * ), INTENT( IN )     :: Input_String
     CHARACTER( LEN( Input_String ) ) :: Output_String
     ! -- Local variables
     INTEGER :: i, n

     ! -- Copy input string
     Output_String = Input_String
     ! -- Loop over string elements
     DO i = 1, LEN( Output_String )
       ! -- Find location of letter in lower case constant string
       n = INDEX( LOWER_CASE, Output_String( i:i ) )
       ! -- If current substring is a lower case letter, make it upper case
       IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
     END DO
   END FUNCTION

   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   PURE FUNCTION lower( Input_String ) RESULT ( Output_String )
   !-------------------------------------------------------------------
     ! -- Argument and result
     CHARACTER( * ), INTENT( IN )     :: Input_String
     CHARACTER( LEN( Input_String ) ) :: Output_String
     ! -- Local variables
     INTEGER :: i, n

     ! -- Copy input string
     Output_String = Input_String
     ! -- Loop over string elements
     DO i = 1, LEN( Output_String )
       ! -- Find location of letter in upper case constant string
       n = INDEX( UPPER_CASE, Output_String( i:i ) )
       ! -- If current substring is an upper case letter, make it lower case
       IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
     END DO
   END FUNCTION

   
!   !-------------------------------------------------------------------
!   ! Find the location of the first .true. element in logicArray.
!   ! * 'ind' is counted from 1 which is the first location in 'logicArray'
!   ! * If no element found, 'ind' returns 0.
!   ! * Fortran2008 has 'findloc()' intrinsic.
!   !-------------------------------------------------------------------
!   PURE FUNCTION findIndex( logicArray ) RESULT ( ind )
!   !-------------------------------------------------------------------
!      logical, intent(in) :: logicArray(:)
!      integer             :: ind
!
!      integer :: i
!
!      if (.not.any(logicArray==.TRUE.)) then      
!         ind = 0         
!      else
!         do i = 1, size(logicArray)
!            if (logicArray(i)==.TRUE.) then
!               ind = i
!               EXIT
!            endif
!         enddo
!      endif
!
!   END FUNCTION   


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    subroutine splitFilePath(fullpath, directory, basename)
!-----------------------------------------------------------------------
        implicit none

        character(len=*), intent(in) :: fullpath
        character(len=*), intent(out) :: basename, directory
        integer :: lastSlash

        basename = repeat(' ', len(basename))
        directory = repeat(' ', len(directory))

        lastSlash = index(fullpath, '/', BACK=.true.)
        if (lastSlash .gt. 0) directory = fullpath(1:lastSlash)
        basename = fullpath(lastSlash + 1:len(trim(fullpath)))
    end subroutine

   
!-----------------------------------------------------------------------
! Find last "word" in "string" and replace it with "newWord"
!-----------------------------------------------------------------------
   logical FUNCTION find_replace_string( string, word, newWord )
!-----------------------------------------------------------------------
      character(*) ,intent(inout) :: string
      character(*) ,intent(in)    :: word
      character(*) ,intent(in)    :: newWord

      integer :: i,j,k,lw
      character(len(string)) :: tempStr
      
      find_replace_string = .FALSE.
      
      i = index( string, word, BACK=.true.)
      if (i>0) then
         lw = len(trim(word))         
         tempStr = string(1:i-1)//trim(newWord)//trim(string(i+lw:))
         string = trim(tempStr)
         find_replace_string = .TRUE.
      endif
      
   END FUNCTION
        
END MODULE

