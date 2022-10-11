!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

Module Module_LineData
   USE Module_ConstParam ,ONLY: FILLINT, FILLREAL

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_LineData, &
             MXBRDMOL, MXBRDMOL3, NLINEREC, &
             INPUT_HEADER, &
             INPUT_BLOCK, &
             LINE_DATA, &
             LINE_SHRINK
             
   
   
   TYPE CLBLM_LineData
      integer                   :: fileLun     =FILLINT   !Line data file Lun
      !
      character(8)              :: rev_num     =''        !version number of the lnfl.f program
      integer                   :: NLTE_flag   =FILLINT   !flag indicating whether the non-LTE info. included in the file
      integer                   :: REJ_flag    =FILLINT   !flag indicating line rejection was selected
      integer                   :: iso_flag    =FILLINT   !flag indicating compatibility with isotopes
      integer                   :: negepp_flag =FILLINT
      !
      integer*4                 :: linmol      =FILLINT   !Number of molecules included in the line data file.
      integer*4                 :: lincnt      =FILLINT   !Total number of lines
      integer*4                 :: ilinlc      =FILLINT   !Number of coupled lines
      integer*4                 :: ilinnl      =FILLINT   !Number of NLTE lines
      real*4                    :: flinlo      =FILLREAL  !Lowest line
      real*4                    :: flinhi      =FILLREAL  !Highest line
      character(8) ,allocatable :: BMOLID(:)              !Molecule ID
      integer*4    ,allocatable :: molcnt(:)              !Number of lines 
      integer*4    ,allocatable :: mcntlc(:)              !Number of coupled lines
      integer*4    ,allocatable :: mcntnl(:)              !Number of NLTE lines
      real*4       ,allocatable :: sumstr(:)              !Sum of line strengths
      integer*4    ,allocatable :: N_NEGEPP(:)            !
      integer*4    ,allocatable :: N_RESETEPP(:)          !
   CONTAINS 
      !procedure ,pass :: init => CLBLM_LineData_init
   END TYPE


   integer ,PARAMETER :: MXBRDMOL  = 7
   integer ,PARAMETER :: MXBRDMOL3 = MXBRDMOL*3
   integer ,PARAMETER :: NLINEREC  = 250
   
   TYPE :: INPUT_HEADER
      SEQUENCE
      REAL*8    :: VMIN
      real*8    :: VMAX
      INTEGER*4 :: NREC
      integer*4 :: NWDS
   END TYPE INPUT_HEADER
   
   TYPE :: INPUT_BLOCK !total length is 39*250 4-byte words
      SEQUENCE
      REAL*8,    DIMENSION(NLINEREC)            :: VNU
      REAL*4,    DIMENSION(NLINEREC)            :: SP, ALFA, EPP
      INTEGER*4, DIMENSION(NLINEREC)            :: MOL
      REAL*4,    DIMENSION(NLINEREC)            :: HWHM, TMPALF, PSHIFT
      INTEGER*4, DIMENSION(NLINEREC)            :: IFLG
      INTEGER*4, DIMENSION(MXBRDMOL,NLINEREC)   :: BRD_MOL_FLG_IN
      REAL*4,    DIMENSION(MXBRDMOL3,NLINEREC)  :: BRD_MOL_DAT
      REAL*4,    DIMENSION(NLINEREC)            :: SPEED_DEP
   END TYPE INPUT_BLOCK
   
   TYPE :: LINE_DATA
      REAL*8,    DIMENSION(NLINEREC)           :: VNU
      REAL,      DIMENSION(NLINEREC)           :: SP, ALFA, EPP
      INTEGER*4, DIMENSION(NLINEREC)           :: MOL
      REAL,      DIMENSION(NLINEREC)           :: HWHM, TMPALF, PSHIFT
      INTEGER,   DIMENSION(NLINEREC)           :: IFLG
      REAL,      DIMENSION(NLINEREC)           :: SPPSP, RECALF, ZETAI
      INTEGER,   DIMENSION(NLINEREC)           :: IZETA
      INTEGER,   DIMENSION(MXBRDMOL,NLINEREC)  :: BRD_MOL_FLG
      REAL,      DIMENSION(MXBRDMOL,NLINEREC)  :: BRD_MOL_HW, BRD_MOL_TMP, BRD_MOL_SHFT
      REAL,      DIMENSION(NLINEREC)           :: SPEED_DEP
   END TYPE LINE_DATA
   
   TYPE :: LINE_SHRINK
      REAL*8, DIMENSION(1250)  :: VNU
      REAL,    DIMENSION(1250)  :: SP, ALFA, EPP
      INTEGER, DIMENSION(1250)  :: MOL
      REAL,    DIMENSION(1250)  :: SPP, SRAD
   END TYPE LINE_SHRINK

   
CONTAINS !===================Module Contains============================


!-----------------------------------------------------------------------
!     Open the line data file and read in the header information 
!-----------------------------------------------------------------------
   SUBROUTINE CLBLM_LineData_init( this,lineFile )
!-----------------------------------------------------------------------
      USE Module_Utility  ,ONLY: getLun,upper
      
      type(CLBLM_LineData)     ,intent(out) :: this
      character(*)             ,intent(in)  :: lineFile

      
      !--- Local variables
      !
      integer :: lun
      logical :: ex,op
      
      integer ,PARAMETER :: MaxNumMol = 64
      
      character(8) :: HLINID(10),HID1(2)
      character(8) :: BMOLID(MaxNumMol)
      integer*4    :: MOLCNT(MaxNumMol)
      integer*4    :: MCNTLC(MaxNumMol)
      integer*4    :: MCNTNL(MaxNumMol)
      real*4       :: SUMSTR(MaxNumMol)
      real*4       :: FLINLO,FLINHI 
      integer*4    :: LINMOL,LINCNT,ILINLC,ILINNL,IREC,IRECTL
      
      CHARACTER(1) :: CNEGEPP(8) 
      integer*4    :: NEGEPP_FLAG
      integer*4    :: N_NEGEPP(MaxNumMol)
      integer*4    :: N_RESETEPP(MaxNumMol)
      real*4       :: XSPACE(4096)
      
      character    :: rev_num*8
      integer      :: NLTE_flag, REJ_flag, iso_flag
      
      
      !--- Open the line data file
      !
      inquire( FILE=trim(lineFile), EXIST=ex, opened=op, number=lun )
      if (.NOT.ex) then
         STOP '--- CLBLM_LineData_init(): Line data file not exist, program stopped.'
      endif
      
      if (.not.op) then
         lun = getlun()
         OPEN (lun,FILE=trim(lineFile),STATUS='OLD',FORM='UNFORMATTED') 
      else
         rewind(lun)
      endif
      
      
      !--- Read in the file header
      READ(lun) HLINID,BMOLID,MOLCNT,MCNTLC, &
                MCNTNL,SUMSTR,LINMOL,FLINLO,FLINHI, &
                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1 

      !---
      ! Test for negative values of ENERGY identified in lnfl             
      ! and read in second header for line information, if needed         
      negepp_flag = 0
      READ (HLINID(7),950) CNEGEPP 
      IF (CNEGEPP(8).eq.'^') THEN 
         negepp_flag = 1 
         READ (lun) n_negepp,n_resetepp,xspace 
      endif 


      !--- Get line data file version 
!180422: Current lnfl.f has a bug, the version number in TAPE3 is wrong. Only the first four digits were saved.      
      rev_num = HLINID(8)

      
      !--- Check NLTE and REJ flags
      ! * REJ_flag is not use in CLBLM. NLTE_flag does used.
!should we add NOCPL and EXBRD flags?
      select case ( upper(trim(adjustl(HLINID(9)))) )
      case ('NLTE')
         NLTE_flag = 1
         REJ_flag =0         
      case ('REJ')
         NLTE_flag = 0
         REJ_flag =1         
      case ('NLTE REJ')
         NLTE_flag = 1
         REJ_flag =1               
      case default
         NLTE_flag = 0
         REJ_flag =0      
         !STOP '--- CLBLM_LineData_init(): Invalid flags in line data file. Program stopped.'
      end select 

      
      !--- Check for isotope inclusion in line data file
      !    170213: iso_flag is not used in CLBLM.
      if ( upper(trim(adjustl(HLINID(10))))=='I' ) then
         iso_flag = 1
      else
         iso_flag = 0
      endif


      !--- Reading header finished. Rewind file for later use.
      rewind(lun)


      
      !--- Load the CLBLM_LineData object
      !
      if (allocated( this%BMOLID     )) deallocate( this%BMOLID ) !Molecule ID
      if (allocated( this%molcnt     )) deallocate( this%molcnt ) !Number of lines 
      if (allocated( this%mcntlc     )) deallocate( this%mcntlc ) !Number of coupled lines
      if (allocated( this%mcntnl     )) deallocate( this%mcntnl ) !Number of NLTE lines
      if (allocated( this%sumstr     )) deallocate( this%sumstr ) !Sum of line strengths
      if (allocated( this%N_NEGEPP   )) deallocate( this%N_NEGEPP )
      if (allocated( this%N_RESETEPP )) deallocate( this%N_RESETEPP )
      
      allocate( this%BMOLID(MaxNumMol) ); this%BMOLID = '' !Molecule ID
      allocate( this%molcnt(MaxNumMol) )                   !Number of lines 
      allocate( this%mcntlc(MaxNumMol) )                   !Number of coupled lines
      allocate( this%mcntnl(MaxNumMol) )                   !Number of NLTE lines
      allocate( this%sumstr(MaxNumMol) )                   !Sum of line strengths
      if (negepp_flag==1) then
         allocate( this%N_NEGEPP(MaxNumMol) )
         allocate( this%N_RESETEPP(MaxNumMol) )
      endif
      
      this%fileLun     = lun         !Line data file Lun
      this%rev_num     = rev_num     !version number of the lnfl.f program
      this%NLTE_flag   = NLTE_flag   !flag indicating whether the non-LTE info. included in the file
      this%REJ_flag    = REJ_flag    !flag for line rejection
      this%iso_flag    = iso_flag    !flag indicating compatibility with isotopes
      this%negepp_flag = negepp_flag
      this%linmol      = linmol      !Number of molecules included in the line data file.
      this%lincnt      = lincnt      !Total number of lines
      this%ilinlc      = ilinlc      !
      this%ilinnl      = ilinnl      !
      this%flinlo      = flinlo      !Lowest line
      this%flinhi      = flinhi      !Highest line
      this%BMOLID(:)   = BMOLID(:)   !Molecule ID
      this%molcnt(:)   = molcnt(:)   !Number of lines 
      this%mcntlc(:)   = mcntlc(:)   !Number of coupled lines
      this%mcntnl(:)   = mcntnl(:)   !Number of NLTE lines
      this%sumstr(:)   = sumstr(:)   !Sum of line strengths
      if (negepp_flag==1) then
         this%N_NEGEPP(:)   = N_NEGEPP(:)   
         this%N_RESETEPP(:) = N_RESETEPP(:)
      endif
      
      
      !--- Print out line file header
      call PRLNHD( HLINID,LINMOL,LINCNT,FLINLO,FLINHI,BMOLID,&
                   MOLCNT,MCNTLC,MCNTNL,SUMSTR,N_NEGEPP,N_RESETEPP,&
                   negepp_flag,HID1 )
      
      
      RETURN
      
  950 FORMAT (8a1) 
      
      END SUBROUTINE
      
!-----------------------------------------------------------------------
!     PRLNHD PRINTS OUT LINE FILE HEADER                                
!-----------------------------------------------------------------------
      SUBROUTINE PRLNHD( HLINID,LINMOL,LINCNT,FLINLO,FLINHI,BMOLID,&
                         MOLCNT,MCNTLC,MCNTNL,SUMSTR,N_NEGEPP,N_RESETEPP,&
                         negepp_flag,HID1 )
!-----------------------------------------------------------------------
      USE Module_Config ,ONLY: IPR
      IMPLICIT NONE

      character(*)   ,intent(in) :: HLINID(10)
      integer*4      ,intent(in) :: LINMOL
      integer*4      ,intent(in) :: LINCNT
      real*4         ,intent(in) :: FLINLO,FLINHI
      character(*)   ,intent(in) :: BMOLID(:)
      integer*4      ,intent(in) :: MOLCNT(:),MCNTLC(:),MCNTNL(:)
      real*4         ,intent(in) :: SUMSTR(:)
      integer*4      ,intent(in) :: N_NEGEPP(:),N_RESETEPP(:)
      integer*4      ,intent(in) :: negepp_flag
      character(*)   ,intent(in) :: HID1(2)

      
      integer     :: M,I
      CHARACTER   :: CHID10*8,CHARID*5,CHARDT*2,CHARI*1,CHTST*1 
  
      DATA CHARI / 'I'/ 
      
 
      WRITE (IPR,900) 
      WRITE (IPR,905) HLINID,HID1 

      ! Output header information regarding lines; if negative values of  
      ! ENERGY were identified in lnfl, output extra header information   

      if (negepp_flag==1) THEN !yma if (CNEGEPP(8).eq.'^') THEN 
         WRITE (IPR,960) 
         WRITE (IPR,965) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),      &
                          N_NEGEPP(I),N_RESETEPP(I), SUMSTR(I),I=1,LINMOL)
      else 
         WRITE (IPR,910) 
         WRITE (IPR,915) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),      &
                          SUMSTR(I),I=1,LINMOL)                                          
      endif 

      WRITE (IPR,920) FLINLO,FLINHI,LINCNT 

      !     When calculating derivative, check make sure the                  
      !     appropriate molecule is included in the linefile.                 
      !     If not, then stop and issue message.                              

      !     CHECK HEADER FOR FLAG INDICATING COMPATIBILITY WITH ISOTOPES      

   30 WRITE (CHID10,925) HLINID(10) 
      READ (CHID10,930) CHARID,CHARDT,CHTST 
      IF (CHTST.NE.CHARI) THEN 
         WRITE (IPR,935) CHARID,CHARDT,CHTST 
         STOP ' PRLNHD - NO ISOTOPE INFO ON LINFIL ' 
      ENDIF 

      RETURN 

  900 FORMAT ('0'/'0',20X,'   LINE FILE INFORMATION ') 
  905 FORMAT ('0',10A8,2X,2(1X,A8,1X)) 
!yma  910 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'SUM LBLRTM ',/,7X,      &
!yma     &        'MOL',5X,'LINES',4X,'LINES',4X,'LINES',4X,'STRENGTHS',/)  
  910 FORMAT ('0',/,24X,'COUPLED',4X,'NLTE',3X,'SUM LBLRTM ',/,7X,      &
     &        'MOL',6X,'LINES',4X,'LINES',4X,'LINES',4X,'STRENGTHS',/)  
!yma  915 FORMAT (' ',4X,A6,' = ',I6,3X,I6,3X,I6,2X,1PE12.4,0P) 
  915 FORMAT (' ',4X,A6,' = ',I7,3X,I6,3X,I6,2X,1PE12.4,0P) 
  920 FORMAT (/,'0 LOWEST LINE = ',F10.3,5X,'  HIGHEST LINE = ',F10.3,  &
     &        5X,' TOTAL NUMBER OF LINES =',I8)                         
  925 FORMAT (A8) 
  930 FORMAT (A5,A2,A1) 
  935 FORMAT (3(/),10X,'LINEFILE PROGRAM: ',A5,3X,'VERSION: ',A2,A1,    &
     &        3(/),3X,52('*'),/,3X,'* THE LINEFILE (TAPE3) IS NOT ',    &
     &        'COMPATIBLE WITH THIS *',/,3X,'* VERSION OF LBLRTM .',    &
     &        '  ISOTOPIC INFORMATION (FROM  *',/,3X,'* HITRAN) ',      &
     &        'MUST BE PRESERVED ON TAPE3.  USE A TAPE3 *',/,3X,        &
     &        '* CREATED WITH THE 91I OR LATER VERSION OF LNFL.   *',   &
     &        /,3X,52('*'))                                             
  940 FORMAT (' Molecule to be retrieved: ',A6,' not in linefile.',/,   &
     &        ' Molecules in linefile: ')                               
  945 FORMAT (24X,A6) 
  950 FORMAT (8a1) 
!yma  960 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'NEGATIVE',3X,           &
!yma     &        'RESET',4X,'SUM LBLRTM',/,7X,'MOL',5X,'LINES',4X,         &
!yma     &        'LINES',4X,'LINES',6X,'EPP',6X,'EPP',                     &
!yma     &        6X,'STRENGTHS',/)                                         
  960 FORMAT ('0',/,24X,'COUPLED',4X,'NLTE',3X,'NEGATIVE',3X,           &
     &        'RESET',4X,'SUM LBLRTM',/,7X,'MOL',6X,'LINES',4X,         &
     &        'LINES',4X,'LINES',6X,'EPP',6X,'EPP',                     &
     &        6X,'STRENGTHS',/)                                         
!yma  965 FORMAT (' ',4X,A6,' = ',I6,                                       &
!yma     &        3X,I6,3X,I6,3X,I6,3X,i6,3X,1PE12.4)                       
  965 FORMAT (' ',4X,A6,' = ',I7,                                       &
     &        3X,I6,3X,I6,3X,I6,3X,i6,3X,1PE12.4)                       

      END SUBROUTINE

END MODULE

