!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
!-----------------------------------------------------------------------
! This module read the TAPE5 into a data structure
!-----------------------------------------------------------------------
MODULE  Module_Scn_TAPE5
   USE Module_ConstParam, ONLY: r8=>kind_r8

   !---
   ! * Type definitions have to visible to module user to access the 
   !   record contents.
   IMPLICIT NONE  
   
   
   TYPE :: Rec_1_1_Type
      character(8) :: XID(10)
   END TYPE
   

   TYPE :: Rec_1_2_Type
      integer :: IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT_tape5,ISCAN,&
                 IFILTR,IPLOT,ITEST, IATM
      integer :: IMRG !character(1) :: CMRG(2)
      integer :: ILAS,IOD,IXSECT,IRAD,MPTS,NPTS!,ISPD
      integer :: ISOTPL, IBRD !LBLRTM v12.7
   END TYPE
   
   
   TYPE :: Rec_1_2a_Type
      real :: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 
   END TYPE
   
   
   TYPE :: Rec_1_2_1_Type
      integer :: INFLAG,IOTFLG,JULDAT
      integer :: ISOLVAR
      real    :: SCON,SOLCYCFRAC,SOLVAR
   END TYPE
   
   
   TYPE :: Rec_1_3_Type
      real(r8) :: V1,V2
      real     :: SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,DPTFAC
      integer  :: ILNFLG
      real     :: DVOUT
      integer  :: nmol_scal, n_xs_scal
   END TYPE
   
   
   TYPE :: Rec_1_3a_and_1_3b_Type
      character(1), allocatable :: hmol_scal(:)
      real,         allocatable :: xmol_scal(:)
      character(1), allocatable :: h_xs_scal(:)
      real,         allocatable :: x_xs_scal(:)
   END TYPE
      
   
   TYPE :: Rec_1_4_Type
      real         :: TMPBND, BNDEMI(3), BNDRFL(3)
      character(1) :: surf_refl
   END TYPE
   
   
   TYPE :: Rec_1_5_Type
      integer :: NSPCRT
   END TYPE
   
   
   TYPE :: Rec_1_6a_Type
      character(55) :: PATH1
      integer       :: laytot
   END TYPE

 
   TYPE :: Rec_2_1_Type
      integer       :: IFORM,NLAYRS,NMOL
      real          :: SECNT0
      character(20) :: HEAD20
      real          :: ZH1
      character(4)  :: HEAD4
      real          :: ZH2
      character(5)  :: HEAD5
      character(8)  :: hol_angle
      character(7)  :: HEAD7
   END TYPE
   
   
   TYPE :: rec_2_1_1_to_2_1_3_Type
      real,         allocatable :: PAVE(:),TAVE(:),SECNTK(:)
      character(3), allocatable :: CINP(:)
      integer,      allocatable :: IPTHRK(:)
      real,         allocatable :: ALTZ(:),PZ(:),TZ(:)
      real,         allocatable :: WKL(:,:), WBRODL(:)
   END TYPE
   
   
   TYPE :: Rec_2_2_Type
      integer :: IXMOLS=0 !default initialize to 0, since scene reader relies on the value of IXMOLS to determine the presence of the x-section data.
      integer :: IXSBIN
   END TYPE


   TYPE :: Rec_2_2_1_Type
      character(10), allocatable :: XSNAME(:)
   END TYPE
   
   
   TYPE :: Rec_2_2_2_Type
      integer       :: IFRMX,NLAYXS,IXMOL
      real          :: SECNTX
      character(60) :: HEDXS
   END TYPE
   
   
   TYPE :: rec_2_2_3_to_2_2_5_Type
      real,         allocatable :: PAVX(:),TAVX(:),SECKXS(:)
      character(3), allocatable :: CINPX(:)
      integer,      allocatable :: IPTHKX(:)
      real,         allocatable :: ALTZX(:),PZX(:),TZX(:)
      real,         allocatable :: XAMNT(:,:),WBRODX(:)
   END TYPE
   
   
   TYPE :: Rec_3_1_Type
      integer       :: MODEL,ITYPE,IBMAX_B,n_zero,NOPRNT,NMOL,IPUNCH,&
                       IFXTYPE,MUNITS
      real          :: RE,HSPACE,XVBAR,dumrd
      character(10) :: sref_lat
   END TYPE
   
   
   TYPE :: Rec_3_2h_Type
      real :: H1F,RANGEF
   END TYPE
   
   
   TYPE :: Rec_3_2_Type
      real    :: H1F,H2F,ANGLEF,RANGEF,BETAF
      integer :: LENF
      real    :: HOBS
   END TYPE


   TYPE :: Rec_3_3a_Type
      real :: AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2
   END TYPE

   
   TYPE :: Rec_3_3b_Type
      real, allocatable :: layerBoundaries(:)
   END TYPE
   
   
   TYPE :: Rec_3_4_Type
      integer      :: IMMAX_B
      character(8) :: HMOD(3)
   END TYPE
   
   
   TYPE :: Rec_3_5_and_3_6_Type
      real,         allocatable :: ZMDL(:),PM(:),TM(:)
      character(1), allocatable :: JCHARP(:),JCHART(:),&
                                   JLONG(:),JCHAR(:,:)
      real,         allocatable :: WMOL(:,:)   
   END TYPE
   
   
   TYPE :: Rec_3_7_Type
      integer :: IXMOLS,IPRFL,IXSBIN
   END TYPE
   
   
   TYPE :: Rec_3_7_1_Type
      character(10), allocatable :: XSNAME(:)
   END TYPE
   
   
   TYPE :: Rec_3_8_Type
      integer       :: LAYX,IZORP
      character(50) :: XTITLE
   END TYPE
   
   
   TYPE :: Rec_3_8_1_Type
      real,         allocatable :: ZX(:)
      real,         allocatable :: PX(:)
      character(1), allocatable :: JCHAR(:,:)
      real,         allocatable :: DTMP(:,:)
   END TYPE

   
   TYPE :: Rec_isotpl_Type !v12.7 
      integer              :: NISOTPL
      integer ,allocatable :: ISOTPL_HCODE(:) ![nISO]
      integer              :: IFRMI, NLAYIS
      real    ,allocatable :: ISOTPL_AMNT_tape5(:,:) ![nISO,nLay] !This reader simply relay the TAPE5 content to the program.
   END TYPE

   
   TYPE :: Rec_6_Type
      real     :: HWHM
      real(r8) :: V1,V2
      integer  :: JEMIT,JFN,JVAR
      real     :: SAMPL
      integer  :: IUNIT,IFILST,NIFILS,JUNIT,NPTS
   END TYPE
   
   
   TYPE :: Rec_6_1_Type
      real :: DIRCOS
   END TYPE
   
   
   TYPE :: Rec_7_1_Type
      real(r8)      :: V1F
      real          :: DVF
      integer       :: NPTF,JEMIT,NNFILE
      character(35) :: HEDDR
   END TYPE
   
   
   TYPE :: Rec_7_2_and_7_3_Type
      character(80)     :: CVAR
      real, allocatable :: XF(:)
   END TYPE

   
   TYPE :: Rec_8_1_Type
      real     :: HWHM
      real(r8) :: V1,V2
      integer  :: JEMIT,JFN,JVAR
      real     :: SAMPL
      integer  :: IUNIT, IFILST,NIFILS,JUNIT,NPTS
   END TYPE
   
   
   TYPE :: Rec_9_1_Type
      real     :: DVO
      real(r8) :: V1,V2
      integer  :: JEMIT,I4PT,IUNIT,IFILST,NIFILS,JUNIT,NPTS
   END TYPE
   
   
   TYPE :: Rec_10_1_Type
      real     :: HWHM
      real(r8) :: V1,V2
      integer  :: JEMIT,JFNIN,MRATIN
      real     :: DVOUT
      integer  :: IUNIT,IFILST,NIFILS,JUNIT,IVX,NOFIX
      real     :: param
   END TYPE
   
   
   TYPE :: Rec_10_2_Type
      real :: PARM1,PARM2,PARM3
   END TYPE
   
   
   TYPE :: rec_11_1_Type
      real(r8)      :: V1F
      real          :: DVF
      integer       :: NPTF,JEMIT,IUNIT,IFILST,NIFILS,junit
      character(35) :: HEDDR
   END TYPE
   

   TYPE :: Rec_11_2_and_11_3_Type
      character(80)     :: CVAR
      real, allocatable :: XF(:)
   END TYPE

   
   TYPE :: CLBLM_TAPE5
      type(Rec_1_1_Type)              :: Rec_1_1
      type(Rec_1_2_Type)              :: Rec_1_2
      type(Rec_1_2a_Type)             :: Rec_1_2a
      type(Rec_1_2_1_Type)            :: Rec_1_2_1
      type(Rec_1_3_Type)              :: Rec_1_3
      type(Rec_1_3a_and_1_3b_Type)    :: Rec_1_3a_and_1_3b
      type(Rec_1_4_Type)              :: Rec_1_4
      type(Rec_1_5_Type)              :: Rec_1_5
      type(Rec_1_6a_Type)             :: Rec_1_6a
      type(Rec_2_1_Type)              :: Rec_2_1
      type(rec_2_1_1_to_2_1_3_Type)   :: rec_2_1_1_to_2_1_3
      type(Rec_2_2_Type)              :: Rec_2_2
      type(Rec_2_2_1_Type)            :: Rec_2_2_1
      type(Rec_2_2_2_Type)            :: Rec_2_2_2
      type(rec_2_2_3_to_2_2_5_Type)   :: rec_2_2_3_to_2_2_5
      type(Rec_3_1_Type)              :: Rec_3_1
      type(Rec_3_2h_Type)             :: Rec_3_2h
      type(Rec_3_2_Type)              :: Rec_3_2
      type(Rec_3_3a_Type)             :: Rec_3_3a
      type(Rec_3_3b_Type)             :: Rec_3_3b
      type(Rec_3_4_Type)              :: Rec_3_4
      type(Rec_3_5_and_3_6_Type)      :: Rec_3_5_and_3_6
      type(Rec_3_7_Type)              :: Rec_3_7
      type(Rec_3_7_1_Type)            :: Rec_3_7_1
      type(Rec_3_8_Type)              :: Rec_3_8
      type(Rec_3_8_1_Type)            :: Rec_3_8_1
      type(Rec_isotpl_Type)           :: Rec_isotpl !v12.7
      type(Rec_6_Type)                :: Rec_6
      type(Rec_6_1_Type)              :: Rec_6_1
      type(Rec_7_1_Type)              :: Rec_7_1
      type(Rec_7_2_and_7_3_Type)      :: Rec_7_2_and_7_3
      type(Rec_8_1_Type)              :: Rec_8_1
      type(Rec_9_1_Type)              :: Rec_9_1
      type(Rec_10_1_Type)             :: Rec_10_1
      type(Rec_10_2_Type)             :: Rec_10_2
      type(rec_11_1_Type)             :: Rec_11_1
      type(Rec_11_2_and_11_3_Type)    :: Rec_11_2_and_11_3   
   END TYPE
   
   !-----------------------------------------
   !type(CLBLM_TAPE5), SAVE :: tape5
   !-----------------------------------------


CONTAINS !====================== MODULE CONTAINS =======================



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE readTAPE5( tape5, tape5FileName )
!-----------------------------------------------------------------------   
      character(*)            ,intent(in)  :: tape5FileName
      type(CLBLM_TAPE5)       ,intent(out) :: tape5
      

      integer   :: IRD,iostat
      character :: fmt*256, iomsg*256
      logical   :: endOfTape5 !end of TAPE5 flag
      
      integer :: IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT_tape5,&
                 ISCAN,IFILTR,IATM,ILAS,IXSECT,&
                 ISOTPL, &
                 nmol_scal, n_xs_scal,&
                 MODEL,ITYPE,IBMAX_B,&
                 IPRFL,&
                 JFNIN
      integer :: IMRG !need for converting CMRG->IMRG


      
      IRD = 55
      endOfTape5 = .FALSE.
      OPEN (IRD,FILE=trim(tape5FileName),STATUS='UNKNOWN',&
            iostat=iostat, iomsg=iomsg )
      
      ! only one cycle is allowed now
      !DO WHILE
         
         !---
         !Record 1.1
         !>>> CXID
         call get_rec_1_1()
         if (endOfTape5) RETURN


         !--- 
         !Record 1.2
         !>>> IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT, &
         !    ITEST,IATM, CMRG,ILAS, IOD,IXSECT, IRAD,MPTS,NPTS,ISOTPL,IBRD !,ISPD                  
         !Or? IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT, &    
         !    ITEST,IATM, CMRG,ILAS, IOD,IXSECT, MPTS,NPTS,ISOTPL,IBRD!,ISPD
         call get_rec_1_2()

         IHIRAC       = tape5%rec_1_2%IHIRAC 
         ILBLF4       = tape5%rec_1_2%ILBLF4 
         ICNTNM       = tape5%rec_1_2%ICNTNM 
         IAERSL       = tape5%rec_1_2%IAERSL 
         IEMIT_tape5  = tape5%rec_1_2%IEMIT_tape5  
         ISCAN        = tape5%rec_1_2%ISCAN  
         IFILTR       = tape5%rec_1_2%IFILTR 
         IATM         = tape5%rec_1_2%IATM
         IMRG         = tape5%rec_1_2%IMRG
         ILAS         = tape5%rec_1_2%ILAS   
         IXSECT       = tape5%rec_1_2%IXSECT
         ISOTPL       = tape5%rec_1_2%ISOTPL !v12.7

         
         !--- 
         !Eecord 1.2a  
         !>>> XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL 
         if (ICNTNM ==6) call get_rec_1_2a() 


         !---
         !Record 1.2.1
         !>>> INFLAG,IOTFLG,JULDAT,ISOLVAR,SCON,SOLCYCFRAC,SOLVAR
         if (IEMIT_tape5 ==2 ) call get_rec_1_2_1()


         !---
         ! Record 1.3
         ! Record 1.3a
         ! Record 1.3b
         if ( (IHIRAC+IAERSL+IEMIT_tape5+IATM+ILAS) >0 ) then

            !Record 1.3
            !>>> V1,V2,SAMPLE,DVSET, &
            !    ALFAL0,AVMASS,DPTMIN,DPTFAC,ILNFLG,DVOUT,&
            !    nmol_scal,n_xs_scal
            call get_rec_1_3()

            nmol_scal = tape5%rec_1_3%nmol_scal
            n_xs_scal = tape5%rec_1_3%n_xs_scal
            

            !Record 1.3a and 1.3b
            !>>> hmol_scal
            !    xmol_scal
            !>>> h_xs_scal
            !    x_xs_scal
            if ( nmol_scal >0 .or. n_xs_scal >0 ) &
                                        call get_rec_1_3a_and_1_3b()
         
         endif



         !---
         ! Record 1.4
         !>>> TMPBND, (BNDEMI(i),i=1,3), (BNDRFL(i),i=1,3), surf_refl
         if ( IEMIT_tape5>0 ) call get_rec_1_4()
         

         
         !---
         !Record 1.5; required for analytic derivative
         !>>> NSPCRT
         if ( IEMIT_tape5==3 .and. any(IMRG==[ 40,41,42,43 ])) then
            call get_rec_1_5()
         endif
         


         !---
         !Record 1.6a         
         !>>> PATH1, laytot
         if (any(IMRG==[ 35,36,40,41,42,43,45,46 ])) then !added 42,43
            call get_rec_1_6a()
         endif
         

         !---
         !%---Record 2.1
         !%---Record 2.1.1
         !%---Record 2.1.2
         !%---Record 2.2
         !%---Record 2.2.1
         !%---Record 2.2.2
         !%---Record 2.2.3
         !%---Record 2.2.4            
         ! user's own path inputs.         
         if (IATM==0) then
 
            !Record 2.1
            !>>> IFORM,NLAYRS,NMOL,SECNT0,HEAD20,ZH1, &
            !    HEAD4,ZH2,HEAD5,hol_angle,HEAD7
            call get_rec_2_1()
 
            !Record 2.1.1 and 2.1.2           
            !>>> PAVE(L),TAVE(L),SECNTK(L),CINP(L),IPTHRK,&
            !    ALTZ(L-1),PZ(L-1),TZ(L-1),&
            !    ALTZ(L),PZ(L),TZ(L)                      
            !>>> WKL(:,:), WBRODL(:)
            call get_rec_2_1_1_to_2_1_3()
             
            if (IXSECT==1) then

               !Record 2.2
               !>>> IXMOLS,IXSBIN
               call get_rec_2_2()

               !Record 2.2.1
               !>>> XSNAME(IXMOLS)
               call get_rec_2_2_1()

               !Record 2.2.2
               !>>> IFRMX,NLAYXS,IXMOL,SECNTX,HEDXS
               call get_rec_2_2_2()

               !Record 2.2.3 and 2.2.4
               !>>> PAVX(L),TAVX(L),SECKXS(L),CINPX(L),IPTHKX,&
               !    ALTZX(L-1),PZX(L-1),TZX(L-1),&
               !    ALTZX(L),PZX(L),TZX(L)
               !>>> XAMNT(:,:),WBRODX(:)
               call get_rec_2_2_3_to_2_2_5()

            endif !if (IXSECT==1) then
         endif !if (IATM==0) then
         
         
         !---
         !Record 3.1            
         !Record 3.2h
         !Record 3.2
         !Record 3.3a
         !Record 3.3b
         !Record 3.4
         !Record 3.5
         !Record 3.6 (3.6.1...3.6.N)
         !Record 3.7
         !Record 3.7.1
         !Record 3.8
         !Record 3.8.1...3.8.n
         if (IATM ==1) then
         
            !Recordd 3.1
            !>>> MODEL,ITYPE,IBMAX_B,n_zero,NOPRNT,NMOL,IPUNCH,  &
            !    IFXTYP,MUNITS,RE,HSPACE,XVBAR,dumrd,sref_lat
            call get_rec_3_1()

            MODEL   = tape5%rec_3_1%MODEL
            ITYPE   = tape5%rec_3_1%ITYPE
            IBMAX_B = tape5%rec_3_1%IBMAX_B

            if (ITYPE==1) then !Horizontal path

               !Recorde 3.2h
               !>>> H1F,RANGEF
               call get_rec_3_2h()

               !Recorde 3.3 not used for this case
               
            elseif (ITYPE==2 .or. ITYPE==3) then

               !Record 3.2
               !>>> H1F,H2F,ANGLEF,RANGEF,BETAF,LENF,HOBS
               call get_rec_3_2()

            endif


            if (IBMAX_B==0) then

               !Record 3.3a
               !>>> AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2
               call get_rec_3_3a()

            else !IBMAX!=0

               !Record 3.3b
               !>>> (layerBoundaries(IB),IB=1,IBMAX)
               call get_rec_3_3b()

            endif


            if (MODEL==0) then !user model atmos.
               
               !Record 3.4
               !>>> IMMAX_B,HMOD
               call get_rec_3_4()

               !Record 3.5
               !Record 3.6 (3.6.1 ... 3.6.N)
               !>>> ZMDL(IM),PM(IM),TM(IM),JCHARP(IM),
               !    JCHART(IM),JLONG(IM),(JCHAR(K,IM),K=1,MXMOL)
               !>>> WMOL(K,IM),K=1,NMOL
               call get_rec_3_5_and_3_6() !USER MODEL 

            endif
      

            if (IXSECT==1) then

               !Record 3.7
               !READ IN THE NUMBER OF MOLECULES IXMOLS, AND THE FLAG IPRFL        
               !INDICATING WHETHER A STANDARD PROFILE (0) OR A USER-INPUT PROFILE 
               !(1) WILL BE USED.
               !>>> IXMOLS,IPRFL,IXSBIN
               call get_rec_3_7()
         
               IPRFL = tape5%rec_3_7%IPRFL
         
               !Record 3.7.1
               !READ IN THE NAMES OF THE MOLECULES                                
               !>>> XSNAME(IXMOLS)
               call get_rec_3_7_1()

               if (IPRFL==0) then !user x-section profile
                  !Record 3.8
                  !A USER-INPUT PROFILE HAS BEEN SELECTED. READ IN THE PROFILES      
                  !AND INTERPOLATE THEM TO THE LEVELS ZMDL.                          
                  !                                                                   
                  !READ IN THE PROFILES. NOTE THAT ZORP CAN BE EITHER ALTITUDE       
                  !OR PRESSURE, DEPENDING UPON THE VALUE OF IZORP: 0 FOR             
                  !ALTITUDE, 1 FOR PRESSURE.
                  !>>> LAYX,IZORP,XTITLE
                  call get_rec_3_8()


                  !Record 3.8.1 .. 3.8.N
                  !>>> ZX(L),(JCHAR(I,L),I=1,IXMOLS) 
                  !    (DTMP(K,L),K=1,IXMOLS) 
                  call get_rec_3_8_1()
         
               endif

            endif
         
         endif !IATM==1
                  

         if (ISOTPL ==1 ) then !v12.7, read isotopologue infomation

            !Record for isotopologue
            !>>> NISOTPL
            !>>> ISOTPL_HCODE(I)
            !>>> IFRMI, NLAYIS
            !>>> ISOTPL_AMNT(K,I,L) !returns ISOTPL_AMNT_tape5(I,L)
            call get_rec_isotpl()
         
         endif
         
         
         if (IAERSL==5) then

            STOP '--- Aerosol and cloud options are not yet implemented.'

            !Record Abs.1
            !Record Abs.2
            !Record Abs.3
            !Record Abs.4.1
            !Record Abs.4.2
         elseif (IAERSL==1 .or. IAERSL==7) then            

            STOP '--- Aerosol and cloud options are not yet implemented.'

            !Record 4.1
            !if any(ICLD==[18,19,20]) then
            !   !Record 4.2
            !end            
            !if (IVSA==1) then
            !   !Record 4.3
            !end
            !if (IAERSL==7) then
            !   !Record 4.4
            !   !Record 4.5
            !end
            !if (IHAZE==7 .or. ICLD==11) then
            !   !Record 4.6.1
            !   !Record 4.6.2
            !   !Record 4.6.3
            !end         
         endif !if (IAERSL==5) then
         
         
                  
         if (ILAS>2) then
            STOP '--- LASER option has been removed.'
            !Record 5.1 ... 5.N
         endif



         !---
         !Record 6
         if ( (IMRG>=13 .and. IMRG<=18) .OR. &
               IMRG==35 .OR. IMRG==36 ) then
            !READ IN CONTROL PARAMETERS for scnmrg
            !>>> HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,IUNIT,&
            !    IFILST,NIFILS,JUNIT,NPTS 
            call get_rec_6()
         endif

         
         !---
         !Record 6.1
         if ( any(IMRG==[ 35,36,45,46 ]) ) then
            !READ IN DIRECTION COSINE 
            !>>> DIRCOS
            call get_rec_6_1()                                           
         endif


         !---
         !Record 7.1
         !Record 7.2 ... 7.N
         if (IMRG>=23 .and. IMRG<=28 ) then
            
            !Record 7.1; for fltmrg
            !>>> V1F,DVF,NPTF,JEMIT,NNFILE,HEDDR
            call get_rec_7_1()
            
            !Record 7.2, 7.3...7.N
            !>>> CVAR
            !>>> XF(I),I=1,NPTS
            call get_rec_7_2_and_7_3()
            
            
         endif
         
         
         !---         
         !Record 8.1
         !Record 9.1
         !Record 10.1
         !Record 10.2
         if (ISCAN==1) then
            
            !Record 8.1
            !READ IN CONTROL PARAMETERS for scanfn
            !>>> HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,IUNIT, &
            !    IFILST,NIFILS,JUNIT,NPTS  
            call get_rec_8_1()

         elseif (ISCAN==2) then

            !Record 9.1
            !for interpolation
            !>>> DVO,V1,V2,JEMIT,I4PT,IUNIT,IFILST,NIFILS,JUNIT,NPTS   
            call get_rec_9_1()

         elseif (ISCAN==3) then
            
            !Record 10.1
            !only one iteration is allowed now.
            !>>> HWHM, V1, V2, JEMIT, JFNIN, MRATIN,          &
            !    DVOUT,IUNIT,IFILST,NIFILS, JUNIT, IVX, NOFIX,param
            call get_rec_10_1()
            
            JFNIN = tape5%rec_10_1%JFNIN
   
            if (abs(JFNIN)>=11) then

               !Record 10.2
               !Functions 11 and above require further parameters 
               !>>> PARM1,PARM2,PARM3 
               call get_rec_10_2()
               
            endif

         endif
         

         !---        
         !Record 11.1         
         !for FILTFN
         if (IFILTR==1) then

            !Record 11.1
            !>>> V1F,DVF,NPTF,JEMIT,IUNIT,IFILST,NIFILS,junit, HEDDR
            call get_rec_11_1()

            !Record 11.2, 11.3...11.N
            !>>> CVAR
            !>>> XF(I),I=1,NPTS
            call get_rec_11_2_and_11_3()

         endif
         
         
         !--- PLOT OPTIONS
         !Record 12.1
         !Record 12.2a
         !Record 12.3a
         !Record 12.2.1a
         !Record 12.2b
         !Record 12.2.1b
         !Record 12.2.2b
         !Reocrd 12.2.3b
         
      !ENDDO !do while
      CLOSE(IRD)


   RETURN
   
   CONTAINS !-------------- Internal Contains ---------------



   !--------------------------------------
   subroutine get_rec_1_1()
   !--------------------------------------
      character(8) :: XID(10)

      do while (.true.)
         read(IRD,'(10(A8))') XID

         if (XID(1)(1:1)/='$') CYCLE !header lines
         if (XID(1)(1:1)=='%') then
            endOfTape5 = .TRUE.
            RETURN
         endif

         EXIT
      enddo
      
      tape5%rec_1_1%XID = XID
   end subroutine

      
   !--------------------------------------
   subroutine get_rec_1_2()
   !--------------------------------------
      character(1) :: CMRG(2)
      integer      :: IMRG !CMRG-->IMRG

      !FORMAT (10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4) !v12.2           
      !fmt = '(10(4X,I1),3X,2A1,3(4X,I1),1X,I4,1X,I4,6X,I4)'
      !fmt = '(10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4)' !v12.2
      !925 FORMAT (10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4,4X,I1,4x,I1,4x,I1) 
      fmt = '(10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4,4X,I1,4x,I1,4x,I1)' !v12.7
      read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) &
         tape5%rec_1_2%IHIRAC ,&
         tape5%rec_1_2%ILBLF4 ,&
         tape5%rec_1_2%ICNTNM ,&
         tape5%rec_1_2%IAERSL ,&
         tape5%rec_1_2%IEMIT_tape5 ,&
         tape5%rec_1_2%ISCAN ,&
         tape5%rec_1_2%IFILTR ,&
         tape5%rec_1_2%IPLOT ,&
         tape5%rec_1_2%ITEST ,&
         tape5%rec_1_2%IATM ,&
         CMRG ,& 
         tape5%rec_1_2%ILAS ,&
         tape5%rec_1_2%IOD ,&
         tape5%rec_1_2%IXSECT ,&
         tape5%rec_1_2%IRAD, &
         tape5%rec_1_2%MPTS ,&
         tape5%rec_1_2%NPTS, &
         ! ,& tape5%rec_1_2%ISPD
         tape5%rec_1_2%ISOTPL, &   !v12.7
         tape5%rec_1_2%IBRD        !v12.7

      IF (CMRG(2).EQ.'A') THEN 
         IMRG = 12 
      ELSEIF (CMRG(2).EQ.'B') THEN 
         IMRG = 22 
      ELSEIF (CMRG(2).EQ.'C') THEN 
         IMRG = 32 
      ELSE 
         READ (CMRG(2),'(I1)') IMRG 
         IF (CMRG(1).EQ.'1') IMRG = IMRG+10 
         IF (CMRG(1).EQ.'2') IMRG = IMRG+20 
         IF (CMRG(1).EQ.'3') IMRG = IMRG+30 
         IF (CMRG(1).EQ.'4') IMRG = IMRG+40 
      ENDIF
      
      tape5%rec_1_2%IMRG = IMRG
      
   end subroutine


   !--------------------------------------
   subroutine get_rec_1_2a()         
   !--------------------------------------
      !real :: XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL 
      read(IRD,*, iostat=iostat, iomsg=iomsg ) tape5%rec_1_2a
   end subroutine                 


   !--------------------------------------
   subroutine get_rec_1_2_1()
   !--------------------------------------
      !!integer :: INFLAG,IOTFLG,JULDAT      
      !
      !!1010 FORMAT (2I5,2X,I3) 
      !fmt = '(2I5,2X,I3)'
      !read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_1_2_1
      
      integer :: INFLAG,IOTFLG,JULDAT,ISOLVAR
      real    :: SCON,SOLCYCFRAC,SOLVAR
      
      !1010 FORMAT (2I5,2X,I3,I5,f10.4,5x,f10.4,5x,2f10.5)
      fmt = '(2I5,2X,I3,I5,f10.4,5x,f10.4,5x,2f10.5)'
      read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_1_2_1            
   end subroutine            


   !--------------------------------------
   subroutine get_rec_1_3()
   !--------------------------------------
      !real(r8) :: V1,V2
      !real     :: SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,DPTFAC
      !integer  :: ILNFLG
      !real     :: DVOUT
      !integer  :: nmol_scal, n_xs_scal
      
      !970 FORMAT (8E10.3,4X,I1,5x,e10.3,3X,i2,3x,i2) 
      fmt = '(8E10.3,4X,I1,5x,e10.3,i5,i5)'
      read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_1_3
   end subroutine


   !--------------------------------------
   subroutine get_rec_1_3a_and_1_3b()
   !--------------------------------------
      character(1), allocatable :: hmol_scal(:)
      real,         allocatable :: xmol_scal(:)
      character(1), allocatable :: h_xs_scal(:)
      real,         allocatable :: x_xs_scal(:)

      integer :: nmol_scal, n_xs_scal
      integer :: i,j,k

      nmol_scal = tape5%rec_1_3%nmol_scal
      n_xs_scal = tape5%rec_1_3%n_xs_scal

      if (nmol_scal>0) then
         allocate(hmol_scal(nmol_scal))
         allocate(xmol_scal(nmol_scal))

         !972 FORMAT (64a1) 
         fmt = '(64a1)'
         read(IRD, trim(fmt),iostat=iostat, iomsg=iomsg ) &
            ( hmol_scal(i),i=1,nmol_scal )

         !973 FORMAT (8e15.7)            
         fmt = '(8e15.7)'
         read(IRD, trim(fmt),iostat=iostat, iomsg=iomsg ) &
            ( xmol_scal(i),i=1,nmol_scal )

         call move_alloc(hmol_scal, tape5%rec_1_3a_and_1_3b%hmol_scal)
         call move_alloc(xmol_scal, tape5%rec_1_3a_and_1_3b%xmol_scal)
      endif
      
      if (n_xs_scal>0) then 
         allocate(h_xs_scal(n_xs_scal))
         allocate(x_xs_scal(n_xs_scal))

         !972 FORMAT (64a1) 
         fmt = '(64a1)'
         read(IRD, trim(fmt),iostat=iostat, iomsg=iomsg ) &
            ( h_xs_scal(i),i=1,n_xs_scal )
   
         !973 FORMAT (8e15.7)
         fmt = '(8e15.7)'
         read(IRD, trim(fmt),iostat=iostat, iomsg=iomsg ) &
            ( x_xs_scal(i),i=1,n_xs_scal )

         call move_alloc(h_xs_scal, tape5%rec_1_3a_and_1_3b%h_xs_scal)
         call move_alloc(x_xs_scal, tape5%rec_1_3a_and_1_3b%x_xs_scal)
      endif
      
   end subroutine
   

   !--------------------------------------
   subroutine get_rec_1_4()
   !--------------------------------------
      !real         :: TMPBND, BNDEMI(3), BNDRFL(3)
      !character(1) :: surf_refl
      
      !971 FORMAT (7E10.3,4X,A1) 
      fmt = '(7E10.3,4X,A1)'
      read(IRD, trim(fmt),iostat=iostat, iomsg=iomsg ) tape5%rec_1_4

   end subroutine         


   !--------------------------------------
   subroutine get_rec_1_5()
   !--------------------------------------
      !integer :: NSPCRT
      
      !1015 FORMAT (3x,I2) 
      fmt = '(3x,I2)'
      read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_1_5
   end subroutine


   !--------------------------------------
   subroutine get_rec_1_6a()            
   !--------------------------------------
      !character(55) :: PATH1
      !integer       :: laytot
      
      !945 FORMAT (A55,1X,I4,2I5)
      fmt = '(A55,1X,I4,2I5)'
      read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_1_6a
   end subroutine


   !--------------------------------------
   subroutine get_rec_2_1()
   !--------------------------------------
      !integer       :: IFORM,NLAYRS,NMOL
      !real          :: SECNT0
      !character(20) :: HEAD20
      !real          :: ZH1
      !character(4)  :: HEAD4
      !real          :: ZH2
      !character(5)  :: HEAD5
      !character(8)  :: hol_angle
      !character(7)  :: HEAD7
      
      !901 FORMAT (1X,I1,I3,I5,F10.2,A20,F8.2,A4,F8.2,A5,A8,A7)
      fmt = '(1X,I1,I3,I5,F10.2,A20,F8.2,A4,F8.2,A5,A8,A7)'
      read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_2_1
   end subroutine


   !--------------------------------------
   subroutine get_rec_2_1_1_to_2_1_3()
   !--------------------------------------
      real,         allocatable :: PAVE(:),TAVE(:),SECNTK(:)
      character(3), allocatable :: CINP(:)
      integer,      allocatable :: IPTHRK(:)
      real,         allocatable :: ALTZ(:),PZ(:),TZ(:)
      real,         allocatable :: WKL(:,:), WBRODL(:)

      integer :: IFORM,NLAYRS, NMOL,NMOL1, L,M
      IFORM  = tape5%rec_2_1%IFORM
      NLAYRS = tape5%rec_2_1%NLAYRS
      NMOL   = tape5%rec_2_1%NMOL

      NMOL1 = NMOL
      if (NMOL1>7) NMOL1 = 7
      
      allocate( PAVE(  NLAYRS),&
                TAVE(  NLAYRS),&
                SECNTK(NLAYRS),&
                CINP(  NLAYRS),&
                IPTHRK(NLAYRS),&
                ALTZ(0:NLAYRS),&
                PZ(  0:NLAYRS),&
                TZ(  0:NLAYRS) )

      allocate( WKL(NMOL,NLAYRS), WBRODL(NLAYRS) )
 
      do L = 1,NLAYRS

          !%---Record 2.1.1
          IF (L==1) THEN 
            IF (IFORM==1) THEN 
               !910 FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2)) 
               fmt='(E15.7,2(F10.4),A3,I2,1X,2(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg )&
                         PAVE(L),TAVE(L),SECNTK(L),CINP(L),IPTHRK(L),&
                         ALTZ(L-1),PZ(L-1),TZ(L-1),&
                         ALTZ(L),PZ(L),TZ(L)                      
            ELSE
               !911 FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))
               fmt='(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg )&
                         PAVE(L),TAVE(L),SECNTK(L),CINP(L),IPTHRK(L),&
                         ALTZ(L-1),PZ(L-1),TZ(L-1),&
                         ALTZ(L),PZ(L),TZ(L)                      
            ENDIF 
         ELSE 
            IF (IFORM==1) THEN 
               !915 FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
               fmt='(E15.7,2(F10.4),A3,I2,23X,(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg )&
                         PAVE(L),TAVE(L),SECNTK(L),CINP(L),IPTHRK(L),&
                         ALTZ(L),PZ(L), TZ(L)                                                
            ELSE 
               !916 FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
               fmt='(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg )&
                         PAVE(L),TAVE(L),SECNTK(L),CINP(L),IPTHRK(L),&
                         ALTZ(L),PZ(L), TZ(L)                                                
            ENDIF
         ENDIF


          !%---Record 2.1.2 and 2.1.3
          IF (IFORM.EQ.1) THEN
            !925 FORMAT (8E15.7)
            fmt='(8E15.7)'
            !READ (IRD,trim(fmt)) (WKL(M,L),M=1,7),WBRODL(L) 
            READ (IRD,trim(fmt)) (WKL(M,L),M=1,NMOL1),WBRODL(L) 
            IF (NMOL.GT.7) &
               READ (IRD,trim(fmt)) (WKL(M,L),M=8,NMOL) 
         ELSE
            !927 FORMAT (8E10.3)
            fmt='(8E10.3)'
            !READ (IRD,trim(fmt)) (WKL(M,L),M=1,7),WBRODL(L) 
            READ (IRD,trim(fmt)) (WKL(M,L),M=1,NMOL1),WBRODL(L) 
            IF (NMOL.GT.7) &
               READ (IRD,trim(fmt)) (WKL(M,L),M=8,NMOL) 
         ENDIF 

      enddo
      
      call move_alloc(PAVE    ,tape5%rec_2_1_1_to_2_1_3%PAVE   )
      call move_alloc(TAVE    ,tape5%rec_2_1_1_to_2_1_3%TAVE   )
      call move_alloc(SECNTK  ,tape5%rec_2_1_1_to_2_1_3%SECNTK )
      call move_alloc(CINP    ,tape5%rec_2_1_1_to_2_1_3%CINP   )
      call move_alloc(IPTHRK  ,tape5%rec_2_1_1_to_2_1_3%IPTHRK )
      call move_alloc(ALTZ    ,tape5%rec_2_1_1_to_2_1_3%ALTZ   )
      call move_alloc(PZ      ,tape5%rec_2_1_1_to_2_1_3%PZ     )
      call move_alloc(TZ      ,tape5%rec_2_1_1_to_2_1_3%TZ     )

      call move_alloc(WKL,   tape5%rec_2_1_1_to_2_1_3%WKL    )
      call move_alloc(WBRODL,tape5%rec_2_1_1_to_2_1_3%WBRODL )

   end subroutine 


   !--------------------------------------
   subroutine get_rec_2_2()
   !--------------------------------------
      !integer :: IXMOLS,IXSBIN
      
      !930 FORMAT (I5,5X,I5)
      fmt='(I5,5X,I5)'
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_2_2
   end subroutine


   !--------------------------------------
   subroutine get_rec_2_2_1()
   !--------------------------------------
      character(10), allocatable :: XSNAME(:)

      integer :: IXMOLS, i,j,k
      IXMOLS = tape5%rec_2_2%IXMOLS

      allocate(XSNAME(IXMOLS))
      
      IF (IXMOLS.GT.7) THEN 
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7) 
         READ (IRD,'(8A10)') (XSNAME(I),I=8,IXMOLS) 
      ELSE 
         READ (IRD,'(7A10)') (XSNAME(I),I=1,IXMOLS) 
      ENDIF

      call move_alloc(XSNAME, tape5%rec_2_2_1%XSNAME )
   end subroutine


   !--------------------------------------
   subroutine get_rec_2_2_2()
   !--------------------------------------
      !integer       :: IFRMX,NLAYXS,IXMOL
      !real          :: SECNTX
      !character(60) :: HEDXS
      
      !900 FORMAT (1X,I1,I3,I5,F10.2,15A4)
      fmt='(1X,I1,I3,I5,F10.2,15A4)'
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_2_2_2
   end subroutine


   !--------------------------------------
   subroutine get_rec_2_2_3_to_2_2_5()
   !--------------------------------------
      !--- ALTZX,PZX,TZX names are different from LBLRTM
      real,         allocatable :: PAVX(:),TAVX(:),SECKXS(:)
      character(3), allocatable :: CINPX(:)
      integer,      allocatable :: IPTHKX(:)
      real,         allocatable :: ALTZX(:),PZX(:),TZX(:)
      real,         allocatable :: XAMNT(:,:),WBRODX(:)

      integer :: IFRMX,NLAYXS, IXMOL,IXMOL1, L,M 
      IFRMX  = tape5%rec_2_2_2%IFRMX
      NLAYXS = tape5%rec_2_2_2%NLAYXS
      IXMOL  = tape5%rec_2_2_2%IXMOL

      IXMOL1 = IXMOL
      if (IXMOL1>7) IXMOL1=7
      
      
      allocate( PAVX(   NLAYXS),&
                TAVX(   NLAYXS),&
                SECKXS( NLAYXS),&
                CINPX(  NLAYXS),&
                IPTHKX( NLAYXS),&
                ALTZX(0:NLAYXS),&
                PZX(  0:NLAYXS),&
                TZX(  0:NLAYXS) )

      allocate( XAMNT(IXMOL,NLAYXS), WBRODX(NLAYXS) )
      
      do L = 1, NLAYXS
      
         !%---Record 2.2.3
         IF (L.EQ.1) THEN 
            IF (IFRMX.EQ.1) THEN
               !910 FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))
               fmt='(E15.7,2(F10.4),A3,I2,1X,2(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt),iostat=iostat,iomsg=iomsg)&
                  PAVX(L),TAVX(L),SECKXS(L),CINPX(L),IPTHKX(L),&
                  ALTZX(L-1),PZX(L-1),TZX(L-1),&
                  ALTZX(L),PZX(L),TZX(L)
            ELSE
               !911 FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))
               fmt='(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt),iostat=iostat,iomsg=iomsg)&
                  PAVX(L),TAVX(L),SECKXS(L),CINPX(L),IPTHKX(L),&
                  ALTZX(L-1),PZX(L-1),TZX(L-1),&
                  ALTZX(L),PZX(L),TZX(L)                             
            ENDIF 
         ELSE 
            IF (IFRMX.EQ.1) THEN
               !915 FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
               fmt='(E15.7,2(F10.4),A3,I2,23X,(F7.2,F8.3,F7.2))'
               READ (IRD,trim(fmt),iostat=iostat,iomsg=iomsg)&
                  PAVX(L),TAVX(L),SECKXS(L),CINPX(L),IPTHKX(L),&
                  ALTZX(L),PZX(L),TZX(L)
            ELSE
               !916 FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
               fmt='(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))' 
               READ (IRD,trim(fmt),iostat=iostat,iomsg=iomsg)&
                  PAVX(L),TAVX(L),SECKXS(L),CINPX(L),IPTHKX(L),&
                  ALTZX(L),PZX(L),TZX(L)
            ENDIF 
         ENDIF 

         !%---Record 2.2.4 and 2.2.5            
         IF (IFRMX.EQ.1) THEN 
            !925 FORMAT (8E15.7) 
            !READ (IRD,'(8E15.7)') (XAMNT(M,L),M=1,7),WBRODX(L) 
            READ (IRD,'(8E15.7)') (XAMNT(M,L),M=1,IXMOL1),WBRODX(L) 
            IF (IXMOL.GT.7)&
               READ (IRD,'(8E15.7)') (XAMNT(M,L),M=8,IXMOL) 
         ELSE
            !927 FORMAT (8E10.3)         
            !READ (IRD,'(8E10.3)') (XAMNT(M,L),M=1,7),WBRODX(L) 
            READ (IRD,'(8E10.3)') (XAMNT(M,L),M=1,IXMOL1),WBRODX(L) 
            IF (IXMOL.GT.7)&
               READ (IRD,'(8E10.3)') (XAMNT(M,L),M=8,IXMOL) 
         ENDIF
      
      enddo
      
      call move_alloc(PAVX    ,tape5%rec_2_2_3_to_2_2_5%PAVX   )
      call move_alloc(TAVX    ,tape5%rec_2_2_3_to_2_2_5%TAVX   )
      call move_alloc(SECKXS  ,tape5%rec_2_2_3_to_2_2_5%SECKXS )
      call move_alloc(CINPX   ,tape5%rec_2_2_3_to_2_2_5%CINPX  )
      call move_alloc(IPTHKX  ,tape5%rec_2_2_3_to_2_2_5%IPTHKX )
      call move_alloc(ALTZX   ,tape5%rec_2_2_3_to_2_2_5%ALTZX  )
      call move_alloc(PZX     ,tape5%rec_2_2_3_to_2_2_5%PZX    )
      call move_alloc(TZX     ,tape5%rec_2_2_3_to_2_2_5%TZX    )

      call move_alloc(XAMNT   ,tape5%rec_2_2_3_to_2_2_5%XAMNT  )
      call move_alloc(WBRODX  ,tape5%rec_2_2_3_to_2_2_5%WBRODX )

   end subroutine
   
   
   
   !--------------------------------------
   subroutine get_rec_3_1()
   !--------------------------------------
      !integer       :: MODEL,ITYPE,IBMAX_B,n_zero,NOPRNT,NMOL,IPUNCH,&
      !                 IFXTYPE,MUNITS
      !real          :: RE,HSPACE,XVBAR,dumrd
      !character(10) :: sref_lat
      
      !900 FORMAT (7I5,I2,1X,I2,4F10.3,A10)
      fmt = '(7I5,I2,1X,I2,4F10.3,A10)'
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_3_1
   end subroutine


   !--------------------------------------
   subroutine get_rec_3_2h()
   !--------------------------------------
      !real :: H1F,RANGEF
      
      !910 FORMAT (F10.3,10X,10X,F10.3)
      READ (IRD,'(F10.3,10X,10X,F10.3)') tape5%rec_3_2h
   end subroutine               


   !--------------------------------------
   subroutine get_rec_3_2()
   !--------------------------------------
      !real    :: H1F,H2F,ANGLEF,RANGEF,BETAF
      !integer :: LENF
      !real    :: HOBS
      
      !932 FORMAT (5F10.4,I5,5X,F10.4)
      fmt='(5F10.4,I5,5X,F10.4)' 
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_3_2
   end subroutine                  


   !--------------------------------------
   subroutine get_rec_3_3a()
   !--------------------------------------
      !real :: AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2
      
      !936 FORMAT (5F10.3)
      READ (IRD,'(5F10.3)', iostat=iostat, iomsg=iomsg ) tape5%rec_3_3a
   end subroutine


   !--------------------------------------
   subroutine get_rec_3_3b()
   !--------------------------------------
      real, allocatable :: BND(:)
      
      integer :: IBMAX_B, IBMAX, IB
      IBMAX_B = tape5%rec_3_1%IBMAX_B
      IBMAX = abs(IBMAX_B)

      allocate(BND(IBMAX))
      !940 FORMAT (8F10.3)
      READ (IRD,'(8F10.3)') (BND(IB),IB=1,IBMAX)
      call move_alloc(BND,tape5%rec_3_3b%layerBoundaries)

   end subroutine


   !--------------------------------------
   subroutine get_rec_3_4()
   !--------------------------------------
      !integer      :: IMMAX_B
      !character(8) :: HMOD(3)
      
      !905 FORMAT (I5,3A8)
      READ (IRD,'(I5,3A8)',iostat=iostat, iomsg=iomsg) tape5%rec_3_4
   end subroutine               


   !--------------------------------------
   subroutine get_rec_3_5_and_3_6() !USER MODEL 
   !--------------------------------------
      USE Module_ConstParam, ONLY: MXMOL 
      
      real,         allocatable :: ZMDL(:),PM(:),TM(:)
      character(1), allocatable :: JCHARP(:),JCHART(:),&
                                   JLONG(:),JCHAR(:,:)
      real,         allocatable :: WMOL(:,:)   
      
     
      integer :: IMMAX_B, IMMAX, NMOL, IM,K
      NMOL    = tape5%rec_3_1%NMOL
      IMMAX_B = tape5%rec_3_4%IMMAX_B
      IMMAX = abs(IMMAX_B)
      
      !--- pls be noted that the data dimensions are different from LBLRTM
      allocate( ZMDL(  IMMAX),&
                PM(    IMMAX),&
                TM(    IMMAX),&
                JCHARP(IMMAX),&
                JCHART(IMMAX),&
                JLONG( IMMAX),&
                JCHAR( MXMOL,IMMAX) );

      allocate(WMOL(NMOL,IMMAX))

      do IM = 1, IMMAX 

         !Record 3.5
         !900 FORMAT (3E10.3,5X,2A1,1X,A1,1X,39A1)         
         !900 FORMAT (3E10.3,5X,2A1,1X,A1,1X,47A1) !v12.7
         fmt = '(3E10.3,5X,2A1,1X,A1,1X,47A1)'
         READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg )&
            ZMDL(IM),PM(IM),TM(IM),JCHARP(IM),&
            JCHART(IM),JLONG(IM),(JCHAR(K,IM),K=1,MXMOL)

         !Record 3.6
         IF (JLONG(IM).EQ.'L') THEN          
            !906 FORMAT (8E15.8)
            READ (IRD,'(8E15.8)') (WMOL(K,IM),K=1,NMOL) 
         ELSEIF (JLONG(IM).EQ.' ') THEN
            !905 FORMAT (8E10.3)         
            READ (IRD,'(8E10.3)') (WMOL(K,IM),K=1,NMOL) 
         ENDIF 
   
      enddo
      
      call move_alloc( ZMDL     ,tape5%rec_3_5_and_3_6%ZMDL  )
      call move_alloc( PM       ,tape5%rec_3_5_and_3_6%PM    )
      call move_alloc( TM       ,tape5%rec_3_5_and_3_6%TM    )
      call move_alloc( JCHARP   ,tape5%rec_3_5_and_3_6%JCHARP)
      call move_alloc( JCHART   ,tape5%rec_3_5_and_3_6%JCHART)
      call move_alloc( JLONG    ,tape5%rec_3_5_and_3_6%JLONG )
      call move_alloc( JCHAR    ,tape5%rec_3_5_and_3_6%JCHAR )

      call move_alloc( WMOL, tape5%rec_3_5_and_3_6%WMOL )

   end subroutine


   !--------------------------------------
   subroutine get_rec_3_7()
   !--------------------------------------
      !integer :: IXMOLS,IPRFL,IXSBIN
      
      !905 FORMAT (3I5)
      READ (IRD,'(3I5)', iostat=iostat, iomsg=iomsg) tape5%rec_3_7
   end subroutine


   !--------------------------------------
   subroutine get_rec_3_7_1()
   !--------------------------------------
      character(10), allocatable :: XSNAME(:)
      
      integer :: IXMOLS, I
      IXMOLS = tape5%rec_3_7%IXMOLS
      
      allocate( XSNAME(IXMOLS) )
      
      IF (IXMOLS.GT.7) THEN 
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7) 
         READ (IRD,'(8A10)') (XSNAME(I),I=8,IXMOLS) 
      ELSE 
         READ (IRD,'(7A10)') (XSNAME(I),I=1,IXMOLS) 
      ENDIF 
      
      call move_alloc(XSNAME, tape5%rec_3_7_1%XSNAME)

   end subroutine


   !--------------------------------------
   subroutine get_rec_3_8()
   !--------------------------------------
      !integer       :: LAYX,IZORP
      !character(50) :: XTITLE
      
      !905 FORMAT (2I5,A)
      READ (IRD, '(2I5,A)', iostat=iostat, iomsg=iomsg) tape5%rec_3_8
   end subroutine                  


   !--------------------------------------
   subroutine get_rec_3_8_1()
   !--------------------------------------
      real,         allocatable :: ZX(:)
      real,         allocatable :: PX(:)
      character(1), allocatable :: JCHAR(:,:)
      real,         allocatable :: DTMP(:,:)

      integer :: IZORP,LAYX,IXMOLS, i,j,k,l
      IZORP  = tape5%rec_3_8%IZORP
      LAYX   = tape5%rec_3_8%LAYX
      IXMOLS = tape5%rec_3_7%IXMOLS
      

      IF (IZORP .EQ. 0) THEN 
         allocate( ZX(LAYX),&
                   JCHAR(IXMOLS,LAYX),&
                   DTMP( IXMOLS,LAYX) )

         DO L = 1, LAYX
            !925 FORMAT (F10.3,5X,38A1)          
            fmt='(F10.3,5X,38A1)'
            READ (IRD,trim(fmt)) ZX(L),(JCHAR(I,L),I=1,IXMOLS)
            !935 FORMAT (8E10.3)
            fmt='(8E10.3)' 
            READ (IRD,trim(fmt)) (DTMP(K,L),K=1,IXMOLS) 
         ENDDO
         
         call move_alloc( ZX    ,tape5%rec_3_8_1%ZX    )
         call move_alloc( JCHAR ,tape5%rec_3_8_1%JCHAR )
         call move_alloc( DTMP  ,tape5%rec_3_8_1%DTMP  )

      ELSE 
         allocate( PX(LAYX),&
                   JCHAR(IXMOLS,LAYX),&
                   DTMP( IXMOLS,LAYX) )

         DO L = 1, LAYX
            !925 FORMAT (F10.3,5X,38A1)         
            fmt='(F10.3,5X,38A1)'
            READ (IRD,trim(fmt)) PX(L),(JCHAR(I,L),I=1,IXMOLS)
            !935 FORMAT (8E10.3)            
            fmt='(8E10.3)' 
            READ (IRD,trim(fmt)) (DTMP(K,L),K=1,IXMOLS) 
         ENDDO

         call move_alloc( PX    ,tape5%rec_3_8_1%PX    )
         call move_alloc( JCHAR ,tape5%rec_3_8_1%JCHAR )
         call move_alloc( DTMP  ,tape5%rec_3_8_1%DTMP  )

      ENDIF 

   end subroutine                  

   
   !--------------------------------------
   subroutine get_rec_isotpl()
   !--------------------------------------
      !USE Module_ConstParam, ONLY: MXMOL, MXISOTPL

      integer              :: I,L
      integer              :: NISOTPL
      integer ,allocatable :: ISOTPL_HCODE(:)
      integer              :: IFRMI,NLAYIS
      real    ,allocatable :: ISOTPL_AMNT_tape5(:,:)
   
        
      fmt='(10I5)' !928 FORMAT (10I5)
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) NISOTPL

      allocate(ISOTPL_HCODE(NISOTPL))     
      fmt='(10I5)' !928 FORMAT (10I5)
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) &
                (ISOTPL_HCODE(I),I=1,NISOTPL)
      
      !900 FORMAT (1X,I1,I3,I5,F10.2,15A4) 
      fmt='(1X,I1,I3,I5,F10.2,15A4)'
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) IFRMI,NLAYIS
            
      !allocate( ISOTPL_AMNT( MXMOL,MXISOTPL,MXLAY ) )
      allocate( ISOTPL_AMNT_tape5( NISOTPL,NLAYIS ) ) !This reader simple transfer the data to the program. 
      do L = 1, NLAYIS
      
         IF (IFRMI.EQ.1) THEN
            !925 FORMAT (8E15.7)
            READ (IRD,'(8E15.7)', iostat=iostat, iomsg=iomsg ) &
                      (ISOTPL_AMNT_tape5(I,L),I=1,NISOTPL)
         ELSE
            !927 FORMAT (8E10.3)
            READ (IRD,'(8E10.3)', iostat=iostat, iomsg=iomsg ) &
                      (ISOTPL_AMNT_tape5(I,L),I=1,NISOTPL)
         ENDIF
                  
      enddo
      

      tape5%rec_isotpl%NISOTPL = NISOTPL
      tape5%rec_isotpl%IFRMI   = IFRMI
      tape5%rec_isotpl%NLAYIS  = NLAYIS
      
      call move_alloc(ISOTPL_HCODE      ,tape5%rec_isotpl%ISOTPL_HCODE )
      call move_alloc(ISOTPL_AMNT_tape5 ,tape5%rec_isotpl%ISOTPL_AMNT_tape5 )
     
   end subroutine
   

   !--------------------------------------
   subroutine get_rec_6()
   !--------------------------------------
      !real    :: HWHM
      !real(r8):: V1,V2
      !integer :: JEMIT,JFN,JVAR
      !real    :: SAMPL
      !integer :: IUNIT,IFILST,NIFILS,JUNIT,NPTS
      
      !900 FORMAT (3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)      
      fmt='(3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)'
      READ (IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_6
   end subroutine


   !--------------------------------------
   subroutine get_rec_6_1()
   !--------------------------------------
      !real :: DIRCOS
      
      !905 FORMAT (F10.8) 
      !READ (IRD, '(F10.8)', iostat=iostat, iomsg=iomsg) tape5%rec_6_1
      READ (IRD, '(F10.7)', iostat=iostat, iomsg=iomsg) tape5%rec_6_1
   end subroutine 


   !--------------------------------------
   subroutine get_rec_7_1()
   !--------------------------------------
      !real          :: V1F
      !real(r8)      :: DVF
      !integer       :: NPTF,JEMIT,NNFILE
      !character(35) :: HEDDR
      
      !900 FORMAT (2F10.4,I5,I5,I5,10X,8A4,A3)
      fmt='(2F10.4,I5,I5,I5,10X,8A4,A3)'
      READ (IRD,trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_7_1
   end subroutine            


   !--------------------------------------
   subroutine get_rec_7_2_and_7_3()
   !--------------------------------------
      character(80)     :: CVAR
      real, allocatable :: XF(:)
      
      integer :: NPTS, I
      NPTS = tape5%rec_7_1%NPTF !NPTF<0 not implemented!
      
      allocate(XF(NPTS))

      READ (IRD, '(A80)') CVAR
      READ (IRD, trim(CVAR)) (XF(I),I=1,NPTS)

      tape5%rec_7_2_and_7_3%CVAR = CVAR
      call move_alloc(XF, tape5%rec_7_2_and_7_3%XF)      

   end subroutine


   !--------------------------------------
   subroutine get_rec_8_1()
   !--------------------------------------
      !real    :: HWHM
      !real(r8):: V1,V2
      !integer :: JEMIT,JFN,JVAR
      !real    :: SAMPL
      !integer :: IUNIT, IFILST,NIFILS,JUNIT,NPTS
      
      !900 FORMAT (3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)
      fmt='(3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)'
      READ (IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_8_1
   end subroutine


   !--------------------------------------
   subroutine get_rec_9_1()
   !--------------------------------------
      !real    :: DVO
      !real(r8):: V1,V2
      !integer :: JEMIT,I4PT,IUNIT,IFILST,NIFILS,JUNIT,NPTS
      
      !900 FORMAT (3F10.3,2I5,15X,5I5)
      fmt='(3F10.3,2I5,15X,5I5)'
      READ (IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_9_1
   end subroutine


   !--------------------------------------
   subroutine get_rec_10_1()
   !--------------------------------------
      !real    :: HWHM
      !real(r8):: V1,V2
      !integer :: JEMIT,JFNIN,MRATIN
      !real    :: DVOUT
      !integer :: IUNIT,IFILST,NIFILS,JUNIT,IVX,NOFIX
      !real    :: param
      
      !10 Format(3F10.3,3I5,F10.5,4I5,I3,I2,f10.6) 
      fmt='(3F10.3,3I5,F10.5,4I5,I3,I2,f10.6)'
      Read(IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_10_1
   end subroutine


   !--------------------------------------
   subroutine get_rec_10_2()
   !--------------------------------------
      !real :: PARM1,PARM2,PARM3
      Read(IRD, '(3F10.4)', iostat=iostat, iomsg=iomsg ) tape5%rec_10_2
   end subroutine


   !--------------------------------------
   subroutine get_rec_11_1()
   !--------------------------------------
      !real(r8)      :: V1F
      !real          :: DVF
      !integer       :: NPTF,JEMIT,IUNIT,IFILST,NIFILS,junit
      !character(35) :: HEDDR
      
      !900 FORMAT (2F10.4,6I5,8A4,A3)
      fmt='(2F10.4,6I5,8A4,A3)'
      READ (IRD, trim(fmt), iostat=iostat, iomsg=iomsg ) tape5%rec_11_1
   end subroutine


   !--------------------------------------
   subroutine get_rec_11_2_and_11_3()
   !--------------------------------------
      character(80)     :: CVAR
      real, allocatable :: XF(:)
      
      integer :: NPTS, I
      NPTS = tape5%rec_11_1%NPTF !NPTF<0 not yet implemented!
      
      allocate(XF(NPTS))

      READ (IRD, '(A80)') CVAR
      READ (IRD, trim(CVAR)) (XF(I),I=1,NPTS)

      tape5%rec_11_2_and_11_3%CVAR = CVAR
      call move_alloc(XF, tape5%rec_11_2_and_11_3%XF)
   end subroutine

   END SUBROUTINE      

   
   
   
   !-----------------------------------------------------------------------
   ! Inquire function to return IATM to determine whether the input atmospheric data
   ! are layer data (path) or level data (profile). 
   !-----------------------------------------------------------------------
   integer FUNCTION get_IATM( tape5 )
   !-----------------------------------------------------------------------           
      type(CLBLM_TAPE5) ,intent(in) :: tape5
     
      get_IATM = tape5%rec_1_2%IATM
         
   END FUNCTION


   !-----------------------------------------------------------------------
   ! Inquire function to check if the isotopologue data present in 
   ! a TAPE5 scene data.
   !-----------------------------------------------------------------------
   logical FUNCTION isotplPresent( tape5 )
   !-----------------------------------------------------------------------      
      type(CLBLM_TAPE5) ,intent(in) :: tape5
              
      if (tape5%rec_1_2%ISOTPL ==1) then
         isotplPresent = .TRUE.
      else
         isotplPresent = .FALSE.
      endif         
   END FUNCTION

   
   !-----------------------------------------------------------------------
   ! Inquire function to check if the cross-section profile is present in 
   ! a TAPE5 scene data.
   !-----------------------------------------------------------------------
   logical FUNCTION separateXSectDataPresent( tape5 )
   !-----------------------------------------------------------------------
      type(CLBLM_TAPE5) ,intent(in) :: tape5
      
      if (tape5%rec_1_2%IXSECT >0) then
         separateXSectDataPresent = .TRUE.
      else
         separateXSectDataPresent = .FALSE.
      endif
         

   END FUNCTION

   !-----------------------------------------------------------------------
   ! Inquire function to check if the user provided layering grid is present in 
   ! a TAPE5 scene data.
   !-----------------------------------------------------------------------
   logical FUNCTION userRTgridPresent( tape5 )
   !-----------------------------------------------------------------------
      type(CLBLM_TAPE5) ,intent(in) :: tape5
         
      if (tape5%rec_3_1%IBMAX_B /=0) then
         userRTgridPresent = .TRUE.
      else
         userRTgridPresent = .FALSE.
      endif
         

   END FUNCTION
   
END MODULE

