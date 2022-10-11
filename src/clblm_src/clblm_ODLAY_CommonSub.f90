!
! CREATION HISTORY:
!       Modified from LBLRTM v12.9
!       Yingtao Ma, AER@NOAA/NESDIS/STAR
!       yma@aer.com; yingtao.ma@noaa.gov
!


! MOLEC,RADFN,RADFNI and line_exception are shared by clblm_ODLAY.f90 and clblm_LineF4.f90
! XINT is shared by clblm_ODLAY.f90 and clblm_Continuum.f90

MODULE Module_ODLAY_CommonSub 

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC :: & !Subroutines and functions
            MOLEC, &
            RADFN, &
            RADFNI, &
            line_exception
            
   PUBLIC :: & !Constants
            AVRAT, CGAUSS, CF1, CF2, CF3, CER

   
   
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   !      BLOCK DATA VOICON 
   !--------------------------------------------------------------------
   !                                                                       
   !     AVRAT CONTAINS THE PARAMTERS AS A FUNCTION OF ZETA USED TO        
   !     OBTAIN THE VOIGTS' WIDTH FROM THE LORENTZ AND DOPPLER WIDTHS.     
   !                                                                       
   !     COMMON /VOICOM/ AVRAT(102),                                       
   !                     CGAUSS(102),CF1(102),CF2(102),CF3(102),CER(102)   
   !                                                                       
   !      COMMON /VOICOM/ AV01(50),AV51(52),CG01(50),CG51(52),CFA01(50),    &
   !     &                CFA51(52),CFB01(50),CFB51(52),CFC01(50),          &
   !     &                CFC51(52),CER01(50),CER51(52)                     
   !                                                                       
   real ,PARAMETER :: AVRAT(1:102) = [&
     &  .10000E+01,  .99535E+00,  .99073E+00,  .98613E+00,  .98155E+00, &
     &  .97700E+00,  .97247E+00,  .96797E+00,  .96350E+00,  .95905E+00, &
     &  .95464E+00,  .95025E+00,  .94589E+00,  .94156E+00,  .93727E+00, &
     &  .93301E+00,  .92879E+00,  .92460E+00,  .92045E+00,  .91634E+00, &
     &  .91227E+00,  .90824E+00,  .90425E+00,  .90031E+00,  .89641E+00, &
     &  .89256E+00,  .88876E+00,  .88501E+00,  .88132E+00,  .87768E+00, &
     &  .87410E+00,  .87058E+00,  .86712E+00,  .86372E+00,  .86039E+00, &
     &  .85713E+00,  .85395E+00,  .85083E+00,  .84780E+00,  .84484E+00, &
     &  .84197E+00,  .83919E+00,  .83650E+00,  .83390E+00,  .83141E+00, &
     &  .82901E+00,  .82672E+00,  .82454E+00,  .82248E+00,  .82053E+00, &
     !  DATA AV51/                                                       &
     &  .81871E+00,  .81702E+00,  .81547E+00,  .81405E+00,  .81278E+00, &
     &  .81166E+00,  .81069E+00,  .80989E+00,  .80925E+00,  .80879E+00, &
     &  .80851E+00,  .80842E+00,  .80852E+00,  .80882E+00,  .80932E+00, &
     &  .81004E+00,  .81098E+00,  .81214E+00,  .81353E+00,  .81516E+00, &
     &  .81704E+00,  .81916E+00,  .82154E+00,  .82418E+00,  .82708E+00, &
     &  .83025E+00,  .83370E+00,  .83742E+00,  .84143E+00,  .84572E+00, &
     &  .85029E+00,  .85515E+00,  .86030E+00,  .86573E+00,  .87146E+00, &
     &  .87747E+00,  .88376E+00,  .89035E+00,  .89721E+00,  .90435E+00, &
     &  .91176E+00,  .91945E+00,  .92741E+00,  .93562E+00,  .94409E+00, &
     &  .95282E+00,  .96179E+00,  .97100E+00,  .98044E+00,  .99011E+00, &
     &  .10000E+01,  .10000E+01]
     
   real ,PARAMETER :: CGAUSS(1:102) = [&
     &  1.00000E+00, 1.01822E+00, 1.03376E+00, 1.04777E+00, 1.06057E+00,&
     &  1.07231E+00, 1.08310E+00, 1.09300E+00, 1.10204E+00, 1.11025E+00,&
     &  1.11766E+00, 1.12429E+00, 1.13014E+00, 1.13523E+00, 1.13955E+00,&
     &  1.14313E+00, 1.14595E+00, 1.14803E+00, 1.14936E+00, 1.14994E+00,&
     &  1.14978E+00, 1.14888E+00, 1.14723E+00, 1.14484E+00, 1.14170E+00,&
     &  1.13782E+00, 1.13319E+00, 1.12782E+00, 1.12171E+00, 1.11485E+00,&
     &  1.10726E+00, 1.09893E+00, 1.08986E+00, 1.08006E+00, 1.06953E+00,&
     &  1.05828E+00, 1.04631E+00, 1.03363E+00, 1.02024E+00, 1.00617E+00,&
     &  9.91403E-01, 9.75966E-01, 9.59868E-01, 9.43123E-01, 9.25745E-01,&
     &  9.07752E-01, 8.89162E-01, 8.69994E-01, 8.50272E-01, 8.30017E-01,&
     ! DATA CG51 /                                                       &
     &  8.09256E-01, 7.88017E-01, 7.66327E-01, 7.44219E-01, 7.21726E-01,&
     &  6.98886E-01, 6.75729E-01, 6.52299E-01, 6.28637E-01, 6.04787E-01,&
     &  5.80794E-01, 5.56704E-01, 5.32565E-01, 5.08428E-01, 4.84343E-01,&
     &  4.60364E-01, 4.36543E-01, 4.12933E-01, 3.89589E-01, 3.66564E-01,&
     &  3.43913E-01, 3.21688E-01, 2.99940E-01, 2.78720E-01, 2.58077E-01,&
     &  2.38056E-01, 2.18701E-01, 2.00053E-01, 1.82148E-01, 1.65021E-01,&
     &  1.48701E-01, 1.33213E-01, 1.18579E-01, 1.04815E-01, 9.19338E-02,&
     &  7.99428E-02, 6.88453E-02, 5.86399E-02, 4.93211E-02, 4.08796E-02,&
     &  3.33018E-02, 2.65710E-02, 2.06669E-02, 1.55667E-02, 1.12449E-02,&
     &  7.67360E-03, 4.82345E-03, 2.66344E-03, 1.16151E-03, 2.84798E-04,&
     &  0.         , 0.         ]
     
   real ,PARAMETER :: CF1(1:102) =[ &
     &  0.         ,-2.56288E-03,-3.05202E-03,-2.50689E-03,-1.18504E-03,&
     &  7.84668E-04, 3.32528E-03, 6.38605E-03, 9.93124E-03, 1.39345E-02,&
     &  1.83758E-02, 2.32392E-02, 2.85120E-02, 3.41837E-02, 4.02454E-02,&
     &  4.66897E-02, 5.35099E-02, 6.07003E-02, 6.82556E-02, 7.61711E-02,&
     &  8.44422E-02, 9.30647E-02, 1.02034E-01, 1.11348E-01, 1.21000E-01,&
     &  1.30988E-01, 1.41307E-01, 1.51952E-01, 1.62921E-01, 1.74208E-01,&
     &  1.85808E-01, 1.97716E-01, 2.09927E-01, 2.22436E-01, 2.35236E-01,&
     &  2.48321E-01, 2.61684E-01, 2.75318E-01, 2.89215E-01, 3.03367E-01,&
     &  3.17764E-01, 3.32399E-01, 3.47260E-01, 3.62338E-01, 3.77620E-01,&
     &  3.93096E-01, 4.08752E-01, 4.24575E-01, 4.40550E-01, 4.56665E-01,&
     ! DATA CFA51 /                                                      &
     &  4.72901E-01, 4.89244E-01, 5.05675E-01, 5.22177E-01, 5.38731E-01,&
     &  5.55315E-01, 5.71913E-01, 5.88502E-01, 6.05059E-01, 6.21561E-01,&
     &  6.37986E-01, 6.54308E-01, 6.70504E-01, 6.86549E-01, 7.02417E-01,&
     &  7.18083E-01, 7.33520E-01, 7.48703E-01, 7.63607E-01, 7.78204E-01,&
     &  7.92472E-01, 8.06384E-01, 8.19918E-01, 8.33050E-01, 8.45759E-01,&
     &  8.58025E-01, 8.69828E-01, 8.81151E-01, 8.91979E-01, 9.02298E-01,&
     &  9.12097E-01, 9.21366E-01, 9.30098E-01, 9.38289E-01, 9.45935E-01,&
     &  9.53036E-01, 9.59594E-01, 9.65613E-01, 9.71101E-01, 9.76064E-01,&
     &  9.80513E-01, 9.84460E-01, 9.87919E-01, 9.90904E-01, 9.93432E-01,&
     &  9.95519E-01, 9.97184E-01, 9.98445E-01, 9.99322E-01, 9.99834E-01,&
     &  1.00000E+00, 1.00000E+00]
     
   real ,PARAMETER :: CF2(1:102) = [&
     &  0.         , 1.15907E-02, 2.32978E-02, 3.51022E-02, 4.69967E-02,&
     &  5.89773E-02, 7.10411E-02, 8.31858E-02, 9.54097E-02, 1.07711E-01,&
     &  1.20089E-01, 1.32541E-01, 1.45066E-01, 1.57663E-01, 1.70330E-01,&
     &  1.83065E-01, 1.95868E-01, 2.08737E-01, 2.21669E-01, 2.34664E-01,&
     &  2.47718E-01, 2.60830E-01, 2.73998E-01, 2.87219E-01, 3.00491E-01,&
     &  3.13812E-01, 3.27178E-01, 3.40587E-01, 3.54035E-01, 3.67520E-01,&
     &  3.81037E-01, 3.94583E-01, 4.08155E-01, 4.21747E-01, 4.35356E-01,&
     &  4.48978E-01, 4.62606E-01, 4.76237E-01, 4.89864E-01, 5.03482E-01,&
     &  5.17086E-01, 5.30669E-01, 5.44225E-01, 5.57746E-01, 5.71226E-01,&
     &  5.84657E-01, 5.98032E-01, 6.11342E-01, 6.24580E-01, 6.37736E-01,&
     ! DATA CFB51 /                                                      &
     &  6.50802E-01, 6.63769E-01, 6.76626E-01, 6.89365E-01, 7.01974E-01,&
     &  7.14444E-01, 7.26764E-01, 7.38924E-01, 7.50912E-01, 7.62717E-01,&
     &  7.74328E-01, 7.85735E-01, 7.96925E-01, 8.07888E-01, 8.18612E-01,&
     &  8.29087E-01, 8.39302E-01, 8.49246E-01, 8.58910E-01, 8.68284E-01,&
     &  8.77358E-01, 8.86125E-01, 8.94577E-01, 9.02706E-01, 9.10506E-01,&
     &  9.17972E-01, 9.25100E-01, 9.31885E-01, 9.38325E-01, 9.44419E-01,&
     &  9.50166E-01, 9.55568E-01, 9.60625E-01, 9.65340E-01, 9.69718E-01,&
     &  9.73763E-01, 9.77481E-01, 9.80878E-01, 9.83962E-01, 9.86741E-01,&
     &  9.89223E-01, 9.91419E-01, 9.93337E-01, 9.94989E-01, 9.96385E-01,&
     &  9.97536E-01, 9.98452E-01, 9.99146E-01, 9.99628E-01, 9.99909E-01,&
     &  1.00000E+00, 1.00000E+00]
     
   real ,PARAMETER :: CF3(1:102) = [&
     &  0.         , 9.88700E-03, 1.98515E-02, 2.99036E-02, 4.00474E-02,&
     &  5.02856E-02, 6.06200E-02, 7.10521E-02, 8.15830E-02, 9.22137E-02,&
     &  1.02945E-01, 1.13778E-01, 1.24712E-01, 1.35749E-01, 1.46889E-01,&
     &  1.58132E-01, 1.69478E-01, 1.80928E-01, 1.92480E-01, 2.04136E-01,&
     &  2.15894E-01, 2.27754E-01, 2.39716E-01, 2.51780E-01, 2.63943E-01,&
     &  2.76205E-01, 2.88564E-01, 3.01020E-01, 3.13571E-01, 3.26214E-01,&
     &  3.38948E-01, 3.51771E-01, 3.64679E-01, 3.77670E-01, 3.90741E-01,&
     &  4.03888E-01, 4.17108E-01, 4.30397E-01, 4.43750E-01, 4.57162E-01,&
     &  4.70628E-01, 4.84142E-01, 4.97700E-01, 5.11293E-01, 5.24915E-01,&
     &  5.38560E-01, 5.52218E-01, 5.65882E-01, 5.79542E-01, 5.93190E-01,&
     ! DATA CFC51 /                                                      &
     &  6.06816E-01, 6.20408E-01, 6.33957E-01, 6.47451E-01, 6.60877E-01,&
     &  6.74223E-01, 6.87477E-01, 7.00624E-01, 7.13651E-01, 7.26544E-01,&
     &  7.39288E-01, 7.51868E-01, 7.64268E-01, 7.76474E-01, 7.88470E-01,&
     &  8.00240E-01, 8.11768E-01, 8.23041E-01, 8.34042E-01, 8.44756E-01,&
     &  8.55171E-01, 8.65271E-01, 8.75044E-01, 8.84478E-01, 8.93562E-01,&
     &  9.02285E-01, 9.10639E-01, 9.18616E-01, 9.26210E-01, 9.33414E-01,&
     &  9.40227E-01, 9.46644E-01, 9.52666E-01, 9.58293E-01, 9.63528E-01,&
     &  9.68373E-01, 9.72833E-01, 9.76915E-01, 9.80625E-01, 9.83973E-01,&
     &  9.86967E-01, 9.89617E-01, 9.91935E-01, 9.93933E-01, 9.95622E-01,&
     &  9.97015E-01, 9.98125E-01, 9.98965E-01, 9.99549E-01, 9.99889E-01,&
     &  1.00000E+00, 1.00000E+00]
     
   real ,PARAMETER :: CER(1:102) = [&
     &  0.         ,-2.11394E-02,-4.08818E-02,-5.97585E-02,-7.79266E-02,&
     & -9.54663E-02,-1.12425E-01,-1.28834E-01,-1.44713E-01,-1.60076E-01,&
     & -1.74933E-01,-1.89289E-01,-2.03149E-01,-2.16515E-01,-2.29388E-01,&
     & -2.41768E-01,-2.53653E-01,-2.65043E-01,-2.75936E-01,-2.86328E-01,&
     & -2.96217E-01,-3.05601E-01,-3.14476E-01,-3.22839E-01,-3.30686E-01,&
     & -3.38015E-01,-3.44822E-01,-3.51105E-01,-3.56859E-01,-3.62083E-01,&
     & -3.66773E-01,-3.70928E-01,-3.74546E-01,-3.77625E-01,-3.80164E-01,&
     & -3.82161E-01,-3.83618E-01,-3.84534E-01,-3.84911E-01,-3.84749E-01,&
     & -3.84051E-01,-3.82821E-01,-3.81062E-01,-3.78778E-01,-3.75976E-01,&
     & -3.72663E-01,-3.68845E-01,-3.64532E-01,-3.59733E-01,-3.54461E-01,&
     ! DATA CER51 /                                                      &
     & -3.48726E-01,-3.42543E-01,-3.35927E-01,-3.28893E-01,-3.21461E-01,&
     & -3.13650E-01,-3.05477E-01,-2.96967E-01,-2.88142E-01,-2.79029E-01,&
     & -2.69652E-01,-2.60040E-01,-2.50221E-01,-2.40225E-01,-2.30084E-01,&
     & -2.19829E-01,-2.09493E-01,-1.99109E-01,-1.88712E-01,-1.78335E-01,&
     & -1.68014E-01,-1.57782E-01,-1.47673E-01,-1.37721E-01,-1.27957E-01,&
     & -1.18414E-01,-1.09120E-01,-1.00105E-01,-9.13939E-02,-8.30122E-02,&
     & -7.49818E-02,-6.73226E-02,-6.00518E-02,-5.31840E-02,-4.67313E-02,&
     & -4.07029E-02,-3.51053E-02,-2.99424E-02,-2.52153E-02,-2.09229E-02,&
     & -1.70614E-02,-1.36249E-02,-1.06056E-02,-7.99360E-03,-5.77750E-03,&
     & -3.94443E-03,-2.48028E-03,-1.36995E-03,-5.97540E-04,-1.46532E-04,&
     &  0.         , 0.         ]
            
   
CONTAINS !=====================Module Contains =========================

!-----------------------------------------------------------------------
!    MOLEC MAKES THE MOLECULAR IDENTIFICATIONS                         
!                                                                      
!    SCOR IS THE FACTOR BY WHICH THE LINE INTENSITY IS CHANGED DUE TO  
!       TEMPERATURE DEPENDENCE OF THE VIB AND ROT PARTITION SUMS       
!                                                                      
!    RHOSLF IS A QUANTITY ('PARTIAL DENSITY') FOR CORRECTING THE       
!       COLLISION WIDTH FOR SELF BROADENING                            
!                                                                      
!    ALFD1 CONTAINS THE DOPPLER WIDTHS AT 1 CM-1                       
!-----------------------------------------------------------------------
      SUBROUTINE MOLEC( IND,SCOR,RHOSLF,ALFD1, &
                        P,TEMP,WK,WTOT,molNames ) 
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: MaxISOTPL_smass, isotpl_mass, &
                                    BOLTZ, AVOGAD, CLIGHT, RADCN2, &
                                    TEMP0=>TEMP296, P0=>Press1013, &
                                    molNum
      USE Module_TIPS        ,ONLY: TIPS_2003

      IMPLICIT NONE !REAL*8           (V) 

      integer     ,intent(in)  :: IND
      real        ,intent(out) :: SCOR(:,:)  !(42,10)
      real        ,intent(out) :: RHOSLF(:)
      real        ,intent(out) :: ALFD1(:,:) !(42,10)
      !           
      real        ,intent(in)  :: P
      real        ,intent(in)  :: TEMP
      real        ,intent(in)  :: WK(:)      !(60)
      real        ,intent(in)  :: WTOT
      character(*),intent(in)  :: molNames(:)

      
      !--- Local variable
      !
      real  ,SAVE :: FAD
      INTEGER     :: ISO, M, imol
      character*8 :: HMOLID(60)
      REAL        :: FLN2, RHORAT, XKT, XKT0

      
      

      IF (IND.EQ.1) THEN 

         FLN2 = LOG(2.) 
         FAD = FLN2*2.*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) 
         XKT0 = TEMP0/RADCN2 

         RETURN 
      ELSE 

         RHORAT = (P/P0)*(TEMP0/TEMP) 
         XKT = TEMP/RADCN2 

!     call tips to get partition sum correction to intensities          
!yma         call tips_2003(nmol,iso_max,temp,scor) 
         call tips_2003( molNames,temp,scor) 

   
         DO 50 M = 1, size(molNames) 
            imol = molNum( molNames(m) )
            DO 40 ISO = 1, MaxISOTPL_smass(imol) 
               RHOSLF(imol) = RHORAT*WK(imol)/WTOT 
               ALFD1(imol,iso) = SQRT(FAD*TEMP/isotpl_mass(imol,iso)) 
   40       CONTINUE 
   50    CONTINUE 
      
      RETURN                                                            
      ENDIF 

  900 FORMAT (A6) 

      END SUBROUTINE

!-----------------------------------------------------------------------                                                                       
!     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE   
!                                                                       !                                                                       
!               LAST MODIFICATION:    12 AUGUST 1991                    
!                                                                       
!                  IMPLEMENTATION:    R.D. WORSHAM                      
!                                                                       
!             ALGORITHM REVISIONS:    S.A. CLOUGH                       
!                                     R.D. WORSHAM                      
!                                     J.L. MONCET                       
!                                                                       
!                                                                       
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.       
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139        
!                                                                       
!                                                                       
!               WORK SUPPORTED BY:    THE ARM PROGRAM                   
!                                     OFFICE OF ENERGY RESEARCH         
!                                     DEPARTMENT OF ENERGY              
!                                                                       
!                                                                       
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL           
!                                                                       
!                                             FASCOD3                   
!                                                                       
!                                                                       
!-----------------------------------------------------------------------                                                                       
      FUNCTION RADFN( VI,XKT ) 
!-----------------------------------------------------------------------                                                                       
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE !REAL*8           (V) 
      
      real                 :: RADFN
      real(r8) ,intent(in) :: VI
      real     ,intent(in) :: XKT
      
      REAL :: EXPVKT, XVI, XVIOKT
      
!      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                       

      XVI = VI 

      IF (XKT.GT.0.0) THEN 

         XVIOKT = XVI/XKT 

         IF (XVIOKT.LE.0.01) THEN 
            RADFN = 0.5*XVIOKT*XVI 

         ELSEIF (XVIOKT.LE.10.0) THEN 
            EXPVKT = EXP(-XVIOKT) 
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT) 

         ELSE 
            RADFN = XVI 
         ENDIF 

      ELSE 
         RADFN = XVI 
      ENDIF 

      RETURN 
      END FUNCTION
      
!-----------------------------------------------------------------------                                                                       
!     FUNCTION RADFNI CALCULATES THE RADIATION TERM FOR THE LINE SHAPE  
!     OVER INTERVAL VI TO VINEW                                         
!                                                                       
! INPUTS: 
!         VI: Starting Wavenumber in cm-1
!        DVI: Wavenumber step in cm-1
!        XKT: Tave/RADCN2 = kb*Tave/(h*c) 
!
! INPUT/OUTPUT:
!     RDLAST: If negative, this is the first call
!             If positive, used as input, should be the radiation term (in cm-1) at VI 
!             In either case, returned value is the radiation term at VINEW
!      VINEW: If negative, this is used as input to determine location of VINEW,
!                  but is adjusted to ensure that VINEW-VI is an integer number of DVI's 
!             If positive, VINEW is initally estimates as 1.003*VI, then adjusted 
!                  to ensure that VINEW-VI is an integer number of DVI's 
!             Adjusted value is returned
!
! OUTPUT:
!     RADFNI: Radiation term (in cm-1) at VI  
!       RDEL: Linear change for each DVI for a straight line connecting R(VI) and R(VINEW)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!               LAST MODIFICATION:    23 AUGUST 1991                    
!                                                                       
!                  IMPLEMENTATION:    R.D. WORSHAM                      
!                                                                       
!             ALGORITHM REVISIONS:    S.A. CLOUGH                       
!                                     R.D. WORSHAM                      
!                                     J.L. MONCET                       
!                                                                       
!                                                                       
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.       
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!               WORK SUPPORTED BY:    THE ARM PROGRAM                   
!                                     OFFICE OF ENERGY RESEARCH         
!                                     DEPARTMENT OF ENERGY              
!                                                                       
!                                                                       
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL           
!                                                                       
!                                             FASCOD3                   
!                                                                       
!-----------------------------------------------------------------------                                                                       
      FUNCTION RADFNI( VI,DVI,XKT,VINEW,RDEL,RDLAST ) 
!-----------------------------------------------------------------------                                                                       
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE !REAL*8           (V) 

      real                    :: RADFNI
      real(r8) ,intent(in)    :: VI
      real     ,intent(in)    :: DVI
      real     ,intent(in)    :: XKT
      real(r8) ,intent(inout) :: VINEW
      real     ,intent(out)   :: RDEL
      real     ,intent(inout) :: RDLAST

      real    ,PARAMETER :: FACT1=3.0E-03
      integer ,PARAMETER :: I_1=1 

      INTEGER :: INTVLS     
      REAL    :: CVIKT,    EXPVKT, EXPVNEWKT
      REAL    :: RDNEXT,   XMINUS, XMINUSNEW, XPLUS
      REAL    :: XPLUSNEW, XVI,    XVINEW,    XVINEWOKT
      REAL    :: XVIOKT    

      
      
!     RADFNI IS COMPUTED AT VI AND AND CALCULATES THE                   
!     WAVENUMBER VALUE (VINEW) FOR NEXT RADFNI CALC.                    

!     IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                        

      XVI = VI 

!     IF FIRST CALL, INITIALIZE RDLAST  (Compare to RADFN)                               

      IF (RDLAST.LT.0.) THEN 
         IF (XKT.GT.0.0) THEN 
            XVIOKT = XVI/XKT 

            IF (XVIOKT.LE.0.01) THEN 
               RDLAST = 0.5*XVIOKT*XVI 

            ELSEIF (XVIOKT.LE.10.0) THEN 
               EXPVKT = EXP(-XVIOKT) 
               RDLAST = XVI*(1.-EXPVKT)/(1.+EXPVKT) 

            ELSE 
               RDLAST = XVI 
            ENDIF 
         ELSE 
            RDLAST = XVI 
         ENDIF 
      ENDIF 

!     SET RADFNI EQUAL TO RADIATION FUNCTION AT VI                      

!     RDLAST IS RADFNI(VI) FOR EACH SUBSEQUENT CALL                     

      RADFNI = RDLAST 

      INTVLS = 1 
      IF (XKT.GT.0.0) THEN 

         XVIOKT = XVI/XKT 

         IF (XVIOKT.LE.0.01) THEN 
            IF (VINEW.GE.0.0) THEN 
               VINEW = VI+FACT1*0.5*XVI 
               INTVLS = (VINEW-VI)/DVI 
               INTVLS = MAX(INTVLS,I_1) 
               VINEW = VI+DVI* REAL(INTVLS) 
            ELSE 
               VINEW = ABS(VINEW) 
               INTVLS = (VINEW-VI)/DVI 
               INTVLS = MAX(INTVLS,I_1) 
            ENDIF 
            XVINEW = VINEW 

!            RDNEXT = 0.5*XVIOKT*XVINEW 
!MJA 20150821, fixing bug in RDLAST reported by Yingtao Ma
            XVINEWOKT = XVINEW/XKT 
            RDNEXT = 0.5*XVINEWOKT*XVINEW 

         ELSEIF (XVIOKT.LE.10.0) THEN 
            EXPVKT = EXP(-XVIOKT) 
            XMINUS = 1.-EXPVKT 
            XPLUS = 1.+EXPVKT 
            IF (VINEW.GE.0.0) THEN 
               CVIKT = XVIOKT*EXPVKT 
               VINEW = VI+FACT1*XVI/(1.+(CVIKT/XMINUS+CVIKT/XPLUS)) 
               INTVLS = (VINEW-VI)/DVI 
               INTVLS = MAX(INTVLS,I_1) 
               VINEW = VI+DVI* REAL(INTVLS) 
            ELSE 
               VINEW = ABS(VINEW) 
               INTVLS = (VINEW-VI)/DVI 
               INTVLS = MAX(INTVLS,I_1) 
            ENDIF 
            XVINEW = VINEW 

!            RDNEXT = XVINEW*XMINUS/XPLUS 
!MJA 20150821, fixing bug in RDLAST reported by Yingtao Ma
            XVINEWOKT = XVINEW/XKT 
            EXPVNEWKT = EXP(-XVINEWOKT) 
            XMINUSNEW = 1.-EXPVNEWKT 
            XPLUSNEW = 1.+EXPVNEWKT                                                             
            RDNEXT = XVINEW*XMINUSNEW/XPLUSNEW 

         ELSE  !XVIOKT > 10.0, assume radiation term is just wavenumber
            IF (VINEW.GE.0.0) THEN 
               VINEW = VI+(FACT1*XVI) 
               INTVLS = (VINEW-VI)/DVI 
               INTVLS = MAX(INTVLS,I_1) 
               VINEW = VI+DVI* REAL(INTVLS) 
            ELSE 
               VINEW = ABS(VINEW) 
               INTVLS = (VINEW-VI)/DVI 
               INTVLS = MAX(INTVLS,I_1) 
            ENDIF 
            XVINEW = VINEW 

            RDNEXT = XVINEW 
         ENDIF 
      ELSE !XKT < 0.0, assume radiation term is just wavenumber
         IF (VINEW.GE.0.0) THEN 
            VINEW = VI+(FACT1*XVI) 
            INTVLS = (VINEW-VI)/DVI 
            INTVLS = MAX(INTVLS,I_1) 
            VINEW = VI+DVI* REAL(INTVLS) 
         ELSE 
            VINEW = ABS(VINEW) 
            INTVLS = (VINEW-VI)/DVI 
            INTVLS = MAX(INTVLS,I_1) 
         ENDIF 
         XVINEW = VINEW 

!         RDNEXT = XVI 
! MJA, 20150821 - fixed bug that used beginning frequency to calculate 
! Radiation term at VINEW
         RDNEXT = XVINEW 
         
      ENDIF 

      RDEL = (RDNEXT-RADFNI)/ REAL(INTVLS) 

      RDLAST = RDNEXT 

      RETURN 
      END FUNCTION
      
      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine line_exception( ind,ipr,h_sub,mol,nmol,iso,iso_max ) 
!-----------------------------------------------------------------------
      IMPLICIT NONE

      integer     ,intent(in) :: ind
      integer     ,intent(in) :: ipr
      character*8 ,intent(in) :: h_sub 
      integer     ,intent(in) :: mol
      integer     ,intent(in) :: nmol
      integer     ,intent(in) :: iso
      integer     ,intent(in) :: iso_max(:)
            
      integer ,SAVE :: mol_max_pr_1=-99
      integer ,SAVE :: iso_max_pr_1=-99

      
      if ((ind.eq.1 .and. mol_max_pr_1.lt.0) .or.                       &
     &    (ind.eq.2 .and. iso_max_pr_1.lt.0)) then                      
         write (*,*) 
         write (*,*) 'Line file exception encountered in', h_sub 
         write (*,*) 'This message only written for first exception',   &
     &               ' for molecule and isotope cases'                  
         write (*,*) 'Other exceptions may exist' 

         write (ipr,*) '****************************************' 
         write (ipr,*) 'Line file exception encountered' 
         write (ipr,*) 'This message only written for first exception' 
         write (ipr,*) 'Other exceptions may exist' 
      endif 

      if (ind .eq. 1) then 
         if (mol_max_pr_1 .lt. 0) then 
            mol_max_pr_1 = 11 
            write (*,*) 
            write (*,*)   ' tape3: molecule number ', mol,             &
     &             ' greater than ', nmol,' encountered and skipped'    
            write (ipr,*) ' tape3: molecule number ', mol,             &
     &             ' greater than ', nmol,' encountered and skipped'    
            write (*,*) 
         endif 
         go to 25 

      else if (ind .eq. 2) then 
      
         if (iso_max_pr_1 .lt. 0) then 
            iso_max_pr_1 = 11 
            write (*,*) 
            write (*,*)   ' tape3: molecule number ', mol 
            write (ipr,*) ' tape3: molecule number ', mol 

            write (*,*)   ' tape3: isotope number ', iso,            &
     &                    ' greater than ', iso_max(mol),            &
     &                    ' encountered and skipped'                 
            write (ipr,*) ' tape3: isotope number ', iso,            &
     &                    ' greater than ', iso_max(mol),            &
     &                    ' encountered and skipped'                 
            write (*,*) 
         endif 
         go to 25 
      endif 

   25 continue 

      return 
      END  SUBROUTINE

      
END MODULE 