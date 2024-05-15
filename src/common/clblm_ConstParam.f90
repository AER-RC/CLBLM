!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

!--------------------------------------------------------------------
! This MODULE contains contents of "BLOCK DATA ATMCON". "BLOCK DATA ATMCON"
! will be removed in the future. Now both are used in the program.
! Just remember don't allow them both appear in a same subroutine/functions.
!--------------------------------------------------------------------
MODULE Module_ConstParam
   USE planet_consts, ONLY: AIRMWT, XMASS_DRY, GRAV_CONST

   IMPLICIT NONE

   integer :: I,J,K


   !integer, PARAMETER :: kind_dp=kind(1.0D0)
   integer, PARAMETER :: kind_r8=selected_real_kind(14)

   real, PARAMETER :: PI = 3.1415926535898 !2.0*ASIN(1.0) !DACOS(-1.D0) !3.1415926535898
   real, PARAMETER :: DEG = 180.0/PI

   ! Parameters for flux calculations
   real, PARAMETER :: EPS = 1.E-4
   !     Here are the weights for the first-order Gaussian quadrature:
! Quadrature weights for angles 1-3, from Radsum (Clough et al. 1992 Table 2)
! 1 angle
   real, PARAMETER :: GWGO1 = 0.50
! 2 angles
   real, PARAMETER :: GWGD1= 0.31804138
   real, PARAMETER :: GWGD2 = 0.18195862
! 3 angles
   real, PARAMETER :: GWGT1 =0.20093191
   real, PARAMETER :: GWGT2 =0.22924111
   real, PARAMETER :: GWGT3 =0.06982698

!HEATFC is the factor one must multiply DELTA-FLUX/DELTA-PRESSURE,
!     with flux in W/M-2 and pressure in Mb, to get the heating rate in
!     units of Degrees/day.  It is equal to
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(3600)(1e-5)/(1.004)
   real, PARAMETER :: HEATFC = 8.4391

   integer ,PARAMETER :: FILLINT = -999
   real    ,PARAMETER :: FILLREAL = -999.9
   real    ,PARAMETER :: ONEPL = 1.01
   real    ,PARAMETER :: ONEMI = 0.999
   real    ,PARAMETER :: ARGMIN = 34.
   real    ,PARAMETER :: EXPMIN = EXP(-ARGMIN)
   !real    ,PARAMETER :: EXPMIN_r8 = tiny(0_kind_r8)
   !real    ,PARAMETER :: EXPMIN_r4 = tiny(0_kind_r4)



   !Constants from NIST May 2010; units are generally cgs
   real, PARAMETER  :: PLANCK = 6.62606876E-27
   real, PARAMETER  :: BOLTZ = 1.3806503E-16
   real, PARAMETER  :: CLIGHT = 2.99792458E+10
   real, PARAMETER  :: AVOGAD = 6.02214199E+23
   real, PARAMETER  :: ALOSMT = 2.6867775E+19  !=2006 CODATA; 2010 CODATA number is 2.6867805E19
   real, PARAMETER  :: GASCON = 8.314472E+07

   !The first and second radiation constants are taken from NIST.
   !They were previously obtained from the relations:
   !RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07
   !RADCN2 = PLANCK*CLIGHT/BOLTZ
   real, PARAMETER  :: RADCN1 = 1.191042722E-12
   real, PARAMETER  :: RADCN2 = 1.4387752

   ! AIRMS1 IS ONE AIRMASS OR THE TOTAL AMOUNT FOR A VERTICAL PATH
   ! FROM GROUND TO SPACE
   real, PARAMETER :: AIRMS1 = 2.153E25

   ! GCAIR IS THE GAS CONSTANT FOR RHO IN MOL CM(-3), P IN MB, AND T IN K
   real, PARAMETER :: GCAIR = 1.0E-3*GASCON/AVOGAD

   ! ADCON IS THE CONSTANT FOR THE DOPPLER HALFWIDTH
   real, PARAMETER :: ADCON = SQRT(2.0* LOG(2.0)*GASCON/CLIGHT**2)


   ! MXFSC  IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO FASE01
   ! MXLAY  IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
   ! MXZMD  IS THE MAXIMUM NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE STORED IN ZMDL (INPUT)
   ! MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH OBTAINED BY MERGING ZMDL AND ZOUT
   ! MXMOL  IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
   integer, parameter :: MXMOL=47, MXTRAC=41, MX_XS=50, MAX_ISO=20, MXISOTPL=10 !MXSPC=5,
   integer, parameter :: MXZMD=6000, MXFSC=600, MXLAY=MXFSC+3
   integer, parameter :: MXPDIM=MXLAY+MXZMD, IM2=MXPDIM-2
   integer, parameter :: NFPTS=2001, NFMX=1.5*NFPTS !yma NFMX=1.3*NFPTS
   integer, parameter :: NMAXCO=4040, NUMZ = 101 !for surface emissivity and reflectivity usage.
   integer, parameter :: IPTS=5050, IPTS2=6000
   integer, parameter :: N_ABSRB=5050!, ~nzeta=101
   integer, parameter :: NT=119, Nmax=MXFSC !Nmax=600
   !integer, parameter :: ~NN_TBL=10000, ~NDIM=2410, ~ND2=5000
   integer, parameter :: MAXSTATE=26
   integer, parameter :: NFLTPT=3001


   ! KMXNOM IS THE DEFAULT number of molecules
   integer, PARAMETER :: KMXNOM = 7

   ! IFXTYP is the flag for fixing the value of ITYL
   !integer, PARAMETER :: IFXTYP = 0

   ! DELTAS IS THE NOMINAL SLANT PATH INCREMENT IN KM.
   real, PARAMETER :: DELTAS = 5.0
   real, PARAMETER :: Press1013 = 1013.25
   real, PARAMETER :: Temp273 = 273.15
   real, PARAMETER :: Temp296 = 296.0

   ! ALZERO IS THE MEAN LORENTZ HALFWIDTH AT PZERO AND 296.0 K.
   ! AVMWT IS THE MEAN MOLECULAR WEIGHT USED TO AUTOMATICALLY
   ! GENERATE THE LBLRTM BOUNDARIES IN AUTLAY
   real, PARAMETER :: ALZERO = 0.04
   real, PARAMETER :: AVMWT  = 36.0

   real, PARAMETER :: secant_diffuse = 1.67






   !**************************************
   !      BLOCK DATA BDMol
   !**************************************
   character(6) ,PARAMETER :: clblmMolNames(1:MXMOL) =[ &   !removed the 0 element ('ALL') from clblmMolNames
        '   H2O','   CO2','    O3','   N2O','    CO','   CH4','    O2',   &
        '    NO','   SO2','   NO2','   NH3','  HNO3','    OH','    HF',   &
        '   HCL','   HBR','    HI','   CLO','   OCS','  H2CO','  HOCL',   &
        '    N2','   HCN',' CH3CL','  H2O2','  C2H2','  C2H6','   PH3',   &
        '  COF2','   SF6','   H2S',' HCOOH','   HO2','     O','CLONO2',   &
        '   NO+','  HOBR','  C2H4',' CH3OH',' CH3BR',' CH3CN','  CF4 ',   &
        ' C4H2 ',' HC3N ','   H2 ','   CS ','  SO3 ']

   ! MOLECULAR WEIGHTS
   ! approx;
   ! MJA, 10/04/2011 fixed ethylene molecular weight (Molecule #38)
   real, PARAMETER :: AMWT(MXMOL) = [ &
          18.015 ,  44.010 , 47.998 , 44.01 , 28.011 ,  16.043 , 31.999 ,&
          30.01  ,  64.06  , 46.01  , 17.03 , 63.01  ,  17.00  , 20.01  ,&
          36.46  ,  80.92  ,127.91  , 51.45 , 60.08  ,  30.03  , 52.46  ,&
          28.014 ,  27.03  , 50.49  , 34.01 , 26.03  ,  30.07  , 34.00  ,&
          66.01  , 146.05  , 34.08  , 46.03 , 33.00  ,  15.99  , 98.    ,&
          30.00  ,  97.    , 28.05  , 32.04 , 94.94  ,  41.05  , 88.0043,&
          50.06  ,  51.05  ,  2.016 , 44.08 , 80.066 ]

   character(6) ,PARAMETER :: CMOL(1:MXMOL) = clblmMolNames(1:MXMOL)
   character(8), PARAMETER :: HMOLC(1:MXMOL) = clblmMolNames(1:MXMOL)


   !  ****************************************
   !   BLOCK DATA BXSECT
   !  ****************************************
   ! ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES
   ! (NOTE: ALL NAMES ARE LEFT-JUSTIFIED)
   character(10) ,PROTECTED :: ALIAS(4,mx_xs)

   DATA (ALIAS(1,I),I=1,mx_xs)/                                      &
       'CLONO2    ', 'HNO4      ', 'CHCL2F    ', 'CCL4      ',       &
       'CCL3F     ', 'CCL2F2    ', 'C2CL2F4   ', 'C2CL3F3   ',       &
       'N2O5      ', 'HNO3      ', 'CF4       ', ' ZZZZZZZZ ',       &
       'CCLF3     ', 'C2CLF5    ', 'NO2       ', 'PAN       ',       &
       'ACET      ', 'CH3CN     ', 'CHF2CF3   ', 'CFH2CF3   ',       &
       'CF3CH3    ', 'CH3CHF2   ', 'CH2F2     ', 'CCl2FCH3  ',       &
       'CH3CClF2  ', 'CHClF2    ', 'CHCl2CF3  ', 'CHCl2C2F5 ',       &
       'C3HCl2F5  ', 'C3HCl2F5  ', 'SO2       ', 'ISOP      ',       &
       'CHF3      ', 'BRO       ', 'HCHO      ', 'CF3CH2CF3 ',       &
       'CHF2CH2CF3', 'NF3       ', 'C2F6      ', 'SF6       ',       &
       'FURAN     ', 'ACETICACI ', 'GLYCOLALD ', 'PROPENE   ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ' /
   DATA (ALIAS(2,I),I=1,mx_xs)/                                      &
       'CLNO3     ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',       &
       'CFCL3     ', 'CF2CL2    ', 'C2F4CL2   ', 'C2F3CL3   ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       'CH3COCH3  ', ' ZZZZZZZZ ', 'HFC-125   ', 'HFC-134a  ',       &
       'HFC-143a  ', 'HFC-152a  ', 'HFC-32    ', 'HCFC-141b ',       &
       'HCFC-142b ', 'HCFC-22   ', 'HCFC-123  ', 'HCFC-124  ',       &
       'HCFC-225ca', 'HCFC-225cb', ' ZZZZZZZZ ', 'C5H8      ',       &
       'HFC-23    ', 'BRO       ', ' ZZZZZZZZ ', 'HFC-236fa ',       &
       'HFC-245fa ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ '/
   DATA (ALIAS(3,I),I=1,mx_xs)/                                      &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',       &
       'CFC11     ', 'CFC12     ', 'CFC114    ', 'CFC113    ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC14     ', 'CFC22     ',       &
       'CFC13     ', 'CFC115    ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       'ACETONE   ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CCL2FCH3  ',       &
       'CH3CCLF2  ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', 'BRO       ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ '/
   DATA (ALIAS(4,I),I=1,mx_xs)/                                      &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F21       ', ' ZZZZZZZZ ',       &
       'F11       ', 'F12       ', 'F114      ', 'F113      ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F14       ', 'F22       ',       &
       'F13       ', 'F115      ',  ' ZZZZZZZZ ',' ZZZZZZZZ ',       &
       'CH3C(O)CH3', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', 'BRO       ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ',       &
       ' ZZZZZZZZ ', ' ZZZZZZZZ '/

   !XSMASS IS MASS OF EACH CROSS-SECTION
   real ,PARAMETER :: XSMASS(MX_XS) = [&
           97.46     ,   79.01     ,  102.92     ,  153.82     ,       &
          137.37     ,  120.91     ,  170.92     ,  187.38     ,       &
          108.01     ,   63.01     ,   88.00     ,   86.47     ,       &
          104.46     ,  154.47     ,   45.99     ,  121.05     ,       &
           58.08     ,   41.05     ,  120.02     ,  102.03     ,       &
           84.04     ,   66.05     ,   52.02     ,  116.95     ,       &
          100.50     ,   86.47     ,  152.93     ,  136.48     ,       &
          202.94     ,  202.94     ,   64.06     ,   68.12     ,       &
           70.01     ,   95.903    ,   30.026    ,  152.039    ,       &
          134.05     ,   71.00     ,  138.01     ,  146.06     ,       &
           68.075    ,   60.052    ,   60.052    ,   42.081    ,       &
           0.        ,   0.        ,   0.        ,   0.        ,       &
           0.       ,    0.]

!  Furan, ACETICACID, GLYCOLALDE, PROPENE  weights are from
!  https://en.wikipedia.org/wiki/Furan
!  https://en.wikipedia.org/wiki/Acetic_acid
!  https://en.wikipedia.org/wiki/Glycolaldehyde
!  https://en.wikipedia.org/wiki/Propene


   !  ****************************************
   !      BLOCK DATA Isotop
   !  ****************************************
   !--- The number of isotopes for a particular molecule:
   integer ,PARAMETER :: MaxISOTPL_smass(MXMOL) = [ &
     &   6,  10,    9,    5,    6,   4,     3,                          &
     &   3,   2,    1,    2,    2,   3,     2,    4,    4,    2,        &
     &   2,   5,    3,    2,    2,   3,     2,    1,    3,    2,    1,  &
     &   2,   1,    3,    1,    1,   1,     2,    1,    2,    2,    1,  &
     &   2,   1,    1,    1,    1,   2,     4,    1 ]
   !   H2O,  CO2,   O3,  N2O,   CO, CH4,    O2,
   !    NO,  SO2,  NO2,  NH3, HNO3,  OH,    HF,  HCl,  HBr,   HI,
   !   ClO,  OCS, H2CO, HOCl,   N2, HCN, CH3Cl, H2O2, C2H2, C2H6,  PH3,
   !  COF2,  SF6,  H2S,HCOOH,  HO2,   O,ClONO2,  NO+, HOBr, C2H4,CH3OH,
   ! CH3Br,CH3CN,  CF4, C4H2, HC3N,  H2,    CS, SO3/


   !--- MOLECULAR MASSES FOR EACH ISOTOPE
   real ,protected :: isotpl_mass(MXMOL,10)
   !
!  H2O:   161,   181,   171,   162,   182,   172
      data (isotpl_mass(1,i),i=1,6)                                           &
     & /  18.01, 20.01, 19.01, 19.01, 21.02, 20.02/
!  CO2:   626,   636,   628,   627,   638,   637,   828,   827,    727,
!         838
      data (isotpl_mass(2,i),i=1,10)                                          &
     & /  43.99, 44.99, 45.99, 44.99, 47.00, 46.00, 48.00, 47.00, 46.00,&
     &    49.00/
!   O3:   666,   668,   686    667    676    886    868    678    768
      data (isotpl_mass(3,i),i=1,9)                                           &
     & /  47.98, 49.99, 49.99, 48.99, 48.99, 51.99, 51.99, 50.99, 50.99/
!  N2O:   446,   456,   546,   448,   447
      data (isotpl_mass(4,i),i=1,5)                                           &
     & /  44.00, 45.00, 45.00, 46.00, 45.00/
!   CO:   26,    36,    28,    27,    38     37
      data (isotpl_mass(5,i),i=1,6)                                           &
     & /  27.99, 28.99, 29.99, 29.00, 31.00, 30.00/
!  CH4:   211,   311,   212,    312
      data (isotpl_mass(6,i),i=1,4)                                           &
     & /  16.03, 17.03, 17.03, 18.03/
!   O2:    66,    68,    67
      data (isotpl_mass(7,i),i=1,3)                                           &
     & /  31.99, 33.99, 32.99/
!   NO:    46,    56,    48,
      data (isotpl_mass(8,i),i=1,3)                                           &
     & /  30.00, 31.00, 32.00/
!  SO2:   626,   646
      data (isotpl_mass(9,i),i=1,2)                                           &
     & /  63.96, 65.96/
!  NO2:   646
      data (isotpl_mass(10,i),i=1,1)                                          &
     & /  45.99/
!  NH3:  4111,  5111;
      data (isotpl_mass(11,i),i=1,2)                                          &
     & /  17.03, 18.02/
! HNO3:
      data (isotpl_mass(12,i),i=1,2)                                          &
     & /  62.99, 63.99/
!   OH:    61,    81,    62
      data (isotpl_mass(13,i),i=1,3)                                          &
     & /   17.00, 19.01, 18.01/
!   HF:    19       29
      data (isotpl_mass(14,i),i=1,2)                                          &
     & /   20.01, 21.01/
!  HCL:   15,    17;       25,    27,
      data (isotpl_mass(15,i),i=1,4)                                          &
     & /  35.98, 37.97, 36.98, 38.97/
!  HBr:   19,    11;     29,     21
      data (isotpl_mass(16,i),i=1,4)                                          &
     & /  79.92, 81.92,80.92, 82.92 /
!   HI:   17, 27
      data (isotpl_mass(17,i),i=1,2)                                          &
     & /  127.91, 128.91/
!  ClO:   56,    76;
      data (isotpl_mass(18,i),i=1,2)                                          &
     & /  50.96, 52.96/
!  OCS:   622,   624,   632,   623,   822
      data (isotpl_mass(19,i),i=1,5)                                          &
     & /  59.97, 61.96, 60.97, 60.97, 61.97/
! H2CO:  126,   136,   128;
      data (isotpl_mass(20,i),i=1,3)                                          &
     & /  30.01, 31.01, 32.01/
! HOCl:  165,   167
      data (isotpl_mass(21,i),i=1,2)                                          &
     & /  51.97, 53.97/
!   N2:   44,      45;
      data (isotpl_mass(22,i),i=1,2)                                          &
     & /  28.01, 29.01/
!  HCN:   124,   134,   125,
      data (isotpl_mass(23,i),i=1,3)                                          &
     & /  27.01, 28.01, 28.01/
! CH3CL:  215,   217;
      data (isotpl_mass(24,i),i=1,2)                                          &
     & /  49.99, 51.99/
!  H2O2:  1661;
      data (isotpl_mass(25,i),i=1,1)                                          &
     & /  34.01/
!  C2H2: 1221,  1231 , 1222
      data (isotpl_mass(26,i),i=1,3)                                          &
     & /  26.01, 27.02, 27.01/
!  C2H6: 1221, , 1231;
      data (isotpl_mass(27,i),i=1,2)                                          &
     & /  30.05, 31.05/
!   PH3:   1111;
      data (isotpl_mass(28,i),i=1,1)                                          &
     & /  34.00/
!  COF2:  269, 369;
      data (isotpl_mass(29,i),i=1,2)                                          &
     & /  65.99, 66.99/
!   SF6:   29
      data (isotpl_mass(30,i),i=1,1)                                          &
     & /  145.96/
!   H2S:  121    141    131;
      data (isotpl_mass(31,i),i=1,3)                                          &
     & /  33.99, 35.98, 34.99/
! HCOOH: 126;
      data (isotpl_mass(32,i),i=1,1)                                          &
     & /  46.01/
!   HO2:   166
      data (isotpl_mass(33,i),i=1,1)                                          &
     & / 33.00/
!     O:   6
      data (isotpl_mass(34,i),i=1,1)                                          &
     & /  15.99/
! ClONO2: 5646   7646;
      data (isotpl_mass(35,i),i=1,2)                                          &
     & /  96.96, 98.95/
!   NO+:   46
      data (isotpl_mass(36,i),i=1,1)                                          &
     & /  30.00 /
!  HOBr:  169,    161
      data (isotpl_mass(37,i),i=1,2)                                          &
     & /  95.92,  97.92/
!   C2H4: 221,   231;
      data (isotpl_mass(38,i),i=1,2)                                          &
     & /  28.05, 29.05/
!   CH3OH: 2161
      data (isotpl_mass(39,i),i=1,1)                                          &
     & /  32.04/
!   CH3Br: 219, 211
      data (isotpl_mass(40,i),i=1,2)                                          &
     & /  93.94, 95.94/
!   CH3CN: 2124
      data (isotpl_mass(41,i),i=1,1)                                          &
     & /  41.05/
!   CF4: 29
      data (isotpl_mass(42,i),i=1,1)                                          &
     & /  88.0043/
!   C4H2: 2211
      data (isotpl_mass(43,i),i=1,1)                                          &
     & /  50.06/
!   HC3N: 1224
      data (isotpl_mass(44,i),i=1,1)                                          &
     & /  51.05/
!   H2: 11, 12
      data (isotpl_mass(45,i),i=1,2)                                          &
     & /  2.016, 3.022/
!   CS: 22, 24, 32, 23
      data (isotpl_mass(46,i),i=1,4)                                          &
     & /  44.08, 46.08, 45.08, 45.08/
!   SO3: 26
      data (isotpl_mass(47,i),i=1,1)                                          &
     & /  80.066/


!------------- Below are from LBLRTM v12.7 lblparams.f90 ---------------
!
! HITRAN 2012 isotopologue information
      integer ,protected :: MaxISOTPL_abd(mxmol)       ! Number of isotopologues for each species
      real    ,protected :: isotpl_abd(mxmol,mxisotpl)  ! Isotopologue abundance by species
      integer ,protected :: isotpl_code(mxmol,mxisotpl) ! Isotopologue code by species

! HITRAN 2012 isotopologue abundance values are based on:
! De Bievre, P., N.E. Holden, and I.L. Barnes, Isotopic Abundances and
! Atomic Weights of the Elements, J. Phys. Chem. Ref. Data, 13,
! 809-891, 1984.

! Molecule names
!      DATA CMOL   /                                                     &
!     &     '  H2O ','  CO2 ','   O3 ','  N2O ','   CO ','  CH4 ',       &
!     &     '   O2 ','   NO ','  SO2 ','  NO2 ','  NH3 ',' HNO3 ',       &
!     &     '   OH ','   HF ','  HCL ','  HBR ','   HI ','  CLO ',       &
!     &     '  OCS ',' H2CO ',' HOCL ','   N2 ','  HCN ','CH3CL ',       &
!     &     ' H2O2 ',' C2H2 ',' C2H6 ','  PH3 ',' COF2 ','  SF6 ',       &
!     &     '  H2S ','HCOOH ','  HO2 ','    O ','ClONO2','  NO+ ',       &
!     &     ' HOBr ',' C2H4 ','CH3OH ',' CH3Br',' CH3CN','  CF4 ',       &
!     &     ' C4H2 ',' HC3N ','   H2 ','   CS ','  SO3 '/

! Number of isotopologues for each species
      data MaxISOTPL_abd(:) /                                              &
     &           6,      10,       5,       5,       6,       4,        &
     &           3,       3,       2,       1,       2,       2,        &
     &           3,       2,       4,       4,       2,       2,        &
     &           5,       3,       2,       2,       3,       2,        &
     &           1,       3,       2,       1,       2,       1,        &
     &           3,       1,       1,       1,       2,       1,        &
     &           2,       2,       1,       2,       1,       1,        &
     &           1,       1,       2,       4,       1/

! HITRAN 2012/AFGL isotopologue codes for each species
! H2O
      data isotpl_code(1,:) /                                           &
     &         161,       181,      171,      162,      182,            &
     &         172,         0,        0,        0,        0/
! CO2
      data isotpl_code(2,:) /                                           &
     &         626,       636,      628,      627,      638,            &
     &         637,       828,      827,      727,      838/
! O3
      data isotpl_code(3,:) /                                           &
     &         666,       668,      686,      667,      676,            &
     &           0,         0,        0,        0,        0/
! N2O
      data isotpl_code(4,:) /                                           &
     &         446,       456,      546,      448,      447,            &
     &           0,         0,        0,        0,        0/
! CO
      data isotpl_code(5,:) /                                           &
     &          26,        36,       28,       27,       38,            &
     &          37,         0,        0,        0,        0/
! CH4
      data isotpl_code(6,:) /                                           &
     &         211,       311,      212,      312,        0,            &
     &           0,         0,        0,        0,        0/
! O2
      data isotpl_code(7,:) /                                           &
     &          66,        68,       67,        0,        0,            &
     &           0,         0,        0,        0,        0/
! NO
      data isotpl_code(8,:) /                                           &
     &          46,        56,       48,        0,        0,            &
     &           0,         0,        0,        0,        0/
! SO2
      data isotpl_code(9,:) /                                           &
     &         626,       646,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! NO2
      data isotpl_code(10,:) /                                          &
     &         646,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! NH3
      data isotpl_code(11,:) /                                          &
     &        4111,      5111,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HNO3
      data isotpl_code(12,:) /                                          &
     &         146,       156,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! OH
      data isotpl_code(13,:) /                                          &
     &          61,        81,       62,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HF
      data isotpl_code(14,:) /                                          &
     &          19,        29,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HCl
      data isotpl_code(15,:) /                                          &
     &          15,        17,       25,       27,        0,            &
     &           0,         0,        0,        0,        0/
! HBr
      data isotpl_code(16,:) /                                          &
     &          19,        11,       29,       21,        0,            &
     &           0,         0,        0,        0,        0/
! HI
      data isotpl_code(17,:) /                                          &
     &          17,        27,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! ClO
      data isotpl_code(18,:) /                                          &
     &          56,        76,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! OCS
      data isotpl_code(19,:) /                                          &
     &         622,       624,      632,      623,      822,            &
     &           0,         0,        0,        0,        0/
! H2CO
      data isotpl_code(20,:) /                                          &
     &         126,       136,      128,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HOCl
      data isotpl_code(21,:) /                                          &
     &         165,       167,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! N2
      data isotpl_code(22,:) /                                          &
     &          44,        45,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HCN
      data isotpl_code(23,:) /                                          &
     &         124,       134,      125,        0,        0,            &
     &           0,         0,        0,        0,        0/
! CH3Cl
      data isotpl_code(24,:) /                                          &
     &         215,       217,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! H2O2
      data isotpl_code(25,:) /                                          &
     &        1661,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! C2H2
      data isotpl_code(26,:) /                                          &
     &        1221,      1231,     1222,        0,        0,            &
     &           0,         0,        0,        0,        0/
! C2H6
      data isotpl_code(27,:) /                                          &
     &        1221,      1231,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! PH3
      data isotpl_code(28,:) /                                          &
     &        1111,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! COF2
      data isotpl_code(29,:) /                                          &
     &         269,       369,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! SF6
      data isotpl_code(30,:) /                                          &
     &          29,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! H2S
      data isotpl_code(31,:) /                                          &
     &         121,       141,      131,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HCOOH
      data isotpl_code(32,:) /                                          &
     &         126,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HO2
      data isotpl_code(33,:) /                                          &
     &         166,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! O
      data isotpl_code(34,:) /                                          &
     &           6,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! ClONO2
      data isotpl_code(35,:) /                                          &
     &        5646,      7646,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! NO+
      data isotpl_code(36,:) /                                          &
     &          46,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HOBr
      data isotpl_code(37,:) /                                          &
     &         169,       161,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! C2H4
      data isotpl_code(38,:) /                                          &
     &         221,       231,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! CH3OH
      data isotpl_code(39,:) /                                          &
     &        2161,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! CH3Br
      data isotpl_code(40,:) /                                          &
     &         219,       211,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! CH3CN
      data isotpl_code(41,:) /                                          &
     &        2124,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! CF4
      data isotpl_code(42,:) /                                          &
     &          29,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! C4H2
      data isotpl_code(43,:) /                                          &
     &        2211,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! HC3N
      data isotpl_code(44,:) /                                          &
     &        1224,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! H2
      data isotpl_code(45,:) /                                          &
     &          11,        12,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/
! CS
      data isotpl_code(46,:) /                                          &
     &          22,        24,       32,       23,        0,            &
     &           0,         0,        0,        0,        0/
! SO3
      data isotpl_code(47,:) /                                          &
     &          26,         0,        0,        0,        0,            &
     &           0,         0,        0,        0,        0/

! HITRAN 2012 isotopologue abundance for each species
! H2O
      data isotpl_abd(1,:) /                                            &
     &     9.973e-1, 1.999e-3, 3.719e-4, 3.107e-4, 6.230e-7,            &
     &     1.158e-7, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CO2
      data isotpl_abd(2,:) /                                            &
     &     9.842e-1, 1.106e-2, 3.947e-3, 7.340e-4, 4.434e-5,            &
     &     8.246e-6, 3.957e-6, 1.472e-6, 1.368e-7, 4.446e-8/
! O3
      data isotpl_abd(3,:) /                                            &
     &     9.929e-1, 3.982e-3, 1.991e-3, 7.405e-4, 3.702e-4,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! N2O
      data isotpl_abd(4,:) /                                            &
     &     9.903e-1, 3.641e-3, 3.641e-3, 1.986e-3, 3.693e-4,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CO
      data isotpl_abd(5,:) /                                            &
     &     9.865e-1, 1.108e-2, 1.978e-3, 3.679e-4, 2.223e-5,            &
     &     4.133e-6, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH4
      data isotpl_abd(6,:) /                                            &
     &     9.883e-1, 1.110e-2, 6.158e-4, 6.918e-6, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! O2
      data isotpl_abd(7,:) /                                            &
     &     9.953e-1, 3.991e-3, 7.422e-4, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NO
      data isotpl_abd(8,:) /                                            &
     &     9.940e-1, 3.654e-3, 1.993e-3, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! SO2
      data isotpl_abd(9,:) /                                            &
     &     9.457e-1, 4.195e-2, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NO2
      data isotpl_abd(10,:) /                                           &
     &     9.916e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NH3
      data isotpl_abd(11,:) /                                           &
     &     9.959e-1, 3.661e-3, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HNO3
      data isotpl_abd(12,:) /                                           &
     &     9.891e-1, 3.636e-3, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! OH
      data isotpl_abd(13,:) /                                           &
     &     9.975e-1, 2.000e-3, 1.554e-4, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HF
      data isotpl_abd(14,:) /                                           &
     &     9.998e-1, 1.557e-4, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HCl
      data isotpl_abd(15,:) /                                           &
     &     7.576e-1, 2.423e-1, 1.180e-4, 3.774e-5, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HBr
      data isotpl_abd(16,:) /                                           &
     &     5.068e-1, 4.931e-1, 7.894e-5, 7.680e-5, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HI
      data isotpl_abd(17,:) /                                           &
     &     9.998e-1, 1.557e-4, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! ClO
      data isotpl_abd(18,:) /                                           &
     &     7.559e-1, 2.417e-1, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! OCS
      data isotpl_abd(19,:) /                                           &
     &     9.374e-1, 4.158e-2, 1.053e-2, 7.399e-3, 1.880e-3,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2CO
      data isotpl_abd(20,:) /                                           &
     &     9.862e-1, 1.108e-2, 1.978e-3, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HOCl
      data isotpl_abd(21,:) /                                           &
     &     7.558e-1, 2.417e-1, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! N2
      data isotpl_abd(22,:) /                                           &
     &     9.927e-1, 7.478e-3, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HCN
      data isotpl_abd(23,:) /                                           &
     &     9.851e-1, 1.107e-2, 3.622e-3, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3Cl
      data isotpl_abd(24,:) /                                           &
     &     7.489e-1, 2.395e-1, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2O2
      data isotpl_abd(25,:) /                                           &
     &     9.950e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C2H2
      data isotpl_abd(26,:) /                                           &
     &     9.776e-1, 2.197e-2, 3.046e-4, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C2H6
      data isotpl_abd(27,:) /                                           &
     &     9.770e-1, 2.195e-2, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! PH3
      data isotpl_abd(28,:) /                                           &
     &     9.995e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! COF2
      data isotpl_abd(29,:) /                                           &
     &     9.865e-1, 1.108e-2, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! SF6
      data isotpl_abd(30,:) /                                           &
     &     9.502e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2S
      data isotpl_abd(31,:) /                                           &
     &     9.499e-1, 4.214e-2, 7.498e-3, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HCOOH
      data isotpl_abd(32,:) /                                           &
     &     9.839e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HO2
      data isotpl_abd(33,:) /                                           &
     &     9.951e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! O
      data isotpl_abd(34,:) /                                           &
     &     9.976e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! ClONO2
      data isotpl_abd(35,:) /                                           &
     &     7.496e-1, 2.397e-1, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NO+
      data isotpl_abd(36,:) /                                           &
     &     9.940e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HOBr
      data isotpl_abd(37,:) /                                           &
     &     5.056e-1, 4.919e-1, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C2H4
      data isotpl_abd(38,:) /                                           &
     &     9.773e-1, 2.196e-2, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3OH
      data isotpl_abd(39,:) /                                           &
     &     9.859e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3Br
      data isotpl_abd(40,:) /                                           &
     &     5.010e-1, 4.874e-1, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3CN
      data isotpl_abd(41,:) /                                           &
     &     9.739e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CF4
      data isotpl_abd(42,:) /                                           &
     &     9.889e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C4H2
      data isotpl_abd(43,:) /                                           &
     &     9.560e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HC3N
      data isotpl_abd(44,:) /                                           &
     &     9.633e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2
      data isotpl_abd(45,:) /                                           &
     &     9.997e-1, 3.114e-4, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CS
      data isotpl_abd(46,:) /                                           &
     &     9.396e-1, 4.168e-2, 1.056e-2, 7.417e-3, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! SO3
      data isotpl_abd(47,:) /                                           &
     &     9.434e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
     &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
!


CONTAINS !=================== MODULE CONTAINS ==========================


   !-----------------------------------------------------------------------
   ! Given a molecule name returns the index number in an array of molecules
   ! If not found returns -1
   !-----------------------------------------------------------------------
   PURE FUNCTION molIndex( molName, molArray )
   !-----------------------------------------------------------------------
      USE Module_Utility, ONLY: upper

      integer                  :: molIndex
      character(*), intent(in) :: molName
      character(*), intent(in) :: molArray(:)

      integer :: im

      molIndex = -1
      do im=1,size(molArray)
         if ( trim(adjustl( molArray(im) )) == &
              trim(adjustl( molName )) )then
            molIndex = im
            EXIT
         endif
      enddo

   END FUNCTION


   !-----------------------------------------------------------------------
   ! To check if a molecule is a member of a set of molecules
   !-----------------------------------------------------------------------
   PURE FUNCTION molIsMember(molName, molSet)
   !-----------------------------------------------------------------------
      logical                  :: molIsMember
      character(*), intent(in) :: molName
      character(*), intent(in) :: molSet(:)

      molIsMember = .FALSE.
      if ( molIndex( molName, molSet ) >0 ) molIsMember = .TRUE.
   END FUNCTION



   !-----------------------------------------------------------------------
   ! Given a molecule name returns the index number in array clblmMolNames (HITRAN molecular code)
   ! If not found returns -1
   !-----------------------------------------------------------------------
   PURE FUNCTION molNum(molName)
   !-----------------------------------------------------------------------
      integer                  :: molNum
      character(*), intent(in) :: molName

      molNum = molIndex( molName, clblmMolNames )
   END FUNCTION



   !-----------------------------------------------------------------------
   ! Given a molecule number (HITRAN mol index) returns mol name string
   !-----------------------------------------------------------------------
   PURE FUNCTION molName( molNum )
   !-----------------------------------------------------------------------
      character(20)       :: molName !In the future, may consider using variable length character, which is a Fortran 2003 feature.
      integer ,intent(in) :: molNum


      !if (molNum<1 .or. molNum>MXMOL) then
      !   STOP ('--- molName(): Invalid molecular number! ') !Fortran pure function can not contains STOP
      !endif

      molName = trim(adjustl( clblmMolNames(molNum) ))

   END FUNCTION

   !-----------------------------------------------------------------------
   ! MATCH THE x-section molecular name AGAINST THE NAMES STORED IN ALIAS
   ! AND DETERMINE THE INDEX VALUE.
   ! If not found returns -1
   !-----------------------------------------------------------------------
   PURE FUNCTION xsMolNum( xsMolName )
   !-----------------------------------------------------------------------
      USE Module_Utility, ONLY: upper

      integer                  :: xsMolNum
      character(*), intent(in) :: xsMolName

      integer :: ind,ia

      xsMolNum = -1
      do ia =1,4
         ind = molIndex( xsMolName, ALIAS(ia,:) )
         if (ind>0) then
            xsMolNum = ind
            EXIT
         endif
      enddo

      ! NO MATCH FOUND
      ! if (xsMolNum==-1) STOP '--- xsMolNum(): Invalid x-section molecular name' !Fortran pure function can not contains STOP

   END FUNCTION


   !-----------------------------------------------------------------------
   ! Given a XS-molecular number, returns the XS molecular name in the ALIAS number ali.
   !-----------------------------------------------------------------------
   PURE FUNCTION xsMolName( xsMolNum, ali )
   !-----------------------------------------------------------------------
      integer          ,intent(in) :: xsMolNum
      integer ,optional,intent(in) :: ali
      character(20)                :: xsMolName

      integer :: ia

      ia = 1
      if (present(ali)) ia=ali

      !if (xsMolNum<1 .or. xsMolNum>MX_XS) then
      !   STOP ('--- xsMolName(): Invalid xsMolNum! ') !Fortran pure function can not contains STOP
      !endif

      !Only the name in ALIAS-I is returned.
      xsMolName = trim(adjustl( ALIAS(ia,xsMolNum) ))

   END FUNCTION


   !-----------------------------------------------------------------------
   ! Comparing the input xs molecular name and the line molecular name, if
   ! there are equal, returns .TRUE. otherwise returns .FALSE.
   !-----------------------------------------------------------------------
   PURE FUNCTION XsName_LnName_areEqual( xsName, lnName )
   !-----------------------------------------------------------------------
      USE Module_Utility, ONLY: upper

      character(*)   ,intent(in) :: xsName, lnName
      logical                    :: XsName_LnName_areEqual

      integer :: xsNo

      xsNo = xsMolNum( xsName )
      if (trim(adjustl( xsMolName( xsNo, 1 ) )) == trim(adjustl( lnName )) .or. &
          trim(adjustl( xsMolName( xsNo, 2 ) )) == trim(adjustl( lnName )) .or. &
          trim(adjustl( xsMolName( xsNo, 3 ) )) == trim(adjustl( lnName )) .or. &
          trim(adjustl( xsMolName( xsNo, 4 ) )) == trim(adjustl( lnName )) ) then

         XsName_LnName_areEqual = .TRUE.
      else
         XsName_LnName_areEqual = .FALSE.
      endif

   END FUNCTION


   !-----------------------------------------------------------------------
   ! Given molecular number and isotopologue number this function
   ! returns "isoName" which is a character string like "CO2_727"
   !-----------------------------------------------------------------------
   PURE FUNCTION molNum_isoNum_to_isoName( iMol, iISO )
   !-----------------------------------------------------------------------
      character(20)       :: molNum_isoNum_to_isoName
      integer, intent(in) :: iMol, iISO

      character(10) :: isoStr

      write(isoStr,'(I0)') isotpl_code( iMol,iISO )
      molNum_isoNum_to_isoName = trim(adjustl(clblmMolNames(iMol))) //'_'// &
                                 trim(adjustl(isoStr))

   END FUNCTION

   !-----------------------------------------------------------------------
   ! Give an iso name, trims the number after "_", returns mol name. &
   ! For example, given "CO2_727" returns "CO2".
   ! If the input molID is not a iso name (doesn't contains "_###") just
   ! return the input molID
   !-----------------------------------------------------------------------
   PURE FUNCTION trimISONum( molID )
   !-----------------------------------------------------------------------
      character(20)            :: trimISONum
      character(*), intent(in) :: molID

      integer :: i, imol, iiso
      character(len(molID)) :: tempID

      tempID = adjustl(molID)
      i=scan( tempID, '_' )
      if (i==0) then
         trimISONum = trim(adjustl(molID))
      else
         trimISONum = trim(adjustl( tempID(1:i-1) ))
      endif

   END FUNCTION


   !-----------------------------------------------------------------------
   ! Given an isotopologue "isoName" returns the molecular code, which is the
   ! HITRAN molecular code or the index number in array clblmMolNames
   ! An isoName is like "CO2_727"
   ! If not found returns -1
   !-----------------------------------------------------------------------
   PURE FUNCTION isoName2molNum( isoName )
   !-----------------------------------------------------------------------
      integer                  :: isoName2molNum
      character(*), intent(in) :: isoName

      integer :: i

      i=scan(isoName,'_')
      isoName2molNum = molNum( isoName(1:i-1) )

   END FUNCTION

   !-----------------------------------------------------------------------
   ! Given an isotopologue "isoName" returns the isotopologue number, which is the
   ! index number in the 2nd dimension of array isotpl_code.
   ! For example, isoName="CO2_727", isoNum=9
   ! If not found returns -1
   !-----------------------------------------------------------------------
   PURE FUNCTION isoName2isoNum( isoName )
   !-----------------------------------------------------------------------
      integer                  :: isoName2isoNum
      character(*), intent(in) :: isoName

      integer :: i, iMol, iISO, is

      i=scan(isoName,'_')
      iMol = molNum( isoName(1:i-1) )
      read( isoName( i+1:len(isoName) ), * ) iISO

      isoName2isoNum = -1
      do is = 1,MaxISOTPL_abd(iMol)
         if ( iISO == isotpl_code(iMol,is) )then
            isoName2isoNum = is
            EXIT
         endif
      enddo
      !if (isoName2isoNum==-1) STOP '--- isoName2isoNum(): Invalid isotopologue ID.' !Fortran pure function can not contains STOP

   END FUNCTION


   !-----------------------------------------------------------------------
   ! Check if a name is a isoName name
   !-----------------------------------------------------------------------
   PURE logical FUNCTION isISOName( molID )
   !-----------------------------------------------------------------------
      character(*), intent(in) :: molID

      integer :: i, imol, iiso

      i=scan( molID, '_' )
      if (i==0) then
         isISOName = .FALSE.
      else
         imol = isoName2molNum( molID )
         iiso = isoName2isoNum( molID )
         if ( (imol>=1 .and. imol<=MXMOL) .and. &
              (iiso>=1 .and. iiso<=MaxISOTPL_abd(imol)) ) then
            isISOName = .TRUE.
         else
            isISOName = .FALSE.
         endif
      endif

   END FUNCTION


END MODULE

