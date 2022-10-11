!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_UnitConv
   USE Module_ConstParam  ,ONLY: MXMOL, &
                                 avogad, alosmt, &
                                 airmwt, &
                                 PZERO=>Press1013,&
                                 TZERO=>Temp273,&
                                 AMWT,HMOLC,molNum
   
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: &
            molDensity2numDen ,& !(numDen, WMOL,JUNIT,molID, P,T,denWat,unitWat)
                 press2mb     ,& !(newP, P,JUNITP)
                  temp2Kelvin ,& !(newT, T,JUNITT)
                 range2Km        !(newR, R,JUNITR)
   PUBLIC :: &
            volMixRat2numDen_dry ,& !(WMOL,P,T,numDenWat); other than water vapor
            masMixRat2numDen_dry ,& !(WMOL,P,T,molWt,numDenWat); other than water vapor
               masDen2numDen     ,& !(WMOL,molWt); for all species
             parPress2numDen     ,& !(WMOL,T); for all species
            volMixRat2numDen_wat ,& !(WMOL,P,T)
            masMixRat2numDen_wat ,& !(WMOL,P,T,molWt)
               dewPtK2numDen_wat ,& !(WMOL,T,molWt)
               dewPtC2numDen_wat ,& !(WMOL,T,molWt)
               relHmd2numDen_wat    !(WMOL,T,molWt)


CONTAINS !================== MODULE CONTAINS ===========================

   
   !------------------------------------------------------------------
   ! THE FUNCTION DENSAT FOR THE SATURATION            
   ! WATER VAPOR DENSITY OVER WATER IS ACCURATE TO BETTER THAN 1       
   ! PERCENT FROM -50 TO +50 DEG C. (SEE THE LOWTRAN3 OR 5 REPORT)     
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION DENSAT(ATEMP,B)
   !------------------------------------------------------------------
      real             :: DENSAT
      real, intent(in) :: ATEMP,B
      
      real, PARAMETER :: C1 = 18.9766
      real, PARAMETER :: C2 = -14.9595
      real, PARAMETER :: C3 = -2.4388

      DENSAT = ATEMP*B*EXP(C1+C2*ATEMP+C3*ATEMP**2)*1.0E-6
   END FUNCTION

   !------------------------------------------------------------------
   ELEMENTAL FUNCTION RHOAIR(P,T)
   !------------------------------------------------------------------
      real             :: RHOAIR
      real, intent(in) :: P,T
      RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T) 
   END FUNCTION

   !------------------------------------------------------------------
   ELEMENTAL FUNCTION AAA(T)
   !------------------------------------------------------------------
      real             :: AAA
      real, intent(in) :: T
      AAA = TZERO/T 
   END FUNCTION

   !------------------------------------------------------------------
   ELEMENTAL FUNCTION BBB(molWt)
   !------------------------------------------------------------------
      real             :: BBB
      real, intent(in) :: molWt
      BBB = AVOGAD/molWt      
   END FUNCTION

   !------------------------------------------------------------------
   ELEMENTAL FUNCTION RRR(molWt)
   !------------------------------------------------------------------
      real             :: RRR
      real, intent(in) :: molWt
      RRR = AIRMWT/molWt
   END FUNCTION

   !------------------------------------------------------------------
   ELEMENTAL FUNCTION DRYAIR(P,T,numDenWat)
   !------------------------------------------------------------------
      real             :: DRYAIR
      real, intent(in) :: P,T,numDenWat
      
      DRYAIR = RHOAIR(P,T) - numDenWat
   END FUNCTION




   !------------------------------------------------------------------
   ! GIVEN VOL. MIXING RATIO                                           
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION volMixRat2numDen_dry(WMOL,P,T,numDenWat)&
             result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,P,T,numDenWat

      DENNUM = WMOL*DRYAIR(P,T,numDenWat)*1.E-6 
   END FUNCTION


   !------------------------------------------------------------------
   ! GIVEN MASS MIXING RATIO (GM KG-1)                                 
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION masMixRat2numDen_dry(WMOL,P,T,molWt,numDenWat) &
             result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,P,T,molWt,numDenWat
      
      DENNUM = RRR(molWt)*WMOL*1.0E-3*DRYAIR(P,T,numDenWat) 
   END FUNCTION


   !------------------------------------------------------------------
   ! GIVEN MASS DENSITY (GM M-3)                                       
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION masDen2numDen(WMOL,molWt) result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,molWt

      DENNUM = BBB(molWt)*WMOL*1.0E-6 
   END FUNCTION


   !------------------------------------------------------------------
   ! GIVEN PARTIAL PRESSURE (MB)                                       
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION parPress2numDen(WMOL,T) result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,T

      DENNUM = ALOSMT*(WMOL/PZERO)*(TZERO/T) 
   END FUNCTION
   !------------------------------------------------------------------


   !------------------------------------------------------------------
   ! Water vapor
   ! GIVEN VOL. MIXING RATIO                                           
   ! Convert using density of dry air.                                                                                                         
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION volMixRat2numDen_wat(inWMOL,P,T) result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: inWMOL,P,T
      
      real :: wMol

      wMOL = inWMOL*1.E-06 
      DENNUM = (wMOL/(1.+wMOL))*RHOAIR(P,T)
   END FUNCTION


   !------------------------------------------------------------------
   ! Water vapor
   ! GIVEN MASS MIXING RATIO (GM KG-1)                                                                        
   ! Convert using density of dry air.  The following quadratic is     
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION masMixRat2numDen_wat(inWMOL,P,T,molWt) &
             result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: inWMOL,P,T,molWt
      
      real :: wMol

      wMOL = inWMOL*RRR(molWt)*1.0E-3 
      DENNUM = (wMOL/(1.+wMOL))*RHOAIR(P,T)
   END FUNCTION


   !------------------------------------------------------------------
   ! Water vapor
   ! GIVEN DEWPOINT (DEG K)
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION dewPtK2numDen_wat(WMOL,T,molWt) result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,T,molWt
      
      real :: ATD

      ATD = TZERO/(WMOL) 
      DENNUM = DENSAT(ATD,BBB(molWt))*(WMOL)/T 
   END FUNCTION


   !------------------------------------------------------------------
   ! Water vapor
   ! GIVEN DEWPOINT (DEG C)                                            
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION dewPtC2numDen_wat(WMOL,T,molWt) result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,T,molWt

      real :: ATD

      ATD = TZERO/(TZERO+WMOL) 
      DENNUM = DENSAT(ATD,BBB(molWt))*(TZERO+WMOL)/T 
   END FUNCTION


   !------------------------------------------------------------------
   ! Water vapor
   ! GIVEN RELATIVE HUMIDITY (PERCENT)                                 
   !------------------------------------------------------------------
   ELEMENTAL FUNCTION relHmd2numDen_wat(WMOL,T,molWt) result(DENNUM)
   !------------------------------------------------------------------
      real             :: DENNUM
      real, intent(in) :: WMOL,T,molWt

      DENNUM = DENSAT( AAA(T),BBB(molWt) )*(WMOL/100.0) 
   END FUNCTION


!-----------------------------------------------------------------------
! SUBROUTINE molDensity2numDen() 
! COMPUTES THE DENSITY (MOL CM-3) 
!                                                                       
!        WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR       
!        UNIFORM DATA INPUT -                                           
!                                                                       
!      JCHAR    JUNIT                                                   
!                                                                       
!    " ",A       10    VOLUME MIXING RATIO (PPMV)                       
!        B       11    NUMBER DENSITY (CM-3)                            
!        C       12    MASS MIXING RATIO (GM(K)/KG(AIR))                
!        D       13    MASS DENSITY (GM M-3)                            
!        E       14    PARTIAL PRESSURE (MB)                            
!        F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY           
!        G       16     "    "     "  (TD IN T(C)) - H2O ONLY           
!        H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY     
!        I       18    AVAILABLE FOR USER DEFINITION                    
!        J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE    
!
!                                                                       
! * WMOL and denWat have to be conformable to each other if input array values
! * If input unit is mixing ration and molID is other than 
!   water vapour, denWat and unitWat must be present. Otherwise
!   denWat and unitWat can be omitted.
! * If molID is water vapour, devWat,unitWat can be omitted.
!-----------------------------------------------------------------------
! non-elemental interface:
!   SUBROUTINE molDensity2numDen( numDen,&
!                        WMOL,JUNIT,molID, P,T,denWat,unitWat)
!-----------------------------------------------------------------------
!      real,              intent(out) :: numDen(:)
!      real,              intent(in)  :: WMOL(:)
!      integer,           intent(in)  :: JUNIT
!      character(*),      intent(in)  :: molID
!      real,              intent(in)  :: P(:) !(mb)
!      real,              intent(in)  :: T(:) !(Kelvin)
!      real,    optional, intent(in)  :: denWat(:)
!      integer, optional, intent(in)  :: unitWat
!-----------------------------------------------------------------------
   ELEMENTAL SUBROUTINE molDensity2numDen( numDen,&
                        WMOL,JUNIT,molID, P,T,denWat,unitWat )
!-----------------------------------------------------------------------
      real,              intent(out) :: numDen
      real,              intent(in)  :: WMOL
      integer,           intent(in)  :: JUNIT
      character(*),      intent(in)  :: molID
      real,              intent(in)  :: P !(mb)
      real,              intent(in)  :: T !(Kelvin)
      real,    optional, intent(in)  :: denWat
      integer, optional, intent(in)  :: unitWat

      real :: numDenWat!(size(WMOL))
      real :: molWeight


      if (JUNIT == 11) then !Already in numDen(CM-3)
         numDen = WMOL
         RETURN
      endif

      molWeight = AMWT( molNum(molID) )

      if ( molNum(molID)==1 ) then !Water vapour

         call watDensity2numDen(numDen,WMOL,JUNIT,P,T )

      else !Other than water vapour

         if (JUNIT==10 .OR. JUNIT==12) then
            ! If unit is mixing ratio and  molecule is not water vapour
            ! denWat and unitWat must be present.
            if (.not.present(denWat) .or. .not.present(unitWat)) then
               numDen = -999.
               RETURN
               !STOP '--- converUnit(): Need water vapour information to convert mixing ratio to number density.'         
            else
               call watDensity2numDen(numDenWat,denWat,unitWat,P,T )
            endif
         endif

         select case (JUNIT)
            case (11); numDen = WMOL
            case (10); numDen = volMixRat2numDen_dry(WMOL,P,T,numDenWat)
            case (12); numDen = masMixRat2numDen_dry(WMOL,P,T,molWeight,numDenWat)
            case (13); numDen =    masDen2numDen    (WMOL,molWeight); 
            case (14); numDen =  parPress2numDen    (WMOL,T)
         end select

      endif !if (molNum(molID)==1) then !Water vapour

   END SUBROUTINE

!-----------------------------------------------------------------------
   ELEMENTAL SUBROUTINE watDensity2numDen(numDen,WMOL,JUNIT,P,T )
!-----------------------------------------------------------------------
      real,     intent(out) :: numDen
      real,     intent(in)  :: WMOL
      integer,  intent(in)  :: JUNIT
      real,     intent(in)  :: P !(mb)
      real,     intent(in)  :: T !(Kelvin)

      real :: watWeight

      watWeight = AMWT(molNum('H2O'))

      select case (JUNIT)
         case (11); numDen = WMOL
         case (10); numDen = volMixRat2numDen_wat(WMOL,P,T)
         case (12); numDen = masMixRat2numDen_wat(WMOL,P,T,watWeight)
         case (13); numDen =    masDen2numDen    (WMOL,watWeight); 
         case (14); numDen =  parPress2numDen    (WMOL,T)
         case (15); numDen =    dewPtK2numDen_wat(WMOL,T,watWeight)
         case (16); numDen =    dewPtC2numDen_wat(WMOL,T,watWeight)
         case (17); numDen =    relHmd2numDen_wat(WMOL,T,watWeight)
      end select

   END SUBROUTINE


!-----------------------------------------------------------------------
! PRESSURE CONVERSIONS                                              
!
!       JCHARP   JUNITP
!                                                                       
!      " ",A     10    PRESSURE IN (MB)                                 
!          B     11       "     "  (ATM)                                
!          C     12       "     "  (TORR)                               
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            
!                                                                                                                                              
!-----------------------------------------------------------------------
   ELEMENTAL SUBROUTINE press2mb(newP, P,JUNITP)
!-----------------------------------------------------------------------
      real,             intent(out) :: newP
      real,             intent(in)  :: P
      integer,          intent(in)  :: JUNITP
                                                                       
      real, PARAMETER :: PMB = 1013.25
      real, PARAMETER :: PTORR = 760.
                                                                       
      IF     (JUNITP.LE.10) THEN; newP = P !Already in mb. <10?
      ELSEIF (JUNITP == 11) THEN; newP = P*PMB         !ATM->mb
      ELSEIF (JUNITP == 12) THEN; newP = P*PMB/PTORR   !TORR->mb
      ELSE;  newP=-999.; RETURN !STOP '--- press2mb(): Invalid P unit.' 
      ENDIF 
            
   END SUBROUTINE

!-----------------------------------------------------------------------
!     TEMPERATURE COMVERSIONS                                           
!
!       JCHART   JUNITT
!                                                                       
!      " ",A     10    AMBIENT TEMPERATURE IN DEG(K)                    
!          B     11       "         "       "  " (C)                    
!          C     12       "         "       "  " (F)                    
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            
!                                                                                                                                              
!-----------------------------------------------------------------------
   ELEMENTAL SUBROUTINE temp2Kelvin(newT, T,JUNITT)
!-----------------------------------------------------------------------
      real,             intent(out) :: newT
      real,             intent(in)  :: T
      integer,          intent(in)  :: JUNITT

      real, PARAMETER :: DEGK = 273.15

      IF     (JUNITT.LE.10) THEN; newT = T !Already in Kelvin. <10?
      ELSEIF (JUNITT == 11) THEN; newT = T+DEGK 
      ELSE;  newT=-999.; RETURN !STOP '--- temp2Kelvin(): Invalid T unit.' 
      ENDIF 

   END SUBROUTINE

!-----------------------------------------------------------------------
! RANGE CONVERSIONS                                                
!
!       JCHARR   JUNITR
!                                                                       
!      " ",A     10    KM
!          B     11    M
!          C     12    CM
!                                                                                                                                              
!-----------------------------------------------------------------------
   ELEMENTAL SUBROUTINE range2Km(newR, R,JUNITR)
!-----------------------------------------------------------------------
      real,             intent(out) :: newR
      real,             intent(in)  :: R
      integer,          intent(in)  :: JUNITR

      IF     (JUNITR.EQ.10) THEN; newR = R !Already in Km
      ELSEIF (JUNITR.EQ.11) THEN; newR = R/1.E3 
      ELSEIF (JUNITR.EQ.12) THEN; newR = R/1.E5 
      ELSE;  newR=-999.; RETURN !STOP '--- range2Km(): Invalid R unit.' 
      ENDIF 

   END SUBROUTINE



END MODULE



   