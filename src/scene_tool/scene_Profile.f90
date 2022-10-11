!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_Scn_Profile

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_Profile, &
             CLBLM_Profile_init, &
             CLBLM_Profile_final, &
             CLBLM_Profile_resizeNumMol, &
             mergeStdAtmos, &
             mergeStdAtmos_xs, &
             convertMolUnit, &
             CMPALT, &
             CMPALT_profile, & !A wrapper of CMPALT accepting prfl structure
             interpAltTempOnPressGrid, &
             interpAltTempOnPressGrid_profile, & !A wrapper of interpAltTempOnPressGrid accepting prfl structures
             interpMolProfile, &
             setTopOfAtmos ,&
             scaleMolProfile, &
             readProfile_tape5, &
             read_XSect_Profile_tape5, &
             CLBLM_LayerGrid, &
             CLBLM_LayerGrid_init, &
             CLBLM_LayerGrid_final, &
             calcAltGridFromPress, &
             readRTgrid_tape5


   
   !* Units also serves as indicator for existence status
   !* Following LBLRTM, if unit<=6, means the profile is missing and 
   !  needs to be filled with the selected std. atmos.
   !* Assuming all levels for a profile has a same unit. i.e 
   !  unit doesn't change from level to level for a single profile.
   !* Q_air is used to derive dry air density when converting 
   !  cross-section mixing ratio to number density.
   TYPE :: CLBLM_Profile 
      integer                    :: nLev =0      !number of levels
      integer                    :: toaLev =0    !top level of the profile. The highest level with Z(toaLev)<HSPACE+0.001. toaLev <= nLev.
      real          ,allocatable :: Z(:)         ![nLev],level heights
      real          ,allocatable :: P(:)         ![nLev],level pressures
      real          ,allocatable :: T(:)         ![nLev],level temperatures
      integer                    :: Z_unit =0    !unit=0 means Z is not present; if=10, in Km; if =210, means Z( Km) is computed from P profile.
      integer                    :: P_unit =0    !unit=0 means P is not present; if=10, in pressure.
      integer                    :: T_unit =0    !
      integer                    :: nMol =0      !num of molecules, array dimmension.
      character(20) ,allocatable :: molID(:)     ![nMol], Molecule identifiers; Molecular name
      integer       ,allocatable :: molUnit(:)   ![nMol], unit=0 means not present
      real          ,allocatable :: Q(:,:)       ![nMol,nLev], level concentration
      real          ,allocatable :: Q_air(:)     ![nLev],air density, in unit of (molecules/cm^3)
   END TYPE   

   
   TYPE CLBLM_layerGrid
      integer           :: nLev =0      !number of grid levels
      real, allocatable :: Z(:)         ![nLev], vertical grid given in level heights
      real, allocatable :: P(:)         ![nLev], vertical grid given in level pressures
      integer           :: Z_unit =0    !unit=0 means Z is not present; if=10, in Km; if =210, means Z( Km) is computed from P profile.
      integer           :: P_unit =0    !unit=0 means P is not present; if=10, in pressure.
   END TYPE

   

CONTAINS !=================== MODULE CONTAINS ==========================

   !-------------------------------------------------------------------
   !* Vertical layering grid. 
   !* The level array index start from 1
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_LayerGrid_init(this, nLev)
   !-------------------------------------------------------------------
      type(CLBLM_LayerGrid)  ,intent(inout) :: this !out
      integer                ,intent(in)    :: nLev
           
      if (allocated(this%Z))  deallocate(this%Z)
      if (allocated(this%P))  deallocate(this%P)
         
      allocate( this%Z(nLev))   ;this%Z(:) = -999.
      allocate( this%P(nLev))   ;this%P(:) = -999.

      this%nLev = nLev
      
   END SUBROUTINE

   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_LayerGrid_final(this)
   !-------------------------------------------------------------------
      type(CLBLM_LayerGrid), intent(inout) :: this
             
      if (allocated(this%Z))  deallocate(this%Z)
      if (allocated(this%P))  deallocate(this%P)

      this%nLev = 0

   END SUBROUTINE

   
   !-------------------------------------------------------------------
   ! * in a profile objcet, prfl%Z, prfl%P and prfl%T array indices are 
   !   start from 1. This follows LBLRTM.
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Profile_init(this, nLev, nMol)
   !-------------------------------------------------------------------
      type(CLBLM_Profile)   ,intent(out) :: this
      integer               ,intent(in)  :: nLev
      integer               ,intent(in)  :: nMol
      

      this%nLev = nLev
      
      if (allocated(this%Z))        deallocate(this%Z)
      if (allocated(this%P))        deallocate(this%P)
      if (allocated(this%T))        deallocate(this%T)                     
      if (allocated(this%molID))    deallocate(this%molID)
      if (allocated(this%molUnit))  deallocate(this%molUnit)
      if (allocated(this%Q))        deallocate(this%Q)
      if (allocated(this%Q_air))    deallocate(this%Q_air)

      allocate( this%Z(nLev)) ;this%Z(:) = 0.
      allocate( this%P(nLev)) ;this%P(:) = 0.
      allocate( this%T(nLev)) ;this%T(:) = 0.
      
      allocate( this%molID( nMol))       ;this%molID(:)  = ''
      allocate( this%molUnit(nMol))      ;this%molUnit(:) = 0
      allocate( this%Q(      nMol,nLev)) ;this%Q(:,:)     = 0.         
      allocate( this%Q_air(       nLev)) ;this%Q_air(:)   = 0.
      
      this%nMol = nMol
      
   END SUBROUTINE


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Profile_final(this)
   !-------------------------------------------------------------------
      type(CLBLM_Profile),  intent(inout) :: this

      if (allocated(this%Z))        deallocate(this%Z)
      if (allocated(this%P))        deallocate(this%P)
      if (allocated(this%T))        deallocate(this%T)      
      if (allocated(this%molID))   deallocate(this%molID)
      if (allocated(this%molUnit))  deallocate(this%molUnit)
      if (allocated(this%Q))        deallocate(this%Q)
      if (allocated(this%Q_air))    deallocate(this%Q_air)

      this%nMol = 0
      this%nLev = 0

   END SUBROUTINE

   
   !-------------------------------------------------------------------
   ! Resize the molecular concentration array, unit array and name array.
   ! If the new size is greater than the old one, copy the old contents to the
   ! first 1:oldSize locations in the new array. If the new size is less 
   ! than the old one, trim the old array to the new size.
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Profile_resizeNumMol(this, newNumMol)
   !-------------------------------------------------------------------
      type(CLBLM_Profile)   ,intent(inout) :: this
      integer               ,intent(in)    :: newNumMol
      
      real          ,allocatable :: newQ(:,:)
      character(20) ,allocatable :: newMolID(:)
      integer       ,allocatable :: newMolUnit(:)
      integer                    :: nLev, nm
      
      
      nLev = this%nLev
      allocate( newQ(      newNumMol, nLev) ) ;newQ(:,:)     = 0.
      allocate( newMolID(  newNumMol) )       ;newMolID(:)  = ''
      allocate( newMolUnit(newNumMol) )       ;newMolUnit(:) = 0
      
      nm = min(this%nMol, newNumMol)
      newQ(      1:nm,:) = this%Q(      1:nm,:)
      newMolID(  1:nm)   = this%molID( 1:nm)
      newMolUnit(1:nm)   = this%molUnit(1:nm)
      
      call move_alloc( newQ,       this%Q )
      call move_alloc( newMolID,   this%molID )
      call move_alloc( newMolUnit, this%molUnit )
      
      this%nMol = newNumMol
      
   END SUBROUTINE



!--------------------------------------------------------------------
! SUBROUTINE mergeStdAtmos()
!
!       JCHAR   JUNIT                                                   
!                                                                       
!     " ",A      10    VOLUME MIXING RATIO (PPMV)                       
!         B      11    NUMBER DENSITY (CM-3)                            
!         C      12    MASS MIXING RATIO (GM(K)/KG(AIR))                
!         D      13    MASS DENSITY (GM M-3)                            
!         E      14    PARTIAL PRESSURE (MB)                            
!         F      15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY           
!         G      16     "    "     "  (TD IN T(C)) - H2O ONLY           
!         H      17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY     
!         I      18    AVAILABLE FOR USER DEFINITION                    
!        1-6    1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            
!                                                   (SEE KEY BELOW)        
!                                                                       
!     ***** OTHER 'JCHAR' SPECIFICATIONS - JCHARP,JCHART                
!                                                                       
!       JCHAR   JUNIT                                                   
!                                                                       
!      " ",A     10    PRESSURE IN (MB)                                 
!          B     11       "     "  (ATM)                                
!          C     12       "     "  (TORR)                               
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            
!                                                                          
!      " ",A     10    AMBIENT TEMPERATURE IN DEG(K)                    
!          B     11       "         "       "  " (C)                    
!          C     12       "         "       "  " (F)                    
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE            
!                                                                       
!     ***** DEFINITION OF "DEFAULT" CHOICES FOR PROFILE SELECTION ***** 
!                                                                       
!      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL          
!      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS   
!      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE                
!                                                                       
!     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)                         
!                                                                       
!      FOR MOLECULES 8-35, ONLY US STD PROFILES ARE AVIALABLE           
!      THEREFORE, WHEN  'JCHAR(K) = 1-5', JCHAR(K) WILL BE RESET TO 6   
!--------------------------------------------------------------------
   SUBROUTINE mergeStdAtmos( userProfile )
!--------------------------------------------------------------------
      !USE Module_Scn_Profile   ,ONLY: CLBLM_Profile
      IMPLICIT NONE

      type(CLBLM_Profile), intent(inout) :: userProfile
      
      integer :: im, il, numLev,numMol
      real    ,allocatable :: WMOL(:)
      integer ,allocatable :: molUNIT(:)

      
      numLev = userProfile%nLev
      numMol = userProfile%nMol

      if ( all([ userProfile%Z_unit ,&
                 userProfile%P_unit ,&
                 userProfile%T_unit ] >6 ) .and. &
           all(  userProfile%molUnit(1:numMol) >6 ) ) &
         RETURN !everything is present. no need for filling.

         
      !---
      allocate( WMOL(numMol) ) 
      allocate( molUnit(numMol) )
      
      
      !--- Fill the missing profiles with standard atmosphere.
      !
      if ( userProfile%Z_unit >6 ) then !Altitude present

         do il = 1,numLev
            molUnit(1:numMol) = userProfile%molUnit(1:numMol) !molUnit will changed in DEFALT, save them in a local array.
            call DEFALT( userProfile%Z(il)   ,& !in
                         userProfile%P(il)   ,& !inout
                         userProfile%T(il)   ,& !inout
                         userProfile%P_unit  ,& !in
                         userProfile%T_unit  ,& !in
                         WMOL(1:numMol)      ,& !inout
                         molUnit(1:numMol)   ,& !inout
                         numMol )               !in
!!test
!print*, '---layer',il
!print*, userProfile%Z(il), userProfile%P(il), userProfile%T(il), userProfile%molUnit(1:numMol)                  
!print*, WMOL(1:numMol)

            do im = 1,numMol
               if (userProfile%molUnit(im)<=6) then
                  userProfile%Q(im,il) = WMOL(im)
               endif
            enddo
         enddo

         userProfile%molUnit(1:numMol) = molUnit(1:numMol) !assuming the unit is the same for all levels.
         
      elseif ( userProfile%Z_unit <=0 .and. &
               userProfile%P_unit >6 ) then !Altitude absent, but pressure present.

         do il = 1,numLev
            molUnit(1:numMol) = userProfile%molUnit(1:numMol) !molUnit will changed in DEFALT, save them in a local array.
            call DEFALT_P( userProfile%P(il)   ,& !in
                           userProfile%T(il)   ,& !inout
                           userProfile%T_unit  ,& !in
                           WMOL(1:numMol)      ,& !inout
                           molUnit(1:numMol)   ,& !inout
                           numMol )               !in
                           
            do im = 1,numMol
               if (userProfile%molUnit(im)<=6) then
                  userProfile%Q(im,il) = WMOL(im)
               endif
            enddo
         enddo

         userProfile%molUnit(1:numMol) = molUnit(1:numMol) !assuming the unit is the same for all levels.

      else
         STOP '--- mergeStdAtmos(): Altitude and pressure can not both be missing.'
      endif
      
         
      ! --- P, T and Mol profiles are all loaded. expect Z, which may be 
      !     still missing. Need to call CMPALT() to get the Z values if it
      !     is missing.
      !
      !     if there are missing profiles:
      !       if Z and P present, call DEFALT  (interp. on Z grid), output: Z P T Q's 
      !       if Z present,       call DEFALT  (interp. on Z grid), output: Z P T Q's 
      !       if P present,       call DEFALT-P(interp. on P grid), output:   P T Q's 
           
!   end associate

      deallocate(WMOL)
      deallocate(molUnit)
      
   END SUBROUTINE

!--------------------------------------------------------------------
!     SUBROUTINE DEFALT 
!                                                                       
!     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES  
!     FROM WHICH IT WILL INTERPOLATE "DEFAULT" VALUES FOR ALTITUDE "Z"  
!                                                                       
!                                                                       
!      ***  THIS SUBROUTINE IS CALLED BY "RDUNIT" WHICH                 
!      ***  READS USER SUPPLIED INPUT PROFILES OR SINGLE VALUES         
!      ***  UNDER "MODEL = 0     " SPECIFICATIONS                       
!                                                                        
!      *** SEE DOCUMENTATION FOR CLARIFICATION ***                      
!                                                                       
!     SUBROUTINE "DEFALT"IS TRIGGERRED WHENEVER ANY ONE OF              
!     THE INPUT PARAMETERS JCHARP, JCART, (JCHAR(K),K=1,NMOL) IS = 1-6  
!                                                                       
!     FOR SIMPLICITY, ALL INTERPOLATIONS ARE DONE AT ONE TIME BECAUSE   
!     THE LAGRANGE WEIGHTS (4PT), BASED ON (ALT-Z), REMAIN UNCHANGED    
!                                                                       
!     JCHARP,JCHART AND JCHAR(K) FOR K<8 ALLOW MODEL-DEPENDENT CHOICES  
!                                                                       
!                   JCHAR=JUNIT                                         
!                                                                       
!                        1       CHOOSES TROPICAL                       
!                        2         "     MID-LATITUDE SUMMER            
!                        3         "     MID-LATITUDE WINTER            
!                        4         "     HIGH-LAT SUMMER                
!                        5         "     HIGH-LAT WINTER                
!                        6         "     US STANDARD                    
!                                                                       
!                                                                       
!     JUNIT(K) FOR K>7 CHOOSES FROM THE SINGLE TRACE CONSTITUENT        
!        PROFILES, ALL APPRORIATE FOR THE US STD ATMOSPHERE             
!                                                                       
!     ***  NOTE ***  T<0 WILL ALSO PRINT OUT A MESSAGE INDICATING       
!     ***  A POSSIBLE MISAPPLICATION OF TEMPERATURE UNITS, (K) VS (C)   
!--------------------------------------------------------------------
      SUBROUTINE DEFALT( Z, P,T,JUNITP,JUNITT, WMOL,JUNIT,NMOL ) 
!--------------------------------------------------------------------
      USE Module_Scn_ModelAtm ,ONLY: NumZMD, ALT, PMATM, TMATM, AMOL, TRAC
      USE Module_Scn_Config   ,ONLY: IPR
      IMPLICIT NONE
      
      real    ,intent(in)    :: Z
      real    ,intent(inout) :: P,T
      integer ,intent(in)    :: JUNITP, JUNITT
      real    ,intent(inout) :: WMOL(:)  !(NMOL)
      integer ,intent(inout) :: JUNIT(:) !(NMOL)
      integer ,intent(in)    :: NMOL
      

      !---Local variables
      !
      integer :: ILOWER, IUPPER
      integer :: I0,I1,I2,I3, IM, IM50 
      real    :: Z0,Z1,Z2,Z3
      real    :: DEN1,DEN2,DEN3,DEN4
      integer :: K,ITR
      integer :: MATM

      !*** 4PT INTERPOLATION FUNCTION                                    
      real :: VAL, A1,A2,A3,A4, X1,X2,X3,X4
      VAL(A1,A2,A3,A4,X1,X2,X3,X4) = A1*X1+A2*X2+A3*X3+A4*X4 

!yma       COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
!yma      &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
!yma      &              NLTEFL,LNFIL4,LNGTH4                                
!yma       COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
!yma      &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
!yma      &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL                     
!yma !                                                                       
!yma       CHARACTER*8      HMOLS 
!yma !                                                                       
!yma       COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
!yma      &               JUNITT                                             
!yma       COMMON /MLATM/ ALT(MXZMD),PMATM(MXZMD,6),TMATM(MXZMD,6),          &
!yma      &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),TST(MXZMD),  &
!yma      &               AMOLS(MXZMD,MXMOL)                                 
!yma       COMMON /MLATMC/ ATMNAM(6) 
!yma       CHARACTER*24 ATMNAM 
!yma       COMMON /TRAC/ TRAC(MXZMD,MXTRAC) 
      
      
      ILOWER = 0 
      IUPPER = 0 
      IM50 = NumZMD !=50 
      DO 10 IM = 2, IM50 
         I2 = IM 
         IF (ALT(IM).GE.Z) GO TO 20 
   10 END DO 
      I2 = IM50 
   20 I1 = I2-1 
      I0 = I2-2 
      I3 = I2+1 
      IF (I0.LT.1) GO TO 30 
      IF (I3.GT.IM50) GO TO 40 

      GO TO 60 

      !     LOWER ENDPOINT CORRECTION 
   30 CONTINUE 
      ILOWER = 1 
      I0 = I1 
      I1 = I2 
      I2 = I3 
      I3 = I3+1 
      GO TO 60 

      !     UPPER ENDPOINT CORRECTION
   40 CONTINUE 
      IUPPER = 1 
      IF (Z.GT.ALT(IM50)) GO TO 50 
      I3 = I2 
      I2 = I1 
      I1 = I0 
      I0 = I1-1 
      GO TO 60 
      
      !      UPPER ENDPOINT EXTRAPOLATION  
   50 CONTINUE 
      Z0 = ALT(I0) 
      Z1 = ALT(I1) 
      Z2 = ALT(I2) 
      Z3 = Z2+2.*(Z-Z2) 
      IUPPER = 2 
      WRITE (IPR,900) Z 

      if (JUNITP.le.6.OR.JUNITT.LE.6) STOP 'DEFAULT Z'
      do k=1,nmol
         IF (JUNIT(K).le.6)  STOP 'DEFAULT Z'
      end do

      !     LAGRANGE CONTINUATION  
   60 CONTINUE 

      !     LAGRANGE COEF DETERMINATION   
      Z1 = ALT(I1) 
      Z2 = ALT(I2) 
      Z0 = ALT(I0) 
      Z3 = ALT(I3) 
      DEN1 = (Z0-Z1)*(Z0-Z2)*(Z0-Z3) 
      DEN2 = (Z1-Z2)*(Z1-Z3)*(Z1-Z0) 
      DEN3 = (Z2-Z3)*(Z2-Z0)*(Z2-Z1) 
      DEN4 = (Z3-Z0)*(Z3-Z1)*(Z3-Z2) 
      A1 = ((Z-Z1)*(Z-Z2)*(Z-Z3))/DEN1 
      A2 = ((Z-Z2)*(Z-Z3)*(Z-Z0))/DEN2 
      A3 = ((Z-Z3)*(Z-Z0)*(Z-Z1))/DEN3 
      A4 = ((Z-Z0)*(Z-Z1)*(Z-Z2))/DEN4 

      !     TEST INPUT PARAMETERS (JUNIT'S) SEQUENTIALLY FOR TRIGGER          
      !      I.E.  JUNIT(P,T,K) = 1-6                                         
      IF (JUNITP.GT.6) GO TO 70 
      MATM = JUNITP 
      
      !     WRITE (IPR,60) Z,MATM                                             
      X1 =  LOG(PMATM(I0,MATM)) 
      X2 =  LOG(PMATM(I1,MATM)) 
      X3 =  LOG(PMATM(I2,MATM)) 
      X4 =  LOG(PMATM(I3,MATM)) 
      IF (IUPPER.EQ.2) X4 = X3+2*(X3-X2) 
      P = VAL(A1,A2,A3,A4,X1,X2,X3,X4) 
      P = EXP(P) 
   70 IF (JUNITT.GT.6) GO TO 80 
      MATM = JUNITT 

      !     WRITE (IPR,65) Z,MATM 
      X1 = TMATM(I0,MATM) 
      X2 = TMATM(I1,MATM) 
      X3 = TMATM(I2,MATM) 
      X4 = TMATM(I3,MATM) 
      T = VAL(A1,A2,A3,A4,X1,X2,X3,X4) 
   80 DO 110 K = 1, NMOL 
         IF (JUNIT(K).GT.6) GO TO 110 

         IF (K.GT.7) GO TO 90 
         MATM = JUNIT(K) 

         !     WRITE (IPR,70) K,HMOLS(K),Z,MATM 
         X1 = AMOL(I0,K,MATM) 
         X2 = AMOL(I1,K,MATM) 
         X3 = AMOL(I2,K,MATM) 
         X4 = AMOL(I3,K,MATM) 
         GO TO 100 
   90    ITR = K-7 
         MATM = 6 

         !     WRITE (IPR,70) K,HMOLS(K),Z,MATM 
         X1 = TRAC(I0,ITR) 
         X2 = TRAC(I1,ITR) 
         X3 = TRAC(I2,ITR) 
         X4 = TRAC(I3,ITR) 

  100    WMOL(K) = VAL(A1,A2,A3,A4,X1,X2,X3,X4) 
         JUNIT(K) = 10 
  110 END DO 

      RETURN 

  900 FORMAT (/,'   *** Z IS GREATER THAN 120 KM ***, Z = ',F10.3) 

      END SUBROUTINE
      
!--------------------------------------------------------------------
!     SUBROUTINE DEFALT_P
!                                                                                                                                              
!     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES  
!     FROM WHICH IT WILL INTERPOLATE "DEFAULT" VALUES FOR PRESSURE "P"  
!                                                                       
!                                                                       
!      ***  THIS SUBROUTINE IS CALLED BY "RDUNIT" WHICH                 
!      ***  READS USER SUPPLIED INPUT PROFILES OR SINGLE VALUES         
!      ***  UNDER "MODEL = 0     " SPECIFICATIONS                       
!                                                                       
!      *** SEE DOCUMENTATION FOR CLARIFICATION ***                      
!                                                                       
!     SUBROUTINE "DEFALT"IS TRIGGERRED WHENEVER ANY ONE OF              
!     THE INPUT PARAMETERS JCHARP, JCART, (JCHAR(K),K=1,NMOL) IS = 1-6  
!                                                                       
!     FOR SIMPLICITY, ALL INTERPOLATIONS ARE DONE AT ONE TIME BECAUSE   
!     THE LAGRANGE WEIGHTS (4PT), BASED ON (ALT-Z), REMAIN UNCHANGED    
!                                                                       
!     JCHARP,JCHART AND JCHAR(K) FOR K<8 ALLOW MODEL-DEPENDENT CHOICES  
!                                                                       
!                   JCHAR=JUNIT                                         
!                                                                       
!                        1       CHOOSES TROPICAL                       
!                        2         "     MID-LATITUDE SUMMER            
!                        3         "     MID-LATITUDE WINTER            
!                        4         "     HIGH-LAT SUMMER                
!                        5         "     HIGH-LAT WINTER                
!                        6         "     US STANDARD                    
!                                                                       
!                                                                       
!     JUNIT(K) FOR K>7 CHOOSES FROM THE SINGLE TRACE CONSTITUENT        
!        PROFILES, ALL APPRORIATE FOR THE US STD ATMOSPHERE             
!                                                                       
!     ***  NOTE ***  T<0 WILL ALSO PRINT OUT A MESSAGE INDICATING       
!     ***  A POSSIBLE MISAPPLICATION OF TEMPERATURE UNITS, (K) VS (C)   
!
!--------------------------------------------------------------------
      SUBROUTINE DEFALT_P( P, T,JUNITT, WMOL,JUNIT,NMOL ) 
!--------------------------------------------------------------------
      USE Module_Scn_ModelAtm ,ONLY: NumZMD, PMATM, TMATM, AMOL, TRAC
      USE Module_Scn_Config   ,ONLY: IPR
      
      IMPLICIT NONE

      real    ,intent(in)     :: P           !input P grid
      real    ,intent(inout)  :: T           !output T on the P grid 
      integer ,intent(in)     :: JUNITT      !intent model# for T
      real    ,intent(inout)  :: WMOL(:)     !(NMOL)  !output WMOL(nmol) on P grid
      integer ,intent(inout)  :: JUNIT(:)    !(NMOL) !intent model# for WMOL(nmol)
      integer ,intent(in)     :: NMOL         

      !---Local variables
      !
      integer :: ILOWER, IUPPER
      integer :: J_MDL, MATM 
      integer :: I0,I1,I2,I3, LVL, LVL_50
      real    :: P_0,P_1,P_2,P_3
      real    :: DEN1,DEN2,DEN3,DEN4
      integer :: K,ITR
      
      ! *** 4PT INTERPOLATION FUNCTION 
      real :: VAL, A1,A2,A3,A4, X1,X2,X3,X4
      VAL(A1,A2,A3,A4,X1,X2,X3,X4) = A1*X1+A2*X2+A3*X3+A4*X4 

      real :: XLOG_P
      XLOG_P =  LOG(P) 

!yma       COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
!yma      &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
!yma      &              NLTEFL,LNFIL4,LNGTH4                                
!yma       COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
!yma      &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
!yma      &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL                     
!yma !                                                                       
!yma       CHARACTER*8      HMOLS 
!yma !                                                                       
!yma       COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
!yma      &               JUNITT                                             
!yma       COMMON /MLATM/ ALT(MXZMD),PMATM(MXZMD,6),TMATM(MXZMD,6),          &
!yma      &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),TST(MXZMD),  &
!yma      &               AMOLS(MXZMD,MXMOL)                                 
!yma       COMMON /MLATMC/ ATMNAM(6) 
!yma       CHARACTER*24 ATMNAM 
!yma       COMMON /TRAC/ TRAC(MXZMD,MXTRAC) 
      
      
      DO 200 J_MDL=1,6 

         ILOWER = 0 
         IUPPER = 0 
         LVL_50 = NumZMD !LVL_50 = 50 
         DO LVL = 2, LVL_50 
            I2 = LVL 
            IF (P .GE. PMATM(LVL,J_MDL)) GO TO 20 
         ENDDO
         I2 = LVL_50 
   20    I1 = I2-1 
         I0 = I2-2 
         I3 = I2+1 
         IF (I0.LT.1) GO TO 30 
         IF (I3.GT.LVL_50) GO TO 40 

         GO TO 60 

         ! LOWER ENDPOINT CORRECTION  
   30    CONTINUE 
         ILOWER = 1 
         I0 = I1 
         I1 = I2 
         I2 = I3 
         I3 = I3+1 
         GO TO 60 

         ! UPPER ENDPOINT CORRECTION 
   40    CONTINUE 
         IUPPER = 1 
         IF (P .LE. PMATM(LVL_50,J_MDL)) GO TO 50 
         I3 = I2 
         I2 = I1 
         I1 = I0 
         I0 = I1-1 
         GO TO 60 


         ! UPPER ENDPOINT EXTRAPOLATION
   50    CONTINUE 
         P_0 =  LOG(PMATM(I0,J_MDL)) 
         P_1 =  LOG(PMATM(I1,J_MDL)) 
         P_2 =  LOG(PMATM(I2,J_MDL)) 
         P_3 = P_2+2.*(XLOG_P-P_2) 
         IUPPER = 2 
         WRITE (IPR,900) P 

         STOP 'DEFAULT P' 

         ! LAGRANGE CONTINUATION  
         !
   60    CONTINUE 

         ! LAGRANGE COEF DETERMINATION 
         P_0 =  LOG(PMATM(I0,J_MDL)) 
         P_1 =  LOG(PMATM(I1,J_MDL)) 
         P_2 =  LOG(PMATM(I2,J_MDL)) 
         P_3 =  LOG(PMATM(I3,J_MDL)) 
         DEN1 = (P_0-P_1)*(P_0-P_2)*(P_0-P_3) 
         DEN2 = (P_1-P_2)*(P_1-P_3)*(P_1-P_0) 
         DEN3 = (P_2-P_3)*(P_2-P_0)*(P_2-P_1) 
         DEN4 = (P_3-P_0)*(P_3-P_1)*(P_3-P_2) 
         A1 = ((XLOG_P-P_1)*(XLOG_P-P_2)*(XLOG_P-P_3))/DEN1 
         A2 = ((XLOG_P-P_2)*(XLOG_P-P_3)*(XLOG_P-P_0))/DEN2 
         A3 = ((XLOG_P-P_3)*(XLOG_P-P_0)*(XLOG_P-P_1))/DEN3 
         A4 = ((XLOG_P-P_0)*(XLOG_P-P_1)*(XLOG_P-P_2))/DEN4 

         ! TEST INPUT PARAMETERS (JUNIT'S) SEQUENTIALLY FOR TRIGGER          
         !  I.E.  JUNIT(P,T,K) = 1-6                                         
         !                                                                   
         ! FOR THIS VERSION OF THE SUBROUTINE DRIVEN BY PRESSURE P           
         ! JUNITP IS THE MODEL ATMOSPHERES TO BE USED FOR THE ALTITUDE       
         !                                                                   
   70    IF (JUNITT.GT.6 .OR. JUNITT.NE.J_MDL) GO TO 80 
         MATM = JUNITT 

         !WRITE (IPR,65) P_,MATM                                            

         X1 = TMATM(I0,MATM) 
         X2 = TMATM(I1,MATM) 
         X3 = TMATM(I2,MATM) 
         X4 = TMATM(I3,MATM) 
         T = VAL(A1,A2,A3,A4,X1,X2,X3,X4) 
         
   80    DO 110 K = 1, NMOL 
            IF (JUNIT(K).GT.6 .OR. JUNIT(K).NE.J_MDL) GO TO 110 
            
            IF (K.GT.7) GO TO 90 
            MATM = JUNIT(K) 
            
            !WRITE (IPR,70) K,HMOLS(K),P_,MATM                                 
                     
            X1 = AMOL(I0,K,MATM) 
            X2 = AMOL(I1,K,MATM) 
            X3 = AMOL(I2,K,MATM) 
            X4 = AMOL(I3,K,MATM) 
            GO TO 100 
   90       ITR = K-7 
            MATM = 6 
            
            !WRITE (IPR,70) K,HMOLS(K),P_,MATM                                 
            
            X1 = TRAC(I0,ITR) 
            X2 = TRAC(I1,ITR) 
            X3 = TRAC(I2,ITR) 
            X4 = TRAC(I3,ITR) 
            
  100       WMOL(K) = VAL(A1,A2,A3,A4,X1,X2,X3,X4) 
            JUNIT(K) = 10 
  110    CONTINUE 

  200 CONTINUE 

      RETURN 

  900 FORMAT (/,'   *** P IS GREATER THAN P(120 KM)  ***, P = ',        &
     &     1PE13.4)                                                     

      END  SUBROUTINE

!-----------------------------------------------------------------------
! Unit conversion for molecular density
!-----------------------------------------------------------------------
      SUBROUTINE convertMolUnit( prfl )
!-----------------------------------------------------------------------
      !USE Module_Scn_Profile     ,ONLY: CLBLM_Profile
      USE Module_UnitConv        ,ONLY: molDensity2numDen
      USE Module_ConstParam      ,ONLY: molNum
      USE Module_Utility         ,ONLY: upper
      IMPLICIT NONE
      
      type(CLBLM_Profile), intent(inout) :: prfl

      integer  :: im,nLev
      integer  :: molUnit
      real    ,allocatable :: numDen(:)

      nLev = prfl%nLev
      allocate(numDen(nLev))

      do im = 1,prfl%nMol

         molUnit = prfl%molUnit(im)
         if ( molUnit ==11 )  CYCLE !Already in number density

         ! If unit is mixing ratio and  molecule is not water vapour
         ! denWat and unitWat must be present.
         if ( (molUnit==10 .OR. molUnit==12) .AND. &
               upper(trim(adjustl(prfl%molID(im))))/='H2O' ) then

            call molDensity2numDen( numDen,&
                                    prfl%Q(im,1:nLev)       ,&
                                    prfl%molUnit(im)        ,&
                                    prfl%molID(im)         ,&
                                    prfl%P(1:nLev)          ,&
                                    prfl%T(1:nLev)          ,&
                                    prfl%Q( molNum('H2O'),1:nLev ),&      !watvap density
                                    prfl%molUnit(molNum('H2O')) )!watvap unit
         else

            call molDensity2numDen( numDen,&
                                    prfl%Q(im,1:nLev)       ,&
                                    prfl%molUnit(im)        ,&
                                    prfl%molID(im)         ,&
                                    prfl%P(1:nLev)          ,&
                                    prfl%T(1:nLev) )
         endif !do im = 1,prfl%nMol

         prfl%Q(im,1:nLev) = numDen(1:nLev)
         prfl%molUnit(im) = 11 ! now in unit of number density


      enddo !do im = 1,prfl%nMol
      
      deallocate(numDen)

      END SUBROUTINE

   
!-----------------------------------------------------------------------
! A wrapper to CMPALT accepting CLBLM_Profile structure as input/output argument.
!-----------------------------------------------------------------------
      SUBROUTINE CMPALT_profile( prfl, latitude, earthRadius )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: molNum
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      IMPLICIT NONE
      
      type(CLBLM_Profile) ,intent(inout) :: prfl
      real                ,intent(in)    :: latitude
      real                ,intent(in)    :: earthRadius
      integer  :: nLev      
      
      nLev = prfl%nLev
      
      if (prfl%Z_unit <=0) then
         call CMPALT( nLev,                                 & !ILVL
                      prfl%P(1:nLev),                       & !PM in (mb)
                      prfl%T(1:nLev),                       & !TM in (K)
                      prfl%Q( molNum('H2O'),1:nLev ),       & !DENW
                      prfl%Z(1),                            & !REF_Z in (Km)
                      latitude,                             & !REF_LAT
                      earthRadius,                          & !RE
                      prfl%Z(1:nLev) )                        !ZMDL in (Km)         

         prfl%Z_unit = 210 !(km) !This is determined by CMPALT. "210" means this is calculated from P grid.
      endif
   
      END SUBROUTINE

      
!-----------------------------------------------------------------------
!
!     AUTHOR: TONY CLOUGH, JENNIFER DELAMERE, JOHN WARDEN               
!             JANUARY 2001                                              
!     PROGRAM TO CALCULATE ALTITUDE LEVEL (ZMDL) GIVEN                  
!     PRESSURE (PM), TEMPERATURE (TM) AND THE NUMBER DENSITY            
!     OF WATER VAPOR (DENW) USING THE HYDROSTATIC EQUATION              
!                                                                       
!     INPUT:                                                            
!      A) PRESSURE (MBAR)                                               
!      B) TEMPERATURE (KELVIN)                                          
!      C) NUMBER DENSITY OF WATER VAPOR                                 
!                                                                       
!     OUTPUT:                                                           
!      IDEAL GAS LAW: P.E.CIDDOR (1996), Refractive index of 
!      air: New equations for the visible and near infrared,
!      Applied Optics, 35(9), 1566-1573.
!-----------------------------------------------------------------------
      SUBROUTINE CMPALT(ILVL,PM,TM,DENW,REF_Z,REF_LAT,RE, ZMDL)
!-----------------------------------------------------------------------                                                                        
      USE Module_ConstParam ,ONLY: BOLTZ, GASCON
      USE Module_ConstParam ,ONLY: XMASS_DRY, GRAV_CONST
      USE Module_Scn_Config ,ONLY: IPR
      
      IMPLICIT NONE
      
      integer ,intent(in)  :: ILVL       !         Number of levels 
      real    ,intent(in)  :: PM(:)      !(ILVL)   (mb), Pressure profile
      real    ,intent(in)  :: TM(:)      !(ILVL)   (K), Temperature profile
      real    ,intent(in)  :: DENW(:)    !(ILVL)   (Molec/cm^3), Water vapor profile
      real    ,intent(in)  :: REF_Z      !         (Km), Reference altitude
      real    ,intent(in)  :: REF_LAT    !         (Deg), Reference latitude
      real    ,intent(in)  :: RE         !         (Km), The Earth radius.
      real    ,intent(out) :: ZMDL(:)    !(ILVL)   (Km), Calculated altitude for output.

           
      real ,PARAMETER :: CA0 = 1.58123E-6, &
                         CA1 = -2.9331E-8, &
                         CA2 = 1.1043E-10 
      real ,PARAMETER :: CB0 = 5.707E-6, &
                         CB1 = -2.051E-8
      real ,PARAMETER :: CC0 = 1.9898E-4, &
                         CC1 = -2.376E-6
      real ,PARAMETER :: CD = 1.83E-11, &
                         CE = -0.0765E-8                                                                        
      real ,PARAMETER :: XMASS_H2O = 0.018015


      INTEGER :: I,J,K
      REAL    :: XMASS_RATIO
      REAL    :: H2O_MIXRAT(ILVL),COMP_FACTOR(ILVL),ZTEMP(ILVL)                 
      REAL    :: G0, GAVE
      REAL    :: Y 
      REAL    :: CHI0, DCHI
      REAL    :: T0,DT
      REAL    :: TOTAL_AIR, DRY_AIR, CHIM
      REAL    :: C1,C2,C3 
      REAL    :: A, B, ALPHA 
      REAL    :: XINT_TOT 

!yma       COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
!yma      &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
!yma      &              NLTEFL,LNFIL4,LNGTH4                                
!yma                                                                         
!yma       COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
!yma      &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
!yma      &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL                     
!yma !                                                                       
!yma       REAL PM(MXZMD),TM(MXZMD),DENW(MXZMD),ZMDL(MXZMD) 
!yma       REAL H2O_MIXRAT(MXZMD),COMP_FACTOR(MXZMD),ZTEMP(MXZMD) 
!yma                                                                         
!yma       REAL Y 
!yma       REAL CHI0 
!yma       REAL T0,DT 
!yma       REAL C1,C2,C3 
!yma       REAL A, B, ALPHA 
!yma !      REAL BTZ                                                         
!yma       REAL XINT_TOT 
!yma                                                                         
!yma       DATA CA0/1.58123E-6/,CA1/-2.9331E-8/,CA2/1.1043E-10/ 
!yma       DATA CB0/5.707E-6/,CB1/-2.051E-8/ 
!yma       DATA CC0/1.9898E-4/,CC1/-2.376E-6/ 
!yma       DATA CD/1.83E-11/,CE/-0.0765E-8/ 
!yma                                                                         
!yma       DATA XMASS_H2O/0.018015/
      
      !--- CALCULATE GRAVITY AT REFERENCE LATITUDE AT SURFACE                                                                        
      G0 = GRAV_CONST(REF_LAT)

      !---
      ! CALCULATE THE NUMBER DENSITY OF TOTAL AIR MOLECULES [MOLEC/CM^3]      
      ! CALCULATE THE COMPRESSIBILITY FACTOR (COMP_FAC) FOR THE               
      ! IDEAL GAS LAW                                                         
      XMASS_RATIO = XMASS_H2O/XMASS_DRY 
      DO J=1,ILVL 
         DT             = TM(J) - 273.15 
         TOTAL_AIR      = PM(J)*1.0E+3/(BOLTZ*TM(J)) 
         DRY_AIR        = TOTAL_AIR - DENW(J) 
         H2O_MIXRAT(J)  = DENW(J)/DRY_AIR 
         CHIM           = XMASS_RATIO*H2O_MIXRAT(J) 
         COMP_FACTOR(J) = 1. - ( PM(J)*100/TM(J) )* &
                          ( CA0 + CA1*DT + CA2*DT**2 + &
                            ( CB0 + CB1*DT )*CHIM + &
                            ( CC0 + CC1*DT )*CHIM**2 ) + &
                          ( CD + CE*CHIM**2 )*( PM(J)*100./TM(J) )**2                        
      ENDDO 

      !--- CONVERT REFERENCE ALTITUDE TO METERS                                  
      ZTEMP(1) = REF_Z*1000.0 
      ZMDL(1) = REF_Z 
      
      DO 20 I=1, ILVL - 1 
         GAVE = G0*(RE/(RE+ZTEMP(I)/1000.0))**2 
         Y = LOG(PM(I+1)/PM(I)) 

         IF (Y .NE. 0.0) THEN 
            CHI0 = H2O_MIXRAT(I) 
            DCHI = (H2O_MIXRAT(I+1)-H2O_MIXRAT(I))/Y 

            T0 = TM(I) 
            DT = (TM(I+1) - TM(I))/Y 

            C1 = T0 + T0*CHI0 
            C2 = T0*DCHI + DT*CHI0 + DT 
            C3 = DT*DCHI 

            B = 1 + XMASS_RATIO*CHI0 
            A = XMASS_RATIO*DCHI 
            ALPHA = A/B 

            IF ( ABS(ALPHA*Y) .GE. 0.01) THEN 
               write(ipr,*) 'LAYER ',I, &
                 ' THICKER THAN IDEAL FOR ALTITUDE CALCULATION'
!v128               PRINT*,'LAYER TOO THICK' 
!v128               STOP 
            ENDIF 

            XINT_TOT = C1*Y + 0.5*(C2-C1*ALPHA)*Y**2 + &
                       0.3333*( C3-C2*ALPHA+C1*ALPHA**2 )*Y**3                                                                       
            XINT_TOT = -XINT_TOT*(GASCON*1.0E-7)/(XMASS_DRY*GAVE*B) 

            ZTEMP(I+1) = ZTEMP(I) + XINT_TOT*COMP_FACTOR(I) 
            ZMDL(I+1) = ZTEMP(I+1)/1000. 
         ELSE 
            ZTEMP(I+1) = ZMDL(I)*1000.0 
            ZMDL(I+1) = ZMDL(I) 
         ENDIF 
   20 ENDDO

      END SUBROUTINE

      
      
!-----------------------------------------------------------------------
! A wrapper to interpAltTempOnPressGrid_profile accepting CLBLM_Profile structures
! as input/output argument. Interpolate the Z and T from xPrfl to prfl grid.
!-----------------------------------------------------------------------
   SUBROUTINE interpAltTempOnPressGrid_profile( xPrfl, prfl, latitude, earthRadius )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: molNum
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      IMPLICIT NONE
      
      type(CLBLM_Profile) ,intent(inout) :: xPrfl
      type(CLBLM_Profile) ,intent(in)    :: prfl
      real                ,intent(in)    :: latitude
      real                ,intent(in)    :: earthRadius
      

      if (xPrfl%Z_unit <=0) then
      
         call interpAltTempOnPressGrid( xPrfl%Z(1:xPrfl%nLev),& !out
                                        xPrfl%T(1:xPrfl%nLev),& !out
                                        xPrfl%P(1:xPrfl%nLev),&
                                        prfl%Z(1:prfl%toaLev),&
                                        prfl%P(1:prfl%toaLev),&
                                        prfl%T(1:prfl%toaLev),&
                                        prfl%Q( molNum('H2O'),1:prfl%toaLev ),&
                                        xPrfl%nLev,&
                                        prfl%toaLev,&
                                        latitude, &
                                        earthRadius, &
                                        prfl%Z_unit==210 ) !if prfl%Z if computed from prfl%P, ignore the interpolated Z in "calcAltTempOnPresGrid()".
      
         xPrfl%Z_unit = 210 !(km)  "210" means Z is computed from P.
         xPrfl%T_unit = 10 !(K)
      endif
      

!yma moved the check out of the if-block      
      !IF (IP .NE. 1) THEN 
      !   IF (ZX(IP).LE.ZX(IP-1)) GO TO 300 
      !ENDIF 
      if (any( xPrfl%Z(2:xPrfl%nLev) <= xPrfl%Z(1:xPrfl%nLev-1) )) then
         STOP '--- interpAltTemp(): Profile altitudes must be in ascending order!'
         !STOP 'ZX IN XPROFL' 
      endif
               
   END SUBROUTINE
    
!-----------------------------------------------------------------------
!
! * TO ENSURE THAT CALCULATED/INPUT ZMDL'S WILL MATCH CALCULATED USER-LEVE
!   ALTITUDES, A COMBINATION OF INTERPOLATION AND HYDROSTATICS ARE USED.  
!   alt = A * F1(P) + (1 - A) * F2(P), WHERE                             
!   F1(P) = INTERPOLATION IN LN(P), F2(P) = HYDROSTATIC CALCULATION
! * v12.7: Ignore exponential interpolation term (zint) if user has input 
!   a profile specified on pressure levels
! * pres and alt are array arguments.                            
! * It returns temperatures at pressure grid as well.
! * "interpAltTempOnPressGrid()" differ from "CMPALT()" in that the latter 
!   assume the available of Temperature and Watervapor at pressure grid.
!   "interpAltTempOnPressGrid()" obtains the Temperature and Watervapor by 
!   interpolation from TM and DENW.
!
!-----------------------------------------------------------------------
      SUBROUTINE interpAltTempOnPressGrid( alt, Temp, &
                        pres,ZMDL,PM,TM,DENW,IBMAX,IMMAX,REF_LAT,RE, &
                        ignoreInterpZ )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out)   :: alt(:)        !(IBMAX)   out
      real    ,intent(out)   :: Temp(:)       !(IBMAX)   out
      real    ,intent(in)    :: pres(:)       !(IBMAX)   boundary level in pressure 
      real    ,intent(in)    :: ZMDL(:)       !(IMMAX)   model atm. Altitude
      real    ,intent(in)    :: PM(:)         !(IMMAX)   model atm. Pressure
      real    ,intent(in)    :: TM(:)         !(IMMAX)   model atm. Temperature
      real    ,intent(in)    :: DENW(:)       !(IMMAX)   model atm. water vapour density.
      integer ,intent(in)    :: IBMAX         !          number of boundary levels
      integer ,intent(in)    :: IMMAX         !          number of model atm. levels
      real    ,intent(in)    :: REF_LAT       !          reference latitude
      real    ,intent(in)    :: RE            !          earth radius
      logical ,intent(in)    :: ignoreInterpZ !          If ZMDL was computed from P profile, Altitude will be calculated from hydrostatic equation only.


      integer :: ISTART, IP, LIP
      real    :: ZTMP(2),PTMP(2),TTMP(2),WVTMP(2)
      real    :: HIP,ZINT
      real    :: TIP
      real    :: WVIP
      real    :: RATP,A

      
      ISTART = 2 

      DO 160 IP=1,IBMAX

         PTMP(1)  = 0.0 
         TTMP(1)  = 0.0 
         WVTMP(1) = 0.0 
         ZTMP(1)  = 0.0 

         PTMP(2)  = 0.0 
         TTMP(2)  = 0.0 
         WVTMP(2) = 0.0 
         ZTMP(2)  = 0.0 

         !what if pres(1) lower than (pres(1)>PM(1))the PM(1)? 

         DO 161 LIP=ISTART,IMMAX 
            IF (pres(IP) .GT. PM(LIP)) GO TO 162 
  161    CONTINUE 
         LIP=IMMAX 
  162    CONTINUE
  
         IF (pres(IP) .EQ. PM(LIP-1)) THEN 
            alt(IP) = ZMDL(LIP-1) 
            Temp(IP) = TM(LIP-1) 
         ELSE 

            IF(pres(IP) .EQ. PM(LIP)) THEN 
               alt(IP) = ZMDL(LIP) 
               Temp(IP) = TM(LIP) 
            ELSE 
               ! PERFORM INTERPOLATION IN LN(PM)                                       
               HIP = (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))                                         
               ZINT = ZMDL(LIP-1)+ HIP* LOG(pres(IP)/PM(LIP-1)) 

               ! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION               
               PTMP(1) = PM(LIP-1) 
               ZTMP(1) = ZMDL(LIP-1) 
               TTMP(1) = TM(LIP-1) 
               WVTMP(1) = DENW(LIP-1) 

               PTMP(2) = pres(IP) 

               TIP = (TM(LIP)-TM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1) )
               TTMP(2) = TM(LIP-1)+ TIP* LOG(pres(IP)/PM(LIP-1) )

               WVIP = (DENW(LIP)-DENW(LIP-1))/ LOG(PM(LIP)/PM(LIP-1)) 
               WVTMP(2) = DENW(LIP-1) + WVIP* LOG(pres(IP)/PM(LIP-1)) 
               CALL CMPALT(2,PTMP,TTMP,WVTMP,ZTMP(1),REF_LAT,RE, ZTMP) 
                                                               
               ! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION
               !
               RATP = LOG(pres(IP)/PM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
               
               !A = RATP**3 
               !v12.7: Ignore exponential interolation term (zint) if user has input a profile specified on pressure levels
               if ( ignoreInterpZ ) then 
                  A =0.
               else 
                  A = RATP**3 
               endif
               
               alt(IP) = A*ZINT + (1-A)*ZTMP(2) 
               Temp(IP) = TTMP(2) 
            ENDIF 
         ENDIF 
         
         ISTART = LIP 

  160 CONTINUE !DO 160 IP=1,IBMAX

      END SUBROUTINE

      
!-----------------------------------------------------------------------
!  Determine the TOA height of the profile based the value of scnGeom%HSPACE
!-----------------------------------------------------------------------  
   SUBROUTINE setTopOfAtmos( prfl, HSPACE )
!-----------------------------------------------------------------------  
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      IMPLICIT NONE
   
      type(CLBLM_Profile) ,intent(inout) :: prfl
      real                ,intent(in)    :: HSPACE

      integer :: ISPACE, il

      !--- 
      ! Check if the input profile altitudes are in a ascending order
      do il = 2,prfl%nLev
         if (prfl%Z(il) <= prfl%Z(il-1)) then
            STOP '--- setTopOfAtmos():INPUT ALTITUDES FOR LBLRTM LAYERS ARE NOT IN ASCENDING ORDER '
         endif
      enddo
      
      !--- Find the top level based on HSPACE. 
      !    The top level is the highest level that has &
      !    Z(toaLev)<HSPACE+0.001. prfl%toaLev <= prfl%nLev
      DO il = 1,prfl%nLev
         IF ( HSPACE+0.001 .GT. prfl%Z(il)) ISPACE = il
      END DO 
      prfl%toaLev = ISPACE 
      
      
      !!--- Zero out molecular amounts above the TOA level
      !if ( prfl%nLev > prfl%toaLev ) then
      !   prfl%Q( :, prfl%toaLev+1 : prfl%nLev ) = 0.
      !   prfl%Q_air(prfl%toaLev+1 : prfl%nLev ) = 0.
      !endif

   END SUBROUTINE
      

      
!-----------------------------------------------------------------------  
! FOR EACH MOLECULE K FOR WHICH JCHAR(K,ILEV) IS '1', THIS SUBROUTIN
! INTERPOLATES THE MIXING RATIO DTMP(K,ILEV) AT THE ALTITUDE Z      
! FROM THE STANDARD PROFILE IN AMOLX ON THE ALTITUDE GRID ALTX.     
! * This subroutine is modified from XTRACT()
!-----------------------------------------------------------------------  
   SUBROUTINE mergeStdAtmos_xs( xPrfl )       
!----------------------------------------------------------------------- 
      USE Module_ConstParam     ,ONLY: xsMolNum
      USE Module_Utility        ,ONLY: EXPINT
      USE Module_Scn_ModelAtm   ,ONLY: NumZMD, ALTX=>ALT, AMOLX
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      ! AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH  
      ! LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL
      
      IMPLICIT NONE

      type(CLBLM_Profile) ,intent(inout) :: xPrfl

      
      integer :: ISTART, iLev, iL, K, IXINDX
      real    :: altZ, A

      
      ISTART = 2
      
      DO iLev = 1,xPrfl%nLev 
         
         altZ = xPrfl%Z(iLev)
      
         ! FIND SMALLEST ALTX(iL) GT altZ  
         DO iL = ISTART, NumZMD 
            IF (altZ.LE.ALTX(iL)) GO TO 20 
         ENDDO 
         iL = NumZMD 
   20    CONTINUE 
         
         DO K = 1, xPrfl%nMol
         
            !IF (JCHAR(K,iLev).EQ.'1') THEN !INTERPOLATE MIXING RATIO FROM STANDARD PROFILE 
            IF (xPrfl%molUnit(K).EQ.1) THEN !INTERPOLATE MIXING RATIO FROM STANDARD PROFILE 
            
               IXINDX = xsMolNum( xPrfl%molID(K) )
            
               A = (altZ-ALTX(iL-1))/(ALTX(iL)-ALTX(iL-1)) 
               CALL EXPINT( xPrfl%Q(K,iLev),&
                            AMOLX(iL,IXINDX),&
                            AMOLX(iL-1, IXINDX), A )
            ENDIF 
            
         ENDDO !DO K = 1, xPrfl%nMol
         
         ISTART = iL
         
      ENDDO !DO iLev = 1,xPrfl%nLev
            
      END SUBROUTINE
         
!-----------------------------------------------------------------------  
!                                                                       
!  This code modified from XINTERP(). Molecular profiles in xPrfl are 
!  interpolated to altitude grid from prfl. 
!
!  THIS SUBROUTINE INTERPLOLATES THE PROFILE DENX ON THE ALTITUDE    
!  GRID ZX INTO DENM ON THE GRID ZMDL.  EXPONENTIAL INTERPOLATION    
!  IS USED. 
!
!-----------------------------------------------------------------------  
      SUBROUTINE interpMolProfile( xPrfl, prfl ) 
!-----------------------------------------------------------------------  
      USE Module_ConstParam     ,ONLY: molNum, &
                                       ALOSMT, PZERO=>Press1013, TZERO=>Temp273
      USE Module_Utility        ,ONLY: EXPINT
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      USE Module_Scn_Config     ,ONLY: IPR      
      IMPLICIT NONE
      
      type(CLBLM_Profile) ,intent(inout) :: xPrfl
      type(CLBLM_Profile) ,intent(in)    :: prfl
      
      
      integer             :: L,IMMAX, LX,LAYX, K,IXMOLS
      real                :: A
      real ,allocatable   :: ZX(:),DENX(:,:)
      real ,allocatable   :: ZMDL(:),DENM(:,:),DRYAIR(:)
      type(CLBLM_Profile) :: newXPrfl
      

      
      IMMAX  = prfl%toaLev
      LAYX   = xPrfl%nLev
      IXMOLS = xPrfl%nMOL
      
      allocate( ZX(          LAYX ) ) 
      allocate( DENX( IXMOLS,LAYX ) ) 
      allocate( ZMDL(        IMMAX) ) 
      allocate( DENM( IXMOLS,IMMAX) ) 
      allocate( DRYAIR(      IMMAX) ) 
            
      ZX(1:LAYX)            = xPrfl%Z(1:LAYX)
      DENX(1:IXMOLS,1:LAYX) = xPrfl%Q(1:IXMOLS,1:LAYX)
      ZMDL(1:IMMAX)         = prfl%Z(1:IMMAX)
      
      !--- Dry air denisty in (molecules/cm^3)
      if ( all( abs(prfl%Q_air(:)) <tiny(0.) ) ) then
         DRYAIR(1:IMMAX) = ALOSMT * (prfl%P(1:IMMAX)/PZERO) / &
                                    (prfl%T(1:IMMAX)/TZERO)
      else
         DRYAIR(1:IMMAX) = prfl%Q_air( 1:IMMAX ) - &
                           prfl%Q( molNum('H2O'),1:IMMAX )
      endif

      
      LX = 2 
      DO 30 L = 1, IMMAX 

         ! FIND THE SMALLEST ZX GE ZMDL(L)
   10    CONTINUE 
         IF (ZMDL(L).LE.ZX(LX).OR.LX.EQ.LAYX) THEN
         
            A = (ZMDL(L)-ZX(LX-1))/(ZX(LX)-ZX(LX-1)) 
            IF (A.LT.0.0 .OR. A.GT.1.0) WRITE (IPR,900) 

            !!IF DRYAIR FOR LAYER NOT CALCULATED PREVIOUSLY 
            !!(USING NORMAL MOLECULES), THEN CALCULATE THE 
            !!NUMBER DENSITY OF AIR
            !IF (DRYAIR(L).EQ.0.) THEN
            !   DRYAIR(L) = ALOSMT*(PM(L)/PZERO)/(TM(L) /TZERO)
            !ENDIF

            DO K = 1, IXMOLS 
               CALL EXPINT (DENM(K,L),DENX(K,LX-1),DENX(K,LX),A) 

               !---(YMa,181015): Extrapolation may result in negative density. 
               if (DENM(K,L)<0.) DENM(K,L)=0.
               
               !CONVERT MIXING RATIO (PPMV) TO NUMBER DENSITY
               DENM(K,L) = DRYAIR(L)*DENM(K,L)*1.0E-6                
            ENDDO
            
            CYCLE !GO TO 30 
         ELSE
         
            LX = LX+1 
         ENDIF
         
         GO TO 10 
   30 ENDDO 
   

      call CLBLM_Profile_init( newXPrfl, prfl%nLev, IXMOLS )
      newXPrfl%toaLev              = prfl%toaLev
      newXPrfl%Z(1:IMMAX)          = prfl%Z(1:IMMAX)
      newXPrfl%P(1:IMMAX)          = prfl%P(1:IMMAX)
      newXPrfl%T(1:IMMAX)          = prfl%T(1:IMMAX)
      newXPrfl%Z_unit              = prfl%Z_unit
      newXPrfl%P_unit              = prfl%P_unit
      newXPrfl%T_unit              = prfl%T_unit
      newXPrfl%molID(1:IXMOLS)     = xPrfl%molID(1:IXMOLS)
      newXPrfl%molUnit(1:IXMOLS)   = 11 !number density
      newXPrfl%Q(1:IXMOLS,1:IMMAX) = DENM(1:IXMOLS,1:IMMAX)
      newXPrfl%Q_air(1:IMMAX)      = prfl%Q_air(1:IMMAX)
      
      xPrfl = newXPrfl
      
      !!--- zero the profiles above IMMAX
      !if (xPrfl%nLev >IMMAX) then
      !   xPrfl%Q(:, xPrfl%nLev+1:IMMAX) = 0.
      !endif
      
      deallocate( ZX,DENX,ZMDL,DENM,DRYAIR ) 
      
  900 FORMAT (//,'  XINTPL: CAUTION- EXTRAPOLATING X-SECTION PROFILE') 
      END SUBROUTINE




!-----------------------------------------------------------------------
! Scale the profile based on user input integrated column amounts
! * User input molecular names and corresponding column amounts that need to be scaled.
! * Only those molecules that to be scaled are needed to be present in molID array (and in columnAmt and scaleMode).
! * Molecular names should be the same as the user used in scene data.
! * Molecular concentration has to be in unit of number density (mol/cm-3), 
!   columnAmt can be in column number density (mol/cm-2), Dobson unit, ppv, column precipitable water (cm), 
!   or it can be a no-unit scaling factor.
! * molID, columAmt and scaleMod are arrays
! * The scaling factor is the same at all levels
! * earthRadius is an optional input.
!
! * The value of scaleMode is interperated ands used as follows:
!
!   'l',or 'L' or 
!    number '1'     scaling factor used directly to scale profile
!   'c' or 'C'      column amount in mol/cm^2 units to which the profile is to be scaled
!   'd' or 'D'      column amount in Dobson units to which the profile is to be scaled
!   'm' or 'M'      volume mixing ratio (ppv) wrt dry air for the total column to which the profile will be scaled!                   
!   'p' or 'P'      value of Precipitable Water Vapor (cm) to which the profile will be scaled (water vapor only)
!
! * Example:
!     molID     = ['H2O', 'CO2', 'O3 ' ]
!     scaleAmt  = [ 0.5 ,   2.,   200. ]
!     scaleMode = [ 'c' ,  '1',   'd'  ]
!
!-----------------------------------------------------------------------  
   SUBROUTINE scaleMolProfile( prfl, molID, scaleAmt, scaleMode, earthRadius )
!-----------------------------------------------------------------------  
      USE Module_Utility        ,ONLY: upper
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      IMPLICIT NONE
   
      type(CLBLM_Profile) ,intent(inout) :: prfl
      character(*)        ,intent(in)    :: molID(:)     !Mol name for molecules to be scaled. Not all mol in prfl, just thost to be scaled.
      real                ,intent(in)    :: scaleAmt(:)   !Column amounts or scale factors for molecules to be scaled.
      character(1)        ,intent(in)    :: scaleMode(:)  !To specify how to interpreted and use the value of columnAmt(:)
      real      ,OPTIONAL ,intent(in)    :: earthRadius   !Needed by profile integration subroutine to determine the minimum thickness of sublayer.

      real              :: RE
      integer           :: im, il, K, Loc, numSclMol
      character(20)     :: molK
      real ,allocatable :: sclAmt(:)
      real              :: totAmt(1), Wair(1), Wh2o(1), Wdry, fac

      
      
      !--- Set RE 
      if (.not.present(earthRadius)) then
         RE = 6371.23
      else
         RE = earthRadius
      endif

      !--- Molecular profiles need to be in number density unit
      if ( any(prfl%molUnit(:) /=11) ) then
         STOP '--- scaleMolProfile(): Gas concentration needs to be in unit of umber density.'
      endif


      numSclMol = size(molID)
      allocate( sclAmt( numSclMol ) )      
      
      !--- Convert the unit of colmumAmt to column number density
      do im =1, numSclMol
         select case (scaleMode(im))
         
            case ('1','l','L') !scaling factor used directly to scale profile
               sclAmt(im) = - scaleAmt(im) !use negative sign as a marker
               
            case ('c','C') !column amount to which the profile is to be scaled (molec/cm^2)            
               sclAmt(im) = scaleAmt(im)
               
            case ('d','D') !column amount in Dobson units to which the profile is to be scaled
               sclAmt(im) = scaleAmt(im)*2.68678e16
               
            case ('m','M') !volume mixing ratio (ppv) wrt dry air for the total column to which the profile will be scaled
   
               ! obtain dry air sum                                             
               !               
               CALL calcTotColumAmount( prfl%Z, &
                                        prfl%Q( 1, 1:prfl%nLev ), & !H2O is assumed to be the first molecule in the array.
                                        Wh2o, &
                                        NMOL=1, &
                                        nLev=prfl%nLev,&
                                        RE=RE )
               CALL calcTotColumAmount( prfl%Z, &
                                        prfl%Q_air( 1:prfl%nLev ), &
                                        Wair, &
                                        NMOL=1, &
                                        nLev=prfl%nLev,&
                                        RE=RE )
               Wdry = Wair(1) - Wh2o(1)
               
               sclAmt(im) = scaleAmt(im)*Wdry
            
            case ('p','P') !value of Precipitable Water Vapor (cm) to which the profile will be scaled (water vapor only)

               if ( trim(adjustl( molID(im) )) /='H2O' ) then
                  STOP '--- scaleMolProfile(): Precipitable water input is allowed for H2O only.'
               endif

               sclAmt(im) = scaleAmt(im)/2.99150e-23; !value from LBLRTM
               
            case default
               STOP '--- scaleMolProfile(): Invalid value for scaleMode.'
         end select
      enddo
      
      

      !--- Start scaling
      !            
      do im = 1, numSclMol
      
         !--- Find location of molecule to be scaled
         Loc = 0         
         do K=1,prfl%nMol
         
            molK = prfl%molID(K)
            
            if ( upper(trim(adjustl( molID(im) ))) == &
                 upper(trim(adjustl( molK )))  ) then
              
               Loc = K
               exit
            endif
         enddo
                  
         !--- No matching found, something is wrong.
         if ( Loc==0 ) then
            STOP '--- scaleMolProfile(): The molecule to be scaled not found!'
         endif
         
         
         !--- Location found, scale the profile
         !
         if ( Loc /=0 ) then
         
            CALL calcTotColumAmount( prfl%Z, prfl%Q(Loc,:), totAmt, &
                                     NMOL=1, nLev=prfl%nLev, RE=RE )
                                     
            if (sclAmt(im)<0) then 
               fac = abs(sclAmt(im))
            else
               fac = sclAmt(im) / totAmt(1)
            endif
            
            do il = 1,prfl%nLev
               prfl%Q( Loc,il ) = fac * prfl%Q( Loc,il )
            enddo
         endif
         
!print*, '---scaleLine',totAmt
      enddo !do im = 1, numSclMol
      
      deallocate( sclAmt )            
   END SUBROUTINE
   

!-----------------------------------------------------------------------
! * Absorbing gas density is assumed to follow an exponential distribution. 
! * Layer temperature is assumed to be isothermal.
! * DENP in unit of mol/cm-3, totAmt in unit of mol/cm-2
!-----------------------------------------------------------------------
      SUBROUTINE calcTotColumAmount( ZPTH, DENP, totAmt, NMOL, nLev, RE )
!-----------------------------------------------------------------------
      USE Module_ConstParam, ONLY: r8=>kind_r8, GCAIR, DELTAS
      IMPLICIT NONE
            
      real     ,intent(in)    :: ZPTH(nLev)  
      real     ,intent(in)    :: DENP(nMol,nLev) !assumed-shape array doesn't work, since user my input one single profile
      real     ,intent(out)   :: totAmt(nMol)
      integer  ,intent(in)    :: nMol
      integer  ,intent(in)    :: nLev
      real     ,intent(in)    :: RE

      
      real ,PARAMETER :: EPSILN=1.0E-5
      
      INTEGER           :: K, J
      REAL(r8)          :: DH, DHMIN, DS, R1
      REAL              :: Z1,Z2, H1,H2,H3, DZ, dsdz
      real ,allocatable :: DENA(:), DENB(:), HDEN(:), AMTP(:,:)


      
      allocate( DENA(NMOL) )
      allocate( DENB(NMOL) )
      allocate( HDEN(NMOL) )
      allocate( AMTP(nMol,nLev) )      
      AMTP(:,:) = 0.0
      
      DO J = 1,nLev-1
      
         !---  INITIALIZE VARIABLES FOR THE CALCULATION OF THE PATH
         !
         Z1    = ZPTH(J) 
         Z2    = ZPTH(J+1) 
         H1    = Z1 
         R1    = RE+H1 
         DHMIN = DELTAS**2/(2.0*R1) 
         DZ    = Z2-Z1
         DO K = 1, NMOL 
            DENA(K) = DENP(K,J) 
            DENB(K) = DENP(K,J+1) 
            IF ((DENA(K).EQ.0.0.OR.DENB(K).EQ.0.0) .OR. &
                (ABS(1.0-DENA(K)/DENB(K)).LE.EPSILN) ) THEN  
                
               ! Use linear interpolation                                           
               HDEN(K) = 0.0
            ELSE  
         
               ! Use exponential interpolation          
               HDEN(K) = -DZ/ LOG(DENB(K)/DENA(K))
            ENDIF 
         END DO
         
         
         ! LOOP THROUGH vertical PATH                                                 
         ! INTEGRATE PATH QUANTITIES USING QUADRATIC INTEGRATION WITH        
         ! UNEQUALLY SPACED POINTS
         !
         DO while(.TRUE.)
         
            DH = DELTAS
            DH = MAX(DH,DHMIN) 
            H3 = H1+DH 
            IF (H3.GT.Z2) H3 = Z2 
            DH = H3-H1 
            DS = DH
            dsdz = 1.                    
            DO K = 1, NMOL 
               IF ((HDEN(K).EQ.0.0).OR. (ABS(DH/HDEN(K)).LT.EPSILN)) THEN 
            
                  ! Linear interpolation                                  
                  ! 1.0E05 factor converts units km to cm                 
                  DENB(K)=DENP(K,J)+(DENP(K,J+1)-DENP(K,J))*(H3-Z1)/DZ 
                  AMTP(K,J) = AMTP(K,J)+0.5*(DENA(K)+DENB(K))*DS*1.0E5 
               ELSE 
            
                  ! Exponential interpolation                             
                  DENB(K) = DENP(K,J)*EXP(-(H3-Z1)/HDEN(K)) 
                  AMTP(K,J) = AMTP(K,J)+DSDZ*HDEN(K) *(DENA(K)-DENB(K))*1.0E5                                                    
               ENDIF 
            ENDDO   
            DO K = 1, NMOL 
               DENA(K) = DENB(K) 
            ENDDO
            
            IF (H3.LT.Z2) THEN 
               H1 = H3 
            ELSE 
               exit
            ENDIF 
            
         ENDDO !end of do while (.true.)
      
      ENDDO !Do J=1,nLev-1
      
      !--- Sum up the layer amounts to get total column.
      do K = 1,nMol
         totAmt(K) = sum( AMTP(K,1:nLev-1) )
      enddo   
         
!print*, AMTP(1,:)
!print*, '---kk', K, nMol,nlev
!print*, sum( AMTP(1,1:nLev-1) )
!   
      deallocate( DENA )
      deallocate( DENB )
      deallocate( HDEN )
      deallocate( AMTP )

   END SUBROUTINE
      

      
      
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   SUBROUTINE readProfile_tape5( prfl, tape5 )
   !-----------------------------------------------------------------------
      USE Module_Scn_TAPE5         ,ONLY: CLBLM_TAPE5
      !USE Module_Scn_Profile       ,ONLY: CLBLM_Profile
      
      type(CLBLM_TAPE5)      ,intent(in)  :: tape5
      type(CLBLM_Profile)    ,intent(out) :: prfl
      
      integer :: MDL, nLnMol
            
     
      if (tape5%rec_3_1%MODEL == 0) then !user provided atmosphere
      
         call getProfile_userAtm_tape5( prfl, tape5 )
         
      else ! standard atmosphere 
      
         MDL    = tape5%rec_3_1%MODEL
         nLnMol = tape5%rec_3_1%NMOL
         !nXsMol = tape5%rec_3_7%IXMOLS
         
         ! This subroutine get profile for line molecules only.
         call getProfile_stdAtm( prfl, MDL, nLnMol )
         
      endif
            
      !--- Check the input
      call grossCheck_profile( prfl )
            
   END SUBROUTINE

   

   !-----------------------------------------------------------------------
   ! * This code assumes that the user input profiles exist for 
   !   all parameters (except altitude). 
   ! * This code assums the units don't change from layer to layer
   !
   !      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL          
   !      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS   
   !      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE                
   !                                                                       
   !     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)                         
   !                                                                       
   !      FOR MOLECULES 8-35, ONLY US STD PROFILES ARE AVIALABLE           
   !      THEREFORE, WHEN  'JCHAR(K) = 1-5', JCHAR(K) WILL BE RESET TO 6   
   !-----------------------------------------------------------------------
   SUBROUTINE getProfile_userAtm_tape5( prfl, tape5 )
   !-----------------------------------------------------------------------
      !USE Module_Scn_Profile     ,ONLY: CLBLM_Profile, &
      !                                  CLBLM_Profile_init
      USE Module_Scn_TAPE5       ,ONLY: CLBLM_TAPE5
      USE Module_ConstParam      ,ONLY: KMXNOM, clblmMolNames, molNum,&
                                        ALOSMT, PZERO=>Press1013, TZERO=>Temp273
      IMPLICIT NONE

      type(CLBLM_Profile)     ,intent(inout) :: prfl !out
      type(CLBLM_TAPE5)       ,intent(in)    :: tape5

      integer :: im
      integer :: nLev, nLnMol

      
      nLev = abs(tape5%rec_3_4%IMMAX_B)
      nLnMol = tape5%rec_3_1%NMOL
      !nXsMol = tape5%rec_3_7%IXMOLS
      IF (nLnMol.EQ.0) nLnMol = KMXNOM

      !--- Allocate memory for the profile structure
      ! * prfl%Z, prfl%P, prfl%T array indices start from 1
      ! * 160504: For this version, the cross-section profile input is 
      !   following the LBLRTM approach from a TAPE5 scene file.
      !   It may need to be interpolated in altitude grid before stored  
      !   together with line-molecule profiles. Here the memory is 
      !   allocated for the final or after-interpolation x-section profiles
      !   which is of dimension (nXsMol, nLev).
      !   This subroutine get profile for line molecules only.
      call CLBLM_Profile_init( prfl, nLev, nLnMol )


      !--- Copy the input profile to a structure.
      !* Assume tape5%rec_3_4%JCHARP(:), 
      !         tape5%rec_3_4%JCHART(:), and
      !         tape5%rec_3_4%JCHAR(:,:)
      !  have a same unit value for each profile.
      !* For boundary altitude ZMDL(km). If IMMAX < 0, altitude levels are
      !  computed from pressure levels PM. If any altitude levels are
      !  provided, they are ignored if IMMAX < 0 (exception: The
      !  first input level must have an accompanying ZM for input
      !  into the hydrostatic equation)
      !
      prfl%molID(1:nLnMol) = clblmMolNames(1:nLnMol)

      prfl%Z(1:nLev) = tape5%rec_3_5_and_3_6%ZMDL(1:nLev) !ZMDL(:) may be missing
      prfl%P(1:nLev) = tape5%rec_3_5_and_3_6%PM(  1:nLev)
      prfl%T(1:nLev) = tape5%rec_3_5_and_3_6%TM(  1:nLev)
      prfl%Q(1:nLnMol,1:nLev) = tape5%rec_3_5_and_3_6%WMOL(1:nLnMol,1:nLev)

      prfl%Q_air(1:nLev) = & !air density in (molecules/cm^3)
                ALOSMT * (prfl%P(1:nLev)/PZERO) * (TZERO/prfl%T(1:nLev))

                
      ! Unit must be the same for whole profile.
      !
      if (any( tape5%rec_3_5_and_3_6%JCHARP(:) /= &
               tape5%rec_3_5_and_3_6%JCHARP(1) )) then
         STOP '--- getProfile_userAtm_tape5(): Pressure unit must be the same for whole profile. Program Stopped.'
      endif
      
      if (any( tape5%rec_3_5_and_3_6%JCHART(:) /= &
               tape5%rec_3_5_and_3_6%JCHART(1) )) then
         STOP '--- getProfile_userAtm_tape5(): Temperature unit must be the same for whole profile. Program Stopped.'
      endif

      do im = 1,nLnMol
         if (any( tape5%rec_3_5_and_3_6%JCHAR(im,:) /= &
                  tape5%rec_3_5_and_3_6%JCHAR(im,1) )) then
            STOP '--- getProfile_userAtm_tape5(): Molecular unit must be the same for whole profile. Program Stopped.'
         endif
      enddo
      
      prfl%P_unit = JOU(tape5%rec_3_5_and_3_6%JCHARP(1)) !assuming the units don't change from layer to layer
      prfl%T_unit = JOU(tape5%rec_3_5_and_3_6%JCHART(1)) !assuming the units don't change from layer to layer
      do im = 1,nLnMol
         prfl%molUnit(im) = JOU(tape5%rec_3_5_and_3_6%JCHAR(im,1)) !assuming the units don't change from layer to layer
      enddo


      ! Boundary altitude (km). If IMMAX < 0, altitude levels are
      ! computed from pressure levels PM. If any altitude levels are
      ! provided, they are ignored if IMMAX < 0 (exception: The
      ! first input level must have an accompanying ZM for input
      ! into the hydrostatic equation)
      if (tape5%rec_3_4%IMMAX_B <0) then
         prfl%Z_unit = 0 !Altitude will be computed from Pressure
      else
         prfl%Z_unit = 10 !Altitude present, in units of 'Km'
      endif
      
      !--- If x-section data present, load the 
      !if (nXsMol>0) then
      !endif
      
   END SUBROUTINE

   !-----------------------------------------------------------------------
   !
   ! "getProfile_stdAtm()" loads the prfl object from standard atmospheres.
   ! It is different from "getProfile_userAtm()" where user has to, at least,
   ! input "Z" or 'P' profiles and the rest profiles can be filled with 
   ! standard atmosphere. Here all profiles are loaded from standard atmosphere.
   !
   ! * assuming the units don't change from layer to layer
   !
   !-----------------------------------------------------------------------
   SUBROUTINE getProfile_stdAtm( prfl, MDL, nLnMol )
   !-----------------------------------------------------------------------
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile, &
      !                                 CLBLM_Profile_init
      USE Module_ConstParam     ,ONLY: clblmMolNames, molNum, &
                                       KMXNOM,MXZMD,MXMOL,MXTRAC
      USE Module_Scn_ModelAtm   ,ONLY: NumZMD, ALT, PMATM, TMATM, AMOL, TRAC
                                    
      IMPLICIT NONE

      type(CLBLM_Profile) ,intent(inout) :: prfl !out
      integer             ,intent(in)    :: MDL  ! atm model#
      integer             ,intent(in)    :: nLnMol !number of line molecules

            
      integer :: I,K,ITR, im
      integer :: nLev, kLnMol
      real    :: DRYAIR
      
      
      nLev = NumZMD !=50
      kLnMol = nLnMol
      IF (kLnMol.EQ.0) kLnMol = KMXNOM

      !--- Allocate memory for the prfl structure
      ! * prfl%Z, prfl%P, prfl%T array indices start from 1
      call CLBLM_Profile_init( prfl, nLev, kLnMol )


      !--- load the prfl structure
      !
      prfl%molID(1:kLnMol) = clblmMolNames(1:kLnMol)

      prfl%Z(1:nLev) = ALT( 1:nLev )
      prfl%P(1:nLev) = PMATM( 1:nLev, MDL )
      prfl%T(1:nLev) = TMATM( 1:nLev, MDL )
      
      do I = 1,nLev
      
         ! Calculate water density and subtract from
         ! total density to obtain dry air density 
         DRYAIR = AMOL(I,8,MDL) - AMOL(I,1,MDL)*AMOL(I,8,MDL)*1.0E-6  !H2O is in the first place.
         
         DO K = 1, 7 
            IF (K.GT.kLnMol) EXIT
            prfl%Q(K,I) = AMOL(I,K,MDL)*1.0E-6*DRYAIR 
         ENDDO

         DO K = 8, 28 
            IF (K.GT.kLnMol) EXIT
            ITR = K-7 
            ! TRAC is the trace constituent information, obtained from LBLLOW 
            prfl%Q(K,I) = TRAC(I,ITR)*1.0E-6*DRYAIR
         ENDDO
         
      enddo
      
      prfl%Q_air(1:nLev) = AMOL(1:nLev,8,MDL) !air density in (molecues/cm^3)
      
      prfl%Z_unit = 10 !(Km)
      prfl%P_unit = 10 !(mb)
      prfl%T_unit = 10 !(K)
      prfl%molUnit(1:kLnMol) = 11 !NUMBER DENSITY (CM-3)
   
   END SUBROUTINE


   
   !----------------------------------------------------------------------------
   ! This subroutine check the before-processing profile. 
   !----------------------------------------------------------------------------
   SUBROUTINE grossCheck_profile( prfl )
   !----------------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: MXZMD,MXMOL,clblmMolNames
      USE Module_Utility        ,ONLY: upper
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile
      
      type(CLBLM_Profile) ,intent(in) :: prfl
      
      character(*) ,parameter :: procedureName = 'grossCheck_profile()'
      integer :: im

      
   
      if ( prfl%nLev <=0 .or. prfl%nLev >MXZMD ) then
         STOP '--- '//procedureName//': Error in number of levels. '//' Program stopped.'
      endif
      
      if ( .NOT.allocated(prfl%Z) .or. &
           .NOT.allocated(prfl%P) .or. &
           .NOT.allocated(prfl%T) ) then
         STOP '--- '//procedureName//': Profile Z(:) or P(:) or T(:) not allocated.'//' Program stopped.'
      endif
      
      if ( .NOT.( (prfl%P_unit >=1 .and. prfl%P_unit <=6) .or. &  !default to value for specified model atmosphere
                  prfl%P_unit ==10 .or. &                        !pressure in (mb)
                  prfl%P_unit ==11 .or. &                        !(atm)
                  prfl%P_unit ==12 ) ) then                      !(torr)
         STOP '--- '//procedureName//': Invalid value for P_unit.'//' Program stopped.'
      endif
      
      if ( prfl%P_unit==0 .and. prfl%nMol>0 ) then
         STOP '--- '//procedureName//': P(:) must be present unless for cross-section profile.'//' Program stopped.'
      endif

      if ( .NOT.( (prfl%T_unit >=1 .and. prfl%T_unit <=6) .or. &  !default to value for specified model atmosphere
                  prfl%T_unit ==10 .or. &                        !(K)
                  prfl%T_unit ==11 ) ) then                      !(C)
         STOP '--- '//procedureName//': Invalid value for T_unit.'//' Program stopped.'
      endif

      if ( prfl%T_unit==0 .and. prfl%nMol>0 ) then
         STOP '--- '//procedureName//': T(:) must be present unless for cross-section profile.'//' Program stopped.'
      endif
      
      if ( .NOT.( prfl%Z_unit ==0 .or. &    !Z(:) not present and will be computed from P(:)
                  prfl%Z_unit ==10 ) ) then  !Z(:) present and is in (Km)
         STOP '--- '//procedureName//': Invalid value for Z_unit.'//' Program stopped.'
      endif

      if ( prfl%Z_unit ==0 .and. &
           .NOT.( prfl%P_unit ==10 .or. (prfl%P_unit>=1 .and. prfl%P_unit<=6)) ) then
         STOP '--- '//procedureName//': Z(:) is to be calculated from P(:), P(:) must be in (mbar) or be one of standard profile.'//' Program stopped.'
      endif
      
      !if ( prfl%Z_unit ==0 .and. prfl%Z(1) ) then
      !   STOP '--- '//procedureName//': Z(:) to be computed from P(:), Z(1) must be present.'//' Program stopped.'
      !endif
      
      if ( prfl%P_unit >6 .and. any(prfl%P(:)<0.) ) then !P(:) is present and not default to standard profile.
         STOP '--- '//procedureName//': Pressure value less then zero.'//' Program stopped.'
      endif
      
      if ( prfl%T_unit==10 .and. any(prfl%T(:)<0.) ) then !input T(:) in (K)
         STOP '--- '//procedureName//': Temperature value less than zero (K).'//' Program stopped.'
      endif


                  
      if ( .NOT.allocated(prfl%Q_air) ) then
         STOP '--- '//procedureName//': Q_air(:) not allocated.'//' Program stopped.'
      endif

      if ( prfl%nMol>0 .and. &
           ( .NOT.allocated(prfl%molID) .or. &
             .NOT.allocated(prfl%molUnit) .or. &
             .NOT.allocated(prfl%Q) ) ) then
         STOP '--- '//procedureName//': Memory not allocated for molecular arrays.'//' Program stopped.'
      endif
      

      
      if ( prfl%nMol >0 ) then 
      
         !do im =1,prfl%nMol
         !   if (.NOT.any( upper(trim(adjustl(clblmMolNames(:)))) == &
         !                 upper(trim(adjustl(prfl%molID(im)))) ) ) then
         !      STOP '--- '//procedureName//': Invalid molecule name.'//' Program stopped.'
         !   endif
         !enddo
      
      
         do im = 1,prfl%nMol
            if ( .NOT.( prfl%molUnit(im) >=1 .or. &
                        prfl%molUnit(im) <=6 .or. &  !default to value for specified model atmosphere
                        prfl%molUnit(im) ==10 .or. & !volume mixing ratio (ppmv):
                        prfl%molUnit(im) ==11 .or. & !number density (cm-3)
                        prfl%molUnit(im) ==12 .or. & !mass mixing ratio (gm/kg)
                        prfl%molUnit(im) ==13 .or. & !mass density (gm m-3)
                        prfl%molUnit(im) ==14 .or. & !partial pressure (mb)
                        (prfl%molUnit(im) ==15 .and. im==1) .or. & !dew point temp (K) *H2O only
                        (prfl%molUnit(im) ==16 .and. im==1) .or. & !dew point temp (C) *H2O only
                        (prfl%molUnit(im) ==17 .and. im==1) & !relative humidity (percent) *H2O only
                        !prfl%molUnit(im) ==18 .or. & ! available for user definition
                        ) ) then
               STOP '--- '//procedureName//': Invalid value for molecular unit.'//' Program stopped.'
            endif
         enddo
         
         
         do im = 1,prfl%nMol
            if ( prfl%molUnit(im) ==10 .AND. &  !volume mixing ratio (PPMV)
                 (any(prfl%Q(im,:)*1e-6 <0.) .or. &
                  any(prfl%Q(im,:)*1e-6 >1.)) ) then  
               STOP '--- '//procedureName//': Column density in PPMV unit, but values*1e-6 <0 or >1.'//' Program stopped.'
            endif
         enddo
         
         
         do im = 1,prfl%nMol
            if ( prfl%molUnit(im) ==12 .AND. &  !mass mixing ratio (gm/kg)
                 (any(prfl%Q(im,:)*1e-3 <0.) .or. &
                  any(prfl%Q(im,:)*1e-3 >1.)) ) then  
               STOP '--- '//procedureName//': Column density in gm/kg unit, but values*1e-3 <0 or >1.'//' Program stopped.'
            endif
         enddo

      endif !if ( prfl%nMol >0 ) then
            
   
   END SUBROUTINE

   

   !-----------------------------------------------------------------------
   ! * This code assums the units don't change from layer to layer
   !
   !      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL          
   !      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS   
   !      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE                
   !                                                                       
   !     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)                         
   !                                                                       
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   SUBROUTINE read_XSect_Profile_tape5( xPrfl, tape5 )
   !-----------------------------------------------------------------------
      !USE Module_Scn_Profile    ,ONLY: CLBLM_Profile, &
      !                                 CLBLM_Profile_init
      USE Module_Scn_TAPE5      ,ONLY: CLBLM_TAPE5
      USE Module_ConstParam     ,ONLY: xsMolNum, MXZMD, MX_XS
      USE Module_Scn_ModelAtm   ,ONLY: NumZMD, ALTX=>ALT, AMOLX
      !AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH  
      !LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL.
      
      IMPLICIT NONE

      type(CLBLM_TAPE5)     ,intent(in)  :: tape5
      type(CLBLM_Profile)   ,intent(out) :: xPrfl

      integer :: im, K
      integer :: nXLev, nXMol
      integer :: stdOrUser  !=1 standard x-section profile; =0 user-input profile.      

   
      stdOrUser = tape5%rec_3_7%IPRFL

      if (stdOrUser == 1) then
         nXLev = NumZMD !=50
      else
         nXLev = tape5%rec_3_8%LAYX
      endif
      
      nXMol = tape5%rec_3_7%IXMOLS
      
      !---
      ! * prfl%Z, prfl%P, prfl%T array indices start from 1
      ! * For a x-section profile object, line molecular components are empty.
      call CLBLM_Profile_init( xPrfl, nXLev,  nXMol ) 

      
      if (stdOrUser ==1) then !standard x-section profile
                       
         xPrfl%molID(1:nXMol) = tape5%rec_3_7_1%XSNAME(1:nXMol)
         
         xPrfl%Z(1:nXLev) = ALTX( 1:nXLev )
         xPrfl%Z_unit = 10 !input altitude must be in (km)
         xPrfl%P_unit = 0
         
         do im = 1,nXMol
            K = xsMolNum( xPrfl%molID(im) )
            xPrfl%Q(im,1:nXLev) = AMOLX(1:nXLev,K) !AMOLX are in unit of VOLUME MIXING RATIO (ppmv)
         enddo         
         xPrfl%molUnit(1:nXMol) = 10 !VOLUME MIXING RATIO (PPMV)

      else
      
         !--- Copy the input profile to "xPrfl" object.
         !* Assuming JCHAR(:,:) have a same unit value for all levels
         !
         xPrfl%molID(1:nXMol) = tape5%rec_3_7_1%XSNAME(1:nXMol)
         
         if (tape5%rec_3_8%IZORP ==0) then !=0 Altitude grid
            xPrfl%Z(1:nXLev) = tape5%rec_3_8_1%ZX(1:nXLev)
            xPrfl%Z_unit = 10 !input is altitude in (km)
            xPrfl%P_unit = 0
         else
            xPrfl%P(1:nXLev) = tape5%rec_3_8_1%PX(1:nXLev)
            xPrfl%P_unit = 10 !input is pressure in (mb), Z will be computed.
            xPrfl%Z_unit = 0
         endif
         
         !prfl%T(1:nLev) = 
         !prfl%T_unit = 
         
         xPrfl%Q(1:nXMol,1:nXLev) = tape5%rec_3_8_1%DTMP(1:nXMol,1:nXLev)                  
         do im = 1,nXMol
            xPrfl%molUnit(im) = JOU(tape5%rec_3_8_1%JCHAR(im,1)) !assuming the units don't change from layer to layer
         enddo

         
      endif !if (stdOrUser ==1)

   END SUBROUTINE


!-----------------------------------------------------------------------
! Check the unit of the input layering grid, if it is given in pressure, 
! calculate the altitude grid in Km.
!-----------------------------------------------------------------------
   SUBROUTINE calcAltGridFromPress( layGrid, prfl, REF_LAT, RE )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: molNum
      !USE Module_Scn_Profile    ,ONLY: CLBLM_LayerGrid, &
      !                                 CLBLM_Profile
      !USE Module_Scn_Profile    ,ONLY: interpAltTempOnPressGrid
      
      IMPLICIT NONE

      type(CLBLM_LayerGrid)  ,intent(inout) :: layGrid
      type(CLBLM_Profile)    ,intent(in)    :: prfl
      real                   ,intent(in)    :: REF_LAT
      real                   ,intent(in)    :: RE !earth radius used by CMPALT()
            
      integer :: nBnd,mxLev
      logical :: ZMDL_fromPM
      real ,allocatable :: dummy(:)

      
      nBnd = layGrid%nLev
      mxLev = prfl%toaLev
      ZMDL_fromPM = (prfl%Z_unit == 210) !v12.7
      allocate(dummy(nBnd))
      
      !---
      ! if boundary is given in pressure, calculate Z and T by 
      ! using hydrostatic equation and interpolation
      if ( layGrid%P_unit==10 ) THEN

         call interpAltTempOnPressGrid( layGrid%Z(1:nBnd),& !out
                                        dummy(    1:nBnd),& !out
                                        layGrid%P(1:nBnd),&
                                        prfl%Z(1:mxLev),&
                                        prfl%P(1:mxLev),&
                                        prfl%T(1:mxLev),&
                                        prfl%Q( molNum('H2O'),1:mxLev ),&
                                        layGrid%nLev,&
                                        mxLev,&
                                        REF_LAT, &
                                        RE, &
                                        ZMDL_fromPM)
         
         layGrid%Z_unit = 210 !(km)  "210" means Z is computed from P grid.
         
      endif

      

      !--- Check to make sure the lower boundary is within profile range.
      IF ( layGrid%Z(1).LT.prfl%Z(1)) THEN 

         IF (ABS(layGrid%Z(1)-prfl%Z(1)).LE.0.0001) THEN 
            layGrid%Z(1) = prfl%Z(1) 
         ELSE 
            !PRINT 946,ZBND(1),ZMDL(1) 
            !IF (NOPRNT.GE.0) WRITE (IPR,946) ZBND(1),ZMDL(1) 
            STOP ' BOUNDARIES OUTSIDE OF ATMOS' 
         ENDIF 

      ENDIF 

   END SUBROUTINE

   

   !-----------------------------------------------------------------------
   ! * Layering grid  ( or boundaries ) must be input in (mb) or (km) unit.
   !-----------------------------------------------------------------------
   SUBROUTINE readRTgrid_tape5( layGrid, tape5 )
   !-----------------------------------------------------------------------
      !USE Module_Scn_Profile      ,ONLY: CLBLM_LayerGrid, &
      !                                   CLBLM_LayerGrid_init
      USE Module_Scn_TAPE5        ,ONLY: CLBLM_TAPE5
      IMPLICIT NONE

      type(CLBLM_LayerGrid) ,intent(inout) :: layGrid !out
      type(CLBLM_TAPE5)     ,intent(in)    :: tape5
            
      integer :: ib, nBnd


      nBnd = abs(tape5%rec_3_1%IBMAX_B)
      if ( nBnd ==0 ) then
         STOP '--- readRTgrid_tape5(): Auto layering (IBMAX=0) is not supported in CLBLM.'
      endif

      !--- call layGrid%init( nBnd )
      ! * layGrid array index start from 1.
      call CLBLM_LayerGrid_init( layGrid, nBND )
      
      !---  load layGrid data from tape5
      if ( tape5%rec_3_1%IBMAX_B <0 ) then
         layGrid%P(1:nBnd) = tape5%rec_3_3b%layerBoundaries(1:nBnd)
         layGrid%P_unit = 10 !must be in (mb)
         layGrid%Z_unit = 0  !IBMAX_B<0, P will be interpolated to Z grid
      else
         layGrid%Z(1:nBnd) = tape5%rec_3_3b%layerBoundaries(1:nBnd)
         layGrid%Z_unit = 10 !must be in (km)
         layGrid%P_unit = 0
      endif


      !--- layer boundaries must be in ascending order
      !
      if ( layGrid%Z_unit >=1 ) then
         do ib = 2,nBnd
           if (layGrid%Z(ib) .LE. layGrid%Z(ib-1)) then
              print*, '--- readRTgrid_tape5(): BOUNDARY ALTITUDES FOR LBLRTM LAYERS ARE NEGATIVE OR NOT IN ASCENDING ORDER'
              STOP !GO TO 300
           endif
         enddo
      endif
      
      if ( layGrid%P_unit >=1 ) then
         do ib = 2,nBnd
            if (layGrid%P(ib) .GE. layGrid%P(ib-1)) then
               print*,'--- readRTgrid_tape5(): BOUNDARY PRESSURES FOR LBLRTM LAYERS ARE POSITIVE OR NOT IN DESCENDING ORDER'
               STOP !GO TO 305 
            endif
         enddo
      endif
      
   END SUBROUTINE
   
   
   
   !--------------------------------------------------------------------
   !     FUNCTION JOU (CHAR) 
   !     Convert profile unit in letters to integer values.
   !--------------------------------------------------------------------
   FUNCTION JOU (CHAR) 
   !--------------------------------------------------------------------                                                                      
      integer                  :: JOU
      character(1), intent(in) :: CHAR

     ! COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     !&              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     !&              NLTEFL,LNFIL4,LNGTH4,IBRD                              
      integer :: INDX,i

      CHARACTER*1 :: HOLVEC(22) 
      DATA (HOLVEC(I),I=1,22) /                                         &
     &                '1','2','3','4','5','6','0','0','0','0',' ','A',  &
     &                'B','C','D','E','F','G','H','I','J','K'/          
      integer :: INDX1(22) 
      DATA (INDX1(I),I=1,22) /                                          &
     &                  1,  2,  3,  4,  5,  6,  0,  0,  0,  0, 10, 10,  &
     &                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20/          


      INDX = 0 
      DO 10 I = 1, 22 
         IF (HOLVEC(I).NE.CHAR) GO TO 10 
         INDX = INDX1(I) 
         GO TO 20 
   10 END DO 
   20 IF (INDX.EQ.0) THEN 
         !WRITE (IPR,900) CHAR 
         STOP ' JOU: BAD PARAM ' 
      ENDIF 
      JOU = INDX 
                                                                       
      RETURN 
!                                                                       
!  900 FORMAT ('0 INVALID PARAMETER :',2X,A1) 
!                                                                       
      END FUNCTION                                       
   
   
!--------------------------------------------------------------------
!      UNITS CONVERSION FOR P AND T                          
!                                                                       
!     A = P OR T     AND  IA =JUNITP(I.E. MB,ATM,TORR)                  
!                            =JUNITT(I.E. DEG K OR C)                   
!                            =JUNITR(I.E. KM,M,OR CM)                   
!--------------------------------------------------------------------
   SUBROUTINE CHECK_P_T_unit(A,IA,KEY) 
!--------------------------------------------------------------------
      IMPLICIT NONE
      
      real,    intent(inout) :: A
      integer, intent(in)    :: IA
      integer, intent(in)    :: KEY
      
      real, PARAMETER :: PMB=1013.25
      real, PARAMETER :: PTORR=760.
      real, PARAMETER :: DEGK=273.15

      
      IF (IA.LE.10) RETURN 

      GO TO (10,20,30) KEY 

      !     PRESSURE CONVERSIONS                                              
   10 IF (IA.EQ.11) THEN 
         A = A*PMB 
         RETURN 
      ELSEIF (IA.EQ.12) THEN 
         A = A*PMB/PTORR 
         RETURN 
      ELSE 
         STOP ' CHECK_P_T_unit(P)' 
      ENDIF 

      !     TEMPERATURE COMVERSIONS                                           
   20 IF (IA.LE.11) THEN 
         A = A+DEGK 
         RETURN 
      ELSE 
         STOP ' CHECK_P_T_unit(T)' 
      ENDIF 

      !      RANGE CONVERSIONS                                                
   30 IF (IA.EQ.11) THEN 
         A = A/1.E3 
         RETURN 
      ELSEIF (IA.EQ.12) THEN 
         A = A/1.E5 
         RETURN 
      ELSE 
         STOP ' CHECK_P_T_unit(R)' 
      ENDIF 

   END SUBROUTINE

   
      
      
END MODULE
