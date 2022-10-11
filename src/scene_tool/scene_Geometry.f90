!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_Scn_Geometry

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_ScnGeom,&
             calcH1H2FromPress, &
             readSceneGeometry_tape5

   TYPE :: CLBLM_ScnGeom
      real    :: Hobs   =0.    !(Km or mb)  observer altitude in Km or mb 
      real    :: Hend   =0.    !(Km or mb)  target altitude in Km or mb
      integer :: H_unit =0     !if <0, H1,H2 are given as pressures in mb; >0 in Km.
      real    :: obsAng =0.    !(deg) path zenith angle at H1
      integer :: LEN    =0     !(NA)  Short/long path through a tangent height. need when H1>H2 and ANGLE>90.
      !real    :: RANGE         !(Km)  Path length  for horizontal path
      !real    :: BETA          !(deg) Earth centerd angle from H1 to H2
      integer :: ITYPE  =-999  != 2  slant path from H1 to H2, use RECORD 3.2; = 3  slant path from H1 to space
         
      real  :: HSPACE      =0.  !altitude definition for space (default = 100 km).
      real  :: earthRadius =0.  !
      real  :: longitude   =0.  !
      real  :: latitude    =0.  !(deg) This is needed when one want to fill up missing with AFGL standard atmos.

      real  :: sunAng =0.   !solar zenith angle
      real  :: sunAzm =0.   !solar azimuth angle
   END TYPE
   
   
             
CONTAINS !=================== MODULE CONTAINS ==========================


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE calcH1H2FromPress( scnGeom, prfl )
!-----------------------------------------------------------------------
      !USE Module_Scn_Geometry   ,ONLY: CLBLM_ScnGeom
      USE Module_ConstParam     ,ONLY: molNum
      USE Module_Scn_Profile    ,ONLY: CLBLM_Profile,&
                                       interpAltTempOnPressGrid
      IMPLICIT NONE

      type(CLBLM_ScnGeom)   ,intent(inout) :: scnGeom
      type(CLBLM_Profile)   ,intent(in)    :: prfl


      logical :: ZMDL_fromPM
      real    :: H1,H2
      integer :: LEN
      real    :: dummyAlt(1), dummyTemp(1)
      integer :: IMMAX

      
      !--- Check input values for Hobs, Hend, obsAng
      if ( scnGeom%Hobs <0 .or. scnGeom%Hend<0 .or. scnGeom%obsAng<0 .or. scnGeom%H_unit==0) then
         STOP '--- calcH1H2FromPress(): Invalid value for Hobs, Hend or obsAng.'          
      endif 

      
      IMMAX = prfl%toaLev
      ZMDL_fromPM = ( prfl%Z_unit == 210 ) !if altitude is computed from pressure profile.
      
                  
      !--- If the input H1 and H2 are pressure values, covert them to Km.
      !
      if ( scnGeom%H_unit <0 ) then

         !--- "interpAltTempOnPressGrid()" output array arguments only. To 
         ! compute scalar H1 and H2, use dummyAlt to receive the output array.
         !
         call interpAltTempOnPressGrid( dummyAlt, dummyTemp, & !intent(out) 
                                        [scnGeom%Hobs], & !H1 in pressure
                                        prfl%Z(1:IMMAX),&
                                        prfl%P(1:IMMAX),&
                                        prfl%T(1:IMMAX),&
                                        prfl%Q( molNum('H2O'),1:IMMAX ),&
                                        1, & !dimension of dummyAlt and dummyTemp
                                        IMMAX,&
                                        scnGeom%latitude, &
                                        scnGeom%earthRadius, &
                                        ZMDL_fromPM )
         H1 = dummyAlt(1)
         

         call interpAltTempOnPressGrid( dummyAlt, dummyTemp, & !intent(out)
                                        [scnGeom%Hend],& !H2 in pressure
                                        prfl%Z(1:IMMAX),&
                                        prfl%P(1:IMMAX),&
                                        prfl%T(1:IMMAX),&
                                        prfl%Q( molNum('H2O'),1:IMMAX ),&
                                        1, & !dimension of dummyAlt and dummyTemp
                                        IMMAX,&
                                        scnGeom%latitude, &
                                        scnGeom%earthRadius, &
                                        ZMDL_fromPM )
         H2 = dummyAlt(1)
         
         scnGeom%Hobs = H1
         scnGeom%Hend = H2
         scnGeom%H_unit = 210 !=210 means the H(Km) is computed from pressure
      endif
      
      
      IF (scnGeom%Hobs.LT.0.0 .OR. scnGeom%Hend.LT.0.0) THEN 
         !PRINT 946, H1,ZTMP(1) 
         !IF (NOPRNT.GE.0) WRITE (IPR,946) H1,ZTMP(1) 
         !STOP ' COMPUTED ALTITUDE VALUE OF H1 or H2 IS NEGATIVE' 
         STOP '--- calcH1H2FromPress(): Altitude value of Hobs or Hend is negative!' 
      ENDIF 

      
      if ( scnGeom%Hobs > scnGeom%HSPACE .and. &
           scnGeom%Hend > scnGeom%HSPACE ) then
         STOP '--- calcH1H2FromPress(): Both H1 and H2 are higher than TOA height (HSPACE)!'
      endif
      
   END SUBROUTINE


   
   
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   SUBROUTINE readSceneGeometry_tape5( scnGeom, tape5)
   !-----------------------------------------------------------------------
      !USE Module_Scn_Geometry    ,ONLY: CLBLM_ScnGeom
      USE Module_Scn_TAPE5       ,ONLY: CLBLM_TAPE5
      IMPLICIT NONE
       
      type(CLBLM_ScnGeom)  ,intent(inout) :: scnGeom !out
      type(CLBLM_TAPE5)    ,intent(in)    :: tape5

      real :: REF_LAT
      real :: RE
      real :: HSPACE
      
      if (tape5%rec_1_2%IATM ==1) then

         IF (tape5%rec_3_1%SREF_LAT .eq. '          ') THEN 
            REF_LAT = 45.0 
         ELSE 
            READ(tape5%rec_3_1%SREF_LAT,"(F10.3)") REF_LAT 
         ENDIF 

         RE = tape5%rec_3_1%RE
         if ( abs(RE) <epsilon(RE) ) then
            select case( tape5%rec_3_1%MODEL )
            case(1);      RE = 6378.39 
            case(4,5);    RE = 6356.91
            case default; RE = 6371.23
            end select
         endif
         
         HSPACE = tape5%rec_3_1%HSPACE
         IF (HSPACE.EQ.0.) HSPACE = 100. 
         
         !only H1,H2,ANGLE are allowed for ITYPE=2,3; warning infor.          
         if (abs(tape5%rec_3_2%RANGEF)>0. .or. abs(tape5%rec_3_2%BETAF)>0.) then
            print*,'--- readSceneGeometry_tape5(): CLBLM uses H1, H2, ANGLE and LEN only'
            STOP
         endif
         
         scnGeom%ITYPE       = tape5%rec_3_1%ITYPE
         scnGeom%Hobs        = tape5%rec_3_2%H1F
         scnGeom%Hend        = tape5%rec_3_2%H2F
         scnGeom%obsAng      = tape5%rec_3_2%ANGLEF
         scnGeom%LEN         = tape5%rec_3_2%LENF
         scnGeom%earthRadius = RE
         scnGeom%latitude    = REF_LAT
         scnGeom%HSPACE      = HSPACE

         if (tape5%rec_3_1%IBMAX_B <0) then; scnGeom%H_unit = -10 !in (mb)
         else;                               scnGeom%H_unit = 10 !in (Km)
         endif

      elseif (tape5%rec_1_2%IATM ==0) then
         !user provided path data, no need to get SceneGeometry.
      endif
      
   END SUBROUTINE

   
      
END MODULE
   