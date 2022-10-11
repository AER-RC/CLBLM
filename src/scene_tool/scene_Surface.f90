!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_Scn_Surface
   USE Module_ConstParam, ONLY: r8=>kind_r8
  
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_Surface, &
             readEmissivityFunc, &
             readReflectivityFunc, &
             readSurface_tape5
             
 
 
   integer, PARAMETER :: NMAXCO=4040
   TYPE CLBLM_EmissivityTable
      logical  :: tableIsReady = .FALSE.
      real(r8) :: V1
      real(r8) :: V2
      real     :: DV
      integer  :: NLIM
      real     :: ZEMIS(NMAXCO) 
   END TYPE

   
   !integer, PARAMETER :: NMAXCO=4040
   TYPE CLBLM_ReflectivityTable
      logical   :: tableIsReady = .FALSE.
      real(r8)  :: V1
      real(r8)  :: V2
      real      :: DV
      integer   :: NLIM
      real      :: ZRFLT(NMAXCO)       
   END TYPE
   

   TYPE :: CLBLM_Surface
      real                          :: Ts =-999.           !surface temperature
      real                          :: emisCoef(3)=0.     !frequency dependent boundary emissivity coefficients 
      real                          :: reflCoef(3)=0.     !frequency dependent boundary reflectivity coefficients 
      type(CLBLM_EmissivityTable)   :: emisTbl
      type(CLBLM_ReflectivityTable) :: reflTbl
      character(1)                  :: surf_refl='s'       !if='s', Specular reflection, default is 's'
                                                           !if='l', Lambertian surface reflection
                                                           !if='0', no surface reflection
   END TYPE
   
      
   !--- Surface emissivity and reflectivity files
   ! This is temporary arrangement. surfEmisFile and surfReflFile should
   ! be input from the scene file.
   character(256) :: surfEmisFile = 'EMISSIVITY'
   character(256) :: surfReflFile = 'REFLECTIVITY'
   

   
!   INTERFACE surfEmisRefl
!      module procedure surfEmisRefl_spect
!      module procedure surfEmisRefl_array   
!   END INTERFACE
   
   
CONTAINS !====================== MODULE CONTAINS =======================   


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Surface_init(this, nPts)
   !-------------------------------------------------------------------
      type(CLBLM_Surface) ,intent(inout) :: this
      integer             ,intent(in)    :: nPts
      
      STOP '--- CLBLM_Surface_init(): not yet implemented.'
      !if (allocated(this%???)   deallocate(this%???)
      !   
      !allocate( this%???(nPts) )   ;this%???(:)   = 0.
      !
      !this%nPts  = nPts
      
   END SUBROUTINE

   
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Surface_final(this)
   !-------------------------------------------------------------------
      type(CLBLM_Surface),  intent(inout) :: this

      STOP '--- CLBLM_Surface_final(): not yet implemented.'
      !if (allocated(this%???))   deallocate(this%???)
      !
      !this%nPts  = 0
      
   END SUBROUTINE


   !--------------------------------------------------------------------
   ! Reads in emission function values directly from file "EMISSIVITY" 
   !--------------------------------------------------------------------
   SUBROUTINE readEmissivityFunc( surf, surfEmisFile ) 
   !--------------------------------------------------------------------
      USE Module_Utility ,ONLY: getlun
      IMPLICIT NONE
      
      type(CLBLM_Surface) ,intent(inout) :: surf
      character(*)        ,intent(in)    :: surfEmisFile   
      
      integer  :: ICOEF,iostat
      real(r8) :: V1EMIS,V2EMIS
      real     :: DVEMIS
      integer  :: NGNU,NLIMEM
      real     :: ZEMIS(NMAXCO)
      

      ICOEF = getlun()
      OPEN (UNIT=ICOEF,FILE=trim(surfEmisFile), STATUS='OLD',IOSTAT=iostat)                                                     
      if ( iostat .gt. 0) &
         stop "FILE 'EMISSIVITY' FOR PATH BOUNDARY NOT FOUND"      
      
      !--- Read header information                                           
      READ (ICOEF,900) V1EMIS,V2EMIS,DVEMIS,NLIMEM 
  900 FORMAT (3E10.3,5X,I5) 

      if (nlimem.gt.nmaxco) then 
         print *, '*********************************************' 
         print *, ' Number of points on EMISSIVITY file > nmaxco' 
         print *, ' Also, check the REFLECTIVITY file' 
         stop 
      endif 

      !--- Read in emissivity values 
      DO NGNU = 1,NLIMEM 
         READ (ICOEF,*) ZEMIS(NGNU) 
      ENDDO
      
      surf%emisTbl%V1              = V1EMIS
      surf%emisTbl%V2              = V2EMIS
      surf%emisTbl%DV              = DVEMIS
      surf%emisTbl%NLIM            = NLIMEM
      surf%emisTbl%ZEMIS(1:NLIMEM) = ZEMIS(1:NLIMEM)
      surf%emisTbl%tableIsReady    = .TRUE.

      CLOSE(ICOEF)
      
   END SUBROUTINE
      
      
      
   !--------------------------------------------------------------------
   ! Reads in reflection function values directly from file "REFLECTIVI
   !--------------------------------------------------------------------
   SUBROUTINE readReflectivityFunc( surf, surfReflFile ) 
   !--------------------------------------------------------------------
      USE Module_Utility ,ONLY: getlun
      IMPLICIT NONE
      
      type(CLBLM_Surface)  ,intent(inout) :: surf
      character(*)         ,intent(in)    :: surfReflFile
      
      integer  :: ICOEF,iostat
      real(r8) :: V1RFLT,V2RFLT
      real     :: DVRFLT
      integer  :: NGNU,NLIMRF
      real     :: ZRFLT(NMAXCO)

      ICOEF = getlun()
      OPEN (UNIT=ICOEF,FILE=trim(surfReflFile), STATUS='OLD',IOSTAT=iostat)
      if ( iostat .gt. 0) &
         stop "FILE 'REFLECTIVITY' FOR PATH BOUNDARY NOT FOUND"                                         
      
      !--- Read header information 
      READ (ICOEF,900) V1RFLT,V2RFLT,DVRFLT,NLIMRF 
  900 FORMAT (3E10.3,5X,I5) 

      if (nlimrf.gt.nmaxco) then 
         print *, '*********************************************' 
         print *, ' Number of points on REFLECTIVITY file > nmaxco' 
         stop 
      endif 

      !--- Read in reflectivity values                                       
      DO NGNU = 1,NLIMRF 
         READ (ICOEF,*) ZRFLT(NGNU) 
      ENDDO
      
      surf%reflTbl%V1              = V1RFLT
      surf%reflTbl%V2              = V2RFLT
      surf%reflTbl%DV              = DVRFLT
      surf%reflTbl%NLIM            = NLIMRF
      surf%reflTbl%ZRFLT(1:NLIMRF) = ZRFLT(1:NLIMRF)
      surf%reflTbl%tableIsReady    = .TRUE.
      
      CLOSE(ICOEF)
      
   END SUBROUTINE

   
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   SUBROUTINE readSurface_tape5( surf, tape5 )
   !-----------------------------------------------------------------------
      USE Module_Utility         ,ONLY: lower
      USE Module_Scn_TAPE5       ,ONLY: CLBLM_TAPE5
      IMPLICIT NONE
            
      type(CLBLM_Surface)     ,intent(inout) :: surf !out
      type(CLBLM_TAPE5)       ,intent(in)    :: tape5
      
      character(1) :: sc

      
      !TYPE :: CLBLM_Surface
      !   real                         :: Ts =0.           !surface temperature
      !   !real                         :: Ps =0.          !surface presure
      !   !real                         :: Zs =0.          !surfave altitude.
      !   real                         :: emisCoef(3)=0.  !frequency dependent boundary emissivity coefficients 
      !   real                         :: reflCoef(3)=0.  !frequency dependent boundary reflectivity coefficients 
      !   type(CLBLM_EmissivityTable   :: emisTbl
      !   type(CLBLM_ReflectivityTable :: reflTbl
      !   character(1)                 :: SURF_REFL='s'     !='s'/'l'/'0'; specular or diffuse or no relection. No bidirectional consideration.
      !END TYPE
      !
      if (tape5%rec_1_4%BNDEMI(1) <0) then
         call readEmissivityFunc( surf, surfEmisFile )
      endif
      
      if (tape5%rec_1_4%BNDRFL(1) <0) then
         call readReflectivityFunc( surf, surfReflFile )
      endif

      surf%Ts             = tape5%rec_1_4%TMPBND
      surf%emisCoef(1:3) = tape5%rec_1_4%BNDEMI(1:3)
      surf%reflCoef(1:3) = tape5%rec_1_4%BNDRFL(1:3)
      
      !--- surf%surf_reflIf takes 'l' or 's'. default value is 's'.
      sc = lower( tape5%rec_1_4%surf_refl )
      if (sc/='l') sc = 's'
      surf%surf_refl = sc

!sfc'0'      !--- added in CLBLM
!sfc'0'      ! If reflCoef(:)==0, 's' or 'l' from surf will be ignored.
!sfc'0'      if (all( surf%reflCoef(:)==0. )) then
!sfc'0'         surf%surf_refl = '0'
!sfc'0'      endif
      
   END SUBROUTINE
   
   
END MODULE
