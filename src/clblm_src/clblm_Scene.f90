!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_Scene
   USE Module_ConstParam ,ONLY: FILLREAL,FILLINT

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_Profile, &
             CLBLM_Profile_init, &
             CLBLM_SceneGeom, &
             CLBLM_ObsGeom, &
             CLBLM_Surface,&
             CLBLM_Surface_init, &
             CLBLM_Scene, &
             CLBLM_Scene_init, &
             chopProfile, &
             readSceneFileHeader, readFileID, readNumScenes, readRTgrid, &
             inputSceneData, readGroupIDs, readSceneNo, &
             check_sceneData, &
             InterpPressToAlt, CMPALT, &
             surfThmEmis, &
             getSurfRefl



   TYPE :: CLBLM_Profile
      integer                    :: nLev    =FILLINT   !number of levels
      real          ,allocatable :: Z(:)               ![nLev],level heights
      real          ,allocatable :: P(:)               ![nLev],level pressures
      real          ,allocatable :: T(:)               ![nLev],level temperatures
      integer                    :: nMol    =FILLINT   !num of molecules, array dimmension.
      character(20) ,allocatable :: molID(:)           ![nMol], Molecule identifiers; Molecular name
      real          ,allocatable :: Q(:,:)             ![nMol,nLev], level concentration
      real                       :: surfAlt =FILLREAL  !First level altitude
   END TYPE


   TYPE :: CLBLM_ObsGeom
      real    :: obsAlt    = FILLREAL   !(Km)  observer level altitude
      real    :: obsPress  = FILLREAL   !(mb)  observer level pressure
      real    :: viewAng   = FILLREAL   !(deg) path zenith angle at H1
      !real    :: targetAlt = FILLREAL   !(Km)  target altitude; If<0, use Hobs and viewAng only (LBLRTM ITYPE=3A), If>0 user H1,H2 and viewAng (LBLRTM ITYPE=2A).
      !integer :: LEN       = 1          !When H1>H2>0 and ANGLE>90, If want a tangent path to space, set Len=1, otherwise set Len=0; ITYPE=3 (Hend<)) will ignore "LEN"
   END TYPE
   TYPE :: CLBLM_SceneGeom
      real                :: earthRadius = 6371.23  !
      real                :: latitude    = 45.      !(deg) This is needed when one want to fill up missing with AFGL standard atmos.
      real                :: sunAng      = 0.       !(deg) Solar zenith angle at TOA
      type(CLBLM_ObsGeom) :: obs                    !a structure containing observing geometry parameters
   END TYPE


   TYPE :: CLBLM_Surface
      real                 :: Tskin       = FILLREAL !surface temperature
      integer              :: ThmReflMode = 1        !(For thermal only) if=0 Lambertian; if=1 Specular
      integer              :: nsf         = FILLINT  !number of hinge points in surface optical property spectral grid
      real    ,allocatable :: surfEmRfGrid(:)     !wavenumber hinge points
      real    ,allocatable :: surfEm(:)           !frequency dependent boundary emissivity coefficients
      real    ,allocatable :: surfRefl(:)         !frequency dependent boundary reflectivity coefficients
   END TYPE


   TYPE :: CLBLM_Scene
      !integer               :: UID
      !integer               :: sceneID
      type(CLBLM_Profile)   :: prfl
      type(CLBLM_SceneGeom) :: geom
      type(CLBLM_Surface)   :: sfc              !a structure containing surface properties
   END TYPE



   INTERFACE surfThmEmis
      module procedure calcSurfThmEmis_spect
      module procedure calcSurfThmEmis_array
   END INTERFACE

   INTERFACE getSurfRefl
      module procedure getSurfRefl_spect
      module procedure getSurfRefl_array
   END INTERFACE

CONTAINS  !===================== Module Contains =======================


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE CLBLM_Profile_init( this, nMol, nLev)
!-----------------------------------------------------------------------
      type(CLBLM_Profile) ,intent(out) :: this
      integer             ,intent(in)  :: nLev
      integer             ,intent(in)  :: nMol

      this%nLev = nLev
      this%nMol = nMol

      if (allocated(this%Z))        deallocate(this%Z)
      if (allocated(this%P))        deallocate(this%P)
      if (allocated(this%T))        deallocate(this%T)
      if (allocated(this%molID))    deallocate(this%molID)
      if (allocated(this%Q))        deallocate(this%Q)

      allocate( this%Z(nLev)) ;this%Z(:) = 0.
      allocate( this%P(nLev)) ;this%P(:) = 0.
      allocate( this%T(nLev)) ;this%T(:) = 0.

      allocate( this%molID( nMol))   ;this%molID(:)  = ''
      allocate( this%Q( nMol,nLev))  ;this%Q(:,:)    = 0.

   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE CLBLM_Surface_init( this, emisPresent, reflPresent, numEmRfFreq )
!-----------------------------------------------------------------------
      type(CLBLM_Surface) ,intent(out) :: this
      logical             ,intent(in)  :: emisPresent   !If=.FALSE. emis not present;
      logical             ,intent(in)  :: reflPresent
      integer             ,intent(in)  :: numEmRfFreq


      if (emisPresent .or. reflPresent) then
         if (allocated(this%surfEmRfGrid)) deallocate(this%surfEmRfGrid)
         allocate( this%surfEmRfGrid( numEmRfFreq ) )

         this%nsf = numEmRfFreq
      endif

      if (emisPresent) then
         if (allocated(this%surfEm)) deallocate(this%surfEm)
         allocate(this%surfEm( numEmRfFreq ))
      endif

      if (reflPresent) then
         if (allocated(this%surfRefl)) deallocate(this%surfRefl)
         allocate(this%surfRefl( numEmRfFreq ))
      endif


   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE CLBLM_Scene_init( this, nMol, nLev, &
                                emisPresent, reflPresent, numEmRfFreq )
!-----------------------------------------------------------------------
      type(CLBLM_Scene) ,intent(out) :: this
      integer           ,intent(in)  :: nLev
      integer           ,intent(in)  :: nMol
      logical           ,intent(in)  :: emisPresent   !If=.FALSE. emis not present;
      logical           ,intent(in)  :: reflPresent
      integer           ,intent(in)  :: numEmRfFreq

      call CLBLM_Profile_init( this%prfl, nMol, nLev )

      call CLBLM_Surface_init( this%sfc, emisPresent, reflPresent,numEmRfFreq )

   END SUBROUTINE

   !--------------------------------------------------------------------
   ! Chop profile, extract profile between levLo and levHi
   !--------------------------------------------------------------------
   SUBROUTINE chopProfile( prfl, levLo, levHi, shortPrfl )
   !--------------------------------------------------------------------
      type(CLBLM_Profile) ,intent(in)   :: prfl
      integer             ,intent(in)   :: levLo, levHi
      type(CLBLM_Profile) ,intent(out)  :: shortPrfl

      call CLBLM_Profile_init( shortPrfl, prfl%nMol, levHi-levLo+1)

      shortPrfl%Z(  :)   = prfl%Z(  levLo:levHi)  ![nLev],level heights
      shortPrfl%P(  :)   = prfl%P(  levLo:levHi)  ![nLev],level pressures
      shortPrfl%T(  :)   = prfl%T(  levLo:levHi)  ![nLev],level temperatures
      shortPrfl%Q(:,:)   = prfl%Q(:,levLo:levHi)  ![nMol,nLev], level concentration
      shortPrfl%molID(:) = prfl%molID(:)          ![nMol], Molecule identifiers; Molecular name

   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readSceneFileHeader( ncID, UID, nScenes,sceneNums, &
                                   pRT,zRT,pRTsize,zRTsize, molID,nMol)
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4)               ,intent(in)  :: ncid
      integer                  ,intent(out) :: UID
      integer                  ,intent(out) :: nScenes
      integer     ,allocatable ,intent(out) :: sceneNums(:)
      real        ,allocatable ,intent(out) :: pRT(:)
      real        ,allocatable ,intent(out) :: zRT(:)
      integer                  ,intent(out) :: pRTsize
      integer                  ,intent(out) :: zRTsize
      character(*),allocatable ,intent(out) :: molID(:)
      integer                  ,intent(out) :: nMol

      integer                 :: i,j,k
      integer(4) ,allocatable :: groupIDs(:)

      call readFileID(    ncID, UID )
      call readNumScenes( ncID, nScenes )
      call readRTgrid(    ncID, pRT, zRT, pRTsize, zRTsize )
      call readMolID(     ncID, molID )
      call readNumMol(    ncID, nMol )

      call readGroupIDs(  ncID, groupIDs )  !one group per scene
      allocate( sceneNums( size(groupIDs) ) )
      do i=1,size(groupIDs)
         call readSceneNo( groupIDs(i), sceneNums(i) )
      enddo

   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readFileID( ncID, UID )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4), intent(in)   :: ncid
      integer,    intent(out)  :: UID

      call check( nf90_get_att( ncid, nf90_global, "UID", UID ) )

   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readNumScenes( ncID, nScn )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4), intent(in)   :: ncid
      integer,    intent(out)  :: nScn

      call check( nf90_get_att( ncid, nf90_global, "numScenes", nScn) )

   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readRTgrid( ncID, pRT, zRT, nRTlev, zRTsize )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check, attIsPresent
      IMPLICIT NONE

      integer(4)         ,intent(in)  :: ncid
      real  ,allocatable ,intent(out) :: pRT(:)
      real  ,allocatable ,intent(out) :: zRT(:)
      integer            ,intent(out) :: nRTlev
      integer            ,intent(out) :: zRTsize

      if ( attIsPresent(ncID, "RT_grid_lev_in_mb") ) then
         call check( nf90_get_att( ncid, nf90_global, "num_RT_grid_lev", nRTlev) )
         allocate( pRT(nRTlev) )
         call check( nf90_get_att( ncid, nf90_global, "RT_grid_lev_in_mb", pRT ) )
      else
         allocate( pRT(0))
         nRTlev = 0
      endif

      if ( attIsPresent(ncID, "RT_grid_lev_in_km") ) then
         call check( nf90_get_att( ncid, nf90_global, "num_RT_grid_lev", zRTsize) )
         allocate( zRT(zRTsize) )
         call check( nf90_get_att( ncid, nf90_global, "RT_grid_lev_in_km", zRT ) )
      else
         zRTsize = 0
      endif

      !integer :: nLev
      !call CLBLM_RTGrid_init( RTgrid, nLev )
      !call check( nf90_get_att( ncid, nf90_global, "RT_grid_lev_in_mb", RTgrid%P ) )
      !call check( nf90_get_att( ncid, nf90_global, "RT_grid_lev_in_km", RTgrid%Z ) )

   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readNumMol( ncID, nMol )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4) ,intent(in)  :: ncid
      integer    ,intent(out) :: nMol

      call check( nf90_get_att( ncid, nf90_global, "numMol", nMol) )

   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readMolID( ncID, molID, nMol )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4)                ,intent(in)  :: ncid
      character(*) ,allocatable ,intent(out) :: molID(:)
      integer         ,optional ,intent(out) :: nMol

      integer :: i,j,k,im,nm
      character :: molIDStr*1000


      call check( nf90_get_att( ncid, nf90_global, "numMol", nm) )
      allocate(molID( nm ))

      call check( nf90_get_att( ncid, nf90_global, "molecules", molIDStr ) )

      !--- Parsing the molIDStr to split the IDs into molID array.
      do im=1,nm-1
         molIDStr = adjustl(molIDStr)
         i=scan( molIDStr,',')
         molID(im) = molIDStr(1:i-1)
         molIDStr(1:i)=''
      enddo
      molID(nm) = trim(adjustl(molIDStr)) !the last moleclue

      if (present(nMol)) nMol = nm

   END SUBROUTINE


   !--------------------------------------------------------------------
   ! get the number of groups and an array of their ids. One group per scene
   !--------------------------------------------------------------------
   SUBROUTINE readGroupIDs( ncID, groupIDs )
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4)              ,intent(in)   :: ncid
      integer(4) ,allocatable ,intent(out)  :: groupIDs(:)

      integer(4) :: nScn, numgrps

      call check( nf90_get_att( ncid, nf90_global, "numScenes", nScn) )
      allocate( groupIDs(nScn) )
      call check( nf90_inq_grps(ncid, numgrps, groupIDs) )

   END SUBROUTINE

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE readSceneNo(grpID, sceneNo )
      USE NetCDF
      USE Module_FileIO ,ONLY: check
      IMPLICIT NONE

      integer(4), intent(in)   :: grpid
      integer,    intent(out)  :: sceneNo

      call check( nf90_get_att( grpID, nf90_global, "sceneNumber", sceneNo) )
   END SUBROUTINE
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! read scene from netCDF scene file
   !--------------------------------------------------------------------
   SUBROUTINE inputSceneData( ncID, sceneNo, scene)
   !--------------------------------------------------------------------
      USE NetCDF
      USE Module_ConstParam ,ONLY: isISOName, FILLREAL
      USE Module_FileIO     ,ONLY: openNetCDFFile, &
                                   varIsPresent, &
                                   dimIsPresent, &
                                   readArray_char, &
                                   readScalar_real, &
                                   readScalar_int, &
                                   readArray_real_1D, &
                                   readArray_real_2D, &
                                   closeNetCDFFile, &
                                   check
      IMPLICIT NONE

      integer(4)          ,intent(in)  :: ncID
      integer             ,intent(in)  :: sceneNo
      type(CLBLM_Scene)   ,intent(out) :: scene

      character(*) ,parameter :: routineName='inputSceneData'

      integer    :: i,j,k, scnNo, status
      integer(4) :: grpID, molDimId, levDimId, EmRfDimID
      integer(4) :: nLev, tempNumMol, EmRfDim
      logical    :: hasISO
      integer    :: sfcPropMode
      real          ,allocatable :: tempQ(:,:), tempQT(:,:)
      character(20) ,allocatable :: tempMolID(:)
      integer(4)    ,allocatable :: groupIds(:)

      !--- Find the groupID corresponding to the sceneNo.
      grpID = -1
      call readGroupIDs( ncID, groupIDs )
      do i=1,size(groupIDs)
         grpID = groupIDs(i)
         call readSceneNo( grpID, scnNo )
         if (scnNo==sceneNo) exit
      enddo
      if (grpID<0) STOP '--- '//routineName//'(): No scene data exist corresponding to the input sceneNo.'

      !--- Read group attribute
      call readMolID( ncID, tempMolID )
      !///to do: check if this scene contains the same moleclues as defined in the header


      !--- Inquire dimensions
      !
      call check( nf90_inq_dimid( grpId, "numPrflLev",levDimId ) )
      call check( nf90_inquire_dimension( grpId, levDimId ,len=nLev ) )
      call check( nf90_inq_dimid( grpId, "numMol", molDimId ) )
      call check( nf90_inquire_dimension( grpId, molDimId, len=tempNumMol ) )
      call check( nf90_inq_dimid( grpId, "numEmisReflPts", EmRfDimID ) )
      call check( nf90_inquire_dimension( grpId, EmRfDimID, len=EmRfDim) )


      call CLBLM_Scene_init( scene, int(tempNumMol), int(nLev), & !int() to default integer type
                             varIsPresent( grpId, "sfcEmis" ), &
                             varIsPresent( grpId, "sfcRefl" ), &
                             int(EmRfDim) )

      allocate(tempQ( nLev,tempNumMol ))
      allocate(tempQT( tempNumMol,nLev ))


      !--- Read variables
      !
      call readArray_real_1D( grpId, "pressure",     scene%prfl%P )
      call readArray_real_1D( grpId, "temperature",  scene%prfl%T )
      call readArray_real_2D( grpId, "molDensities", tempQ ); !Isotopologues are input as fractions
      call readScalar_real(   grpId, "viewAngle",    scene%geom%obs%viewAng )

      if (varIsPresent( grpId, "altitude") )         call readArray_real_1D( grpId, "altitude",          scene%prfl%Z )
      if (varIsPresent( grpId, "sfcAltitude" ))      call readScalar_real(   grpId, "sfcAltitude",       scene%prfl%surfAlt )
      if (varIsPresent( grpId, "earthRadius" ))      call readScalar_real(   grpId, "earthRadius",       scene%geom%earthRadius )
      if (varIsPresent( grpId, "latitude" ))         call readScalar_real(   grpId, "latitude",          scene%geom%latitude )
      if (varIsPresent( grpId, "sfcSkinTemp" ))      call readScalar_real(   grpId, "sfcSkinTemp",       scene%sfc%Tskin)
      if (varIsPresent( grpId, "sfcPropInputMode" )) call readScalar_int(    grpId, "sfcPropInputMode",  sfcPropMode)
      if (varIsPresent( grpID, "sfcPropSpectrGrid" ))call readArray_real_1D( grpId, "sfcPropSpectrGrid", scene%sfc%surfEmRfGrid )
      if (varIsPresent( grpID, "sfcEmis" ))          call readArray_real_1D( grpId, "sfcEmis",           scene%sfc%surfEm )
      if (varIsPresent( grpID, "sfcRefl" ))          call readArray_real_1D( grpId, "sfcRefl",           scene%sfc%surfRefl )

      if (varIsPresent( grpId, "obsAltitude" ))      call readScalar_real(   grpId, "obsAltitude", scene%geom%obs%obsAlt )
      if (varIsPresent( grpId, "obsPressure" ))      call readScalar_real(   grpId, "obsPressure", scene%geom%obs%obsPress )


      !--- Set surface thermal reflection mode
      if     (sfcPropMode==1) then; scene%sfc%ThmReflMode = 1;  !Specular thermal reflection
      elseif (sfcPropMode==2) then; scene%sfc%ThmReflMode = 0;  !Lambertian thermal reflection
      endif

      print *, 'sfcPropMode', sfcPropMode
      print *, 'ThmReflMode', scene%sfc%ThmReflMode

      !--- Isotopologue are input as fractions. Convert to number density.
      ! Whenever a molecule has isotopologule fractions input, fractions for
      ! all of isotopologules for that molecule will be adjusted to keep the
      ! total density unchanged. Therefore the output prfoile will
      ! contain densities for all isotopologules for that molecule. Total
      ! number of molecule will be changed.
      !
      !-Check if there is isotopologule included
      hasISO = .FALSE.
      do i=1,size(tempMolID)
         if (isISOName( tempMolID(i) )) then
            hasISO = .TRUE.
            exit
         endif
      enddo

      tempQT = transpose(tempQ) !->[nMol,nLev]
      if (hasISO) then
         call intervalSub_ISOTPL_fraction2density( tempQT, tempMolID )
      endif
      scene%prfl%Q     = tempQT
      scene%prfl%molID = tempMolID
      scene%prfl%nMol  = size(tempMolID)

   CONTAINS !----------------- Internal subroutine ---------------------


   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
      subroutine  intervalSub_ISOTPL_fraction2density( prflQ, molID )
   !-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: MXMOL, MXISOTPL, MXLAY, &
                                     ISOTPL_ABD, MaxISOTPL_abd, &
                                     molIndex, getMolNum=>molnum, &
                                     isISOName, &
                                     trimISONum, &
                                     isoName2molNum, &
                                     isoName2isoNum, &
                                     molNum_isoNum_to_isoName
      USE Module_Config       ,ONLY: IPR
      IMPLICIT NONE

      real         ,allocatable ,intent(inout) :: prflQ(:,:) ![nMol,nLay]
      character(*) ,allocatable ,intent(inout) :: molID(:)


      character(*) ,parameter :: routineName = 'inputSceneData'

      integer               :: i,j,k,l,m, im,imm, nMol,nLay, pMol, iMol, molNo, isoNo, thisMolNo
      integer               :: NISOTPL, NLAYIS
!      integer               :: ISOTPL_HCODE(MXMOL*MXISOTPL)
      integer               :: MOLCOUNT
      integer               :: MOLNUM(MXMOL*MXISOTPL)
      integer               :: INPTYP(MXMOL*MXISOTPL)
      integer               :: ISOTPLNUM(MXMOL*MXISOTPL)
      integer               :: ISOTPL_FLAG(MXMOL,MXISOTPL)
 !     integer               :: ISOTPL_MAIN_FLAG(MXMOL)
      integer               :: MOLLST(MXMOL)
      character(2)          :: MOLSTR
      real                  :: ISOTPL_ABD_SUM(MXMOL), ISOTPL_ABD_SUBSUM(MXMOL)
      real                  :: ISOTPL_AMNT_SUM(MXMOL)
      real                  :: ISOTPL_AMNT(MXMOL,MXISOTPL,MXLAY)
!      real                  :: WKL(path%lnMol%nMol,path%nLay)
      real                  :: WKL_ISOTPL(MXMOL,MXISOTPL,MXLAY)
      real                  :: ADJ_FRAC
      integer               :: ISOCOUNT!, IDXCNT
!      integer               :: IDXM(MXMOL*MXISOTPL),IDXJ(MXMOL*MXISOTPL)
      real                  :: QKL(MXMOL,MXLAY)
      character(len(molID)) :: molNm
      integer               :: newNumMol,nNonIso
      integer               ,allocatable :: isoMark(:)
      real                  ,allocatable :: newPrflQ(:,:)
      character(len(molID)) ,allocatable :: newMolID(:)



      nMol = size(molID)
      nLay = size(prflQ,2) ![nMol,nLay]

      NLAYIS  = nLay

      allocate( isoMark(nMol) )
      isoMark(:) = 0

      NISOTPL = 0
      do im = 1,nMol

         if (isISOName(molID(im))) then

            molNm = trimISONum( molID(im) )
            molNo = isoName2molNum( molID(im) )
            isoNo = isoName2isoNum( molID(im) )

            ISOTPL_AMNT( molNo,isoNo,1:nLay ) = prflQ(im,1:nLay)

            pMol = molIndex( molNm, molID ) !find the location of the parent molecule in molID array
            if (pMol<0) STOP '---'//routineName//'(): Concentration of parent molecule must be present when input isotopologue fractions.'
            QKL( molNo,1:nLay ) = prflQ( pMol,1:nLay )

            NISOTPL = NISOTPL+1
            MOLNUM(    NISOTPL ) = molNo
            ISOTPLNUM( NISOTPL ) = isoNo

            !--- Mark the isotopologues and their parent molecule in molID list
            ! isoMark(i)==0 means there is no isotopologules for molecule i exist in the list.
            ! isoMark(i)/=0 means the molecule i is either a parent molecule or a isotopologule.
            ! e.g. molID = [ H2O,CO2,O3,CO2_727,H2O_181,H2O_161]
            !    isoMark = [  3,  1,  0,   1,      2       3 ]
            isoMark( pMol ) = NISOTPL
            isoMark( im   ) = NISOTPL
         endif
      enddo

      ISOTPL_FLAG(:,:) = 0


      MOLCOUNT = 1
      DO I = 1,NISOTPL

         ISOTPL_FLAG( MOLNUM(I),ISOTPLNUM(I) ) = 1

         IF (I.EQ.1) THEN
            MOLLST(MOLCOUNT) = MOLNUM(I)
         ELSE
            IF (MOLNUM(I).NE.MOLNUM(I-1)) THEN !yma160317, what if ISOTPL_HCODE is not in monotonic order?
               MOLCOUNT = MOLCOUNT + 1
               IF (MOLCOUNT.GT.MXMOL) THEN
                  STOP 'MOLCOUNT GREATER THAN MXMOL IN PATH'
               ENDIF
               MOLLST(MOLCOUNT) = MOLNUM(I)
            ENDIF
         ENDIF
      ENDDO


      DO L = 1,NLAYIS

         !     Check for negative input, which is permitted to identify when amounts
         !     for isotopologues not directly specified as column amount are to be
         !     set to their original default values (i.e. WKL * HITRAN ratio).
         DO I = 1,NISOTPL
            IF (ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L).LT.0.) THEN

                ISOTPL_AMNT( MOLNUM(I),ISOTPLNUM(I),L ) = &
                QKL(MOLNUM(I),L) * ISOTPL_ABD( MOLNUM(I),ISOTPLNUM(I) )
            ENDIF
         ENDDO

         !     Check input for consistent type for each specified molecular species;
         DO I = 1,NISOTPL
            INPTYP(I) = 0
            IF (ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L).LE.1.) THEN
               INPTYP(I) = 1
            ENDIF
         ENDDO

         DO I = 1,NISOTPL-1
            DO J = I+1,NISOTPL
              IF (MOLNUM(I).EQ.MOLNUM(J).AND.INPTYP(I).NE.INPTYP(J)) THEN
                WRITE(IPR,'(I3,2E15.7)') L,ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L), &
                                 ISOTPL_AMNT(MOLNUM(J),ISOTPLNUM(J),L)
                WRITE(*,'(I3,2E15.7)')&
                               L,ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L), &
                               ISOTPL_AMNT(MOLNUM(J),ISOTPLNUM(J),L)
                WRITE(MOLSTR,'(I2)') MOLNUM(I)
                STOP 'ISOTPL_AMNT INPUT TYPES DIFFERENT'
             ENDIF
            ENDDO
         ENDDO
      ENDDO !DO L = 1,NLAYIS

      !     Get sum of isotopologue Hitran ratios for each molecule
      DO I=1,MXMOL
         ISOTPL_ABD_SUM(I) = 0.0
         DO J=1,MaxISOTPL_abd(I)
            ISOTPL_ABD_SUM(I) = ISOTPL_ABD_SUM(I) + ISOTPL_ABD(I,J)
         ENDDO
      ENDDO



      !----------------------------------------------------------
      !
      !         CONVERSION FROM FRACTION TO COLUMN AMOUNT
      !
      ! Isotopologue amounts are specified as fractions (in terms
      ! of moleulces, not mass) for that species.  Fractions of all
      ! isotopologues are adjusted to force the total density for
      ! that species to remain fixed, and all isotopologues are
      ! activated for that species. Adjusted fractions are used to
      ! define the density for each isotopologue from the density of
      ! the parent molecule density input earlier.
      !
      DO L = 1,NLAYIS

         DO I = 1,MOLCOUNT
            M = MOLLST(I)

            ! Generate sums needed for fractional input
            ISOTPL_ABD_SUBSUM(M) = ISOTPL_ABD_SUM(M)  !sum of un-selected fractions.
            ISOTPL_AMNT_SUM(M) = 0.0                  !sum of slected fractions.
            DO J = 1, MaxISOTPL_abd(M)
               IF (ISOTPL_FLAG(M,J).EQ.1 .AND. &
                   ISOTPL_AMNT(M,J,L).LE.1.) THEN

                  ISOTPL_ABD_SUBSUM(M) = ISOTPL_ABD_SUBSUM(M) - &
                                         ISOTPL_ABD(M,J)

                  ISOTPL_AMNT_SUM(M) = ISOTPL_AMNT_SUM(M) + &
                                       ISOTPL_AMNT(M,J,L)
               ENDIF
            ENDDO

            ! Test to check sum of specified fractions less than sum of Hitran ratios
            IF (ISOTPL_AMNT_SUM(M) .GT. ISOTPL_ABD_SUM(M)) THEN
               WRITE(IPR,'(2I3,2E15.7)') L, I, ISOTPL_AMNT_SUM(M), &
                             ISOTPL_ABD_SUM(M)
               WRITE(*,'(2I3,2E15.7)') L, I, ISOTPL_AMNT_SUM(M), &
                           ISOTPL_ABD_SUM(M)
               STOP 'ISOTPL_AMNT NOT PROPERLY SPECIFIED IN PATH; ' &
                   //'SUM EXCEEDS SUM OF HITRAN RATIOS'
            ENDIF


            DO J = 1, MaxISOTPL_abd(M)

               ! Isotopologue fraction input
               ! Total density for parent species remains constant
               ! Isotopologue fraction scaled by Hitran ratio and converted to density

               IF (ISOTPL_AMNT(M,J,L) >1.) STOP '--- '//routineName//'(): Isotopologue profile must be input as fractions.'

               ! Selected molecules/isotopologues; unmodified fractions
               IF (ISOTPL_FLAG(M,J).EQ.1) THEN

                  ADJ_FRAC = ISOTPL_AMNT(M,J,L)
               ENDIF

               ! Unselected molecules/isotopologues; modified fractions
               IF (ISOTPL_FLAG(M,J).EQ.0) THEN
                  ! Test to prevent divide by zero, though code is not likely to reach this
                  ! point in that condition
                  IF (ISOTPL_ABD_SUBSUM(M) .EQ. 0.0) THEN
                     WRITE(IPR,'(2I3,E15.7)') L, I, &
                           ISOTPL_ABD_SUBSUM(M)
                     WRITE(*,'(2I3,E15.7)') L, I, &
                           ISOTPL_ABD_SUBSUM(M)
                     STOP 'ISOTPL_ABD_SUBSUM IS ZERO IN PATH'
                  ENDIF

                  ADJ_FRAC = &
                      (ISOTPL_ABD_SUM(M) - ISOTPL_AMNT_SUM(M)) &
                    * (ISOTPL_ABD(M,J) / ISOTPL_ABD_SUBSUM(M))
               ENDIF

               ! Generate isotopologue density with fractions, from
               ! density for parent molecule input earlier
               WKL_ISOTPL(M,J,L) = QKL(M,L) * ADJ_FRAC

            ENDDO !iso
         ENDDO !mol
      ENDDO !layer

      ! Activate all isotopologues for molecules with selected isotopologues
      ! provided as fractional input
      DO I = 1,NISOTPL
         M = MOLNUM(I)
         IF (MAXVAL(ISOTPL_FLAG(M,1:MXISOTPL)).GT.0.AND.&
             INPTYP(I).EQ.1) THEN

            DO J = 1, MaxISOTPL_abd(M)
               ISOTPL_FLAG(M,J) = 1
            ENDDO
         ENDIF
      ENDDO


      ! Get total number of activated isotoplogues
      ISOCOUNT = SUM(ISOTPL_FLAG)


      newNumMol = nMol-NISOTPL + ISOCOUNT
      allocate( newPrflQ( newNumMol, nLay ))
      allocate( newMolID( newNumMol ))

      !--- Count the number of non-iso molecues (not iso and not their parent)
      nNonIso = 0
      do im = 1,nMol
         if (isoMark(im)==0) nNonIso=nNonIso+1
      enddo

      !--- Load the new prfile array
      ! Prfl contains densities for both isotopologules and their parents, although
      ! all isotopologules are activated. Continuum and RHOSLF calculations still
      ! use total amounts (amounts for parent molecules).
      !
      iMol = 0
!      iIso = nNonIso
      do im = 1,nMol

         if ( isoMark(im)>0 ) then !is a iso molecule or its parent molecule

            thisMolNo = getMolNum( trimISONum( molID(im) ))

            !--- Find and load parent molecule
            do imm = im,nMol
               if ( getMolNum(trimISONum(molID(imm))) == ThisMolNo ) then
                  if ( .not.isISOName( molID(imm) ) ) then
                     iMol = iMol+1
                     newPrflQ( iMol,1:nLay ) = prflQ(imm,1:nLay)
                     newMolID( iMol ) = adjustl( molID(imm) )
                  endif

                  !At this stage, parent molecule has been loaded, All isotopologules
                  !will be loaded in the next step. Marks are not needed anymore.
                  !Set all marks for this molecule including parent and isotopologules
                  !to a negative value.
                  isoMark(imm) = -99
                endif
            enddo


            !--- Load isotopologules, all of them.
            do isoNo = 1, MaxISOTPL_abd( thisMolNo )

               iMol = iMol+1
               newPrflQ( iMol,1:nLay ) = WKL_ISOTPL( thisMolNo,isoNo,1:nLay)
               newMolID( iMol ) = adjustl( molNum_isoNum_to_isoName( thisMolNo,isoNo ) )
            enddo

         elseif (isoMark(im)==0) then !non-iso molecule

            iMol = iMol+1
            newPrflQ( iMol,1:nLay ) = prflQ(im,1:nLay)
            newMolID( iMol ) = adjustl(molID(im))
         endif

      enddo !do im = 1,nMol

      call move_alloc( newPrflQ, prflQ )
      call move_alloc( newMolID, molID )

      end subroutine

   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE check_sceneData( scn, &
                               molID, numMol, &
                               pRT, pRTsize, zRt, zRTsize, RTgridFlag, &
                               ovrdObsAlt, ovrdViewAng, outCtrl )
!-----------------------------------------------------------------------
      USE Module_Utility    ,ONLY: upper
      USE Module_ConstParam ,ONLY: FILLREAL, molIndex
      USE Module_Config     ,ONLY: CLBLM_Output_Ctrl

      type(CLBLM_Scene)       ,intent(inout) :: scn
      character(*)            ,intent(in)    :: molID(:)
      integer                 ,intent(in)    :: numMol
      real                    ,intent(in)    :: pRT(:) !RT grid will be convert to altitude grid.
      integer                 ,intent(in)    :: pRTsize
      real       ,allocatable ,intent(inout) :: zRT(:)
      integer                 ,intent(inout) :: zRTsize
      integer                 ,intent(in)    :: RTgridFlag
      real                    ,intent(in)    :: ovrdObsAlt
      real                    ,intent(in)    :: ovrdViewAng
      type(CLBLM_Output_Ctrl) ,intent(in)    :: outctrl

      character(*) ,parameter :: routineName = 'check_sceneData'
      logical :: altIsPresent, ignoreInterpZ, ODonly
      integer :: il,im,nLev
      real    :: tempAlt(1)



      !--- Check if all scenes have the same molecules
      if (scn%prfl%nMol /= numMol) STOP '--- '//routineName//'(): All scenes in a file should have the same number of molecules.'
      do im=1,numMol
         if ( upper(trim(adjustl(scn%prfl%molID(im)))) /= upper(trim(adjustl(molID(im)))) ) then
            STOP '--- '//routineName//'(): All scenes in a file must consist of same molecules.'
         endif
      enddo

      !///to do: if (speciesBroad), co2 needs h2o+co2, o2 needs h2o, h2o needs co2

      !--- If to use generic P-grid or Z-grid from the header, hdrRTgrid must be exist
      if ( RTgridFlag==1 .and. pRTsize==0 .and. zRTsize==0) &
         STOP '--- '//routineName//'(): Choose to use the generic P-grid or Z-grid from header but the P-grid/Z-grid doesnt exist.'


      !--- Detect the presence of altitude profile before it is calcualted from pressure grid.
      altIsPresent = .not. all(scn%prfl%Z(2:) <=0.)


      !--- Check if altitude exist or not, if not, calcualte from pressure
      ignoreInterpZ = .FALSE.
      if ( all(scn%prfl%Z <=0.) ) then

         if (scn%prfl%surfAlt <=FILLREAL) &
            STOP '--- '//routineName//'(): Altitude (Km) for the first pressure level must be present when altitude profile is not present.'

         nLev = scn%prfl%nLev
         call CMPALT( nLev,                                                & !ILVL
                      scn%prfl%P(1:nLev),                                  & !PM in (mb)
                      scn%prfl%T(1:nLev),                                  & !TM in (K)
                      scn%prfl%Q( molIndex('H2O',scn%prfl%molID),1:nLev ), & !DENW
                      scn%prfl%surfAlt,                                    & !REF_Z in (Km)
                      scn%geom%latitude,                                   & !REF_LAT
                      scn%geom%earthRadius,                                & !RE
                      scn%prfl%Z(1:nLev) )                                  !ZMDL in (Km)

         ignoreInterpZ = .TRUE.
      endif

      !--- Check if the input profile altitudes are in a ascending order
      do il = 2,scn%prfl%nLev
         if (scn%prfl%Z(il) <= scn%prfl%Z(il-1)) then
            STOP '--- '//routineName//'(): INPUT ALTITUDES ARE NOT IN ASCENDING ORDER.'//' Program stopped.'
         endif
      enddo


      !--- Check the input of obsAlt and obsPress
      !
      if (altIsPresent .and. RTgridFlag/=1) then

         if (scn%geom%obs%obsAlt <=FILLREAL) &
            STOP '---'//routineName//'(): Observer height should be given in altitude (Km) when altitude profile is present.'

      else !If input altitudes are not supplied or if generic P-grid in file header is used.

         if ( scn%geom%obs%obsPress >0. ) then
            call InterpPressToAlt( tempAlt, &
                                   [scn%geom%obs%obsPress], &
                                   scn%prfl%Z, &
                                   scn%prfl%P, &
                                   scn%prfl%T, &
                                   scn%prfl%Q( molIndex('H2O',scn%prfl%molID),: ),&
                                   scn%geom%latitude,  &
                                   scn%geom%earthRadius, &
                                   ignoreInterpZ )
            scn%geom%obs%obsAlt= tempAlt(1)

         elseif ( scn%geom%obs%obsAlt < scn%prfl%Z( scn%prfl%nLev )) then
           ! STOP '---'//routineName//'(): When using pressure RT grid, obsAltitude is only allowed for higher than TOA cases.'
         endif
      endif


      !--- If override values exist, override the values from the scene file.
      if (.not. ovrdObsAlt<=FILLREAL) scn%geom%obs%obsAlt =ovrdObsAlt
      if (.not. ovrdViewAng<=FILLREAL) scn%geom%obs%viewAng =ovrdViewAng


      !--- If OD only calculation, viewing angle reset to vertical direction.'
      !
      ODonly = (outCtrl%Rad==0 .and. outCtrl%Tx==0 .and. outCtrl%OD/=0 )
      if (ODonly) then
         if (scn%geom%obs%viewAng <=90.) scn%geom%obs%viewAng = 0.
         if (scn%geom%obs%viewAng  >90.) scn%geom%obs%viewAng = 180.
         print*, '--- '//routineName//'(): OD-only calculation, viewing angle reset to vertical direction.'
      endif


      !--- Hobs and Hend has to be >0.
      !IF ( scn%geom%obs%obsAlt.LT.0.0 ) then !.OR. scn%geom%obs%targetAlt.LT.0.0) THEN
      !   STOP '--- '//routineName//'(): Hobs or Hend is negative.'//' Program stopped.'
      !ENDIF


      !--- Check Hobs, if it is given higher then HSPACE, change H1 to HSPACE and
      !   view angle to be the value at HSPPCE
      if ( scn%geom%obs%obsAlt > scn%prfl%Z( scn%prfl%nLev ) ) then
         scn%geom%obs%viewAng = viewAngAtTOA( scn%geom%obs%viewAng, &
                                              scn%geom%obs%obsAlt, &
                                              scn%prfl%Z( scn%prfl%nLev ), &
                                              scn%geom%earthRadius )
         scn%geom%obs%obsAlt = scn%prfl%Z( scn%prfl%nLev)
         print*, '--- '//routineName//'(): Hobs higher than top of profile, Hobs and obsAng adjusted to the values at top of profile.'
      endif


      print *,'zRT', zRT
      print *,'pRTsize', pRTsize
      !--- Calculate generic RT grid in (Km)
      if (pRTsize>0) then
         print *,'Calculate generic RT grid in (Km)'
         zRTsize = pRTsize
         allocate( zRT(zRTsize) )
         call InterpPressToAlt( zRT(1:zRTsize), & !out
                                pRT(1:pRTsize), &
                                scn%prfl%Z,&
                                scn%prfl%P,&
                                scn%prfl%T,&
                                scn%prfl%Q( molIndex('H2O',scn%prfl%molID),: ),&
                                scn%geom%latitude, &
                                scn%geom%earthRadius, &
                                ignoreInterpZ ) !ignoreInterpZ=.TRUE., use hydrostatic result only.
      endif

   END SUBROUTINE

!-----------------------------------------------------------------------
! * Given view angle at observer height, calculate the view angle at TOA height.
!   Mainly used for satellite application when user inputs view angle at satellite
!   height and wants view angle at TOA level.
! * Input viewAngAtObs can be either zenith or nadir angle. If input zenith
!   angle, the output zenith angle too, and vise versa.
! * Hobs can lower or higher than Htoa
!-----------------------------------------------------------------------
   real FUNCTION viewAngAtTOA( viewAngAtObs, Hobs, Htoa, RE )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: pi
      IMPLICIT NONE

      real ,intent(in) :: viewAngAtObs !(Degree)
      real ,intent(in) :: Hobs, Htoa   !Altitudes in (Km)
      real ,intent(in) :: RE           !The earth radius (Km)

      real :: Robs, Rtoa
      real :: sinTOA, angTOA


      Robs = RE + Hobs
      Rtoa = RE + Htoa

      sinTOA = Robs/Rtoa*sin( viewAngAtObs *pi/180. )
      angTOA = asin( sinTOA ) *180./pi

      if (viewAngAtObs>90.)  angTOA = 180.-angTOA

      viewAngAtTOA = angTOA ! in (Degree)

   END FUNCTION


!-----------------------------------------------------------------------
! * TO ENSURE THAT CALCULATED/INPUT ZMDL'S WILL MATCH CALCULATED USER-LEVE
!   ALTITUDES, A COMBINATION OF INTERPOLATION AND HYDROSTATICS ARE USED.
!   alt = A * F1(P) + (1 - A) * F2(P), WHERE
!   F1(P) = INTERPOLATION IN LN(P), F2(P) = HYDROSTATIC CALCULATION
! * v12.7: Ignore exponential interpolation term (zint) if user has input
!   a profile specified on pressure levels
! * pres and alt are array arguments.
! * It returns temperatures at pressure grid as well.
! * "InterpPressToAlt()" differ from "CMPALT()" in that the latter
!   assume the available of Temperature and Watervapor at pressure grid.
!   "InterpPressToAlt()" obtains the Temperature and Watervapor by
!   interpolation from TM and DENW.
!-----------------------------------------------------------------------
   SUBROUTINE InterpPressToAlt( alt, pres, &
                                ZMDL,PM,TM,DENW, REF_LAT, RE, ignoreInterpZ )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: alt(:)        ! out
      real    ,intent(in)  :: pres(:)       ! boundary level in pressure
      real    ,intent(in)  :: ZMDL(:)       ! model atm. Altitude
      real    ,intent(in)  :: PM(:)         ! model atm. Pressure
      real    ,intent(in)  :: TM(:)         ! model atm. Temperature
      real    ,intent(in)  :: DENW(:)       ! model atm. water vapour density.
      real    ,intent(in)  :: REF_LAT       ! reference latitude
      real    ,intent(in)  :: RE            ! earth radius
      logical ,intent(in)  :: ignoreInterpZ ! If ZMDL was computed from P profile, Altitude will be calculated from hydrostatic equation only.


      integer           :: IBMAX  ! number of boundary levels
      integer           :: IMMAX  ! number of model atm. levels
      real ,allocatable :: Temp(:)
      integer           :: ISTART, IP, LIP
      real              :: ZTMP(2),PTMP(2),TTMP(2),WVTMP(2)
      real              :: HIP,ZINT
      real              :: TIP
      real              :: WVIP
      real              :: RATP,A


      IBMAX = size(pres)
      IMMAX = size(ZMDL)
      allocate( Temp(IBMAX) )

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
      USE Module_Config     ,ONLY: IPR

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

   !--------------------------------------------------------------------
   ! Returns surface emitted radiance and emissivity. This is a wrapper for calcSurfRadEmis_array.
   !--------------------------------------------------------------------
   SUBROUTINE calcSurfThmEmis_spect( surf, V1,DV,NLIM, sfcRad, sfcEmis )
   !--------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init

      type(CLBLM_Surface)           ,intent(in)  :: surf
      real(r8)                      ,intent(in)  :: V1
      real                          ,intent(in)  :: DV
      integer                       ,intent(in)  :: NLIM
      type(CLBLM_Spectrum) ,optional,intent(out) :: sfcRad
      type(CLBLM_Spectrum) ,optional,intent(out) :: sfcEmis

      real(r8) :: dummyV2
      real,allocatable :: sRad(:)
      real,allocatable :: sEmis(:)


      !---
      ! * Please be noted that sfcRad%indV1 = 1
      if ( present(sfcRad) .and. present(sfcEmis) ) then
      print *,'returns sfc rad and emis from calcSurfThmEmis_spect in clblm_scene'
         allocate(sRad(NLIM))
         allocate(sEmis(NLIM))

         call calcSurfThmEmis_array( surf, V1,dummyV2,DV,NLIM, sRad, sEmis )

         call CLBLM_Spectrum_init( sfcRad,  V1,DV,NLIM )
         call CLBLM_Spectrum_init( sfcEmis, V1,DV,NLIM )
         call move_alloc( sRad,   sfcRad%spect )
         call move_alloc( sEmis,  sfcEmis%spect )

      elseif (present(sfcRad)) then
         allocate(sRad(NLIM))
         call calcSurfThmEmis_array( surf, V1,dummyV2,DV,NLIM, sRad )
         call CLBLM_Spectrum_init( sfcRad, V1,DV,NLIM )
         call move_alloc( sRad, sfcRad%spect )

      elseif (present(sfcEmis)) then

         allocate(sEmis(NLIM))
         call calcSurfThmEmis_array( surf, V1,dummyV2,DV,NLIM, sEmis=sEmis )
         call CLBLM_Spectrum_init( sfcEmis, V1,DV,NLIM )
         call move_alloc( sEmis,  sfcEmis%spect )
      endif

   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE calcSurfThmEmis_array( surf, V1PO,V2PO,DVPO,NLIMO, sRad, sEmis )
   !--------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: RADCN2, r8=>kind_r8
      USE Module_EMLAY       ,ONLY: planck

      type(CLBLM_Surface) ,intent(in)  :: surf
      real(r8)            ,intent(in)  :: V1PO
      real(r8)            ,intent(in)  :: V2PO !used to distinguish between calcSurfThmEmis_array and calcSurfThmEmis_spect
      real                ,intent(in)  :: DVPO
      integer             ,intent(in)  :: NLIMO
      real      ,OPTIONAL ,intent(out) :: sRad(:)  !(NLIMO)
      real      ,OPTIONAL ,intent(out) :: sEmis(:) !(NLIMO)


      integer  :: I,J,K
      real     :: XKTBND, BB
      real     :: EMISIV, EMDEL, EMLAST
      real(r8) :: V, VI, VIDV
      integer  :: NLIM1, NLIM2


      XKTBND  = surf%Tskin / RADCN2

      VI = V1PO-DVPO
      VIDV = VI
!      VIDVEM = VI

      EMLAST = -1.
!      EMDUM  = 0.

      NLIM1 = 0
      NLIM2 = 0


!      EMISIV = EMISFN( VI,DVPO,VIDVEM,EMDEL,EMDUM, surf )
!      IEMBB = 1

!      if ( present(sRefl) ) then
!         VIDVRF = VI
!         RFLAST = -1.
!         RFDUM = 0.
!         REFLCT = REFLFN( VI,DVPO,VIDVRF,RFDEL,RFDUM, surf )
!         if ( VIDVRF.LE.VIDVEM)  then
!            IEMBB = 2
!         endif
!      endif


   60 NLIM1 = NLIM2+1

         VI = V1PO+ REAL(NLIM1-1)*DVPO
!         IF (IEMBB.EQ.1) THEN
            EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST, surf)
!            if (present(sRefl)) then
!               VIDVRF = -VIDV
!               REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFLAST, surf)
!            endif
!         ELSE
!            REFLCT = REFLFN(VI,DVPO,VIDV,RFDEL,RFLAST, surf)
!            VIDVEM = -VIDV
!            EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST, surf)
!         ENDIF

         IF (VIDV.GE.9.E+4) THEN
            NLIM2 = NLIMO+1
         ELSE
            NLIM2 = (VIDV-V1PO)/DVPO+1.001
         ENDIF
         NLIM2 = MIN(NLIM2,NLIMO)

         if ( present(sRad) .and. present(sEmis) ) then

            DO J = NLIM1, NLIM2
               V=V1PO+ REAL(J-1)*DVPO
               BB = planck(V,XKTBND)
               sRad(J) = BB*EMISIV
               !print *,'sRad(J)', sRad(J)
               !print *,'emisiv', emisiv
               sEmis(J) = EMISIV
               EMISIV = EMISIV+EMDEL
            ENDDO

         elseif (present(sRad)) then

            DO J = NLIM1, NLIM2
               V=V1PO+ REAL(J-1)*DVPO
               BB = planck(V,XKTBND)
               sRad(J) = BB*EMISIV
               EMISIV = EMISIV+EMDEL
            ENDDO

         elseif (present(sEmis)) then

            DO J = NLIM1, NLIM2
               sEmis(J) = EMISIV
               EMISIV = EMISIV+EMDEL
            ENDDO

         endif

      IF (NLIM2.LT.NLIMO) GO TO 60

   END SUBROUTINE

   !--------------------------------------------------------------------
   !     FUNCTION EMISFN CALCULATES BOUNDARY EMISSIVITY FOR WAVE NUMBER
   !     VALUE CORRESPONDING TO VI AND VINEM, AND THEN CALCULATES THE
   !     LINEAR CHANGE BETWEEN THE EMISSIVITY VALUES AT VI AND VINEM
   !                                                                          !
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
   !
   !               WORK SUPPORTED BY:    THE ARM PROGRAM
   !                                     OFFICE OF ENERGY RESEARCH
   !                                     DEPARTMENT OF ENERGY
   !
   !
   !      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
   !
   !                                             FASCOD3
   ! --------------------------------
   !
   !   EQUIVALENCE (BNDEMI(1),A) , (BNDEMI(2),B) , (BNDEMI(3),C)
   !
   !--------------------------------------------------------------------
   real FUNCTION EMISFN( VI,DVI,VINEM,EMDEL,EMLAST, surf )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real(r8)            ,intent(in)    :: VI
      real                ,intent(in)    :: DVI
      real(r8)            ,intent(inout) :: VINEM
      real                ,intent(out)   :: EMDEL
      real                ,intent(inout) :: EMLAST
      type(CLBLM_Surface) ,intent(in)    :: surf


      character(*) ,parameter :: routineName = 'EMISFN'
      integer  ,PARAMETER :: I_1 = 1
      real     ,PARAMETER :: FACTOR = 0.001

      real      :: A,B,C
      real(r8)  :: V1EMIS,V2EMIS
      real      :: DVEMIS
      integer   :: NLIMEM
      !real      :: ZEMIS(NMAXCO)

      integer  :: NELMNT
      real(r8) :: V1A,V1B
      real     :: ZINT,ZDEL
      real     :: XVI
      real     :: XVNEXT
      integer  :: INTVLS
      real     :: EMNEXT



      if ( surf%nsf>1 ) then !IF (A.LT.0.) THEN !---Surface emissivity input as a look-up-table

         !--->>> 150922 ytma
         V1EMIS = surf%surfEmRfGrid( 1 )
         !print *,'V1EMIS from look-up-table', V1EMIS
         V2EMIS = surf%surfEmRfGrid( surf%nsf )
         DVEMIS = surf%surfEmRfGrid(2) - surf%surfEmRfGrid(1)
         NLIMEM = surf%nsf
         if ( any( [surf%surfEmRfGrid(2:NLIMEM) - surf%surfEmRfGrid(1:NLIMEM-1)] - DVEMIS >1e-6 ) ) then
            STOP '--- '//routineName//'(): Emissivity table needs on a uniform grid.'
         endif
         !ZEMIS(:) = surf%emisTbl%ZEMIS(:)
         !---<<< 150922 ytma


         !--- Determine elements of EMISSION function to use with
         !    input frequency
         NELMNT = INT((VI-V1EMIS)/DVEMIS)

         !--- Test for bounds on EMISSION function
         IF ((NELMNT.LT.0).OR.(NELMNT.GE.NLIMEM-1)) THEN
            WRITE(*,*) 'Frequency range of calculation exceeded',       &
            ' emissivity input.'
            WRITE(*,*) ' VI = ',VI,' V1EMIS = ',V1EMIS,' V2EMIS = ',    &
            V2EMIS
            STOP 'ERROR IN EMISFN'
         ENDIF

         !--- Interpolate to obtain appropriate EMISSION value
         V1A = V1EMIS+DVEMIS*NELMNT
         V1B = V1EMIS+DVEMIS*(NELMNT+1)
         CALL LINTCO(V1A,surf%surfEm(NELMNT+1),&
                     V1B,surf%surfEm(NELMNT+2), VI,ZINT,ZDEL)
         EMISFN = ZINT
         VINEM = V1B
         EMDEL = ZDEL*DVI
         EMLAST = surf%surfEm(NELMNT+1)

      elseif (surf%nsf==1) then !--- Emissivity input as a constant value

         A = surf%surfEm(1)
         print *,'constant emis from emisfn in clblm_scene', A
         B = 0. !surf%surfEm(2)
         C = 0. !surf%surfEm(3)

         !--- CHECK FOR CONSTANT E (INDEPENDENT OF VI)
         !    IF CONSTANT RETURN LARGE VALUE FOR VINEM
         IF (B.EQ.0..AND.C.EQ.0.) THEN
            EMISFN = A
            VINEM = 9.99E+9
            EMDEL = 0.0
            EMLAST = EMISFN

            RETURN
         ENDIF

         XVI = VI
         IF (EMLAST.LT.0.) THEN
            EMLAST = A+B*XVI+C*XVI*XVI
         ENDIF

         !--- SET EMISFN EQUAL TO EMISSIVITY AT VI
         !    EMLAST IS EMISFN(VI) FOR EACH SUBSEQUENT CALL
         EMISFN = EMLAST

         IF (VINEM.GE.0.0) THEN
            XVNEXT = XVI+FACTOR/ABS((B+2.*C*XVI))
            XVNEXT = MIN(XVNEXT,(XVI+DVI*2400))
            INTVLS = (XVNEXT-XVI)/DVI
            INTVLS = MAX(INTVLS,I_1)
            XVNEXT = XVI+DVI* REAL(INTVLS)
         ELSE
            XVNEXT = ABS(VINEM)
            INTVLS = (XVNEXT-XVI)/DVI
            INTVLS = MAX(INTVLS,I_1)
         ENDIF

         EMNEXT = A+B*XVNEXT+C*XVNEXT*XVNEXT

         EMDEL = (EMNEXT-EMISFN)/ REAL(INTVLS)

         VINEM = XVNEXT
         EMLAST = EMNEXT

      endif !if ( surf%nsf>1 )

   END FUNCTION


   !--------------------------------------------------------------------
   ! Returns surface reflectance. This is a wrapper for getSurfRefl_array.
   ! * Please be noted that sfcRad%indV1 = 1
   !--------------------------------------------------------------------
   SUBROUTINE getSurfRefl_spect( surf, V1,DV,NLIM, sfcRefl )
   !--------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, &
                                    CLBLM_Spectrum_init

      type(CLBLM_Surface)  ,intent(in)  :: surf
      real(r8)             ,intent(in)  :: V1
      real                 ,intent(in)  :: DV
      integer              ,intent(in)  :: NLIM
      type(CLBLM_Spectrum) ,intent(out) :: sfcRefl

      real,allocatable :: sRefl(:)


      allocate( sRefl( NLIM ))

      call getSurfRefl_array( surf, V1,DV,NLIM, sRefl )

      call CLBLM_Spectrum_init( sfcRefl, V1,DV,NLIM )

      call move_alloc( sRefl, sfcRefl%spect )

   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE getSurfRefl_array( surf, V1PO,DVPO,NLIMO, sRefl )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8

      type(CLBLM_Surface) ,intent(in)  :: surf
      real(r8)            ,intent(in)  :: V1PO
      real                ,intent(in)  :: DVPO
      integer             ,intent(in)  :: NLIMO
      real                ,intent(out) :: sRefl(:) !(NLIMO)

      integer  :: I,J,K
      real     :: REFLCT, RFDEL, RFLAST
      real(r8) :: V, VI, VIDV
      integer  :: NLIM1, NLIM2


      VI     = V1PO-DVPO
      VIDV   = VI
      RFLAST = -1.
      NLIM1  = 0
      NLIM2  = 0

   60 NLIM1 = NLIM2+1

         VI = V1PO+ REAL(NLIM1-1)*DVPO
         REFLCT = REFLFN(VI,DVPO,VIDV,RFDEL,RFLAST, surf)

         IF (VIDV.GE.9.E+4) THEN
            NLIM2 = NLIMO+1
         ELSE
            NLIM2 = (VIDV-V1PO)/DVPO+1.001
         ENDIF
         NLIM2 = MIN(NLIM2,NLIMO)

         DO J = NLIM1, NLIM2
            sRefl(J) = REFLCT
            REFLCT = REFLCT+RFDEL
         ENDDO

      IF (NLIM2.LT.NLIMO) GO TO 60

   END SUBROUTINE



   !--------------------------------------------------------------------
   ! FUNCTION REFLFN CALCULATES BOUNDARY REFLECTIVITY FOR WAVE NUMBER
   ! VALUE CORRESPONDING TO VI AND VINRF, AND THEN CALCULATES THE
   ! LINEAR CHANGE BETWEEN THE REFLECTIVITY VALUES AT VI AND VINRF
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
   !--------------------------------------------------------------------
   real FUNCTION REFLFN (VI,DVI,VINRF,RFDEL,RFLAST, surf)
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real(r8)            ,intent(in)    :: VI
      real                ,intent(in)    :: DVI
      real(r8)            ,intent(inout) :: VINRF
      real                ,intent(out)   :: RFDEL
      real                ,intent(inout) :: RFLAST
      type(CLBLM_Surface) ,intent(in)    :: surf


      character(*) ,parameter :: routineName = 'REFLFN'
      real    ,PARAMETER :: FACTOR = 0.001
      integer ,PARAMETER :: I_1 = 1

      real     :: A,B,C
      real(r8) :: V1RFLT,V2RFLT
      real     :: DVRFLT
      integer  :: NLIMRF
      !real     :: ZRFLT(NMAXCO)

      integer  :: NELMNT
      real(r8) :: V1A,V1B
      real     :: ZINT,ZDEL
      real     :: XVI
      real     :: XVNEXT
      integer  :: INTVLS
      real     :: RFNEXT



      if ( surf%nsf>1 ) then !IF (A.LT.0.) THEN  !!--- Surface reflectance input as a look-up-table

         !--->>> 150922 ytma
         V1RFLT = surf%surfEmRfGrid( 1 )
         V2RFLT = surf%surfEmRfGrid( surf%nsf )
         DVRFLT = surf%surfEmRfGrid(2) - surf%surfEmRfGrid(1)
         NLIMRF = surf%nsf
         if ( any( [surf%surfEmRfGrid(2:NLIMRF) - surf%surfEmRfGrid(1:NLIMRF-1)] - DVRFLT >1e-6 ) ) then
            STOP '--- '//routineName//'(): surfRefl table needs to be on uniform grid.'
         endif
         !ZRFLT(:) = surf%reflTbl%ZRFLT(:)
         !---<<< 150922 ytma


         !--- Determine elements of REFLECTION function to use with
         !    input frequency
         NELMNT = INT((VI-V1RFLT)/DVRFLT)

         !--- Test for bounds on REFLECTION function
         IF ((NELMNT.LT.0).OR.(NELMNT.GE.NLIMRF-1)) THEN
            WRITE(*,*) 'Frequency range of calculation exceeded',       &
            ' reflectivity input.'
            WRITE(*,*) ' VI = ',VI,' V1RFLT = ',V1RFLT,' V2RFLT = ',    &
            V2RFLT
            STOP 'ERROR IN REFLFN'
         ENDIF


         !--- Interpolate to obtain appropriate reflection value
         V1A = V1RFLT+DVRFLT*NELMNT
         V1B = V1RFLT+DVRFLT*(NELMNT+1)
         CALL LINTCO(V1A,surf%surfRefl(NELMNT+1),&
                     V1B,surf%surfRefl(NELMNT+2), VI,ZINT,ZDEL)
         REFLFN = ZINT
         VINRF = V1B
         RFDEL = ZDEL*DVI
         RFLAST = surf%surfRefl(NELMNT+1)

      elseif (surf%nsf==1) then !--- Surface reflectane input as a quadratic formula

         A = surf%surfRefl(1)
         B = 0. !surf%surfRefl(2)
         C = 0. !surf%surfRefl(3)

         !--- CHECK FOR CONSTANT R (INDEPENDENT OF VI)
         !    IF CONSTANT RETURN LARGE VALUE FOR VINRF
         IF (B.EQ.0..AND.C.EQ.0.) THEN
            REFLFN = A
            VINRF = 9.99E+9
            RFDEL = 0.0
            RFLAST = REFLFN

            RETURN
         ENDIF


         XVI = VI
         IF (RFLAST.LT.0.) THEN
            RFLAST = A+B*XVI+C*XVI*XVI
         ENDIF

         !--- SET REFLFN EQUAL TO REFLECTIVITY AT VI
         !    RFLAST IS REFLFN(VI) FOR EACH SUBSEQUENT CALL
         REFLFN = RFLAST

         IF (VINRF.GE.0.0) THEN
            XVNEXT = XVI+FACTOR/ABS((B+2.*C*XVI))
            XVNEXT = MIN(XVNEXT,(XVI+DVI*2400))
            INTVLS = (XVNEXT-XVI)/DVI
            INTVLS = MAX(INTVLS,I_1)
            XVNEXT = XVI+DVI* REAL(INTVLS)
         ELSE
            XVNEXT = ABS(VINRF)
            INTVLS = (XVNEXT-XVI)/DVI
            INTVLS = MAX(INTVLS,I_1)
         ENDIF

         RFNEXT = A+B*XVNEXT+C*XVNEXT*XVNEXT

         RFDEL = (RFNEXT-REFLFN)/ REAL(INTVLS)

         VINRF = XVNEXT
         RFLAST = RFNEXT

      endif !if ( surf%nsf>1 ) then

   END FUNCTION


   !--------------------------------------------------------------------
   ! Linearly interpolates emission and reflection values which
   ! are directly read in from ASCII files
   !--------------------------------------------------------------------
    SUBROUTINE LINTCO(V1,Z1,V2,Z2,VINT,ZINT,ZDEL)
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real(r8) ,intent(in)  :: V1,V2,VINT
      real     ,intent(in)  :: Z1,Z2
      real     ,intent(out) :: ZINT, ZDEL

      real :: ZCEPT

      ! ZDEL is the slope of the line
      ZDEL = (Z2-Z1)/(V2-V1)

      ! ZCEPT is the intercept for V = 0.0
      ZCEPT = Z1 - ZDEL*V1

      ! Calculate ZINT value at VINT
      ZINT = ZDEL*VINT + ZCEPT

   END SUBROUTINE



END MODULE
