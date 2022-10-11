!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_AtmPath
   USE Module_ConstParam ,ONLY: FILLINT, FILLREAL

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: &
             CLBLM_PathCompound, &
             CLBLM_PathCompound_init, &
             CLBLM_Path, &
             CLBLM_Path_init, &
             CLBLM_Layer, &
             CLBLM_Layer_init, &
             buildRTgrid, &
             calculatePath, &
             calcCLBLMPaths, &
             loadLayerFromPath, &
             extendPath, &
             densityOnPathLevel, &
             scalePathAmount_interface, &
             CLBLM_PathAmountModifier, &
             autoLayering



   TYPE :: CLBLM_PathLevelData
      integer           :: nRTlev =0  !number of levels
      real, allocatable :: zRT(:)     ![nLev],level heights
      real, allocatable :: pRT(:)     ![nLev],level pressures
      real, allocatable :: T(:)       ![nLev],level temperatures
      real, allocatable :: Q(:,:)     ![nLev,nMol],level concentration
   END TYPE


   TYPE :: CLBLM_PathGeom !LBLRTM standard path parameters.
      real    :: Hobs   = FILLREAL !(Km) Observer hight
      real    :: Hend   = FILLREAL !(Km) final hight
      real    :: ANGLE  = FILLREAL !(Deg) path zenith angle at Hobs
      real    :: PHI    = FILLREAL !(Deg) final (Hend) zenith angle along the path
      real    :: HMIN   = FILLREAL !(Km) The minimum height along the path
      integer :: LEN    = FILLINT  !(NA) Short/long path through a tangent height.
   END TYPE


   TYPE ,EXTENDS(CLBLM_PathLevelData) :: CLBLM_Path
      integer                    :: nLay =0       !number of layers
      integer                    :: nMol =0       !num of molecules, array dimmension.
      character(20) ,allocatable :: molID(:)      ![nMol], Molecule identifiers; Molecular name; Isotopologue tag is like "CO2_727"
      real          ,allocatable :: W(:,:)        !(cm-2), [nMol,nLay],layer column amounts
      real          ,allocatable :: Wtot(:)       !(cm-2) [nLay],layer column amounts for broadening gases.
      real          ,allocatable :: Pave(:)       !(hPa)  [nLay],average pressures for each layer
      real          ,allocatable :: Tave(:)       !(K)    [nLay],average temporatures for each layer
      !real          ,allocatable :: airMass(:)    !(n/a) [nLay] air mass factor = slantPathAmount/verticalPathAmount
      integer       ,allocatable :: IPATH(:)      !(n/a) [nLay],relationship between the layer and the rest of the path.
      type(CLBLM_PathGeom)       :: geom          ! Store the final path geometry information.

      real          ,allocatable :: dTave_dTup(:) !(K/K)       [nLay]      Tave derivative wrt upper level temperature
      real          ,allocatable :: dTave_dTlo(:) !(K/K)       [nLay]      Tave derivative wrt lower level temperature
      real          ,allocatable :: dW_dQup(:)    !(cm-3/cm-2) [nLay,nMol] W derivative wrt upper level concentration
      real          ,allocatable :: dW_dQlo(:)    !(cm-3/cm-2) [nLay,nMol] W derivative wrt lower level concentration
   END TYPE


   TYPE :: CLBLM_PathCompound
      integer          :: obsLev           = FILLINT !observer level
      logical          :: AMscalingThermal = .FALSE. !
      logical          :: AMscalingSolar   = .FALSE. !
      integer          :: refPath          = 1       !if=0, vertical; if=1, up path
      type(CLBLM_Path) :: vert
      type(CLBLM_Path) :: view
      type(CLBLM_Path) :: down
      type(CLBLM_Path) :: sun
   END TYPE



   TYPE :: CLBLM_Layer
      integer                   :: layNo  =0      !layer number
      integer                   :: nMol =0        !num of molecules
      character(20),allocatable :: molID(:)      !molecular name; Isotopologue tag is like "CO2_727"
      real         ,allocatable :: W(:)           !(cm-2), [nMol],molecular layer amounts
      real                      :: Wtot           !Wbroad =0.    !(cm-2) layer column amounts for broadening gases.
      real                      :: Pave =0.       !(hPa)  aveage pressures for the layer
      real                      :: Tave =0.       !(K)    average temporatures for the layer

      real                      :: Ztop =0.       !top level height
      real                      :: Zbot =0.       !bottom level height
      real                      :: Ptop =0.       !top level pressures
      real                      :: Pbot =0.       !bottom level pressures
      real                      :: Ttop =0.       !top level temporatures
      real                      :: Tbot =0.       !bottom level temporatures
      real         ,allocatable :: Qtop(:)        !top level concentration
      real         ,allocatable :: Qbot(:)        !bottom level concentration

      real                      :: airMass = 0.   !airMass factor
      integer                   :: IPATH = 0
   END TYPE




   !--- TYPE CLBLM_PathAmountModifier
   ! This type is solely used for uses to implement their own path amount scaling
   ! routine. If no user-provided scaling subroutine, the default scaling subroutine
   ! contained here will be called. To implement a user subroutine, the user
   ! need to extend the "CLBLM_PathAmountModifier" type and write a subroutine to override
   ! the "scalePathAmount()" contained here.
   !
   ! For example,
   !
   ! TYPE, EXTENDS(CLBLM_PathAmountModifier) :: userPathAmountModifier
   ! CONTAINS
   !    procedure ,nopass :: scalePathAmount => userSubroutineToScalePathAmount
   ! END TYPE
   !
   ! INTERFACE
   !    subroutine userSubroutineToScalePathAmount( wkl, wbrodl, nMol, nLnMol, nLay, molID )
   !        real         ,intent(inout) :: wkl(nMol,nLay) !(cm-2), Absorber column amount
   !        real         ,intent(inout) :: wbrodl(nLay)   !(cm-2)Column amount other than line molecules
   !        integer      ,intent(in)    :: nMol           !Number of absorbing molecules, = nLnMOl + nXsMol
   !        integer      ,intent(in)    :: nLnMol         !Number of line molecules contained in wkl
   !        integer      ,intent(in)    :: nLay           !Number of layers
   !        character(*) ,intent(in)    :: molID(nMol)   !Molecular name
   !    end subroutine
   ! END INTERFACE
   !
   ! Then in user's program, user define a type(userPathAmountModifier) and
   ! call testscal%scalePathAmount().
   !
   ! type(userPathAmountModifier) :: testscal
   ! real                         :: testW(1,1)=[1.0], testWb(1)=1.0
   !
   ! call testscal%scalePathAmount( testW,testWb,1,1,1,'H2O')
   !
   TYPE CLBLM_PathAmountModifier
   CONTAINS
      procedure, nopass :: scalePathAmount
   END TYPE


CONTAINS !=================== MODULE CONTAINS ==========================


   !-------------------------------------------------------------------
   !* Boundary array subscript may start from 0 or 1, depending on the value of the
   !  optional argument 'levelStartFrom'.
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_PathLevelData_init(this, nMol, nLev, levelStartFrom)
   !-------------------------------------------------------------------
      type(CLBLM_PathLevelData)   ,intent(inout) :: this !out
      integer                     ,intent(in)    :: nMol
      integer                     ,intent(in)    :: nLev
      integer           ,optional ,intent(in)    :: levelStartFrom

      integer :: lev1, lev2

      if (present(levelStartFrom)) then
         lev1 = levelStartFrom
         lev2 = lev1 + nLev -1
      else
         lev1 = 1
         lev2 = nLev
      endif

      if (allocated(this%zRT))  deallocate(this%zRT)
      if (allocated(this%pRT))  deallocate(this%pRT)
      if (allocated(this%T))    deallocate(this%T)
      if (allocated(this%Q))    deallocate(this%Q)

      allocate( this%zRT(    lev1:lev2))   ;this%zRT(:) = 0.
      allocate( this%pRT(    lev1:lev2))   ;this%pRT(:) = 0.
      allocate( this%T(      lev1:lev2))   ;this%T(:) = 0.
      allocate( this%Q(nMol, lev1:lev2))   ;this%Q(:,:) = 0.

      this%nRTlev = nLev

   END SUBROUTINE


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Path_init( this, nMol, nLay )
   !-------------------------------------------------------------------
      type(CLBLM_Path)  ,intent(out) :: this
      integer           ,intent(in)  :: nMol
      integer           ,intent(in)  :: nLay


      !--- parent type initialization
      call CLBLM_PathLevelData_init( this%CLBLM_PathLevelData, nMol, nLay+1, &
                                     levelStartFrom =1 ) !Level data subscript start from 1

      if (allocated(this%molID))    deallocate(this%molID)
      if (allocated(this%W))        deallocate(this%W)
      if (allocated(this%Wtot))     deallocate(this%Wtot)
      if (allocated(this%Pave))     deallocate(this%Pave)
      if (allocated(this%Tave))     deallocate(this%Tave)
      !if (allocated(this%airMass))  deallocate(this%airMass)
      if (allocated(this%IPATH))    deallocate(this%IPATH)
      !if (allocated(this%SECANT))   deallocate(this%SECANT)
      !if (allocated(this%Wbroad))   deallocate(this%Wbroad)

      allocate( this%molID(  nMol))  ;this%molID(:)   = ''
      allocate( this%W( nMol,nLay))  ;this%W(:,:)     = 0.
      allocate( this%Wtot(   nLay))  ;this%Wtot(:)    = 0.
      allocate( this%Pave(   nLay))  ;this%Pave(:)    = 0.
      allocate( this%Tave(   nLay))  ;this%Tave(:)    = 0.
      !allocate( this%airMass(nLay))  ;this%airMass(:) = 0.
      allocate( this%IPATH(  nLay))  ;this%IPATH(:)   = 0
      !allocate( this%SECANT( nLay))  ;this%SECANT(:)  = 0.
      !allocate( this%Wbroad(nLay))   ;this%Wbroad(:)  = 0.

      !if (allocated(this%dTave_dTup) deallocate(this%dTave_dTup)
      !if (allocated(this%dTave_dTlo) deallocate(this%dTave_dTlo)
      !if (allocated(this%dW_dQup)    deallocate(this%dW_dQup)
      !if (allocated(this%dW_dQlo)    deallocate(this%dW_dQlo)
      !allocate( this%dTave_dTup(   nLay))  ;this%dTave_dTup(:) = 0.
      !allocate( this%dTave_dTlo(   nLay))  ;this%dTave_dTlo(:) = 0.
      !allocate( this%dW_dQup( nMol,nLay))  ;this%dW_dQup(:,:)  = 0.
      !allocate( this%dW_dQlo( nMol,nLay))  ;this%dW_dQlo(:,:)  = 0.

      this%nMol = nMol
      this%nLay = nLay

   END SUBROUTINE


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_PathCompound_init( this, nMol, nLay )
   !-------------------------------------------------------------------
      type(CLBLM_PathCompound)  ,intent(out) :: this
      integer                   ,intent(in)  :: nMol
      integer                   ,intent(in)  :: nLay

      call CLBLM_Path_init( this%vert, nMol, nLay )
      call CLBLM_Path_init( this%view, nMol, nLay )
      call CLBLM_Path_init( this%down, nMol, nLay )
      call CLBLM_Path_init( this%sun,  nMol, nLay )

   END SUBROUTINE


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   SUBROUTINE CLBLM_Layer_init(this, nMol, layerNum)
   !-------------------------------------------------------------------
      type(CLBLM_Layer) ,intent(inout) :: this
      integer           ,intent(in)    :: nMol
      integer           ,intent(in)    :: layerNum
      !integer           ,intent(in)    :: nISO   !number of isotopologues

      if (allocated(this%molID))   deallocate(this%molID)
      if (allocated(this%W))       deallocate(this%W)
      if (allocated(this%Qtop))    deallocate(this%Qtop)
      if (allocated(this%Qbot))    deallocate(this%Qbot)

      allocate( this%molID(  nMol)) ;this%molID(:)  = ''
      allocate( this%W(      nMol)) ;this%W(:)       = 0.
      allocate( this%Qtop(   nMol)) ;this%Qtop(:)    = 0.
      allocate( this%Qbot(   nMol)) ;this%Qbot(:)    = 0.

      this%nMol  = nMol
      this%layNo = layerNum

   END SUBROUTINE



   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE loadLayerFromPath( aLayer, path, iLay )
   !--------------------------------------------------------------------
      type(CLBLM_Layer),  intent(inout) :: aLayer
      type(CLBLM_Path),   intent(in)    :: path
      integer,            intent(in)    :: iLay

      integer :: nMol


      nMol = path%nMol

      call CLBLM_Layer_init( aLayer, nMol, iLay )

      aLayer%molID(1:nMol) = path%molID(1:nMol)
      aLayer%W(    1:nMol) = path%W(    1:nMol,iLay)

      aLayer%Wtot    = path%Wtot(iLay)
      aLayer%Pave    = path%Pave(iLay)
      aLayer%Tave    = path%Tave(iLay)

      aLayer%Ztop    = path%zRT(iLay+1)
      aLayer%Zbot    = path%zRT(iLay)
      aLayer%Ptop    = path%pRT(iLay+1)
      aLayer%Pbot    = path%pRT(iLay)
      aLayer%Ttop    = path%T(iLay+1)
      aLayer%Tbot    = path%T(iLay)
      aLayer%Qtop(:) = path%Q(:,iLay+1)
      aLayer%Qbot(:) = path%Q(:,iLay)

      !aLayer%airMass = path%airMass(iLay)   !airMass factor
      aLayer%IPATH   = path%IPATH(iLay)
      !aLayer%SECANT  = path%SECANT(iLay)    !"effective" secant of zenith angle, SECNTA(L) = SOUT(L)/(ZOUT(L+1)-ZOUT(L))
      !aLayer%Wbroad  = path%Wbroad(iLay)

   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE scalePathAmount_interface( path )
!-----------------------------------------------------------------------
      type(CLBLM_Path)  ,intent(inout) :: path

      type(CLBLM_PathAmountModifier) :: amtScl
      integer                        :: il
      integer                        :: nMol          !Number of absorbing molecules
      integer                        :: nLay          !Number of layers
      character(20) ,allocatable     :: molID(:)      !Molecule identifiers;
      real          ,allocatable     :: wkl(:,:)      !(cm-2), Absorber column amount
      real          ,allocatable     :: wtot(:)       !(cm-2), Total column amount


      nMol = path%nMol
      nLay = path%nLay

      allocate( molID( nMol), &
                wkl( nMol,nLay), &
                wtot( nLay) )

      molID( 1:nMol )      = path%molID( 1:nMol )  !Molecule names
      wkl( 1:nMol,1:nLay ) = path%W( 1:nMol,1:nLay) !(cm-2), Absorber column amount
      wtot( 1:nLay )       = path%Wtot( 1:nLay)     !(cm-2)Total Column amount

      call amtScl % scalePathAmount( wkl, wtot, nMol, nLay, molID )

      path%W( 1:nMol,1:nLay ) = wkl( 1:nMol,1:nLay )
      path%Wtot(     1:nLay ) = wtot( 1:nLay )


      deallocate(molID)
      deallocate(wkl)
      deallocate(wtot)

   END SUBROUTINE



!-----------------------------------------------------------------------
! * The input column amounts are the amounts along the viewing path
!   for each layer.
! * "wkl(nMol,nLay)" contains both line molecules and cross-section molecules.
!   This example assumes no isotopologues included in the molecular list.
!   If isotopologues are included, one needs to be careful to keep the
!   iso fraction not changed.
!-----------------------------------------------------------------------
   SUBROUTINE  scalePathAmount( wkl, wtot, &
                                nMol, nLay, molID )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: molNum, xsMolNum, molIndex
      IMPLICIT NONE

      real         ,intent(inout) :: wkl(nMol,nLay) !(cm-2), Absorber column amount
      real         ,intent(inout) :: wtot(nLay)     !(cm-2), Total column amount
      integer      ,intent(in)    :: nMol           !Number of absorbing molecules, = nLnMOl + nXsMol
      integer      ,intent(in)    :: nLay           !Number of layers
      character(*) ,intent(in)    :: molID(nMol)   !Molecular name



      !--------- Below is an example of path scaling code ---------
      !------------------------------------------------------------
      !  *** It should be noted that no attempt has been made to keep the
      !      mass in a given layer constant, i.e. level pressure nor retained
      !
      integer      :: m,l,k, iH2O
      logical      :: isH2O, isXSMol
      character(1) :: hmol_scal(nMol)
      real         :: xmol_scal(nMol)
      real         :: xmol_scal_m
      real         :: wmt(nMol), wbroad(nLay), wsum_brod, wsum_drair


      if (nMol<=0) RETURN !No need to do path scaling, return

      hmol_scal(1:nMol) = ''
      xmol_scal(1:nMol) = 1.


      do m = 1, nMol
         wmt(m) = 0.
      enddo

      do m = 1, nmol
         do l = 1, nLay
            wmt(m) = wmt(m) + wkl(m,l)
         enddo
      enddo

      iH2O = molIndex('H2O', molID)
      if ( iH2O<0 ) then
         STOP '--- scalePathAmount(): H2O is not present.'
      endif
      wsum_drair = sum( wtot(1:nLay) ) - wmt( iH2O )

      !---Save the amounts before scaling
      ! This example assumes no isotopologues included in the molecular list.
      ! If isotopologues are included, one needs to be careful to keep the
      ! iso fraction not changed.
      do l = 1,nLay
         wbroad(l) = wtot(l) - sum(wkl(:,l))
      enddo

!if (nXsMol>0) hmol_scal(nLnMol+1:nMol) = '1' !Since in in LBLRTM x_xs_scal(m) = x_xs_scal_m

      do m = 1, nMol

         isXSMol = (xsMolNum(molID(m)) > 0)
         if (isXSMol) hmol_scal(m)= '1'

         isH2O = (molNum( molID(m) ) == molNum('H2O'))

         xmol_scal_m = xmol_scal(m)
         xmol_scal(m) = -999.9

         if (hmol_scal(m).eq.' ') xmol_scal(m) = 1.
         if (hmol_scal(m).eq.'0') xmol_scal(m) = 0.
         if (hmol_scal(m).eq.'1') xmol_scal(m) = xmol_scal_m

         if (hmol_scal(m).eq.'C' .or. hmol_scal(m).eq.'c')           &
              xmol_scal(m) = xmol_scal_m/wmt(m)

         if (hmol_scal(m).eq.'M' .or. hmol_scal(m).eq.'m')           &
              xmol_scal(m) = xmol_scal_m/(wmt(m)/wsum_drair)

         if ((hmol_scal(m).eq.'P' .or. hmol_scal(m).eq.'p') .and. isH2O)  &
              xmol_scal(m) = (xmol_scal_m/2.99150e-23)/wmt(m)
!             value from vpayne 2006/07/24

         if (hmol_scal(m).eq.'D' .or. hmol_scal(m).eq.'d')           &
              xmol_scal(m) =  (xmol_scal_m*2.68678e16)/wmt(m)

         if ((hmol_scal(m).eq.'P' .or. hmol_scal(m).eq.'p') .and. .not.isH2O) then
            write (*,*) 'm = ', m !yma write (ipr,*) 'm = ', m
            stop ' (hmol_scal(m).eq."P" .and. m.ne.1) '
         endif


         if (xmol_scal(m) .lt. -998.) then
            !write (ipr,*) 'm = ', m,' h_mol_scal(m) not valid '
            write (  *,*) 'm = ', m,' h_mol_scal(m) not valid '
            stop
         endif

         do l = 1, nLay
            wkl(m,l) = wkl(m,l) * xmol_scal(m)
         enddo

         wmt(m) = wmt(m)*xmol_scal(m)

      enddo !do m = 1, nMol


      !--- Recalculate the total amounts
      do l=1,nLay
         wtot(l) = wbroad(l) + sum(wkl(:,l))
      enddo

   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE calcCLBLMPaths( scene, pathCtrl, paths, hdrRTgrid,hdrRTsize,flux_Flags )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: PI, secant_diffuse, molIndex
      USE Module_Config     ,ONLY: CLBLM_Flux_Ctrl
      USE Module_Config     ,ONLY: CLBLM_Path_Ctrl
      USE Module_Scene      ,ONLY: CLBLM_scene, &
                                   CLBLM_SceneGeom
      IMPLICIT NONE

      type(CLBLM_Scene)            ,intent(in)    :: scene
      type(CLBLM_Path_Ctrl)        ,intent(in)    :: pathCtrl
      type(CLBLM_PathCompound)     ,intent(inout) :: paths
      real               ,OPTIONAL ,intent(in)    :: hdrRTgrid(:)
      integer            ,OPTIONAL ,intent(in)    :: hdrRTsize
      type(CLBLM_Flux_Ctrl)        ,intent(in)    :: Flux_flags


      character(*) ,parameter :: routineName = 'calcCLBLMPaths'
      real         ,parameter :: TOL = 5.e-4

      integer               :: il,nBnd,mxLev, obsLevel, nRTlev
      real     ,allocatable :: zRT(:)
      real                  :: obsAlt
      type(CLBLM_Path)      :: path, downPath, solPath, vertPath, extPath
      type(CLBLM_SceneGeom) :: tempGeom
      logical               :: SpaceLook, reflectingSurf
      logical               :: flux_flag




      !--- Build RT grid
      ! * The returned zRT is based on user input pRT or scene%prfl (autolay/profile level),
      !   it is not necessarily include the observer level. The following
      !   calculatePath() for viewing path will include the observer level into viewpath%zrt.
      CALL buildRTgrid( scene%prfl, scene%geom, pathCtrl, zRT,nRTlev, hdrRTgrid,hdrRTsize )


      obsAlt = scene%geom%obs%obsAlt

      !--- Calculate viewing path
      ! * If the obsAlt is in the middle of atmopshere, viewing path is extended in
      !   the backward direction to conver whole atmosphere, unless it is a space-looking path.
      !
      call calculatePath( path, scene%prfl, scene%geom, pathCtrl, zRT,nRTlev )

      SpaceLook = path%geom%Hend > path%geom%HMIN

      !---Check to see if zRT covers the whole atmosphere defined by profile
      if ( .not.SpaceLook ) then
         reflectingSurf = ( allocated(scene%sfc%surfEm)   .and. any(scene%sfc%surfEm <1.)) .OR. &
                          ( allocated(scene%sfc%surfRefl) .and. any(scene%sfc%surfRefl >0.))
         if ( reflectingSurf .and. &
            ( zRT(nRTLev) < scene%prfl%Z( scene%prfl%nLev ) .or. &
              zRT(1)      < scene%prfl%Z( 1 )  )) then
            STOP '--- '//routineName//'(): For down looking sensor over reflecting surface, the RT grid has to cover from surface to TOA.'
         endif
      endif

      !--- For the viewing path, if Hobs is in the middle of atmosphere,
      !  and look down, extends the viewing path to whole atmosphere
      if ( obsAlt < scene%prfl%Z( scene%prfl%nLev ) .and. &
           obsAlt > scene%prfl%Z( 1 ) .and. &
           .not.SpaceLook ) then

         tempGeom = scene%geom
         tempGeom%obs%viewAng = 180. - scene%geom%obs%viewAng
         CALL calculatePath( extPath, scene%prfl, tempGeom, pathCtrl, zRT,nRTlev )
         CALL extendPath( path, extPath, above=.TRUE.) !add extPath above path
      endif



      !--- Calculate downwelling thermal path
      !
      if ( SpaceLook ) then !Up looking case
         downPath = path
      else

         tempGeom = scene%geom

         print *, 'path flux_flag', flux_flags%flux_flag

         print *, 'ThmReflMode', scene%sfc%ThmReflMode

         if ( scene%sfc%ThmReflMode == 1) then !specular reflection
            print *, 'Specular Reflection in calcCLBLMPaths'
            print *, 'Specular PHI', path%geom%PHI
            tempGeom%obs%obsAlt  = path%geom%Hend
            tempGeom%obs%viewAng = path%geom%PHI
            call calculatePath( downPath, scene%prfl, tempGeom, pathCtrl, path%zRT,path%nRTlev )

         elseif ( scene%sfc%ThmReflMode ==0) then !Lambertian reflection, use diffusivity angle
            print *, 'Lambertian Reflection in calcCLBLMPaths'
            print *, 'secant_diffuse', secant_diffuse
            print *, 'Lambertian PHI', path%geom%PHI
            tempGeom%obs%obsAlt  = path%geom%Hend
            if (flux_flags%flux_flag .eqv. .false.) then
               print*, 'use secant of diffusivity angle for path'
               tempGeom%obs%viewAng = acos(1./secant_diffuse)*180./PI
            elseif (flux_flags%flux_flag .eqv. .true.) then
               print*, 'use vertical path'
               tempGeom%obs%viewAng = path%geom%PHI
            endif
            call calculatePath( downPath, scene%prfl, tempGeom, pathCtrl, path%zRT,path%nRTlev )
         endif
      endif


      !--- Calculate solar path
      !
      if ( SpaceLook ) then
         solPath = path
      else
         tempGeom = scene%geom
         tempGeom%obs%obsAlt  = path%geom%Hend    !scene%prfl%Z( scene%prfl%nLev )
         tempGeom%obs%viewAng = scene%geom%sunAng !180. - scene%geom%sunAng
         call calculatePath( solPath, scene%prfl, tempGeom, pathCtrl, path%zRT,path%nRTlev )
      endif



      !--- Calculate vertical path
      !
      !tempGeom = scene%geom
      !
      !if  (path%geom%ANGLE <=90. ) then !Up looking cases ! elseif (path%geom%obsAlt==path%geom%HMIN) then
      !
      !   tempGeom%obs%viewAng = 0.
      !
      !elseif (abs(path%geom%Hend - path%geom%HMIN) <TOL) then !Down looking and line of sight ended at surface
      !
      !   tempGeom%obs%targetAlt = path%geom%HMIN
      !   tempGeom%obs%viewAng = 180.
      !
      !elseif ( path%geom%Hend > path%geom%HMIN ) then !Down looking but tangent line of sight
      !
      !   tempGeom%obs%obsAlt = path%geom%HMIN
      !   tempGeom%obs%viewAng = 0.
      !endif
      !
      tempGeom = scene%geom
      tempGeom%obs%obsAlt = path%geom%HMIN
      tempGeom%obs%viewAng = 0.
      call calculatePath( vertPath, scene%prfl, tempGeom, pathCtrl, path%zRT,path%nRTlev )
      !--- If limb path, keep the IPATH flags. IPATH=2 indicates the layer is passed twice.
      if (SpaceLook) then
         vertPath%IPATH(:) = path%IPATH(:)
      endif


      !--- Make sure the three paths are on the same vertical grid.
      if ( any( abs(path%zRT - downPath%zRT) >TOL ) .or. &
           any( abs(path%zRT - solPath%zRT) >TOL ) .or. &
           any( abs(path%zRT - vertPath%zRT) >TOL ) ) then
         STOP '--- '//routineName//'(): Altitude grids for path, downPath and solPath are not consistent.'
      endif


      !--- Determine oberserv level#
      ! The observer altitude must be one of the layer boundary.
      obsLevel = -1
      do il =1,path%nRTlev
         if ( abs( path%zRT(il) - obsAlt ) <TOL ) then
            obsLevel = il
            EXIT
         endif
      enddo
      if (obsLevel<0) STOP '--- '//routineName//'(): Observer height is not on a layer boundary.'



      !--- Collect the paths in the compound structure
      !
      paths%view   = path
      paths%down   = downPath
      paths%sun    = solPath
      paths%vert   = vertPath
      paths%obsLev = obsLevel

      paths%AMscalingThermal = pathCtrl%AMscalingThermal
      paths%AMScalingSolar   = pathCtrl%AMScalingSolar
      paths%refPath          = pathCtrl%refPath


   END SUBROUTINE



   !-------------------------------------------------------------------
   ! Combines two paths into a single path. The 'extPath' will be added to 'path'.
   ! To the above or below of the 'path' depends on the value of "above", if=.T.
   ! 'extPath' will be added to the above of 'path', if=.F., 'extPath' will be added
   ! to the below of 'path'.
   !-------------------------------------------------------------------
   SUBROUTINE extendPath( path, extPath, above )
   !-------------------------------------------------------------------
      type(CLBLM_Path)  ,TARGET,intent(inout) :: path
      type(CLBLM_Path)  ,TARGET,intent(in)    :: extPath
      logical                  ,intent(in)    :: above

      character(*) ,parameter :: routineName = 'extendPath'

      type(CLBLM_Path)          :: tPath
      integer                   :: nLayLo, nLayHi, nLevLo, frstHi
      integer                   :: nMol
      type(CLBLM_Path) ,POINTER :: ptrLo=>null()
      type(CLBLM_Path) ,POINTER :: ptrHi=>null()


      if (path%nMol /= extPath%nMol) then
         STOP '--- '//routineName//'(): The two paths need to have same molecules.'
      endif


      !--- "path" is always the major one. All level values from "path" are
      ! kept unchanged, that is, the values for the overlap level come from 'path'.
      if (above) then !add extPath above path
         ptrLo  => path
         ptrHi  => extPath
         nLayLo =  path%nLay
         nLayHi =  extPath%nLay
         nLevLo =  nLayLo+1
         frstHi =  2
      else                    !add extPath below path
         ptrLo  => extPath
         ptrHi  => path
         nLayLo =  extPath%nLay
         nLayHi =  path%nLay
         nLevLo =  nLayLo
         frstHi =  1
      endif


      nMol = path%nMol
      CALL CLBLM_Path_init( tPath, nMol, nLayLo+nLayHi )


      !--- The lower part of tPath contains ptrLo
      !
      tPath%zRT(     1:nLevLo) = ptrLo%zRT(     1:nLevLo)
      tPath%pRT(     1:nLevLo) = ptrLo%pRT(     1:nLevLo)
      tPath%T(       1:nLevLo) = ptrLo%T(       1:nLevLo)
      tPath%Q(1:nMol,1:nLevLo) = ptrLo%Q(1:nMol,1:nLevLo)

      tPath%W(1:nMol,1:nLayLo) = ptrLo%W(1:nMol,1:nLayLo)
      tPath%Wtot(    1:nLayLo) = ptrLo%Wtot(    1:nLayLo)
      tPath%Pave(    1:nLayLo) = ptrLo%Pave(    1:nLayLo)
      tPath%Tave(    1:nLayLo) = ptrLo%Tave(    1:nLayLo)
      tPath%IPATH(   1:nLayLo) = ptrLo%IPATH(   1:nLayLo)
      !tPath%airMass( 1:nLayLo) = ptrLo%airMass( 1:nLayLo)

      tPath%molID(1:nMol) = ptrLo%molID(1:nMol)


      !--- The upper part of tPath contains ptrHi
      !
      tPath%zRT(      nLevLo+1 : nLayLo+nLayHi+1) = ptrHi%zRT(     frstHi:nLayHi+1)
      tPath%pRT(      nLevLo+1 : nLayLo+nLayHi+1) = ptrHi%pRT(     frstHi:nLayHi+1)
      tPath%T(        nLevLo+1 : nLayLo+nLayHi+1) = ptrHi%T(       frstHi:nLayHi+1)
      tPath%Q(1:nMol, nLevLo+1 : nLayLo+nLayHi+1) = ptrHi%Q(1:nMol,frstHi:nLayHi+1)

      tPath%W(1:nMol, nLayLo+1 : nLayLo+nLayHi) = ptrHi%W(1:nMol,1:nLayHi)
      tPath%Wtot(     nLayLo+1 : nLayLo+nLayHi) = ptrHi%Wtot(    1:nLayHi)
      tPath%Pave(     nLayLo+1 : nLayLo+nLayHi) = ptrHi%Pave(    1:nLayHi)
      tPath%Tave(     nLayLo+1 : nLayLo+nLayHi) = ptrHi%Tave(    1:nLayHi)
      tPath%IPATH(    nLayLo+1 : nLayLo+nLayHi) = ptrHi%IPATH(   1:nLayHi)
      !tPath%airMass(  nLayLo+1 : nLayLo+nLayHi) = ptrHi%airMass( 1:nLayHi)


      !--- Keep path geometry information (for path only).
      tPath%geom = path%geom

      !--- Return the new extended path
      path = tPath

   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   SUBROUTINE buildRTgrid( prfl, geom, pathCtrl, zRT,zRTsize, hdrRTgrid,hdrRTsize )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: molIndex
      USE Module_Config     ,ONLY: CLBLM_Path_Ctrl
      USE Module_Scene      ,ONLY: CLBLM_Profile, &
                                   CLBLM_SceneGeom, &
                                   chopProfile, &
                                   InterpPressToAlt
      IMPLICIT NONE

      type(CLBLM_Profile)          ,intent(in)    :: prfl
      type(CLBLM_SceneGeom)        ,intent(in)    :: geom
      type(CLBLM_Path_Ctrl)        ,intent(in)    :: pathCtrl
      real            ,allocatable ,intent(out)   :: zRT(:)
      integer                      ,intent(out)   :: zRTsize
      real               ,OPTIONAL ,intent(in)    :: hdrRTgrid(:)
      integer            ,OPTIONAL ,intent(in)    :: hdrRTsize

      character(*) ,parameter :: routineName = 'buildRTgrid'

      integer               :: il, nBnd, mxLev, nExtRTLev
      logical               :: autoLayFlag
      real     ,allocatable :: extZRT(:), combineZRT(:)
      type(CLBLM_Profile)   :: tempPrfl



      !--- Build RT grid
      ! * rtGrid (layer boundaries) should be in altitude unit (Km)
      ! * obsAlt should be on a level of RTgrid
      !
      if ( pathCtrl%RTgridFlag == 1 ) then !RT grid from pRT which is input from the scene file header

         if (.not.present(hdrRTgrid) .or. .not.present(hdrRTsize)) &
            STOP '--- '//routineName//'(): The generic RT grid is not present. '

         ! hdrRTgrid is alread in Km unit. simply copy into zRT array.
         allocate( zRT( hdrRTsize ) )
         zRT(1:hdrRTsize) = hdrRTgrid(1:hdrRTsize)
         zRTsize = hdrRTsize

      elseif ( pathCtrl%RTgridFlag == 2 ) then !auto-layering; Calculate RTgrid internally

         call autoLayering( prfl, pathCtrl, zRT )
         zRTsize = size(zRT)

      elseif ( pathCtrl%RTgridFlag == 3 ) then !Use profile grid as RT grid.

         allocate( zRT(prfl%nLev) )
         zRT(:) = prfl%Z(:)
         zRTsize = size(zRT)
      endif



      !!--- If the zRT doesn't cover the whole profile, calculate the RTgrid for
      !! atmosphere between zRT top and profile top.
      !if ( zRT(zRTsize) < prfl%Z( prfl%nLev )) then
      !
      !   !--- Find the lowest level above zRT
      !   do il=1,prfl%nLev
      !      if ( prfl%Z(il) > zRT(zRTsize) ) exit
      !   enddo
      !
      !   !--- Get prfl segment above zRT top
      !   call chopProfile( prfl, il, prfl%nLev, tempPrfl )
      !
      !   !--- Calculate RT grid for the atmopsher above zRT top
      !   !tempGeom = geom
      !   !tempGeom%obs%viewAng = 180. - geom%obs%viewAng
      !   !CALL buildRTgrid( tempPrfl, tempGeom, pathCtrl, extZRT, nExtRTlev )
      !   call autoLayering( tempPrfl, pathCtrl, extZRT )
      !   nExtRTLev = size(extZRT)
      !
      !   !--- If autolayering gives too many layers, use the prfile levels instead
      !   if (nExtRTlev > tempPrfl%nLev) then
      !      call move_alloc( tempPrfl%Z, extZRT )
      !      nExtRTLev = tempPrfl%nLev
      !   endif
      !
      !   !---Combine the user zRT and the extZRT and update the zRT
      !   allocate(combineZRT( zRTsize+nExtRTlev ))
      !   combineZRT = [zRT,extZRT]
      !   call move_alloc( combineZRT, zRT )
      !   zRTsize = zRTsize+nExtRTlev
      !endif



      !--- Check the RTgrid is in ascending order
      do il = 2,zRTsize
         if (zRT(il) <= zRT(il-1)) then
            STOP '--- '//routineName//'(): Layer boundaries are not in ascending order.'
         endif
      enddo


      !--- Check the RTgrid is within atmosphere
      !///


      !!--- If Hobs is in the middle of RTgrid, insert Hobs in RTgrid to
      !! make sure obsLevel is one of the boundaries in viewPath, downPaht and solarPath.
      !! And the layering results for viewing, down and solar are the same.
      !!--- Insert Hobs into RTgrid
      !if ( obsAlt > zRT(1) .and. obsAlt < zRT( zRTsize ) ) then
      !
      !   CALL insertObsAltInRTgrid( zRT,zRTsize, obsAlt )
      !endif

   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE calculatePath( path, scnPrfl, scnGeom, pathCtrl, zRT,nRTlev )
!-----------------------------------------------------------------------
      USE Module_Config       ,ONLY: CLBLM_Path_Ctrl
      USE Module_Scene        ,ONLY: CLBLM_Profile, &
                                     CLBLM_SceneGeom
      USE Module_ConstParam   ,ONLY: molIndex
      IMPLICIT NONE

      type(CLBLM_Path)       ,intent(out)   :: path    !out
      type(CLBLM_Profile)    ,intent(in)    :: scnPrfl
      type(CLBLM_SceneGeom)  ,intent(in)    :: scnGeom
      type(CLBLM_Path_Ctrl)  ,intent(in)    :: pathCtrl
      real                   ,intent(in)    :: zRT(:)
      integer                ,intent(in)    :: nRTlev


      integer               :: N,M,iAng,il, mxLev
      real ,allocatable     :: refracIndex(:)
      type(CLBLM_PathGeom)  :: pathGeo
      type(CLBLM_SceneGeom) :: vertGeom
      type(CLBLM_Path)      :: vertPath
      real                  :: saveViewAng



      !--- COMPUTE THE REFRACTIVE INDEX PROFILE
      mxLev = scnPrfl%nLev  !prfl%toaLev
      allocate( refracIndex(mxLev) ) ;refracIndex(:)=0.
      CALL refractiveIndex( refracIndex(1:mxLev), &
                            scnPrfl%P(1:mxLev), &
                            scnPrfl%T(1:mxLev), &
                            scnPrfl%Q( molIndex('H2O',scnPrfl%molID),1:mxLev ),&
                            pathCtrl%V_refrac, &
                            mxLev)


      !--- REDUCE SLANT PATH PARAMETERS TO STANDARD FORM
      !* The standard form: [H1,H2,Angle,Phi, Hmin,Len]
      CALL standardPathGeometry( pathGeo, scnPrfl, scnGeom, refracIndex, zRT,nRTlev )


      !--- Merge profile and user input layer boundaries, perform ray
      ! tracing, and pack the final result to LBLRTM computational layers
      CALL rayTraceAndLayerPack( path, &                         !out
                                 scnPrfl, &                      !in
                                 refracIndex(1:scnPrfl%nLev), & !in, refracIndex(1:prfl%toaLev), &
                                 zRT,nRTlev, &                   !in
                                 pathGeo, &                      !in
                                 scnGeom%earthRadius, &          !in
                                 pathCtrl%zeroAbsAmt )           !in


!      !--- Calculate air-mass-factor
!      !
!      if ( path%geom%ANGLE ==0. .or. path%geom%ANGLE ==180.) then
!         path%airmass(:) = 1.
!      else
!
!         vertGeom = scnGeom
!
!         if  (path%geom%ANGLE <=90. ) then !Up looking cases ! elseif (path%geom%obsAlt==path%geom%HMIN) then
!
!            vertGeom%obs%viewAng = 0.
!
!         elseif (path%geom%Hend == path%geom%HMIN) then !Down looking and line of sight ended at surface
!
!            vertGeom%obs%targetAlt = path%geom%HMIN
!            vertGeom%obs%viewAng = 180.
!
!         elseif ( path%geom%Hend > path%geom%HMIN ) then !Down looking but tangent line of sight
!
!            vertGeom%obs%obsAlt = path%geom%HMIN
!            vertGeom%obs%viewAng = 0.
!         endif
!
!         call standardPathGeometry( pathGeo, scnPrfl, vertGeom, refracIndex, zRT,nRTlev )
!         call rayTraceAndLayerPack( vertPath, &
!                                    scnPrfl, &
!                                    refracIndex(1:scnPrfl%nLev), & !refracIndex(1:prfl%toaLev), &
!                                    zRT,nRTlev, &
!                                    pathGeo, &
!                                    vertGeom%earthRadius, &
!                                    pathCtrl%zeroAbsAmt )
!         do il = 1,path%nLay
!            path%airMass(il) = path%Wtot(il) / vertPath%Wtot(il)
!         enddo
!
!      endif


      !---
      deallocate(refracIndex)

   END SUBROUTINE


!-----------------------------------------------------------------------
! SUBROUTINE refractiveIndex()
!
! COMPUTE THE REFRACTIVE INDEX PROFILE
! EQUATION FOR RFNDXM IS FROM LOWTRAN (REF 3)
!        RFNDXM(IM) = ((77.46+0.459E-8*XVBAR**2)*PM(IM)/TM(IM)-
!    *                (PPH2O/1013.0)*(43.49-0.347E-8*XVBAR**2))*
!    *                1.0E-6
! RFNDXM IS 1.0-INDEX
!-----------------------------------------------------------------------
   SUBROUTINE refractiveIndex(RFNDXM,PM,TM,DENW,XVBAR,IMMAX)
!-----------------------------------------------------------------------
      USE Module_ConstParam,  ONLY: ALOSMT,&
                                    PZERO=>Press1013, TZERO=>Temp273
      IMPLICIT NONE

      real    ,intent(out) :: RFNDXM(:)    !(IMMAX)
      real    ,intent(in)  :: PM(:)        !(IMMAX)
      real    ,intent(in)  :: TM(:)        !(IMMAX)
      real    ,intent(in)  :: DENW(:)      !(IMMAX)  water vapour density
      real    ,intent(in)  :: XVBAR        !         average wavenumber
      integer ,intent(in)  :: IMMAX

      integer :: im
      real :: PPH2O

      DO 170 IM = 1, IMMAX
         !PPH2O = DENM(1,IM)*PZERO*TM(IM)/(TZERO*ALOSMT)
         PPH2O = DENW(IM)*PZERO*TM(IM)/(TZERO*ALOSMT)

         !	 Approximation to refraction index (from LOWTRAN6)
         RFNDXM(IM)=((83.42+(185.08/(1.0-(XVBAR/1.14E+5)**2))+    &
         (4.11/(1.0-(XVBAR/6.24E+4)**2)))*(PM(IM)*288.15)/        &
         (1013.25*TM(IM))-(43.49-(XVBAR/1.7E+4)**2)*(PPH2O/       &
         1013.25)) *1.0E-06
  170 CONTINUE

   END SUBROUTINE



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE standardPathGeometry( pathGeo, &
                                    scnPrfl, scnGeom, refracIndex, zRT,nRTlev )
!-----------------------------------------------------------------------
      USE Module_Scene  ,ONLY: CLBLM_Scene, &
                               CLBLM_Profile, &
                               CLBLM_SceneGeom
      IMPLICIT NONE

      type(CLBLM_PathGeom)  ,intent(inout) :: pathGeo !out
      type(CLBLM_Profile)   ,intent(in)    :: scnPrfl
      type(CLBLM_SceneGeom) ,intent(in)    :: scnGeom
      real                  ,intent(in)    :: refracIndex(:) !refractive index profile
      real                  ,intent(in)    :: zRT(:)
      integer               ,intent(in)    :: nRTlev


      integer :: M,N
      integer :: NMOL
      integer :: ITYPE, IERROR
      !logical :: ZMDL_fromPM
      real    :: H1,H2,ANGLE,RANGE,BETA,PHI,HMIN !,HOBS
      integer :: LEN
      real    :: dummyAlt(1), dummyTemp(1)
      real    :: ZMAX,ZMIN,RE
      integer :: IBMAX,IMMAX

      real ,allocatable :: ZMDL(:),PM(:),TM(:),RFNDXM(:),ZBND(:)
      real ,allocatable :: DENM(:,:)


      RE = scnGeom%earthRadius

      IMMAX = scnPrfl%nLev     !prfl%toaLev
      ZMAX  = scnPrfl%Z(IMMAX) !prfl%Z( prfl%toaLev )
      ZMIN  = scnPrfl%Z(1)

      IBMAX = nRTlev
      allocate( ZBND(IBMAX) )
      ZBND(1:IBMAX) = zRT(1:IBMAX)

      allocate( ZMDL(   IMMAX) )
      allocate( PM(     IMMAX) )
      allocate( TM(     IMMAX) )
      allocate( RFNDXM( IMMAX) )
      ZMDL(  1:IMMAX) = scnPrfl%Z(1:IMMAX)
      PM  (  1:IMMAX) = scnPrfl%P(1:IMMAX)
      TM  (  1:IMMAX) = scnPrfl%T(1:IMMAX)
      RFNDXM(1:IMMAX) = refracIndex(1:IMMAX)


      !---
      ! * FSCGEO does not compute path amount, so it actually does NOT use DENM
      !
      NMOL  = scnPrfl%nMol  !nLnMol + nXsMol
      allocate( DENM(NMOL,IMMAX) )
      DENM(1:NMOL,1:IMMAX) = scnPrfl%Q(1:nMol,1:IMMAX)


      !---
      ! CLBLM input geometry parameters: H2 and ANGLE only.
      ! LBLRMT input geometry parameters: H1,H2,ANGLE,RANGE,BETA
      ! internal used parameters: H1,H2,ANGLE,PHI,HMIN,LEN
      !
      RANGE = 0.
      BETA  = 0.

      ITYPE = 2  !LBLRTM case2A, input H1,H2 and ANGLE
      LEN   = 0
      H1    = scnGeom%obs%obsAlt   !(Km)
      ANGLE = scnGeom%obs%viewAng  !(Deg)
      if (ANGLE <=90.) then
         H2 = ZMAX
      else
         H2 = ZMIN
      endif

      !--- Call LBLRTM subroutines
      !* FSCGEO needs to call RFPATH to perform reduction
      IERROR = 0
      CALL FSCGEO ( H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,IERROR, & !HOBS, &
                    ZBND,ZMAX,ZMIN,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX,RE, &
                    DENM,NMOL )

      !--- If view angle >90 and the path is a tanget path to space (not reach the ZMIN)
      ! recalcualte the standard (internal) path parameters.
      if (IERROR==2) then

         RANGE = 0.
         BETA  = 0.

         ITYPE = 2  !LBLRTM case2A, input H1,H2 and ANGLE
         H1    = scnGeom%obs%obsAlt   !(Km)
         ANGLE = scnGeom%obs%viewAng  !(Deg)
         H2    = ZMAX

         IERROR=0
         CALL FSCGEO ( H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,IERROR, & !HOBS, &
                       ZBND,ZMAX,ZMIN,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX,RE, &
                       DENM,NMOL )
      endif

      if (IERROR.NE.0) then !GO TO 310
         print*, '--- standartPathGeometry(): ',&
                 'ERROR FLAG RETURNED FROM FSCGEO:  AN ERROR OCCURED ',   &
                 'IN PROCESSING THE SLANT PATH PARAMETERS.'
         STOP
      endif

      !--- Save the results into pathGeo structure.
      pathGeo%Hobs  = H1 !in (Km)
      pathGeo%Hend  = H2 !in (Km)
      pathGeo%ANGLE = ANGLE
      pathGeO%PHI   = PHI
      pathGeo%HMIN  = HMIN
      pathGeo%Len   = LEN


      !---
      deallocate( ZMDL,PM,TM,RFNDXM,ZBND )
      deallocate( DENM)

   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE rayTraceAndLayerPack( path, &
                      scnPrfl, refracIndex, zRT,nRTlev, pathGeo, earthRadius, &
                      zeroAbsAmt )
!-----------------------------------------------------------------------
      USE Module_Scene ,ONLY: CLBLM_Scene, &
                              CLBLM_Profile
      IMPLICIT NONE

      type(CLBLM_Path)       ,intent(out) :: path !out
      type(CLBLM_Profile)    ,intent(in)  :: scnPrfl
      real                   ,intent(in)  :: refracIndex(:) !(prfl%toaLev)
      real                   ,intent(in)  :: zRT(:)
      integer                ,intent(in)  :: nRTlev
      type(CLBLM_PathGeom)   ,intent(in)  :: pathGeo
      real                   ,intent(in)  :: earthRadius
      integer                ,intent(in)  :: zeroAbsAmt


      integer :: i,j,k, M,N
      integer :: NMOL
      real    :: H1,H2,ANGLE,PHI,HMIN
      integer :: LEN, IAMT
      real    :: RANGE,BETA,BENDNG
      real    :: FAC
      integer :: IBMAX,IMMAX
      real    :: RE

      integer              :: IOUTMX
      integer              :: IPMAX
      real    ,allocatable :: ZBND(:),ZMDL(:),PM(:),TM(:),RFNDXM(:)
      real    ,allocatable :: ZOUT(:),ZPTH(:),PP(:),TP(:)
      real    ,allocatable :: SP(:),PPSUM(:),TPSUM(:),RHOPSM(:),DENM(:,:),AMTP(:,:)
      real    ,allocatable :: ALTZ(:),PZ(:),TZ(:)
      real    ,allocatable :: SOUT(:),PBAR(:),TBAR(:),RHOSUM(:),AMOUNT(:,:),WN2L(:)
      real    ,allocatable :: SECNTA(:)
      integer ,allocatable :: IPATH(:)



      H1    = pathGeo%Hobs
      H2    = pathGeo%Hend
      ANGLE = pathGeo%ANGLE
      PHI   = pathGeo%PHI
      HMIN  = pathGeo%HMIN
      LEN   = pathGeo%LEN

      IAMT = 1 !IAMT = 1: CALCULATE AMOUNTS, IAMT = 2: DO NOT CALCULATE AMOUNTS

      IBMAX  = nRTlev
      IMMAX  = scnPrfl%nLev !prfl%toaLev
      RE     = earthRadius
      NMOL   = scnPrfl%nMol !nLnMol + nXsMol

      allocate( ZBND(IBMAX) )
      ZBND(1:IBMAX) = zRT(1:IBMAX)

      allocate( ZMDL(   IMMAX) )
      allocate( PM(     IMMAX) )
      allocate( TM(     IMMAX) )
      allocate( RFNDXM( IMMAX) )
      ZMDL(  1:IMMAX) = scnPrfl%Z(1:IMMAX)
      PM(    1:IMMAX) = scnPrfl%P(1:IMMAX)
      TM(    1:IMMAX) = scnPrfl%T(1:IMMAX)
      RFNDXM(1:IMMAX) = refracIndex(1:IMMAX)


      !--- Allocate output arrays
      !
      allocate( ZOUT(         IBMAX+3) )
      allocate( ZPTH(   IMMAX+IBMAX+3) )
      allocate( PP(     IMMAX+IBMAX+3) )
      allocate( TP(     IMMAX+IBMAX+3) )
      allocate( SP(     IMMAX+IBMAX+3) )
      allocate( PPSUM(  IMMAX+IBMAX+3) )
      allocate( TPSUM(  IMMAX+IBMAX+3) )
      allocate( RHOPSM( IMMAX+IBMAX+3) )

      allocate( DENM( NMOL,IMMAX        ) ) !this is for input
      allocate( AMTP( NMOL,IMMAX+IBMAX+3) )

      DENM(1:NMOL,1:IMMAX) = scnPrfl%Q(1:nMol,1:IMMAX)

      !---
      CALL RFPATH( H1,                                                  & !inout
                   H2,                                                  & !in
                   ANGLE,                                               & !in
                   PHI,                                                 & !inout
                   LEN,                                                 & !in
                   HMIN,                                                & !in
                   IAMT,                                                & !in
                   RANGE,BETA,BENDNG,                                   & !out
                   ZBND,                                                & !in
                   ZMDL,                                                & !inout
                   PM,TM,RFNDXM,IBMAX,IMMAX, RE,                        & !in
                   ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX,  & !out
                   DENM,                                                & !in
                   AMTP,                                                & !out
                   NMOL )                                                 !in


      !--- Call FPACK
      ! * Following LBLRTM, WN2L is the layer amount excluding line molecules.
      ! ? What about return RHOSUM, the total layer air amount, instead of WN2L?
      !
      allocate( ALTZ(      0:IOUTMX-1) )
      allocate( PZ(        0:IOUTMX-1) )
      allocate( TZ(        0:IOUTMX-1) )
      allocate( SOUT(        IOUTMX-1) )
      allocate( PBAR(        IOUTMX-1) )
      allocate( TBAR(        IOUTMX-1) )
      allocate( RHOSUM(      IOUTMX-1) )
      allocate( AMOUNT( NMOL,IOUTMX-1) )
      !allocate( WN2L(        IOUTMX-1) )
      allocate( SECNTA(      IOUTMX-1) )
      allocate( IPATH(       IOUTMX-1) )

      CALL FPACK( H1,H2,LEN,zeroAbsAmt, &                     !intent(in)
                  ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM, &    !intent(in)
                  IOUTMX, &                                   !intent(inout)
                  IPMAX,&                                     !intent(in)
                  AMTP,nMol,scnPrfl%molID, &                  !intent(in)
                  ALTZ,PZ,TZ,SOUT,PBAR,TBAR,RHOSUM,AMOUNT, &  !intent(out)
                  SECNTA,IPATH )                              !intent(out)


      !--- Load the path object
      !
      call CLBLM_Path_init( path, nMol, IOUTMX-1 )

      path%molID(1:nMol)          = scnPrfl%molID(1:nMol)
      path%W( 1:nMol,1:path%nLay) = AMOUNT(1:nMol,1:path%nLay)
      path%Wtot(     1:path%nLay) = RHOSUM(1:path%nLay)  !WN2L(  1:path%nLay) + sum( AMOUNT, 1 )
      path%Pave(     1:path%nLay) = PBAR(  1:path%nLay)
      path%Tave(     1:path%nLay) = TBAR(  1:path%nLay)
      path%IPATH(    1:path%nLay) = IPATH( 1:path%nLay)
      !path%Wbroad(   1:path%nLay) = WN2L(  1:path%nLay)
      !path%SECANT(   1:path%nLay) = SECNTA(1:path%nLay)

      path%zRT( 1:path%nRTlev) = ALTZ(  0:path%nLay)
      path%pRT( 1:path%nRTlev) = PZ(    0:path%nLay)
      path%T(   1:path%nRTlev) = TZ(    0:path%nLay)

      !--- Level concentration on RTgrid.
      call densityOnPathLevel( path%Q, & ![nMol,1:nLev], output
                               path%zRT, &
                               path%nRTlev, &
                               path%molID, &
                               scnPrfl )

      !path%airMass will be calculated later when both vertical and slant
      !             path amounts are available.
      !path%isoMol%nMol,path%isoMol%molID,path%isoMol%W will be added later after.

      !TBOUND = path%T(1) will be taken care out of the CLBLM-ATM module.

      !--- Store the final path geometry.
      path%geom = pathGeo


      !--- Deallocate memory
      deallocate( ZBND,ZMDL,PM,TM,RFNDXM )
      deallocate( ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM )
      deallocate( DENM,AMTP )
      deallocate( ALTZ,PZ,TZ,SOUT,PBAR,TBAR,RHOSUM,AMOUNT, &
                  SECNTA,IPATH )

   END SUBROUTINE

   !--------------------------------------------------------------------
   !Interpolate molecular density onto RT grid (rtZ).
   !--------------------------------------------------------------------
   SUBROUTINE densityOnPathLevel( Q_rtLev, Z_rtLev, nRTLev, ajMolNames, prfl )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: molIndex
      USE Module_Utility    ,ONLY: EXPINT
      USE Module_Scene      ,ONLY: CLBLM_Profile
      IMPLICIT NONE

      real                ,intent(out) :: Q_rtLev(:,:)   ![nAJMol,nLev], density on RTgrid
      real                ,intent(in)  :: Z_rtLev(:)     ![nLev], Level altitude on RTgrid
      integer             ,intent(in)  :: nRTLev         !number of RTgrid leves
      character(*)        ,intent(in)  :: ajMolNames(:)  !molecular names for A.J. calculation
      type(CLBLM_Profile) ,intent(in)  :: prfl

      real ,PARAMETER :: TOL=5.E-4
      integer :: lp,lr,k,im
      real    :: A


      lp = 2
      DO lr = 1, nRTLev

         do ! FIND THE SMALLEST ZX GE rtZ(L) and do interpolation

            if ( Z_rtLev(lr) <= prfl%Z(lp) .OR. lp==prfl%nLev ) then !or prfl%ToaLev?

               if ( abs(Z_rtLev(lr)-prfl%Z(lp))<TOL ) then
                  A = 1.0
               elseif ( abs(Z_rtLev(lr)-prfl%Z(lp-1))<TOL ) then
                  A = 0.0
               else
                  A = ( Z_rtLev(lr) - prfl%Z(lp-1) ) / &
                      ( prfl%Z(lp)  - prfl%Z(lp-1) )
               endif
               if ( A < 0.0 .OR. A > 1.0 ) print*, ('--- densityOnPathLevel(): Extrapolating RT altitude grid ')

               do k = 1,size(ajMolNames)
                  im = molIndex( ajMolNames(k), prfl%molID )
                  call EXPINT( Q_rtLev(k,lr), prfl%Q(im,lp-1), prfl%Q(im,lp), A )
               enddo

               EXIT !do next RT level
            else
               lp = lp+1
            endif

         enddo

      ENDDO  !DO lr = 1, nRTLev

   END SUBROUTINE


!-----------------------------------------------------------------------
!   SUBROUTINE FSCGEO ()
!
!     -------------------------------------------------------------
!     This routine was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     It was also modified to eliminate GOTO statements in order to
!     make the program easier to understand.
!     These changes were obtained from H. Snell (March, 1996).
!     -------------------------------------------------------------
!
!     FSCGEO INTERPRETS THE ALLOWABLE COMBINATIONS OF INPUT PATH
!     PARAMETERS INTO THE STANDARD SET H1,H2,ANGLE,PHI,HMIN, AND LEN.
!     THE ALLOWABLE COMBINATIONS OF INPUT PARAMETERS ARE-
!
!     FOR ITYPE = 2 (SLANT PATH H1 TO H2)
!        A. H1, H2, AND ANGLE
!        B. H1, ANGLE, AND RANGE
!        C. H1, H2, AND RANGE
!        D. H1, H2, AND BETA
!     FOR ITYPE = 3 (SLANT PATH H1 TO SPACE, H2 = ZMAX(=100 KM,M=1 TO 6
!        A. H1 AND ANGLE,
!        B. H1 AND HMIN (INPUT AS H2).
!
!     THE SUBROUTINE ALSO DETECTS BAD INPUT (IMPOSSIBLE GEOMETRY) AND
!     ITYPE = 2 CASES WHICH INTERSECT THE EARTH, AND RETURNS THESE
!     CASES WITH ERROR FLAGS.
!     THE SUBROUTINE FNDHMN IS CALLED TO CALCULATE HMIN, THE MINIMUM
!     HEIGHT ALONG THE PATH, AND PHI, THE ZENITH ANGLE AT H2, USING THE
!     ATMOSPHERIC PROFILE STORED IN /MDATA/
!
!-----------------------------------------------------------------------
      SUBROUTINE FSCGEO( H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,IERROR, & !HOBS,&
                         inZBND,ZMAX,ZMIN,inZMDL,PM,TM,RFNDXM,inIBMAX,IMMAX,RE, &
                         DENM,NMOL )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: pi, DEG
      USE Module_Config       ,ONLY: IPR,noPrnt
      IMPLICIT NONE

      real     ,intent(inout) :: H1,H2
      real     ,intent(inout) :: ANGLE
      real     ,intent(inout) :: RANGE
      real     ,intent(inout) :: BETA
      integer  ,intent(in)    :: ITYPE
      integer  ,intent(inout) :: LEN
      real     ,intent(out)   :: HMIN
      real     ,intent(out)   :: PHI
      integer  ,intent(inout) :: IERROR
      !
      real     ,intent(in)    :: inZBND(:)
      real     ,intent(in)    :: ZMAX,ZMIN
      real     ,intent(in)    :: inZMDL(:)
      real     ,intent(in)    :: PM(:)
      real     ,intent(in)    :: TM(:)
      real     ,intent(in)    :: RFNDXM(:)
      integer  ,intent(in)    :: inIBMAX
      integer  ,intent(in)    :: IMMAX
      real     ,intent(in)    :: RE
      !
      real     ,intent(in)    :: DENM(:,:) !Not really used
      integer  ,intent(in)    :: NMOL


      ! ZMDL may be changed by RFPATH/AMERGE
      ! ZBND may be changed by FDBETA
      ! IBMAX may be changed by FDBETA
      real, allocatable :: ZBND(:)
      real, allocatable :: ZMDL(:)
      integer           :: IBMAX

      INTEGER ::  ISELCT,     ITER
      REAL    ::  DIFFANGLE,  ERARG2,   ERARG3,  H2ST
      REAL    ::  H_TOA,      R1,      R2
      REAL    ::  RADCONV,    SINANGLE, SINPHI,  SINTOA
      REAL    ::  SINTOA_SAT,           TOA_ANG, TOA_SAT
      REAL    ::  ZARG2,      ZARG3


      allocate(ZBND(size(inZBND)))
      allocate(ZMDL(size(inZMDL)))
      ZBND(:) = inZBND(:)
      ZMDL(:) = inZMDL(:)
      IBMAX   = inIBMAX

      ITER = 0

      ! Check for error

      IF ((ITYPE.NE.3).AND.(ITYPE.NE.2)) GOTO 90

      IF (ITYPE.EQ.3) THEN

         ! Slant path to space
         ! NOTE: If both HMIN and ANGLE are zero, then ANGLE is
         !       assumed specified

         IF (H2.EQ.0) THEN

            ! Case 3A: H1,SPACE,ANGLE

            WRITE (IPR,900)
            H2 = ZMAX
            CALL FNDHMN( H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR, &
                         ZMDL, RFNDXM, IMMAX, RE )

         ELSE

            ! Case 3B: H1,HMIN,SPACE

            WRITE (IPR,905)
            HMIN = H2
            H2 = ZMAX
            IF (H1.LT.HMIN) GO TO 80
            CALL FNDHMN( HMIN,90.0,H1,LEN,ITER,HMIN,ANGLE,IERROR, &
                         ZMDL, RFNDXM, IMMAX, RE )
            CALL FNDHMN( HMIN,90.0,H2,LEN,ITER,HMIN,PHI,IERROR, &
                         ZMDL, RFNDXM, IMMAX, RE )
            IF (HMIN.LT.H1) LEN = 1
         ENDIF
      ENDIF

      IF (ITYPE.EQ.2) THEN

         ! Assign the variable ISELCT to the following cases
         ! (depending on input parameters):
         !
         ! -----------------------------------------------
         ! H1   H2   ANGLE  RANGE  BETA  =>   CASE  ISELCT
         ! -----------------------------------------------
         ! X    X      X                       2A     21
         ! X           X      X                2B     22
         ! X    X             X                2C     23
         ! X    X                   X          2D     24
         ! -----------------------------------------------

         IF (RANGE.GT.0.0) THEN

            !Must be Case 2B or Case 2C

            IF (H2.GT.0.0) THEN

               !Case 2C

               ISELCT=23
            ELSEIF (ANGLE.EQ.0.0) THEN
               WRITE(IPR,1000)
               WRITE(*,1000)
               ISELCT=23
            ELSE

               ! Case 2B

               ISELCT=22
            ENDIF
         ELSEIF (BETA.GT.0.0) THEN

            ! Case 2D (beta cannot be zero)

            ISELCT=24
         ELSE

            !Case 2A, since RANGE and BETA are both zero

            ISELCT=21
         ENDIF

         IF (ISELCT.EQ.21) THEN

            !Case 2A: H1, H2, ANGLE

            if (noprnt .ge.0) WRITE (IPR,910)
            IF (H1.GE.H2.AND.ANGLE.LE.90.0) GO TO 110
            IF (H1.EQ.0.0.AND.ANGLE.GT.90.0) GO TO 120
            IF (H2.LT.H1.AND.ANGLE.GT.90.0) WRITE (IPR,915) LEN
            H2ST = H2
            CALL FNDHMN( H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR, &
                         ZMDL, RFNDXM, IMMAX, RE )
            IF (H2.NE.H2ST) GO TO 120
         ENDIF

         IF (ISELCT.EQ.22) THEN

             ! Case 2B: H1, ANGLE, RANGE
             ! Assume refraction

            if (noprnt .ge.0) WRITE (IPR,920)
            CALL NEWH2( H1,H2,ANGLE,RANGE,BETA,LEN,HMIN,PHI, &
                        ZMDL, RFNDXM, ZMAX, RE, IMMAX )
         ENDIF

         IF (ISELCT.EQ.23) THEN

            ! Case 2C: H1, H2, RANGE

            if (noprnt .ge.0) WRITE (IPR,930)
            IF (ABS(H1-H2).GT.RANGE) GO TO 100
            R1 = H1+RE
            R2 = H2+RE

            ZARG2 = (H1**2-H2**2+RANGE**2+2.0*RE*(H1-H2)) /(2.0*R1*     &
            RANGE)
            ERARG2 = ABS(ZARG2)-1.0
            IF ((ERARG2.LE.1.0E-6).AND.(ERARG2.GE.0.0)) THEN
               IF (ZARG2.LT.0.0) THEN
                  ZARG2 = -1.0
               ELSE
                  ZARG2 = 1.0
               ENDIF
            ENDIF
            ANGLE = 180.0-ACOS(ZARG2)*DEG
            ZARG3 = (H2**2-H1**2+RANGE**2+2*RE*(H2-H1))/(2.0*R2*RANGE)
            ERARG3 = ABS(ZARG3)-1.0
            IF ((ERARG3.LE.1.0E-6).AND.(ERARG3.GE.0.0)) THEN
               IF (ZARG3.LT.0.0) THEN
                  ZARG3 = -1.0
               ELSE
                  ZARG3 = 1.0
               ENDIF
            ENDIF
            PHI = 180.0-ACOS(ZARG3)*DEG
            BETA = PHI+ANGLE-180.

            IF (RANGE.GT.2.0.AND.BETA.GT.0) THEN
               CALL FDBETA (H1,H2,BETA,ANGLE,PHI,LEN,HMIN,IERROR, &
                            ZMDL,PM,TM,RFNDXM,ZBND,ZMAX,ZMIN,IBMAX,IMMAX,RE,&
                            DENM,NMOL )
            ELSE
               LEN = 0
               IF (ANGLE.GT.90.0.AND.PHI.GT.90.0) LEN = 1
               CALL FNDHMN( H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR, &
                            ZMDL, RFNDXM, IMMAX, RE )
            ENDIF
         ENDIF

         IF (ISELCT.EQ.24) THEN

            ! Case 2D: H1, H2, BETA

            CALL FDBETA (H1,H2,BETA,ANGLE,PHI,LEN,HMIN,IERROR, &
                         ZMDL,PM,TM,RFNDXM,ZBND,ZMAX,ZMIN,IBMAX,IMMAX,RE,&
                         DENM,NMOL )
         ENDIF
      ENDIF
      ! End of allowed cases


      ! Test IERROR and recheck LEN
      !
      IF (IERROR.NE.0) RETURN

      LEN = 0
      IF (HMIN.LT.  MIN(H1,H2)) LEN = 1

      ! Reduce path endpoints above ZMAX to ZMAX
      IF (HMIN.GE.ZMAX) GO TO 130
      IF (H1.GT.ZMAX.OR.H2.GT.ZMAX) &
         CALL REDUCE( H1,H2,ANGLE,PHI,ITER, ZMDL,RFNDXM,ZMAX,IMMAX,RE )


      ! At this point the following parameters are defined-
      !     H1,H2,ANGLE,PHI,HMIN,LEN


      ! Calculate sin(PHI) and sin(ANGLE) and output
      radconv = 2.*pi/360. !yma 151225, Pi was not defined before. Corrected.
      sinphi = sin(radconv*phi)
      sinangle = sin(radconv*angle)
      if (noprnt .ge. 0) WRITE (IPR,935)                                &
       &                   H1,H2,ANGLE,sinangle,PHI,sinphi,HMIN,LEN

      RETURN

      ! Error messages
   80 CONTINUE
      WRITE (IPR,940) H1,HMIN
      GO TO 140
   90 WRITE (IPR,945) ITYPE,ITYPE
      GO TO 140
  100 WRITE (IPR,950) H1,H2,RANGE
      GO TO 140
  110 CONTINUE
      WRITE (IPR,955) H1,H2,ANGLE
      GO TO 140
  120 WRITE (IPR,960)
      GO TO 140
  130 WRITE (IPR,965) ZMAX,H1,H2,HMIN
  140 IERROR = 1

      RETURN

  900 FORMAT (//,' CASE 3A: GIVEN H1,H2=SPACE,ANGLE')
  905 FORMAT (//,' CASE 3B: GIVEN H1, HMIN, H2=SPACE')
  910 FORMAT (//,' CASE 2A: GIVEN H1, H2, ANGLE')
  915 FORMAT (//,' EITHER A SHORT PATH (LEN=0) OR A LONG PATH ',        &
     &        'THROUGH A TANGENT HEIGHT (LEN=1) IS POSSIBLE: LEN = ',   &
     &        I3)
  920 FORMAT (//,' CASE 2B:, GIVEN H1, ANGLE, RANGE',//,10X,            &
     &        'NOTE: H2 IS COMPUTED FROM H1, ANGLE, AND RANGE ',        &
     &        'ASSUMING REFRACTION')
  925 FORMAT (//,10X,'CALCULATED H2 IS LESS THAN ZERO:',/,10X,          &
     &        'RESET H2 = 0.0 AND RANGE = ',F10.3)
  930 FORMAT (//,' CASE 2C: GIVEN H1, H2, RANGE',//,10X,                &
     &        'NOTE: ANGLE IS COMPUTED FROM H1, H2, AND RANGE ',        &
     &        'ASSUMING NO REFRACTION')
  935 FORMAT (///,' SLANT PATH PARAMETERS IN STANDARD FORM',/           &
     &        /,10X,'H1         = ',F12.6,' KM',                        &
     &        /,10X,'H2         = ',F12.6,' KM',                        &
     &        /,10X,'ANGLE      = ',F12.6,' DEG',                       &
     &        /,10X,'sin(ANGLE) = ',F12.6,                              &
     &        /,10X,'PHI        = ',F12.6,' DEG',                       &
     &        /,10X,'sin(PHI)   = ',F12.6,                              &
     &        /,10X,'HMIN       = ',F12.6,' KM',                        &
     &        /,10X,'LEN        = ',I10)
  937 FORMAT (///,' SLANT PATH PARAMETERS AT SATELLITE',/               &
     &        /,10X,'H_SAT        = ',F12.6,' KM',                      &
     &        /,10X,'PHI_SAT      = ',F12.6,' DEG'                      &
     &        /,10X,'sin(PHI_SAT) = ',F12.6,                            &
     &        /,10X,'PHI_SAT-PHI  = ',F12.6,' DEG')
  940 FORMAT ('0FSCGEO: CASE 3B (H1,HMIN,SPACE): ERROR IN INPUT DATA',  &
     &        //,10X,'H1 = ',F12.6,'    IS LESS THAN HMIN = ',F12.6)
  945 FORMAT ('0FSCGEO: ERROR IN INPUT DATA, ITYPE NOT EQUAL TO ',      &
     &        ' 2, OR 3.   ITYPE = ',I10,E23.14)
  950 FORMAT ('0FSCGEO: CASE 2C (H1,H2,RANGE): ERROR IN INPUT DATA',    &
     &        //,10X,'ABS(H1-H2) GT RANGE;  H1 = ',F12.6,'    H2 = ',   &
     &        F12.6,'    RANGE = ',F12.6)
  955 FORMAT ('0FSCGEO: CASE 2A (H1,H2,ANGLE): ERROR IN INPUT DATA',    &
     &        //,10X,'H1 = ',F12.6,'    IS GREATER THAN OR EQUAL TO',   &
     &        ' H2 = ',F12.6,/,10X,'AND ANGLE = ',F12.6,'    IS LESS',  &
     &        ' THAN OR EQUAL TO 90.0')
  960 FORMAT ('0FSCGEO: ITYPE = 2: SLANT PATH INTERSECTS THE EARTH',    &
     &        ' AND CANNOT REACH H2')
  965 FORMAT (' FSCGEO:  THE ENTIRE PATH LIES ABOVE THE TOP ZMAX ',     &
     &        'OF THE ATMOSPHERIC PROFILE',//,10X,'ZMAX = ',G13.6,5X,   &
     &        '  H1 = ',G13.6,5X,'  H2 = ',G13.6,'  HMIN = ',G13.6)

 1000 FORMAT (/3X, 'Ambiguous Inputs:',/3X,'H1 and RANGE are both > 0', &
     &    /3X,'but H2 and ANGLE = 0',//3X,'Path could be 2B or 2C',     &
     &    //5X,'will assume 2C',//3X,                                   &
     &    'change in FSCGEO if 2B is desired')

     END SUBROUTINE


!----------------------------------------------------------------------
!     SUBROUTINE FNDHMN ()
!
!     THIS SUBROUTINE CALCULATES THE MINIMUM ALTITUDE HMIN ALONG
!     THE REFRACTED PATH AND THE FINAL ZENITH ANGLE PHI.
!     THE PARAMETER LEN INDICATES WHETHER THE PATH GOES THROUGH
!     A TANGENT HEIGHT (LEN=1) OR NOT (LEN=0).  IF ANGLE > 90 AND
!     H1 > H2, THEN LEN CAN EITHER BE 1 OR 0, AND THE CHOICE IS
!     LEFT TO THE USER.
!     THE (INDEX OF REFRACTION - 1.0) IS MODELED AS AN EXPONENTIAL
!     BETWEEN THE LAYER BOUNDARIES, WITH A SCALE HEIGHT SH AND AN
!     AMOUNT AT THE GROUND GAMMA.
!     CPATH IS THE REFRACTIVE CONSTANT FOR THIS PATH AND
!     EQUALS  INDEX(H1)*(RE+H1)*SIN(ANGLE).
!----------------------------------------------------------------------
      SUBROUTINE FNDHMN( H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR, &
                         ZMDL, RFNDXM, IMMAX, RE )
!----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, DEG
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      real    ,intent(in)    :: H1
      real    ,intent(in)    :: ANGLE
      real    ,intent(inout) :: H2
      integer ,intent(inout) :: LEN
      integer ,intent(in)    :: ITER
      real    ,intent(out)   :: HMIN
      real    ,intent(out)   :: PHI
      integer ,intent(inout) :: IERROR
      real    ,intent(in)    :: ZMDL(:)
      real    ,intent(in)    :: RFNDXM(:)
      integer ,intent(in)    :: IMMAX
      real    ,intent(in)    :: RE


      ! ETA MAY BE TOO SMALL FOR SOME COMPUTERS. TRY 1.0E-7 FOR 32 BIT
      ! WORD MACHINES
      real ,PARAMETER :: DH=0.2
      real ,PARAMETER :: ETA=5.0E-7

      INTEGER  :: N
      REAL(r8) :: CH2,  CMIN,  CPATH, CT1
      REAL(r8) :: CTP
      REAL     :: DC,   DERIV
      REAL(r8) :: GAMMA
      REAL     :: H,    HT1,   HTP
      REAL(r8) :: SH


      N = 0
      CALL FNDSHD( H1,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CPATH = CRFRCT(H1)*SIN(ANGLE/DEG)
      CALL FNDSHD( H2,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CH2 = CRFRCT(H2)
      IF (ABS(CPATH/CH2).GT.1.0) GO TO 70
      IF (ANGLE.LE.90.0) THEN
         LEN = 0
         HMIN = H1
         GO TO 60
      ENDIF
      IF (H1.LE.H2) LEN = 1
      IF (LEN.NE.1) THEN
         LEN = 0
         HMIN = H2
         GO TO 60
      ENDIF

      ! LONG PATH THROUGH A TANGENT HEIGHT.
      ! SOLVE ITERATIVELY FOR THE TANGENT HEIGHT HT.
      ! HT IS THE HEIGHT FOR WHICH  INDEX(HT)*(RE+HT) = CPATH.
      !
      CALL FNDSHD( 0.0,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CMIN = CRFRCT(0.0)

      !FOR BETA CASES (ITER>0), ALLOW FOR HT < 0.0
      IF (ITER.EQ.0.AND.CPATH.LT.CMIN) GO TO 50
      HT1 = H1*SIN(ANGLE/DEG)+(SIN(ANGLE/DEG)-1.0)*RE

      !ITERATE TO FIND HT
   30 CONTINUE
      N = N+1
      CALL FNDSHD( HT1,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CT1 = CRFRCT(HT1)
      IF (ABS((CPATH-CT1)/CPATH).LT.ETA) GO TO 40
      IF (N.GT.15) GO TO 80
      HTP = HT1-DH
      CALL FNDSHD( HTP,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CTP = CRFRCT(HTP)
      DERIV=(CT1-CTP)/DH
      HT1=HT1+(CPATH-CT1)/DERIV
      GO TO 30
   40 CONTINUE
      HMIN=HT1
      GO TO 60
   50 CONTINUE

      ! TANGENT PATH INTERSECTS EARTH
      H2 = 0.0
      HMIN = 0.0
      LEN = 0
      CH2 = CMIN
      WRITE (IPR,900) H1,ANGLE
   60 CONTINUE

      !CALCULATE THE ZENITH ANGLE PHI AT H2
      PHI = ASIN(CPATH/CH2)*DEG
      IF (ANGLE.LE.90.0.OR.LEN.EQ.1) PHI = 180.0-PHI

      RETURN

      ! H2 LT TANGENT HEIGHT FOR THIS H1 AND ANGLE
   70 CONTINUE
      WRITE (IPR,905)
      IERROR = 2

      RETURN

   80 CONTINUE
      DC = CPATH-CT1
      WRITE (IPR,910) N,CPATH,CT1,DC,HT1

      STOP ' FNDHMN '

  900 FORMAT (///,' TANGENT PATH WITH H1 = ',F10.3,' AND ANGLE = ',     &
     &        F10.3,' INTERSECTS THE EARTH',//,10X,                     &
     &        'H2 HAS BEEN RESET TO 0.0 AND LEN TO 0')
  905 FORMAT ('0H2 IS LESS THAN THE TANGENT HEIGHT FOR THIS PATH ',     &
     &        'AND CANNOT BE REACHED')
  910 FORMAT (///,'0FROM SUBROUTINE FNDHMN :',//,10X,                   &
     &        'THE PROCEEDURE TO FIND THE TANGENT HEIGHT DID NOT ',     &
     &        'CONVERG AFTER ',I3,'  ITERATIONS',//,10X,'CPATH   = ',   &
     &        F12.5,' KM',//,10X,'CT1     = ',F12.5,' KM',//,10X,       &
     &        'DC      = ',E12.3,' KM',//,10X,'HT1     = ',F12.5,' KM')

      CONTAINS !---------------- Internal function -------------------

      FUNCTION CRFRCT(H) !H is local. RE,SH and GAMMA are visible from host.
         REAL(r8)         :: CRFRCT
         real ,intent(in) :: H

         CRFRCT = (RE+H)*ANDEXD(H,SH,GAMMA) !ANDEXD() is a module procedure
      END FUNCTION !internal function

      END  SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE NEWH2()
!
!     Changed for LBLRTM to correct geometry problems
!
!     THIS ROUTINE DETERMINES H2,BETA, TANGENT HEIGHT AND LEN.
!     ADOPTED FROM THE MODTRAN2 GEOMETRY PACKAGE
!
!     INPUTS ARE: H1, ZENTIH ANGLE (ANGLE) AND RANGE.
!     LEN = 1 IF THE PATH GOES THROUGH HTAN.
!
!     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO FASE01
!     MXLAY IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
!     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
!         OBTAINED BY MERGING ZMDL AND ZOUT
!     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
!-----------------------------------------------------------------------
      SUBROUTINE NEWH2( H1,H2,ANGLE,RANGE,BETA,LEN,HTAN,PHI, &
                        ZMDL, RFNDXM, ZMAX, RE, IMMAX)
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, DEG
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      real     ,intent(in)    :: H1
      real     ,intent(inout) :: H2
      real     ,intent(in)    :: ANGLE
      real     ,intent(inout) :: RANGE
      real     ,intent(inout) :: BETA
      integer  ,intent(inout) :: LEN
      real     ,intent(out)   :: HTAN
      real     ,intent(out)   :: PHI
      real     ,intent(in)    :: ZMDL(:)
      real     ,intent(in)    :: RFNDXM(:)
      real     ,intent(in)    :: ZMAX
      real     ,intent(in)    :: RE
      integer  ,intent(in)    :: IMMAX


      INTEGER  :: J,     JMAX
      REAL(r8) :: CPATH,  CPJ,   CPJ1, GAMMA
      REAL     :: H
      REAL(r8) :: RE2,    SH
      REAL     :: ZJ,     ZJ1


      RE2=RE
      ! COMPUTE CPATH OR PATH CONSTANT
      CALL FNDSHD( H1,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CPATH = CRFRCT(H1)*SIN(ANGLE/DEG)

      ! ANGLE = 90 at H1 implies that H1 = tangent height

      IF (ANGLE.EQ.90.0) THEN
         HTAN=H1
      ELSE
         DO 100 J=1,IMMAX
            IF (H1.GE.ZMDL(J)) JMAX=J
  100    CONTINUE
         JMAX=JMAX+1
         ZJ1=ZMDL(JMAX)
         CPJ1=CRFRCT(ZJ1)
         HTAN=-1.0
         DO 200 J=JMAX,1,-1
            IF (HTAN.LT.0.0) THEN
               IF (J.EQ.1) THEN
                  HTAN=0.0
               ELSE
                  CPJ=CPJ1
                  ZJ=ZJ1
                  ZJ1=ZMDL(J-1)
                  CPJ1=CRFRCT(ZJ1)
                  IF ((CPATH.LE.CPJ).AND.(CPATH.GE.CPJ1)) THEN
                     HTAN=RTBIS( ZJ1,CPJ1,ZJ,CPJ,CPATH, ZMDL,RFNDXM,RE,IMMAX )
                  ENDIF
               ENDIF
            ENDIF
  200    CONTINUE
      ENDIF

      ! Find H2, BETA AND LEN
      !
      CALL FNDPTH( CPATH,H1,HTAN,H2,RANGE,BETA,LEN,ANGLE,PHI, &
                   ZMDL, RFNDXM, DEG, RE, ZMAX, IMMAX )


      ! Ensure LEN is not reset in FSCGEO if direct path
      IF (LEN.EQ.0) HTAN=H2

      ! IF (ANGLE .LE. 90.0) HTAN CARRIES HMIN NOT HTAN
      IF (ANGLE .LE. 90.0) HTAN = MIN(H1,H2)

      RETURN

      CONTAINS !----------------- Internal function ----------------

      FUNCTION CRFRCT(H)  !H is local. RE2,SH and GAMMA are visible from host.
         REAL(r8)         :: CRFRCT
         real ,intent(in) :: H

         CRFRCT =(RE2+H)*ANDEXD(H,SH,GAMMA) !ANDEXD() is a module procedure
      END FUNCTION

      END  SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE FDBETA ()
!
!     GIVEN H1,H2,AND BETA (THE EARTH CENTERED ANGLE) THIS SUBROUTINE
!     CALCULATES THE INITIAL ZENITH ANGLE AT H1 THROUGH AN ITERATIVE
!     PROCEDURE
!-----------------------------------------------------------------------
      SUBROUTINE FDBETA( H1,H2,BETAS,ANGLE,PHI,LEN,HMIN,IERROR, &
                         ZMDL,PM,TM,RFNDXM,ZBND,ZMAX,ZMIN,IBMAX,IMMAX,RE,&
                         DENM,NMOL )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, DEG
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      real     ,intent(in)    :: H1,H2
      real     ,intent(inout) :: BETAS
      real     ,intent(inout) :: ANGLE
      real     ,intent(out)   :: PHI
      integer  ,intent(inout) :: LEN
      real     ,intent(out)   :: HMIN
      integer  ,intent(inout) :: IERROR
      !
      real     ,intent(inout) :: ZMDL(:)
      real     ,intent(in)    :: PM(:)
      real     ,intent(in)    :: TM(:)
      real     ,intent(in)    :: RFNDXM(:)
      real     ,intent(inout) :: ZBND(:)
      real     ,intent(in)    :: ZMAX,ZMIN
      integer  ,intent(inout) :: IBMAX
      integer  ,intent(in)    :: IMMAX
      real     ,intent(in)    :: RE
      !
      real     ,intent(in)    :: DENM(:,:) !not really used for FDBETA
      integer  ,intent(in)    :: NMOL


      !--- Local variables
      !
      real ,PARAMETER :: TOLRNC=5.0E-3
      real ,PARAMETER :: ITERMX=10
      real ,PARAMETER :: BETD=0.04
      real ,PARAMETER :: ZER=0.


      INTEGER  :: IAMTB,  IBMSAV, IFLAG,  IORDER, ITER
      REAL     :: ANGLEP, ANGLS1, ANGLS2, BENDNG
      REAL     :: BETA1,  BETA2,  BETAP
      REAL     :: DANG, DC, DERIV
      REAL     :: HA, HB, HMING
      REAL     :: RANGE
      REAL     :: TEMP

      REAL(r8) :: RA,RB,SG,ANGLE1,ANGLE2,BETA,DBETA


      !--- These variables are really used for FDBETA.
      integer  :: IOUTMX
      integer  :: IPMAX
      real    ,allocatable :: ZOUT(:),ZPTH(:),PP(:),TP(:)
      real    ,allocatable :: SP(:),PPSUM(:),TPSUM(:),RHOPSM(:)
      real    ,allocatable :: AMTP(:,:)


      allocate( ZOUT(            IBMAX+3) )
      allocate( ZPTH(      IMMAX+IBMAX+3) )
      allocate( PP(        IMMAX+IBMAX+3) )
      allocate( TP(        IMMAX+IBMAX+3) )
      allocate( SP(        IMMAX+IBMAX+3) )
      allocate( PPSUM(     IMMAX+IBMAX+3) )
      allocate( TPSUM(     IMMAX+IBMAX+3) )
      allocate( RHOPSM(    IMMAX+IBMAX+3) )
      allocate( AMTP( NMOL,IMMAX+IBMAX+3) )


      BETA = BETAS
      IFLAG = 0
      IF (H1.LE.H2) THEN
         IORDER = 1
         HA = H1
         HB = H2
      ELSE
         IORDER = -1
         HA = H2
         HB = H1
      ENDIF

      !IF AUTOLAYERING SELECTED(IBMAX = 0) THEN SET UP DUMMY
      !LBLRTM OUTPUT LAYERS
      !
      IBMSAV = IBMAX
      IF (IBMAX.EQ.0) THEN
         IBMAX = 2
         ZBND(1) = ZMIN
         ZBND(2) = ZMAX
      ENDIF

      !SET PARAMETER TO SUPRESS CALCULATION OF AMOUNTS

      IAMTB = 2

      !GUESS AT ANGLE, INTEGRATE TO FIND BETA, TEST FOR
      !CONVERGENCE, AND ITERATE
      !FIRST GUESS AT ANGLE: USE THE GEOMETRIC SOLUTION (NO REFRACTION)

      WRITE (IPR,900)
      ITER = 0
      RA = RE+HA
      RB = RE+HB
      SG = SQRT((HA-HB)**2+4.0*RA*RB*(SIN(BETA/(2.0*DEG)))**2)
      ANGLE1 = 180.0-ACOS((HA**2-HB**2+2.0*RE*(HA-HB)+SG**2)            &
     &         /(2.0*RA*SG))*DEG
      HMIN = HA
      IF (ANGLE1.GT.90.0) HMIN = RA*SIN(ANGLE1/DEG)-RE
      HMING = HMIN
      ANGLS1 = ANGLE1
      CALL FNDHMN( HA,ANGLS1,HB,LEN,ITER,HMIN,PHI,IERROR, &
                   ZMDL, RFNDXM, IMMAX, RE )
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH( HA,HB,ANGLS1,PHI,LEN,HMIN,IAMTB,RANGE,BETA1,BENDNG, &
                   ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, RE, &
                   ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX, &
                   DENM,AMTP,NMOL )
      WRITE (IPR,905) ITER,ANGLS1,BETA,ZER,SG,HMING,ZER,ZER

      !OBTAIN DERIVATIVE

      SG = SQRT((HA-HB)**2+4.0*RA*RB*(SIN((BETA+BETD)/(2.0*DEG)))**2)
      ANGLEP = 180.0-ACOS((HA**2-HB**2+2.0*RE*(HA-HB)+SG**2)            &
     &         /(2.0*RA*SG))*DEG
      DANG = ANGLE1-ANGLEP
      IF (HMIN.LT.0.0) THEN
         IFLAG = 1
         HMIN = 0.0
         CALL FNDHMN( HMIN,90.0,HA,LEN,ITER,HMIN,ANGLS1,IERROR, &
                      ZMDL, RFNDXM, IMMAX, RE )
      ENDIF
      ITER = 1
      LEN = 0
      IF (ANGLE1.GT.90.0) LEN = 1
      CALL FNDHMN( HA,ANGLS1,HB,LEN,ITER,HMIN,PHI,IERROR, &
                   ZMDL, RFNDXM, IMMAX, RE )
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH( HA,HB,ANGLS1,PHI,LEN,HMIN,IAMTB,RANGE,BETA1,BENDNG, &
                   ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, RE, &
                   ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX, &
                   DENM,AMTP,NMOL )
      DBETA = BETA-BETA1
      WRITE (IPR,905) ITER,ANGLS1,BETA1,DBETA,RANGE,HMIN,PHI,BENDNG
      IF (IFLAG.EQ.1.AND.BETA1.LT.BETA) GO TO 90
   50 CONTINUE
      ANGLEP = ANGLE1-DANG
      LEN = 0
      IF (ANGLEP.GT.90.0) LEN = 1
      CALL FNDHMN( HA,ANGLEP,HB,LEN,ITER,HMIN,PHI,IERROR, &
                   ZMDL, RFNDXM, IMMAX, RE )
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH( HA,HB,ANGLEP,PHI,LEN,HMIN,IAMTB,RANGE,BETAP,BENDNG, &
                   ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, RE, &
                   ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX, &
                   DENM,AMTP,NMOL )
      IF (ABS(BETA1-BETAP).LT.TOLRNC) GO TO 60
      ITER = ITER+1
      DC = BETAP-BETA1
      DERIV = -DC/BETD
      ANGLE2 = ANGLE1+(ANGLE1-ANGLEP)*(BETA-BETA1)/(BETA1-BETAP)
      ANGLS2 = ANGLE2
      LEN = 0
      IF (ANGLE2.GT.90.0) LEN = 1
      CALL FNDHMN( HA,ANGLS2,HB,LEN,ITER,HMIN,PHI,IERROR, &
                   ZMDL, RFNDXM, IMMAX, RE )
      LEN = 0
      IF (HMIN.LT.HA) LEN = 1
      CALL RFPATH( HA,HB,ANGLS2,PHI,LEN,HMIN,IAMTB,RANGE,BETA2,BENDNG, &
                   ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, RE, &
                   ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX, &
                   DENM,AMTP,NMOL )
      DBETA = BETA-BETA2
      WRITE (IPR,905) ITER,ANGLS2,BETA2,DBETA,RANGE,HMIN,PHI,BENDNG
      IF (BETA2.LT.BETA.AND.HMIN.LT.0.0) GO TO 90
      ANGLE1 = ANGLE2
      ANGLS1 = ANGLE1
      BETA1 = BETA2
      IF (ABS(BETA-BETA2).LT.TOLRNC) GO TO 70
      IF (ITER.GT.ITERMX) GO TO 100
      GO TO 50
   60 ANGLE2 = ANGLEP
      ANGLS2 = ANGLE2
      BETA = BETAP
   70 CONTINUE
      IF (HMIN.LT.0.0) GO TO 90

      !CONVERGED TO A SOLUTION

      ANGLE = ANGLE2
      BETA = BETA2

      !ASSIGN ANGLE AND PHI TO PROPER H1 AND H2

      IF (IORDER.NE.1) THEN
         TEMP = PHI
         PHI = ANGLE
         ANGLE = TEMP
      ENDIF
      IBMAX = IBMSAV
      BETAS = BETA

      RETURN

      !ERROR MESSAGES

   90 CONTINUE
      WRITE (IPR,910)
      GO TO 110
  100 CONTINUE
      WRITE (IPR,915) H1,H2,BETA,ITER,ANGLE1,BETA1,ANGLE2,BETA2

  110 IERROR = 1


      deallocate( ZOUT,ZPTH,PP,TP )
      deallocate( SP,PPSUM,TPSUM,RHOPSM )
      deallocate( AMTP )

      RETURN

  900 FORMAT (///,' CASE 2D: GIVEN H1, H2,  BETA:',//,                  &
     &        ' ITERATE AROUND ANGLE UNTIL BETA CONVERGES',//,          &
     &        ' ITER    ANGLE',T21,'BETA',T30,'DBETA',T40,'RANGE',      &
     &        T51,'HMIN',T61,'PHI',T70,'BENDING',/,T10,'(DEG)',T21,     &
     &        '(DEG)',T30,'(DEG)',T41,'(KM)',T51,'(KM)',T60,'(DEG)',    &
     &        T71,'(DEG)',/)
  905 FORMAT (I5,3F10.4,2F10.3,2F10.4)
  910 FORMAT ('0FDBETA, CASE 2D(H1,H2,BETA): REFRACTED TANGENT ',       &
     &        'HEIGHT IS LESS THAN ZERO-PATH INTERSECTS THE EARTH',     &
     &        //,10X,'BETA IS TOO LARGE FOR THIS H1 AND H2')
  915 FORMAT ('0FDBETA, CASE 2D (H1,H2,BETA): SOLUTION DID NOT ',       &
     &        ' CONVERGE',//,10X,'H1 = ',F12.6,'    H2 = ',F12.6,       &
     &        '    BETA = ',F12.6,'    ITERATIONS = ',I4,//,10X,        &
     &        'LAST THREE ITERATIONS ',//,(10X,'ANGLE = ',F15.9,        &
     &        '    BETA = ',F15.9))

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE REDUCE ()
!
!     ZMAX IS THE HIGHEST LEVEL IN THE ATMOSPHERIC PROFILE STORED IN
!     COMMON /MDATA/.  IF H1 AND/OR H2 ARE GREATER THAN ZMAX, THIS
!     SUBROUTINE REDUCES THEM TO ZMAX AND RESETS ANGLE AND/OR PHI
!     AS NECESSARY. THIS REDUCTION IS NECESSARY,FOR EXAMPLE FOR
!     SATELLITE ALTITUDES, BECAUSE (1) THE DENSITY PROFILES ARE
!     POORLY DEFINED ABOVE ZMAX AND (2) THE CALCULATION TIME FOR
!     PATHS ABOVE ZMAX CAN BE EXCESSIVE ( EG. FOR GEOSYNCRONOUS
!     ALTITUDES)
!-----------------------------------------------------------------------
      SUBROUTINE REDUCE( H1,H2,ANGLE,PHI,ITER, &
                         ZMDL,RFNDXM,ZMAX,IMMAX,RE )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: DEG
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      real    ,intent(inout) :: H1
      real    ,intent(inout) :: H2
      real    ,intent(inout) :: ANGLE
      real    ,intent(inout) :: PHI
      integer ,intent(in)    :: ITER
      !
      real    ,intent(in)    :: ZMDL(:)
      real    ,intent(in)    :: RFNDXM(:)
      real    ,intent(in)    :: ZMAX
      integer ,intent(in)    :: IMMAX
      real    ,intent(in)    :: RE


      real :: CPATH, CZMAX, ANGMAX
      real :: SH, GAMMA



      IF (H1.LE.ZMAX.AND.H2.LE.ZMAX) RETURN

      CALL FINDSH( H1,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CPATH = ANDEX(H1,SH,GAMMA)*(RE+H1)*SIN(ANGLE/DEG)

      CALL FINDSH( ZMAX,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
      CZMAX = ANDEX(ZMAX,SH,GAMMA)*(RE+ZMAX)

      ANGMAX = 180.0-ASIN(CPATH/CZMAX)*DEG

      IF (H1.LE.ZMAX) GO TO 10
      H1 = ZMAX
      ANGLE = ANGMAX
   10 CONTINUE
      IF (H2.LE.ZMAX) GO TO 20
      H2 = ZMAX
      PHI = ANGMAX
   20 CONTINUE

      IF (ITER.EQ.0) WRITE (IPR,900) ZMAX,ANGMAX

      RETURN

  900 FORMAT (///,' FROM SUBROUTINE REDUCE : ',/,10X,'ONE OR BOTH OF',  &
     &        ' H1 AND H2 ARE ABOVE THE TOP OF THE ATMOSPHERIC ',       &
     &        'PROFILE ZMAX = ',F10.3,'  AND HAVE BEEN RESET TO ZMAX.', &
     &        /,10X,'ANGLE AND/OR PHI HAVE ALSO BEEN RESET TO THE ',    &
     &        'ZENITH ANGLE AT ZMAX = ',F10.3,' DEG')

      END  SUBROUTINE


!-----------------------------------------------------------------------
!     SUBROUTINE FNDPTH()
!
!     THIS ROUTINE DETERMINES H2, BETA AND LEN.
!     INPUTS ARE H1, HTAN (TANGENT HEIGHT), RANGE (RANGEI) AND
!     THE PATH CONSTANT, CPATH.
!     RANGEO IS THE OUTPUT RANGE WHICH SHOULD EQUAL THE INPUT RANGE.
!-----------------------------------------------------------------------
      SUBROUTINE FNDPTH( CPATH,H1,HTAN,H2,RANGEI,BETA,LEN,ANGLE,PHI, &
                         ZMDL, RFNDXM, DEG, RE, ZMAX, IMMAX)
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real(r8) ,intent(in)    :: CPATH
      real     ,intent(in)    :: H1
      real     ,intent(in)    :: HTAN
      real     ,intent(inout) :: H2
      real     ,intent(inout) :: RANGEI
      real     ,intent(inout) :: BETA
      integer  ,intent(inout) :: LEN
      real     ,intent(in)    :: ANGLE
      real     ,intent(out)   :: PHI
      real     ,intent(in)    :: ZMDL(:)
      real     ,intent(in)    :: RFNDXM(:)
      real     ,intent(in)    :: DEG
      real     ,intent(in)    :: RE
      real     ,intent(in)    :: ZMAX
      integer  ,intent(in)    :: IMMAX


      real ,PARAMETER :: DR=0.005

      INTEGER  :: I
      REAL     :: BASE
      REAL(r8) :: CAPRJ,  CTHET1, CTHETA, DBETA
      REAL(r8) :: DIFF
      REAL(r8) :: DRNG,   DX
      REAL     :: DZ
      REAL(r8) :: GAMMA
      REAL     :: PERP
      REAL(r8) :: PNTGRN, R
      REAL     :: R1,     R2,     RANGEO
      REAL(r8) :: RATIO,  RPLDR,  RX,     SAVE
      REAL(r8) :: SH,     STHETA
      REAL     :: Z,      Z2



      IF (RANGEI .LT. DR) STOP 'STOPPED IN FNDPTH'

      !(RANGEI .LT. DR) SHOULD NOT HAPPEN; SO THIS CHECK IS REDUNDANT.
      RANGEO = 0
      !BETA = 0  !v12.7 removed this
      DO 200 I = 1, 2

         IF (ANGLE .LE. 90.0000 .AND. I .EQ. 1) GO TO 200

         !IF (ANGLE .LE. 90.0000) THE PATH DOES NOT GO THROUGH HTAN.
         !IF (ANGLE .LE. 90.0000) THE I = 1 CALCULATION SHOULD NOT BE DON
         !IF (ANGLE .LE. 90.0000) FOR I = 2, R1 = H1

         IF (I .EQ. 1) THEN
            R1 = H1
            R2 = HTAN
         ELSEIF (I .EQ. 2) THEN
            IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) GO TO 200

            !IF (HTAN APPROXIMATELY 0) THEN YOU ARE ABOUT TO HIT THE EART

            R2 = ZMAX
            IF (ANGLE .LE. 90.0000) THEN
               R1 = H1
            ELSE
               R1 =HTAN
            ENDIF
         ENDIF
         IF (R2 .LT. R1) THEN
            DZ = -DR
         ELSE
            DZ = DR
         ENDIF


         z = r1
         DO 100 while (z.lt.r2)
         Z2=Z
         R=Z+RE
         CALL FNDSHD( Z2,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
         RX=ANDEXD(Z2,SH,GAMMA)
         STHETA = CPATH/(RX*R)
         IF (STHETA .GT. 1.0) STHETA = 1.
         IF (STHETA .LT.-1.0) STHETA =-1.
         SAVE = STHETA
         CTHETA = SQRT(1.0-STHETA**2)
         IF (R1 .GT. R2) CTHETA = -CTHETA

         !IF (R1 .GT. R2) THEN CTHETA IS NEGATIVE BECAUSE THETA .GT. 9

         RATIO=-(RX*SH)/(RX-1.0)
         CAPRJ = -R/RATIO
         PNTGRN = 1.0/(1.0-CAPRJ*STHETA*STHETA)
         RPLDR = R+DZ
         Z2 = Z+DZ
         CALL FNDSHD( Z2,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
         RX=ANDEXD(Z2,SH,GAMMA)
         STHETA = CPATH/(RX*RPLDR)
         CTHET1 = CTHETA
         CTHETA = SQRT(1.0-STHETA**2)
         IF (R1 .GT. R2) CTHETA = -CTHETA
         DX=CTHETA*DZ+(CTHETA-CTHET1)*R
         DRNG = PNTGRN*DX
         RANGEO = RANGEO + DRNG

         DBETA = (((SAVE+STHETA)*0.5) * (PNTGRN*DX)) / (Z-0.5*DZ+RE)
         BETA = BETA+DBETA
         IF (RANGEO .GE. RANGEI) THEN
            DIFF = (RANGEI-(RANGEO-DRNG))
            H2 = Z + (DZ/DRNG)*DIFF
            BETA = BETA*DEG
            IF (I .EQ. 2) THEN
               LEN = 1
               IF (ANGLE .LE. 90.0000) LEN = 0
               IF (H2 .LT. HTAN) THEN

                  !  THIS WILL BE THE CASE IF I = 2, AND YOU HAVE
                  !  GONE THROUGH THE R-LOOP BARELY (ONLY) ONCE.

                  H2 = HTAN
                  LEN = 0
               ENDIF
            ELSE
               LEN = 0
            ENDIF

            !CORRECTION FOR VERY SHORT PATHS; HERE IT IS ABOUT 5 KM

            IF (RANGEI .LT. 5.0 .AND. RANGEO/RANGEI .GT. 1.05) THEN

               !CALCULATE BETA BY STARIGHT LINE GEOMETRY.

               PERP = SIN(ANGLE/DEG)*RANGEI
               BASE = COS(ANGLE/DEG)*RANGEI + RE+H1
               BETA = ATAN(PERP/BASE)*DEG
               RANGEO = RANGEI

               !H2 = BASE - RE

               H2 = COS(ANGLE/DEG)*RANGEI+H1
            ENDIF
            PHI = 180.0 - ACOS(CTHETA)*DEG
            RETURN
         ENDIF
         z=z+dz
  100    CONTINUE
  200 END DO

      ! COMES HERE IF YOU HAVE REACHED ZMAX, BUT YOUR RANGEI IS STILL
      ! NOT EQUAL TO OUTPUT VALUE.
      ! IN THIS CASE DO THE FOLLOWING.

      RANGEI = RANGEO
      H2 = ZMAX
      IF (ANGLE .LE. 90) THEN
         LEN = 0
      ELSE
         LEN = 1
      ENDIF
      IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) THEN

         !YOU HAVE HIT THE EARTH IF YOU ARE AT THIS POINT OF THE CODE

         LEN = 0
         H2 = 0
      ENDIF
      BETA = BETA*DEG
      PHI = 180.0 - ACOS(CTHETA)*DEG

      RETURN
      END  SUBROUTINE

!-----------------------------------------------------------------------
!     THIS FUNCTION FINDS THE ROOT OF
!            FUNC(X) = X*REFRACTIVE INDEX - CPA
!
!     THE ROOT IS ACTUALLY THE TANGENT HEIGHT, BETWEEN X1 AND X2.
!     THIS ROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS, ET AL.
!-----------------------------------------------------------------------
      FUNCTION RTBIS( X1,CX1,X2,CX2,CPATH, ZMDL,RFNDXM,RE,IMMAX )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real                 :: RTBIS
      real     ,intent(in) :: X1,X2
      real(r8) ,intent(in) :: CX1,CX2
      real(r8) ,intent(in) :: CPATH
      real     ,intent(in) :: ZMDL(:)
      real     ,intent(in) :: RFNDXM(:)
      real     ,intent(in) :: RE
      integer  ,intent(in) :: IMMAX


      integer ,PARAMETER :: JMAX=40
      real    ,PARAMETER :: XACC=1E-5

      integer  :: J
      real(r8) :: F, FMID, SH, GAMMA
      real     :: DX, XMID


      FMID=CX2-CPATH
      F=CX1-CPATH
      IF(F*FMID.GE.0.) STOP 'ROOT MUST BE BRACKETED FOR BISECTION.'
      IF(F.LT.0.)THEN
         RTBIS=X1
         DX=X2-X1
      ELSE
         RTBIS=X2
         DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
         DX=DX*.5
         XMID=RTBIS+DX
         CALL FNDSHD( XMID,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
         FMID=ANDEXD(XMID,SH,GAMMA)*(XMID+RE)-CPATH
         IF(FMID.LE.0.)RTBIS=XMID
         IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
   11 END DO

      !     COMES HERE IF UNABLE TO SOLVE.
      IF (ABS(CX2) .LT. ABS(CX1)) THEN
         RTBIS = X2
      ELSE
         RTBIS = X1
      ENDIF

      RETURN
      END FUNCTION

!-----------------------------------------------------------------------
!     SUBROUTINE RFPATH ()
!
!     -------------------------------------------------------------
!     This routine was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     It was also modified to eliminate GOTO statements in order to
!     make the program easier to understand.
!     These changes were obtained from H. Snell (March, 1996).
!     -------------------------------------------------------------
!
!     THIS SUBROUTINE TRACES THE REFRACTED RAY FROM H1 WITH AN
!     INITIAL ZENITH ANGLE ANGLE TO H2 WHERE THE ZENITH ANGLE IS PHI,
!     AND CALCULATES THE ABSORBER AMOUNTS (IF IAMT.EQ.1) ALONG
!     THE PATH.  IT STARTS FROM THE LOWEST POINT ALONG THE PATH
!     (THE TANGENT HEIGHT HMIN IF LEN = 1 OR HA = MIN(H1,H2) IF LEN = 0
!     AND PROCEEDS TO THE HIGHEST POINT.  BETA AND RANGE ARE THE
!     EARTH CENTERED ANGLE AND THE TOTAL DISTANCE RESPECTIVELY
!     FOR THE REFRACTED PATH FROM H1 TO H2, AND BENDNG IS THE TOTAL
!     BENDING ALONG THE PATH
!
!-----------------------------------------------------------------------
      SUBROUTINE RFPATH( H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA,BENDNG,&
                         ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, RE, &
                         ZOUT,ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX, &
                         DENM,AMTP,NMOL )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, DEG
      USE Module_Config      ,ONLY: IPR,noPrnt
      IMPLICIT NONE

      real    ,intent(inout) :: H1
      real    ,intent(in)    :: H2
      real    ,intent(in)    :: ANGLE
      real    ,intent(inout) :: PHI
      integer ,intent(in)    :: LEN
      real    ,intent(in)    :: HMIN
      integer ,intent(in)    :: IAMT
      real    ,intent(out)   :: RANGE
      real    ,intent(out)   :: BETA
      real    ,intent(out)   :: BENDNG
      !
      real    ,intent(in)    :: ZBND(:)
      real    ,intent(inout) :: ZMDL(:)
      real    ,intent(in)    :: PM(:)
      real    ,intent(in)    :: TM(:)
      real    ,intent(in)    :: RFNDXM(:)
      integer ,intent(in)    :: IBMAX
      integer ,intent(in)    :: IMMAX
      real    ,intent(in)    :: RE
      !
      real    ,intent(out)   :: ZOUT(:)
      real    ,intent(out)   :: ZPTH(:)
      real    ,intent(out)   :: PP(:)
      real    ,intent(out)   :: TP(:)
      real    ,intent(out)   :: SP(:)
      real    ,intent(out)   :: PPSUM(:)
      real    ,intent(out)   :: TPSUM(:)
      real    ,intent(out)   :: RHOPSM(:)
      integer ,intent(out)   :: IOUTMX
      integer ,intent(out)   :: IPMAX
      !
      real    ,intent(in)    :: DENM(:,:)
      real    ,intent(out)   :: AMTP(:,:)
      integer ,intent(in)    :: NMOL


      !---Local variables
      !
      character*2 ,PARAMETER :: HLOW(2)=['H1','H2']
      integer     ,PARAMETER :: I_2=2

      INTEGER  :: IHIGH,  IHLOW,   IORDER
      INTEGER  :: J,      J2
      REAL     :: ANGLEA
      REAL(r8) :: COSAI,  CPATH,   DBEND
      REAL     :: DBETA
      REAL(r8) :: DS,     GAMMA
      REAL     :: HA,     HB,      PBAR,    RHOBAR
      REAL(r8) :: S,      SH,      SINAI
      REAL     :: TBAR,   THETA
      integer  :: IPHMID

      real  ,allocatable :: RFNDXP(:) !RFNDXP(IM2)
      real  ,allocatable :: DENP(:,:) !DENP(MXMOL+MX_XS, MXPDIM)



      !--- Zero out the common block variables
      !DO N = 1, IPDIM
      !   IF (N.LE.IPDIM-2) THEN
      !      ZPTH(N) = 0.
      !      PP(N) = 0.
      !      TP(N) = 0.
      !      !RFNDXP(N) = 0.
      !      SP(N) = 0.
      !      PPSUM(N) = 0.
      !      TPSUM(N) = 0.
      !      RHOPSM(N) = 0.
      !   ENDIF
      !   DO M = 1, (MXMOL + MX_XS)
      !      !DENP(M,N) = 0.
      !      AMTP(M,N) = 0.
      !   ENDDO
      !END DO

      ZOUT(:) = 0.
      ZPTH(:) = 0.
      PP(:) = 0.
      TP(:) = 0.
      SP(:) = 0.
      PPSUM(:) = 0.
      TPSUM(:) = 0.
      RHOPSM(:) = 0.
      AMTP(:,:) = 0.

      allocate( RFNDXP(    IMMAX+IBMAX+3) ); RFNDXP(:)=0.
      allocate( DENP( NMOL,IMMAX+IBMAX+3) ); DENP(:,:)=0.



      ! REORDER H1 AND H2 TO HA AND HB (HA .LE. HB)

      IF (H1.LE.H2) THEN
         IORDER = 1
         HA = H1
         HB = H2
         ANGLEA = ANGLE
      ELSE
         IORDER = -1
         HA = H2
         HB = H1
         ANGLEA = PHI
      ENDIF

      ! MERGE THE ATMOSPHERIC PROFILE STORED IN ZMDL WITH H1,H2,(HMIN) AN
      ! THE BOUNDARIES ZBND
      CALL AMERGE( H1,H2,HMIN,LEN, &
                   ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, &
                   ZOUT,ZPTH,PP,TP,RFNDXP,IOUTMX,IPMAX,IPHMID, &
                   DENM,DENP,NMOL )
      IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,900)

      ! CALCULATE CPATH SEPERATELY FOR LEN = 0,1

      IF (LEN.EQ.0) THEN
         CALL FNDSHD( HA,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
         CPATH = (RE+HA)*ANDEXD(HA,SH,GAMMA)*SIN(ANGLEA/DEG)
      ELSE
         CALL FNDSHD( HMIN,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
         CPATH = (RE+HMIN)*ANDEXD(HMIN,SH,GAMMA)
      ENDIF

      BETA = 0.0
      S = 0.0
      BENDNG = 0.0
      IF (LEN.EQ.1) THEN

         ! TANGENT PATH

         IF (IORDER.EQ.-1) THEN
            IHLOW = 2
         ELSE
            IHLOW = 1
         ENDIF
         IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,905) HLOW(IHLOW)
         SINAI = 1.0
         COSAI = 0.0
         THETA = 90.0
      ELSE

         !     SHORT PATH
         !
         !     ANGLEA IS THE ZENITH ANGLE AT HA IN DEG
         !     SINAI IS SIN OF THE INCIDENCE ANGLE
         !     COSAI IS CARRIED SEPERATELY TO AVOID A PRECISION PROBLEM
         !     WHEN SINAI IS CLOSE TO 1.0

         THETA = ANGLEA
         IF (ANGLEA.LE.45.0) THEN
            SINAI = SIN(ANGLEA/DEG)
            COSAI = -COS(ANGLEA/DEG)
         ELSE
            SINAI = COS((90.0-ANGLEA)/DEG)
            COSAI = -SIN((90.0-ANGLEA)/DEG)
         ENDIF
         IF (IORDER.EQ.-1) THEN
            IHLOW = 2
         ELSE
            IHLOW = 1
         ENDIF
         IHIGH = MOD(IHLOW,I_2)+1
         IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,910) HLOW(IHLOW),    &
         HLOW(IHIGH)
      ENDIF

      ! LOOP OVER THE LAYERS

      J2 = IPMAX-1
      DO 100 J = 1, J2
         CALL SCLHTD (ZPTH(J),ZPTH(J+1),RFNDXP(J),RFNDXP(J+1),SH,GAMMA)
         CALL ALAYER( J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,DS,DBEND, &
                      ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,&
                      DENP,AMTP,NMOL,RE )
         DBEND = DBEND*DEG
         PHI = ASIN(SINAI)*DEG
         DBETA = THETA-PHI+DBEND
         PHI = 180.0-PHI
         S = S+DS
         BENDNG = BENDNG+DBEND
         BETA = BETA+DBETA
         IF (IAMT.EQ.1) THEN
            PBAR = PPSUM(J)/RHOPSM(J)
            TBAR = TPSUM(J)/RHOPSM(J)
            RHOBAR = RHOPSM(J)/DS
            IF (NOPRNT.ge.0) WRITE (IPR,915) J,ZPTH(J),ZPTH(J+1),       &
            THETA,DS,S,DBETA,BETA,PHI,DBEND,BENDNG,PBAR, TBAR,RHOBAR
         ENDIF
         THETA = 180.0-PHI

         IF (LEN.EQ.1) THEN

            ! For tangent paths, double the quantities BENDNG,BETA,
            ! and S for the symmetric part of the path

            IF ((J+1).EQ.IPHMID) THEN
               BENDNG = 2.0*BENDNG
               BETA = 2.0*BETA
               S = 2.0*S
               IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,920) S,BETA,   &
               BENDNG
               IF (IPHMID.NE.IPMAX) THEN
                  IF (IORDER.EQ.-1) THEN
                     IHLOW = 2
                  ELSE
                     IHLOW = 1
                  ENDIF
                  IHIGH = MOD(IHLOW,I_2)+1
                  IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,910) HLOW(  &
                  IHLOW),HLOW(IHIGH)
               ENDIF
            ENDIF
         ENDIF
  100 END DO
      IF (IORDER.EQ.-1) PHI = ANGLEA
      RANGE = S


      deallocate(RFNDXP)
      deallocate(DENP)

      RETURN

  900 FORMAT ('1CALCULATION OF THE REFRACTED PATH THROUGH THE ',        &
     &        'ATMOSPHERE',///,T5,'I',T14,'ALTITUDE',T30,'THETA',T38,   &
     &        'DRANGE',T47,'RANGE',T57,'DBETA',T65,'BETA',T76,'PHI',    &
     &        T84,'DBEND',T91,'BENDING',T102,'PBAR',T111,'TBAR',T119,   &
     &        'RHOBAR',/,T11,'FROM',T22,'TO',/,T11,'(KM)',T21,'(KM)',   &
     &        T30,'(DEG)',T39,'(KM)',T48,'(KM)',T57,'(DEG)',T65,        &
     &        '(DEG)',T75,'(DEG)',T84,'(DEG)',T92,'(DEG)',T102,'(MB)',  &
     &        T112,'(K)',T117,'(MOL CM-3)',/)
  905 FORMAT (' ',T10,'TANGENT',T20,A2,/,T10,'HEIGHT',/)
  910 FORMAT (' ',T14,A2,' TO ',A2,/)
  915 FORMAT (' ',I4,2F10.3,10F9.3,1PE9.2)
  920 FORMAT ('0',T10,'DOUBLE RANGE, BETA, BENDING',/,T10,              &
     &        'FOR SYMMETRIC PART OF PATH',T44,F9.3,T62,F9.3,T89,       &
     &        F9.3,/)

      END  SUBROUTINE


!-----------------------------------------------------------------------
!     SUBROUTINE AMERGE ()
!
!     AMERGE CREATES A SET OF LAYER BOUNDARIES ZOUT WHICH INCLUDES
!     HMIN, (HMID), HMAX AND ALL OF ZBND BETWEEN HMIN AND HAMX.
!     ZOUT DEFINES THE LAYERS FOR THE LBLRTM CALCULATION.
!     ZOUT IS THEN MERGED WITH THE ATMOSPHERIC PROFILE IN ZMDL INTO ZPT
!     INTERPOLATING TO THE LEVELS ZOUT WHEN NECESSARY.  THE RAY
!     TRACE IS CALCULATED USING THE PROFILE IN ZPTH.
!-----------------------------------------------------------------------
      SUBROUTINE AMERGE( H1,H2,HMIN,LEN, &
                         ZBND,ZMDL,PM,TM,RFNDXM,IBMAX,IMMAX, &
                         ZOUT,ZPTH,PP,TP,RFNDXP,IOUTMX,IPMAX,IPHMID, &
                         DENM,DENP,NMOL )
!-----------------------------------------------------------------------
      USE Module_Utility   ,ONLY: EXPINT
      USE Module_Config    ,ONLY: IPR
      IMPLICIT NONE


      real     ,intent(inout) :: H1
      real     ,intent(in)    :: H2
      real     ,intent(in)    :: HMIN
      integer  ,intent(in)    :: LEN
      !
      real     ,intent(in)    :: ZBND(:)
      real     ,intent(inout) :: ZMDL(:)
      real     ,intent(in)    :: PM(:)
      real     ,intent(in)    :: TM(:)
      real     ,intent(in)    :: RFNDXM(:)
      integer  ,intent(in)    :: IBMAX
      integer  ,intent(in)    :: IMMAX
      !
      real     ,intent(out)   :: ZOUT(:)
      real     ,intent(out)   :: ZPTH(:)
      real     ,intent(out)   :: PP(:)
      real     ,intent(out)   :: TP(:)
      real     ,intent(out)   :: RFNDXP(:)
      integer  ,intent(out)   :: IOUTMX
      integer  ,intent(out)   :: IPMAX
      integer  ,intent(out)   :: IPHMID
      !
      real     ,intent(in)    :: DENM(:,:)
      real     ,intent(out)   :: DENP(:,:)
      integer  ,intent(in)    :: NMOL


      !---Local variables
      !
      real    ,PARAMETER :: TOL=5.E-4
      integer ,PARAMETER :: I_2=2

      INTEGER  :: I1,    IB,    IH,     IHMAX
      INTEGER  :: IM,    IOUT,   IP
      INTEGER  :: JM,    K
      REAL     :: A,     HMAX,  HMID
      REAL     :: ZH(3)



      ! HMID .EQ. MINIMUM OF H1, H2

      HMID =   MIN(H1,H2)
      HMAX =   MAX(H1,H2)
      IHMAX = 2
      ZH(1) = HMIN
      IF (LEN.EQ.0) THEN
         ZH(2) = HMAX
      ELSE
         ZH(2) = HMID
         IF (ABS(H1-H2).LT.TOL) H1 = H2
         IF (H1.NE.H2) THEN
            IHMAX = 3
            ZH(3) = HMAX
         ENDIF
      ENDIF

      ! MERGE ZH AND ZBND BETWEEN ZH(1) AND ZH(IHMAX) TO CREAT ZOUT

      ZOUT(1) = ZH(1)
      DO 30 I1 = 1, IBMAX
         IF (ABS(ZBND(I1)-ZH(1)).LT.TOL) ZH(1) = ZBND(I1)
         IF (ZBND(I1).GT.ZH(1)) GO TO 40
   30 END DO
      I1 = IBMAX
   40 CONTINUE

      ! ZBND(I1) IS SMALLEST ZBND .GT. ZH(1)

      IOUT = 1
      IB = I1
      IH = 2
   50 CONTINUE
         IOUT = IOUT+1
         IF (IB.GT.IBMAX) GO TO 60
         IF (ABS(ZBND(IB)-ZH(IH)).LT.TOL) ZH(IH) = ZBND(IB)
         IF (ZBND(IB).LT.ZH(IH)) GO TO 70
         IF (ZBND(IB).EQ.ZH(IH)) IB = IB+1

         ! INSERT ZH(IH)

   60    CONTINUE
         ZOUT(IOUT) = ZH(IH)
         IH = IH+1
         IF (IH.GT.IHMAX) GO TO 80
         GO TO 50

         ! INSERT ZBND(IB)

   70    CONTINUE
         ZOUT(IOUT) = ZBND(IB)
         IB = IB+1
      GO TO 50
   80 CONTINUE
      IOUTMX = IOUT

      ! NOW MERGE ZOUT AND ZMDL INTO ZPTH (FROM ZOUT(1) TO ZOUT(IOUTMX))
      ! AND INTERPOLATE PRESSURE, TEMPERATURE, AND DENSITY WHEN
      ! NECESSARY
      !
      ! FIND SMALLEST ZMDL .GT. HMIN
      !
      DO 90 IM = 1, IMMAX
         IF (ZMDL(IM).GE.HMIN) GO TO 100
   90 END DO
      WRITE (IPR,900) HMIN
      STOP ' AMERGE - HMIN '
  100 CONTINUE
      IPHMID = 0
      IP = 0
      IOUT = 1
  110 CONTINUE
         IP = IP+1
         !yma IF (IP.GT.IPDIM) THEN
         IF (IP.GT.(IMMAX+IBMAX+3)) THEN
            WRITE (IPR,905) (IMMAX+IBMAX+3)  !yma WRITE (IPR,905) IPDIM
            STOP ' AMERGE - IPDIM '
         ENDIF
         IF (IM.GT.IMMAX) GO TO 130
         IF (ABS(ZOUT(IOUT)-ZMDL(IM)).LT.TOL) ZMDL(IM) = ZOUT(IOUT)
         IF (ZOUT(IOUT).LT.ZMDL(IM)) GO TO 130
         IF (ZOUT(IOUT).EQ.ZMDL(IM)) IOUT = IOUT+1

         ! INSERT ZMDL(IM)

         ZPTH(IP) = ZMDL(IM)
         PP(IP) = PM(IM)
         TP(IP) = TM(IM)
         RFNDXP(IP) = RFNDXM(IM)
         DO 120 K = 1, NMOL
            DENP(K,IP) = DENM(K,IM)
  120    END DO
         IM = IM+1
         IF (ABS(ZPTH(IP)-HMID).LT.TOL) HMID = ZPTH(IP)
         IF (ZPTH(IP).EQ.HMID) IPHMID = IP
         IF (ABS(ZPTH(IP)-ZOUT(IOUTMX)).LT.TOL) ZOUT(IOUTMX) = ZPTH(IP)
         IF (ZPTH(IP).EQ.ZOUT(IOUTMX)) GO TO 150
         GO TO 110

         !INSERT LEVEL FROM ZOUT(IOUT) AND INTERPOLATE

  130    CONTINUE
         ZPTH(IP) = ZOUT(IOUT)
         JM = IM
         JM = MAX(JM,I_2)
         A = (ZOUT(IOUT)-ZMDL(JM-1))/(ZMDL(JM)-ZMDL(JM-1))
         CALL EXPINT (PP(IP),PM(JM-1),PM(JM),A)
         TP(IP) = TM(JM-1)+(TM(JM)-TM(JM-1))*A
         CALL EXPINT (RFNDXP(IP),RFNDXM(JM-1),RFNDXM(JM),A)
         DO 140 K = 1, NMOL
            CALL EXPINT (DENP(K,IP),DENM(K,JM-1),DENM(K,JM),A)
  140    END DO
         IF (ABS(ZPTH(IP)-HMID).LT.TOL) ZPTH(IP) = HMID
         IF (ZPTH(IP).EQ.HMID) IPHMID = IP
         IOUT = IOUT+1
         IF (ABS(ZPTH(IP)-ZOUT(IOUTMX)).LT.TOL) ZPTH(IP) = ZOUT(IOUTMX)
         IF (ZPTH(IP).EQ.ZOUT(IOUTMX)) GO TO 150
         GO TO 110
  150 CONTINUE
      IPMAX = IP

      RETURN

  900 FORMAT ('0FROM AMERGE- ATMOSPHERIC PROFILE IN ZMDL DOES NOT',     &
     &        ' EXTEND UP TO HMIN = ',E12.5)
  905 FORMAT ('0FROM AMERGE- MERGING THE ATMOSPHERIC PROFILE AND THE ', &
     &        'LBLRTM BOUNDARIES INTO ZPTH(IPDIM) EXCEEDS THE ',        &
     &        'DIMENSION IPDIM = ',I5)

      END SUBROUTINE


!-----------------------------------------------------------------------
!     SUBROUTINE ALAYER ()
!
!     -------------------------------------------------------------
!     This routine was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     It was also modified to eliminate GOTO statements in order to
!     make the program easier to understand.
!     These changes were obtained from H. Snell (March, 1996).
!     -------------------------------------------------------------
!
!     THIS SUBROUTINE TRACES THE OPTICAL RAY THROUGH ONE LAYER FROM
!     Z1 TO Z2 AND IF IAMT.NE.2 CALCULATES THE INTEGRATED ABSORBER
!     AMOUNTS FOR THE LAYER. SINAI IS THE SIN OF THE INITIAL INCIDENCE
!     ANGLE (= 180 - ZENITH ANGLE). COSAI IS CARRIED SEPERATELY TO
!     AVOID A PRECISION PROBLEM NEAR SINAI = 1. CPATH IS THE CONSTANT
!     OF REFRACTION FOR THE PATH = INDEX*RADIUS*SINAI, SH AND GAMMA ARE
!     THE SCALE HEIGHT AND THE AMOUNT AT THE GROUND FOR THE REFRACTIVIT
!     (= 1-INDEX OF REFRACTION), S IS THE REFRACTED PATH LENGTH THROUGH
!     THE LAYER, BETA IS THE EARTH CENTERED ANGLE, AND BEND IS THE
!     BENDING THROUGH THE LAYER. IAMT CONTROLS WHETHER AMOUNTS ARE
!     CALCULATED OR NOT.
!-----------------------------------------------------------------------
      SUBROUTINE ALAYER( J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,S,BEND, &
                         ZPTH,PP,TP,SP,PPSUM,TPSUM,RHOPSM,&
                         DENP,AMTP,NMOL,RE )
!-----------------------------------------------------------------------
      USE Module_ConstParam, ONLY: r8=>kind_r8, GCAIR, DELTAS
      IMPLICIT NONE

      integer  ,intent(in)    :: J
      real(r8) ,intent(inout) :: SINAI
      real(r8) ,intent(inout) :: COSAI
      real(r8) ,intent(in)    :: CPATH
      real(r8) ,intent(in)    :: SH, GAMMA
      integer  ,intent(in)    :: IAMT
      real(r8) ,intent(out)   :: S
      real(r8) ,intent(out)   :: BEND
      !
      real     ,intent(in)    :: ZPTH(:)
      real     ,intent(in)    :: PP(:)
      real     ,intent(in)    :: TP(:)
      real     ,intent(inout) :: SP(:)
      real     ,intent(inout) :: PPSUM(:)
      real     ,intent(inout) :: TPSUM(:)
      real     ,intent(inout) :: RHOPSM(:)
      real     ,intent(in)    :: DENP(:,:)
      real     ,intent(inout) :: AMTP(:,:)
      integer  ,intent(in)    :: NMOL
      real     ,intent(in)    :: RE


      !--- Local variables
      !
      real ,PARAMETER :: EPSILN=1.0E-5

      INTEGER  :: K, N
      REAL(r8) :: COSAI1,    COSAI2,    COSAI3,    D21
      REAL(r8) :: D31,       D32,       DBEND,     DBNDX1
      REAL(r8) :: DBNDX2,    DBNDX3
      REAL(r8) :: DH,        DHMIN,     DS,        DSDX1
      REAL(r8) :: DSDX2,     DSDX3
      REAL     :: DSDZ
      REAL(r8) :: DX
      REAL     :: DZ,        H1,        H2,        H3,    HP
      REAL     :: HRHO,      PA,        PB
      REAL(r8) :: R1,        R2,        R3,        RATIO1
      REAL(r8) :: RATIO2,    RATIO3
      REAL     :: RHOA,      RHOB
      REAL(r8) :: SINAI1,    SINAI2,    SINAI3
      REAL     :: TA,        TB
      REAL(r8) :: W1,        W2,        W3,        X1
      REAL(r8) :: X2,        X3,        Y1,        Y3
      REAL     :: Z1,        Z2


      real ,allocatable :: DENA(:), DENB(:), HDEN(:)



      allocate( DENA(NMOL) )
      allocate( DENB(NMOL) )
      allocate( HDEN(NMOL) )

      ! INITIALIZE VARIABLES FOR THE CALCULATION OF THE PATH

      N = 0
      Z1 = ZPTH(J)
      Z2 = ZPTH(J+1)
      H1 = Z1
      R1 = RE+H1
      DHMIN = DELTAS**2/(2.0*R1)
      SINAI1 = SINAI
      COSAI1 = COSAI
      IF ((1.0-SINAI).LT.EPSILN)                                        &
     &     Y1 = COSAI1**2/2.0+COSAI1**4/8.0+COSAI1**6*3.0/48.0
      Y3 = 0.0
      X1 = -R1*COSAI1
      RATIO1 = R1/RADRFD(H1,SH,GAMMA)
      DSDX1 = 1.0/(1.0-RATIO1*SINAI1**2)
      DBNDX1 = DSDX1*SINAI1*RATIO1/R1
      S = 0.0
      BEND = 0.0
      IF (IAMT.NE.2) THEN

         ! Initialize the variables for the calculation of the
         ! absorber amounts

         PA = PP(J)
         PB = PP(J+1)
         IF (PB.EQ.PA) THEN
            WRITE(*,*) PB
            STOP 'LBLATM: PRESSURES IN ADJOINING LAYERS MUST DIFFER'
         ENDIF
         TA = TP(J)
         TB = TP(J+1)
         RHOA = PA/(GCAIR*TA)
         RHOB = PB/(GCAIR*TB)
         DZ = ZPTH(J+1)-ZPTH(J)
         HP = -DZ/ LOG(PB/PA)
         IF (ABS(RHOB/RHOA-1.0).GE.EPSILN) THEN
            HRHO = -DZ/ LOG(RHOB/RHOA)
         ELSE
            HRHO = 1.0E30
         ENDIF
         DO 40 K = 1, NMOL
            DENA(K) = DENP(K,J)
            DENB(K) = DENP(K,J+1)
            IF ((DENA(K).EQ.0.0.OR.DENB(K).EQ.0.0).OR. (ABS(1.0-DENA(K)/&
            DENB(K)).LE.EPSILN)) THEN

               ! Use linear interpolation

               HDEN(K) = 0.0
            ELSE

               ! Use exponential interpolation

               HDEN(K) = -DZ/ LOG(DENB(K)/DENA(K))
            ENDIF
   40    CONTINUE
      ENDIF

      ! LOOP THROUGH PATH
      ! INTEGRATE PATH QUANTITIES USING QUADRATIC INTEGRATION WITH
      ! UNEQUALLY SPACED POINTS

   60 CONTINUE
      N = N+1
      DH = -DELTAS*COSAI1
      DH = MAX(DH,DHMIN)
      H3 = H1+DH
      IF (H3.GT.Z2) H3 = Z2
      DH = H3-H1
      R3 = RE+H3
      H2 = H1+DH/2.0
      R2 = RE+H2
      SINAI2 = CPATH/(ANDEXD(H2,SH,GAMMA)*R2)
      SINAI3 = CPATH/(ANDEXD(H3,SH,GAMMA)*R3)
      RATIO2 = R2/RADRFD(H2,SH,GAMMA)
      RATIO3 = R3/RADRFD(H3,SH,GAMMA)
      IF ((1.0-SINAI2).LE.EPSILN) THEN

         ! Near a tangent height, COSAI = -SQRT(1-SINAI**2) loses
         ! precision. use the following algorithm to get COSAI.

         Y3 = Y1+(SINAI1*(1.0-RATIO1)/R1+4.0*SINAI2*(1.0-RATIO2)/R2+    &
         SINAI3*(1.0-RATIO3)/R3)*DH/6.0
         COSAI3 = -SQRT(2.0*Y3-Y3**2)
         X3 = -R3*COSAI3
         DX = X3-X1
         W1 = 0.5*DX
         W2 = 0.0
         W3 = 0.5*DX
      ELSE
         COSAI2 = -SQRT(1.0-SINAI2**2)
         COSAI3 = -SQRT(1.0-SINAI3**2)
         X2 = -R2*COSAI2
         X3 = -R3*COSAI3

         ! Calculate weights

         D31 = X3-X1
         D32 = X3-X2
         D21 = X2-X1
         IF (D32.EQ.0.0.OR.D21.EQ.0.0) THEN
            W1 = 0.5*D31
            W2 = 0.0
            W3 = 0.5*D31
         ELSE
            W1 = (2.0-D32/D21)*D31/6.0
            W2 = D31**3/(D32*D21*6.0)
            W3 = (2.0-D21/D32)*D31/6.0
         ENDIF
      ENDIF
      DSDX2 = 1.0/(1.0-RATIO2*SINAI2**2)
      DSDX3 = 1.0/(1.0-RATIO3*SINAI3**2)
      DBNDX2 = DSDX2*SINAI2*RATIO2/R2
      DBNDX3 = DSDX3*SINAI3*RATIO3/R3

      ! INTEGRATE

      DS = W1*DSDX1+W2*DSDX2+W3*DSDX3
      S = S+DS
      DBEND = W1*DBNDX1+W2*DBNDX2+W3*DBNDX3
      BEND = BEND+DBEND
      IF (IAMT.NE.2) THEN

         ! Calculate amounts

         DSDZ = DS/DH
         PB = PA*EXP(-DH/HP)
         RHOB = RHOA*EXP(-DH/HRHO)
         IF ((DH/HRHO).GE.EPSILN) THEN
            PPSUM(J) = PPSUM(J)+DSDZ*(HP/(1.0+HP/HRHO)) *(PA*RHOA-PB*   &
            RHOB)
            TPSUM(J) = TPSUM(J)+DSDZ*HP*(PA-PB)/GCAIR
            RHOPSM(J) = RHOPSM(J)+DSDZ*HRHO*(RHOA-RHOB)
         ELSE
            PPSUM(J) = PPSUM(J)+0.5*DS*(PA*RHOA+PB*RHOB)
            TPSUM(J) = TPSUM(J)+0.5*DS*(PA+PB)/GCAIR
            RHOPSM(J) = RHOPSM(J)+0.5*DS*(RHOA+RHOB)
         ENDIF
         DO 130 K = 1, NMOL
            IF ((HDEN(K).EQ.0.0).OR. (ABS(DH/HDEN(K)).LT.EPSILN)) THEN

               ! Linear interpolation
               ! 1.0E05 factor converts units km to cm

               DENB(K)=DENP(K,J)+(DENP(K,J+1)-DENP(K,J))*(H3-Z1)/DZ
               AMTP(K,J) = AMTP(K,J)+0.5*(DENA(K)+DENB(K))*DS*1.0E5
            ELSE

               ! Exponential interpolation

               DENB(K) = DENP(K,J)*EXP(-(H3-Z1)/HDEN(K))
               AMTP(K,J) = AMTP(K,J)+DSDZ*HDEN(K) *(DENA(K)-DENB(K))*   &
               1.0E5
            ENDIF
  130    CONTINUE
         PA = PB
         RHOA = RHOB
         DO 140 K = 1, NMOL
            DENA(K) = DENB(K)
  140    CONTINUE
      ENDIF

      IF (H3.LT.Z2) THEN
         H1 = H3
         R1 = R3
         SINAI1 = SINAI3
         RATIO1 = RATIO3
         Y1 = Y3
         COSAI1 = COSAI3
         X1 = X3
         DSDX1 = DSDX3
         DBNDX1 = DBNDX3
      ELSE
         SINAI = SINAI3
         COSAI = COSAI3
         SP(J) = S
         RETURN
      ENDIF

      GO TO 60


      !deallocate( DENA, DENB, HDEN )

      END  SUBROUTINE

!-----------------------------------------------------------------------
!
! FPACK TAKES THE AMOUNTS STORED IN THE LAYERS DEFINED BY ZPTH AND
! PACKS THEM INTO THE LAYERS DEFINED BY ZOUT.  IT ALSO ZEROS OUT
! LAYER AMOUNTS IF THE AMOUNT FOR THAT LAYER AND ABOVE IS LESS
! THAN 0.1 PERCENT OF THE TOTAL FOR THAT MOLECULE if THE
! n_zero OPTION IS SELECTED.
!
! * In CLBLM/FPACK, cross-section molecules are treated the same
!   as line molecules. they all counted when computing PBAR, TBAR
!   and layer amounts and may be zeroed out. Please be noted that
!   in LBLRTM/FPACK, PBAR, TBAR and layer amounts are computed
!   without taking into account cross-section molecules and
!   cross-section can not being zeroed out either.
!
!         I2 = IPMAX-1
!         IOUT = 1
!         DO 100 IP = 1, I2
!            DO 90 K = 1, IXMOLS
!               XAMNT(K,IOUT) = XAMNT(K,IOUT)+AMTP(K,IP)
!   90       CONTINUE
!            IF (ZPTH(IP+1).EQ.ZOUT(IOUT+1)) IOUT = IOUT+1
!  100    CONTINUE
!         IF (IOUT.NE.IOUTMX) THEN
!            WRITE (IPR,935) IOUT,IOUTMX
!            STOP 'STOPPED IN XAMNTS, IOUT .NE. IOUTMX'
!         ENDIF
!         IOUTMX = IOMXSV
!         LMAX = IOUTMX-1
!
!-----------------------------------------------------------------------
      SUBROUTINE FPACK ( H1,H2,LEN,n_zero, &
                         ZOUT,ZPTH,PP,TP,SP,  &
                         PPSUM,TPSUM,RHOPSM,IOUTMX,IPMAX,AMTP,nMol,molID, &
                         ALTZ,PZ,TZ,SOUT,PBAR, TBAR, RHOSUM,   AMOUNT,&
                         SECNTA,IPATH )
!-----------------------------------------------------------------------
      USE Module_Config     ,ONLY: IPR
      USE Module_ConstParam ,ONLY: molIndex
      IMPLICIT NONE

      real        ,intent(in)    :: H1, H2
      integer     ,intent(in)    :: LEN
      integer     ,intent(in)    :: n_zero
      !
      real        ,intent(in)    :: ZOUT(:)
      real        ,intent(in)    :: ZPTH(:)
      real        ,intent(in)    :: PP(:)
      real        ,intent(in)    :: TP(:)
      real        ,intent(in)    :: SP(:)
      real        ,intent(in)    :: PPSUM(:)
      real        ,intent(in)    :: TPSUM(:)
      real        ,intent(in)    :: RHOPSM(:)
      integer     ,intent(inout) :: IOUTMX
      integer     ,intent(in)    :: IPMAX
      real        ,intent(in)    :: AMTP(:,:)
      integer     ,intent(in)    :: nMol
      character(*),intent(in)    :: molID(:)
      !
      real        ,intent(out)   :: ALTZ(0:)
      real        ,intent(out)   :: PZ(0:)
      real        ,intent(out)   :: TZ(0:)
      real        ,intent(out)   :: SOUT(:)
      real        ,intent(out)   :: PBAR(:)
      real        ,intent(out)   :: TBAR(:)
      real        ,intent(out)   :: RHOSUM(:)
      real        ,intent(out)   :: AMOUNT(:,:)
      real        ,intent(out)   :: SECNTA(:)
      integer     ,intent(out)   :: IPATH(:)


      !---Local variables
      !
      integer :: O2Ind, I, I2, IP, IOUT, K, L2, L, LMAX
      real    :: HMID
      integer :: ISKPT, nmol_max
      real    :: FAC

      real     ,allocatable :: AMTTOT(:)
      real     ,allocatable :: AMTCUM(:)
      integer  ,allocatable :: ISKIP(:)




      HMID = MIN(H1,H2)

      allocate( AMTTOT(nMOL) )
      allocate( AMTCUM(nMOL) )
      allocate( ISKIP(nMOL) )

      AMTTOT(:) = 0.0
      DO I = 1, IPMAX-1
         FAC = 1.0
         IF (LEN.EQ.1.AND.ZPTH(I+1).LE.HMID) FAC = 2.0
         DO K = 1, nMol
            AMTTOT(K) = AMTTOT(K)+FAC*AMTP(K,I)
         ENDDO
      ENDDO

      SOUT(:)     = 0.
      PBAR(:)     = 0.
      TBAR(:)     = 0.
      RHOSUM(:)   = 0.
      AMOUNT(:,:) = 0.


      I2    = IPMAX-1
      IOUT  = 1
      PZ(0) = PP(1)
      TZ(0) = TP(1)

!yma160504      ! If entry in TAPE5 for TBOUND < 0, use TZ(O) as boundary temperature
!yma160504      IF (TBOUND.LT.0.) TBOUND = TZ(0)

      DO 20 IP = 1, I2
         PBAR(IOUT)   = PBAR(IOUT)+PPSUM(IP)
         TBAR(IOUT)   = TBAR(IOUT)+TPSUM(IP)
         RHOSUM(IOUT) = RHOSUM(IOUT)+RHOPSM(IP)
         SOUT(IOUT)   = SOUT(IOUT)+SP(IP)
         DO 10 K = 1, NMOL
            AMOUNT(K,IOUT) = AMOUNT(K,IOUT)+AMTP(K,IP)
   10    CONTINUE
         IF (ZPTH(IP+1).EQ.ZOUT(IOUT+1)) THEN
            PZ(IOUT) = PP(IP+1)
            TZ(IOUT) = TP(IP+1)
            IOUT = IOUT+1
         ENDIF
   20 END DO
      IF (IOUT.NE.IOUTMX) GO TO 110

      ! CALCULATE THE DENSITY WEIGHTED PRESSURE AND TEMPERATURE AND
      ! ZERO OUT LAYER AMOUNTS AFTER 99.9 PERCENT OF THE TOTAL
      !
!      !iskip(7) = 0
!      if ( nLnMol>=7 ) iskip(7) = 0 !yma1604504
      O2Ind = molIndex('O2',molID) !O2 is a member of molID, O2Ind>0
      if (O2Ind>0) iskip( O2Ind ) = 0 !yma 180520

      DO 30 K = 1, NMOL
         AMTCUM(K) = 0.0
         ISKIP(K) = 0
         IF (AMTTOT(K).EQ.0.0) ISKIP(K) = 1
   30 END DO

      L2 = IOUTMX-1
      LMAX = L2
      DO 90 L = 1, L2
         PBAR(L) = PBAR(L)/RHOSUM(L)
         TBAR(L) = TBAR(L)/RHOSUM(L)

         ! ADJUST RHOSUM FOR THE PATH LENGTH IN CM NOT KM

         RHOSUM(L) = RHOSUM(L)*1.0E+5


         ! CALCULATE 'EFFECTIVE SECANT' SECNTA

         SECNTA(L) = SOUT(L)/(ZOUT(L+1)-ZOUT(L))
         IF (L.EQ.1) ALTZ(0) = ZOUT(1)
         ALTZ(L) = ZOUT(L+1)

         ! SET  IPATH

         IF (LEN.EQ.1) GO TO 50
         IF (H1.LT.H2) IPATH(L) = 3
2         IF (H1.GT.H2) IPATH(L) = 1
         GO TO 60
   50    CONTINUE
         IF (ZOUT(L).LT.HMID) IPATH(L) = 2
         IF (ZOUT(L).GE.HMID.AND.H1.GT.H2) IPATH(L) = 1
         IF (ZOUT(L).GE.HMID.AND.H1.LT.H2) IPATH(L) = 3
   60    CONTINUE

         ! TEST FOR ZEROING OF AMOUNTS

         ISKPT = 0
         nmol_max = nmol

         !!IF (ISKIP(7).EQ.1) nmol_max = nmol - 1
         !if ( nLnMol>=7 ) then
         !   IF (ISKIP(7).EQ.1 ) nmol_max = nmol - 1 !yma160504
         !endif
         if ( O2Ind>0 ) then
            if ( iskip(O2Ind) ==1 ) nmol_max = nmol - 1 !yma180520
         endif

         FAC = 1.0
         IF (IPATH(L).EQ.2) FAC = 2.0

         DO 80 K = 1, NMOL

            !---If n_zero==1, zero the layer amount if <0.1%
            IF (n_zero.ne.1) go to 70  !yma IF (n_zero.ne.2) go to 70

            IF (ISKIP(K).NE.1) THEN
               !IF (K.EQ.7 .OR. (IEMIT.EQ.1.AND.IPATH(L).NE.3)) GO TO 70
               !IF ( (nLnMol>=7.and.K.EQ.7) .OR. (IEMIT.EQ.1.AND.IPATH(L).NE.3)) GO TO 70
               if ( (O2Ind>0 .and. K==O2Ind) .or. (IPATH(L).NE.3) ) GO TO 70
               IF ( ((AMTTOT(K)-AMTCUM(K))/AMTTOT(K)).GT.0.001 ) GO TO 70
            ENDIF

            ! ZERO OUT THIS AMOUNT
            ISKIP(K) = 1
            AMOUNT(K,L) = 0.0
            ISKPT = ISKPT+1

            ! IF ALL BUT O2 ARE ZEROED, ELIMINATE ALL HIGHER LAYERS
            IF (ISKPT.GE.(NMOL_max)) GO TO 100

   70       CONTINUE  !Do not zero the layer amount
            AMTCUM(K) = AMTCUM(K)+FAC*AMOUNT(K,L)
   80    CONTINUE

         LMAX = L
   90 END DO

  100 CONTINUE
      IOUTMX = LMAX+1


      deallocate( AMTTOT ) !yma160504
      deallocate( AMTCUM ) !yma160504
      deallocate( ISKIP ) !yma160504

      RETURN

  110 CONTINUE
      WRITE (IPR,900) IOUT,IOUTMX

      STOP ' ERROR FPACK '

  900 FORMAT ('0FROM FPACK-  ERROR, IOUT = ',I5,'  DOES NOT MATCH ',    &
     &        'IOUTMX = ',I5)

      END SUBROUTINE




!-----------------------------------------------------------------------
!     FUNCTION ANDEXD ()
!
!     Double precision version of ANDEX - needed for improved geometry
!
!     COMPUTES THE INDEX OF REFRACTION AT HEIGHT H, SH IS THE
!     SCALE HEIGHT, GAMMA IS THE VALUE AT H=0 OF THE REFRACTIVITY =
!     INDEX-1
!-----------------------------------------------------------------------
      FUNCTION ANDEXD( H,SH,GAMMA )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real(r8)             :: ANDEXD
      real     ,intent(in) :: H
      real(r8) ,intent(in) :: SH,GAMMA


      IF (SH.EQ.0.0) THEN
         ANDEXD = 1.0+GAMMA
      ELSE
         ANDEXD = 1.0+GAMMA*EXP(-H/SH)
      ENDIF

      RETURN
      END  FUNCTION
!
!-----------------------------------------------------------------------
!     FUNCTION RADRFD ()
!
!     Double precision version of RADREF - needed for improved geometry
!
!     COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FOR
!     A HORIZONTAL PATH.  RADREF = ANDEX/ D(ANDEX)/D(RADIUS)
!-----------------------------------------------------------------------
      FUNCTION RADRFD (H,SH,GAMMA)
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real(r8)             :: RADRFD
      real     ,intent(in) :: H
      real(r8) ,intent(in) :: SH,GAMMA

      real(r8) ,PARAMETER :: BIGNUM=1.0D36


      IF (SH.EQ.0.0) GO TO 10
      RADRFD = SH*(1.0+EXP(H/SH)/GAMMA)

      RETURN

   10 RADRFD = BIGNUM

      RETURN

      END FUNCTION

!-----------------------------------------------------------------------
!     FUNCTION ANDEX ()
!
!     COMPUTES THE INDEX OF REFRACTION AT HEIGHT H, SH IS THE
!     SCALE HEIGHT, GAMMA IS THE VALUE AT H=0 OF THE REFRACTIVITY =
!     INDEX-1
!-----------------------------------------------------------------------
      FUNCTION ANDEX( H,SH,GAMMA )
!-----------------------------------------------------------------------
      IMPLICIT NONE
      real             :: ANDEX
      real ,intent(in) :: H,SH,GAMMA

      IF (SH.EQ.0.0) THEN
         ANDEX = 1.0+GAMMA
      ELSE
         ANDEX = 1.0+GAMMA*EXP(-H/SH)
      ENDIF

      RETURN
      END  FUNCTION

!-----------------------------------------------------------------------
!     FUNCTION RADREF ()
!
!     COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FOR
!     A HORIZONTAL PATH.  RADREF = ANDEX/ D(ANDEX)/D(RADIUS)
!-----------------------------------------------------------------------
      FUNCTION RADREF (H,SH,GAMMA)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      real             :: RADREF
      real ,intent(in) :: H,SH,GAMMA

      real, PARAMETER :: BIGNUM=1.0E36


      IF (SH.EQ.0.0) GO TO 10
      RADREF = SH*(1.0+EXP(H/SH)/GAMMA)

      RETURN

   10 RADREF = BIGNUM

      RETURN

      END   FUNCTION


!-----------------------------------------------------------------------
!     SUBROUTINE FNDSHD ()
!
!     Double precision version of FINDSH - needed for improved geometry
!
!     GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARIES
!     Z(I1) AND Z(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE
!     HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE
!     REFRACTIVITY (INDEX OF REFRACTION -1)
!-----------------------------------------------------------------------
      SUBROUTINE FNDSHD( H,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
!-----------------------------------------------------------------------
      USE Module_ConstParam, ONLY: r8=>kind_r8
      IMPLICIT NONE

      real     ,intent(in)  :: H
      real(r8) ,intent(out) :: SH
      real(r8) ,intent(out) :: GAMMA
      real     ,intent(in)  :: ZMDL(:)
      real     ,intent(in)  :: RFNDXM(:)
      integer  ,intent(in)  :: IMMAX

      integer :: I1,I2,IM


      DO 10 IM = 2, IMMAX
         I2 = IM
         IF (ZMDL(IM).GE.H) GO TO 20
   10 END DO
      I2 = IMMAX
   20 CONTINUE
      I1 = I2-1
      CALL SCLHTD (ZMDL(I1),ZMDL(I2),RFNDXM(I1),RFNDXM(I2),SH,GAMMA)

      RETURN
      END  SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE FINDSH ()
!
!     GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARIES
!     Z(I1) AND Z(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE
!     HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE
!     REFRACTIVITY (INDEX OF REFRACTION -1)
!-----------------------------------------------------------------------
      SUBROUTINE FINDSH( H,SH,GAMMA, ZMDL,RFNDXM,IMMAX )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(in)  :: H
      real    ,intent(out) :: SH
      real    ,intent(out) :: GAMMA
      real    ,intent(in)  :: ZMDL(:)
      real    ,intent(in)  :: RFNDXM(:)
      integer ,intent(in)  :: IMMAX

      integer :: I1,I2, IM


      DO 10 IM = 2, IMMAX
         I2 = IM
         IF (ZMDL(IM).GE.H) GO TO 20
   10 END DO
      I2 = IMMAX
   20 CONTINUE
      I1 = I2-1
      CALL SCALHT (ZMDL(I1),ZMDL(I2),RFNDXM(I1),RFNDXM(I2),SH,GAMMA)

      RETURN
      END  SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE SCLHTD ()
!
!     Double precision version of SCALHT - needed for improved geometry
!
!     THIS SUBROUTINE CALCULATES THE SCALE HEIGHT SH OF THE (INDEX OF
!     REFRACTION-1.0) FROM THE VALUES OF THE INDEX AT THE ALTITUDES Z1
!     AND Z2 ( Z1 < Z2). IT ALSO CALCULATES THE EXTRAPOLATED VALUE
!     GAMMA OF THE (INDEX-1.0) AT Z = 0.0
!-----------------------------------------------------------------------
      SUBROUTINE SCLHTD( Z1,Z2,RFNDX1,RFNDX2,SH,GAMMA )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real     ,intent(in)  :: Z1,Z2
      real     ,intent(in)  :: RFNDX1,RFNDX2
      REAL(r8) ,intent(out) :: SH,GAMMA

      real :: RATIO,RF1,RF2


      RF1 = RFNDX1+1.0E-20
      RF2 = RFNDX2+1.0E-20
      RATIO = RF1/RF2
      IF (ABS(RATIO-1.0).LT.1.0E-05) GO TO 10
      SH = (Z2-Z1)/ LOG(RATIO)
      GAMMA = RF1*(RF2/RF1)**(-Z1/(Z2-Z1))
      GO TO 20
   10 CONTINUE

      !     THE VARIATION IN THE INDEX OF REFRACTION WITH HEIGHT IS
      !     INSIGNIFICANT OR ZERO
      SH = 0.0
      GAMMA = RFNDX1
   20 CONTINUE

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE SCALHT ()
!
!     THIS SUBROUTINE CALCULATES THE SCALE HEIGHT SH OF THE (INDEX OF
!     REFRACTION-1.0) FROM THE VALUES OF THE INDEX AT THE ALTITUDES Z1
!     AND Z2 ( Z1 < Z2). IT ALSO CALCULATES THE EXTRAPOLATED VALUE
!     GAMMA OF THE (INDEX-1.0) AT Z = 0.0
!-----------------------------------------------------------------------
      SUBROUTINE SCALHT( Z1,Z2,RFNDX1,RFNDX2,SH,GAMMA )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real ,intent(in)  :: Z1,Z2
      real ,intent(in)  :: RFNDX1, RFNDX2
      real ,intent(out) :: SH,GAMMA

      real :: RATIO,RF1,RF2

      RF1 = RFNDX1+1.0E-20
      RF2 = RFNDX2+1.0E-20
      RATIO = RF1/RF2
      IF (ABS(RATIO-1.0).LT.1.0E-05) GO TO 10
      SH = (Z2-Z1)/ LOG(RATIO)
      GAMMA = RF1*(RF2/RF1)**(-Z1/(Z2-Z1))
      GO TO 20
   10 CONTINUE

      ! THE VARIATION IN THE INDEX OF REFRACTION WITH HEIGHT IS
      ! INSIGNIFICANT OR ZERO

      SH = 0.0
      GAMMA = RFNDX1
   20 CONTINUE

      RETURN
      END  SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE autoLayering( prfl, pathCtrl, zRT )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: MXFSC, IBDIM=>MXFSC
      USE Module_Scene       ,ONLY: CLBLM_Profile
      USE Module_Config      ,ONLY: CLBLM_Path_Ctrl
      IMPLICIT NONE

      type(CLBLM_Profile)   ,intent(in)  :: prfl
      type(CLBLM_Path_Ctrl) ,intent(in)  :: pathCtrl
      real     ,allocatable ,intent(out) :: zRT(:)


      integer  :: IERROR
      integer  :: nLev, IBMAX
      real     :: XVBAR
      real     :: ZMAX, HMIN, HMAX
      real     :: AVTRAT, TDIFF1, TDIFF2, ALTD1, ALTD2
      real     :: ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC)
      real     :: ALORNZ(MXFSC),ADOPP(MXFSC),AVOIGT(MXFSC)


      nLev   = prfl%nLev
      ZMAX   = prfl%Z( nLev )
      HMIN   = prfl%Z( 1 )
      HMAX   = ZMAX
      XVBAR  = pathCtrl%V_refrac
      AVTRAT = pathCtrl%maxAlphaRat
      TDIFF1 = pathCtrl%maxTmpDiff(1)
      TDIFF2 = pathCtrl%maxTmpDiff(2)
      ALTD1  = pathCtrl%maxTmpRefAlt(1)
      ALTD2  = pathCtrl%maxTmpRefAlt(2)

      call AUTLAY( HMIN,HMAX,XVBAR,AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2, &  !in
                   prfl%Z(1:nLev),prfl%P(1:nLev),prfl%T(1:nLev),ZMAX, & !in
                   ZBND,PBND,TBND,ALORNZ,ADOPP,AVOIGT,IBMAX,&           !out
                   nLev, IERROR )                                       !in
      if (IERROR/=0) STOP '--- autoLayering(): Error in auto layering, the number of generated layer boundaries exceeds the maximum dimension.'


      allocate( zRT(IBMAX) )
      zRT(1:IBMAX) = ZBND(1:IBMAX)

   END SUBROUTINE

!-----------------------------------------------------------------------
!     THIS SUBROUTINE AUTOMATICALLY SELECTS A SET OF LBLRTM BOUNDARY
!     LEVELS WHICH SATISFY THE FOLLOWING TWO TESTS:
!          1. THE RATIO OF THE VOIGT HALFWIDTHS BETWEEN BOUNDARIES
!             IS LESS THAN OR EQUAL TO AVTRAT, AND
!          2. THE TEMPERATURE DIFFERENCE BETWEEN BOUNDARIES IS
!             LESS THAN OR EQUAL TO TDIFF
!     TDIFF VARIES FROM TDIFF1 AT HMIN TO TDIFF2 AT HMAX,
!     WITH EXPONENTIAL INTERPOLATION BETWEEN
!     THESE BOUNDARIES ARE ROUNDED DOWN TO THE NEAREST TENTH KM
!     NOTE THAT THESE TESTS APPLY TO THE LAYER BOUNDARIES
!     NOT TO THE AVERAGE VALUES FROM ONE LAYER TO THE NEXT.
!-----------------------------------------------------------------------
   SUBROUTINE AUTLAY( minH,maxH,XVBAR,AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2, &
                      ZMDL,PM,TM,ZMAX, &
                      ZBND,PBND,TBND,ALORNZ,ADOPP,AVOIGT,IBMAX, &
                      IMMAX, IERROR )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: MXZMD, IBDIM=>MXFSC
      USE Module_Config     ,ONLY: IPR
      USE Module_Utility    ,ONLY: EXPINT

      real     ,intent(in)  :: minH, maxH
      real     ,intent(in)  :: XVBAR
      real     ,intent(in)  :: AVTRAT
      real     ,intent(in)  :: TDIFF1, TDIFF2
      real     ,intent(in)  :: ALTD1, ALTD2
      integer  ,intent(out) :: IERROR
      !
      real     ,intent(in)  :: ZMDL(:), PM(:), TM(:)
      real     ,intent(in)  :: ZMAX
      real     ,intent(out) :: ZBND(:), PBND(:), TBND(:)
      real     ,intent(out) :: ALORNZ(:), ADOPP(:), AVOIGT(:)
      integer  ,intent(out) :: IBMAX
      integer  ,intent(in)  :: IMMAX


      integer  :: IM, IM2, IB, IBM1, IND, IHMIN, IPASS
      real     :: ZROUND,ZX
      real     :: HMIN, HMAX, HTOP, ZZ
      real     :: P,T,AL,AD, AVTM(MXZMD)
      real     :: TMIN,TMAX
      real     :: ZBNDTI, FAC, TDIFF
      real     :: X, ALOGX, Y, ALOGY



      !--- Statement functioin ZROUND ROUNDS THE ALTITUDE Z DOWN TO THE NEAREST TENTH KM
      ZROUND(ZX) = 0.1* REAL( INT(10.0*ZX))

      IERROR = 0
      HMIN = minH
      HMAX = maxH

      HMIN = MAX(HMIN,ZMDL(1))

      DO 10 IM = 2, IMMAX
         IHMIN = IM
         IF (ZMDL(IM).GT.HMIN) GO TO 20
   10 END DO
   20 CONTINUE
      HTOP = HMAX
      HTOP = MIN(HTOP,ZMAX)
      IM = IHMIN-1
      ZZ = ZMDL(IM)
      CALL HALFWD (ZZ,XVBAR,P,T,AL,AD,AVTM(IM), ZMDL,PM,TM,IMMAX)
      IB = 1
      ZBND(IB) = HMIN
      IM = IHMIN
      CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),         &
     &             ADOPP(IB),AVOIGT(IB), ZMDL,PM,TM,IMMAX)

      !--- BEGIN IB LOOP. Loop over boundary levels
      !
   30 CONTINUE
         IB = IB+1
         IF (IB.GT.IBDIM) GO TO 90
         IBM1 = IB-1
         TMIN = TBND(IBM1)
         TMAX = TBND(IBM1)
         IND = 0

         !--- BEGIN IM LOOP. Loop over profile levels
         !
   40    CONTINUE
               IPASS = 0
               ZBND(IB) = ZMDL(IM)
               ZBNDTI = ZMDL(IM)
               IF (ZBND(IB).GE.HTOP) ZBND(IB) = HTOP
               CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),  &
                            ADOPP(IB),AVOIGT(IB), ZMDL,PM,TM,IMMAX)
               AVTM(IM) = AVOIGT(IB)

               !--- TEST THE RATIO OF THE VOIGT WIDTHS AGAINST AVTRAT
               IF ((AVOIGT(IB-1)/AVOIGT(IB)).LT.AVTRAT) GO TO 50

            !--- ZMDL(IM) FAILS THE HALFWIDTH RATIO TEST
            IPASS = 1
            AVOIGT(IB) = AVOIGT(IB-1)/AVTRAT
            X = AVTM(IM)/AVTM(IM-1)
            ALOGX = 1.-X
            IF (ABS(ALOGX).LT.0.001) THEN
               ZBND(IB) = (ZMDL(IM)+ZMDL(IM-1))/2.
               GO TO 50
            ELSE
               ALOGX = LOG(X)
            ENDIF
            Y = AVOIGT(IB)/AVTM(IM-1)
            ALOGY = 1.-Y
            IF (ABS(ALOGY).GT.0.001) ALOGY =  LOG(Y)
            ZBND(IB) = ZMDL(IM-1)+(ZMDL(IM)-ZMDL(IM-1))*ALOGY/ALOGX

   50          CONTINUE
               !--- TEST THE TEMPERATURE DIFFERENCE AGAINST TDIFF
               FAC = (ZBND(IB-1)-ALTD1)/(ALTD2-ALTD1)
               CALL EXPINT (TDIFF,TDIFF1,TDIFF2,FAC)
               IF (TM(IM).GT.TMAX) THEN
                  IND = 1
                  TMAX = TM(IM)
               ENDIF
               IF (TM(IM).LT.TMIN) THEN
                  IND = 2
                  TMIN = TM(IM)
               ENDIF

               IF (TMAX-TMIN.LE.TDIFF) GO TO 60

            IF (IND.EQ.1) TBND(IB) = TMIN+TDIFF
            IF (IND.EQ.2) TBND(IB) = TMAX-TDIFF

            !--- ZBND(IB) FAILS THE TEMPERATURE DIFFERENCE TEST
            IPASS = 2
            IF (ABS(TM(IM)-TM(IM-1)).LT.0.0001) THEN
               ZBNDTI = (ZMDL(IM)+ZMDL(IM-1))/2.
            ELSE
               ZBNDTI = ZMDL(IM-1)+(ZMDL(IM)-ZMDL(IM-1))* (TBND(IB)-TM(IM-1))/&
               (TM(IM)-TM(IM-1))
            ENDIF

   60          CONTINUE
               IF (ZBNDTI.LT.ZBND(IB)) ZBND(IB) = ZBNDTI
               IF (ZBND(IB).GE.HTOP) THEN
                  ZBND(IB) = HTOP
                  IF (ZBND(IB)-ZBND(IB-1).LE.0.1) THEN
                     IB = IB-1
                     ZBND(IB) = HTOP
                     CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),   &
                                  ADOPP(IB),AVOIGT(IB), ZMDL,PM,TM,IMMAX)
                  ENDIF
                  GO TO 80
               ENDIF
               IF (IPASS.NE.0) GO TO 70

               !--- BOTH HALFWIDTH AND TEMPERATURE TEST PASS FOR ZBND(IB) = ZMDL(IM),
               !  NOW TRY ZBND(IB) = ZMDL(IM+1)
               IM = IM+1
         GO TO 40
   70    CONTINUE

         !--- ONE OF THE TESTS FAILED AND A NEW BOUNDRY ZBND WAS PRODUCED
         ZBND(IB) = ZROUND(ZBND(IB))
         CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB), &
                      ADOPP(IB),AVOIGT(IB), ZMDL,PM,TM,IMMAX)
      GO TO 30
   80 CONTINUE
      IBMAX = IB
      WRITE (IPR,900) AVTRAT,TDIFF1,HMIN,TDIFF2,HMAX

      RETURN

   90 CONTINUE
      WRITE (IPR,905) IBDIM
      IBMAX = IBDIM
      IERROR = 5

      RETURN

  900 FORMAT (///,                                                      &
     &        ' LBLRTM LAYER BOUNDARIES PRODUCED BY THE AUTOMATIC ',    &
     &        'LAYERING ROUTINE AUTLAY',/,' THE USER SHOULD EXAMINE ',  &
     &        'THESE BOUNDARIES AND MODIFY THEM IF APPROPRIATE',/,      &
     &        ' THE FOLLOWING PARAMETERS ARE USED:',//,10X,             &
     &        'AVTRAT    = ',F8.2,'       = MAX RATIO OF VOIGT WIDTHS', &
     &        /,10X,'TDIFF1    = ',F8.2,'       = MAX TEMP DIFF AT ',   &
     &        F4.0,' KM',/10X,'TDIFF2    = ',F8.2,                      &
     &        '       = MAX TEMP DIFF AT ',F4.0,' KM')
  905 FORMAT (///,' ERROR IN AUTLAY:',/,5X,'THE NUMBER OF ',            &
     &        'GENERATED LAYER BOUNDARIES EXCEEDS THE DIMENSION IBDIM', &
     &        ' OF THE ARRAY ZBND.  IBDIM = ',I5,/,5X,'PROBABLE CAUSE', &
     &        ': EITHER AVTRAT AND/OF TDIFF ARE TOO SMALL',/,5X,        &
     &        'THE GENERATED LAYERS FOLLOW')

      END SUBROUTINE

!-----------------------------------------------------------------------
!     GIVEN AN ALTITUDE Z AND AN AVERAGE WAVENUMBER VBAR, THIS
!     SUBROUTINE INTERPOLATES P AND T FROM THE PROFILE IN ZMDL  AND
!     CALCULATES THE LORENTZ, THE DOPPLER, AND THE VOIGT HALFWIDTHS
!     (AT HALFHEIGHT) ALORNZ, ADOPP, AND AVOIGT RESPECTIVELY FOR
!     THE ALTITUDE Z
!     AN AVERAGE LORENTZ WIDTH ALZERO AND AN AVERAGE MOLECULAR
!     WEIGHT AVMWT ARE ASSUMED
!-----------------------------------------------------------------------
   SUBROUTINE HALFWD (Z,XVBAR,P,T,ALORNZ,ADOPP,AVOIGT, ZMDL,PM,TM,IMMAX)
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: ALZERO, ADCON, AVMWT, &
                                   PZERO=>Press1013, TZERO=>Temp273
      USE Module_Utility    ,ONLY: EXPINT

      real    ,intent(in)  :: Z
      real    ,intent(in)  :: XVBAR
      real    ,intent(out) :: P,T
      real    ,intent(out) :: ALORNZ
      real    ,intent(out) :: ADOPP
      real    ,intent(out) :: AVOIGT
      !
      real    ,intent(in)  :: ZMDL(:),PM(:),TM(:)
      integer ,intent(in)  :: IMMAX

      integer :: I2, IM
      real    :: ALPHAL,ALPHAD,ALPHAV
      real    :: FAC,AL,AD,V


      !--- Statement functions
      ! * ALZERO IS AT 1013.25 MB AND 296.0 K
      ALPHAL(P,T) = ALZERO*(P/PZERO)*SQRT(296.0/T)
      ALPHAD(T,V) = ADCON*V*SQRT(T/AVMWT)
      ALPHAV(AL,AD) = 0.5*(AL+SQRT(AL**2+4.0*AD**2))

      DO 10 I2 = 2, IMMAX
         IM = I2
         IF (ZMDL(IM).GE.Z) GO TO 20
   10 END DO
      IM = IMMAX
   20 CONTINUE
      FAC = (Z-ZMDL(IM-1))/(ZMDL(IM)-ZMDL(IM-1))
      CALL EXPINT (P,PM(IM-1),PM(IM),FAC)
      T = TM(IM-1)+(TM(IM)-TM(IM-1))*FAC
      ALORNZ = ALPHAL(P,T)
      ADOPP = ALPHAD(T,XVBAR)
      AVOIGT = ALPHAV(ALORNZ,ADOPP)

      RETURN
      END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE insertObsAltInRTgrid( zRT,nRTlev, Hobs )
!-----------------------------------------------------------------------
      real ,allocatable,intent(inout) :: zRT(:)
      integer          ,intent(inout) :: nRTlev
      real             ,intent(in)    :: Hobs

      real ,parameter    :: TOL = 5.e-4
      real ,allocatable  :: newRTgrid(:)
      integer            :: il
      real               :: f


      !--- Insert Hobs into RTgrid
      if ( Hobs > zRT(1) .and. Hobs < zRT( nRTlev ) ) then

         do il = 1,nRTlev

            if ( abs( zRT(il)-Hobs ) < TOL ) then !There is a level that is very close to Hobs in the array.
               zRT(il) = Hobs
               EXIT
            endif

            if ( Hobs < zRT(il) ) then

               allocate( newRTgrid( nRTlev+1 ))

               !--- Insert Z
               newRTgrid(1:il-1)        = zRT(1:il-1)
               newRTgrid(il)            = Hobs
               newRTgrid(il+1:nRTlev+1) = zRT(il:nRTlev)

               !!--- Interpolate P
               !newRTgrid%P(1:il-1)              = RTgrid%P(1:il-1)
               !newRTgrid%P(il+1:newRTgrid%nLev) = RTgrid%P(il:RTgrid%nLev)
               !f = (Hobs - RTgrid%Z(il-1)) / (RTgrid%Z(il) - RTgrid%Z(il-1))
               !CALL EXPINT (newRTgrid%P(il), RTgrid%P(il-1), RTgrid%P(il), f )

               !--- Replace the old RTgrid and exit the loop.
               call move_alloc( newRTgrid, zRT )
               nRTlev = nRTlev+1
               EXIT
            endif
         enddo
      endif

   END SUBROUTINE



END MODULE
