!
! CREATION HISTORY:
!       Modified from LBLRTM v12.9
!       Yingtao Ma, AER@NOAA/NESDIS/STAR
!       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_ODLAY

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: ODLAY


   INTERFACE ODLAY
      module procedure ODLAY_struct
      module procedure ODLAY_nonStruct
   END INTERFACE

CONTAINS !===================== Module Contains ========================


!-----------------------------------------------------------------------
! A wrapper subroutine to call non-structure arguments version of ODLAY
!-----------------------------------------------------------------------
   SUBROUTINE ODLAY_struct( aLayer, odCtrl, spectGrid, dvCtrl, &
                            layOD, NLTE,layNLTEEmisFac, JacMolName )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8
      USE Module_Config         ,ONLY: CLBLM_OD_Ctrl,&
                                       CLBLM_DV_Ctrl
      USE Module_DV             ,ONLY: CLBLM_SpectGrid
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum, &
                                       CLBLM_Spectrum_init
      USE Module_AtmPath        ,ONLY: CLBLM_Layer
      IMPLICIT NONE

      type(CLBLM_Layer)               ,intent(in)  :: aLayer
      type(CLBLM_OD_Ctrl)             ,intent(in)  :: odCtrl
      type(CLBLM_SpectGrid)           ,intent(in)  :: spectGrid
      type(CLBLM_DV_Ctrl)             ,intent(in)  :: dvCtrl
      type(CLBLM_Spectrum)            ,intent(out) :: layOD          !optical depth spectrum range from V1 to V2
      logical               ,optional ,intent(in)  :: NLTE
      type(CLBLM_Spectrum)  ,optional ,intent(out) :: layNLTEEmisFac !spectrum of the NLTE emission correction factor
      character(*)          ,optional ,intent(in)  :: JacMolName     !molecular name for radiance derivative


      logical                   :: NLTE_flag
      integer                   :: nMol
      character(20),allocatable :: molID(:)
      real                      :: Pav
      real                      :: Tav
      real         ,allocatable :: W(:) ![nMol]
      real                      :: Wtot
      real                      :: Zbot
      real                      :: Ztop
      real(r8)                  :: V1,V2
      real                      :: DVnormal
      real                      :: DVnarrow
      real                      :: SAMPLE
      real                      :: ALFAL0
      character(20)             :: ajMolName
      integer                   :: IPATH
      integer                   :: layNo
      integer                   :: nSamp
      real         ,allocatable :: monoOD(:)
      real         ,allocatable :: nlteEmisFac(:)
      real                      :: DVout


      allocate( molID( aLayer%nMol) )
      allocate( W(     aLayer%nMol) )

      !--- Arguments to array-version of ODLAY
      nMol        = aLayer%nMol
      molID(:)    = aLayer%molID
      Pav         = aLayer%Pave
      Tav         = aLayer%Tave
      W(:)        = aLayer%W
      Wtot        = aLayer%Wtot
      Zbot        = aLayer%Zbot
      Ztop        = aLayer%Ztop
      !odCtrl_arg  = odCtrl
      V1          = spectGrid%V1
      V2          = spectGrid%V2
      DVnormal    = spectGrid%DVnormal( aLayer%layNo )
      DVnarrow    = spectGrid%DVnarrow( aLayer%layNo )
      SAMPLE      = dvCtrl%SAMPLE
      ALFAL0      = dvCtrl%ALFAL0
      layNo       = aLayer%layNo
      !nSamp       =
      !monoOD      =
      !NLTE        =
      !nlteEmisFac =
      DVout       = spectGrid%DVout
      ajMolName   = JacMolName


      !--- Call array-argument version of ODLAY
      CALL ODLAY_nonStruct( nMol, molID, Pav, Tav, W, Wtot, Zbot,Ztop, odCtrl, &
                            V1, V2, DVnormal, DVnarrow, &
                            SAMPLE, ALFAL0, IPATH, LayNo, &
                            nSamp, monoOD, NLTE,nlteEmisFac, DVout, ajMolName )


      call CLBLM_Spectrum_init( layOD, V1,DVout,nSamp )
      if (size(monoOD)==nSamp) then
         call move_alloc( monoOD, layOD%spect )
      else
         layOD%spect(1:nSamp) = monoOD(1:nSamp)
      endif

      if (present(NLTE)) then
      if (NLTE) then
         call CLBLM_Spectrum_init( layNLTEEmisFac, V1,DVout,nSamp )
         if (size(nlteEmisFac)==nSamp) then
            call move_alloc( nlteEmisFac, layNLTEEmisFac%spect )
         else
            layNLTEEmisFac%spect(1:nSamp) = nlteEmisFac(1:nSamp)
         endif
      endif
      endif

   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE ODLAY_nonStruct( nMol, molID, Pav, Tav, W, Wtot, Zbot,Ztop, odCtrl, &
                               V1, V2, DVnormal, DVnarrow, &
                               SAMPLE, ALFAL0, IPATH, LayNo, &
                               nSamp, monoOD, NLTE,nlteEmisFac, DVout, ajMolName )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8, FILLINT,&
                                       MXMOL,MXISOTPL,MX_XS, &
                                       ISOTPL_ABD, &
                                       isoName2molNum, isoName2isoNum, &
                                       molNum, molIndex
      USE Module_Utility        ,ONLY: getLun
      USE Module_Config         ,ONLY: CLBLM_OD_Ctrl,&
                                       CLBLM_DV_Ctrl,&
                                       ioFiles, IPR
      USE Module_DV             ,ONLY: CLBLM_SpectGrid
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum, &
                                       CLBLM_Spectrum_init, &
                                       interpSpectrum
      USE Module_FileIO         ,ONLY: NWDL
      USE Module_AtmPath        ,ONLY: CLBLM_Layer
      USE Module_LineData       ,ONLY: CLBLM_LineData
      USE Module_XSect          ,ONLY: xsTbl,&
                                       getXsInfo, &
                                       readXsMasterFile
      IMPLICIT NONE

      integer                     ,intent(in)    :: nMol
      character(*)                ,intent(in)    :: molID(:)
      real                        ,intent(in)    :: Pav
      real                        ,intent(in)    :: Tav
      real                        ,intent(in)    :: W(:) ![nMol]
      real                        ,intent(in)    :: Wtot
      real                        ,intent(in)    :: Zbot
      real                        ,intent(in)    :: Ztop
      type(CLBLM_OD_Ctrl)         ,intent(in)    :: odCtrl
      real(r8)                    ,intent(in)    :: V1,V2
      real                        ,intent(in)    :: DVnormal
      real                        ,intent(in)    :: DVnarrow
      real                        ,intent(in)    :: SAMPLE
      real                        ,intent(in)    :: ALFAL0
      integer                     ,intent(in)    :: IPATH
      integer                     ,intent(in)    :: layNo
      integer                     ,intent(out)   :: nSamp
      real ,allocatable           ,intent(out)   :: monoOD(:)
      logical           ,optional ,intent(in)    :: NLTE
      real ,allocatable ,optional ,intent(out)   :: nlteEmisFac(:)
      real              ,optional ,intent(inout) :: DVout    !output the real DV used for OD calculation
      character(*)      ,optional ,intent(in)    :: ajMolName



      !--- Local variables
      !
logical, SAVE :: firstPassJRAD=.TRUE. !this is to make the setting of JRAD to be consistent with LBLRTM
integer, SAVE :: JRAD
      character(*) ,parameter :: routineName = 'ODLAY'

      real, PARAMETER :: HZ=6.5 !scale height

      integer  :: i,j,k,iv, iStat
      logical  :: OP,exist
      integer  :: lineFileLun, lineF4FileLun, nlteStatPopFileLun, narrowLineFile
      real     :: outputDV
      real     :: ALFAV
!jrad      integer  :: JRAD
      logical  :: NLTE_flag

      integer                              :: nLnMol, nXsMol, nISO, iMol, iISO
      character( len(molID) ) ,allocatable :: lnMolNames(:)
      character( len(molID) ) ,allocatable :: xsMolNames(:)
      character( len(molID) ) ,allocatable :: isoMolNames(:)

      character*10 :: XSFILE(6,5,MX_XS)
      integer      :: IXFORM(5,MX_XS),NSPECR(MX_XS),NTEMPF(5,MX_XS), NUMXS
      real         :: XDOPLR(5,MX_XS)
      integer      :: LnOrXs(5,MX_XS)
      integer      :: ILBLF4,IXSCNT,IPATHL,IAERSL
      integer      :: JCNVF4
      real         :: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
      integer      :: ILNFLG
      real         :: DPTMIN,DPTFAC,ALTAV !WTOT,SAMPLE,DPTMIN,DPTFAC,ALTAV
      real         :: WKI(MXMOL,MXISOTPL)
      integer      :: ISOTPL_FLAG(MXMOL,MXISOTPL)
      integer      :: IXSBIN
      real         :: WXM(MX_XS)
      integer      :: NOPR,LINFIL,LNFIL4,NFHDRF,NPHDRF,NPHDRL,NLNGTH
      integer      :: icflg
      logical      :: speciesBroad
      logical      :: XsectOff
      logical      :: contnmOff

      real    ,ALLOCATABLE :: r1BuffNorm(:)
      real    ,ALLOCATABLE :: r1BuffFine(:)
      real    ,ALLOCATABLE :: rr1BuffNorm(:)
      real    ,ALLOCATABLE :: rr1BuffFine(:)
      real    ,ALLOCATABLE :: tempBuff(:)
      integer              :: sizeOfBuffer
      integer              :: numNormPoints, numFinePoints
      integer              :: lengthOfSpect
      real                 :: DVnorm
      integer              :: numNarrowLines
      real                 :: narrowWidth
! MJI - add small epsilon value for sizeOfBuffer check
      real, parameter      :: v2eps = 1.e-7



      !--- Local common blocks.
      ! Used only for determining the size of data records in memory
      !
      REAL        :: PAVE,TAVE,WK,PZL,PZU,TZL,TZU,WBROAD,DV,TBOUND,EMISIV,FSCDID,YI1
      REAL*8      :: SECANT,XALTZ,V1com,V2com
      INTEGER     :: NMOLcom,LAYER,LSTWDF
      CHARACTER*8 :: XID,HMOLID,YID
      COMMON /FILHDR_0/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4), &
                      WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1com ,V2com ,TBOUND, &
                      EMISIV,FSCDID(17),NMOLcom,LAYER ,YI1,YID(10),LSTWDF

      REAL    :: DVD
      REAL*8  :: V1D,V2D
      INTEGER :: NLND,IWLD
      COMMON /HDRF_0/ V1D,V2D,DVD,NLND,IWLD !used for counting the record length

      REAL    :: SD,AD,EPD,HWHD,TMPD,PSHD,FLGD
      REAL*8  :: VD
      INTEGER :: MOLD,ILS2D
      COMMON /NGTH_0/ VD,SD,AD,EPD,MOLD,HWHD,TMPD,PSHD,FLGD,ILS2D  !used for counting the record length

      REAL*8   :: V1LD,VL2D
      INTEGER  :: NLD,NWDS,ILST3D
      COMMON /HDRL_0/ V1LD,VL2D,NLD,NWDS,ILST3D

      INTEGER :: IWD(2),IWD2(2),IWD3(2),IWD4(2)  !used for counting the record length
      EQUIVALENCE (IWD(1),XID(1)),(IWD2(1),V1D),(IWD3(1),VD),(IWD4(1),V1LD)  !used for counting the record length



   !==================== Executable code start here ====================

      !--- Some options are now fixed in CLBLM
      ILBLF4 = 1
      NOPR   = 0
      IAERSL = 0 !used in continuum
      NLTE_flag = .FALSE.
      if (present(NLTE)) NLTE_flag = NLTE

      !--- Open line data file (TAPE3)
      inquire( FILE=trim(ioFiles%lineFile), EXIST=exist, opened=op, number=lineFileLun )
      if (.not.exist)  STOP '--- '//routineName//'(): Line data file not exist, program stopped.'
      if (.not.op) then
         lineFileLun = getlun()
         OPEN (lineFileLun,FILE=trim(ioFiles%lineFile),STATUS='OLD',FORM='UNFORMATTED')
      endif
      ! Rewind the line data file
      REWIND( lineFileLun )


      !--- Read the XS master file
      if ( xsTbl%numXS ==0 ) then !if table is not yet readed
         CALL readXsMasterFile( trim(ioFiles%xsFilePath)//'FSCDXS' )
      endif


      !--- Temporary file to save the shrunk lines for line function 4 (TAPE9)
      if (ILBLF4 >0) then
         inquire( file=trim(ioFiles%lineF4File), opened=OP, number=lineF4FileLun )
         if (.NOT.OP) then
            lineF4FileLun = getLun()
            OPEN( lineF4FileLun,FILE=trim(ioFiles%lineF4File), &
                                STATUS='UNKNOWN',FORM='UNFORMATTED')
         else
            rewind( lineF4FileLun )
         endif
      endif


      !--- Separate line/XS/isotopologue
      CALL seperateLineXsISO( molID, nMol, &
                              V1, V2,&
                              lnMolNames, xsMolNames, isoMolNames, &
                              nLnMol, nXsMol, nISO )


      ! Call internal subroutine to check with the line data base
      !call internalSub_checkLineData()


      !---If x-section, read xs information
      if (.not.odCtrl%XSectOff) then
         CALL getXsInfo( V1, V2, &
                         nXsMol, xsMolNames, &
                         XSFILE, IXFORM, NSPECR, NTEMPF, XDOPLR, LnOrXs )
      endif


      !--- NLTE vibration state population file (TAPE4)
      if (NLTE_flag) then
         inquire( file=trim(ioFiles%nlteStatPopFile), opened=OP, number=nlteStatPopFileLun )
         if (.NOT.OP) then
            nlteStatPopFileLun = getLun()
            OPEN( nlteStatPopFileLun,FILE=trim(ioFiles%nlteStatPopFile),&
                                     STATUS='OLD')
         else
            rewind( nlteStatPopFileLun )
         endif
      endif


      !---If doing double convolution, open a scratch file to save narrow lines
      !if ( odCtrl%doubleConv ) then
         narrowLineFile = getlun()
         OPEN (narrowLineFile, FILE=trim(ioFiles%narrowLineFile), STATUS='UNKNOWN', FORM='unformatted')
      !endif



      !--- Initialize COMMON BLOCK variables ---
      !
      !EQUIVALENCE FSCDID YI1
      IPATHL = IPATH !for RT and "IF (IEMST.EQ.0.AND.IPATHL.EQ.2) DPTMIN = 2.*DPTMST"
if (firstPassJRAD) then !to be consistent with LBLRTM
      JRAD   = odCtrl%JRAD !may be modified in ODLAY subs
firstPassJRAD=.FALSE.
endif

      !COMMON /CONVF/ !use: (block data popdpt) JCNVF4
      JCNVF4 = 0

      !COMMON /CNTSCL/ !use: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
      XSELF = odCtrl%XSELF
      XFRGN = odCtrl%XFRGN
      XCO2C = odCtrl%XCO2C
      XO3CN = odCtrl%XO3CN
      XO2CN = odCtrl%XO2CN
      XN2CN = odCtrl%XN2CN
      XRAYL = odCtrl%XRAYL

      !COMMON /RCNTRL/ !use: ILNFLG
      !=-1  deativate line rejection (This is added option, not exist in LBLRTM)
      != 0  line rejection information not recorded
      != 1  write line rejection information to REJ1, REJ4
      != 2  read line rejection information from REJ1, REJ4
      if ( odCtrl%lineRejec ) then
         ILNFLG = 0 !line rejection information not recorded
      else
         ILNFLG = -1 !deativate line rejection
      endif

      !COMMON /MANE/ !use: P0,TEMP0,WTOT, SAMPLE,ALFAL0,DPTMIN,DPTFAC,-(opdepth)IPFLAG, ALTAV (for NLTE)
      DPTMIN = odCtrl%DPTMIN
      DPTFAC = odCtrl%DPTFAC
      IF (ZBot >= 0.) THEN
         ALTAV = ZBot - HZ*LOG(.5*(1.+EXP(-(ZTop-ZBot)/HZ)))
      ELSE
         ALTAV = ZTop
      ENDIF


      !COMMON /FILHDR/ !use: XID(10),PAVE,TAVE,*HMOLID(60),WK(60),(continuum)WBROAD,
      !                      DV ,V1 ,V2,FSCDID(17), NMOL,(opdpth)LAYRS,
      !                      (IOD)YI1,YID(10)
      LAYER = layNo
      PAVE  = Pav
      TAVE  = Tav
      DV    = DVnormal !First do normal line convolution. If there are narrow lines exist, a second round convolution will be done.
      WK(:) = 0.;
      do i = 1,nLnMol
         imol = molNum( lnMolNames(i) )
         WK(imol) = W( molIndex(lnMolNames(i), molID) )
      enddo


      !COMMON /PATH_ISOTPL/ use: ISOTPL_FLAG(MXMOL,MXISOTPL), WKI(MXMOL,MXISOTPL)  !v12.7
      WKI(:,:) = 0.
      ISOTPL_FLAG(:,:) = 0
      do i = 1,nISO
         iMOL = isoName2molNum( isoMolNames(i) )
         iISO = isoName2isoNum( isoMolNames(i) )
         !--- Scale the isotopologule amounts by HITRAN abundance ratio.
         WKI(iMOL,iISO) = W( molIndex(isoMolNames(i), molID) ) / ISOTPL_ABD(iMOL,iISO)
         ISOTPL_FLAG(iMOL,iISO) = 1
      enddo


      !COMMON /XSECTR/ and /XSECTF/ !use: WXM,NUMXS,IXSBIN...
      NUMXS  = nXsMol !xsMast%NUMXS
      WXM(:) = 0.
      do i =1,NUMXS
         WXM(i) = W( molIndex(xsMolNames(i),molID) )
      enddo
      IXSBIN = 0
      if (.not.odCtrl%xsPressConv) IXSBIN=1 !IXSBIN = odCtrl%IXSBIN


      !COMMON /IFIL/ !use: (xsread)IRD,IPR,NOPR,NFHDRF,NPHDRF,NLNGTH,
      !                  LINFIL,LNFIL4,speciesBroad
      LINFIL = lineFileLun
      if (ILBLF4 >0) then
         LNFIL4 = lineF4FileLun
      endif
      speciesBroad   = odCtrl%speciesBroad !v12.7

      LSTWDF = -654321 ;NFHDRF = NWDL(IWD,LSTWDF)  !used by XS>
      IWLD   = -654321 ;NPHDRF = NWDL(IWD2,IWLD)   !used by XS>
      ILS2D  = -654321 ;NLNGTH = NWDL(IWD3,ILS2D)  !line record length
      ILST3D = -654321 ;NPHDRL = NWDL(IWD4,ILST3D) !used by LIN4>

      !COMMON /CDERIV/ !use: icflg, !outOnlyButNoUseOutside: v1absc,v2absc,dvabsc,nptabsc
      icflg = FILLINT !=nspcrt !for derivative use in continuum
      if (present(ajMolName)) icflg = molNum( ajMolName )


      XsectOff    = odCtrl%XsectOff
      contnmOff   = odCtrl%contnmOff
      ALFAV       = SAMPLE * max(DVnormal,DVnarrow)

      if ( abs(DVnormal-DVnarrow)/DVnormal >0.2 ) then
         !--- If very much different, do two runs.
         ! In the first pass, check if there are narrow lines exist.
         ! numNarrowLines <0 means unknown, need to check, =0 means none, >0 number of narrow lines
         numNarrowLines = -1
         narrowWidth = ALFAV !ALFAL0
      else
         !--- If not very different, just use the DVnarrow
         numNarrowLines = 0
         DV = DVnarrow
      endif
! Skip narrow lines code   KCP !
      numNarrowLines = 0


      !--- Allocate buffer(s) for OD spectrum/spectra
      !
      ! R1PRNT() PRINTS THE FIRST NPTS=MPTS VALUES STARTING AT JLO
      ! AND THE LAST NPTS=MPTS VALUES ENDING AT NLIM OF THE R1 ARRAY
      !
      sizeOfBuffer = ceiling( (V2-V1)/DV  + 1. )
! MJI - Revised sizeOfBuffer check by adding v2eps to allow check to pass and add extra point beyond numerical precision
      if (V1+real(sizeOfBuffer-1)*DV <=V2+v2eps) sizeOfBuffer=sizeOfBuffer+1  ! "+1" to ensure the last point is beyond V2

      allocate( r1BuffNorm( sizeOfBuffer ), STAT=iStat )
      !if (NLTE==.TRUE.) allocate( rr1BuffNorm( sizeOfBuffer ), STAT=iStat )
      allocate( rr1BuffNorm( sizeOfBuffer ), STAT=iStat ) !This can be conditional on IHRAC later when all sub's are placed in Fortran module. Cannot be done now because optional argument need explicit interface and allocatable argument needs allocated real argument.

      !--- Calculate OD for this layer ---
      CALL OPDPTH( r1BuffNorm, rr1BuffNorm, numNormPoints, &
                   LAYER,ALTAV,PAVE,TAVE,lnMolNames(1:nLnMol),WK,WXM,WTOT,WKI,ISOTPL_FLAG, &
                   DV,V1,V2,NLTE_flag,speciesBroad,SAMPLE,ALFAV,DPTFAC,DPTMIN,IPATHL,JRAD,ILNFLG,LINFIL, &
                   ILBLF4,LNFIL4,JCNVF4,&
                   XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL,IAERSL,icflg,&
                   odCtrl%lineOff, XsectOff, contnmOff,IXSBIN,XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                   NOPR,NFHDRF,NPHDRF,NPHDRL,NLNGTH, &
                   numNarrowLines, narrowWidth)


      !--- For very narrow lines, do a second round convolution on finer spectral grid.
      if ( numNarrowLines>0 ) then

         !--- Do a second round convolution using selected narrow lines
         LINFIL = narrowLineFile
         numNarrowLines = 0

         if (NLTE_flag) then
            rewind( nlteStatPopFileLun )
         endif

         XsectOff  = .TRUE.
         contnmOff = .TRUE.
         JRAD      = odCtrl%JRAD
         !JCNVF4    = 0
         !WK(:)     = !may be zeroed when lineOff=.TRUE.
         DVnorm    = DV
         DV        = DVnarrow !Do narrow line convolution using finer DV

         sizeOfBuffer = ceiling( (V2-V1)/DV  + 1. )
         if (V1+real(sizeOfBuffer-1)*DV <=V2) sizeOfBuffer=sizeOfBuffer+1  ! "+1" to ensure the last point is beyond V2

         allocate( r1BuffFine( sizeOfBuffer ), STAT=iStat )
         allocate( rr1BuffFine( sizeOfBuffer ), STAT=iStat ) !This can be conditional on IHRAC later when all sub's are placed in Fortran module. Cannot be done now because optional argument need explicit interface and allocatable argument needs allocated real argument.

         CALL OPDPTH( r1BuffFine, rr1BuffFine, numFinePoints, &
                      LAYER,ALTAV,PAVE,TAVE,lnMolNames(1:nLnMol),WK,WXM,WTOT,WKI,ISOTPL_FLAG, &
                      DV,V1,V2,NLTE_flag,speciesBroad,SAMPLE,ALFAV,DPTFAC,DPTMIN,IPATHL,JRAD,ILNFLG,LINFIL, &
                      ILBLF4,LNFIL4,JCNVF4,&
                      XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL,IAERSL,icflg,&
                      odCtrl%lineOff, XsectOff, contnmOff,IXSBIN,XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                      NOPR,NFHDRF,NPHDRF,NPHDRL,NLNGTH, &
                      numNarrowLines, narrowWidth)


         allocate( tempBuff( size(r1BuffFine) ))

         !--- Combine fine line OD and normal line OD
         call interpSpectrum( r1BuffNorm, DVnorm, V1, numNormPoints, &
                              tempBuff,   DV,     V1, numFinePoints )
         do iv=1,numFinePoints
            r1BuffFine(iv) = r1BuffFine(iv) + tempBuff(iv)
         enddo

         !--- Combine file line NLTE factor and normal line NLTE factor
         if (NLTE_flag) then
            call interpSpectrum( rr1BuffNorm, DVnorm, V1, numNormPoints, &
                                 tempBuff,    DV,     V1, numFinePoints )
            do iv=1,numFinePoints
               rr1BuffFine(iv) = rr1BuffFine(iv) + tempBuff(iv)
            enddo
         endif

         !---
         deallocate( tempBuff )
         deallocate( r1BuffNorm )
         deallocate( rr1BuffNorm )

         !--- Delete the scratch file
         !CLOSE( narrowLineFile,  STATUS='delete')
      else

         call move_alloc( r1BuffNorm, r1BuffFine )
         call move_alloc( rr1BuffNorm, rr1BuffFine )
         numFinePoints = numNormPoints

      endif !if ( numNarrowLines>0 )

      !--- Delete the scratch file
      CLOSE( narrowLineFile,  STATUS='delete')


      !--- If user request uniform DV out, interpolate the OD spectra.
      outputDV = DV
      if (present(DVout)) then
         if ( DVout >0. )  outputDV = DVout
      endif

      if ( abs(outputDV-DV) < epsilon(DV) ) then !use exact DV

         if (size(r1BuffFine)==numFinePoints) then
            call move_alloc( r1BuffFine, monoOD )
         else
            allocate( monoOD(numFinePoints) )
            monoOD(1:numFinePoints) = r1BuffFine(1:numFinePoints)
         endif

         if (NLTE_flag) then
            if (size(rr1BuffFine)==numFinePoints) then
               call move_alloc( rr1BuffFine, nlteEmisFac )
            else
               allocate(nlteEmisFac(numFinePoints))
               nlteEmisFac(1:numFinePoints) = rr1BuffFine(1:numFinePoints)
            endif
         endif

         nSamp = numFinePoints

      else !Interpolate from DV to outputDV

         if (outputDV > DV) Print*, '--- '//routineName//'(): Optical depth is interpolated to coarser grid.'

         lengthOfSpect = ceiling( (V2-V1)/outputDV  + 1. )
         allocate( monoOD(lengthOfSpect) )
         call interpSpectrum( r1BuffFine, DV, V1, numFinePoints, &
                              monoOD, outputDV, V1, lengthOfSpect )

         if (NLTE_flag) then
            allocate( nlteEmisFac(lengthOfSpect) )
            call interpSpectrum( rr1BuffFine, DV, V1, numFinePoints, &
                                 nlteEmisFac, outputDV, V1, lengthOfSpect )
         endif

         nSamp = lengthOfSpect
      endif

      !--- Output the actual DV used
      if (present(DVout)) DVout = outputDV


      !--- Deallocate the buffers
      if (allocated(r1BuffNorm))  deallocate(r1BuffNorm)
      if (allocated(r1BuffFine))  deallocate(r1BuffFine)
      if (allocated(rr1BuffNorm)) deallocate(rr1BuffNorm)
      if (allocated(rr1BuffFine)) deallocate(rr1BuffFine)
      if (allocated(tempBuff))    deallocate(tempBuff)

   END SUBROUTINE


   !-----------------------------------------------------------------
   ! Check each selected molecule, to see if it is a isotopologue, a line molecule, a
   ! cross-section (xs) molecule or a molecule with both line and xs. If the molecule has both line data and
   ! xs data in the interval (V1~V2), check further for each sub interval
   ! where xs data are available to see which one is preferred, line or xs?
   ! The selected molecules are grouped into three sub groups, isotoplogue group,
   ! line group and xs group. Those that has both line and xs data in
   ! V1~V2 are assigned to both line group and xs group.
   !-----------------------------------------------------------------
   subroutine seperateLineXsISO( molNames, nMol, V1, V2,&
                                 lnMolNames, xsMolNames, isoMolNames, &
                                 nLnMol, nXsMol, nISO )
   !-----------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, isISOName, molNum, xsMolNum
      USE Module_XSect       ,ONLY: xsTbl

      character(*)              ,intent(in)  :: molNames(:)
      integer                   ,intent(in)  :: nMol
      real(r8)                  ,intent(in)  :: V1,V2
      character(*) ,allocatable ,intent(out) :: lnMolNames(:)
      character(*) ,allocatable ,intent(out) :: xsMolNames(:)
      character(*) ,allocatable ,intent(out) :: isoMolNames(:)
      integer                   ,intent(out) :: nLnMol, nXsMol, nISO

      character(*),parameter :: routineName = 'seperateLineXsISO'
      integer     ,parameter :: MaxNumXsSpect = 15
      ! increase number of allowed xs molecules to 15 (LBLRTM allows 50)
      ! integer     ,parameter :: MaxNumXsSpect = 6

      integer  :: im,is, kln,kxs,kis
      integer  :: mnum,xnum
      logical  :: useLine, useXS
      real(r8) :: v1x,v2x
      integer  ,allocatable :: molMark(:)



      !--- Mark the type of molecule
      ! mark value =1 use line data
      !            =2 use xs data
      !            =3 isotopologue
      !            =4 use both line and xs in V1~V2 region
      !
      allocate( molMark(nMol) )
      nLnMol = 0
      nXsMol = 0
      nISO = 0
      do im = 1,nMol

         if ( isISOName( molNames(im) ) ) then ! molecule is a isotopologue

            molMark(im) = 3
            nISO = nISO+1

         else !line molecules or/and xs molecules

            mnum = molNum( molNames(im) )
            xnum = xsMolNum( molNames(im) )
            useLine = .FALSE.
            useXS   = .FALSE.

            if (mnum>0 .and. xnum>0) then !both line data and xs exist

               !--- check if there are cross-sections exist in V1~V2 interval and the preference
               ! FSCDXS value for preference
               ! 0: only XS exist for this molecule in this spectral range
               ! 1: both line parameters and XS exist, line parameters are preferable
               ! 2: both line parameters and XS exist, XS are preferable
               do is = 1,xsTbl%numSpect(xnum)

                  v1x = xsTbl%V1FX(is,xnum)
                  v2x = xsTbl%V2FX(is,xnum)

                  if ( (v2x.GT.V1.AND.v1x.LT.V2) ) then !The molecule has xs in V1~V2
                     if ( xsTbl%lineOrXs(is,xnum)==1 ) then !preference for line over xs
                        useLine = .TRUE.
                     else
                        useXS = .TRUE.
                     endif
                  else ! xs not in V1~V2
                     useLine = .TRUE.
                  endif
               enddo

            elseif (mnum>0 ) then !is a line molecule
               useLine = .TRUE.
            elseif (xnum>0) then !is a xs molecule
               useXS = .TRUE.
            else
               STOP '--- '//routineName//'(): Unknown molecular name.'
            endif

            if (useLine .and. useXS) then !use line
               molMark(im) = 4
               nLnMol = nLnMol+1
               nXsMol = nXsMol+1
            elseif (useLine) then
               molMark(im) = 1
               nLnMol = nLnMol+1
            elseif (useXS) then
               molMark(im) = 2
               nXsMol = nXsMol+1
            endif

         endif !( isISOName( molNames(im) ) ) then
      enddo !im = 1,nMol


      !--- Assign molecules to sub groups. Molecules with both line and
      ! xs data available are assigned to both line and xs sub groups.
      !
      allocate(lnMolNames(nLnMol))
      allocate(xsMolNames(nXsMol))
      allocate(isoMolNames(nISO))
      kln=0
      kxs=0
      kis=0
      do im = 1,nMol
         select case ( molMark(im) )
         case(1)
            kln = kln+1
            lnMolNames( kln ) = molNames(im)
         case(2)
            kxs = kxs+1
            xsMolNames( kxs ) = molNames(im)
         case(3)
            kis = kis+1
            isoMolNames( kis ) = molNames(im)
         case(4)
            kln = kln+1
            kxs = kxs+1
            lnMolNames( kln ) = molNames(im)
            xsMolNames( kxs ) = molNames(im)
         end select
      enddo

   END SUBROUTINE




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE OPDPTH( r1Buffer, rr1Buffer, numTotPoints, &
                         LAYER,ALTAV,PAVE,TAVE,lnMolNames,WK,WXM,WTOT,WKI,ISOTPL_FLAG,&
                         DV,V1,V2,NLTE,speciesBroad,SAMPLE,ALFAV,DPTFAC,DPTMIN,IPATHL,JRAD,ILNFLG,LINFIL, &
                         ILBLF4,LNFIL4,JCNVF4,&
                         XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL,IAERSL,icflg,&
                         lineOff,XsectOff,contnmOff,IXSBIN,XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                         NOPR,NFHDRF,NPHDRF,NPHDRL,NLNGTH, &
                         numNarrowLines, narrowWidth)
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, RADCN2, &
                                    N_ABSRB, MAXSTATE, Max_ISO, MXMOL
      USE Module_LineF4      ,ONLY: LINF4Q
      USE Module_Continuum   ,ONLY: CONTNM
      USE Module_Config      ,ONLY: IPR

      IMPLICIT NONE !REAL*8           (V)

      real         ,intent(out)   :: r1Buffer(:)
      real         ,intent(out)   :: rr1Buffer(:)
      integer      ,intent(out)   :: numTotPoints
      !
      integer      ,intent(in)    :: LAYER
      real         ,intent(in)    :: ALTAV !used in non-LTE
      real         ,intent(in)    :: PAVE
      real         ,intent(in)    :: TAVE
      character(*) ,intent(in)    :: lnMolNames(:)
      real         ,intent(inout) :: WK(:)   !(60), may be set to zero when lineOff=.TRUE.
      real         ,intent(in)    :: WXM(:)
      real         ,intent(in)    :: WTOT
      real         ,intent(inout) :: WKI(:,:) !like WK, WKI may be changed by lineOff
      integer      ,intent(in)    :: ISOTPL_FLAG(:,:)
      !
      real         ,intent(inout) :: DV    !may be modified by READ from from REJ file
      real(r8)     ,intent(in)    :: V1, V2
      logical      ,intent(in)    :: NLTE
      logical      ,intent(in)    :: speciesBroad
      real         ,intent(in)    :: SAMPLE, ALFAV
      real         ,intent(in)    :: DPTFAC
      real         ,intent(inout) :: DPTMIN !will be modified and restored
      integer      ,intent(in)    :: IPATHL
      integer      ,intent(inout) :: JRAD
      integer      ,intent(in)    :: ILNFLG
      integer      ,intent(in)    :: LINFIL
      !
      integer      ,intent(inout) :: ILBLF4 !may be changed by lineOff
      integer      ,intent(in)    :: LNFIL4
      integer      ,intent(inout) :: JCNVF4
      !
      real         ,intent(in)    :: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
      integer      ,intent(in)    :: IAERSL !used in the continuum module
      integer      ,intent(in)    :: icflg  !used in the continuum module for derivative
      !
      logical      ,intent(in)    :: lineOff,XsectOff,contnmOff
      integer      ,intent(in)    :: IXSBIN
      character*10 ,intent(in)    :: XSFILE(:,:,:) !(6,5,MX_XS)
      integer      ,intent(in)    :: IXFORM(:,:)   !(5,MX_XS)
      integer      ,intent(in)    :: NSPECR(:)     !(MX_XS)
      integer      ,intent(in)    :: NTEMPF(:,:)   !(5,MX_XS)
      real         ,intent(in)    :: XDOPLR(:,:)   !(5,MX_XS)
      integer      ,intent(in)    :: LnOrXs(:,:)   !(5,MX_XS)
      integer      ,intent(in)    :: NUMXS
      !
      integer      ,intent(in)    :: NOPR
      integer      ,intent(in)    :: NFHDRF,NPHDRF,NPHDRL,NLNGTH
      integer      ,intent(inout) :: numNarrowLines
      real         ,intent(in)    :: narrowWidth


      !--- Local variables
      !
      integer   :: NPTABS                   ! ABSORB
      real      :: ABSRB(N_ABSRB)           ! ABSORB
      real      :: DVABS                    ! ABSORB
      real(r8)  :: V1ABS, V2ABS             ! ABSORB
      real      :: BOUND4
      real      :: DVR4
      real(r8)  :: V1R4, V2R4
      integer   :: IPFLAG
      integer   :: rejLAYER

      integer   :: NEGEPP_FLAG
      real      :: ratState(MAXSTATE*Max_ISO,MXMOL)
      integer   :: numState(MXMOL)

      integer ,PARAMETER :: I_10=10

      INTEGER  :: I,       IEMST,  IPTS4
      INTEGER  :: M
      REAL     :: ALFAV4,  DPTMST
      REAL(r8) :: V1L4,    V2L4
      REAL     :: XKT
!yma       REAL     :: TIME0,   TIME1,   TIMEO



      IPFLAG = 0

      DPTMST = DPTMIN
      IF (IPATHL.EQ.2) DPTMIN = 2.*DPTMST

      if (.not.contnmOff .OR. NLTE) then !if continuum or NLTE
         IF (PAVE.LE.0.5) IPFLAG = 1
      endif

      !      PRINT LAYER INFORMATION
      IF (NOPR.EQ.0) THEN
         !IF (IMRG.LE.10) WRITE (IPR,900)
         WRITE (IPR,905) LAYER
         !IF (ILAS.GT.0) WRITE (IPR,910) VLAS,V1,V2
         !WRITE (IPR,915) XID,(YID(M),M=1,2),TIME0  !141129 ytma
      ENDIF

      if (NLTE) then
         call nonLTE_RatState( RATSTATE, NUMSTATE, ALTAV, TAVE )
      endif


      !     JRAD= -1  NO RADIATION TERM IN ABSORPTION COEFFICIENTS
      !     JRAD=  0  RADIATION TERM PUT IN BY PANEL
      !     JRAD=  1  RADIATION TERM INCLUDED IN LINE STRENGTHS
      XKT = TAVE/RADCN2
      IF (((V1/XKT).LT. 5.).AND.(JRAD.NE.-1)) JRAD = 0

      DVABS = 0.
      IF (.not.contnmOff) THEN
         DVABS = 1.
         V1ABS = INT(V1)
         IF (V1.LT.0.) V1ABS = V1ABS-1.
         V1ABS = V1ABS-3.*DVABS
         V2ABS = INT(V2+3.*DVABS+0.5)
         NPTABS = (V2ABS-V1ABS)/DVABS+1.5
         !yma IF (PAVE.LE.0.5) IPFLAG = 1
         DO 10 I = 1, n_absrb
            ABSRB(I) = 0.
   10    CONTINUE
         CALL CONTNM( JRAD, &
                      V1ABS,V2ABS,DVABS,NPTABS,ABSRB, &
                      XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL, &
                      V1,V2,LAYER,PAVE,TAVE,Wtot,WK,IAERSL,ICFLG )
      ENDIF


      DVR4 = 0.
      !yma IF (ILBLF4.GE.1) THEN
      IF (ILBLF4.GE.1 .or. numNarrowLines<0) THEN
         !yma ALFAV = SAMPLE*DV
         !yma ALFAV4 = 64.*ALFAV
         ALFAV4 = 64.*SAMPLE*DV
         !     Read in DVR4 from REJ file
         IF (ILNFLG.EQ.2) THEN
            READ(16) rejLAYER, DVR4
         ELSE !yma ILNFLG==-1,0,1
         !     Compute DVR4
            DVR4 = ALFAV4/SAMPLE
            IF (ILNFLG.EQ.1) WRITE(16) LAYER, DVR4
         ENDIF
         BOUND4 = 25.
         IF (ILBLF4.EQ.2.AND.IPFLAG.EQ.1) BOUND4 = 5.
         IPTS4 = BOUND4/DVR4

         IF (NOPR.EQ.0) WRITE (IPR,920) IPTS4,DVR4,BOUND4

         REWIND LINFIL
         REWIND LNFIL4
         V1R4 = V1-2.*DVR4
         V2R4 = V2+2.*DVR4
         V1L4 = V1R4-BOUND4-DVR4
         V2L4 = V2R4+BOUND4+2*DVR4
         CALL LINF4Q( V1L4,V2L4, &
                      LINFIL,LNFIL4,NEGEPP_FLAG,DVR4,V1,V2,&
                      PAVE,TAVE,lnMolNames,WK,WTOT,WKI,ISOTPL_FLAG,&
                      DPTFAC,DPTMIN,NLTE,speciesBroad,NUMSTATE,RATSTATE,&
                      NPHDRL,NLNGTH,NOPR, &
                      numNarrowLines, narrowWidth)
      ENDIF


      !    Write out DV to REJ1 file
      IF (ILNFLG.EQ.1) WRITE(15) LAYER, DV
      !    Read in DV from  REJ1 file
      IF (ILNFLG.EQ.2) READ(15) rejLAYER,DV

      CALL HIRACQ( r1Buffer, rr1Buffer, numTotPoints,&
                   LAYER,PAVE,TAVE,lnMolNames,WK,WXM,WTOT,WKI,ISOTPL_FLAG,&
                   NUMSTATE,RATSTATE,&
                   DV,V1,V2,NLTE,speciesBroad,SAMPLE,ALFAV,DPTFAC,DPTMIN,JRAD,ILNFLG,LINFIL,&
                   ILBLF4,LNFIL4,NEGEPP_FLAG,DVR4,V1R4,V2R4,BOUND4,JCNVF4,&
                   ABSRB,DVABS,V1ABS,V2ABS,&
                   lineOff,XsectOff,contnmOff,IXSBIN,XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                   NOPR,NFHDRF,NPHDRF,NLNGTH, &
                   numNarrowLines, narrowWidth)

      DPTMIN = DPTMST
      !CALL CPUTIM (TIME1)
      !TIMEO = TIME1-TIME0
      !WRITE (IPR,925) TIME1,TIMEO

      RETURN

  900 FORMAT ('1')
  905 FORMAT (/'0 LAYER = ',I8)
  910 FORMAT ('0 VLAS  ',F20.8,8X,'V1 RESET ',F12.5,8X,'V2 RESET ',     &
     &        F12.5)
  915 FORMAT ('0',10A8,2X,2(1X,A8,1X),/,'0 TIME ENTERING OPDPTH ',      &
     &        F15.3)
  920 FORMAT ('0  IPTS4 FOR LINF4 = ',I5,3X,' DV FOR LINF4 = ',F10.5,   &
     &        5X,'BOUND FOR LINF4 =',F10.4)
  925 FORMAT ('0 TIME LEAVING OPDPTH ',F15.3,'  TOTAL FOR LAYER ',      &
     &        F15.3)

      END  SUBROUTINE


!-----------------------------------------------------------------------
!*
!*    CALCULATES MONOCHROMATIC ABSORPTION COEFFICIENT FOR SINGLE LAYER
!*
!*
!*            USES APPROXIMATE VOIGT ALGORITHM
!*
!*
!*              VAN VLECK WEISSKOPF LINE SHAPE
!*
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
!                                     M.W.SHEPHARD
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     131 Hartwell Ave,  Lexington,  MA   02421
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
!-----------------------------------------------------------------------
      SUBROUTINE HIRACQ( r1Buffer, rr1Buffer, numTotPoints,&
                         LAYER,PAVE,TAVE,lnMolNames,WK,WXM,WTOT,WKI,ISOTPL_FLAG,&
                         NUMSTATE,RATSTATE,&
                         DV,V1,V2,NLTE,speciesBroad,SAMPLE,ALFAV,DPTFAC,DPTMIN,JRAD,ILNFLG,LINFIL,&
                         ILBLF4,LNFIL4,NEGEPP_FLAG,DVR4,V1R4,V2R4,BOUND4,JCNVF4,&
                         ABSRB,DVABS,V1ABS,V2ABS,&
                         lineOff,XsectOff,contnmOff,IXSBIN,XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                         NOPR,NFHDRF,NPHDRF,NLNGTH, &
                         numNarrowLines, narrowWidth)
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8,RADCN2, &
                                         MXMOL,MXISOTPL,NFPTS,NFMX, &
                                         molNum
      USE Module_LineData         ,ONLY: MXBRDMOL,NLINEREC
      USE Module_ODLAY_CommonSub  ,ONLY: MOLEC,RADFN
      USE Module_Spectrum         ,ONLY: XINT
      USE Module_LineF4           ,ONLY: LBLF4Q
      USE Module_XSect            ,ONLY: XSECTM
      USE Module_Config           ,ONLY: IPR

      IMPLICIT NONE !REAL*8           (V)

      real         ,intent(out)   :: r1Buffer(:) !may change to r1Buffer(:) in the future
      real         ,intent(out)   :: rr1Buffer(:) !may change to rr1Buffer(:) in the future
      integer      ,intent(out)   :: numTotPoints
      !
      integer      ,intent(in)    :: LAYER
      real         ,intent(in)    :: PAVE
      real         ,intent(in)    :: TAVE                          ! FILHDR
      character(*) ,intent(in)    :: lnMolNames(:)
      real         ,intent(inout) :: WK(:)            !(60), may be set to zero when lineOff=.TRUE.
      real         ,intent(in)    :: WXM(:)           !(MX_XS)
      real         ,intent(in)    :: WTOT
      real         ,intent(inout) :: WKI(:,:)         !(MXMOL,MXISOTPL) may be set to zero when lineOff=.TRUE.
      integer      ,intent(in)    :: ISOTPL_FLAG(:,:) !(MXMOL,MXISOTPL)
      !
      integer      ,intent(in)    :: NUMSTATE(:)      !(MXMOL)
      real         ,intent(in)    :: RATSTATE(:,:)    !(MAXSTATE*Max_ISO,MXMOL)
      !
      real         ,intent(in)    :: DV                            ! FILHDR
      real(r8)     ,intent(in)    :: V1, V2                        ! FILHDR
      logical      ,intent(in)    :: NLTE
      logical      ,intent(in)    :: speciesBroad
      real         ,intent(in)    :: SAMPLE, ALFAV                 ! MANE
      real         ,intent(in)    :: DPTFAC, DPTMIN                ! MANE
      integer      ,intent(in)    :: JRAD
      integer      ,intent(in)    :: ILNFLG
      integer      ,intent(in)    :: LINFIL
      !
      integer      ,intent(inout) :: ILBLF4 !may be changed by lineOff
      integer      ,intent(in)    :: LNFIL4
      integer      ,intent(in)    :: NEGEPP_FLAG
      real         ,intent(in)    :: DVR4                          ! LBLF
      real(r8)     ,intent(inout) :: V1R4
      real(r8)     ,intent(inout) :: V2R4                          ! LBLF
      real         ,intent(in)    :: BOUND4
      integer      ,intent(inout) :: JCNVF4
      !
      real         ,intent(in)    :: ABSRB(:) !(N_ABSRB)           ! ABSORB
      real         ,intent(in)    :: DVABS                         ! ABSORB
      real(r8)     ,intent(in)    :: V1ABS, V2ABS                  ! ABSORB
      !
      logical      ,intent(in)    :: lineOff,XsectOff,contnmOff
      integer      ,intent(in)    :: IXSBIN
      character*10 ,intent(in)    :: XSFILE(:,:,:) !(6,5,MX_XS)
      integer      ,intent(in)    :: IXFORM(:,:)   !(5,MX_XS)
      integer      ,intent(in)    :: NSPECR(:)     !(MX_XS)
      integer      ,intent(in)    :: NTEMPF(:,:)   !(5,MX_XS)
      real         ,intent(in)    :: XDOPLR(:,:)   !(5,MX_XS)
      integer      ,intent(in)    :: LnOrXs(:,:)   !(5,MX_XS)
      integer      ,intent(in)    :: NUMXS
      !
      integer      ,intent(in)    :: NOPR                          ! IFIL
      integer      ,intent(in)    :: NFHDRF, NPHDRF, NLNGTH        ! IFIL
      !
      integer      ,intent(in)    :: numNarrowLines
      real         ,intent(in)    :: narrowWidth


      !--- Local variables
      !
      integer           :: lbR1 ,ubR1
      integer           :: lbR2 ,ubR2
      integer           :: lbR3 ,ubR3
      integer           :: lbR4 ,ubR4
      real ,allocatable :: R1(:) ,RR1(:)
      real ,allocatable :: R2(:) ,RR2(:)
      real ,allocatable :: R3(:) ,RR3(:)
      real ,allocatable :: R4(:) ,RR4(:)

      real(r8)    :: VNU(   NLINEREC)                      ! %BLANK
      real        :: SP(    NLINEREC)                      ! %BLANK
      real        :: ALFA0( NLINEREC)
      real        :: EPP(   NLINEREC)
      integer*4   :: MOL(   NLINEREC)
      real        :: HWHMS( NLINEREC)
      real        :: TMPALF(NLINEREC)
      real        :: PSHIFT(NLINEREC)
      integer     :: IFLG(  NLINEREC)
      real        :: SPPSP( NLINEREC)                    ! %BLANK
      real        :: RECALF(NLINEREC)                    ! %BLANK
      real        :: ZETAI( NLINEREC)                    ! %BLANK
      integer     :: IZETA( NLINEREC)                    ! %BLANK
      integer*4   :: BRD_MOL_FLG( MXBRDMOL,NLINEREC)
      real        :: BRD_MOL_HW(  MXBRDMOL,NLINEREC)
      real        :: BRD_MOL_TMP( MXBRDMOL,NLINEREC)
      real        :: BRD_MOL_SHFT(MXBRDMOL,NLINEREC)
      integer     :: IOUT(NLINEREC)

      integer     :: N1MAX, N2MAX, N3MAX
      integer     :: NX1,   NX2,   NX3                ! CMSHAP
      real        :: DXF1,  DXF2,  DXF3
      real        :: HWF1,  HWF2,  HWF3               ! CMSHAP

      integer     :: MAX1,  MAX2,  MAX3
      integer     :: N1R1,  N1R2,  N1R3               ! SUB1
      integer     :: N2R1,  N2R2,  N2R3
      integer     :: NLIM1, NLIM2, NLIM3              ! SUB1
      integer     :: NLO, NHI                         ! SUB1
      real        :: DVR2, DVR3                       ! SUB1
      real        :: DVP                              ! XPANEL
      integer     :: NSHIFT                           ! XPANEL
      real(r8)    :: VBOT, VTOP                       ! XSUB
      real(r8)    :: VFT                              ! XSUB
      integer     :: LIMIN                            ! XSUB
      integer     :: ILO, IHI                         ! XSUB
      integer     :: IEOF, IPANEL, ISTOP, IDATA
      integer     :: NPTR4                            !LBLF

      integer     :: LINCNT, NCHNG, NLIN, NMINUS, NPLUS   ! LNC1
      real        :: SCOR(MXMOL,MXISOTPL)  !(42,10)       ! LNC1  !for dummy only
      real        :: RHOSLF(MXMOL)                        ! LNC1  !for dummy only
      real        :: ALFD1(MXMOL,MXISOTPL) !(42,10)       ! LNC1  !for dummy argument only

      real        :: ALFMAX
      real        :: DPTFC                      ! LNC1
      real        :: DPTMN

      integer     :: nALFMAX

      real        :: SUMALF, SUMZET             ! LNC1
      integer     :: ILIN4, ILIN4T              ! R4SUB
      integer     :: NPTS                       ! XPANEL     !no use

!yma       integer     :: L4NLN,   L4NLS    ! L4TIMG
!yma       real        :: L4TIM,   L4TMR    ! L4TIMG
!yma       real        :: L4TMS,   LOTHER   ! L4TIMG
!yma       real        :: TF4,     TF4CNV   ! XTIME
!yma       real        :: TF4PNL,  TF4RDF   ! XTIME
!yma       real        :: TIMCNV,  TIME     ! XTIME
!yma       real        :: TIMPNL,  TIMRDF   ! XTIME
!yma       real        :: TXS,     TXSCNV   ! XTIME
!yma       real        :: TXSPNL,  TXSRDF   ! XTIME

      !--- Saved local variables
      integer  ,SAVE :: IFN=0                  ! FNSHQ
      real     ,SAVE :: F1(NFMX),  F2(NFMX)    ! FNSHQ
      real     ,SAVE :: F3(NFMX),  FG(NFMX)    ! FNSHQ
      real     ,SAVE :: XVER(NFMX)             ! FNSHQ


      ! NOTE that DXFF1 = (HWFF1/(NFPTS-1))
      ! and       DXFF2 = (HWFF2/(NFPTS-1))
      ! and       DXFF3 = (HWFF3/(NFPTS-1))
      real    ,PARAMETER :: HWFF1=4.    ,HWFF2=16.   ,HWFF3=64.
      real    ,PARAMETER :: DXFF1=0.002 ,DXFF2=0.008 ,DXFF3=0.032
      integer ,PARAMETER :: NXF1=NFPTS  ,NXF2=NFPTS  ,NXF3=NFPTS
      integer ,PARAMETER :: NF1MAX=NFMX ,NF2MAX=NFMX ,NF3MAX=NFMX

      integer ,PARAMETER :: I_10=10

      integer :: MEFDP(64)=0

      integer      :: imol
      INTEGER      :: I,         IENTER,    IFPAN
      INTEGER      :: IFST,      IR4
      INTEGER      :: I_TIME_LAY,           LTIME,     M
      INTEGER      :: NBOUND,    NLNCR
      REAL         :: AVALF,     AVZETA,    BOUND,     BOUNF3
      REAL         :: DVSAV,     DV_LBL
      REAL(r8)     :: V1R4ST,    V2R4ST,    VF1,       VF2



!     SET INPUT FLAG FOR USE BY X-SECTIONS

      IFST = -99
      IR4 = 0
      IENTER = 0

!     SET COMMON BLOCK CMSHAP

      HWF1 = HWFF1
      DXF1 = DXFF1
      NX1 = NXF1
      N1MAX = NF1MAX
      HWF2 = HWFF2
      DXF2 = DXFF2
      NX2 = NXF2
      N2MAX = NF2MAX
      HWF3 = HWFF3
      DXF3 = DXFF3
      NX3 = NXF3
      N3MAX = NF3MAX

      DPTMN = DPTMIN
      IF (JRAD.NE.1) DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)
      DPTFC = DPTFAC
      ILIN4 = 0
      ILIN4T = 0
      LIMIN = 250
      NSHIFT = 32

!v128--->
!     SAMPLE IS AVERAGE ALPHA / DV
! "*0.04/ALFAL0" is a temporary solution to keep ALFMAX to be close to the
! original default value. It may enlarge the ALFMAX in regions dominated
! Doppler broadening where ALFAL0 is supposed to have no effect if user
! sets ALFAL0 to a value less than 0.04.
!    previous expression:  NBOUND = 4.*(2.*HWF3)*SAMPLE+0.01
!
      ALFMAX = 4.*ALFAV !4*SAMPLE*DV * 0.04/ALFAL0
      nALFMAX = ALFMAX/DV
      NBOUND = (2.*HWF3)*nALFMAX+0.01

      NLIM1 = 2401
      NLIM2 = (NLIM1/4)+1
      NLIM3 = (NLIM2/4)+1

      IF (IFN.EQ.0) THEN
         !CALL CPUTIM(TPAT0)
         CALL SHAPEL( F1,F2,F3, &
                      HWF1,HWF2,HWF3,DXF1,DXF2,DXF3,NX1,NX2,NX3,N1MAX,N2MAX,N3MAX )
         CALL SHAPEG( FG,DXF1,NX1,N1MAX )
         CALL VERFN( XVER,DXF1,N1MAX )
         IFN = IFN+1
         !CALL CPUTIM(TPAT1)
         !TSHAPE = TSHAPE+TPAT1-TPAT0
      ENDIF

      !CALL CPUTIM(TPAT0)
      CALL MOLEC( 1,SCOR,RHOSLF,ALFD1, &
                  PAVE,TAVE,WK,WTOT,lnMolNames )

      !CALL CPUTIM(TPAT1)
      !TMOLEC = TMOLEC+TPAT1-TPAT0
      REWIND LINFIL
      IEOF = 0
      ILO = 0
      IHI = -999
      NMINUS = 0
      NPLUS = 0

!     NOTE (DXF3/DXF1) IS 16 AND (DXF3/DXF2) IS 4

      DVP = DV
      DVR2 = (DXF2/DXF1)*DV
      DVR3 = (DXF3/DXF1)*DV
      MAX1 = NSHIFT+NLIM1+NSHIFT+NBOUND/2
      MAX2 = MAX1/4
      MAX3 = MAX1/16
      MAX1 = MAX1 +  4*nALFMAX
      MAX2 = MAX2 + 16*nALFMAX/4 + 1  !+1 to account for rounding error
      MAX3 = MAX3 + 64*nALFMAX/16 + 1 !+1 to account for rounding error

      !CALL CPUTIM(TPAT0)
      BOUND =  REAL(NBOUND)*DV/2.
      BOUNF3 = BOUND/2.
      NLO = NSHIFT+1
      NHI = NLIM1+NSHIFT-1

      lbR1 = -ceiling(  real( 4*(HWF3+4)  )*ALFAV/DV +1.)
      lbR2 = -ceiling( (real( 4*(HWF3+16) )*ALFAV/DV)/4. +1.)
      lbR3 = -ceiling( (real( 4*(HWF3+64) )*ALFAV/DV)/16. +1.)
      ubR1 = ceiling(  real(2*NSHIFT+NLIM1) + real( 4*(HWF3+4 ) )*ALFAV/DV +1.)
      ubR2 = ceiling( (real(2*NSHIFT+NLIM1) + real( 4*(HWF3+16) )*ALFAV/DV)/4. +1.)
      ubR3 = ceiling( (real(2*NSHIFT+NLIM1) + real( 4*(HWF3+64) )*ALFAV/DV)/16. +1.)
      lbR4 = 1
      ubR4 = 2502

      allocate( R1(lbR1:ubR1) )
      allocate( R2(lbR2:ubR2) )
      allocate( R3(lbR3:ubR3) )
      allocate( R4(lbR4:ubR4) )

      DO 10 I = 1, MAX1 !yma,151201
         R1(I) = 0.
   10 END DO
      DO 20 I = 1, MAX2
         R2(I) = 0.
   20 END DO
      DO 30 I = 1, MAX3
         R3(I) = 0.
   30 END DO
      IF (ILBLF4.EQ.0) THEN
         DO 40 I = 1, ubR4!2502
            R4(I) = 0.
   40    CONTINUE
      ENDIF
      if ( NLTE ) then !yma,151201
         DO I = 1, MAX1
            RR1(I) = 0.
         END DO
         DO I = 1, MAX2
            RR2(I) = 0.
         END DO
         DO I = 1, MAX3
            RR3(I) = 0.
         END DO
         IF (ILBLF4.EQ.0) THEN
            DO I = 1, ubR4!2502
               RR4(I) = 0.
            ENDDO
         ENDIF
      endif

      !CALL CPUTIM(TPAT1)
      !TLOOPS = TLOOPS + TPAT1-TPAT0

      !CALL CPUTIM(TPAT1)
      !TODFIL = TODFIL + TPAT1-TPAT0

      IF (lineOff) THEN
!          DO 50 M = 1, NMOL
!             WK(M) = 0.
!    50    CONTINUE
         WK(:) = 0.
         WKI(:,:) = 0.
         ILBLF4 = 0
      ENDIF


      VFT = V1- REAL(NSHIFT)*DV
      VBOT = V1-BOUND
      VTOP = V2+BOUND

      LINCNT = 0
      NLIN = 0
      AVALF = 0.
      AVZETA = 0.
      SUMALF = 0.
      SUMZET = 0.
      NCHNG = 0
      NLNCR = 0

      numTotPoints = 0 !yma

      V1R4ST = V1R4
      V2R4ST = V2R4 !V2R4 will be modified by LBLF4Q
      IF (ILBLF4.GE.1) CALL LBLF4Q( JRAD,V1R4,V2R4, &
                                    R4,RR4,NPTR4,DVR4,BOUND4,TAVE,DPTFAC,DPTMIN,NLTE,&
                                    LNFIL4,NEGEPP_FLAG,ILNFLG,JCNVF4,ILIN4,ILIN4T )
      IFPAN = 1

   60 CONTINUE

         !CALL CPUTIM (TIME0)
         IF (IEOF.NE.0) GO TO 80

         ! THERE ARE (LIMIN * 9) QUANTITIES READ IN:
         ! VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG

         CALL RDLIN( VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG, &
                     BRD_MOL_FLG,BRD_MOL_HW,BRD_MOL_TMP,BRD_MOL_SHFT, &
                     IOUT, LINFIL,NLNGTH,ILO,IHI,LIMIN,VBOT,VTOP, &
                     IDATA,IEOF,NOPR )

         !CALL CPUTIM (TIME)
         !TIMRDF = TIMRDF+TIME-TIME0

         IF (IEOF.NE.0) GO TO 80

         ! MODIFY LINE DATA FOR TEMPERATURE, PRESSURE, AND COLUMN DENSITY
         !CALL CPUTIM(TPAT0)
         CALL LNCORQ( NLNCR,IHI,ILO,MEFDP, &
                      VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG,&
                      SPPSP,RECALF,ZETAI,IZETA,&
                      BRD_MOL_FLG,BRD_MOL_HW,BRD_MOL_TMP,BRD_MOL_SHFT, &
                      IOUT,DV,ALFMAX,ISOTPL_FLAG,&
                      PAVE,TAVE,WK,WKI,WTOT,NUMSTATE,RATSTATE,&
                      VFT,HWF3,ILNFLG,DVR4,V1R4,NPTR4,R4,DPTFC,DPTMN,&
                      NCHNG,NLIN,NMINUS,NPLUS,SUMALF,SUMZET,LINCNT,&
                      lnMolNames,NLTE,speciesBroad,JRAD, &
                      numNarrowLines, narrowWidth )
         !CALL CPUTIM(TPAT1)
         !TLNCOR = TLNCOR+TPAT1-TPAT0

   70    CONTINUE


         CALL CNVFNQ( VNU,SP,EPP,SPPSP,RECALF, & ! SABS=>SP,SRAD=>EPP
                      R1,R2,R3,RR1,RR2,RR3,lbR1,lbR2,lbR3, F1,F2,F3,FG,XVER, &
                      ZETAI,IZETA, &
                      DV,DVR2,DVR3,HWF1,DXF1,NX1,NX2,NX3,&
                      VFT,ILO,IHI,MAX1,IOUT,IPANEL,IDATA )


      IF (IPANEL.EQ.0) GO TO 60 !the record of lines used up, go to read more lines

   80 CONTINUE

      !    FOR FIRST PANEL     N1R1=   1    N1R2=  1    N1R3=  1
      ! FOR SUBSEQUENT PANELS  N1R1=  33   *N1R2= 13   *N1R3=  6
      !     FOR ALL PANELS     N2R1=2432   *N2R2=612   *N2R3=155
      !
      !        NOTE: THE VALUES FOR N1R2, N1R3, N2R2 AND N2R3 WHICH
      !              ARE MARKED WITH AN ASTERISK, CONTAIN A 4 POINT
      !              OFFSET WHICH PROVIDES THE NECESSARY OVERLAP FOR
      !              THE INTERPOLATION OF R3 INTO R2, AND R2 INTO R1.

      IF (IFPAN.EQ.1) THEN
         IFPAN = 0
         N1R1 = 1
         N1R2 = 1
         N1R3 = 1
      ELSE
         N1R1 = NSHIFT+1
         N1R2 = (NSHIFT/4)+1+4
         N1R3 = (NSHIFT/16)+1+3
      ENDIF
      N2R1 = NLIM1+NSHIFT-1
      N2R2 = NLIM2+(NSHIFT/4)-1+4
      N2R3 = NLIM3+(NSHIFT/16)-1+3

      IF (VFT.LE.0.) THEN
         CALL RSYM (R1,lbR1,DV,VFT)
         CALL RSYM (R2,lbR2,DVR2,VFT)
         CALL RSYM (R3,lbR3,DVR3,VFT)
         if ( NLTE ) then
            CALL RSYM (RR1,lbR1,DV,VFT)
            CALL RSYM (RR2,lbR2,DVR2,VFT)
            CALL RSYM (RR3,lbR3,DVR3,VFT)
         endif
      ENDIF

      IF (.not.XsectOff.AND.IR4.EQ.0) THEN
         !CALL CPUTIM (TIME0)
         CALL XSECTM( IFST,IR4,&
                      R1,R2,R3,R4,lbR1,lbR2,lbR3,lbR4,&
                      V1,V2,V1R4,V2R4,DV,DVR2,DVR3,DVR4,&
                      N1R1,N2R1,N1R2,N2R2,N1R3,N2R3,NPTR4,NHI,VFT, &
                      PAVE,TAVE,WXM,&
                      XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                      IXSBIN,JRAD,ILBLF4,DPTMIN,NFHDRF,NPHDRF )
         !CALL CPUTIM (TIME)
         !TXS = TXS+TIME-TIME0
      ENDIF

      !CALL CPUTIM(TPAT0)
      IF (ILBLF4.GE.1) THEN
         CALL XINT (V1R4,V2R4,DVR4,R4,1.0,VFT,DVR3,R3(1:),N1R3,N2R3)
         if ( NLTE ) then
            CALL XINT (V1R4,V2R4,DVR4,RR4,1.0,VFT,DVR3,RR3(1:),N1R3,N2R3)
         endif
      ENDIF

      IF (.not.contnmOff)                                                  &
     &    CALL XINT (V1ABS,V2ABS,DVABS,ABSRB,1.,VFT,DVR3,R3(1:),N1R3,N2R3)
      !CALL CPUTIM(TPAT1)
      !TXINT = TXINT + TPAT1-TPAT0

      CALL panelOut( R1,R2,R3,RR1,RR2,RR3,lbR1,lbR2,lbR3, JRAD, &
                     NLTE, r1Buffer, rr1Buffer, numTotPoints, &
                     VFT,V2,DV,DVP,NLO,NHI,&
                     NLIM1,NLIM2,NLIM3,MAX1,MAX2,MAX3,&
                     N1R2,NSHIFT,TAVE,ISTOP )

      IF (ISTOP.NE.1) THEN
         IF (ILBLF4.GE.1) THEN

            VF1 = VFT-2.*DVR4
            VF1 = V1+floor((VF1-V1)/DVR4)*DVR4
            VF2 = VFT+2.*DVR4+ REAL(N2R3+4)*DVR3

            IF (VF2.GT.V2R4.AND.V2R4.NE.V2R4ST) THEN

               V1R4 = VF1
               V2R4 = V2R4ST !V2R4 will modified in LBLF4Q
               CALL LBLF4Q( JRAD,V1R4,V2R4, &
                            R4,RR4,NPTR4,DVR4,BOUND4,TAVE,DPTFAC,DPTMIN,NLTE,&
                            LNFIL4,NEGEPP_FLAG,ILNFLG,JCNVF4,ILIN4,ILIN4T )

               IF (.not.XsectOff.AND.IR4.EQ.1) THEN
                  !CALL CPUTIM (TIME0)
                  CALL XSECTM( IFST,IR4,&
                               R1,R2,R3,R4,lbR1,lbR2,lbR3,lbR4,&
                               V1,V2,V1R4,V2R4,DV,DVR2,DVR3,DVR4,&
                               N1R1,N2R1,N1R2,N2R2,N1R3,N2R3,NPTR4,NHI,VFT, &
                               PAVE,TAVE,WXM,&
                               XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                               IXSBIN,JRAD,ILBLF4,DPTMIN,NFHDRF,NPHDRF )
                  !CALL CPUTIM (TIME)
                  !TXS = TXS+TIME-TIME0
               ENDIF
            ENDIF
         ENDIF
         GO TO 70
      ENDIF

      !CALL CPUTIM (TIMEH1)
      !TIME = TIMEH1-TIMEH0-TF4-TXS

      IF (NOPR.NE.1) THEN
         IF (ILBLF4.GE.1) WRITE (IPR,905) DVR4,BOUND4
         IF (NMINUS.GT.0) WRITE (IPR,910) NMINUS
         IF (NPLUS.GT.0) WRITE (IPR,915) NPLUS
         WRITE (IPR,*) "NLIN,LINCNT,NCHNG=", NLIN,LINCNT,NCHNG

         IF (LINCNT.GE.1) THEN
            AVALF = SUMALF/ REAL(LINCNT)
            AVZETA = SUMZET/ REAL(LINCNT)
         ENDIF
         WRITE (IPR,925) AVALF,AVZETA

         DO 90 M = 1, size(lnMolNames)
            !IF (MEFDP(M).GT.0) WRITE (IPR,930) MEFDP(M),M
            imol = molNum( lnMolNames(m) )
            IF (MEFDP(imol).GT.0) WRITE (IPR,930) MEFDP(imol),imol
   90    CONTINUE
      ENDIF

      deallocate( R1 )
      deallocate( R2 )
      deallocate( R3 )
      deallocate( R4 )

      RETURN

  900 FORMAT ('0  * HIRAC1 *  OUTPUT ON FILE ',I5,10X,' DV = ',F12.8,   &
     &        10X,' BOUNDF3(CM-1) = ',F8.4)
  905 FORMAT ('0 DV FOR LBLF4 = ',F10.5,5X,'BOUND FOR LBLF4 =',F10.4)
  910 FORMAT ('0 -------------------------',I5,' HALF WIDTH CHANGES')
  915 FORMAT ('0 +++++++++++++++++++++++++',I5,' HALF WIDTH CHANGES')
  920 FORMAT ('0',20X,'TIME',11X,'READ',4X,'CONVOLUTION',10X,'PANEL',   &
     &        9X,'OTHER+',                                              &
     &        6X,'NO. LINES',3X,'AFTER REJECT',5X,'HW CHANGES',/,       &
     &        2x,'LINF4',3X,2F15.3,15X,2F15.3,2I15,/,                   &
     &        2X,'XSECT ',2X,4F15.3,/,2X,'LBLF4 ',2X,4F15.3,15X,2I15,/, &
     &        2X,'HIRAC1',2X,5F15.3,3I15)
  921 FORMAT (2x,'LINF4',3X,2F15.3,15X,2F15.3,                          &
     &        2X,'XSECT ',2X,4F15.3,                                    &
     &        2X,'LBLF4 ',2X,4F15.3,                                    &
     &        2X,'HIRAC1',2X,5F15.3)
  922 FORMAT ('0',20X,'TIME',11X,'READ',4X,'CONVOLUTION',10X,'PANEL',   &
     &        9X,'OTHER+',/,                                            &
     &        2x,'LINF4',3X,2F15.3,15X,2F15.3,/,                        &
     &        2X,'XSECT ',2X,4F15.3,/,2X,'LBLF4 ',2X,4F15.3,15X,/,      &
     &        2X,'HIRAC1',2X,5F15.3)
  925 FORMAT ('0  * HIRAC1 *  AVERAGE WIDTH = ',F8.6,                   &
     &        ',  AVERAGE ZETA = ',F8.6)
  930 FORMAT ('0 ********  HIRAC1  ********',I5,' STRENGTHS FOR',       &
     &        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,            &
     &        ' SET TO ZERO')
  935 FORMAT (/,'0     + OTHER timing includes:',/,                     &
     &          '0             In LINF4:  MOLEC, BUFIN, BUFOUT, ',      &
     &          'NWDL, ENDFIL, and SHRINQ',/,                           &
     &          '0             In HIRAC:  LNCOR, XINT, SHAPEL, ',       &
     &          'SHAPEG, VERFN, MOLEC, and other loops and ',           &
     &          'file maintenance within HIRAC',/)

      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE LNCORQ( NLNCR,IHI,ILO,MEFDP, &
                         VNU,SP,ALFA0,SRAD,MOL,HWHMS,TMPALF,PSHIFT,IFLG,&
                         SPPSP,RECALF,ZETAI,IZETA,&
                         BRD_MOL_FLG,BRD_MOL_HW,BRD_MOL_TMP,BRD_MOL_SHFT, &
                         IOUT,DV,ALFMAX,ISOTPL_FLAG,&
                         PAVE,TAVE,WK,WKI,WTOT,NUMSTATE,RATSTATE,&
                         VFT,HWF3,ILNFLG,DVR4,V1R4,NPTR4,R4,DPTFC,DPTMN,&
                         NCHNG,NLIN,NMINUS,NPLUS,SUMALF,SUMZET,LINCNT,&
                         molNames,NLTE,speciesBroad,JRAD,&
                         numNarrowLines, narrowWidth )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, RADCN2, ONEPL, &
                                         MXMOL,MAX_ISO,MXISOTPL,MaxISOTPL_smass,MAXSTATE,&
                                         P0=>Press1013, TEMP0=>Temp296, &
                                         molNum, xsMolNum
      USE Module_LineData         ,ONLY: MXBRDMOL, NLINEREC
      USE Module_ODLAY_CommonSub  ,ONLY: AVRAT, MOLEC, RADFN, line_exception
      USE Module_Config           ,ONLY: IPR
      USE Module_XSect            ,ONLY: xsTbl

      IMPLICIT NONE !REAL*8           (V)

      integer      ,intent(inout) :: NLNCR
      integer      ,intent(in)    :: IHI
      integer      ,intent(in)    :: ILO
      integer      ,intent(inout) :: MEFDP(:)         !(64)
      !
      real(r8)     ,intent(inout) :: VNU(:)           !(NLINEREC)
      real         ,intent(inout) :: SP(:)            !(NLINEREC) !carry in S and then store values for SP
      real         ,intent(in)    :: ALFA0(:)         !(NLINEREC)
      real         ,intent(inout) :: SRAD(:)          !(NLINEREC) !carry in EPP and then store values for SRAD
      integer*4    ,intent(inout) :: MOL(:)           !(NLINEREC)
      real         ,intent(in)    :: HWHMS(:)         !(NLINEREC)
      real         ,intent(in)    :: TMPALF(:)        !(NLINEREC)
      real         ,intent(in)    :: PSHIFT(:)        !(NLINEREC)
      integer      ,intent(in)    :: IFLG(:)          !(NLINEREC)
      real         ,intent(out)   :: SPPSP(:)         !(NLINEREC)
      real         ,intent(out)   :: RECALF(:)        !(NLINEREC)
      real         ,intent(out)   :: ZETAI(:)         !(NLINEREC)
      integer      ,intent(out)   :: IZETA(:)         !(NLINEREC)
      integer*4    ,intent(in)    :: BRD_MOL_FLG(:,:) !(MXBRDMOL,NLINEREC)
      real         ,intent(in)    :: BRD_MOL_HW(:,:)  !(MXBRDMOL,NLINEREC)
      real         ,intent(in)    :: BRD_MOL_TMP(:,:) !(MXBRDMOL,NLINEREC)
      real         ,intent(in)    :: BRD_MOL_SHFT(:,:)!(MXBRDMOL,NLINEREC)
      integer      ,intent(in)    :: IOUT(:)          !(NLINEREC)
      real         ,intent(in)    :: DV
      real         ,intent(in)    :: ALFMAX
      integer      ,intent(in)    :: ISOTPL_FLAG(:,:) !(MXMOL,MXISOTPL)
      real         ,intent(in)    :: PAVE, TAVE
      real         ,intent(in)    :: WK(:)            !(60)
      real         ,intent(in)    :: WKI(:,:)         !(MXMOL,MXISOTPL)
      real         ,intent(in)    :: WTOT
      integer      ,intent(in)    :: NUMSTATE(:)      !(MXMOL)
      real         ,intent(in)    :: RATSTATE(:,:)    !(MAXSTATE*MAX_ISO,MXMOL)
      real(r8)     ,intent(in)    :: VFT
      real         ,intent(in)    :: HWF3
      integer      ,intent(in)    :: ILNFLG
      real         ,intent(in)    :: DVR4
      real(r8)     ,intent(in)    :: V1R4
      integer      ,intent(in)    :: NPTR4
      real         ,intent(in)    :: R4(:)
      real         ,intent(in)    :: DPTFC
      real         ,intent(in)    :: DPTMN
      integer      ,intent(inout) :: NCHNG
      integer      ,intent(inout) :: NLIN
      integer      ,intent(inout) :: NMINUS
      integer      ,intent(inout) :: NPLUS
      real         ,intent(inout) :: SUMALF
      real         ,intent(inout) :: SUMZET
      integer      ,intent(inout) :: LINCNT
      character(*) ,intent(in)    :: molNames(:)
      logical      ,intent(in)    :: NLTE
      logical      ,intent(in)    :: speciesBroad
      integer      ,intent(in)    :: JRAD
      integer      ,intent(in)    :: numNarrowLines !number of narrow lines counted in LINF4. used here as an indicator to turn on/off skipping narrow lines
      real         ,intent(in)    :: narrowWidth


      !--- Local variables
      !
      character(*) ,parameter :: routineName = 'LNCORQ'

      real*4        :: AMOL
      real    ,SAVE :: XKT
      real    ,SAVE :: DELTMP !no use
      real    ,SAVE :: BETACR
      real    ,SAVE :: SCOR(MXMOL,MXISOTPL) !(42,10)
      real    ,SAVE :: RHOSLF(MXMOL)
      real    ,SAVE :: ALFD1(MXMOL,MXISOTPL)!(42,10)
      real    ,SAVE :: TRATIO
      real    ,SAVE :: RHORAT
      real    ,SAVE :: PAVP0
      real    ,SAVE :: PAVP2
      integer ,SAVE :: ILC
      real    ,SAVE :: RECTLC
      real    ,SAVE :: TMPDIF

      ! TEMPERATURES FOR LINE COUPLING COEFFICIENTS
      real        ,PARAMETER :: TEMPLC(4)=[ 200.0, 250.0, 296.0, 340.0 ]
      character   ,PARAMETER :: HREJ='0', HNOREJ='1'
      integer     ,PARAMETER :: NWDTH=0
      integer     ,PARAMETER :: I_1=1, I_100=100, I_1000=1000
      character*8 ,PARAMETER :: h_lncor1=' lncor1 '

      CHARACTER*1  :: FREJ(nlinerec)
      real         :: TMPCOR_ARR(MXBRDMOL), ALFA_TMP(MXBRDMOL)
      Real(r8)     :: ALFSUM
      real         :: A(4), B(4)

      INTEGER     :: I,        IFLAG,     IL,       INDLOW
      INTEGER     :: INDUPP,   ISO,       IZ
      INTEGER     :: J,        JJ
      INTEGER     :: M,        MFULL
      INTEGER     :: NLOW,     NMINAD,    NPLSAD,   NUPP
      REAL        :: ALFA0I,   ALFAD,    ALFL
      REAL        :: ALFV,     DELTA,    FNLTE
      REAL        :: FREQ,     FZETA,     GAMMA1,   GAMMA2
      REAL        :: GI,       HWHMSI,    RLOW
      REAL        :: RUPP,     SLFABS,    SLOPEA,   SLOPEB
      REAL        :: SPEAK,    SPPI,      SUI
      REAL        :: TEST,     TMPCOR,    XKT0,     YI
      REAL        :: ZETA,     ZETDIF

      integer :: im, is, ind
      real    :: v1x,v2x
      integer :: moNum( size(molNames) )
      integer :: xsNum( size(molNames) )



      NLNCR = NLNCR+1
      IF (NLNCR.EQ.1) THEN

         XKT0 = TEMP0/RADCN2
         XKT = TAVE/RADCN2
         DELTMP = ABS(TAVE-TEMP0)
         BETACR = (1./XKT)-(1./XKT0)
         CALL MOLEC( 2,SCOR,RHOSLF,ALFD1, &
                     PAVE,TAVE,WK,WTOT,molNames )

         TRATIO = TAVE/TEMP0
         RHORAT = (PAVE/P0)*(TEMP0/TAVE)

         PAVP0 = PAVE/P0
         PAVP2 = PAVP0*PAVP0

         ! FIND CORRECT TEMPERATURE AND INTERPOLATE FOR Y AND G

         DO IL = 1, 3
            ILC = IL
            IF (TAVE.LT.TEMPLC(ILC+1)) GO TO 20
         ENDDO
   20    IF (ILC.EQ.4) ILC = 3

         RECTLC = 1.0/(TEMPLC(ILC+1)-TEMPLC(ILC))
         TMPDIF = TAVE-TEMPLC(ILC)

      ENDIF

      IF (ILNFLG.EQ.2) READ(15)(FREJ(J),J=ILO,IHI)


      do im=1,size(molNames)
         moNum(im) = molNum( molNames(im) )
         xsNum(im) = xsMolNum( molNames(im) )
      enddo


      DO 30 J = ILO, IHI
         YI = 0.
         GI = 0.
         GAMMA1 = 0.
         GAMMA2 = 0.
         I = IOUT(J)
         IFLAG = IFLG(I)
         MFULL=MOL(I)
         ! Molecule number for this line, last 2 digits of MFULL
         M = MOD(MOL(I),I_100)
         ! ISO=(MOD(MOL(I),1000)-M)/100   IS PROGRAMMED AS:
         ! Isotope number for this line, 3rd digit from right of MFULL
         ISO = MOD(MOL(I),I_1000)/100
         if (iso.eq.0) iso=10

         !--- Check if this line is from a selected molecule and when
         ! both line data and xs data are available if this line prefer xs over line,
         ! if so, skip this line.
         !
         if ( .not.any(moNum==m) ) then !This line is not from selected molecules
            GO TO 25
         else  !line is from a selected mol, check if it is also a xs molecule.

            ind = minloc( abs(moNum-m), 1 )
            if ( xsNum(ind)>0 ) then !this molecule is also a xs molecule and has xs data in V1~V2

               !This line is from a xs molecule, unless the lineOrXs marker ==1, the absorption
               !is handled using x-section data. This assums xs data covers all possible intervals.
               if( all( xsTbl%lineOrXs(:,xsNum(ind)) /=1) ) GO TO 25

               !At least one xs spectral interval is with lineOrXs marker==1
               do is = 1,xsTbl%numSpect( xsNum(ind) )
                  v1x = xsTbl%V1FX(is,xsNum(ind))
                  v2x = xsTbl%V2FX(is,xsNum(ind))
                  if ( v1x < VNU(I).and.VNU(I) < v2x ) then
                     if( xsTbl%lineOrXs(is, xsNum(ind)) /=1 ) then
                        GO TO 25 !This interval handeled using xs section
                     else
                        EXIT !this interval handeled as line absorption.
                     endif
                  endif
               enddo
               if (is > xsTbl%numSpect( xsNum(ind) ) ) then
                  STOP '---'//routineName//'(): Found a line from xs molecule but not in any xs intervals listed in FSCDXS.'
               endif
            endif

            !if ( xsNum(ind)>0 ) then !this molecule is also a xs molecule and has xs data in V1~V2
            !   do is = 1,xsTbl%numSpect( xsNum(ind) )
            !      if ( xsTbl%lineOrXs(is, xsNum(ind)) /=2 ) CYCLE  !this xs sub interval doesnot have line data (lineOrXs=0), or it does and treaded line-by-line (lineOrXs=1)
            !      v1x = xsTbl%V1FX(is,xsNum(ind))
            !      v2x = xsTbl%V2FX(is,xsNum(ind))
            !      if ( v1x < VNU(I).and.VNU(I) < v2x ) then
            !         GO TO 25 !this line lies in a xs sub interval and is treaded by xs
            !      endif
            !   enddo
            !endif
         endif

         ! check if lines are within allowed molecular and isotopic limits
         if (m.gt.mxmol .or. m.lt. 1) then
            call line_exception (1,ipr,h_lncor1,m,mxmol,iso,MaxISOTPL_smass)
            go to 25
         else if (iso .gt. MaxISOTPL_smass(m)) then
            call line_exception (2,ipr,h_lncor1,m,mxmol,iso,MaxISOTPL_smass)
            go to 25
         endif

         MOL(I) = M

         IF (ISOTPL_FLAG(M,ISO).EQ.0) THEN
            SUI = SP(I)*WK(M)  !yma SUI = S(I)*WK(M)
         ELSE
            SUI = SP(I)*WKI(M,ISO)  !yma SUI = S(I)*WKI(M,ISO)
         ENDIF

         !MJA, 20150821
         ! using the VNU(I) approximation for the radiation term
         ! can cause issues for wavenumbers around 1000 cm-1,
         ! so use the full rad term instead
         ! IF (JRAD.EQ.1) SUI = SUI*VNU(I)
         IF (JRAD.EQ.1) SUI = SUI*RADFN(VNU(I),XKT)

         IF (SUI.EQ.0.) GO TO 25

         NLIN = NLIN+1

         ! Y'S AND G'S ARE STORED IN I+1 POSTION OF VNU,S,ALFA0,EPP...
         ! A(1-4),  B(1-4) CORRESPOND TO TEMPERATURES TEMPLC(1-4) ABOVE

         IF (IFLAG.EQ.1.OR.IFLAG.EQ.3) THEN
            A(1) = VNU(I+1)
            B(1) = SP(I+1)                   !yma B(1) = S(I+1)
            A(2) = ALFA0(I+1)
            B(2) = SRAD(I+1)                 !yma B(2) = EPP(I+1)
            A(3) = transfer( MOL(I+1),AMOL ) !yma A(3) = AMOL(I+1)
            B(3) = HWHMS(I+1)
            A(4) = TMPALF(I+1)
            B(4) = PSHIFT(I+1)

            ! CALCULATE SLOPE AND EVALUATE

            SLOPEA = (A(ILC+1)-A(ILC))*RECTLC
            SLOPEB = (B(ILC+1)-B(ILC))*RECTLC

            IF (IFLAG.EQ.1) THEN
               YI = A(ILC)+SLOPEA*TMPDIF
               GI = B(ILC)+SLOPEB*TMPDIF
            ELSE
               GAMMA1 = A(ILC)+SLOPEA*TMPDIF
               GAMMA2 = B(ILC)+SLOPEB*TMPDIF
            ENDIF
         ENDIF

         ! IFLAG = 2 IS RESERVED FOR LINE COUPLING COEFFICIENTS ASSOCIATED
         !           WITH AN EXACT TREATMENT (NUMERICAL DIAGONALIZATION)
         ! IFLAG = 3 TREATS LINE COUPLING IN TERMS OF REDUCED WIDTHS

         VNU(I) = VNU(I)+RHORAT*PSHIFT(I)
         if(sum(brd_mol_flg(:,i)).gt.0.AND. speciesBroad) then
            vnu(i) = vnu(i)+sum( rhoslf(1:mxbrdmol)*brd_mol_flg(:,i)* &
     &                           (brd_mol_shft(:,i)-pshift(i)) )
         endif

         ! TEMPERATURE CORRECTION OF THE HALFWIDTH
         ! SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN

         TMPCOR = TRATIO**TMPALF(I)
         ALFA0I = ALFA0(I)*TMPCOR
         HWHMSI = HWHMS(I)*TMPCOR
         ALFL = ALFA0I*(RHORAT-RHOSLF(m))+HWHMSI*RHOSLF(m)

         if(sum(brd_mol_flg(:,i)).gt.0.AND. speciesBroad) then
            tmpcor_arr = tratio**brd_mol_tmp(:,i)
            alfa_tmp = brd_mol_hw(:,i)*tmpcor_arr
            alfsum = sum( rhoslf(1:mxbrdmol)*brd_mol_flg(:,i)*alfa_tmp )
            alfl = ( rhorat-sum(rhoslf(1:mxbrdmol)*brd_mol_flg(:,i)) )&
                   *alfa0i + alfsum
            if(brd_mol_flg(m,i).eq.0)   &
     &           alfl = alfl + rhoslf(m)*(hwhmsi-alfa0i)
         end if

! mji - Revise code to skip broadening for incoming lines with flag = -1
!      if(brd_mol_flg(m,i).eq.-1.and.ibrd.gt.0) then
! bb - I think ibrd in LBL is speciesBroad in clblm
         if (m.le.mxbrdmol) then
            if(brd_mol_flg(m,i).eq.-1.and.speciesBroad) then
               alfl = hwhmsi
            end if
         endif

         IF (IFLAG.EQ.3) ALFL = ALFL*(1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)

         ALFAD = VNU(I)*ALFD1(m,iso)
         ZETA = ALFL/(ALFL+ALFAD)
         ZETAI(I) = ZETA
         FZETA = 100.*ZETA
         IZ = FZETA + ONEPL
         IZETA(I) = IZ
         ZETDIF = FZETA - REAL(IZ-1)

         ALFV = (AVRAT(IZ)+ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)))*(ALFL+ALFAD)

         IF (ALFV.LT.DV) THEN
            ALFV = DV
            NMINAD = 1
         ELSE
            NMINAD = 0
         ENDIF
         IF (ALFV.GT.ALFMAX) THEN
            ALFV = ALFMAX
            NPLSAD = 1
         ELSE
            NPLSAD = 0
         ENDIF

!------------------------>>>
         !--- Narrow lines will be treated in a second round
         !  convolution, skip them for this round
         if ( numNarrowLines>0 .and. ALFV<narrowWidth ) GO TO 25
!------------------------>>><<<

         IF (HWF3*ALFV+VNU(I) .LT. VFT) GO TO 25

         RECALF(I) = 1./ALFV

         !   TREAT TRANSITIONS WITH negative EPP AS SPECIAL CASE
         !>> an epp value between -0.9999 and 0.  cm-1 is taken as valid
         !>> an epp value of -1. is assumed set by hitran indicating an unknown
         !   value: no temperature correction is performed
         !>> for an epp value of less than -1., it is assumed that value has
         !   been provided as a reasonable value to be used for purposes of
         !   temperature correction.  epp is set positive

         if (SRAD(i).le.-1.001) SRAD(i) = abs(SRAD(i))   !yma if (epp(i).le.-1.001) epp(i) = abs(epp(i))

         if (SRAD(i).le.-0.999) MEFDP(M) = MEFDP(M)+1    !yma if (epp(i).le.-0.999) MEFDP(M) = MEFDP(M)+1

         ! temperature correction:
         if (SRAD(i) .gt. -0.999) then                   !yma if (epp(i) .gt. -0.999) then
            SUI = SUI*SCOR(m,iso)* EXP(-SRAD(I)*BETACR)*(1.+EXP(-VNU(I)/XKT) )  !yma SUI = SUI*SCOR(m,iso)* EXP(-EPP(I)*BETACR)*(1.+EXP(-VNU(I)/XKT) )
         endif


         SP(I) = SUI*(1.+GI*PAVP2)
         SPPI = SUI*YI*PAVP0
         SPPSP(I) = SPPI/SP(I)
         SRAD(I)=0.0

! ---from nlte:
!yma         IF (MFULL.GE.1000) THEN
         IF ( NLTE .and. MFULL.GE.1000) THEN
            FREQ=VNU(I)
            ! NLOW is 4th and 5th digit from right of MFULL
            NLOW=MOD(MFULL/1000,100)
            ! NUPP is 6th and 7th digit from right of MFULL
            NUPP=MFULL/100000
            RLOW=1.0
            RUPP=1.0
            DELTA=EXP(-FREQ/XKT)

            ! PICK OUT MOLECULAR TYPES WITH VIBRATIONAL STATES

            IF(NUMSTATE(M).GT.0) THEN
             IF (NLOW.GT.NUMSTATE(M)) STOP 'NLOW GT NUMSTATE IN LNCORQ'
               INDLOW=NLOW + (ISO-1)*MAXSTATE
               IF (NLOW.GT.0) RLOW=RATSTATE(INDLOW,M)
             IF (NUPP.GT.NUMSTATE(M)) STOP 'NUPP GT NUMSTATE IN LNCORQ'
               INDUPP=NUPP + (ISO-1)*MAXSTATE
               IF (NUPP.GT.0) RUPP=RATSTATE(INDUPP,M)
            ELSE
               PRINT 900,M
  900        FORMAT('LNCORQ: MOL IN TROUBLE',I10)
               SP(I)=0.
               SRAD(I)=0.
               SPPSP(I)=0.
               GO TO 30
            END IF

            ! RLOW AND RUPP NOW SET

            FNLTE=SP(I)/(1.0-DELTA)
            SP(I)=FNLTE*(RLOW-RUPP*DELTA)
            SRAD(I)=FNLTE*(RLOW-RUPP)

            IF (IFLAG .EQ. 0) THEN
               SPEAK = SP(I)*ABS(RECALF(I))
               SLFABS = SPEAK
               IF(SPEAK.GE.5.) THEN
                  SLFABS = 1.
               ELSE
                  IF(SPEAK.GT.0.01) SLFABS = 1.-EXP(-SPEAK)
               ENDIF
               TEST = SLFABS *(1.-SRAD(I)/SP(I))
               IF(TEST.LE.DPTMN) GOTO 25
            END IF

         ! --- above from nlte
         ELSE

            IF (IFLAG.EQ.0) THEN
               IF (ILNFLG==0 .or. ILNFLG==1) THEN  !!IF (ILNFLG.LE.1) THEN
                  FREJ(J) = HNOREJ
                  SPEAK = SUI*RECALF(I)
                  IF (DVR4.LE.0.) THEN
                     IF (SPEAK.LE.DPTMN) THEN
                        FREJ(J) = HREJ
                        GO TO 25
                     ENDIF
                  ELSE
                     JJ = (VNU(I)-V1R4)/DVR4+1.
                     JJ = MAX(JJ,1)
                     JJ = MIN(JJ,NPTR4)
!                         IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ))) THEN
                     IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) .and. sppsp(i)   &
                     .eq.0.) THEN
                        FREJ(J) = HREJ
                        GO TO 25
                     ENDIF
                  ENDIF
               ELSE if (ILNFLG==2) then
               ! "ELSE" IS TRUE WHEN "ILNFLG" EQUALS 2

                  IF (FREJ(J).EQ.HREJ) GO TO 25
               ENDIF
            ENDIF
         ENDIF

         NMINUS = NMINUS+NMINAD
         NPLUS = NPLUS+NPLSAD
         SUMALF = SUMALF+ALFV
         SUMZET = SUMZET+ZETA
         LINCNT = LINCNT+1

         GO TO 30

   25    SP(I)=0.0
         SRAD(I)=0.0
         SPPSP(I) = 0.0

   30 END DO

      NCHNG = NMINUS+NPLUS
      IF (ILNFLG.EQ.1) WRITE(15)(FREJ(J),J=ILO,IHI)


      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE CNVFNV PERFORMS THE CONVOLUTION OF THE LINE DATA WITH
!     THE VOIGT LINE SHAPE (APPROXIMATED)
!
!     IMPLEMENTATION:    R.D. WORSHAM
!
!     ALGORITHM REVISIONS:    S.A. CLOUGH
!     R.D. WORSHAM
!     J.L. MONCET
!
!
!     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!     WORK SUPPORTED BY:    THE ARM PROGRAM
!     OFFICE OF ENERGY RESEARCH
!     DEPARTMENT OF ENERGY
!
!
!     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!     FASCOD3
!
!-----------------------------------------------------------------------
      SUBROUTINE CNVFNQ( VNU,SABS,SRAD,SPPSP,RECALF,&
                         R1,R2,R3,RR1,RR2,RR3,lbR1,lbR2,lbR3,F1,F2,F3,FG,XVER,&
                         ZETAI,IZETA, &
                         DV,DVR2,DVR3,HWF1,DXF1,NX1,NX2,NX3,&
                         VFT,ILO,IHI,MAX1,IOUT,IPANEL,IDATA )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8
      USE Module_ODLAY_CommonSub  ,ONLY: CF1,CF2,CF3,CER,CGAUSS
      IMPLICIT NONE !REAL*8           (V)

      real(r8)  ,intent(in)    :: VNU(:)
      real      ,intent(in)    :: SABS(:)
      real      ,intent(in)    :: SRAD(:)
      real      ,intent(in)    :: SPPSP(:)
      real      ,intent(in)    :: RECALF(:)
      real      ,intent(inout) :: R1(lbR1:), R2(lbR2:), R3(lbR3:)
      real      ,intent(inout) :: RR1(lbR1:), RR2(lbR2:), RR3(lbR3:)
      integer   ,intent(in)    :: lbR1,lbR2,lbR3
      real      ,intent(in)    :: F1(:), F2(:), F3(:), FG(:)
      real      ,intent(in)    :: XVER(:)
      real      ,intent(in)    :: ZETAI(:)
      integer   ,intent(in)    :: IZETA(:)
      !
      real      ,intent(in)    :: DV, DVR2, DVR3
      real      ,intent(in)    :: HWF1
      real      ,intent(in)    :: DXF1
      integer   ,intent(in)    :: NX1, NX2, NX3
      real(r8)  ,intent(in)    :: VFT
      integer   ,intent(inout) :: ILO,IHI
      integer   ,intent(in)    :: MAX1
      integer   ,intent(in)    :: IOUT(250)
      integer   ,intent(out)   :: IPANEL
      integer   ,intent(in)    :: IDATA


      !--- Local variables
      INTEGER  :: I,      ILAST,  IZ1,    IZ2
      INTEGER  :: IZ3,    IZM,    J,      J1
      INTEGER  :: J2,     J2SHFT, J3,     J3SHFT
      INTEGER  :: JMAX1,  JMIN1,  JMIN2,  JMIN3
      REAL     :: BHWDXF, CLC1,   CLC2,   CLC3
      REAL     :: CONF2,  CONF3,  DEPTHA, DEPTHR
      REAL     :: DPTRAT, HWDXF,  RSHFT,  STRDA
      REAL     :: STRDR,  STRFA1, STRFA2, STRFA3
      REAL     :: STRFR1, STRFR2, STRFR3, STRVRA
      REAL     :: STRVRR, WAVDXF, ZETDIF
      REAL     :: ZF1,    ZF1L,   ZF2,    ZF2L
      REAL     :: ZF3,    ZF3L,   ZINT,   ZSLOPE


      CLC1 = 4./( REAL(NX1-1))
      CLC2 = 16./( REAL(NX2-1))
      CLC3 = 64./( REAL(NX3-1))
      WAVDXF = DV/DXF1
      HWDXF = HWF1/DXF1
      CONF2 = DV/DVR2
      CONF3 = DV/DVR3
      ILAST = ILO-1

      IF (ILO.LE.IHI) THEN
         DO 30 J = ILO, IHI
            I = IOUT(J)
            IF (SABS(I).NE.0.) THEN
               DEPTHA = SABS(I)*RECALF(I)
               DEPTHR = SRAD(I)*RECALF(I)
               IZM = IZETA(I)

               ZETDIF = 100.*ZETAI(I)- REAL(IZM-1)
               STRFA1 = DEPTHA*(CF1(IZM)+ZETDIF*(CF1(IZM+1)-CF1(IZM)))
               STRFA2 = DEPTHA*(CF2(IZM)+ZETDIF*(CF2(IZM+1)-CF2(IZM)))
               STRFA3 = DEPTHA*(CF3(IZM)+ZETDIF*(CF3(IZM+1)-CF3(IZM)))
               STRDA = DEPTHA*(CGAUSS(IZM)+ZETDIF*(CGAUSS(IZM+1)-       &
               CGAUSS(IZM)))
               STRVRA = DEPTHA*(CER(IZM)+ZETDIF*(CER(IZM+1)-CER(IZM)) )

               if ( abs(DEPTHR) >0. ) then !yma,151201; if DEPTHR/=0.
                  STRFR1 = DEPTHR*(CF1(IZM)+ZETDIF*(CF1(IZM+1)-CF1(IZM)))
                  STRFR2 = DEPTHR*(CF2(IZM)+ZETDIF*(CF2(IZM+1)-CF2(IZM)))
                  STRFR3 = DEPTHR*(CF3(IZM)+ZETDIF*(CF3(IZM+1)-CF3(IZM)))
                  STRDR = DEPTHR*(CGAUSS(IZM)+ZETDIF*(CGAUSS(IZM+1)-       &
                  CGAUSS(IZM)))
                  STRVRR = DEPTHR*(CER(IZM)+ZETDIF*(CER(IZM+1)-CER(IZM)) )
               endif

               ZSLOPE = RECALF(I)*WAVDXF
               ZINT = (VNU(I)-VFT)/DV
               BHWDXF = HWDXF/ZSLOPE
               JMAX1 = ZINT+BHWDXF+1.5
               IF (JMAX1.GT.MAX1) THEN
                  ILAST = J-1
                  IPANEL = 1
                  GO TO 40
               ENDIF
               JMIN1 = ZINT-BHWDXF+1.5
               RSHFT = 0.5
               IF (ZINT.LT.0.0) RSHFT = -RSHFT
               J2SHFT = ZINT*(1.-CONF2)+RSHFT
               J3SHFT = ZINT*(1.-CONF3)+RSHFT
               JMIN2 = JMIN1-J2SHFT
               JMIN3 = JMIN1-J3SHFT
               ZF1L = ( REAL(JMIN1-2)-ZINT)*ZSLOPE
               ZF2L = ( REAL(JMIN2-2)-ZINT*CONF2)*ZSLOPE
               ZF3L = ( REAL(JMIN3-2)-ZINT*CONF3)*ZSLOPE
               ZF1 = ZF1L
               ZF2 = ZF2L
               ZF3 = ZF3L
!
               DO 10 J1 = JMIN1, JMAX1
                  J2 = J1-J2SHFT
                  J3 = J1-J3SHFT
                  ZF3 = ZF3+ZSLOPE
                  ZF2 = ZF2+ZSLOPE
                  ZF1 = ZF1+ZSLOPE
                  IZ3 = ABS(ZF3)+1.5
                  IZ2 = ABS(ZF2)+1.5
                  IZ1 = ABS(ZF1)+1.5
                  R3(J3) = R3(J3)+STRFA3*F3(IZ3)
                  R2(J2) = R2(J2)+STRFA2*F2(IZ2)
                  R1(J1) = R1(J1)+STRFA1*F1(IZ1)+STRDA*FG(IZ1) +STRVRA* &
                  XVER(IZ1)
                  IF (DEPTHR.NE.0) THEN
                     RR3(J3) = RR3(J3)+STRFR3*F3(IZ3)
                     RR2(J2) = RR2(J2)+STRFR2*F2(IZ2)
                     RR1(J1) = RR1(J1)+STRFR1*F1(IZ1)+STRDR*FG(IZ1)     &
                     +STRVRR*XVER(IZ1)
                  ENDIF
   10          CONTINUE

               IF (SPPSP(I).NE.0.) THEN

!                 THE FOLLOWING DOES LINE COUPLING

!                 SPPSP(I) = SPP(I)/SP(I)

                  DPTRAT = SPPSP(I)
                  STRFA3 = STRFA3*CLC3*DPTRAT
                  STRFA2 = STRFA2*CLC2*DPTRAT
                  STRFA1 = STRFA1*CLC1*DPTRAT
                  STRDA = STRDA*CLC1*DPTRAT
                  STRVRA = STRVRA*CLC1*DPTRAT
                  if ( abs(DEPTHR) >0. ) then !yma, 151201; if (DEPTHR /=0)
                     STRFR3 = STRFR3*CLC3*DPTRAT
                     STRFR2 = STRFR2*CLC2*DPTRAT
                     STRFR1 = STRFR1*CLC1*DPTRAT
                     STRDR = STRDR*CLC1*DPTRAT
                     STRVRR = STRVRR*CLC1*DPTRAT
                  endif

                  DO 20 J1 = JMIN1, JMAX1

                     J2 = J1-J2SHFT
                     J3 = J1-J3SHFT
                     ZF3L = ZF3L+ZSLOPE
                     ZF2L = ZF2L+ZSLOPE
                     ZF1L = ZF1L+ZSLOPE
                     IZ3 = ABS(ZF3L)+1.5
                     IZ2 = ABS(ZF2L)+1.5
                     IZ1 = ABS(ZF1L)+1.5
                     R3(J3) = R3(J3)+STRFA3*F3(IZ3)*ZF3L
                     R2(J2) = R2(J2)+STRFA2*F2(IZ2)*ZF2L
                     R1(J1) = R1(J1)+(STRFA1*F1(IZ1)+STRDA*FG(IZ1)      &
                     +STRVRA*XVER(IZ1))*ZF1L

                     IF (DEPTHR.NE.0) THEN
                        RR3(J3) = RR3(J3)+STRFR3*F3(IZ3)*ZF3L
                        RR2(J2) = RR2(J2)+STRFR2*F2(IZ2)*ZF2L
                        RR1(J1) = RR1(J1)+(STRFR1*F1(IZ1)+ STRDR*FG(IZ1)&
                        +STRVRR*XVER(IZ1))*ZF1L
                     ENDIF

   20             CONTINUE

               ENDIF
            ENDIF

   30    CONTINUE
         ILAST = IHI

!        IDATA=0 FOR MORE DATA REQUIRED
!        IDATA=1 IF NO MORE DATA REQUIRED

         IPANEL = IDATA
      ELSE
         IPANEL = 1
      ENDIF

   40 ILO = ILAST+1
      !CALL CPUTIM (TIME)
      !TIMCNV = TIMCNV+TIME-TIME0
      RETURN

      END  SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE PANEL COMBINES RESULTS OF R3, R2, AND R1 INTO R1 ARRAY
!     AND OUTPUTS THE ARRAY R1
!
!               LAST MODIFICATION:    28 AUGUST 1992
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
      !SUBROUTINE PANELQ (R1,R2,R3,RR1,RR2,RR3,KFILE,JRAD,IENTER)
      SUBROUTINE panelOut( R1,R2,R3,RR1,RR2,RR3,lbR1,lbR2,lbR3, JRAD, &
                           NLTE, r1Buffer, rr1Buffer, numTotPoints, &
                           VFT,V2,DV,DVP,NLO,NHI,&
                           NLIM1,NLIM2,NLIM3,MAX1,MAX2,MAX3,&
                           N1R2,NSHIFT,TAVE,ISTOP )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, RADCN2
      USE Module_ODLAY_CommonSub  ,ONLY: RADFNI

      IMPLICIT NONE !REAL*8           (V)

      real      ,intent(inout) :: R1(lbR1:),R2(lbR2:),R3(lbR3:)
      real      ,intent(inout) :: RR1(lbR1:),RR2(lbR2:),RR3(lbR3:)
      integer   ,intent(in)    :: lbR1,lbR2,lbR3
      integer   ,intent(in)    :: JRAD
      logical   ,intent(in)    :: NLTE
      real      ,intent(inout) :: r1Buffer(:)
      real      ,intent(inout) :: rr1Buffer(:)
      integer   ,intent(inout) :: numTotPoints
      !
      real(r8)  ,intent(inout) :: VFT
      real(r8)  ,intent(in)    :: V2
      real      ,intent(in)    :: DV
      real      ,intent(in)    :: DVP
      integer   ,intent(inout) :: NLO
      integer   ,intent(inout) :: NHI
      integer   ,intent(in)    :: NLIM1, NLIM2, NLIM3
      integer   ,intent(in)    :: MAX1, MAX2, MAX3
      integer   ,intent(in)    :: N1R2
      integer   ,intent(in)    :: NSHIFT
      real      ,intent(in)    :: TAVE
      integer   ,intent(out)   :: ISTOP


      !--- Internal variables
      real     :: tmpRADVI
      real(r8) :: V1P, V2P
      integer  :: NLIM
      INTEGER  :: I,      J,     J2,      J3
      INTEGER  :: LIMHI,  LIMLO, NPTSI1,  NPTSI2
      REAL     :: RADVI, RDEL,    RDLAST
      REAL(r8) :: VI,     VITST
      REAL     :: X00,    X01,   X02,     X03
      REAL     :: X10,    X11,   XKT



      !CALL CPUTIM (TIME0)
      X00 = -7./128.
      X01 = 105./128.
      X02 = 35./128.
      X03 = -5./128.
      X10 = -1./16.
      X11 = 9./16.
      ISTOP = 0

!     Test for last panel.  If last, set the last point to one point
!     greater than V1 specified on TAPE5 (to ensure last point for
!     every layer is the same)

      IF ((VFT+(NHI-1)*DVP).GT.V2) THEN
         NHI = (V2-VFT)/DVP + 1.
         V2P = VFT+ REAL(NHI-1)*DVP
         IF (V2P.LT.V2) THEN
            V2P = V2P+DVP
            NHI = NHI+1
         ENDIF
         ISTOP = 1
      ELSE
         V2P = VFT+ REAL(NHI-1)*DV
      ENDIF
      NLIM = NHI-NLO+1
      V1P = VFT+ REAL(NLO-1)*DV

      LIMLO = N1R2
      IF (N1R2.EQ.1) LIMLO = LIMLO+4
      LIMHI = (NHI/4)+1


      DO 10 J = LIMLO, LIMHI, 4
         J3 = (J-1)/4+1
         R2(J) = R2(J)+R3(J3)
         R2(J+1) = R2(J+1)+X00*R3(J3-1)+ &
                           X01*R3(J3)+ &
                           X02*R3(J3+1)+ &
                           X03*R3(J3+2)
         R2(J+2) = R2(J+2)+X10*(R3(J3-1)+R3(J3+2))+ &
                           X11*(R3(J3)+R3(J3+1))
         R2(J+3) = R2(J+3)+X03*R3(J3-1)+ &
                           X02*R3(J3)+ &
                           X01*R3(J3+1)+ &
                           X00*R3(J3+2)

   10 END DO
      if ( NLTE ) then !NLTE case
         DO J = LIMLO, LIMHI, 4
            J3 = (J-1)/4+1
            RR2(J) = RR2(J)+RR3(J3)
            RR2(J+1) = RR2(J+1)+X00*RR3(J3-1)+ &
                                X01*RR3(J3)+ &
                                X02*RR3(J3+1)+ &
                                X03*RR3(J3+2)
            RR2(J+2) = RR2(J+2)+X10*(RR3(J3-1)+RR3(J3+2))+ &
                                X11*(RR3(J3)+RR3(J3+1))
            RR2(J+3) = RR2(J+3)+X03*RR3(J3-1)+ &
                                X02*RR3(J3)+ &
                                X01*RR3(J3+1)+ &
                                X00*RR3(J3+2)
         END DO
      endif

      !--- If the last panel, interpolate the first point of the next 4-DV2 segment
      ! The first point is exactly aligned, so interpolation is simply taking the
      ! corresponding R3 value.
      if (LIMLO<=LIMHI .and. ISTOP==1) then
         J3 = J3 + 1
         R2(J) = R2(J)+R3(J3)
         R2(J+1) = R2(J+1)+X00*R3(J3-1)+X01*R3(J3)+X02*R3(J3+1)+X03*R3(J3+2)
         if ( NLTE ) RR2(J) = RR2(J) + RR3(J3)
      endif

      DO 20 J = NLO, NHI, 4
         J2 = (J-1)/4+1
         R1(J) = R1(J)+R2(J2)
         R1(J+1) = R1(J+1)+X00*R2(J2-1)+ &
                           X01*R2(J2)+ &
                           X02*R2(J2+1)+ &
                           X03*R2(J2+2)
         R1(J+2) = R1(J+2)+X10*(R2(J2-1)+ &
                                R2(J2+2))+ &
                           X11*(R2(J2)+&
                                R2(J2+1))
         R1(J+3) = R1(J+3)+X03*R2(J2-1)+ &
                           X02*R2(J2)+ &
                           X01*R2(J2+1)+ &
                           X00*R2(J2+2)

   20 END DO
      if ( NLTE ) then !NLTE case
         DO J = NLO, NHI, 4
            J2 = (J-1)/4+1
            RR1(J) = RR1(J)+RR2(J2)
            RR1(J+1) = RR1(J+1)+X00*RR2(J2-1)+ &
                                X01*RR2(J2)+ &
                                X02*RR2(J2+1)+ &
                                X03*RR2(J2+2)
            RR1(J+2) = RR1(J+2)+X10*(RR2(J2-1)+ &
                                     RR2(J2+2))+ &
                                X11*(RR2(J2)+&
                                     RR2(J2+1))
            RR1(J+3) = RR1(J+3)+X03*RR2(J2-1)+ &
                                X02*RR2(J2)+ &
                                X01*RR2(J2+1)+ &
                                X00*RR2(J2+2)
         END DO
      endif

!     IN THE FOLLOWING SECTION THE ABSORPTION COEFFICIENT IS MULTIPIIED
!     BY THE RADIATION FIELD

      IF (JRAD.EQ.0) THEN

         XKT = TAVE/RADCN2
         VI = V1P-DV
         VITST = VI
         RDLAST = -1.
         NPTSI1 = NLO-1
         NPTSI2 = NLO-1

   30    NPTSI1 = NPTSI2+1

         VI = VFT+ REAL(NPTSI1-1)*DV
         RADVI = RADFNI(VI,DV,XKT,VITST,RDEL,RDLAST)

!         NPTSI2 = (VITST-VFT)/DV+1.001
!MJA 20150819 Implementing Yingtao Fix
         NPTSI2 = (VITST-VFT)/DV+0.001
         NPTSI2 = MIN(NPTSI2,NHI)

         tmpRADVI = RADVI
         DO 40 I = NPTSI1, NPTSI2
            R1(I) = R1(I)*tmpRADVI
            tmpRADVI = tmpRADVI+RDEL
   40    CONTINUE
         if ( NLTE ) then !NLTE case
            tmpRADVI = RADVI
            DO I = NPTSI1, NPTSI2
               rr1(i) = rr1(i)*tmpRADVI
               tmpRADVI = tmpRADVI+RDEL
            ENDDO
         endif

         IF (NPTSI2.LT.NHI) GO TO 30

      ENDIF

!     V1P IS FIRST FREQ OF PANEL
!     V2P IS LAST FREQ OF PANEL
      r1Buffer( numTotPoints+1:&
                numTotPoints+NLIM ) = R1(NLO:NHI)
      if ( NLTE ) then
         rr1Buffer( numTotPoints+1:&
                    numTotPoints+NLIM ) = RR1(NLO:NHI)
      endif
      numTotPoints = numTotPoints+NLIM


      VFT = VFT+ REAL(NLIM1-1)*DV
      IF (ISTOP.NE.1) THEN

         DO 50 J = NLIM1, MAX1
            R1(J-NLIM1+1) = R1(J)
   50    CONTINUE
         DO 60 J = MAX1-NLIM1+2, MAX1
            R1(J) = 0.
   60    CONTINUE
         DO 70 J = NLIM2, MAX2
            R2(J-NLIM2+1) = R2(J)
   70    CONTINUE
         DO 80 J = MAX2-NLIM2+2, MAX2
            R2(J) = 0.
   80    CONTINUE
         DO 90 J = NLIM3, MAX3
            R3(J-NLIM3+1) = R3(J)
   90    CONTINUE
         DO 100 J = MAX3-NLIM3+2, MAX3
            R3(J) = 0.
  100    CONTINUE

         if ( NLTE ) then
            DO J = NLIM1, MAX1
               RR1(J-NLIM1+1) = RR1(J)
            ENDDO
            DO J = MAX1-NLIM1+2, MAX1
               RR1(J) = 0.
            ENDDO
            DO J = NLIM2, MAX2
               RR2(J-NLIM2+1) = RR2(J)
            ENDDO
            DO J = MAX2-NLIM2+2, MAX2
               RR2(J) = 0.
            ENDDO
            DO J = NLIM3, MAX3
               RR3(J-NLIM3+1) = RR3(J)
            ENDDO
            DO J = MAX3-NLIM3+2, MAX3
               RR3(J) = 0.
            ENDDO
         endif

         NLO = NSHIFT+1
      ENDIF
      !CALL CPUTIM (TIME)
      !TIMPNL = TIMPNL+TIME-TIME0

      RETURN
      END SUBROUTINE






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!   Below are from OPRAP.f90   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!-----------------------------------------------------------------------
!
!     SUBROUTINE RDLIN INPUTS LINE DATA FROM FILE LINFIL
!
!-----------------------------------------------------------------------
      SUBROUTINE RDLIN( VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG, &
                        BRD_MOL_FLG,BRD_MOL_HW,BRD_MOL_TMP,BRD_MOL_SHFT, &
                        IOUT, LINFIL,NLNGTH,ILO,IHI,LIMIN,VBOT,VTOP, &
                        IDATA,IEOF,NOPR )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_LineData    ,ONLY: MXBRDMOL, NLINEREC, &
                                    INPUT_HEADER, INPUT_BLOCK
      USE Module_Config      ,ONLY: IPR

      IMPLICIT NONE !REAL*8           (V)

      real(r8)  ,intent(out)   :: VNU(:)
      real      ,intent(out)   :: SP(:)
      real      ,intent(out)   :: ALFA0(:)
      real      ,intent(out)   :: EPP(:)
      integer*4 ,intent(out)   :: MOL(:)
      real      ,intent(out)   :: HWHMS(:)
      real      ,intent(out)   :: TMPALF(:)
      real      ,intent(out)   :: PSHIFT(:)
      integer   ,intent(out)   :: IFLG(:)
      integer*4 ,intent(out)   :: BRD_MOL_FLG(:,:)  !(MXBRDMOL,NLINEREC)
      real      ,intent(out)   :: BRD_MOL_HW(:,:)   !(MXBRDMOL,NLINEREC)
      real      ,intent(out)   :: BRD_MOL_TMP(:,:)  !(MXBRDMOL,NLINEREC)
      real      ,intent(out)   :: BRD_MOL_SHFT(:,:) !(MXBRDMOL,NLINEREC)
      integer   ,intent(out)   :: IOUT(:) !=NLINEREC=250
      integer   ,intent(in)    :: LINFIL
      integer   ,intent(in)    :: NLNGTH
      integer   ,intent(inout) :: ILO,IHI
      integer   ,intent(in)    :: LIMIN
      real(r8)  ,intent(in)    :: VBOT, VTOP
      integer   ,intent(out)   :: IDATA
      integer   ,intent(inout) :: IEOF
      integer   ,intent(in)    :: NOPR


      !--- Local variables
      !
      character*8 :: HLINID(10)
      real        :: SDEP(NLINEREC)
      integer*4   :: N_NEGEPP(64)
      integer*4   :: N_RESETEPP(64)
      real*4      :: XSPACE(4096)

      integer*4 ,PARAMETER :: i_1=1
      integer*4 ,PARAMETER :: npnlhd=6

      type(input_header) :: rdlnpnl
      type(input_block)  :: rdlnbuf, dumbuf

      integer*4   :: molcnt(64), mcntlc(64), mcntnl(64)
      integer*4   :: linmol, lincnt, ilinlc, ilinnl, irec, irectl
      real*4      :: sumstr(64), flinlo, flinhi
      integer*4   :: nrec, nwds, lnfl,leof
      real*4      :: str, hw_f, e_low, xmol, hw_s, hw_T, shft
      integer*4   :: mol_id, jflg
      CHARACTER*1 :: CNEGEPP(8)

      INTEGER  :: M, I, J, J1, IJ, ILNGTH
      REAL     :: AMOL(NLINEREC)
      REAL     :: RVMR
      REAL(r8) :: VMIN




      lnfl = linfil

!     THERE ARE (LIMIN * 9) QUANTITIES READ IN:
!     VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG

      ILNGTH = NLNGTH*LIMIN
      IDATA = 0

!     BUFFER PAST FILE HEADER if necessary

      IF (ILO.LE.0) THEN
         REWIND LNFL
         read (lnfl) HLINID
         READ (HLINID(7),950) CNEGEPP
         IF (CNEGEPP(8).eq.'^') THEN
            read (lnfl) n_negepp,n_resetepp,xspace
         endif
      ENDIF

      !print *,'will read pnlhd'
   10 CALL BUFIN_sgl(Lnfl,LEOF,rdlnpnl,npnlhd)
      IF (LEOF.EQ.0) THEN
         IF (NOPR.EQ.0) WRITE (IPR,900)
         IEOF = 1
         RETURN
      ENDIF

      IF (rdlnpnl%nrec.GT.LIMIN) STOP 'RDLIN; NREC GT LIMIN'

      !print *,'will read to advance'
      IF (rdlnpnl%VMAX.LT.VBOT) THEN
         CALL BUFIN_sgl(Lnfl,LEOF,dumbuf,i_1)
         GO TO 10
      ENDIF

      !print *,'will read panel '
      CALL BUFIN_sgl(Lnfl,LEOF,rdlnbuf,rdlnpnl%NWDS)

      ! precision conversion occurs here:
      ! incoming on right: vlin is real*8, others real*4 and integer*4

      do 15 i=1,rdlnpnl%nrec

         IFLG(i) = rdlnbuf%iflg(i)
         VNU(i) = rdlnbuf%vnu(i)
         SP(i) = rdlnbuf%sp(i)
         ALFA0(i) = rdlnbuf%alfa(i)
         EPP(i) = rdlnbuf%epp(i)

         !yma if (iflg(i) .ge. 0) then
         !yma    MOL(i) = rdlnbuf%mol(i) ! int*4 to int*8 (if compiled as r8)
         !yma else
         !yma   xmol = transfer(rdlnbuf%mol(i),xmol) !int*4 to real*4)
         !yma   amol(i) = xmol   ! real*4 to real*8 (if compiled as r8)
         !yma   mol(i) = transfer (amol(i), mol(i))  ! real*8 to int*8  (if compiled as r8)
         !yma endif
         MOL(i) = rdlnbuf%mol(i) !The above transfer has problem when compile with i4 and r8 options. So force MOL to be integer*4

         HWHMS(i) = rdlnbuf%hwhm(i)
         TMPALF(i)= rdlnbuf%tmpalf(i)
         PSHIFT(i)= rdlnbuf%pshift(i)

         do j=1,mxbrdmol
            brd_mol_flg(j,i)=rdlnbuf%brd_mol_flg_in(j,i)
         end do
         j = 1
         do j1 = 1,mxbrdmol
            brd_mol_hw(j1,i) = rdlnbuf%brd_mol_dat(j,i)
            brd_mol_tmp(j1,i) = rdlnbuf%brd_mol_dat(j+1,i)
            brd_mol_shft(j1,i) = rdlnbuf%brd_mol_dat(j+2,i)
            j = j+3
         end do
         sdep(i) = rdlnbuf%speed_dep(i)


         !  MJA 20140909
         !  HITRAN provides widths for broadening by air; LBLRTM and MONORTM have always treated these widths as foreign
         !  This assumption is valid for most species, but not for N2 or O2. We now adjust the HITRAN widths to obtain
         !  true foreign widths. Similar ajdustment is applied if self shift information is available.
!v128         M = MOD(MOL(I),100)
!v128         !WRITE(*,*) M
!v128         if (M.eq.7 .AND. IFLG(i).ge.0) then
!v128            !WRITE(*,*) M, ALFA0(i),HWHMS(i)
!v128            rvmr = 0.21
!v128            ALFA0(i) = ( ALFA0(i)-rvmr*HWHMS(i))/(1.0-rvmr)
!v128            !WRITE(*,*) M, ALFA0(i),HWHMS(i)
!v128         endif
!v128         if (M.eq.22 .AND. IFLG(i).ge.0) then
!v128             rvmr = 0.79
!v128             ALFA0(i) = ( ALFA0(i)-rvmr*HWHMS(i))/(1.0-rvmr)
!v128         endif
         M = MOD(MOL(I),100)
         !WRITE(*,*) M
         if (M.eq.7 .AND. IFLG(i).ge.0) then
            !WRITE(*,*) M, ALFA0(i),HWHMS(i)
            rvmr = 0.21
            ALFA0(i) = ( ALFA0(i)-rvmr*HWHMS(i))/(1.0-rvmr)
            if (brd_mol_flg(m,i).gt.0) then
               pshift(i) = (pshift(i)-rvmr*brd_mol_shft(m,i))/(1.0-rvmr)
            endif
            !WRITE(*,*) M, ALFA0(i),HWHMS(i)
         endif
         if (M.eq.22 .AND. IFLG(i).ge.0) then
             rvmr = 0.79
             ALFA0(i) = ( ALFA0(i)-rvmr*HWHMS(i))/(1.0-rvmr)
             ! Currently SBS broadening is only code for the first seven HITRAN species
             ! When it becomes available for N2, the next three lines should become executable
             !if (brd_mol_flg(m,i).gt.0) then
             !   pshift(i) = (pshift(i)-rvmr*brd_mol_shft(m,i))/(1.0-rvmr)
             !endif
         endif


   15 continue

!yma      IF ((ILO.EQ.0).AND.(VMIN.GT.VBOT)) WRITE (IPR,905)
      IF ((ILO.EQ.0).AND.(rdlnpnl%VMIN.GT.VBOT)) WRITE (IPR,905)
      ILO = 1

      IJ = 0
      DO 20 I = 1, rdlnpnl%NREC
         IF (IFLG(I).GE.0) THEN
            IJ = IJ+1
            IOUT(IJ) = I
         ENDIF
   20 END DO

      DO 30 I = IJ+1, 250
         IOUT(I) = rdlnpnl%NREC
   30 END DO

!yma      IF (VMIN.LT.VBOT) THEN
      IF (rdlnpnl%VMIN.LT.VBOT) THEN
         DO 40 J = 1, IJ
            I = IOUT(J)
            ILO = J
            IF (VNU(I).GE.VBOT) GO TO 50
   40    CONTINUE
      ENDIF

   50 CONTINUE
      DO 60 J = ILO, IJ
         I = IOUT(J)
         IF (MOL(I).GT.0) THEN
            IHI = J
            IF (VNU(I).GT.VTOP) GO TO 70
         ENDIF
   60 END DO

      ! the following test is to see if more data is required
      ! idata = 1 means data requirements have been met

   70 IF (IHI.LT.IJ) IDATA = 1

      RETURN

  900 FORMAT ('  EOF ON LINFIL (MORE LINES MAY BE REQUIRED) ')
  905 FORMAT (                                                          &
     &   ' FIRST LINE ON LINFIL USED (MORE LINES MAY BE REQUIRED) ')
  950 FORMAT (8a1)

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE SHAPEL CONSTRUCTS THE SUB-FUNCTIONS FOR THE
!     LORENTZ LINE SHAPE
!-----------------------------------------------------------------------
      SUBROUTINE SHAPEL( F1,F2,F3, &
                         HWF1,HWF2,HWF3,DXF1,DXF2,DXF3,NX1,NX2,NX3,N1MAX,N2MAX,N3MAX )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: F1(:), F2(:), F3(:)
      real    ,intent(in)  :: HWF1,  HWF2,  HWF3
      real    ,intent(in)  :: DXF1,  DXF2,  DXF3
      integer ,intent(in)  :: NX1,   NX2,   NX3
      integer ,intent(in)  :: N1MAX, N2MAX, N3MAX

      !--- Local variables
      INTEGER :: I,      J1LIM, J1LIMP, J2LIM
      INTEGER :: J2LIMP, JJ
      REAL    :: A1,     A2,    A3,     B1
      REAL    :: B2,     B3,    RECPI,  SUM
      REAL    :: TOTAL,  X,     XSQ,    Z0

      !--- Statement Functions
      !
      REAL :: A, B
      REAL :: XLORNZ, Q1FN, Q2FN, Q3FN
      !
      XLORNZ(XSQ) = 1./(1.+XSQ)
      Q1FN(XSQ) = A1+B1*XSQ
      Q2FN(XSQ) = A2+B2*XSQ
      Q3FN(XSQ) = A3+B3*XSQ
      !
      A(Z0) = (1.+2.*Z0*Z0)/(1.+Z0*Z0)**2
      B(Z0) = -1./(1.+Z0*Z0)**2


      RECPI = 1./(2.*ASIN(1.))
      TOTAL = 0.

      A1 = A(HWF1)
      B1 = B(HWF1)

      A2 = A(HWF2)
      B2 = B(HWF2)

      A3 = A(HWF3)
      B3 = B(HWF3)

      DO 10 I = 1, N1MAX
         F1(I) = 0.
   10 END DO
      F1(1) = RECPI*(XLORNZ(0.)-Q1FN(0.))
      SUM = F1(1)
      DO 20 JJ = 2, NX1
         X = REAL(JJ-1)*DXF1
         XSQ = X*X
         F1(JJ) = RECPI*(XLORNZ(XSQ)-Q1FN(XSQ))
         SUM = SUM+F1(JJ)*2.
   20 END DO
      F1(NX1) = 0.
      SUM = SUM*DXF1
      TOTAL = TOTAL+SUM

      DO 30 I = 1, N2MAX
         F2(I) = 0.
   30 END DO
      F2(1) = RECPI*(Q1FN(0.)-Q2FN(0.))
      SUM = F2(1)
      J1LIM = HWF1/DXF2+1.001
      DO 40 JJ = 2, J1LIM
         X = REAL(JJ-1)*DXF2
         XSQ = X*X
         F2(JJ) = RECPI*(Q1FN(XSQ)-Q2FN(XSQ))
         SUM = SUM+F2(JJ)*2.
   40 END DO
      J1LIMP = J1LIM+1
      DO 50 JJ = J1LIMP, NX2
         X = REAL(JJ-1)*DXF2
         XSQ = X*X
         F2(JJ) = RECPI*(XLORNZ(XSQ)-Q2FN(XSQ))
         SUM = SUM+F2(JJ)*2.
   50 END DO
      F2(NX2) = 0.
      SUM = SUM*DXF2
      TOTAL = TOTAL+SUM

      DO 60 I = 1, N3MAX
         F3(I) = 0.
   60 END DO
      F3(1) = RECPI*(Q2FN(0.)-Q3FN(0.))
      SUM = F3(1)
      J2LIM = HWF2/DXF3+1.001
      DO 70 JJ = 2, J2LIM
         X = REAL(JJ-1)*DXF3
         XSQ = X*X
         F3(JJ) = RECPI*(Q2FN(XSQ)-Q3FN(XSQ))
         SUM = SUM+F3(JJ)*2.
   70 END DO
      J2LIMP = J2LIM+1
      DO 80 JJ = J2LIMP, NX3
         X = REAL(JJ-1)*DXF3
         XSQ = X*X
         F3(JJ) = RECPI*(XLORNZ(XSQ)-Q3FN(XSQ))
         SUM = SUM+F3(JJ)*2.
   80 END DO
      SUM = SUM*DXF3
      TOTAL = TOTAL+SUM

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE SHAPEG CONSTRUCTS THE FUNCTION FOR THE DOPPLER PROFILE
!-----------------------------------------------------------------------
      SUBROUTINE SHAPEG( FG,DXF1,NX1,N1MAX )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: FG(:)
      real    ,intent(in)  :: DXF1
      integer ,intent(in)  :: NX1
      integer ,intent(in)  :: N1MAX


      !--- Local variabls
      INTEGER :: I,      JJ
      REAL    :: FGNORM, FLN2,  RECPI,  SUM
      REAL    :: TOTAL,  X,     XSQ

      real :: FGAUSS
      FGAUSS(XSQ) = EXP(-FLN2*XSQ)


      FLN2 =  LOG(2.)
      RECPI = 1./(2.*ASIN(1.))
      FGNORM = SQRT(FLN2*RECPI)
      TOTAL = 0.
      DO 10 I = 1, N1MAX
         FG(I) = 0.
   10 END DO
      FG(1) = FGNORM*FGAUSS(0.)
      SUM = FG(1)
      DO 20 JJ = 2, NX1
         X = REAL(JJ-1)*DXF1
         XSQ = X*X
         FG(JJ) = FGNORM*FGAUSS(XSQ)
         SUM = SUM+FG(JJ)*2.
   20 END DO
      FG(NX1) = 0.
      SUM = SUM*DXF1
      TOTAL = TOTAL+SUM

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     VERFN IS A FUNCTION USED TO IMPROVE THE ACCURACY OF THE
!     VOIGT APPROXIMATION IN THE DOMANE 0 - 4 HALFWIDTHS.
!-----------------------------------------------------------------------
      SUBROUTINE VERFN( XVER,DXF1,N1MAX )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: pi
      IMPLICIT NONE

      real     ,intent(out) :: XVER(:)
      real     ,intent(in)  :: DXF1
      integer  ,intent(in)  :: N1MAX

      !--- Local variables
      INTEGER :: I
      REAL    :: AE0,  AE2,   AE4,  CE0
      REAL    :: CE2,  CE4,   CEXP, FACTOR
      REAL    :: SE0,   SUM0, SUM2
      REAL    :: SUM4, SUMER, XE0,  XE2
      REAL    :: XE4,  Z,     Z2

!     FOR ZETA = 0.3
      DATA CEXP,CE0,CE2,CE4 / 0.45,1.,-.20737285249,-.00872684335747 /
      DATA SUM0,SUM2,SUM4,SUMER / 4*0. /

      !--- Statement Function
      real :: ERFN
      ERFN(Z2) = (1./(CE0+CE2+CE4))*(CE0+CE2*AE2*Z2+CE4*AE4*Z2*Z2)*XE0


      IF (SUMER.NE.0.) RETURN
      SE0 = SQRT(CEXP/PI)
      AE0 = 1.
      AE2 = 2.*CEXP
      AE4 = AE2*AE2/3.
      FACTOR = 1.

      DO 10 I = 1, N1MAX
         XVER(I) = 0.
   10 END DO

      DO 20 I = 1, N1MAX
         Z = DXF1* REAL(I-1)
         Z2 = Z*Z
         XE0 = SE0*EXP(-CEXP*Z2)
         XE2 = AE2*Z2*XE0
         XE4 = AE4*Z2*Z2*XE0
         XVER(I) = ERFN(Z2)
         SUM0 = SUM0+FACTOR*DXF1*XE0
         SUM2 = SUM2+FACTOR*DXF1*XE2
         SUM4 = SUM4+FACTOR*DXF1*XE4
         SUMER = SUMER+FACTOR*DXF1*XVER(I)
         FACTOR = 2.
   20 END DO

!PRT  WRITE (IPR,900) Z,SUM0,SUM2,SUM4,SUMER

      RETURN

!  900 FORMAT (F10.3,6F15.10)

      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE RSYM( R,lb,DV,VFT )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE !REAL*8           (V)

      real     ,intent(inout) :: R(lb:)
      integer  ,intent(in)    :: lb    !The lower bound of R(:)
      real     ,intent(in)    :: DV
      real(r8) ,intent(in)    :: VFT

      !--- Local variables
      INTEGER :: I,  IP,  IPMAX, K
      REAL    :: B1, B2,  C1,    C2
      REAL    :: P,  PST, W0,    W1
      REAL    :: W2, WN1


      IP = (-VFT/DV)+1.-.000001
      IP = IP+1
      P = ( REAL(IP-1)+VFT/DV)*2.
      PST = P
      IF (P.GT.1.) P = P-1.

!     VFT/DV- INT(VFT/DV)= 0. TO 0.5

      WN1 = -P*(P-1.)*(P-2.)/6.
      W0 = (P*P-1.)*(P-2.)/2.
      W1 = -P*(P+1.)*(P-2.)/2.
      W2 = P*(P*P-1.)/6.
      K = IP
      IPMAX = IP+IP-1
      IF (PST.LE.1.) GO TO 20
      B1 = R(IP-2)
      B2 = R(IP-1)
      C1 = R(K)
      DO 10 I = IP, IPMAX
         K = K-1
         C2 = C1
         IF (K.LT.1) GO TO 40
         C1 = R(K)
         R(K) = R(K)+WN1*R(I+1)+W0*R(I)+W1*B2+W2*B1
         B1 = B2
         B2 = R(I)
         IF (K.LE.2) GO TO 10
         R(I) = R(I)+WN1*C2+W0*C1+W1*R(K-1)+W2*R(K-2)
   10 END DO
      GO TO 40

!    VFT/DV- INT(VFT/DV) = 0.5 TO 1.0

   20 C1 = R(IP)
      C2 = R(IP+1)
      B2 = R(IP-1)
      DO 30 I = IP, IPMAX
         K = K-1
         B1 = B2
         B2 = R(I)
         IF (K.LE.1) GO TO 40
         R(I) = R(I)+WN1*C2+W0*C1+W1*R(K)+W2*R(K-1)
         C2 = C1
         C1 = R(K)
         R(K) = R(K)+WN1*R(I+2)+W0*R(I+1)+W1*B2+W2*B1
   30 END DO

   40 RETURN
      END SUBROUTINE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
     SUBROUTINE nonLTE_RatState( RATSTATE, NUMSTATE, ALTAV, TAVE )
!-----------------------------------------------------------------------
      USE Module_Config      ,ONLY: IPR, ioFiles
      USE Module_ConstParam  ,ONLY: CMOL, MAXSTATE, Max_ISO, MXMOL

      IMPLICIT NONE !REAL*8           (V)

      real     ,intent(out) :: RATSTATE(:,:)!(MAXSTATE*MAX_ISO,MXMOL)
      integer  ,intent(out) :: NUMSTATE(:)  !(MXMOL)
      real     ,intent(in)  :: ALTAV
      real     ,intent(in)  :: TAVE


      !--- Local variables
      !
      integer      :: NLTEFL
      character*5  :: IDSTATE(MAXSTATE,MXMOL)
      real         :: EESTATE(MAXSTATE*MAX_ISO,MXMOL)
      integer      :: NDGSTATE(MAXSTATE*MAX_ISO,MXMOL)
      logical      :: ISODATA
      character*80 :: TEXTLINE
      integer      :: ISORDER(MAXSTATE*Max_ISO)
      character*6  :: TXTISO

      INTEGER :: I, ID, IVIB
      INTEGER :: MOLNEQ
      REAL    :: XKT




      CALL DEFNLTEDAT( NUMSTATE,IDSTATE,EESTATE,NDGSTATE,RATSTATE )

      ISODATA=.FALSE.

      inquire( file=ioFiles%nlteStatPopFile, number=NLTEFL )
      REWIND NLTEFL

      ! READ UP TO 20 LINES OF TEXT AT BEGINNING OF TAPE4
      WRITE(IPR,890)
  890 FORMAT(/2X,'TAPE4 HEADER:')
                   ! 20 MAX TEXT LINES AT BEGINNING OF TAPE4
      DO I=1,20
         READ(NLTEFL,900) TEXTLINE
  900    FORMAT(A80)
         IF(TEXTLINE(1:1).NE.'!') GO TO 915
         WRITE(IPR,910) TEXTLINE
  910    FORMAT(2X,A80)
      END DO

  915 READ(TEXTLINE,920) IVIB,MOLNEQ
  920 FORMAT(2I5)
      WRITE(IPR,921) IVIB,ALTAV,TAVE
  921 FORMAT(/' IVIB =',I5,/'  ALT = ',F10.3,'  TEMP =',F10.3)

      READ(NLTEFL,900) TEXTLINE
      write(ipr,940) textline
  940 FORMAT(A80)

      DO I=1,74
         IF(TEXTLINE(I:I+6).EQ.'VIBRATI') ISODATA=.TRUE.
      END DO

      IF(ISODATA) THEN
   10    CALL GETINDEX(TEXTLINE,CMOL,MXMOL,ID,TXTISO)
!           END OF DATA ENCOUNTERED
            IF(ID.EQ.0) GO TO 30

            !yma CALL RDNLTE( NLTEFL,TEXTLINE,TXTISO,NUMSTATE(ID), IDSTATE(1,ID),&
            !yma              EESTATE(1,ID),NDGSTATE(1,ID), ISORDER)
            CALL RDNLTE( NLTEFL,TEXTLINE,TXTISO,NUMSTATE(ID), IDSTATE(:,ID),&
                         EESTATE(:,ID),NDGSTATE(:,ID), ISORDER)

            IF(IVIB.EQ.1) THEN
               !yma CALL VIBTMP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),   &
               !yma              NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,        &
               !yma              TEXTLINE,ISORDER, TAVE )
               CALL VIBTMP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(:,ID),   &
                            NDGSTATE(:,ID),EESTATE(:,ID), RATSTATE(:,ID),TXTISO,        &
                            TEXTLINE,ISORDER, TAVE )
            ELSE
               !yma CALL VIBPOP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),   &
               !yma              NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,        &
               !yma              TEXTLINE,ISORDER, TAVE )
               CALL VIBPOP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(:,ID),   &
                            NDGSTATE(:,ID),EESTATE(:,ID), RATSTATE(:,ID),TXTISO,        &
                            TEXTLINE,ISORDER, TAVE )
            END IF
!                    IF END OF DATA ENCOUNTERED
            DO I=1,70
               IF(TEXTLINE(I:I+10).EQ.'END OF DATA') GO TO 30
            END DO
!           READ DATA FOR NEXT SPECIE
         GO TO 10
      ELSE
         DO ID=1,MXMOL
         IF(NUMSTATE(ID).GT.0) THEN
            TXTISO = adjustl(CMOL(ID)) !yma CALL DROPSPACE(CMOL(ID),TXTISO)

            IF(IVIB.EQ.1) THEN
               !yma CALL VIBTMP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),&
               !yma              NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,     &
               !yma              TEXTLINE,ISORDER, TAVE )
               CALL VIBTMP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(:,ID),&
                            NDGSTATE(:,ID),EESTATE(:,ID), RATSTATE(:,ID),TXTISO,     &
                            TEXTLINE,ISORDER, TAVE )
            ELSE
               !yma CALL VIBPOP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),&
               !yma              NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,     &
               !yma              TEXTLINE,ISORDER, TAVE )
               CALL VIBPOP( XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(:,ID),&
                            NDGSTATE(:,ID),EESTATE(:,ID), RATSTATE(:,ID),TXTISO,     &
                            TEXTLINE,ISORDER, TAVE )
            END IF
         END IF
         END DO
      END IF

   30 continue

      WRITE(IPR,930)
  930 FORMAT(//)

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE DEFNLTEDAT( NUMSTATE,IDSTATE,EESTATE,NDGSTATE,RATSTATE )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: CMOL, MAXSTATE, MAX_ISO, MXMOL
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      integer     ,intent(out) :: NUMSTATE(:)   !(MXMOL)
      character*5 ,intent(out) :: IDSTATE(:,:)  !(MAXSTATE,MXMOL)
      real        ,intent(out) :: EESTATE(:,:)  !(MAXSTATE*MAX_ISO,MXMOL)
      integer     ,intent(out) :: NDGSTATE(:,:) !(MAXSTATE*MAX_ISO,MXMOL)
      real        ,intent(out) :: RATSTATE(:,:) !(MAXSTATE*MAX_ISO,MXMOL)


      INTEGER      :: I,     ID,    INDEX,   ISO
      INTEGER      :: ISO1,  ISO2,  ISOTOPE, J
      INTEGER      :: J0
      INTEGER      :: K
      integer      :: NUMISO(MXMOL)
      CHARACTER*6  :: TXTISO

      !gFortran CHARACTER*80 ,PARAMETER :: ISOMOL(5) = ['H2O','CO2','O3','CO','NO']
      character(80),SAVE :: ISOMOL(5)
      data ISOMOL / 'H2O', 'CO2', 'O3', 'CO', 'NO' /
      integer      ,PARAMETER :: MAXH2O=8, MAXCO2=26, MAXO3=18, MAXCO=3, MAXNO=3

      integer     :: JSOH2O(MAXH2O), JSOCO2(MAXCO2), JSOO3(MAXO3), JSOCO(MAXCO), JSONO(MAXNO)
      character*5 ::   AH2O(MAXH2O),   ACO2(MAXCO2),   AO3(MAXO3),   ACO(MAXCO),   ANO(MAXNO)
      real        ::   FH2O(MAXH2O),   FCO2(MAXCO2),   FO3(MAXO3),   FCO(MAXCO),   FNO(MAXNO)
      integer     :: MDGH2O(MAXH2O), MDGCO2(MAXCO2), MDGO3(MAXO3), MDGCO(MAXCO), MDGNO(MAXNO)


      DATA  (JSOH2O(I),AH2O(I),FH2O(I),MDGH2O(I),I=1,8)/                &
     &     1, '000' ,     0.   , 1,                                     &
     &     1, '010' ,  1594.750, 1,                                     &
     &     1, '020' ,  3151.630, 1,                                     &
     &     1, '100' ,  3657.053, 1,                                     &
     &     1, '001' ,  3755.930, 1,                                     &
     &     1, '030' ,  4666.793, 1,                                     &
     &     1, '110' ,  5234.977, 1,                                     &
     &     1, '011' ,  5333.269, 1/

      DATA  (JSOCO2(I),ACO2(I),FCO2(I),MDGCO2(I),I=1, 9)/               &
     &     1, '00001' ,    0.   , 1 ,                                   &
     &     1, '01101' ,  667.380, 2 ,                                   &
     &     1, '10002' , 1285.409, 1 ,                                   &
     &     1, '02201' , 1335.132, 2 ,                                   &
     &     1, '10001' , 1388.185, 1 ,                                   &
     &     1, '11102' , 1932.470, 2 ,                                   &
     &     1, '03301' , 2003.246, 2 ,                                   &
     &     1, '11101' , 2076.856, 2 ,                                   &
     &     1, '00011' , 2349.143, 1 /
      DATA  (JSOCO2(I),ACO2(I),FCO2(I),MDGCO2(I),I=10,26)/              &
     &     1, '20003' , 2548.366, 1 ,                                   &
     &     1, '12202' , 2585.022, 2 ,                                   &
     &     1, '20002' , 2671.143, 1 ,                                   &
     &     1, '04401' , 2671.717, 2 ,                                   &
     &     1, '12201' , 2760.725, 2 ,                                   &
     &     1, '20001' , 2797.135, 1 ,                                   &
     &     1, '01111' , 3004.012, 2 ,                                   &
     &     1, '10012' , 3612.842, 1 ,                                   &
     &     1, '02211' , 3659.273, 2 ,                                   &
     &     1, '10011' , 3714.783, 1 ,                                   &
     &     1, '11112' , 4247.706, 2 ,                                   &
     &     1, '03311' , 4314.914, 2 ,                                   &
     &     1, '11111' , 4390.629, 2 ,                                   &
     &     1, '20013' , 4853.623, 1 ,                                   &
     &     1, '04411' , 4970.931, 1 ,                                   &
     &     1, '20012' , 4977.834, 1 ,                                   &
     &     1, '20011' , 5099.660, 1 /

      DATA (JSOO3(I),AO3(I),FO3(I),MDGO3(I),I=1,18)/                    &
     &     1, '000' ,    0.   , 1,                                      &
     &     1, '010' ,  700.931, 1,                                      &
     &     1, '001' , 1042.084, 1,                                      &
     &     1, '100' , 1103.140, 1,                                      &
     &     1, '020' , 1399.275, 1,                                      &
     &     1, '011' , 1726.528, 1,                                      &
     &     1, '110' , 1796.261, 1,                                      &
     &     1, '002' , 2057.892, 1,                                      &
     &     1, '101' , 2110.785, 1,                                      &
     &     1, '200' , 2201.157, 1,                                      &
     &     1, '111' , 2785.245, 1,                                      &
     &     1, '003' , 3041.200, 1,                                      &
     &     1, '004' , 3988.   , 1,                                      &
     &     1, '005' , 4910.   , 1,                                      &
     &     1, '006' , 5803.   , 1,                                      &
     &     1, '007' , 6665.   , 1,                                      &
     &     1, '008' , 7497.   , 1,                                      &
     &     1, '009' , 8299.   , 1/

      DATA (JSOCO(I),ACO(I),FCO(I),MDGCO(I),I=1,3)/                     &
     &     1, '0' ,    0.   , 1,                                        &
     &     1, '1' , 2143.272, 1,                                        &
     &     1, '2' , 4260.063, 1/

      DATA (JSONO(I),ANO(I),FNO(I),MDGNO(I),I=1,3)/                     &
     &     1, '0' ,    0.   , 1,                                        &
     &     1, '1' , 1878.077, 1,                                        &
     &     1, '2' , 3724.067, 1/

      DO J=1,MXMOL
         NUMSTATE(J)=0
         NUMISO(J)=0
         DO I=1,MAXSTATE
            IDSTATE(I,J)='     '
            DO ISO=1,Max_ISO
               INDEX=(ISO-1)*MAXSTATE + I
               EESTATE(INDEX,J)=0.
               NDGSTATE(INDEX,J)=0
               RATSTATE(INDEX,J)=1.
            END DO
         END DO
      END DO

!---H2O
      CALL GETINDEX(ISOMOL(1),CMOL,MXMOL,ID,TXTISO)
      DO I=1,MAXH2O
         ISOTOPE=JSOH2O(I)
         IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
     &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
         IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
         IF(ISOTOPE.EQ.1) THEN
            NUMSTATE(ID) = NUMSTATE(ID)+1
            INDEX=NUMSTATE(ID)
            IDSTATE(INDEX,ID)=AH2O(I)
            ISO1=1
            ISO2=Max_ISO
         ELSE
            INDEX=0
            DO J=1,NUMSTATE(ID)
               IF(AH2O(I).EQ.IDSTATE(J,ID)) INDEX=J
            END DO
            IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR H2O ISOTOPE>1'
            ISO1=ISOTOPE
            ISO2=ISOTOPE
         END IF
         DO J=ISO1,ISO2
            K=(J-1)*MAXSTATE + INDEX
            NDGSTATE(K,ID)=MDGH2O(I)
            EESTATE(K,ID)=FH2O(I)
         END DO
      END DO

!---CO2
      CALL GETINDEX(ISOMOL(2),CMOL,MXMOL,ID,TXTISO)
      DO I=1,MAXCO2
         ISOTOPE=JSOCO2(I)
         IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
     &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
         IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
         IF(ISOTOPE.EQ.1) THEN
            NUMSTATE(ID) = NUMSTATE(ID)+1
            INDEX=NUMSTATE(ID)
            IDSTATE(INDEX,ID)=ACO2(I)
            ISO1=1
            ISO2=Max_ISO
         ELSE
            INDEX=0
            DO J=1,NUMSTATE(ID)
            IF(ACO2(I).EQ.IDSTATE(J,ID)) INDEX=J
            END DO
            IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR CO2 ISOTOPE>1'
            ISO1=ISOTOPE
            ISO2=ISOTOPE
         END IF
         DO J=ISO1,ISO2
            K=(J-1)*MAXSTATE + INDEX
            NDGSTATE(K,ID)=MDGCO2(I)
            EESTATE(K,ID)=FCO2(I)
         END DO
      END DO

!---O3
      CALL GETINDEX(ISOMOL(3),CMOL,MXMOL,ID,TXTISO)
      DO I=1,MAXO3
         ISOTOPE=JSOO3(I)
         IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
     &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
         IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
         IF(ISOTOPE.EQ.1) THEN
            NUMSTATE(ID) = NUMSTATE(ID)+1
            INDEX=NUMSTATE(ID)
            IDSTATE(INDEX,ID)=AO3(I)
            ISO1=1
            ISO2=Max_ISO
         ELSE
            INDEX=0
            DO J=1,NUMSTATE(ID)
            IF(AO3(I).EQ.IDSTATE(J,ID)) INDEX=J
            END DO
            IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR O3 ISOTOPE>1'
            ISO1=ISOTOPE
            ISO2=ISOTOPE
         END IF
         DO J=ISO1,ISO2
            K=(J-1)*MAXSTATE + INDEX
            NDGSTATE(K,ID)=MDGO3(I)
            EESTATE(K,ID)=FO3(I)
         END DO
      END DO

!---CO
      CALL GETINDEX(ISOMOL(4),CMOL,MXMOL,ID,TXTISO)
      DO I=1,MAXCO
         ISOTOPE=JSOCO(I)
         IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
     &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
         IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
         IF(ISOTOPE.EQ.1) THEN
            NUMSTATE(ID) = NUMSTATE(ID)+1
            INDEX=NUMSTATE(ID)
            IDSTATE(INDEX,ID)=ACO(I)
            ISO1=1
            ISO2=Max_ISO
         ELSE
            INDEX=0
            DO J=1,NUMSTATE(ID)
            IF(ACO(I).EQ.IDSTATE(J,ID)) INDEX=J
            END DO
            IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR CO ISOTOPE>1'
            ISO1=ISOTOPE
            ISO2=ISOTOPE
         END IF
         DO J=ISO1,ISO2
            K=(J-1)*MAXSTATE + INDEX
            NDGSTATE(K,ID)=MDGCO(I)
            EESTATE(K,ID)=FCO(I)
         END DO
      END DO

!---NO
      CALL GETINDEX(ISOMOL(5),CMOL,MXMOL,ID,TXTISO)
      DO I=1,MAXNO
         ISOTOPE=JSONO(I)
         IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
     &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
         IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
         IF(ISOTOPE.EQ.1) THEN
            NUMSTATE(ID) = NUMSTATE(ID)+1
            INDEX=NUMSTATE(ID)
            IDSTATE(INDEX,ID)=ANO(I)
            ISO1=1
            ISO2=Max_ISO
         ELSE
            INDEX=0
            DO J=1,NUMSTATE(ID)
            IF(ANO(I).EQ.IDSTATE(J,ID)) INDEX=J
            END DO
            IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR NO ISOTOPE>1'
            ISO1=ISOTOPE
            ISO2=ISOTOPE
         END IF
         DO J=ISO1,ISO2
            K=(J-1)*MAXSTATE + INDEX
            NDGSTATE(K,ID)=MDGNO(I)
            EESTATE(K,ID)=FNO(I)
         END DO
      END DO

      DO I=1,MXMOL
         IF(NUMSTATE(I).GT.MAXSTATE) THEN
         WRITE(IPR,*) 'ERROR IN DEFNLTEDAT: MAXSTATE NEEDS TO BE '                       !WRITE(IPR,*) 'ERROR IN DEFNLTEDAT: MAXSTATE NEEDS TO BE '
         WRITE(IPR,*) 'INCREASED TO ',I,' TO ACCOMODATE DEFAULT ','ISOTOPE STATE DATA'   !WRITE(IPR,*) 'INCREASED TO ',I,' TO ACCOMODATE DEFAULT ','ISOTOPE STATE DATA'
         STOP 'ERROR IN DEFNLTEDAT'
         END IF
      END DO

      DO ID=1,MXMOL
         IF(NUMSTATE(ID).GT.0) THEN
            write(ipr,*) 'Defaults for Molecule',cmol(Id),'  id=',id,   &
     &           '  NUMSTATE= ',NUMSTATE(ID),'  NUMISO=',NUMISO(ID)
            write(ipr,901) (idstate(j,id),j=1,numstate(id))
  901       format('IDSTATE= ',26a6)
            do iso=1,NUMISO(ID)
               j0=(iso-1)*maxstate
               write(ipr,902) (ndgstate(j+j0,id),j=1,numstate(id))
  902          format('NDGSTATE= ',26I3)
            end do
            do iso=1,NUMISO(ID)
               j0=(iso-1)*maxstate
               write(ipr,903) (eestate(j+j0,id),j=1,numstate(id))
  903          FORMAT('EESTATE=',26f9.3)
            end do
         END IF
      end do

      RETURN
      END SUBROUTINE


!-----------------------------------------------------------------------
!   GET MOLECULE INDEX FROM CMOL ARRAY FOR VIBRATIONAL DATA
!-----------------------------------------------------------------------
      SUBROUTINE GETINDEX( TEXTLINE,CMOL,MXMOL,ID,TXTISO )
!-----------------------------------------------------------------------
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      character*80 ,intent(in)  :: TEXTLINE
      character*6  ,intent(in)  :: CMOL(:) !(MXMOL)
      integer      ,intent(in)  :: MXMOL
      integer      ,intent(out) :: ID
      character*6  ,intent(out) :: TXTISO

      CHARACTER*6 ,PARAMETER :: BLANKS='      '
      CHARACTER*6  :: TXTMOL
      INTEGER      :: I, I1, I2, ISP


      I1=0
      I2=0
      DO I=1,80
         IF(TEXTLINE(I:I).NE.'-' .AND. TEXTLINE(I:I).NE.' ') THEN
         IF(I1.EQ.0) I1=I
         END IF
         IF(I1.GT.0 .AND. TEXTLINE(I:I).EQ.' ') THEN
         IF(I2.EQ.0) I2=I-1
         END IF
      END DO
      TXTISO=TEXTLINE(I1:I2)//BLANKS
      ID=0
      DO ISP=1,MXMOL
         TXTMOL = adjustl(CMOL(ISP)) !yma CALL DROPSPACE(CMOL(ISP),TXTMOL)
         IF(TXTISO.EQ.TXTMOL) THEN
         ID=ISP
            WRITE(IPR,*) 'GETINDEX: ',TXTISO,'  INDEX=',ID
         RETURN
         END IF
      END DO
      WRITE(IPR,*) 'READING NONLTE DATA HEADER FROM TAPE4'
      WRITE(IPR,*) 'BUT THIS SPECIE ',TXTISO,' NOT FOUND IN MOLECULE LIST'
!##     STOP 'ERROR READING NLTE DATA FROM TAPE4'

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE RDNLTE( NLTEFL,TEXTLINE,TXTISO,IMAX,                    &
     &                   IDX,EEX,NDG,ISORDER )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: MAXSTATE, MAX_ISO
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      integer      ,intent(in)    :: NLTEFL
      character*80 ,intent(out)   :: TEXTLINE
      character*6  ,intent(in)    :: TXTISO
      integer      ,intent(in)    :: IMAX
      character*5  ,intent(in)    :: IDX(:)     !(MAXSTATE)
      real         ,intent(inout) :: EEX(:)     !(MAXSTATE*MAX_ISO)
      integer      ,intent(inout) :: NDG(:)     !(MAXSTATE*MAX_ISO)
      integer      ,intent(out)   :: ISORDER(:) !(MAXSTATE*MAX_ISO)


      !--- Local variables
      !
      CHARACTER*5 ,PARAMETER :: BLANKS=' '

      CHARACTER*1 :: QUOTE
      CHARACTER*5 :: IDENT
      INTEGER     :: I,        I1,      I2,       INDEX
      INTEGER     :: INDEX2,   INUM,    ISOTOPE
      INTEGER     :: NN



      ! QUOTE = '
      QUOTE=CHAR(39)

      !       INITIALIZE ARRAY TO HOLD ISOTOPE ORDER INFO
      DO I=1,MAXSTATE*Max_ISO
         ISORDER(I)=0
      END DO

      INUM=0
   50 READ(NLTEFL,940) TEXTLINE
  940 FORMAT(A80)

         IF(TEXTLINE(1:2).EQ.'--') RETURN

         INUM=INUM+1
         I1=0
         I2=0
         DO I=1,80
            IF(TEXTLINE(I:I).EQ.QUOTE) THEN
            IF(I1.EQ.0) THEN
               I1=I
            ELSE IF(I2.EQ.0) THEN
               I2=I
            END IF
            END IF
         END DO

         READ(TEXTLINE(1:I1-1),*) NN,ISOTOPE
         IF(NN.NE.INUM .AND. ISOTOPE.EQ.1) THEN
            WRITE(IPR,*) 'WARNING IN TAPE 4: ISOTOPE DATA FOR ',TXTISO
            WRITE(IPR,*) 'LINE ',INUM,' NOT IN ORDER'
!            STOP 'NN.NE.INUM IN NONLTE DATA FROM TAPE4'
         END IF
         IF(ISOTOPE.EQ.1 .AND. INUM.GT.IMAX) THEN
            WRITE(IPR,*) 'READING NONLTE DATA FROM TAPE4'
            WRITE(IPR,*) 'ISOTOPE DATA FOR ',TXTISO, &
     &           ' CONTAINS MORE STATES THAN DEFAULT ALLOWANCE'
            STOP 'NUMBER OF NONLTE STATES IN TAPE4 EXCEEDS DEFAULT'
         END IF
         IF(ISOTOPE.GT.Max_ISO) THEN
            WRITE(IPR,*) 'ERROR IN TAPE4:  ISOTOPE NUMBER ',ISOTOPE,&
     &           ' GREATER THAN Max_ISO=',Max_ISO
            STOP 'TAPE4 ERROR: ISOTOPE NUMBER TOO LARGE'
         END IF

         IDENT=TEXTLINE(I1+1:I2-1)//BLANKS
         INDEX=0
         DO I=1,MAXSTATE
            IF(IDENT.EQ.IDX(I)) INDEX=I
         END DO
         IF(INDEX.EQ.0) THEN
            WRITE(IPR,*) 'ERROR IN TAPE 4 ISOTOPE DATA FOR ',TXTISO, &
            ' ON LINE ',INUM
            WRITE(IPR,*) '   ISOTOPE NUMBER ',ISOTOPE,' IDENTIFIER ',IDENT,&
     &          ' NOT CONSISTENT WITH EXPECTED INPUT'
            STOP 'ERROR IN ISOTOPE INPUT IN NONLTE DATA FROM TAPE4'
         END IF

         INDEX2=(ISOTOPE-1)*MAXSTATE + INDEX
         ISORDER(INUM)=INDEX2

         READ(TEXTLINE(I2+1:80),*) EEX(INDEX2),NDG(INDEX2)
         IF(INUM.EQ.1 .AND. ABS(EEX(INDEX2)).GE.1.E-25) THEN
            WRITE(IPR,*) 'RDNLTE: GROUND STATE NOT SPECIFIED ', &
     &           'FOR MOLECULE ',TXTISO
            STOP 'TAPE4 GROUND STATE MISSING'
         END IF
!         WRITE(IPR,*) 'RDNLTE, LINE ',inum,' state=',index,'  index=',  &
!           index2,'  isorder=',isorder(inum)

      GO TO 50

      END SUBROUTINE


!-----------------------------------------------------------------------
!     SUBROUTINE VIBPOP USES THE NON-LTE POPULATION DATA FROM
!     TAPE4 TO CALCULATE THE VIBRATIONAL POPULATION ENHANCEMENT
!     RATIOS FOR SELECTED VIBRATIONAL STATES OF H2O,CO2,NO AND O3.
!-----------------------------------------------------------------------
      SUBROUTINE VIBPOP( XKT,HT,NLTEFLAG,NUM,IDX,NDEG,EH,RAT,  &
                         TITMOL,TEXTLINE,ISORDER, TAVE )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, RADCN2, MAXSTATE, MAX_ISO
      USE Module_Utility     ,ONLY: EXPINT_real8
      USE Module_Config      ,ONLY: IPR

      IMPLICIT NONE !REAL*8           (V)

      real         ,intent(out)   :: XKT
      real         ,intent(in)    :: HT
      integer      ,intent(in)    :: NLTEFLAG
      integer      ,intent(in)    :: NUM                    !for information
      character*5  ,intent(in)    :: IDX(:)     !(MAXSTATE) !for information
      integer      ,intent(in)    :: NDEG(:)    !(MAXSTATE*MAX_ISO)
      real         ,intent(in)    :: EH(:)      !(MAXSTATE*MAX_ISO)
      real         ,intent(out)   :: RAT(:)     !(MAXSTATE*MAX_ISO)
      character*6  ,intent(in)    :: TITMOL
      character*80 ,intent(inout) :: TEXTLINE
      integer      ,intent(in)    :: ISORDER(:) !MAXSTATE*MAX_ISO)
      !
      real         ,intent(in)    :: TAVE


      !--- Internal Variables
      !
      character*6 ,PARAMETER :: BLANKS=' '

      character*6 :: HMNLTE
      real(r8)    :: VQNE(   MAXSTATE*Max_ISO), &
                     VQEQ(   MAXSTATE*Max_ISO), &
                     VPNE1(  MAXSTATE*Max_ISO), &
                     VPNE2(  MAXSTATE*Max_ISO), &
                     VQNEST( MAXSTATE*Max_ISO), &
                     VQNEIN( MAXSTATE*Max_ISO)
      real        :: TNE(    MAXSTATE*Max_ISO)

      INTEGER  :: I,        I1,      I2,     INDEX
      INTEGER  :: ISOTOPE,  J,       K,      LVL
      INTEGER  :: NUMIN
      REAL     :: A,        ALT1,    ALT2,   DEN
      REAL     :: POPEQ,    POPNE
      REAL     :: TOTMOL
      REAL(r8) :: VST1




!   READ NLTE VIB POPULATIONS

      XKT= TAVE/RADCN2
      DO I=1,MAXSTATE*Max_ISO
         IF(ISORDER(I).GT.0) NUMIN=I
      end do
!      write(ipr,*) 'called vibpop with txtmol=',titmol,'  alt=',ht
!      write(ipr,*) 'num=',num,'  numin=',numin
!      write(ipr,*) (isorder(i),i=1,numin)

      I1=0
      I2=0
      DO I=1,80
         IF(TEXTLINE(I:I).NE.'-' .AND. TEXTLINE(I:I).NE.' ') THEN
            IF(I1.EQ.0) I1=I
         END IF
         IF(I1.GT.0 .AND. TEXTLINE(I:I).EQ.' ') THEN
            IF(I2.EQ.0) I2=I-1
         END IF
      END DO
      HMNLTE=TEXTLINE(I1:I2)//BLANKS
      IF(HMNLTE.NE.TITMOL) THEN
         READ (NLTEFLAG,902) HMNLTE
         i1=0
         DO I=1,6
         IF(HMNLTE(I:I).EQ.' ') I1=I
         END DO
         HMNLTE=HMNLTE(I1+1:6)//BLANKS
      END IF
      IF(HMNLTE.NE.TITMOL) THEN
         WRITE(IPR,*) 'READING VIBPOP DATA FROM TAPE4'
         WRITE(IPR,*) 'EXPECTED PROFILE DATA FOR ',TOTMOL
         WRITE(IPR,*) 'BUT READ SPECIE ',HMNLTE
         STOP 'ERROR READING NLTE DATA FROM TAPE4'
      END IF
      READ (NLTEFLAG,*)  ALT1,(VPNE1(I),I=1,NUMIN)
   10 READ (NLTEFLAG,*)  ALT2,(VPNE2(I),I=1,NUMIN)

!*****WOG, 11/06/2000: ALT1 -> AL2:
      !     IF( ALT1.LT.HT) THEN
      IF( ALT2.LE.HT) THEN
         ALT1 = ALT2
         DO 20 I=1,NUMIN
            VPNE1(I) = VPNE2(I)
   20    CONTINUE
         GO TO 10
      ENDIF
      !     CALL LININT(HT,ALT1,ALT2,NUMIN,VPNE1,VPNE2,VQNEIN)
      A = (HT-ALT1)/(ALT2-ALT1)
      DO 30 I=1,NUMIN
         !CALL EXPINT(VQNEIN(I),VPNE1(I),VPNE2(I),A )
         CALL EXPINT_real8(VQNEIN(I),VPNE1(I),VPNE2(I),A )
   30 END DO
      DO I=1,NUMIN
         INDEX=ISORDER(I)
         VQNE(INDEX)=VQNEIN(I)
      END DO

      POPEQ=0.0
      POPNE=0.0

      DO 50 LVL=1,MAXSTATE*Max_ISO
         VQEQ(LVL)=NDEG (LVL)*EXP(-EH(LVL)/XKT)
         VQNEST(LVL)=VQNE(LVL)
         POPEQ=POPEQ + VQEQ(LVL)
         POPNE=POPNE + VQNE (LVL)
   50 END DO

!    NORMALIZE POPULATIONS AND CALCULATE RATIOS

      DO 100 LVL=1,MAXSTATE*Max_ISO
         I=LVL
         VQEQ(LVL)=VQEQ(LVL)/POPEQ
         VQNE (LVL)=VQNE(LVL)/POPNE
         RAT(I)=VQNE(I)/VQEQ(I)
         IF(LVL.EQ.1) THEN
            VST1=VQNE(1)
            TNE(I)=TAVE
         ELSE
            DEN=(NDEG(1)*VQNE(LVL)/(NDEG(LVL)*VST1))
            TNE(I)=-RADCN2*EH(LVL)/ LOG(DEN)
         END IF
  100 END DO

      WRITE(IPR,906) TITMOL
      WRITE(IPR,935)
      DO J=1,NUMIN
         I=ISORDER(J)
         ISOTOPE=(I-1)/MAXSTATE + 1
         K=I-(ISOTOPE-1)*MAXSTATE
         WRITE(IPR,920) ISOTOPE,IDX(K),EH(I),VQEQ(I),VQNE(I),RAT(I),TNE(I),VQNEST(I)
      END DO

!     READ TO THE END OF THE VIBRATIONAL DATA

      CALL RDSKIP(NLTEFLAG,TEXTLINE)
      RETURN

  902 FORMAT(4x,A6)
  904 FORMAT (F7.0,1P,7E11.4,     /(18X, 6E11.4))
  906 FORMAT(//,A10,'  ENERGY LEVELS',10(/,20X,1PE11.4))
  920 FORMAT(I7,2X,A10,4G15.5,F10.2,G15.5)
  935 FORMAT ('ISOTOPE',2X,'VIB',10X,'E(CM-1)',11X,'POP LTE',7X,        &
     & 'POP NLTE',6X,'NLTE/LTE',7X,'NLTE TMP',7X,'NLTE POP ORIG')

      END SUBROUTINE


!-----------------------------------------------------------------------
!     SUBROUTINE VIBTMP USES THE NON-LTE TEMPERATURE DATA FROM
!     TAPE4 TO CALCULATE THE VIBRATIONAL POPULATION ENHANCEMENT
!     RATIOS FOR SELECTED VIBRATIONAL STATES OF H2O, CO2, NO AND
!     O3.  THE NLTE VIBRATIONAL TEMPERATURES ARE WITH RESPECT TO
!     THE GROUND VIBRATIONAL STATE.
!-----------------------------------------------------------------------
      SUBROUTINE VIBTMP( XKT,HT,NLTEFLAG,NUM,IDX,NDEG,EH,RAT,  &
                         TITMOL,TEXTLINE,ISORDER, TAVE )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8, RADCN2, MAXSTATE, MAX_ISO
      USE Module_Config       ,ONLY: IPR
      IMPLICIT NONE !REAL*8           (V)

      real         ,intent(out)   :: XKT
      real         ,intent(in)    :: HT
      integer      ,intent(in)    :: NLTEFLAG
      integer      ,intent(in)    :: NUM                    !for information
      character*5  ,intent(in)    :: IDX(:)     !(MAXSTATE) !for information
      integer      ,intent(in)    :: NDEG(:)    !(MAXSTATE*MAX_ISO)
      real         ,intent(in)    :: EH(:)      !(MAXSTATE*MAX_ISO)
      real         ,intent(out)   :: RAT(:)     !(MAXSTATE*MAX_ISO)
      character*6  ,intent(in)    :: TITMOL
      character*80 ,intent(inout) :: TEXTLINE
      integer      ,intent(in)    :: ISORDER(:) !(MAXSTATE*MAX_ISO)
      !
      real         ,intent(in)    :: TAVE


      !--- Local variables
      !
      character*6 ,PARAMETER ::  BLANKS=' '

      real(r8)    :: VQNE(   MAXSTATE*Max_ISO)
      real(r8)    :: VQEQ(   MAXSTATE*Max_ISO)
      real        :: TNE(    MAXSTATE*Max_ISO), &
                     TNESAV( MAXSTATE*Max_ISO), &
                     TEM1(   MAXSTATE*Max_ISO), &
                     TEM2(   MAXSTATE*Max_ISO), &
                     TNEIN(  MAXSTATE*Max_ISO)
      character*6 :: HMNLTE

      INTEGER  :: I,        I1,     I2,      INDEX
      INTEGER  :: ISOTOPE,  J,      K
      INTEGER  :: NUMIN,  NUMISO
      REAL     :: ALT1,     ALT2,   RATTV
      REAL     :: SUMNQ,    SUMQ




!     READ NLTE VIB TEMPERATURE
!     If this routine is called repeatedly for different altitudes,
!     this code is very inefficient, reading the input files every tim
!     to search for the proper altitude

      XKT= TAVE/RADCN2
      DO I=1,MAXSTATE*Max_ISO
         IF(ISORDER(I).GT.0) NUMIN=I
      end do
      write(ipr,*) 'called vibtmp with txtmol=',titmol,'  alt=',ht
      write(ipr,*) 'num=',num,'  numin=',numin
      write(ipr,*) 'isorder=',(isorder(i),i=1,numin)

      I1=0
      I2=0
      DO I=1,80
         IF(TEXTLINE(I:I).NE.'-' .AND. TEXTLINE(I:I).NE.' ') THEN
         IF(I1.EQ.0) I1=I
         END IF
         IF(I1.GT.0 .AND. TEXTLINE(I:I).EQ.' ') THEN
         IF(I2.EQ.0) I2=I-1
         END IF
      END DO
      HMNLTE=TEXTLINE(I1:I2)//BLANKS
      IF(HMNLTE.NE.TITMOL) THEN
         READ (NLTEFLAG,902) HMNLTE
         i1=0
         DO I=1,6
         IF(HMNLTE(I:I).EQ.' ') I1=I
         END DO
         HMNLTE=HMNLTE(I1+1:6)//BLANKS
      END IF
      IF(HMNLTE.NE.TITMOL) THEN
         WRITE(IPR,*) 'READING VIBTMP DATA FROM TAPE4'
         WRITE(IPR,*) 'EXPECTED PROFILE DATA FOR ',TITMOL
         WRITE(IPR,*) 'BUT READ SPECIE ',HMNLTE
         STOP 'ERROR READING NLTE DATA FROM TAPE4'
      END IF
      READ (NLTEFLAG,*)  ALT1,(TEM1(I),I=1,NUMIN)
   10 READ (NLTEFLAG,*)  ALT2,(TEM2(I),I=1,NUMIN)
      IF(TEM1(1).LE.1.E-12 .OR. TEM2(1).LE.1.E-12) THEN
         WRITE(IPR,*) 'VIBRATIONAL TEMPRATURE BASELINE FOR ',TITMOL,' MISSSING'  !WRITE(IPR,*) 'VIBRATIONAL TEMPRATURE BASELINE FOR ',TITMOL,' MISSSING'
         STOP 'KINETIC BASELINE TEMPERATURE MISSING'
      END IF
!*****WOG, 11/06/2000: ALT1 -> AL2:
!     IF( ALT1.LT.HT) THEN
      IF( ALT2.LE.HT) THEN
         ALT1 = ALT2
         DO 20 I=1,NUMIN
            TEM1(I) = TEM2(I)
   20    CONTINUE
         GO TO 10
      ENDIF

!     SET ZERO TEMP TO AMBIENT TEMP (NOW DONE AFTER INTERPOLATION)

!      write(ipr,941) 'tem1',alt1,(tem1(i),i=1,numin)
!      write(ipr,941) 'tem2',alt2,(tem2(i),i=1,numin)
      CALL LININT(HT,ALT1,ALT2,NUMIN,TEM1,TEM2,TNEIN)
!      write(ipr,941) 'TNEIN',ht,(tnein(i),i=1,numin)
      DO I=1,MAXSTATE
                              ! this rescales to ambient
         TNESAV(I)=TNEIN(1)
      END DO
      DO I=1,NUMIN
         INDEX=ISORDER(I)
                                     ! FOR ISOTOPE #1
         IF(INDEX.LE.MAXSTATE) THEN
         IF(TNEIN(I).GT.1.E-25) TNESAV(INDEX)=TNEIN(I)
         END IF
      END DO

!       FOR ISOTOPES>1, DEFAULT TO ISOTOPE #1 TEMPERATURE DATA

      DO ISOTOPE=2,Max_ISO
         DO K=1,MAXSTATE
            I= K + (ISOTOPE-1)*MAXSTATE
            TNESAV(I)=TNESAV(K)
         END DO
      END DO
      NUMISO=0
      DO I=1,NUMIN
         INDEX=ISORDER(I)
         ISOTOPE=(INDEX-1)/MAXSTATE +1
         IF(ISOTOPE.GT.NUMISO) NUMISO=ISOTOPE
         IF(ISOTOPE.GT.1) THEN
         IF(TNEIN(I).GT.1.E-25) THEN
            TNESAV(INDEX)=TNEIN(I)
         END IF
         END IF
      END DO
      WRITE(IPR,*) 'NUMBER OF ISOTOPES=',NUMISO
!      write(ipr,941) 'TNESAV',ht,(tnesav(i),i=1,MAXSTATE*NUMISO)
  941 format(a6,2x,27f8.3/26F8.3/26F8.3)

!     CORRECT TEMP TO ATMOSPHERIC

      RATTV = TAVE /TNESAV(1)
!     loop 40 corrects the input temperatures when the input
!     level temperature does not match the computed layer
!     temperature
      !write(*,*) ' temperature not corrected'

      DO 40 I = 1,MAXSTATE*Max_ISO
         TNE(I) = TNESAV(I) * RATTV
   40 END DO
!      write(ipr,941) 'TNE',ht,(tne(i),i=1,MAXSTATE*NUMISO)

      SUMQ=0
      SUMNQ=0
      DO 50 I=1,MAXSTATE*Max_ISO
         VQNE(I)=1.
         IF (TNE(I).GT.0.0) VQNE(I)=NDEG(I)*EXP(-RADCN2*EH(I)/TNE(I))
         VQEQ(I)=NDEG(I)*EXP(-EH(I)/XKT)
         SUMQ=SUMQ+VQEQ(I)
         SUMNQ=SUMNQ+VQNE(I)
   50 END DO
      DO 100 I=1,MAXSTATE*Max_ISO
         VQNE(I)=VQNE(I)/SUMNQ
         VQEQ(I)=VQEQ(I)/SUMQ
         IF(VQNE(I).GT.0.) THEN
            RAT(I)=VQNE(I)/VQEQ(I)
         ELSE
            RAT(I)=1.0
         END IF
  100 END DO

      WRITE(IPR,906) TITMOL
      WRITE(IPR,935)
      DO J=1,NUMIN
         I=ISORDER(J)
         ISOTOPE=(I-1)/MAXSTATE + 1
         K=I-(ISOTOPE-1)*MAXSTATE
         WRITE(IPR,920)ISOTOPE,IDX(K),EH(I),VQEQ(I),VQNE(I),RAT(I),TNESAV(I), TNE(I)
      END DO

!     READ TO THE END OF THE VIBRATIONAL DATA

      CALL RDSKIP(NLTEFLAG,TEXTLINE)
      RETURN

  902 FORMAT(4x,A6)
  904 FORMAT(F7.0,7F11.3/(18X,6F11.3))
  906 FORMAT(//,5X,A10,'  ENERGY LEVELS',10(/,20X,1PE11.4))
  920 FORMAT(I7,2X,A10,4G12.5,2F9.3)
  935 FORMAT ('ISOTOPE',9X,'VIB E(CM-1)',9X,'POP LTE    POP NLTE ',     &
     &     'NLTE/LTE    NLTE TMP 2-STATE NLTE TMP')

      END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE RDSKIP( NTAPE,TEXTLINE )
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer      ,intent(in)  :: NTAPE
      CHARACTER*80 ,intent(out) :: TEXTLINE

      !CHARACTER *1 HRD,HMINUS
      !DATA HMINUS /'-'/

   10 READ(NTAPE,900,END=50) TEXTLINE
  900 FORMAT(A80)
      IF(TEXTLINE(1:2).EQ.'--') RETURN
      GO TO 10
   50 RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE LININT( HT,ALT1,ALT2,NUM,T1,T2,TNE )
!-----------------------------------------------------------------------
      USE Module_Config ,ONLY:IPR
      IMPLICIT NONE

      integer ,intent(in)  :: NUM
      real    ,intent(in)  :: ALT1,ALT2
      real    ,intent(in)  :: HT
      real    ,intent(in)  :: T1(:)
      real    ,intent(in)  :: T2(:)
      real    ,intent(out) :: TNE(:)

      INTEGER :: I
      REAL    :: AM



      !*****WOG 11/03/2000
      !*****Correct for divide by zero if two altitudes are the same
      !*****0.001 = 1 meter, small enough to use average, large enough to
      !     prevent numerical error.
      IF (ABS(ALT2-ALT1) .LE. 0.001) THEN
         AM = 0.5
      ELSE
         AM = (HT-ALT1)/(ALT2-ALT1)
      ENDIF

      DO 10 I=1,NUM
         IF(ABS(T1(I)).LE.1.E-25 .AND. ABS(T2(I)).GT.1.E-25) GOTO 50
         IF(ABS(T1(I)).GT.1.E-25 .AND. ABS(T2(I)).LE.1.E-25) GOTO 50
         TNE(I)=T1(I)+AM*(T2(I)-T1(I))
   10 END DO

!     DO 10 I=1,NUM
!         AM = (T2(I)-T1(I))/(ALT2-ALT1)
!         C=T1(I)-AM*ALT1
!         TNE(I)=AM*HT+C
!  10  CONTINUE
      RETURN

   50 WRITE(IPR,*) 'ERROR IN LININT FOR Tvib INPUT'
      WRITE(IPR,*) 'Tvib=0 AT SOME ALTITUDES AND NOT OTHERS'
      write(ipr,*) 'ht=',ht,'  alt1,alt2=',alt1,alt2
      write(ipr,*) 'T1=',(t1(i),i=1,num)
      write(ipr,*) 'T2=',(t2(i),i=1,num)
      STOP 'ERROR IN LININT FOR Tvib INPUT'

      END SUBROUTINE



END MODULE
