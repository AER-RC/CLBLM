!
! CREATION HISTORY:
!       Modified from LBLRTM v12.9
!       Yingtao Ma, AER@NOAA/NESDIS/STAR
!       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_XSect
   USE Module_ConstParam ,ONLY: r8=>kind_r8,MX_XS, FILLINT

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: CLBLM_FSCDXS, xsTbl, &
             readXsMasterFile, &
             getXsInfo, &
             XSECTM


   !---
   ! NUMXS IS THE NUMBER OF 'CROSS SECTION' MOLECULES TO BE USED
   !
   ! XSFILE(ITEMP,ISPEC,NXS) IS THE NAME OF THE FILE CONTAINING THE
   !                         'CROSS SECTION' DATA.  THE THREE INDICES
   !                         ARE DEFINED AS FOLLOWS:
   !
   !                         ITEMP - DENOTES THE TEMPERATURE FOR WHICH
   !                                 THE 'CROSS SECTION' IS VALID
   !                                 (IMPLEMENTED FOR HITRAN 91 TAPE)
   !                         ISPEC - DENOTES THE SECTRAL REGION FOR
   !                                 WHICH THE FILE PERTAINS
   !                         NXS   - IS THE INCREMENT FOR THE 'CROSS
   !                                 SECTION' INDEX
   !
   ! NTEMPF(ISPEC,NXS) IS THE NUMBER OF TEMPERATURE FILES TO BE USED
   !                   FOR EACH SPECTRAL REGION OF EACH MOLECULE
   !
   ! NSPECR(NXS) IS THE NUMBER OF SPECTRAL REGIONS FOR THE MOLECULE NX
   !
   integer ,parameter :: MaxNumTemp = 6
   integer ,parameter :: MaxNumSpec = 6
   TYPE :: CLBLM_FSCDXS
      integer       :: numXs                                  =0
      character(10) :: xsName(                         MX_XS) =''
      integer       :: numSpect(                       MX_XS) =0    !NSPECR
      real(r8)      :: V1FX(               MaxNumSpec, MX_XS) =0.0
      real(r8)      :: V2FX(               MaxNumSpec, MX_XS) =0.0
      real          :: DVFX(               MaxNumSpec, MX_XS) =0.0
      integer       :: numTemp(            MaxNumSpec, MX_XS) =0    !NTEMPF
      integer       :: xsForm(             MaxNumSpec, MX_XS) =0    !IXFORM
      character(10) :: xsFile( MaxNumTemp, MaxNumSpec, MX_XS) =''
      real          :: XDOPLR(             MaxNumSpec, MX_XS) =0.0  !XDOPLR
      integer       :: lineOrXs(           MaxNumSpec, MX_XS) =FILLINT ! 0: only XS exist for this molecule in this spectral range
                                                                       ! 1: both line parameters and XS exist, line parameters are preferable
                                                                       ! 2: both line parameters and XS exist, XS are preferable
   END TYPE

   type( CLBLM_FSCDXS ), PROTECTED, SAVE :: xsTbl



CONTAINS !=====================Module Contians==========================



!-----------------------------------------------------------------------
! Read in the contents of FSCDXS and save them into a table. The order
! or xs molecules in the table follows the order of molecules in ALIAS arrays.
!-----------------------------------------------------------------------
   SUBROUTINE readXsMasterFile( xsMastFile )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: xsMolNum, XSMASS, MX_XS
      IMPLICIT NONE

      character(*) ,intent(in)  :: xsMastFile

      character(*) ,parameter :: routineName='readXsMasterFile'

      ! T296 IS TEMPERATURE FOR INITAL CALCULATIN OF DOPPLER WIDTHS
      real      ,PARAMETER :: T296=296.0
      character ,PARAMETER :: CASTSK='*'
      character ,PARAMETER :: CPRCNT='%'
      character ,PARAMETER :: CN    ='N'
      character ,PARAMETER :: CF    ='F'

      integer       :: ix, J, numXs, IOSTAT, IXFIL
      character     :: XSREC*130
      CHARACTER(10) :: XNAME, oldXNAME
      REAL(r8)      :: V1X, V2X
      REAL          :: DVX
      INTEGER       :: NTEMP
      integer       :: IFRM
      CHARACTER     :: CFRM
      character(10) :: XFILS(MaxNumTemp)
      integer       :: lnOxs
      integer       :: NSPECR(MX_XS)
      integer       :: IXFORM(MaxNumSpec, MX_XS)


      ! READ IN "CROSS SECTION" MASTER FILE FSCDXS
      IXFIL = 8
      OPEN (IXFIL,FILE=trim(xsMastFile), STATUS='OLD',FORM='FORMATTED',IOSTAT=iostat)
      if (IOSTAT.gt.0) STOP '--- '//routineName//'(): FSCDXS does not exist'

      !--- Read the header line
      REWIND IXFIL
      READ (IXFIL,905)

      !--- Loop over the records in FSCDXS
      oldXNAME = ''
      NSPECR(:) = 0
      numXs = 0
   50 READ (IXFIL,'(A)',END=80) XSREC

         IF (XSREC(1:1).EQ.CASTSK) GO TO 50 !IF (CFLG.EQ.CASTSK) GO TO 50
         IF (XSREC(1:1).EQ.CPRCNT) GO TO 80 !IF (CFLG.EQ.CPRCNT) GO TO 80

         READ (XSREC,915) XNAME,V1X,V2X,DVX,NTEMP,IFRM,CFRM, (XFILS(J),J=1,MaxNumTemp), lnOxs

         XNAME=adjustl(XNAME)
         if ( trim(adjustl(XNAME)) /= trim(adjustl(oldXNAME)) ) then
            ix = xsMolNum( XNAME )
            if (ix<0) then !no matching found in ALIAS(1:4), ignore it
               !print*, ('--- '//routineName//'(): '//XNAME//'from FSCDXS not found in XS ALIAS tables.')
               GO TO 50
            endif
            xsTbl%xsName(ix) = XNAME
            oldXNAME = XNAME
            numXs = numXs+1
         endif

         NSPECR(ix) = NSPECR(ix)+1
         IF (NSPECR(ix).GT.6) THEN
            STOP '--- '//routineName//'(): NSPECR .GT. 6'
         ENDIF
         xsTbl%numSpect(ix) = NSPECR(ix)

         xsTbl%V1FX(NSPECR(ix),ix) = V1X
         xsTbl%V2FX(NSPECR(ix),ix) = V2X
         xsTbl%DVFX(NSPECR(ix),ix) = DVX
         xsTbl%numTemp(NSPECR(ix),ix) = NTEMP

         IXFORM(NSPECR(ix),ix) = 91
         IF (IFRM.EQ.86) IXFORM(NSPECR(ix),ix) = IFRM
         IF (CFRM.NE.CN) IXFORM(NSPECR(ix),ix) = IXFORM(NSPECR(ix),ix)+100
         IF (CFRM.EQ.CF) IXFORM(NSPECR(ix),ix) = -IXFORM(NSPECR(ix),ix)
         xsTbl%xsForm(NSPECR(ix),ix) = IXFORM(NSPECR(ix),ix)

         DO J = 1, NTEMP
            xsTbl%xsFile(J,NSPECR(ix),ix) = XFILS(J)
         ENDDO

         !3.58115E-07 = SQRT( 2.* LOG(2.)*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) )
         !XDOPLR(NSPECR(ix),ix) = 3.58115E-07 * ( 0.5*(V1X+V2X) ) * SQRT( T296/XSMASS(IXINDX(ix)) )
         xsTbl%XDOPLR(NSPECR(ix),ix) = 3.58115E-07 * ( 0.5*(V1X+V2X) ) * SQRT( T296/XSMASS( xsMolNum(XNAME) ) )
         xsTbl%lineOrXs(NSPECR(ix),ix) = lnOxs

      GO TO 50

   80 continue
      xsTbl%numXs = numXs

  905 FORMAT (/)
!----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3
!F11         830.0010  859.9960 .01482699    6          N    xs/F11AT6 xs/F11AT5 xs/F11AT4 xs/F11AT3 xs/F11AT2 xs/F11AT1     0
  915 FORMAT (A10,2F10.4,F10.8,I5,5X,I5,A1,4X,6A10,I5)
  920 FORMAT (/,'******* ERROR IN XSREAD ** MOLECULE SECLECTED -',A10,  &
     &        '- HAS ',I2,' SPECTRAL REGIONS ON FILE FSCDXS, BUT THE',  &
     &        ' MAXIMUM ALLOWED IS 6 *******',/)
   END SUBROUTINE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE getXsInfo( XV1, XV2, IXMOLS, xsNames, &
                            XSFILE, IXFORM, NSPECR, NTEMPF, XDOPLR, LnOrXs )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8, ALIAS, MX_XS, molIndex, xsMolNum
      USE Module_Config       ,ONLY: IPR
      IMPLICIT NONE

      real(r8)           ,intent(in)  :: XV1,XV2
      integer            ,intent(in)  :: IXMOLS        !NUMBER OF THESE MOLECULES SELECTED,
      character(*)       ,intent(in)  :: xsNames(:)    !(MX_XS)
      character(10)      ,intent(out) :: XSFILE(:,:,:) !(6,5,MX_XS)
      integer            ,intent(out) :: IXFORM(:,:)   !(5,MX_XS)
      !real(r8)           ,intent(out) :: V1FX(:,:)     !(5,MX_XS)
      !real(r8)           ,intent(out) :: V2FX(:,:)     !(5,MX_XS)
      integer            ,intent(out) :: NSPECR(:)     !(MX_XS)
      integer            ,intent(out) :: NTEMPF(:,:)   !(5,MX_XS)
      real               ,intent(out) :: XDOPLR(:,:)   !(5,MX_XS)
      integer            ,intent(out) :: LnOrXs(:,:)   !(5,MX_XS)


      !--- Local variables
      INTEGER   :: I,J, ind,ia, ix,is,nspc,nt
      integer   :: IXINDX(MX_XS)    !INDEX VALUES OF SELECTED MOLECULES (E.G. 1=CLONO2),
      REAL(r8)  :: V1X, V2X



      ! MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS
      ! AND DETERMINE THE INDEX VALUE.  STOP IF NO MATCH IS FOUND.
      ! NAME MUST BE ALL IN CAPS.
      DO I = 1, IXMOLS
         IXINDX(I) = xsMolNum( xsNames(I) )
         if ( IXINDX(I) <0 ) then  ! NO MATCH FOUND
            WRITE (IPR,900) xsNames(I)
            STOP 'STOPPED IN XSREAD'
         endif
      ENDDO


      !--- Look up the xsTbl
      do ix = 1,IXMOLS

         !--- Search the xs in xsTbl
         do ia = 1,4 !4 aliases
            ind = molIndex( ALIAS(ia,IXINDX(ix)), xsTbl%xsName )
            if (ind>0) EXIT
         enddo

         !--- If not found in the xsTbl, stop
         if (ind<0) then
            WRITE (*,925) xsNames(ix)
            STOP ' IXFLAG - XSREAD '
         endif

         !--- Found in the xsTbl, read the xs information.
         ! Only xs intervals cover XV1~XV2 are output
         !
         nspc = 0
         do is = 1, xsTbl%numSpect( ind )

            V1X = xsTbl%V1FX(is,ind)
            V2X = xsTbl%V2FX(is,ind)

            IF (V2X.GT.XV1.AND.V1X.LT.XV2) THEN

               nspc = nspc+1
               nt = xsTbl%numTemp(is,ind)

               XSFILE(1:nt,nspc,ix) = xsTbl%xsFile(1:nt,is,ind)
               IXFORM(     nspc,ix) = xsTbl%xsForm(is,ind)
               !V1FX(       nspc,ix) = V1X
               !V2FX(       nspc,ix) = V2X
               NTEMPF(     nspc,ix) = nt
               XDOPLR(     nspc,ix) = xsTbl%XDOPLR(is,ind)
               LnOrXs(     nspc,ix) = xsTbl%lineOrXs(is,ind)
            endif
         enddo

         NSPECR(ix) = nspc !xsTbl%numSpect(ind)
         IF (NSPECR(ix).GT.6) THEN
            WRITE (*,920) ix,xsNames(ix),NSPECR(ix)
            STOP ' XSREAD - NSPECR .GT. 6'
         ENDIF

      enddo !do ix = 1,IXMOLS

      RETURN

  900 FORMAT (/,'  THE NAME: ',A10, ' IS NOT ONE OF THE ',              &
     &        'CROSS SECTION MOLECULES. CHECK THE SPELLING.')
  920 FORMAT (/,'******* ERROR IN XSREAD ** MOLECULE SECLECTED -',A10,  &
     &        '- HAS ',I2,' SPECTRAL REGIONS ON FILE FSCDXS, BUT THE',  &
     &        ' MAXIMUM ALLOWED IS 6 *******',/)
  925 FORMAT (/,'******* MOLECULE SELECTED -',A10,'- IS NOT FOUND ON',  &
     &        ' FILE FSCDXS *******',/)

      END SUBROUTINE


!-----------------------------------------------------------------------
!     THIS SUBROUTINE MOVES THE CROSS SECTIONS INTO
!     THE APPROPRIATE ARRAY R1, R2, R3, R4, OR ABSRB
!-----------------------------------------------------------------------
      SUBROUTINE XSECTM( IFST,IR4,&
                         R1,R2,R3,R4,lbR1,lbR2,lbR3,lbR4,&
                         V1,V2,V1R4,V2R4,DV,DVR2,DVR3,DVR4,&
                         N1R1,N2R1,N1R2,N2R2,N1R3,N2R3,NPTR4,NHI,VFT, &
                         PAVE,TAVE,WXM,&
                         XSFILE,IXFORM,NSPECR,NTEMPF,XDOPLR,LnOrXs,NUMXS,&
                         IXSBIN,JRAD,ILBLF4,DPTMIN,NFHDRF,NPHDRF )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, MX_XS, RADCN2
      USE Module_ODLAY_CommonSub  ,ONLY: RADFNI
      USE Module_Spectrum         ,ONLY: XINT
      USE Module_Config           ,ONLY: IPR
      IMPLICIT NONE !REAL*8           (V)

      integer       ,intent(inout) :: IFST
      integer       ,intent(inout) :: IR4
      !
      integer       ,intent(in)    :: lbR1,lbR2,lbR3,lbR4 !lower bounds of R1,R2,R3,and R4
      real          ,intent(inout) :: R1(lbR1:),R2(lbR2:),R3(lbR3:),R4(lbR4:)
      real(r8)      ,intent(in)    :: V1,V2
      real(r8)      ,intent(in)    :: V1R4,V2R4
      real          ,intent(in)    :: DV
      real          ,intent(in)    :: DVR2
      real          ,intent(in)    :: DVR3
      real          ,intent(in)    :: DVR4
      integer       ,intent(in)    :: N1R1,N2R1
      integer       ,intent(in)    :: N1R2,N2R2
      integer       ,intent(in)    :: N1R3,N2R3
      integer       ,intent(in)    :: NPTR4
      integer       ,intent(in)    :: NHI
      real(r8)      ,intent(in)    :: VFT
      real          ,intent(in)    :: PAVE
      real          ,intent(in)    :: TAVE
      real          ,intent(in)    :: WXM(:)        !(MX_XS)
      character(10) ,intent(in)    :: XSFILE(:,:,:) !(6,5,MX_XS)
      integer       ,intent(in)    :: IXFORM(:,:)   !(5,MX_XS)
      integer       ,intent(in)    :: NSPECR(:)     !(MX_XS)
      integer       ,intent(in)    :: NTEMPF(:,:)   !(5,MX_XS)
      real          ,intent(in)    :: XDOPLR(:,:)   !(5,MX_XS)
      integer       ,intent(in)    :: LnOrXs(:,:)   !(5,MX_XS)
      integer       ,intent(in)    :: NUMXS
      integer       ,intent(in)    :: IXSBIN
      integer       ,intent(in)    :: JRAD
      integer       ,intent(in)    :: ILBLF4
      real          ,intent(in)    :: DPTMIN
      integer       ,intent(in)    :: NFHDRF,NPHDRF



      !---Saved local variables
      !
      integer        ,SAVE :: NFILEX(5,MX_XS)   =0 !yma added
      real           ,SAVE :: XSTEMP(6,5,MX_XS)
      real           ,SAVE :: XSMAX(6,5,MX_XS)  =0.0
      real           ,SAVE :: PDX(6,5,MX_XS)
      real(r8)       ,SAVE :: V1FX(5,MX_XS)     =0.0
      real(r8)       ,SAVE :: V2FX(5,MX_XS)     =0.0
      real           ,SAVE :: DVFX(5,MX_XS)     =0.0
      integer        ,SAVE :: NPTSFX(5,MX_XS)
      real(r8)       ,SAVE :: V1XS              =0.0
      real(r8)       ,SAVE :: V2XS              =0.0
      real           ,SAVE :: DVXS              =0.0
      integer        ,SAVE :: NPTSXS            =0
      real           ,SAVE :: RX(13000)
      real(r8)       ,SAVE :: V1X
      real(r8)       ,SAVE :: V2X
      real           ,SAVE :: DVX
      integer        ,SAVE :: NPTSX
      integer        ,SAVE :: NMODES            =1
      integer        ,SAVE :: NPANEL            =0
      integer        ,SAVE :: JINPUT
      integer        ,SAVE :: IXBIN(5,MX_XS)
      integer        ,SAVE :: IXSBN(5,MX_XS)
      real           ,SAVE :: DVXPR(5,MX_XS)
      character(100) ,SAVE :: HEADT1(MX_XS)


      !--- Non-saved local variables
      !
      integer ,PARAMETER :: IFILE=115
      integer ,PARAMETER :: JFILE=116
      integer ,PARAMETER :: I_0=0

      integer   :: NLIMX
      real(r8)  :: V1DX, V2DX
      real      :: DVDX
      integer   :: NPTSDX
      real      :: RDX1(520)
      real      :: RDX2(520)
      real      :: PF, TF

      INTEGER  :: I,      IAFORM,  IFL
      INTEGER  :: IMAX,   IMAXSV, IRPEAT
      INTEGER  :: JI,     JMAX,   KI
      INTEGER  :: LIMOUT, N1RX,   N2RX
      INTEGER  :: NBSKIP, NEOF,   NFILET, NI
      INTEGER  :: NMAX,   NMODE,  NNSKIP, NP
      INTEGER  :: NPAN,   NPTSI1, NPTSI2, NPTST
      INTEGER  :: NRSKIP, NS,     NSKIP,  NT1
      INTEGER  :: NT2
      LOGICAL  :: OPCL
      REAL     :: DVFXX,  DVXMIN, RADVI
      REAL     :: RDEL,   RDLAST, TFACT1, TFACT2
      REAL(r8) :: V1FP,   V1XMIN, V1XT,   V2FP
      REAL(r8) :: V2XT,   VFX2,   VI,     VITST
      REAL(r8) :: VREF
      REAL     :: WXM1,   WXM2,   XKT




      NLIMX = 510
      LIMOUT = 13000
      IRPEAT = 0

      PF = PAVE
      TF = TAVE

!     IF FIRST ENTRANCE, RESET QUANTITES

      IF (IFST.EQ.-99) THEN !IFST set to -99 for each layer
         IFST = 0
         JINPUT = -1
         NMODES = 1

         V1X = V1XS
         V2X = V2XS
         DVX = DVXS
         NPTSX = NPTSXS

         DO 10 NI = 1, NUMXS
            DO 9 NS = 1, NSPECR(NI)
               IXBIN(NS,NI) = 1
               IXSBN(NS,NI) = 0
               NFILEX(NS,NI) = ABS(NFILEX(NS,NI))
    9       continue
   10    CONTINUE
      ENDIF

!     CHECK V1X FOR INPUT

   20 VFX2 = VFT+2.*DVX+ REAL(NHI)*DV  !yma: wanted end point of the XS range for this panel
      IF (IR4.EQ.1) VFX2 = V2R4+2.*DVX
      IF (V1X.GT.VFX2) GO TO 140

      IF (JINPUT.EQ.-1) THEN
         JINPUT = 1
      ELSE
         VFX2 = MIN(VFX2,V2+2.*DVX)
         IF (VFX2.GT.V2X) THEN  !yma: V2X is available end point of the XS range
            JINPUT = 1
            IF (IRPEAT.EQ.1) THEN
               V1X = V2X-2.*DVX
            ELSE
               V1X = VFT-2.*DVX
!Matt Alvarado 20150819 Added Yingtao Ma's fix for V1X
               Vref = V1-32*DV
               V1X = Vref + floor((V1X-Vref)/DVx)*DVX
               IF (IR4.EQ.1) V1X = V1R4-2.*DVX
            ENDIF
            V2X = V1X+ REAL(LIMOUT-1)*DVX
            IF (V2X.GT.V2) V2X = V1X+ REAL(INT((V2-V1X)/DVX)+3)*DVX
!Matt Alvarado 20150819 Added Yingtao Ma's fix for NPTSX
            NPTSX = ceiling((V2X-V1X)/DVX+1.)
         ENDIF

         !--- Check if there are xs lies within the interval V1XT ~V2XT
         IFL = 0
         V1XT = V1X+2.*DVX
         V2XT = V2X-2.*DVX
         IF (V1XT.GT.V2XT) GO TO 140
         DO 30 NI = 1, NUMXS
            DO 29 NS = 1, NSPECR(NI)
               IF (NFILEX(NS,NI).EQ.0) GO TO 30
               IF (V1FX(NS,NI).LE.V1XT.AND.V2FX(NS,NI).GE.V1XT) IFL = 1
               IF (V1FX(NS,NI).LE.V2XT.AND.V2FX(NS,NI).GE.V2XT) IFL = 1
               IF (V1FX(NS,NI).GE.V1XT.AND.V2FX(NS,NI).LE.V2XT) IFL = 1
   29       CONTINUE
   30    CONTINUE

         !--- No xs found in this interval
         IF (IFL.EQ.0) GO TO 140
      ENDIF


      !--- JINPUT==1, need to READ IN CROSS SECTION and calculate xs OD. (xs OD in RX array.)
      !
      IF (JINPUT.EQ.1) THEN
         JINPUT = 0
         DO 40 I = 1, LIMOUT
            RX(I) = 0.0
   40    CONTINUE
         DO 50 I = 1, NLIMX+10
            RDX1(I) = 0.0
            RDX2(I) = 0.0
   50    CONTINUE

!     FOR NPANEL = 0, READ IN FILE HEADERS

         IF (NPANEL.EQ.0) THEN
            DVXMIN = V2-V1
            V1XMIN = V2
            IMAX = 0
            NT2 = 0
            NMODE = 0
            DO 60 NI = 1, NUMXS
               DO 59 NS = 1, NSPECR(NI)
                  DO 58 NT1 = 1, NTEMPF(NS,NI)
                     NFILEX(NS,NI) = 1
                     !CALL CPUTIM (TIME0)
                     CALL XSECIN( NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,IMAX,NEOF, &
                                  XSFILE,IXFORM,IXSBN,HEADT1,V1DX,V2DX,DVDX,NPTSDX, &
                                  XSTEMP,XSMAX,PDX,RDX1(5:),RDX2(5:),NFHDRF,NPHDRF )
                     !CALL CPUTIM (TIME)
                     !TXSRDF = TXSRDF+TIME-TIME0

!     CHECK FOR WAVENUMBER BOUNDS AND SMALLEST DV

                     IF (V1DX.GT.V2.OR.V2DX.LT.V1.OR.NEOF.EQ.1) THEN !xs interval is out of V1~V2 range
                        NFILEX(NS,NI) = 0
                     elseif ( LnOrXs(NS,NI)==1 ) then !both line data and xs are available for this interval, line data are used.
!ab.OrderMol
                        NFILEX(NS,NI) = 0
                     ELSE
                        DVXMIN = MIN(DVXMIN,DVDX)
                        V1XMIN = MIN(V1XMIN,V1DX)
                        V1FX(NS,NI) = V1DX
                        V2FX(NS,NI) = V2DX
                        DVFX(NS,NI) = DVDX
                        NPTSFX(NS,NI) = NPTSDX
                     ENDIF

!mji  CHECK FOR TEMPERATURES; MUST BE IN ASCENDING ORDER

                     IF (NT1.GT.1) then
                        if (XSTEMP(NT1,NS,NI).LT.XSTEMP(NT1-1,NS,NI)) THEN
                           WRITE(IPR,900)
                           STOP 'XSTEMP - XSECTM'
                        endif
                     ENDIF
   58             CONTINUE
   59          CONTINUE
   60       CONTINUE
            DVX = DVXMIN
            V1X = MAX(VFT,V1XMIN)
            V1X = V1X-2.*DVX
            V2X = V1X+ REAL(LIMOUT-1)*DVX
            IF (V2X.GT.V2) V2X = V1X+ REAL(INT((V2-V1X)/DVX)+2)*DVX
            NPTSX = ceiling((V2X-V1X)/DVX+1.)
            V1XS = V1X
            V2XS = V2X
            DVXS = DVX
            NPTSXS = NPTSX
            IF (V1X.GT.VFX2) THEN
               JINPUT = 1
               NPANEL = -1
               GO TO 140
            ENDIF
         ENDIF !IF (NPANEL.EQ.0) THEN

         NFILET = 0
         NMODES = 0

         DO 110 NI = 1, NUMXS
            DO 109 NS = 1, NSPECR(NI)
               NPANEL = -1
               IF (NFILEX(NS,NI).LE.0) GO TO 105
               IF (V1FX(NS,NI).GT.V2X) GO TO 105
               IF (V2FX(NS,NI).LT.V1X) THEN
                  NFILEX(NS,NI) = -NFILEX(NS,NI)
                  GO TO 105
               ENDIF

!     DETERMINE TEMPERATURE FILES AND TEST ON DPTMIN

               CALL XSNTMP( NI,NS,NT1,NT2,NMODE, &
                            TAVE,XSTEMP,NTEMPF,V1FX,WXM,XSMAX,DPTMIN,JRAD )

!     DPTMIN TEST - IF NMODE = 0, SKIP CROSS SECTION

               IF (NMODE.EQ.0) GO TO 105
               NMODES = NMODES+NMODE

!     FOR PRESSURE BROADENED CROSS-SECTION
!     CREATE TEMPERATURE AVERAGED BINARY FILE

               IF (IXSBIN.EQ.0.AND.IXBIN(NS,NI).EQ.1) THEN
                  !CALL CPUTIM (TIME0)
                  CALL XSBINF( NI,NS,NT1,NT2,NMODE, &
                               IXSBN,XSFILE,IXFORM,V1FX,V2FX,DVFX,NPTSFX,&
                               XDOPLR,PF,TF,TAVE,XSTEMP,PDX,HEADT1,NFHDRF,NPHDRF,DVXPR )
                  !CALL CPUTIM (TIME)
                  !TXSCNV = TXSCNV+TIME-TIME0
                  IXBIN(NS,NI) = 0
               ENDIF
               IF (IXSBN(NS,NI).EQ.1) THEN
                  DVFXX = DVXPR(NS,NI)
               ELSE
                  DVFXX = DVFX(NS,NI)
               ENDIF
               NFILET = NFILET+NFILEX(NS,NI)

               NNSKIP = (V1X-V1FX(NS,NI))/DVFXX
               NSKIP = (NNSKIP-3)/10
               NSKIP = MAX(NSKIP,I_0)
               NRSKIP = NSKIP*10
               NBSKIP = NSKIP

!     FOR BLOCKED DATA, V1FP MUST REFLECT SHORT RECORD

               IAFORM = ABS(IXFORM(NS,NI))
               IF (IXSBN(NS,NI).EQ.1) IAFORM = IAFORM+100
               IF (IAFORM.GT.100) THEN
                  NBSKIP = NSKIP/51
                  NRSKIP = (NBSKIP-1)*510+500
                  NRSKIP = MAX(NRSKIP,I_0)
               ENDIF
               V1FP = V1FX(NS,NI)+ REAL(NRSKIP)*DVFXX
               V2FP = V2X+2.0*DVFXX
               V2FP = MIN(V2FP,V2FX(NS,NI))
               NMAX = (V2FP-V1FP)/DVFXX+1.
               NPAN = (NMAX+NLIMX-1)/NLIMX
               IF (IAFORM.GT.100.AND.NPANEL.LE.0.AND.NBSKIP.EQ.0) THEN
                  NPTST = NMAX-500-(NPAN-1)*NLIMX
                  IF (NPTST.GT.0) NPAN = NPAN+1
                  IF (NMAX.GT.500) NMAX = NMAX+10
               ENDIF
               N2RX = ((V1FP-4.*DVFXX-V1X)/DVX+0.999)-1.
               N2RX = MAX(N2RX,I_0)

!     IMAX = -4 TO PLACE THE FIRST PANEL V1 AT ARRAY LOCATION 1

               IMAX = -4
               DO 100 NP = 1, NPAN
                  V1FP = V1FP+ REAL(IMAX)*DVFXX
                  IMAX = NMAX-(NP-1)*NLIMX
                  IF (IAFORM.GT.100 .AND. NPANEL.LE.0 .AND. NBSKIP.EQ.0 &
                 & .AND. IMAX.GT.500) IMAX = 500
                  IMAX = MIN(IMAX,NLIMX)

!     FOR V2FP IMAX + 3 GIVES US ARRAY LOCATION 514
!             (504 FOR FIRST PANEL OF BLOCKED DATA)

                  V2FP = V1FP+ REAL(IMAX+3)*DVFXX

                  IF (NP.GT.1) THEN
                     DO 70 JI = 1, 4
                        RDX1(JI) = RDX1(JI+IMAXSV)
                        RDX2(JI) = RDX2(JI+IMAXSV)
   70                CONTINUE
                  ENDIF
                  IMAXSV = IMAX

                  !CALL CPUTIM (TIME0)
                  CALL XSECIN( NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,IMAX,NEOF, &
                               XSFILE,IXFORM,IXSBN,HEADT1,V1DX,V2DX,DVDX,NPTSDX,&
                               XSTEMP,XSMAX,PDX,RDX1(5:),RDX2(5:),NFHDRF,NPHDRF )
                  !CALL CPUTIM (TIME)
                  !TXSRDF = TXSRDF+TIME-TIME0

                  IF (NP.EQ.1) THEN
                     DO 80 JI = 1, 4
                        RDX1(JI) = RDX1(5)
                        RDX2(JI) = RDX2(5)
   80                CONTINUE
                     NPANEL = ABS(NPANEL)
                     NSKIP = 0
                  ENDIF

!     IF LAST PANEL OF FILE, FILL IN ADDITIONAL POINTS
!     TO ENSURE ENDPOINT CAPTURE

                  IF (NP.EQ.NPAN) THEN
                     JMAX = IMAX+4
                     DO 90 JI = 1, 4
                        KI = JMAX+JI
                        RDX1(KI) = RDX1(JMAX)
                        RDX2(KI) = RDX2(JMAX)
   90                CONTINUE

                     V2FP = V2FP+3.*DVFXX
                  ENDIF

                  N1RX = MAX(1,N2RX+1)
                  N2RX = (V2FP-DVFXX-V1X)/DVX+.999
                  N2RX = MIN(N2RX,LIMOUT)

                  WXM1 = WXM(NI)

!     FOR TWO TEMPERATURES LINEARLY INTERPOLATE FACTOR

                  !CALL CPUTIM (TIME0)
                  IF (NMODE.EQ.2.AND.IXBIN(NS,NI).EQ.1) THEN
                     TFACT2 = (TAVE-XSTEMP(NT1,NS,NI))/ (XSTEMP(NT2,NS, &
                     NI)-XSTEMP(NT1,NS,NI))
                     TFACT1 = 1.-TFACT2
                     WXM1 = WXM(NI)*TFACT1
                     WXM2 = WXM(NI)*TFACT2
                     CALL XINT (V1FP,V2FP,DVFXX,RDX2,WXM2,V1X,DVX,RX,N1RX,N2RX)
                  ENDIF
                  CALL XINT (V1FP,V2FP,DVFXX,RDX1,WXM1,V1X,DVX,RX,N1RX,N2RX)
                  !CALL CPUTIM (TIME)
                  !TXSPNL = TXSPNL+TIME-TIME0

  100          CONTINUE

!              Continue for GOTO statements at E04540, E04550, E04580, &
!              E04670.

  105          CONTINUE

  109       CONTINUE
  110    CONTINUE
         IF (NFILET.EQ.0) GO TO 140

!     FACTOR OUT RADIATION FIELD IF REQUIRED

         IF (JRAD.EQ.0) THEN
            !CALL CPUTIM (TIME0)
            XKT = TAVE/RADCN2
            VI = V1X-DVX
            VITST = VI
            RDLAST = -1.
            NPTSI1 = 0
            NPTSI2 = 0

  120       NPTSI1 = NPTSI2+1
            NPTSX = ceiling((V2X-V1X)/DVX+1.)
            NPTSX = min(NPTSX,LIMOUT) !yma 171220; NPTSX may reach to 13001

            VI = V1X+ REAL(NPTSI1-1)*DVX
            RADVI = RADFNI(VI,DVX,XKT,VITST,RDEL,RDLAST)

!v128            NPTSI2 = (VITST-V1X)/DVX+1.001
            NPTSI2 = (VITST-V1X)/DVX+0.001
            NPTSI2 = MIN(NPTSI2,NPTSX)

            DO 130 I = NPTSI1, NPTSI2
               VI = VI+DVX
               RX(I) = RX(I)/RADVI
               RADVI = RADVI+RDEL
  130       CONTINUE

            IF (NPTSI2.LT.NPTSX) GO TO 120
            !CALL CPUTIM (TIME)
            !TXSPNL = TXSPNL+TIME-TIME0

         ENDIF
      ENDIF  !IF (JINPUT.EQ.1) THEN

      IF (NMODES.EQ.0) GO TO 140


!     DETERMINE TARGET ARRAY

!      ===> R1

      !CALL CPUTIM (TIME0)
      IF (DVX.LT.DVR2) THEN
         CALL  XINT( V1X,V2X,DVX,RX,1.0,VFT,DV,R1(1:),N1R1,N2R1 ) !XINT( V1X,V2X,DVX,RX,1.0,VFT,DV,R1,N1R1,N2R1)
         IR4 = 0

!      ===> R2

      ELSEIF (DVX.LT.DVR3) THEN
         CALL XINT( V1X,V2X,DVX,RX,1.0,VFT,DVR2,R2(1:),N1R2,N2R2 ) !XINT (V1X,V2X,DVX,RX,1.0,VFT,DVR2,R2,N1R2,N2R2)
         IR4 = 0

!      ===> R3

      ELSEIF (DVX.LT.DVR4.OR.ILBLF4.EQ.0) THEN
         CALL XINT( V1X,V2X,DVX,RX,1.0,VFT,DVR3,R3(1:),N1R3,N2R3 ) !XINT (V1X,V2X,DVX,RX,1.0,VFT,DVR3,R3,N1R3,N2R3)
         IR4 = 0

!      ===> R4

      ELSE
         CALL XINT (V1X,V2X,DVX,RX,1.0,V1R4,DVR4,R4(1:),1,NPTR4) !XINT (V1X,V2X,DVX,RX,1.0,V1R4,DVR4,R4,1,NPTR4)
         IF (IR4.EQ.0) VFX2 = V2R4+2.*DVX
         IR4 = 1
      ENDIF
      !CALL CPUTIM (TIME)
      !TXSPNL = TXSPNL+TIME-TIME0

      IRPEAT = 1
      IF (VFX2.GT.V2X) GO TO 20

  140 INQUIRE (UNIT=IFILE,OPENED=OPCL)
      IF (OPCL) CLOSE (IFILE)
      INQUIRE (UNIT=JFILE,OPENED=OPCL)
      IF (OPCL) CLOSE (JFILE)
      IF (NMODES.EQ.0) IR4 = 1

  900 FORMAT(/,'******* ERROR IN XSECTM *******',/                      &
     &         'CROSS-SECTION FILES MUST BE IN ASCENDING ORDER ',       &
     &         'BY TEMPERATURE IN FSCDXS.')
      RETURN

      END SUBROUTINE

!-----------------------------------------------------------------------
!     THIS SUBROUTINE READS IN THE DESIRED CROSS SECTIONS
!-----------------------------------------------------------------------
      SUBROUTINE XSECIN( NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,NMAX,IEOF, &
                         XSFILE,IXFORM,IXSBN,HEADT1,V1DX,V2DX,DVDX,NPTSDX,&
                         XSTEMP,XSMAX,PDX,RDX1,RDX2,NFHDRF,NPHDRF )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Config      ,ONLY: ioFiles
      IMPLICIT NONE !REAL*8           (V)

      integer        ,intent(in)    :: NPANEL
      integer        ,intent(in)    :: NI,NS
      integer        ,intent(in)    :: NT1,NT2
      integer        ,intent(in)    :: NMODE
      integer        ,intent(inout) :: NSKIP
      integer        ,intent(in)    :: NMAX
      integer        ,intent(out)   :: IEOF
      !
      character(10)  ,intent(in)    :: XSFILE(:,:,:) !(6,5,MX_XS)
      integer        ,intent(in)    :: IXFORM(:,:)   !(5,MX_XS)
      integer        ,intent(in)    :: IXSBN(:,:)    !(5,MX_XS)
      character(100) ,intent(inout) :: HEADT1(:)     !(MX_XS)
      real(r8)       ,intent(inout) :: V1DX, V2DX
      real           ,intent(inout) :: DVDX
      integer        ,intent(inout) :: NPTSDX
      real           ,intent(inout) :: XSTEMP(:,:,:) !(6,5,MX_XS)
      real           ,intent(inout) :: XSMAX(:,:,:)  !(6,5,MX_XS)
      real           ,intent(inout) :: PDX(:,:,:)    !(6,5,MX_XS)
      real           ,intent(inout) :: RDX1(:)
      real           ,intent(inout) :: RDX2(:)
      integer        ,intent(in)    :: NFHDRF, NPHDRF


      !--- Locale variables
      !
      character(4)   ,PARAMETER :: XSTMP='TMPX'
      integer        ,PARAMETER :: LIMXX=516
      integer        ,PARAMETER :: IFILE=181
      integer        ,PARAMETER :: JFILE=182
      character(10)  ,PARAMETER :: UNBFRM='(10E10.3)'
      character(10)  ,PARAMETER :: BLKFRM='(510E10.3)'
      character(10)  ,PARAMETER :: CTORR='      TORR'
      integer        ,PARAMETER :: I_100=100

      CHARACTER     :: AMOL*8, BMOL*6, CI, BFRM*10, HEADER*100
      character     :: XSFIL1*256, XSFIL2*256, XSNUM*3
      CHARACTER(10) :: SOURCE(3)
      INTEGER       :: I,         IAFORM,  ICM
      INTEGER       :: IMFORM,    IOSTAT,  ISFORM,  ITEMP
      INTEGER       :: J,       JEOF
      INTEGER       :: NSTRT,     NXMODE
      LOGICAL       :: OP,        OPCL
      REAL          :: PATMMB,  PRES,    PTORMB
      REAL          :: SMAX,      TEMP


      !--- Common blocks are for local use only
      !
      REAL         :: PAV1,TAV1,W1,P1L,P1U,T1L,T1U,WBROA1,DVB,TBOUN1,EMISI1,FSCDI1,Y11
      REAL*8       :: V1B,V2B,SECAN1,XALT1
      INTEGER      :: NMO1,LAYER1,LSTWD1
      CHARACTER*8  :: XI1,HMOLI1,YID1
      COMMON /FXSHDR_xsecin/ XI1(10),SECAN1,PAV1,TAV1,HMOLI1(60),XALT1(4), &
     &                W1(60),P1L,P1U,T1L,T1U,WBROA1,DVB,V1B,V2B,TBOUN1, &
     &                EMISI1,FSCDI1(17),NMO1,LAYER1,Y11,YID1(10),LSTWD1

      REAL :: FILHDR(2)
      EQUIVALENCE (XI1(1), FILHDR(1))

      REAL    :: DVPX,RBX
      REAL*8  :: V1PX,V2PX
      INTEGER :: NLIMPX
      COMMON /PXSHDR_xsecin/ V1PX,V2PX,DVPX,NLIMPX,RBX(2050)

      REAL     :: PNLHDR(2),DUM(2)
      EQUIVALENCE (V1PX, PNLHDR(1))


      !%%%%%LINUX_PGI90 (-i8)%%%%%      integer*4 iostat
      !real :: RDXX1(516),RDXX2(516), &
      !        RDXA1(510),RDXA2(510), &
      !        RDXH1(500),RDXH2(500)
      !EQUIVALENCE (RDX1(5),RDXX1(1),RDXA1(1),RDXH1(1))
      !EQUIVALENCE (RDX2(5),RDXX2(1),RDXA2(1),RDXH2(1))




!     DEFINE PRESSURE CONVERSIONS
!        PTORMB = 1013. MB / 760. TORR  (TORR TO MILLIBARS)
!        PATMMB = 1013. MB / 1.0  ATM   (ATMOPHERES TO MILLIBARS)
      PTORMB = 1013./760.
      PATMMB = 1013.

      IEOF = 0
      ISFORM = IXFORM(NS,NI)
      NXMODE = NMODE
      IF (IXSBN(NS,NI).EQ.1) THEN
         IF (ABS(ISFORM).LT.100) ISFORM = ABS(ISFORM)+100
         ISFORM = -ISFORM
         NXMODE = 1
      ENDIF
      IAFORM = ABS(ISFORM)
      IMFORM = MOD(IAFORM,I_100)

!     IF NPANEL <= 0, OPEN FILE AND READ HEADER

      IF (NPANEL.LE.0) THEN
         IF (IXSBN(NS,NI).EQ.0) THEN
            XSFIL1 = trim(ioFiles%xsFilePath)//trim(XSFILE(NT1,NS,NI))  !yma XSFIL1 = XSFILE(NT1,NS,NI)
         ELSE
            WRITE (XSNUM,'(I1,I2.2)') NS,NI
            XSFIL1 = trim(ioFiles%scratchPath)//trim(XSTMP)//XSNUM
         ENDIF

         INQUIRE (FILE=trim(XSFIL1),OPENED=OP)
         INQUIRE (UNIT=IFILE,OPENED=OPCL)
         IF (.NOT.OP.AND.OPCL) CLOSE (IFILE)
         IF (.NOT.OP) THEN
            IF (ISFORM.GT.0) THEN
               OPEN (IFILE,FILE=trim(XSFIL1),STATUS='OLD',FORM='FORMATTED', IOSTAT=iostat)
               if (IOSTAT.gt.0) then
                  print*,('in oprop (a) - No file XSFIL1: '//trim(xsfil1))
                  STOP
               endif
            ELSE
               OPEN (IFILE,FILE=trim(XSFIL1),STATUS='OLD',FORM='UNFORMATTED', IOSTAT=iostat)
               if (IOSTAT.gt.0) then
                  print*,('in oprop (b) - No file XSFIL1: '//trim(xsfil1))
                  STOP
               endif
            ENDIF
         ENDIF
         REWIND IFILE
         IF (NXMODE.EQ.2) THEN
            XSFIL2 = trim(ioFiles%xsFilePath)//trim(XSFILE(NT2,NS,NI))  !yma XSFIL2 = XSFILE(NT2,NS,NI)

            INQUIRE (FILE=trim(XSFIL2),OPENED=OP)
            INQUIRE (UNIT=JFILE,OPENED=OPCL)
            IF (.NOT.OP.AND.OPCL) CLOSE (JFILE)
            IF (.NOT.OP) THEN
               IF (ISFORM.GT.0) THEN
                  OPEN (JFILE,FILE=trim(XSFIL2),STATUS='OLD',FORM='FORMATTED', IOSTAT=iostat)
                  if (IOSTAT.gt.0) then
                     print*,('in oprop - No file XSFIL2: '//trim(xsfil2))
                     STOP
                  endif
               ELSE
                  OPEN (JFILE,FILE=trim(XSFIL2),STATUS='OLD',FORM='UNFORMATTED', IOSTAT=iostat)
                  if (IOSTAT.gt.0) then
                     print*,('in oprop - No file XSFIL2: '//trim(xsfil2))
                     STOP
                  endif
               ENDIF
            ENDIF
            REWIND JFILE
         ENDIF

!     HEADER: 86 FORMAT
!
!             AMOL,V1,V2,NPTS,BMOL,PRES,ICM,ITEMP,SOURCE
!
!     HEADER: 91 FORMAT
!
!             AMOL,V1,V2,NPTS,TEMP,PRES,SMAX,SOURCE
!
!
!     IAFORM < 100, UNBLOCKED DATA (100 CHARACTERS/RECORD)
!
         IF (IAFORM.LT.100) THEN
            READ (IFILE,900,END=30) HEADER
            HEADT1(NI) = HEADER
            IF (IMFORM.EQ.86) THEN
               READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,ICM,   &
               ITEMP,SOURCE
               IF (NPANEL.EQ.0) THEN
                  XSTEMP(NT1,NS,NI) = REAL(ITEMP)+273.15
                  XSMAX(NT1,NS,NI) = 0.0
                  PDX(NT1,NS,NI) = PRES*PTORMB
               ENDIF
            ELSE
               READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,SMAX,  &
               SOURCE
               IF (NPANEL.EQ.0) THEN
                  XSTEMP(NT1,NS,NI) = TEMP
                  XSMAX(NT1,NS,NI) = SMAX
                  IF (SOURCE(3).EQ.CTORR) THEN
                     PDX(NT1,NS,NI) = PRES*PTORMB
                  ELSE
                     PDX(NT1,NS,NI) = PRES
                  ENDIF
               ENDIF
            ENDIF
            IF (NXMODE.EQ.2) THEN
               READ (JFILE,900,END=30) HEADER
               IF (IMFORM.EQ.86) THEN
                  READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,    &
                  ICM,ITEMP,SOURCE
                  IF (NPANEL.EQ.0) THEN
                     XSTEMP(NT2,NS,NI) = REAL(ITEMP)+273.15
                     XSMAX(NT2,NS,NI) = 0.0
                     PDX(NT2,NS,NI) = PRES*PTORMB
                  ENDIF
               ELSE
                  READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,    &
                  SMAX,SOURCE
                  IF (NPANEL.EQ.0) THEN
                     XSTEMP(NT2,NS,NI) = TEMP
                     XSMAX(NT2,NS,NI) = SMAX
                     IF (SOURCE(3).EQ.CTORR) THEN
                        PDX(NT2,NS,NI) = PRES*PTORMB
                     ELSE
                        PDX(NT2,NS,NI) = PRES
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ELSE

!     IAFORM > 100, BLOCKED DATA (51*100 CHARACTERS/RECORD)

            IF (ISFORM.GT.0) THEN
               READ (IFILE,915,END=30) HEADER,(RDX1(J),J=1,500) !(RDXX1(J),J=1,500)
            ELSE
               CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
               IF (IEOF.LE.0) GO TO 30
               WRITE (HEADER,'(10A8)') XI1
               CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)
               IF (IEOF.LE.0) GO TO 30
               CALL BUFIN( IFILE,IEOF,RDX1(1),NLIMPX) !(IFILE,IEOF,RDXH1(1),NLIMPX)
            ENDIF
            HEADT1(NI) = HEADER
            IF (IMFORM.EQ.86) THEN
               READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,ICM,   &
               ITEMP,SOURCE
               IF (NPANEL.EQ.0) THEN
                  XSTEMP(NT1,NS,NI) = REAL(ITEMP)+273.15
                  XSMAX(NT1,NS,NI) = 0.0
                  PDX(NT1,NS,NI) = PRES*PTORMB
               ENDIF
            ELSE
               READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,SMAX,  &
               SOURCE
               IF (NPANEL.EQ.0) THEN
                  XSTEMP(NT1,NS,NI) = TEMP
                  XSMAX(NT1,NS,NI) = SMAX
                  IF (SOURCE(3).EQ.CTORR) THEN
                     PDX(NT1,NS,NI) = PRES*PTORMB
                  ELSE
                     PDX(NT1,NS,NI) = PRES
                  ENDIF
               ENDIF
            ENDIF
            IF (NXMODE.EQ.2) THEN
               IF (ISFORM.GT.0) THEN
                  READ (JFILE,915,END=30) HEADER,(RDX2(J),J=1,500) !(RDXX2(J),J=1,500)
               ELSE
                  CALL BUFIN (JFILE,JEOF,FILHDR(1),NFHDRF)
                  IF (JEOF.LE.0) GO TO 30
                  WRITE (HEADER,'(10A8)') XI1
                  CALL BUFIN (JFILE,JEOF,PNLHDR(1),NPHDRF)
                  IF (JEOF.LE.0) GO TO 30
                  CALL BUFIN( JFILE,JEOF,RDX2(1),NLIMPX) !(JFILE,JEOF,RDXH2(1),NLIMPX)
               ENDIF
               IF (IMFORM.EQ.86) THEN
                  READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,    &
                  ICM,ITEMP,SOURCE
                  IF (NPANEL.EQ.0) THEN
                     XSTEMP(NT2,NS,NI) = REAL(ITEMP)+273.15
                     XSMAX(NT2,NS,NI) = 0.0
                     PDX(NT2,NS,NI) = PRES*PTORMB
                  ENDIF
               ELSE
                  READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,    &
                  SMAX,SOURCE
                  IF (NPANEL.EQ.0) THEN
                     XSTEMP(NT2,NS,NI) = TEMP
                     XSMAX(NT2,NS,NI) = SMAX
                     IF (SOURCE(3).EQ.CTORR) THEN
                        PDX(NT2,NS,NI) = PRES*PTORMB
                     ELSE
                        PDX(NT2,NS,NI) = PRES
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         DVDX = (V2DX-V1DX)/ REAL(NPTSDX-1)

!     FOR NPANEL = -1, SKIP REQUIRED NUMBER OF RECORDS

         IF (NPANEL.EQ.-1) THEN
            NSTRT = 1
            IF (IAFORM.GT.100) THEN
               NSKIP = NSKIP/51
               NSTRT = 2
               IF (NSKIP.EQ.0) RETURN
            ENDIF
!
            DO 10 I = NSTRT, NSKIP
               IF (ISFORM.GT.0) THEN
                  READ (IFILE,920,END=30) CI
                  IF (NXMODE.EQ.2) READ (JFILE,920,END=30) CI
               ELSE
                  CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)
                  IF (IEOF.LE.0) GO TO 30
                  CALL BUFIN (IFILE,IEOF,DUM(1),1)
                  IF (NXMODE.EQ.2) THEN
                     CALL BUFIN (JFILE,JEOF,PNLHDR(1),NPHDRF)
                     IF (JEOF.LE.0) GO TO 30
                     CALL BUFIN (JFILE,JEOF,DUM(1),1)
                  ENDIF
               ENDIF
   10       CONTINUE
         ENDIF
      ENDIF

!     FOR ABS(NPANEL) > 0, READ IN MORE DATA

      IF (ABS(NPANEL).GT.0) THEN
         BFRM = UNBFRM
         IF (IAFORM.GT.100) BFRM = BLKFRM
         IF (ISFORM.GT.0) THEN
            READ (IFILE,BFRM,END=30) (RDX1(J),J=1,NMAX) !(RDXX1(J),J=1,NMAX)
            IF (NXMODE.EQ.2) READ (JFILE,BFRM,END=30) (RDX2(J),J=1,NMAX) !(RDXX2(J),J=1,NMAX)
         ELSE
            CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)
            IF (IEOF.LE.0) GO TO 30
            CALL BUFIN( IFILE,IEOF,RDX1(1),NLIMPX) !(IFILE,IEOF,RDXA1(1),NLIMPX)
            IF (NXMODE.EQ.2) THEN
               CALL BUFIN (JFILE,JEOF,PNLHDR(1),NPHDRF)
               IF (JEOF.LE.0) GO TO 30
               CALL BUFIN( JFILE,JEOF,RDX2(1),NLIMPX) !(JFILE,JEOF,RDXA2(1),NLIMPX)
            ENDIF
         ENDIF
      ENDIF

      DO 20 I = NMAX+1, LIMXX
         RDX1(I) = 0.0 !RDXX1(I) = 0.0
         IF (NXMODE.EQ.2) RDX2(I) = 0.0 !RDXX2(I) = 0.0
   20 END DO

      RETURN

   30 IEOF = 1

      DO 40 I = NMAX+1, LIMXX
         RDX1(I) = 0.0 !RDXX1(I) = 0.0
         IF (NXMODE.EQ.2) RDX2(I) = 0.0  !RDXX2(I) = 0.0
   40 END DO

      RETURN

  900 FORMAT (A100)
  905 FORMAT (A8,2F10.4,I10,1X,A6,F4.2,5X,I4,3X,I5,3A10)
  910 FORMAT (A10,2F10.4,I10,3G10.3,3A10)
  915 FORMAT (A100,50(10E10.3))
  920 FORMAT (A1)

      END SUBROUTINE

!-----------------------------------------------------------------------
!     THIS SUBROUTINE DETERMINES THE CORRECT MODE
!     AND BRACKETS THE LAYER TEMPERATURE
!-----------------------------------------------------------------------
      SUBROUTINE XSNTMP( NI,NS,NT1,NT2,NMODE, &
                         TAVE,XSTEMP,NTEMPF,V1FX,WXM,XSMAX,DPTMIN,JRAD )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, RADCN2
      USE Module_ODLAY_CommonSub  ,ONLY: RADFN

      IMPLICIT NONE !REAL*8           (V)

      integer  ,intent(in)  :: NS, NI
      integer  ,intent(out) :: NMODE
      integer  ,intent(out) :: NT1, NT2
      !
      real     ,intent(in)  :: TAVE
      real     ,intent(in)  :: XSTEMP(:,:,:) !(6,5,MX_XS)
      integer  ,intent(in)  :: NTEMPF(:,:)   !(5,MX_XS)
      real(r8) ,intent(in)  :: V1FX(:,:)     !(5,MX_XS)
      real     ,intent(in)  :: WXM(:)        !(MX_XS)
      real     ,intent(in)  :: XSMAX(:,:,:)  !(6,5,MX_XS)
      real     ,intent(in)  :: DPTMIN
      integer  ,intent(in)  :: JRAD

      !--- Local variables
      INTEGER  :: I,       IDPTMN
      REAL     :: RADVI1,  RADVI2
      REAL(r8) :: VI
      REAL     :: WXM1,    WXM2,    XKT1,   XKT2



      NT1 = 0
      NT2 = 0
      IDPTMN = 0 !yma

      IF (NTEMPF(NS,NI).LE.1) THEN
         NMODE = 1
         NT1 = 1
         NT2 = 1
      ELSE
         NMODE = 2
         DO 10 I = 2, NTEMPF(NS,NI)
            IF (TAVE.LT.XSTEMP(I,NS,NI)) THEN
               NT1 = I-1
               NT2 = I
               GO TO 20
            ENDIF
   10    CONTINUE
      ENDIF

   20 IF (NT1.EQ.0) THEN
         NT2 = NTEMPF(NS,NI)
         NT1 = NT2-1
      ENDIF

!     CHECK VERSUS DPTMIN

      IF (XSMAX(NT1,NS,NI).NE.0.0) THEN
         WXM1 = WXM(NI)*XSMAX(NT1,NS,NI)
         IF (WXM1.LT.DPTMIN) IDPTMN = IDPTMN+1
         WXM2 = 0.0
         IF (NMODE.EQ.2) THEN
            WXM2 = WXM(NI)*XSMAX(NT2,NS,NI)
            IF (WXM2.LT.DPTMIN) IDPTMN = IDPTMN+1
         ENDIF
         IF (JRAD.EQ.0) THEN
            XKT1 = XSTEMP(NT1,NS,NI)/RADCN2
            IF (NMODE.EQ.2) XKT2 = XSTEMP(NT2,NS,NI)/RADCN2
            VI = V1FX(NS,NI)

            RADVI1 = RADFN(VI,XKT1)
            IF (NMODE.EQ.2) RADVI2 = RADFN(VI,XKT2)
            WXM1 = WXM(NI)*XSMAX(NT1,NS,NI)/RADVI1
            IF (NMODE.EQ.2) WXM2 = WXM(NI)*XSMAX(NT2,NS,NI)/RADVI2
         ENDIF

!     DETERMINE IDPTMN --- IF IDPTMN = NMODE  ==> BELOW THRESHOLD


         IDPTMN = 0
         IF (IDPTMN.EQ.NMODE) NMODE = 0
      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     THIS SUBROUTINE PERFORMS A TEMPERATURE DEPENDENT CONVOLUTION
!     ON THE CROSS-SECTIONS PRODUCING A BINARY INTERMEDIATE FILE
!-----------------------------------------------------------------------
      SUBROUTINE XSBINF( NI,NS,NT1,NT2,NMODE, &
                         IXSBN,XSFILE,IXFORM,V1FX,V2FX,DVFX,NPTSFX,&
                         XDOPLR,PF,TF,argTAVE,XSTEMP,PDX,HEADT1,NFHDRF,NPHDRF,DVXPR )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8, MX_XS, ONEPL
      USE Module_Config      ,ONLY: ioFiles
      IMPLICIT NONE !REAL*8           (V)

      integer       ,intent(in)    :: NS, NI
      integer       ,intent(in)    :: NMODE
      integer       ,intent(in)    :: NT1, NT2
      !
      integer       ,intent(inout) :: IXSBN(:,:)    !(5,MX_XS)
      character(10) ,intent(in)    :: XSFILE(:,:,:) !(6,5,mx_xs)
      integer       ,intent(in)    :: IXFORM(:,:)   !(5,MX_XS)
      real(r8)      ,intent(in)    :: V1FX(:,:)     !(5,MX_XS)
      real(r8)      ,intent(in)    :: V2FX(:,:)     !(5,MX_XS)
      real          ,intent(in)    :: DVFX(:,:)     !(5,MX_XS)
      integer       ,intent(in)    :: NPTSFX(:,:)   !(5,MX_XS)
      real          ,intent(in)    :: XDOPLR(:,:)   !(5,MX_XS)
      real          ,intent(in)    :: PF
      real          ,intent(in)    :: TF
      real          ,intent(in)    :: argTAVE
      real          ,intent(inout) :: XSTEMP(:,:,:) !(6,5,MX_XS)
      real          ,intent(inout) :: PDX(:,:,:)    !(6,5,MX_XS)
      character*100 ,intent(inout) :: HEADT1(:)     !(MX_XS)
      integer       ,intent(in)    :: NFHDRF, NPHDRF
      real          ,intent(inout) :: DVXPR(:,:)


      !--- Local variables
      !
      real    ,PARAMETER :: HWJ=16.
      real    ,PARAMETER :: DXJ=0.02
      integer ,PARAMETER :: NJ=801
      integer ,PARAMETER :: NJMX=851
      real    ,PARAMETER :: SMPLJ=4.
      real    ,PARAMETER :: XSCAL=0.

!     STANDARD PRESSURE AND TEMPERATURE
      real ,PARAMETER :: P0=1013.
      real ,PARAMETER :: T0=273.15

!     ASSUMED MEAN HALFWIDTH AT P0 IS 0.10
!     DOPPLER VALUES ARE INITIALIZED AT T296.
      real ,PARAMETER :: HWHM0=0.10
      real ,PARAMETER :: T296=296.0

      integer   :: NF
      integer   :: NFMAX
      real      :: DXF
      real      :: HWF
      integer   :: JABS
      integer   :: NSHIFT
      real      :: DVO
      real(r8)  :: V1, V2
      integer   :: JEMIT
      real      :: HWHM
      real      :: SAMPLE
      real      :: XF(851)
      integer   :: NLIMX
      real(r8)  :: V1DX,V2DX
      real      :: DVDX
      integer   :: NPTSDX
      real      :: XSMAX(6,5,MX_XS)
      real      :: RDX1(520)
      real      :: RDX2(520)

      character*4 :: XSTMP='TMPX'
      integer     :: IFILEO=193, JFILEO=191

      CHARACTER    :: XSFIL*256, tmpxbinFile*256
      CHARACTER*3  :: XSNUM
      INTEGER      :: IEOF,  IMAX,   IOTPAN
      INTEGER      :: JI,      JJ,     LPMAX
      INTEGER      :: NEOF,    NFMXSV, NFSV
      INTEGER      :: NMAX,   NP
      INTEGER      :: NPAN,    NPANEL,  NSKIP
      LOGICAL      :: OP,      OPCL
      REAL         :: DXFSV,     FAC1,   FAC2
      REAL         :: HDOPLR,  HWFSV,  HWHMD
      REAL         :: HWHMF,   HWHMSC
      REAL         :: PD


      !--- Local common blocks
      REAL         :: PAV1,TAV1,W1,P1L,P1U,T1L,T1U,WBROA1,DVB,TBOUN1,EMISI1,FSCDI1,Y11
      REAL*8       :: V1B,V2B,SECAN1,XALT1
      INTEGER      :: NMO1,LAYER1,LSTWD1
      CHARACTER*8  :: XI1,HMOLI1,YID1
      REAL         :: PAVE,TAVE,WK,PZL,PZU,TZL,TZU,WBROAD,DV0,TBOUND,EMISIV,FSCDID,YI1
      REAL*8       :: V10,V20,SECANT,XALTZ
      INTEGER      :: NMOL,LAYER,LSTWDF
      CHARACTER*8  :: XID,HMOLID,YID
      COMMON /FXSHDR_xsbinf/ XI1(10),SECAN1,PAV1,TAV1,HMOLI1(60),XALT1(4),     &
     &                W1(60),P1L,P1U,T1L,T1U,WBROA1,DVB,V1B,V2B,TBOUN1, &
     &                EMISI1,FSCDI1(17),NMO1,LAYER1,Y11,YID1(10),LSTWD1
      COMMON /FILHDR_xsbinf/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV0,V10,V20,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF

      REAL    :: FILHDR(2),FILHDS(2)
      real    :: XSCID,XHWHM
      integer :: IEMIT,ISCHDR
      EQUIVALENCE (IEMIT,  FSCDI1(5) )
      EQUIVALENCE (ISCHDR, FSCDI1(6) )
      EQUIVALENCE (XSCID,  FSCDI1(12))
      EQUIVALENCE (XHWHM,  FSCDI1(13))
      EQUIVALENCE (XI1(1), FILHDR(1) )
      EQUIVALENCE (XID(1), FILHDS(1) )

      REAL     :: DVPX,RBX
      REAL*8   :: V1PX, V2PX
      INTEGER  :: NLIMPX
      COMMON /PXSHDR_xsbinf/ V1PX,V2PX,DVPX,NLIMPX,RBX(2050)

      REAL :: PNLHDR(2)
      EQUIVALENCE (V1PX,PNLHDR(1))



      TAVE=argTAVE !yma160616; added to make /FILHDR/ a local block


!     INITIALIZE IXSBN AND TEMPERATURE RATIO

      IXSBN(NS,NI) = 0

      FAC1 = 1.0
      FAC2 = 0.0
      IF (NMODE.EQ.2) THEN
         FAC2 = (TAVE-XSTEMP(NT1,NS,NI))/ (XSTEMP(NT2,NS,NI)-XSTEMP(NT1,&
         NS,NI))
         FAC1 = 1.0-FAC2
      ENDIF
      PD = PDX(NT1,NS,NI)*FAC1+PDX(NT2,NS,NI)*FAC2

!     NOTE THAT AT THIS POINT, THE CROSS-SECTIONS HAVE BEEN
!     LINEARLY INTERPOLATED IN TEMPERATURE (HENCE TD=TF),
!     AND THE CONVOLUTION WILL BE DONE ONLY FOR PRESSURE

      HWHMF = HWHM0*(PF/P0)*(T0/TF)
      HWHMD = HWHM0*(PD/P0)*(T0/TF)

!     SET MINIMUM HALF-WIDTH TO DOPPLER

      HDOPLR=XDOPLR(NS,NI)*SQRT(TF/T296)
      HWHMD=MAX(HDOPLR,HWHMD)
      HWHMSC = HWHMF-HWHMD
      IF (HWHMSC/HWHMD.LT.0.1) GO TO 30
      HWHM = HWHMSC

!     BOUND AT THIS POINT IS THE WAVENUMBER VALUE
!     OF HALF THE SCANNING FUNCTION

      DVO = HWHMF/2.0
      IF (HWHMSC/2.0.LT.DVFX(NS,NI)) GO TO 30
      SAMPLE = HWHM/DVO
      XHWHM = HWHM

!     OPEN FILE AND SET FILHDR

      tmpxbinFile = trim(ioFiles%scratchPath)//'TMPXBIN'
      INQUIRE (FILE=trim(tmpxbinFile),OPENED=OP)
      INQUIRE (UNIT=IFILEO,OPENED=OPCL)
      IF (.NOT.OP) THEN
         IF (OPCL) CLOSE (IFILEO)
         OPEN (IFILEO,FILE=trim(tmpxbinFile),STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND IFILEO
         CALL BUFOUT (IFILEO,FILHDS(1),NFHDRF)
         REWIND IFILEO
         CALL BUFIN (IFILEO,IEOF,FILHDR(1),NFHDRF)
         READ (HEADT1(NI),'(10A8)') XI1
      ENDIF
      REWIND IFILEO

      NLIMX = 510
      IOTPAN = 1
      LPMAX = 0

      NPAN = (NPTSFX(NS,NI)+9)/NLIMX+1
      V1PX = V1FX(NS,NI)
      NPANEL = -1
      NSKIP = 0

      DO 20 NP = 1, NPAN
         IF (NP.NE.1) NPANEL = 1
         NMAX = 510
         IF (NP.EQ.1) NMAX = 500
         V2PX = V1PX+ REAL(LPMAX+NMAX-1)*DVFX(NS,NI)
         V2PX = MIN(V2PX,V2FX(NS,NI))
         NMAX = ((V2PX-V1PX)/DVFX(NS,NI)+ONEPL)-LPMAX

         IMAX = MIN(NMAX,NLIMX)

         !CALL CPUTIM (TIME0)
         CALL XSECIN( NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,IMAX,NEOF, &
                      XSFILE,IXFORM,IXSBN,HEADT1,V1DX,V2DX,DVDX,NPTSDX, &
                      XSTEMP,XSMAX,PDX,RDX1(5:),RDX2(5:),NFHDRF,NPHDRF )
         !CALL CPUTIM (TIME)
         !TXSRDF = TXSRDF+TIME-TIME0
         !TXSCNV = TXSCNV-TIME+TIME0

         DO 10 JI = 1, NMAX
            JJ = JI+4
            RBX(JI+LPMAX) = FAC1*RDX1(JJ)
            IF (NMODE.EQ.2) RBX(JI+LPMAX) = RBX(JI+LPMAX)+FAC2*RDX2(JJ)
   10    CONTINUE

         LPMAX = LPMAX+NMAX
         IOTPAN = IOTPAN+1

         IF (NP.EQ.1) THEN
            V1B = V1FX(NS,NI)
            V2B = V2FX(NS,NI)
            DVB = DVFX(NS,NI)
            XSCID = -99
            ISCHDR = 0
            IEMIT = 0
            CALL BUFOUT (IFILEO,FILHDR(1),NFHDRF)
         ENDIF

         IF (IOTPAN.EQ.5.OR.NP.EQ.NPAN) THEN
            IOTPAN = 1
            DVPX = DVFX(NS,NI)
            NLIMPX = LPMAX
            CALL BUFOUT (IFILEO,PNLHDR(1),NPHDRF)
            CALL BUFOUT (IFILEO,RBX(1),NLIMPX)
            LPMAX = 0
            V1PX = V2PX+DVFX(NS,NI)
         ENDIF

   20 END DO

      NSHIFT = 0

      V1 = V1FX(NS,NI)
      V2 = V2FX(NS,NI)

      JEMIT = 0
      JABS = 0

      HWF = HWJ
      DXF = DXJ
      NF = NJ
      NFMAX = NJMX
!CP   XSCALE = XSCAL
      CALL SLRENZ( XF, NF,NFMAX,DXF )

      WRITE (XSNUM,'(I1,I2.2)') NS,NI
      XSFIL = trim(ioFiles%scratchPath)//trim(XSTMP)//XSNUM

      INQUIRE (FILE=trim(XSFIL),OPENED=OP)
      INQUIRE (UNIT=JFILEO,OPENED=OPCL)
      IF (.NOT.OP.AND.OPCL) CLOSE (JFILEO)
      IF (.NOT.OP) OPEN (JFILEO,FILE=trim(XSFIL),STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND JFILEO

      CALL XSCNVN( IFILEO,JFILEO,NS,NI, &
                   V1,V2,DVO,NSHIFT,HWF,DXF,HWHM,SAMPLE,XF, &
                   JEMIT,JABS,NFHDRF,NPHDRF,DVXPR )

      CLOSE (IFILEO)
      CLOSE (JFILEO)

      IXSBN(NS,NI) = 1

   30 RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!
!     DRIVER FOR CONVOLVING SPECTRUM WITH INSTRUMENTAL SCANNING FUNCTIO
!
!-----------------------------------------------------------------------
      SUBROUTINE XSCNVN( IFILE,JFILE,NS,NI, &
                         V1,V2,DVO,NSHIFT,HWF,DXF,HWHM,SAMPLE,XF, &
                         JEMIT,JABS,NFHDRF,NPHDRF,DVXPR )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8
      USE Module_ScanFilter     ,ONLY: CONVSC, SHRKSC
      USE Module_FileIO      ,ONLY: SKIPFL !yma 151225

      IMPLICIT NONE !REAL*8           (V)

      integer   ,intent(in)    :: IFILE
      integer   ,intent(in)    :: JFILE
      integer   ,intent(in)    :: NS,NI
      !
      real(r8)  ,intent(in)    :: V1,V2
      real      ,intent(in)    :: DVO
      integer   ,intent(in)    :: NSHIFT
      real      ,intent(in)    :: HWF
      real      ,intent(in)    :: DXF
      real      ,intent(in)    :: HWHM
      real      ,intent(in)    :: SAMPLE
      real      ,intent(in)    :: XF(:) !(851)
      integer   ,intent(inout) :: JEMIT
      integer   ,intent(inout) :: JABS
      integer   ,intent(in)    :: NFHDRF
      integer   ,intent(in)    :: NPHDRF
      real      ,intent(inout) :: DVXPR(:,:)


      !---Local variables
      !
      real(r8)  :: V1I,V2I
      real      :: DVI
      integer   :: NNI
      real      :: S(2050),R1(1025),SS(200)
      real      :: SUMIN
      real(r8)  :: VFT,VBOT,VTOP
      integer   :: NREN
      integer   :: JVAR
      integer   :: IPRT
      integer   :: ILO,IHI
      integer   :: NLO,NHI
      integer   :: NLIMF
      integer   :: MAXF
      integer   :: IEOFSC
      integer   :: IPANEL
      integer   :: IDATA
      integer   :: ISTOP


      INTEGER  :: I,       IDABST, IEOF,    IEOFT
      INTEGER  :: IFILST,  INIT,   ISCAN
      INTEGER  :: IUNIT,   JTREM,  JUNIT,   NBOUND
      INTEGER  :: NIFILS,  NXPAN
      REAL     :: BOUND,   DVOSAV, DVSAV,   SMAX
      REAL     :: SMIN,    SUMOUT, SUMR(4)


      !--- Common block /SCNHDR/ for internal use only. Not for passing variables between subroutines.
      !
      REAL         :: PAVE,TAVE,WK,PZL,PZU,TZL,TZU,WBROAD,DV,TBOUND,EMISIV,FSCDID,YI1
      REAL*8       :: V1C,V2C,SECANT,XALTZ
      INTEGER      :: NMOL,LAYER,LSTWDF
      CHARACTER*8  :: XID,HMOLID,YID
      COMMON /SCNHDR_xscnvn/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV,V1C,V2C,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF

      real :: FILHDR(2),XHWHM
      integer :: IEMIT,ISCHDR,IDABS
      EQUIVALENCE (FILHDR(1)  ,XID(1))
      EQUIVALENCE (FSCDID(5)  ,IEMIT )
      EQUIVALENCE (FSCDID(6)  ,ISCHDR)
      EQUIVALENCE (FSCDID(13) ,XHWHM )
      EQUIVALENCE (FSCDID(14) ,IDABS )



!     IUNIT INPUT FILE
!     JUNIT OUTPUT FILE
      IUNIT = IFILE
      JUNIT = JFILE
      NREN = 0
      IDABS = 0
      JVAR = 0
      IPRT = 0

      IF (JEMIT.LT.0) THEN
         JABS = 1
         JEMIT = 0
         IDABS = -1
      ENDIF
      IDABST = IDABS
      IFILST = 1
      NIFILS = 9999

      SUMOUT = 0.
      SMIN = 999999.
      SMAX = -99999.
      DVOSAV = 0.
      SUMR(1) = SUMOUT
      SUMR(2) = SMIN
      SUMR(3) = SMAX
      SUMR(4) = DVOSAV

      REWIND IUNIT
      CALL BUFIN (IUNIT,IEOF,FILHDR(1),NFHDRF)
      IF (IEOF.EQ.0) GO TO 50

      DVSAV = DV
      IDABS = IDABST

      ISCAN = ISCHDR
      JTREM = 3

!     JTREM=3   SCANFN CONVOLVED WITH OPTICAL DEPTH

      DVI = DV

!     BOUND AT THIS POINT IS THE WAVENUMBER VALUE
!     OF HALF THE SCANNING FUNCTION

      BOUND = HWF*HWHM
      DV = DVO
      V1C = V1
      V2C = V2
      XHWHM = HWHM
      IEMIT = 0
      CALL BUFOUT (JUNIT,FILHDR(1),NFHDRF)
      NBOUND = (2.*HWF)*SAMPLE+0.01

!     BOUND AT THIS POINT IS THE WAVENUMBER VALUE
!     OF THE FULL SCANNING FUNCTION

      BOUND =  REAL(NBOUND)*DVO/2.

      NXPAN = 500
      NLO = NBOUND+1
      NLIMF = NLO+NXPAN-NSHIFT
      NHI = NLIMF+NSHIFT-1
      MAXF = NLIMF+2*NBOUND

      !TIMRDF = 0.
      !TIMCNV = 0.
      !TIMPNL = 0.
      IEOFSC = 1
      SUMIN = 0.
      DO 10 I = 1, MAXF
         R1(I) = 0.
   10 END DO
      INIT = 0
      IDATA = -1
      VFT = V1-2.*BOUND
      VBOT = V1-BOUND
      VTOP = V2+BOUND

   20 continue  !---start loop over panels

      !CALL CPUTIM (TIME0)
      IF (IEOFSC.LE.0) GO TO 40
      CALL RDSCAN( S,JTREM,IUNIT,ISCAN,IPRT, &
                   V1I,V2I,DVI,NNI,ILO,IHI,VBOT,VTOP,NREN,&
                   IDATA,IEOFSC,JABS,NPHDRF )

      !CALL CPUTIM (TIME)
      !TIMRDF = TIMRDF+TIME-TIME0

      IF (IEOFSC.LE.0) GO TO 40
      CALL SHRKSC( INIT,HWHM, &
                   S,SS,V1I,V2I,DVI,NNI,VBOT,VTOP,ILO,IHI,NREN )

!     SHRKSC MAY SHRINK (COMPRESS) THE DATA;
!     DVI IS MODIFIED ACCORDINGLY

   30 CONTINUE
      CALL CONVSC( S,HWHM,R1,XF, &
                   HWF,DXF,V1I,DVI,ILO,IHI,MAXF,VFT,V1,DVO,&
                   JVAR,IDATA,IPANEL )

      IF (IPANEL.EQ.0) GO TO 20

   40 CALL PNLCNV( R1,JUNIT,SUMR,NS,NI, &
                   NLO,NHI,NLIMF,MAXF,NSHIFT,VFT,V2,DVO,ISTOP,NPHDRF,DVXPR )
      IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 30
      IF (ISTOP.NE.1) GO TO 40
      !CALL CPUTIM (TIME)

      SUMIN = SUMIN*DVSAV

      IF (IEOFSC.EQ.1) CALL SKIPFL (1,IUNIT,IEOFSC)

!yma      IEOFT = IEOFT+1

      SUMOUT = SUMR(1)
      SMIN = SUMR(2)
      SMAX = SUMR(3)
      DVOSAV = SUMR(4)

      SUMOUT = SUMOUT*DVOSAV

   50 RETURN

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE PNLCNV OUTPUTS THE RESULTS OF THE CONVOLUTION
!     TO FILE JFILE
!-----------------------------------------------------------------------
      SUBROUTINE PNLCNV( R1,JFILE,SUMR,NS,NI, &
                         NLO,NHI,NLIMF,MAXF,NSHIFT,VFT,V2,DVO,ISTOP,NPHDRF,DVXPR )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE !REAL*8           (V)

      real     ,intent(inout) :: R1(:)
      integer  ,intent(in)    :: JFILE
      real     ,intent(inout) :: SUMR(:)
      integer  ,intent(in)    :: NS,NI
      !
      integer  ,intent(inout) :: NLO,NHI
      integer  ,intent(inout) :: NLIMF
      integer  ,intent(in)    :: MAXF
      integer  ,intent(in)    :: NSHIFT
      real(r8) ,intent(inout) :: VFT
      real(r8) ,intent(in)    :: V2
      real     ,intent(in)    :: DVO
      integer  ,intent(out)   :: ISTOP
      integer  ,intent(in)    :: NPHDRF
      real     ,intent(inout) :: DVXPR(:,:) !(5,MX_XS)


      !--- Local variables
      !
      INTEGER :: I,      J,    JF
      INTEGER :: NLIMHI, NNHI
      REAL    :: SMAX,   SMIN, SUMOUT

      !--- Common block /SPANEL/ is for internal use only. not used to share variable with any external subroutines.
      real(r8) :: V1P,V2P
      integer  :: NLIM
      real     :: DV
      COMMON /SPANEL_pnlcnv/ V1P,V2P,DV,NLIM

      real :: PNLHDR(2)
      EQUIVALENCE (PNLHDR(1),V1P)




      !CALL CPUTIM (TIME0)

      SUMOUT = SUMR(1)
      SMIN = SUMR(2)
      SMAX = SUMR(3)
      DV = DVO
      ISTOP = 0
      NNHI = (V2-VFT)/DV+1.5
      IF (NHI.GE.NNHI) THEN
         ISTOP = 1
         NHI = NNHI
      ENDIF
      NLIM = NHI-NLO+1
      V1P = VFT+ REAL(NLO-1)*DV
      V2P = VFT+ REAL(NHI-1)*DV

!     V1P IS FIRST FREQ OF PANEL
!     V2P IS LAST  FREQ OF PANEL
!
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)
      CALL BUFOUT (JFILE,R1(NLO),NLIM)

      VFT = VFT+ REAL(NLIMF-1)*DV
      DVXPR(NS,NI) = DV
      NLIMHI = NLIM+NLO-1
      DO 10 I = NLO, NLIMHI
         SMIN = MIN(SMIN,R1(I))
         SMAX = MAX(SMAX,R1(I))
         SUMOUT = SUMOUT+R1(I)
   10 END DO
      IF (ISTOP.EQ.1) GO TO 40
      JF = 1
      DO 20 J = NLIMF, MAXF
         R1(JF) = R1(J)
         JF = JF+1
   20 END DO
      DO 30 J = JF, MAXF
         R1(J) = 0.
   30 END DO
      NLIMF = 511
      NLO = NSHIFT+1
      NHI = NLIMF+NSHIFT-1
   40 SUMR(1) = SUMOUT
      SUMR(2) = SMIN
      SUMR(3) = SMAX
      SUMR(4) = DVO
      !CALL CPUTIM (TIME)
      !TIMPNL = TIMPNL+TIME-TIME0

      RETURN

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE SLRENZ SETS UP THE LORENZ SCANNING FUNCTION
!-----------------------------------------------------------------------
      SUBROUTINE SLRENZ( XF, NF,NFMAX,DXF )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: pi
      IMPLICIT NONE

      real    ,intent(out) :: XF(:)
      integer ,intent(in)  :: NF
      integer ,intent(in)  :: NFMAX
      real    ,intent(in)  :: DXF

      integer :: I
      real    :: SUM, X, XNORM


      XNORM = 1.0/PI
      DO 10 I = 1, NFMAX
         XF(I) = 0.
   10 END DO
      XF(1) = XNORM
      SUM = XF(1)
      DO 20 I = 2, NF
         X = REAL(I-1)*DXF
         XF(I) = XNORM*(1./(1.+X**2))
         SUM = SUM+2.*XF(I)
   20 END DO
      SUM = SUM*DXF

!     RENORMALIZE

      XNORM = 1.0/SUM
      DO 30 I = 1, NF
         XF(I) = XNORM*XF(I)
   30 END DO

!     WRITE(IPR,900) NF,DXF,SUM

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE RDSCAN INPUTS PANELS FROM IFILE RESULTING
!     FROM THE LBLRTM CALCULATION FOR CONVOLUTION
!     WITH THE SELECTED SCANNING FUNCTION
!-----------------------------------------------------------------------
      SUBROUTINE RDSCAN( S,JTREM,IFILE,ISCAN,IPRT, &
                         VMIN,VMAX,DVI,NNI,ILO,IHI,VBOT,VTOP,NREN,&
                         IDATA,IEOFSC,JABS,NPHDRF )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8, ARGMIN, EXPMIN
      USE Module_Config       ,ONLY: IPR
      IMPLICIT NONE !REAL*8          (V)

      real     ,intent(inout) :: S(:)
      integer  ,intent(in)    :: JTREM
      integer  ,intent(in)    :: IFILE
      integer  ,intent(in)    :: ISCAN
      integer  ,intent(in)    :: IPRT
      !
      real(r8) ,intent(out)   :: VMIN,VMAX
      real     ,intent(out)   :: DVI
      integer  ,intent(out)   :: NNI
      integer  ,intent(out)   :: ILO, IHI
      real(r8) ,intent(in)    :: VBOT,VTOP
      integer  ,intent(in)    :: NREN
      integer  ,intent(inout) :: IDATA
      integer  ,intent(inout) :: IEOFSC
      integer  ,intent(in)    :: JABS
      integer  ,intent(in)    :: NPHDRF


      integer ,PARAMETER :: I_1000=1000

      INTEGER :: I,      IDUM1,    IDUM2,  ISCANT
      INTEGER :: NLOW,   NNB
      REAL    :: DIF,    SI, DUMMY(2)


      !--- Common block /RSCAN_RDSCAN/ for internal use only. To buffin the panel header.
      real(r8) :: cmVMIN,cmVMAX
      real     :: cmDVI
      integer  :: cmNNI
      COMMON /RSCAN_rdscan/ cmVMIN,cmVMAX,cmDVI,cmNNI

      real :: PNLHDR(2)
      EQUIVALENCE (PNLHDR(1),cmVMIN)




!PRT  WRITE(IPR,900) VBOT,VTOP

      IDUM1 = 0
      IDUM2 = 0
      ISCANT = MOD(ISCAN,I_1000)
      IF(JTREM.EQ.0.AND.ISCANT.GE.1) GO TO 60
      IF (ISCAN.LT.1) THEN
         IF (JTREM.EQ.1) IDUM1 = 1
         IF (JTREM.EQ.2) IDUM2 = 1
      ENDIF

   10 CALL BUFIN (IFILE,IEOFSC,PNLHDR(1),NPHDRF)
         IF (IEOFSC.LE.0) GO TO 50
!------------->>>yma 160613
         VMIN = cmVMIN
         VMAX = cmVMAX
         DVI  = cmDVI
         NNI  = cmNNI
!------------->>><<<yma 160613

         NLOW = NREN+1
         IF (NREN.LE.0) NLOW = 1
         VMIN = VMIN-(NLOW-1)*DVI
         NNB = NNI
         NNI = NNI+NLOW-1
         IF ((IDATA.EQ.-1).AND.(VMIN.GT.VBOT).AND.(IPRT.EQ.1))             &
     &        WRITE (IPR,905)
         IDATA = 0
         IF (VMAX.GE.VBOT) GO TO 20

         IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
         CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
         IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
      GO TO 10

   20 IF (JTREM.EQ.0 .OR. JTREM.EQ.4 ) THEN
         CALL BUFIN (IFILE,IEOFSC,S(NLOW),NNB)
         DO 30 I = NLOW, NNI
            SI = S(I)
            S(I) = 1.
            IF (SI.GT.1.0E-04) THEN
               IF (SI.LT.ARGMIN) THEN
                  S(I) = EXP(-SI)
               ELSE
                  S(I) = EXPMIN
               ENDIF
            ELSE
               S(I) = 1.-SI
            ENDIF
   30    CONTINUE
      ELSE

         IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
         CALL BUFIN (IFILE,IEOFSC,S(NLOW),NNB)
         IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
      ENDIF

!PRT  WRITE(IPR,910) VMIN,VMAX,DVI,NLOW,NNI

      IF (JABS.NE.0) THEN
         DO 40 I = NLOW, NNI
            S(I) = 1.-S(I)
   40    CONTINUE
      ENDIF

      ILO = 1
      IHI = NNI
      DIF = (VMIN-VBOT)/DVI
      IF (DIF.LT.0.) ILO = -DIF+1.5
      IF (VMAX.LE.VTOP) RETURN
      IHI = (VTOP-VMIN)/DVI+1.5
      IDATA = 1
      RETURN

   50 IF (IPRT.EQ.1) WRITE (IPR,915)
      RETURN

   60 WRITE(IPR,920) JTREM,ISCAN
      RETURN

  900 FORMAT ('0',/,'0   READING SPECTRUM, VBOT =',F10.3,', VTOP =',    &
     &        F10.3)
  905 FORMAT ('0 ********** FIRST VALUE USED ON IFILE; CHECK IFILE ')
  910 FORMAT (10X,'VMIN =',F10.3,',  VMAX =',F10.3,',  DVI=',F7.5,',    &
     &        NLOW=',I4,',  NNI=',I4)
  915 FORMAT ('0 ********** END OF FILE ENCOUNTERED; CHECK IFILE ')
  920 FORMAT(' ERROR IN INPUT',/,'  JTREM =',I2,'  ISCAN=',I5)

      END SUBROUTINE


END MODULE
