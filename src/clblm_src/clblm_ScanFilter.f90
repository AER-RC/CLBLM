!
! CREATION HISTORY:
!       Modified from LBLRTM v12.9
!       Yingtao Ma, AER@NOAA/NESDIS/STAR
!       yma@aer.com; yingtao.ma@noaa.gov
!


MODULE Module_ScanFilter

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: clblm_SCANFN, &
             clblm_INTERP, &
             clblm_FLTRFN, &
             SHRKSC, &
             CONVSC

   INTERFACE clblm_SCANFN
      module procedure clblm_SCANFN_array
      module procedure clblm_SCANFN_struct
   END INTERFACE
             
CONTAINS !=======================Module Contains========================


!-----------------------------------------------------------------------
! DRIVER FOR CONVOLVING INSTRUMENTAL SCANNING FUNCTION
! WITH SPECTRUM
! * The input functID is in CLBLM definition, which starts from 1, so it is
!   1 larger than JFN used in program.
!-----------------------------------------------------------------------
   SUBROUTINE clblm_SCANFN_array( spData, spV1, spV2, spDV, spNP, &
                                  dataOut, V1out,V2out,DVout,NPout, functID,functHWHM )                               
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      USE Module_Config     ,ONLY: IPR
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum, &
                                   CLBLM_Spectrum_init
      IMPLICIT NONE
      
      real                   ,intent(in)    :: spData(:)
      real(r8)               ,intent(in)    :: spV1,spV2
      real                   ,intent(in)    :: spDV
      integer                ,intent(in)    :: spNP
      real      ,allocatable ,intent(out)   :: dataOut(:)
      real(r8)               ,intent(inout) :: V1out, V2out
      real                   ,intent(inout) :: DVout
      integer                ,intent(inout) :: NPout
      integer                ,intent(in)    :: functID
      real                   ,intent(in)    :: functHWHM
      
      type(CLBLM_Spectrum)  :: spectObj
      
      
      call CLBLM_Spectrum_init( spectObj, spV1,spDV,spNP )
      spectObj%spect(1:spNP) = spData(1:spNP)
      
      call clblm_SCANFN_struct( spectObj, V1out,V2out,DVout,functID,functHWHM )      
                         
      call move_alloc( spectObj%spect, dataOut )
      V1out = spectObj%V1
      V2out = spectObj%V2
      DVout = spectObj%DV
      NPout = spectObj%NLIM
      
   END SUBROUTINE
   
!-----------------------------------------------------------------------
   SUBROUTINE clblm_SCANFN_struct( spect, V1out,V2out,DVout, functID,functHWHM )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8
      USE Module_Spectrum     ,ONLY: CLBLM_Spectrum,&
                                     CLBLM_Spectrum_init
      USE Module_Config       ,ONLY: IPR
      IMPLICIT NONE !REAL*8          (V)


      type(CLBLM_Spectrum)  ,intent(inout) :: spect !spect may be modified by both INTERP and SCANFN.
      real(r8)              ,intent(inout) :: V1out
      real(r8)              ,intent(inout) :: V2out
      real                  ,intent(inout) :: DVout
      integer               ,intent(in)    :: functID
      real                  ,intent(in)    :: functHWHM
      !integer               ,intent(in)    :: argJVAR
      !real                  ,intent(in)    :: argSAMPL


      real    ,allocatable :: oBuffer(:)
      integer              :: processedInPoints
      integer              :: processedOutPoints
real(r8) :: processedInV2 !for mimicing LBLRTM step scheme.

      integer  ,PARAMETER :: I_1=1, I_1000=1000

      integer  :: I,M
      integer  :: INIT
      logical  :: OP
      real     :: PI
      real     :: S(3850),R1(5000),N1(5000)
      integer  :: NFNMAX
      integer  :: NLIMF,NSHIFT,NREN
      integer  :: JFN
      real     :: HWHM
      integer  :: JVAR
      real     :: SAMPL
      real(r8) :: V1,V2
      integer  :: IFN
      real     :: HWF,DXF
      integer  :: NF,NFMAX
      real     :: SAMPLE, XSCALE
      real     :: DVINT
      real     :: XF(6018)
      real     :: DV
      real(r8) :: V1I,V2I
      real     :: DVI
      integer  :: NNI
      real     :: DVO
      integer  :: IRATIO
      real     :: BOUND
      integer  :: NBOUND
      integer  :: MAXF
      integer  :: ILO,IHI,NLO,NHI
      integer  :: IEOFSC,IPANEL,ISTOP,IDATA
      real(r8) :: VFT,VBOT,VTOP
      integer  :: JFLG



      !-----------------------------------------------------------------
      !
      !    ADDITIONAL SCANNING FUNCTIONS MAY READILY BE ADDED TO THOSE
      !      CURRENTLY IMPLEMENTED IN THIS VERSION OF LBLRTM:
      !
      !    A SHAPE SUBROUTINE FOR THE DESIRED FUNCTION MUST BE CREATED-
      !     THIS SUBROUTINE PRECALCULATES THE FUNCTION FOR SUBSEQUENT
      !      LOOKUP.  SEE FOR EXAMPLE SUBROUTINE SHAPEG FOR THE GAUSSIAN
      !
      !    THE SHAPE SUBROUTINE SETS UP THE SYMMETRIC FUNCTION IN ARRAY FG
      !     AT EQUAL INCREMENTS OF THE HALFWIDTH, 'DXF'. THE VALUE OF 'DXF'
      !     IS SET IN THIS SUBROUTINE BY THE VALUE OF 'DXJ(?)'
      !
      !    A DATA CARD MUST BE CREATED FOR EACH SCANNING FUNCTION DEFINING
      !    THE FOLLOWING QUANTITIES:
      !
      !    HWJ(?)   EXTENT OF THE FUNCTION (BOUND) FROM THE CENTER IN UNITS
      !               OF HALFWIDTH
      !
      !    DXJ(?)   INCREMENT AT WHICH THE FUNCTION IS STORED IN UNITS
      !               OF HALFWIDTH
      !
      !    NJ(?)    THE NUMBER OF POINTS FROM THE CENTER TO THE FUNCTION
      !               BOUND
      !
      !    NJMAX(?) SIZE OF THE ARRAY IN WHICH THE FUNCTION IS STORED
      !               FUNCTION VALUES BETWEEN NJ AND NJMAX ARE ZERO
      !
      !    SMPL(?)  DEFAULT VALUE OF THE SAMPLING INCREMENT IN RECIPRICAL
      !               HALFWIDTH UNITS: E.G. A VALUE OF FOUR MEANS THAT THE
      !               OUTPUT SPACING, 'DV', IN WAVENUMBERS WILL BE 1/4 THE
      !               HALFWIDTH VALUE IN WAVENUMBERS, 'HWHM'.
      !
      !    XSCAL(?) REQUIRED FOR PERIODIC FUNTIONS. THE VALUE OF THE
      !               FUNCTION ARGUMENT IN RADIANS FOR WHICH THE
      !               FUNCTION VALUE IS 0.5, E.G.
      !                   SINX/X = 0.5 FOR X = 1.89549425, XSCAL(4)
      !
      !    CONSIDERATION MUST BE GIVEN TO THE ISSUE OF FUNCTION
      !      NORMALIZATION FOR FUNCTIONS THAT DO NOT HAVE RAPID
      !      CONVERGENCE TO ZERO (SINX/X)
      !
      !                                                               SAC
      !
      !-----------------------------------------------------------------
      character*8  :: HSCNID(0:6), SCANID
      real         :: HWJ(0:6)
      real         :: DXJ(0:6)
      integer      :: NJ(0:6)
      integer      :: NJMX(0:6)
      real         :: SMPLJ(0:6)
      real         :: XSCAL(0:6)

      DATA HSCNID(0) / 'RECTANGL'/,HWJ(0) / 1.         /,               &
     &     DXJ(0) / 0.0  /,NJ(0) / 0    /,NJMX(0) / 0    /,             &
     &     SMPLJ(0) / .5 /,XSCAL(0) / 0.          /
      DATA HSCNID(1) / 'TRIANGLE'/,HWJ(1) / 2.         /,               &
     &     DXJ(1) / 0.02 /,NJ(1) / 101  /,NJMX(1) / 251  /,             &
     &     SMPLJ(1) / 2. /,XSCAL(1) / 0.          /
      DATA HSCNID(2) / 'GAUSS   '/,HWJ(2) / 4.         /,               &
     &     DXJ(2) / 0.02 /,NJ(2) / 201  /,NJMX(2) / 251  /,             &
     &     SMPLJ(2) / 4. /,XSCAL(2) / 0.          /

      ! SINCSQ: 54.18 HALFWIDTHS CORRESPONDS TO 24 ZERO CROSSINGS
      ! PI CORRESPONDS TO X=2.257609141
      DATA HSCNID(3) / 'SINCSQ  '/,HWJ(3) / 54.1826    /,               &
     &     DXJ(3) / 0.02 /,NJ(3) / 2710 /,NJMX(3) / 2760 /,             &
     &     SMPLJ(3) / 4. /,XSCAL(3) / 1.391557377 /

     ! SINC: 119.33 HALFWIDTHS CORRESPONDS TO 72 ZERO CROSSINGS
     ! PI CORRESPONDS TO X=1.657400255
      DATA HSCNID(4) / 'SINC    '/,HWJ(4) / 119.332818 /,               &
     &     DXJ(4) / 0.02 /,NJ(4) / 5968 /,NJMX(4) / 6018 /,             &
     &     SMPLJ(4) / 4. /,XSCAL(4) / 1.89549425  /
      DATA HSCNID(5) / 'VRCTCENT'/,HWJ(5) / 1.         /,               &
     &     DXJ(5) / 0.0  /,NJ(5) / 0    /,NJMX(5) / 0    /,             &
     &     SMPLJ(5) / .5 /,XSCAL(5) / 0.          /
      DATA HSCNID(6) / 'VRCTLEFT'/,HWJ(6) / 1.         /,               &
     &     DXJ(6) / 0.0  /,NJ(6) / 0    /,NJMX(6) / 0    /,             &
     &     SMPLJ(6) / .5 /,XSCAL(6) / 0.          /




      !--- ASSIGN CVS VERSION NUMBER TO MODULE
      !yma HVRPST = '$Revision: 16421 $'

      PI = 2.*ASIN(1.)

      !--- SET THE MAXIMIM NUMBER OF AVAILABLE FUNCTIONS:
      !yma NFNMAX = 6 !function 5,6 (FOV correction) are no longer available in CLBLM.
      NFNMAX = 4

      !--- NLIMF IS ONE MORE THAN THE SIZE OF OUTPUT (CONVOLVED) ARRAY
      NLIMF = 2401
      NREN = 0
      NSHIFT = 32

      !yma SUMOUT = 0.
      !yma SMIN = 999999.
      !yma SMAX = -99999.
      !yma DVOSAV = 0.
      !yma SUMR(1) = SUMOUT
      !yma SUMR(2) = SMIN
      !yma SUMR(3) = SMAX
      !yma SUMR(4) = DVOSAV

      JFN    = abs(functID)-1; JFN=sign(JFN, functID) !CLBLM functID start from 1, so minus 1 here.
      HWHM   = functHWHM
      V1     = V1out
      V2     = V2out
      SAMPL  = - abs(DVout) !SAMPL<0, user is supplying DV
      JVAR   = 0             !argJVAR
      !NPTS   = argNPTS

      IF (HWHM.LE.0.) GO TO 70

      ! JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)
      ! THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2
      ! HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED

      IFN = ABS(JFN)
      IF (IFN.GT.NFNMAX) THEN
         WRITE (IPR,*) 'SCANF; JFN GT LIMIT'
         STOP 'SCANF; JFN GT LIMIT'
      ENDIF

      READ (HSCNID(IFN),905) SCANID

      HWF = HWJ(IFN)
      DXF = DXJ(IFN)
      NF = NJ(IFN)
      NFMAX = NJMX(IFN)
      SAMPLE = SMPLJ(IFN) !if SAMPL==0
      XSCALE = XSCAL(IFN)


      ! CHECK FOR NEGATIVE JFN OR NEGATIVE SAMPL
      !
      ! FOR NEGATIVE JFN, USER IS SUPPLYING FIRST ZERO CROSSING FOR THE
      ! PERIODIC FUNCTION IN HWHM.  SET HWHM=(FIRST ZERO)/(PI/XSCALE)
      !
      ! For JFN=5,6 user is supplying instrument field of view half angle
      ! in degrees in HWHM.  Trap if JFN=-5,-6.
      !
      IF (JFN.LT.0) THEN
         JFN = ABS(JFN)
         IF ((JFN.EQ.3).OR.(JFN.EQ.4)) THEN
            HWHM = HWHM/(PI/XSCALE)
         ELSE
            WRITE (IPR,910) JFN
            STOP 'SCANFN; INVALID JFN'
         ENDIF
      ENDIF


      ! SET DVINT TO DETERMINE IF INTERPOLATION IS NECESSARY
      ! - For JFN = 5,6, set DVINT to 1/12 the width of the first box.
      !   HWHM should carry the value of the field of view half angle
      !   (in degrees).  This is converted to radians.  The box width
      !   formula is
      !
      !              width = V1*(1/2 angle FOV)**2/2
      !
      !   and the degrees-to-radians formula is
      !
      !              rad = deg*3.141592654/180.
      !
      ! - For JFN not equal to 5 or 6, set DVINT to 1/12 the value of
      !   HWHM.  HWHM should carry the true value of the Half Width
      !   at Half Maximum of the scanning function at this point.
      !
      IF ((JFN.EQ.5).OR.(JFN.EQ.6)) THEN
         DVINT = V1*(HWHM*3.141592654/180.)**2/24
      ELSE
         DVINT = HWHM/12.
      ENDIF

      ! - For positive SAMPL, set SAMPLE equal to SAMPL (the number
      !   of points per half width).
      ! - For negative SAMPL, user is supplying desired DELVO
      !   (outgoing spectral spacing).  SAMPLE (the number of sample
      !   points per half width) is set such that SAMPLE=HWHM/DELVO
      !   (Half Width at Half Max over user input outgoing spectral
      !   spacing), and the outgoing spectral spacing DVO will be
      !   recalculated using HWHM and SAMPLE below.
      !
      IF (SAMPL.LT.0.) THEN
         SAMPLE = HWHM/(-SAMPL)
      ELSEIF (SAMPL.GT.0.) THEN
         SAMPLE = SAMPL
      ENDIF

      ! SET UP SELECTED SCANNING FUNCTION:
      IF (JFN.EQ.1) CALL SHAPET( XF, NF,NFMAX,DXF )
      IF (JFN.EQ.2) CALL SHAPEG( XF, NF,NFMAX,DXF )
      IF (JFN.EQ.3) CALL SINCSQ( XF,XSCALE, NF,NFMAX,DXF )
      IF (JFN.EQ.4) CALL SINC(   XF,XSCALE, NF,NFMAX,DXF )


      ! IF DV NOT FINE ENOUGH, FIRST INTERPOLATE
      ! The spect object is used for input and output, i.e, the
      ! interpolated results will be stored in spect too.
      DV = spect%DV
      IF (DV.GT.DVINT) THEN
         call clblm_INTERP( spect, spect%V1, spect%V2, DVINT )
      ENDIF


      DV = spect%DV
      DVI = DV
      !yma DVSAV = DVI

      ! Compute output spectral spacing.  For JFN not 5 or 6, at this
      ! point HWHM always contains the value of the Half Width at Half
      ! Maximum of the scanning function, and SAMPLE always contains
      ! the number of points per half width of the scanning function.
      !
      ! For JFN = 5,6 at this point, HWHM contains the value of the
      ! field of view half angle (in degrees), and SAMPLE contains
      ! the ratio of the field of view half angle to the specified
      ! output spectral spacing (the quotient of HWHM and SAMPLE
      ! results in the circuitous calculation of the previously input
      ! DVO).
      !
      DVO = HWHM/SAMPLE
      IF (JFN.EQ.0) THEN
         IRATIO = DVO/DVI+0.5
         DVO = REAL(IRATIO)*DVI
         IF (IRATIO.LT.2) THEN
            WRITE (IPR,950)
            RETURN !yma151017 GO TO 10
         ENDIF
      ENDIF

      ! BOUND AT THIS POINT IS THE WAVENUMBER VALUE
      ! OF HALF THE SCANNING FUNCTION
      BOUND = HWF*HWHM
      DV = DVO
      !yma V1C = V1
      !yma V2C = V2
      !yma SCIND = JVAR+10*(JFN+10*(JEMIT))
      !yma XSCID = SCIND+0.01
      !yma XHWHM = HWHM
      !yma CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)
      WRITE (IPR,955) HWHM,BOUND,-99,V1,V2,DVO  !WRITE (IPR,955) HWHM,BOUND,JFILE,V1,V2,DVO


      ! NBOUND IS THE NUMBER OF SPECTRAL VALUES SPANNED
      ! BY THE FULL SCANNING FUNCTION
      !
      ! RESET BOUND BASED ON NBOUND
      !
      NBOUND = (2.*HWF)*SAMPLE+0.01

      BOUND =  REAL(NBOUND)*DVO/2.
      MAXF = NLIMF+2*NBOUND+NSHIFT

      !yma SUMIN = 0.

      IEOFSC = 1
      NLO = NSHIFT+1
      NHI = NLIMF+NSHIFT-1  !NLIMF is 2401 not 2400
      DO I = 1, MAXF
         N1(I) = 0.
         R1(I) = 0.
      ENDDO
      INIT = 0
      IDATA = -1  !flag to indicate data present or not in the in-buffer
      IPANEL = -1 !flag to indicate out-buffer full or not
      JFLG = -1   !flag to control the increment or not of NB
      VFT = V1- REAL(NSHIFT)*DV
      VBOT = V1-BOUND
      VTOP = V2+BOUND

      !--- Allocate the output buffer. The size is an approximate at this time.
      allocate( oBuffer( ceiling((VTOP-VBOT)/DVO + 1.) ))
      processedInPoints = 0
      processedOutPoints = 0
processedInV2 = spect%V1-spect%DV !for mimicing LBLRTM step scheme


   40 continue !loop until the output panels is full

      IF (IEOFSC.LE.0) GO TO 60 !end of panel loop


      ! READ DATA TO BE CONVOLVE FROM the spectral object AND PUT INTO ARRAY S
      !yma CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)
      call panelInput( S,&
                       NREN,&
                       V1I,V2I,DVI,NNI, & !out
                       VBOT, VTOP,&
                       ILO,IHI, & !out
                       IEOFSC,IDATA,& !out
                       spect, processedInPoints,&
                       processedInV2 )


      IF (IEOFSC.LE.0) GO TO 60 !end of panel loop


      ! SHRKSC MAY SHRINK (COMPRESS) THE DATA; DVI IS MODIFIED ACCORDINGL
      IF ((JFN.NE.0).AND.(JFN.NE.5).AND.(JFN.NE.6)) THEN
         CALL SHRKSC( INIT,HWHM,&
                      S,N1,V1I,V2I,DVI,NNI,VBOT,VTOP,ILO,IHI,NREN )

      ENDIF

   50 CONTINUE

      ! PERFORM THE CONVOLUTION OF XF ON S TO GIVE R1
      IF (JFN.EQ.0) THEN

         CALL CNVRCT( S,HWHM,R1, &!XF)
                      V1I,V2I,DVI,NNI,NLO,NHI,V1,V2,DVO, &
                      IDATA,IPANEL,JFLG )

      ELSEIF (JFN.EQ.5) THEN

         CALL CNVVRC( S,HWHM,R1, &!XF)
                      V1I,V2I,DVI,NNI,NLO,NHI,V1,V2,DVO,&
                      IDATA,IPANEL,JFLG )

      ELSEIF (JFN.EQ.6) THEN

         CALL CNVVRL( S,HWHM,R1, &!XF)
                      V1I,V2I,DVI,NNI,NLO,NHI,V1,V2,DVO,&
                      IDATA,IPANEL,JFLG )
      ELSE

         CALL CONVSC( S,HWHM,R1,XF, &
                      HWF,DXF,V1I,DVI,ILO,IHI,MAXF,VFT,V1,DVO,&
                      JVAR,IDATA,IPANEL )

      ENDIF

      IF (IPANEL.EQ.0) GO TO 40

   60 CONTINUE


      ! OUTPUT PANEL TO a output buffer, NPTS VALUES OF R1
      !yma IF (JFN.EQ.0.OR.JFN.EQ.5.OR.JFN.EQ.6) THEN
      !yma    CALL PNLRCT (R1,JFILE,SUMR,NPTS)
      !yma ELSE
      !yma    CALL PANLSC (R1,JFILE,SUMR,NPTS)
      !yma ENDIF
      call panelOutput( oBuffer, processedOutPoints,&
                        R1,DVO,NLO,NHI,NSHIFT,NLIMF,MAXF,&
                        VFT,V2,ISTOP )
      IF (JFN.EQ.0.OR.JFN.EQ.5.OR.JFN.EQ.6) IPANEL=-1


      IF ((ISTOP.NE.1).AND.(IEOFSC.LT.0)) GO TO 60
      IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 50


      !--- Fill the output spectrum object. The input spect is
      !    deallocated and reused for the output spect.
      call CLBLM_Spectrum_init( spect, V1,DVO,processedOutPoints )
      spect%spect( spect%indV1:&
                   spect%indV2 ) = oBuffer(1:processedOutPoints)
      !yma CALL ENDFIL (JFILE)

      V1out = spect%V1
      V2out = spect%V2
      DVout = spect%DV
      
      deallocate(oBuffer)


      !yma SUMIN = SUMIN*DVSAV
      !yma WRITE (IPR,975) SUMIN

      !yma SUMOUT = SUMR(1)
      !yma SMIN = SUMR(2)
      !yma SMAX = SUMR(3)
      !yma DVOSAV = SUMR(4)
      !yma SUMOUT = SUMOUT*DVOSAV
      !yma WRITE (IPR,980) SUMOUT,SMIN,SMAX


   70 CONTINUE

   80 RETURN

  900 FORMAT (3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)
  905 FORMAT (A8)
  910 FORMAT (//,' *****  INVALID VALUE FOR JFN = ',I2,'  *****',/)
  915 FORMAT ('1',' **SCANFN** ',/,'0',10A8,2X,2(1X,A8,1X))
  920 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
  925 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',  &
     &        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/     &
     &        '0 V2(CM-1) = ',F12.6)
  930 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',     &
     &        1PE10.3,/(5X,A6,' = ',1PE10.3))
  935 FORMAT ('0',' **SCANFN** ',/)
  940 FORMAT (A4,I2.2)
  945 FORMAT ('0','***',A8,'***',//6X,'INPUT FILE NUMBER =',I3,         &
     &        ' ,IFILST = ',I5,' ,NIFILS = ',I5,',JEMIT =',I2,          &
     &        ' ,JFN =',I2,' ,JVAR =',I2,'  ,JABS =',I2)
  950 FORMAT ('0',60X,'****** IRATIO LESS THAN 2, NO SCANFN ******')
  955 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,   &
     &        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,    &
     &        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',    &
     &        F12.5,5X,' DV OUT',F12.8)
  960 FORMAT (///,'0',5X,A12,/)
  965 FORMAT ('0',5X,'IEOFSC =',I3,'  IDATA =',I3,'  IPANEL =',I3,/)
  970 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3,  &
     &        ',  PANEL =',F6.3)
  975 FORMAT ('0    SUMIN  =',1P,E16.9)
  980 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)

      END SUBROUTINE

!-----------------------------------------------------------------------
!
!     SUBROUTINE SHAPET SETS UP THE TRIANGULAR SCANNING FUNCTION
!
!-----------------------------------------------------------------------
      SUBROUTINE SHAPET( XF, NF,NFMAX,DXF )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: XF(:)
      integer ,intent(in)  :: NF,NFMAX
      real    ,intent(in)  :: DXF

      integer :: I
      real    :: SUM,X
      real    :: XTRIAN
      XTRIAN(X) = 1.-0.5*X



      DO 10 I = 1, NFMAX
         XF(I) = 0.
   10 END DO
      XF(1) = 0.5
      SUM = XF(1)
      DO 20 I = 2, NF
         X = REAL(I-1)*DXF
         XF(I) = 0.5*XTRIAN(X)
         SUM = SUM+2.*XF(I)
   20 END DO
      SUM = SUM*DXF

      !PRT  WRITE(IPR,900) NF,DXF,SUM

      RETURN

  900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE SHAPEG CONSTRUCTS THE FUNCTION FOR THE DOPPLER PROFILE
!
! 160608: This is a copy of SHAPEG from clblm_odlay.f90. This duplication is solely
!         for resolving mutually module dependencies:
!         clblm_Post use clblm_odlay/SHAPEG  and clblm_odlay use clblm_Post/CONVSC,SHRKSC
!-----------------------------------------------------------------------
      SUBROUTINE SHAPEG( FG,NX1,N1MAX,DXF1 )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: FG(:)
      integer ,intent(in)  :: NX1,N1MAX
      real    ,intent(in)  :: DXF1

      integer :: I,JJ
      real    :: FGNORM,FLN2,RECPI,SUM
      real    :: TOTAL,X,XSQ
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
!
!     SUBROUTINE SINCSQ SETS UP THE SINCSQ SCANNING FUNCTION
!
!-----------------------------------------------------------------------
      SUBROUTINE SINCSQ( XF,XSCALE, NF,NFMAX,DXF )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: XF(:)
      real    ,intent(in)  :: XSCALE
      integer ,intent(in)  :: NF, NFMAX
      real    ,intent(in)  :: DXF

      integer :: I
      real    :: PI,SUM,X,XNORM
      !DATA XSCALE / 1.391557377 /
      real :: XSINC2
      XSINC2(X) = (SIN(X)/X)**2


      PI = 2.*ASIN(1.)

      !PI CORRESPONDS TO X=2.257609141

      XNORM = XSCALE/PI
      DO 10 I = 1, NFMAX
         XF(I) = 0.
   10 END DO
      XF(1) = XNORM
      SUM = XF(1)
      DO 20 I = 2, NF
         X = REAL(I-1)*DXF
         XF(I) = XNORM*XSINC2(X*XSCALE)
         SUM = SUM+2.*XF(I)
   20 END DO
      SUM = SUM*DXF

      !PRT  WRITE(IPR,900) NF,DXF,SUM

      RETURN

  900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)

      END SUBROUTINE


!-----------------------------------------------------------------------
!
!     SUBROUTINE SINC SETS UP THE SINC SCANNING FUNCTION
!
!-----------------------------------------------------------------------
      SUBROUTINE SINC( XF,XSCALE, NF,NFMAX,DXF )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real    ,intent(out) :: XF(:)
      real    ,intent(in)  :: XSCALE
      integer ,intent(in)  :: NF, NFMAX
      real    ,intent(in)  :: DXF

      integer :: I
      real    :: PI,SUM,X,XNORM
      !DATA XSCALE / 1.89549425  /
      real :: XSINC
      XSINC(X) = (SIN(X)/X)


      PI = 2.*ASIN(1.)

      ! PI CORRESPONDS TO X=1.657400255

      XNORM = XSCALE/PI
      DO 10 I = 1, NFMAX
         XF(I) = 0.
   10 END DO
      XF(1) = XNORM
      SUM = XF(1)
      DO 20 I = 2, NF
         X = REAL(I-1)*DXF
         XF(I) = XNORM*XSINC(X*XSCALE)
         SUM = SUM+2.*XF(I)
   20 END DO
      SUM = SUM*DXF

      !PRT  WRITE(IPR,900) NF,DXF,SUM

      RETURN

  900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)

      END SUBROUTINE


!-----------------------------------------------------------------------
!
!     THIS SUBROUTINE COMPRESSES (SHRINKS) THE INPUT TO THE CONVOLUTION
!     ROUTINE FOR THE SCANNING FUNCTION TO ACCELERATE THE CALCULATION
!
!-----------------------------------------------------------------------
      SUBROUTINE SHRKSC( INIT,HWHM,&
                         S,SS,V1I,V2I,DVI,NLIM,VBOT,VTOP,ILO,IHI,NREN )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE !REAL*8          (V)


      integer   ,intent(inout) :: INIT
      real      ,intent(in)    :: HWHM

      real      ,intent(inout) :: S(:) !S(3850)
      real      ,intent(inout) :: SS(:) !SS(200)
      real(r8)  ,intent(inout) :: V1I
      real(r8)  ,intent(out)   :: V2I
      real      ,intent(inout) :: DVI
      integer   ,intent(inout) :: NLIM
      real(r8)  ,intent(in)    :: VBOT
      real(r8)  ,intent(in)    :: VTOP
      integer   ,intent(out)   :: ILO
      integer   ,intent(out)   :: IHI
      integer   ,intent(inout) :: NREN


      !--- Local variables
      !
      integer ,PARAMETER :: I_1=1
      integer ,PARAMETER :: JRATIO(24)=(/ 1,2,3,4,5,6,8,10,12,15,16,20,&
                                          24,25,30,32,40,48,50,  &
                                          60,75,80,100,120 /)

      real      ,SAVE :: DVSC
      real      ,SAVE :: SRATIO
      real(r8)  ,SAVE :: V1SHFT
      integer   ,SAVE :: IRATM1
      integer   ,SAVE :: IRATSH


      INTEGER  :: I,    II,   IMAX, IMIN
      INTEGER  :: J,    JHI
      INTEGER  :: K,    NLIMS
      REAL     :: SUMK




!yma      CALL CPUTIM (TIME0)
      NLIMS = NLIM
      IF (NREN.GT.0) THEN
         DO 10 I = 1, NREN
            S(I) = SS(I)
   10    CONTINUE
      ENDIF
      NREN = 0
      IF (INIT.EQ.0) THEN
         DVSC = HWHM/12.
         IRATSH = DVSC/DVI+0.5
         DO 20 I = 2, 24
            IF (JRATIO(I).GT.IRATSH) THEN
               IRATSH = JRATIO(I-1)
               GO TO 30
            ENDIF
   20    CONTINUE
   30    IF (IRATSH.GT.JRATIO(24)) IRATSH = JRATIO(24)
         IF (IRATSH.LE.1) RETURN
         DVSC = REAL(IRATSH)*DVI
         V1SHFT = REAL(IRATSH-1)*DVI/2.
         WRITE (IPR,900) IRATSH
         SRATIO = IRATSH
         IRATM1 = IRATSH-1
         INIT = 1
      ENDIF
      IF (IRATSH.LE.1) RETURN
      NREN = NLIM-(NLIM/IRATSH)*IRATSH

      !PRT  WRITE(IPR,905) V1I,V1SHFT,DVSC,NREN

      V1I = V1I+V1SHFT
      IMIN = 1
      IMAX = NLIM-IRATM1-NREN

      K = 0
      DO 50 I = IMIN, IMAX, IRATSH
         SUMK = 0.
         JHI = I+IRATM1
         K = K+1
         DO 40 J = I, JHI
            SUMK = SUMK+S(J)
   40    CONTINUE
         S(K) = SUMK/SRATIO
   50 END DO

      V2I = V1I+DVSC* REAL(K-1)
      NLIM = K
      DVI = DVSC
      ILO = ((VBOT-V1I)/DVI)+1.5
      ILO = MAX(ILO,I_1)
      IHI = ((VTOP-V1I)/DVI)+1.5
      IHI = MIN(IHI,NLIM)

      !PRT  WRITE(IPR,910) ILO,IHI

      IF (NREN.GT.0) THEN
         DO 60 I = 1, NREN
            II = NLIMS-NREN+I
            SS(I) = S(II)
   60    CONTINUE
      ENDIF
!yma      CALL CPUTIM (TIME)
!yma      TIMCNV = TIMCNV+TIME-TIME0

      RETURN

  900 FORMAT ('   SHRINK RATIO = ',I5)
  905 FORMAT ('   V1I =',F10.3,'  V1SHFT =',F10.3,'  DVSC =',F12.5,     &
     &        '   NREN =',I4)
  910 FORMAT ('   ILO =',I4,'  IHI =',I4)

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SUBROUTINE CONVSC PERFORMS THE CONVOLUTION WITH THE SELECTED
!     SCANNING FUNCTION
!-----------------------------------------------------------------------
      SUBROUTINE CONVSC( S,HWHMV1,R1,XF ,&
                         HWF,DXF,V1I,DVI,ILO,IHI,MAXF,VFT,V1,DVO,&
                         JVAR,IDATA,IPANEL )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE !REAL*8          (V)

      real     ,intent(in)    :: S(:)
      real     ,intent(in)    :: HWHMV1
      real     ,intent(inout) :: R1(:)
      real     ,intent(in)    :: XF(:)
      !
      real     ,intent(in)    :: HWF
      real     ,intent(in)    :: DXF
      real(r8) ,intent(in)    :: V1I
      real     ,intent(in)    :: DVI
      integer  ,intent(inout) :: ILO
      integer  ,intent(in)    :: IHI
      integer  ,intent(in)    :: MAXF
      real(r8) ,intent(in)    :: VFT
      real(r8) ,intent(in)    :: V1
      real     ,intent(in)    :: DVO
      integer  ,intent(in)    :: JVAR
      integer  ,intent(in)    :: IDATA
      integer  ,intent(out)   :: IPANEL

      !---Local variables
      !
      real :: RATIO
      real :: SUMIN=0.

      INTEGER  :: I,     ILAST,  IT,   ITST
      INTEGER  :: JF,    JMAX,   JMIN
      REAL     :: DVODX, HWBND,  HWHM, SI
      REAL(r8) :: VI
      REAL     :: XNORM, ZBOUND, ZF,   ZINT
      REAL     :: ZPEAK, ZSLOPE

      integer ,PARAMETER :: I_1=1


      !CALL CPUTIM (TIME0)
      IF (ILO.GT.IHI) GO TO 60
      RATIO = DVI/DVO
      DVODX = DVO/DXF
      HWBND = HWF/DVO
      ZINT = ((V1I-VFT)/DVO)
      HWHM = HWHMV1
      ITST = -1

      DO 50 I = ILO, IHI
         IF (S(I).EQ.0.) GO TO 50
         IF (I.LT.ITST) GO TO 20
         ITST = 9999
         IF (JVAR.EQ.0) GO TO 10
         VI = REAL(I-1)*DVI+V1I
         HWHM = HWHMV1*(VI/V1)**2
         ITST = I+ INT(1./DVI)
   10    CONTINUE
         ZSLOPE = DVODX/HWHM
         ZBOUND = HWBND*HWHM
         XNORM = DVI/HWHM

!PRT         WRITE(IPR,900) VI,HWHM

   20    CONTINUE
         ZPEAK = REAL(I-1)*RATIO+ZINT
         JMAX = ZPEAK+ZBOUND+1.5
         IF (JMAX.LE.MAXF) GO TO 30
         ILAST = I-1
         GO TO 60

   30    JMIN = ZPEAK-ZBOUND+1.5
         JMIN = MAX(JMIN,I_1)
         SUMIN = SUMIN+S(I)
         SI = XNORM*S(I)
         ZF = ( REAL(JMIN-1)-ZPEAK)*ZSLOPE
         DO 40 JF = JMIN, JMAX
            IT = ABS(ZF)+1.5
            R1(JF) = R1(JF)+SI*XF(IT)
            ZF = ZF+ZSLOPE
   40    CONTINUE

   50 END DO
      ILAST = IHI
      IPANEL = IDATA
      GO TO 70

   60 IPANEL = 1
   70 continue
!yma      CALL CPUTIM (TIME)
!yma      TIMCNV = TIMCNV+TIME-TIME0
      ILO = ILAST+1

      RETURN

  900 FORMAT ('0 AVE PANEL WAVENUMBER = ',F12.4,5X,'HWHM = ',F10.5)

      END SUBROUTINE

!-----------------------------------------------------------------------
!
!     SUBROUTINE CNVVRC PERFORMS THE CONVOLUTION WITH A RECTANGULAR
!     SCANNING FUNCTION OF VARIABLE SIZE, WHERE THE BOX SIZE IS
!     WAVENUMBER DEPENDENT.  V1, V2, AND DVO ARE USED TO DEFINE THE
!     CENTER OF THE OUTPUT BOXES.  BOXES OVERLAP WHERE NECESSARY TO
!     INSURE A CONSTANT DVO.
!
!     BFOV is used to determine the resolution (box size), which is
!     spectrally variable.
!
!     AFOV is passed in from calls to CNVVRC from SCANFN and SCNMRG
!     as HWHM, since the value of HWHM on Record 8.1 on TAPE5 holds
!     the place of the half angle of the instrument field of view in
!     degrees.
!
!     BOX WIDTH EQUALS V*B**2/2, AND THE SHIFT EQUALS HALF THE BOX WIDTH
!     V*(1-B**2/4), WHERE B IS THE HALF ANGLE OF THE INSTRUMENT FIELD
!     OF VIEW IN RADIANS.
!
!     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT
!     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHI
!     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
!     BY THE SUM OF THE WEIGHTS.
!
!-----------------------------------------------------------------------
      SUBROUTINE CNVVRC( S,AFOV,R1, &!XF)
                         V1I,V2I,DVI,NNI,NLO,NHI,V1,V2,DVO, &
                         IDATA,IPANEL,JFLG )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8

      IMPLICIT NONE !REAL*8          (V)

      real    ,intent(in)    :: S(:)
      real    ,intent(in)    :: AFOV
      real    ,intent(inout) :: R1(:)
      !real    ,intent() :: XF(:)
      !
      real(r8) ,intent(in)    :: V1I, V2I
      real     ,intent(in)    :: DVI
      integer  ,intent(inout) :: NNI
      integer  ,intent(in)    :: NLO, NHI
      real(r8) ,intent(in)    :: V1, V2
      real     ,intent(in)    :: DVO
      !
      integer  ,intent(in)    :: IDATA
      integer  ,intent(inout) :: IPANEL
      integer  ,intent(inout) :: JFLG


      !---Local variabls
      !
      real     ,SAVE :: SUMJ
      real     ,SAVE :: RNJ
      integer  ,SAVE :: JN
      integer  ,SAVE :: NB
      integer  ,SAVE :: IPC

      real(r8) :: VLFT,VCNT,VRGT
      real     :: WGTL,WGTR
      real     :: RATIO
!yma      real     :: SUMIN

      INTEGER  :: I,     IH,   IL, NEP
      INTEGER  :: NNIV2
      REAL     :: BFOV,  RL,   RR
      REAL(r8) :: VLBLL, VLBLR


      ! LBLRTM flags
      ! JFLG = -1:  first time through; increment NB when box is full
      ! JFLG =  0:  subsequent calls:increment NB when box full
      ! JFLG =  1:  out of data, return for more; do not increment NB
      ! IDATA =-1:  first time through; or, need more data
      ! IDATA = 0:  data present
      ! IDATA = 1:  no data present
      ! IPANEL=-1:  first time through; or, after panel written
      ! IPANEL= 0:  panel not full
      ! IPANEL= 1:  panel is full

      !CALL CPUTIM (TIME0)

      ! Convert AFOV to BFOV, the half angle of the instrument field of
      ! view in radians. (For IRIS-D, AFOV equals 2.5 degrees)
      BFOV = AFOV*3.141592654/180.

      RATIO = DVO/DVI

      ! During first call or if entering after writing a panel,
      ! initialize: SUMJ (radiance sum), and
      !             RNJ (accumulator for number of input points in
      !                  current box, i.e. the sum of the weights)
      !             JN (box counter from 1 at VFT)
      ! During first call only,
      ! initialize: NB (box counter from 1 at V1),
      !             IPC (output panel counter).

      IF (IPANEL.EQ.-1) THEN
         SUMJ = 0.
         RNJ = 0.
         JN = NLO
         IF (JFLG.EQ.-1) THEN
            NB = 1
            IPC = 1
         ENDIF
      ENDIF

      ! Check that number of points in current panel, NNI, is correct.

      NNIV2 = (V2I-V1I)/DVI+1.0001
      IF (NNI.GT.NNIV2) NNI = NNIV2

      ! Top of loop over NB boxes

   10 IF (NLO.LE.NHI) THEN

         ! For current box find wavenumber at center and left/right edges.
         ! For first box, VCNT equals V1.  When current box exceeds V2,
         ! then exit.

         VCNT = V1+(NB-1)*DVO
         IF (VCNT.GT.V2) THEN
            IPANEL = 1
            RETURN
         ENDIF

         VLFT = VCNT*(1-BFOV**2/4)
         VRGT = VCNT*(1+BFOV**2/4)

         ! Find lbl panel indices for points which fall within current box.

         RL = (VLFT-V1I)/DVI+1
         RR = (VRGT-V1I)/DVI+1

         IL = INT(RL+0.5)
         IH = INT(RR+0.5)

         ! Calculate weight for each end point, inner points weighted as 1.
         ! NEP is the number of endpoints in use.

         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = (VLBLR-VLFT)/DVI
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = (VRGT-VLBLL)/DVI
         NEP = 2

         ! Set flag if last data point on current input panel reached

         IF (IH.GT.NNI) THEN
            IH = NNI
            JFLG = 1

            ! If retrieving next panel while box sum is in progress, then
            ! use weight of 1. for temp. right endpoint at IH = NNI = 2400
            ! calculate partial sum below, then return.  If only one point
            ! is included in this sum (IL = IH = NNI), use weight of 0
            ! for right point, and add only left endpoint to sum.
            VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
            WGTR = 1.
            IF (IL.EQ.IH) THEN
               WGTR = 0.
               NEP = 1
            ENDIF
         ENDIF

         ! If returning with new panel to partially summed box, then set
         ! weight for temporary left endpoint to 1.  If it's the last
         ! point going into the box, then count it as final right endpoint
         ! and use weight of 0 for left point (since IL = IH = 1).

         IF (IL.LE.1) THEN
            IL = 1
            VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
            WGTL = 1.
            IF (IL.GE.IH) THEN
               WGTL = 0.
               NEP = 1
            ENDIF
         ENDIF

         ! If retrieving next panel while box sum is not progress, then
         ! check that left edge of current output box is beyond last data
         ! panel (IL.GT.NNI), if so, go back for new panel without summing

         IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
            IPANEL = 0
            JFLG = 0
            RETURN
         ENDIF

         ! If last point on current input panel is reached, and there is
         ! no more data to retrieve, then return

         IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN

         ! Compute sum for current box number NB, for all points but
         ! the end points

         DO 20 I = IL+1, IH-1
            SUMJ = SUMJ+S(I)
   20    CONTINUE

         ! Add weighted end points to sum
         SUMJ = SUMJ+S(IL)*WGTL
         SUMJ = SUMJ+S(IH)*WGTR

         !  Define sum of the weights, where all inner points are weighted
         !  as 1, and the end points are weighted with the fraction that
         !  occurs within the box.

         RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR

         ! If out of data on current input panel, go back for more;
         ! partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV

         IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN
            IPANEL = 0
            JFLG = 0
            RETURN
         ENDIF

         ! IPANEL=IDATA

!yma         SUMIN = SUMIN+SUMJ


         ! Compute average radiance for completed box

         R1(JN) = SUMJ/RNJ

         ! Increment current box counters, initialize SUMJ and RNJ
         JN = JN+1
         SUMJ = 0.
         RNJ = 0.

         ! Output panel when number of boxes, NB, reaches a multiple of
         ! 2400, using then incrementing current output panel number, IPC.

         IF (NB.EQ.IPC*(NHI-NLO+1)) THEN
            IPANEL = 1
            IPC = IPC+1
            NB = NB+1
            RETURN
         ENDIF

         ! Increment NB
         NB = NB+1

         ! Go back to top of loop over NB boxes

         GO TO 10

      ENDIF

      RETURN
      END SUBROUTINE


!-----------------------------------------------------------------------
!
!     SUBROUTINE CNVVRL PERFORMS THE CONVOLUTION WITH A RECTANGULAR
!     SCANNING FUNCTION OF VARIABLE SIZE, WHERE THE BOX SIZE IS
!     WAVENUMBER DEPENDANT.  V1, V2, AND DVO ARE USED TO DEFINE THE
!     LEFT EDGE OF THE OUTPUT BOXES.  BOXES OVERLAP WHERE NECESSARY
!     TO INSURE A CONSTANT DVO.
!
!     BFOV is used to determine the resolution (box size), which is
!     spectrally variable.
!
!     AFOV is passed in from calls to CNVVRL from SCANFN and SCNMRG
!     as HWHM, since the value of HWHM on Record 8.1 on TAPE5 holds
!     the place of the half angle of the instrument field of view in
!     degrees.
!
!     BOX WIDTH EQUALS V*B**2/2, AND THE SHIFT EQUALS HALF THE BOX WIDTH
!     V*(1-B**2/4), WHERE B IS THE HALF ANGLE OF THE INSTRUMENT FIELD
!     OF VIEW IN RADIANS.
!
!     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT
!     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHI
!     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
!     BY THE SUM OF THE WEIGHTS.
!
!-----------------------------------------------------------------------
      SUBROUTINE CNVVRL( S,AFOV,R1, &!XF)
                         V1I,V2I,DVI,NNI,NLO,NHI,V1,V2,DVO,&
                         IDATA,IPANEL,JFLG )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8

      IMPLICIT NONE !REAL*8          (V)

      real     ,intent(in)    :: S(:)
      real     ,intent(in)    :: AFOV
      real     ,intent(inout) :: R1(:)
      !real     ,intent() :: XF(:)
      !
      real(r8) ,intent(in)    ::  V1I, V2I
      real     ,intent(in)    ::  DVI
      integer  ,intent(inout) ::  NNI
      integer  ,intent(in)    ::  NLO,NHI
      real(r8) ,intent(in)    ::  V1,V2
      real     ,intent(in)    ::  DVO
      !
      integer  ,intent(in)    ::  IDATA
      integer  ,intent(inout) ::  IPANEL
      integer  ,intent(inout) ::  JFLG


      !---Local variables
      !
      real     ,SAVE ::  SUMJ
      real     ,SAVE ::  RNJ
      integer  ,SAVE ::  JN
      integer  ,SAVE ::  NB
      integer  ,SAVE ::  IPC

      real(r8) :: VLFT, VCNT, VRGT
      real     :: WGTL, WGTR
      real     :: RATIO
!yma      real     :: SUMIN

      INTEGER  :: I,     IH,   IL, NEP
      INTEGER  :: NNIV2
      REAL     :: BFOV,  RL,   RR
      REAL(r8) :: VLBLL, VLBLR


      ! LBLRTM flags
      ! JFLG = -1:  first time through; increment NB when box is full
      ! JFLG =  0:  subsequent calls:increment NB when box full
      ! JFLG =  1:  out of data, return for more; do not increment NB
      ! IDATA =-1:  first time through; or, need more data
      ! IDATA = 0:  data present
      ! IDATA = 1:  no data present
      ! IPANEL=-1:  first time through; or, after panel written
      ! IPANEL= 0:  panel not full
      ! IPANEL= 1:  panel is full

      ! CALL CPUTIM (TIME0)

      ! Convert AFOV to BFOV, the half angle of the instrument field of
      ! view in radians. (For IRIS-D, AFOV equals 2.5 degrees)
      BFOV = AFOV*3.141592654/180.

      RATIO = DVO/DVI

      ! During first call or if entering after writing a panel,
      ! initialize: SUMJ (radiance sum), and
      !             RNJ (accumulator for number of input points in
      !                  current box, i.e. the sum of the weights)
      !             JN (box counter from 1 at VFT)
      ! During first call only,
      ! initialize: NB (box counter from 1 at V1),
      !             IPC (output panel counter).
      !
      IF (IPANEL.EQ.-1) THEN
         SUMJ = 0.
         RNJ = 0.
         JN = NLO
         IF (JFLG.EQ.-1) THEN
            NB = 1
            IPC = 1
         ENDIF
      ENDIF

      ! Check that number of points in current panel, NNI, is correct.
      NNIV2 = (V2I-V1I)/DVI+1.0001
      IF (NNI.GT.NNIV2) NNI = NNIV2

      ! Top of loop over NB boxes
   10 IF (NLO.LE.NHI) THEN

         ! For current box find wavenumber at the left and right edges.
         ! For first box, VLFT equals V1.  When current box exceeds V2,
         ! then exit.
         !
         VLFT = V1+(NB-1)*DVO
         IF (VLFT.GT.V2) THEN
            IPANEL = 1
            RETURN
         ENDIF

         VCNT = VLFT*(1+BFOV**2/4)
         VRGT = VLFT*(1+BFOV**2/2)

         ! Find lbl panel indices for points which fall within current box.
         RL = (VLFT-V1I)/DVI+1
         RR = (VRGT-V1I)/DVI+1

         IL = INT(RL+0.5)
         IH = INT(RR+0.5)

         ! Calculate weight for each end point, inner points weighted as 1,
         ! NEP is the number of endpoints in use.
         !
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = (VLBLR-VLFT)/DVI
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = (VRGT-VLBLL)/DVI
         NEP = 2

         ! Set flag if last data point on current input panel reached
         IF (IH.GT.NNI) THEN
            IH = NNI
            JFLG = 1

            ! If retrieving next panel while box sum is in progress, then
            ! use weight of 1. for temp. right endpoint at IH = NNI,
            ! calculate partial sum below, then return.  If only one point
            ! is included in this sum (IL = IH = NNI), use weight of 0
            ! for right point, and add only left endpoint to sum.
            VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
            WGTR = 1.
            IF (IL.EQ.IH) THEN
               WGTR = 0.
               NEP = 1
            ENDIF
         ENDIF

         ! If returning with new panel to partially summed box, then set
         ! weight for temporary left endpoint to 1.  If it's the last
         ! point going into the box, then count it as final right endpoint
         ! and use weight of 0 for left point (since IL = IH = 1).
         !
         IF (IL.LE.1) THEN
            IL = 1
            VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
            WGTL = 1.
            IF (IL.GE.IH) THEN
               WGTL = 0.
               NEP = 1
            ENDIF
         ENDIF

         ! If retrieving next panel while box sum is not in progress, then
         ! check that left edge of current output box is beyond last data
         ! panel (IL.GT.NNI), if so, go back for new panel without summing
         !
         IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
            IPANEL = 0
            JFLG = 0
            RETURN
         ENDIF

         ! If last point on current input panel is reached, and there is
         ! no more data to retrieve, then return
         !
         IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN

         ! Compute sum for current box number NB, using a weight of 1.0
         ! for all points but the end points, which use WGTL and WGTR
         !
         DO 20 I = IL+1, IH-1
            SUMJ = SUMJ+S(I)
   20    CONTINUE

         ! Add weighted end points to sum
         SUMJ = SUMJ+S(IL)*WGTL
         SUMJ = SUMJ+S(IH)*WGTR

         ! Define sum of the weights, where all inner points are weighted
         ! as 1, and the end points are weighted with the fraction that
         ! occurs within the box.
         !
         RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR

         ! If out of data on current input panel, go back for more;
         ! partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV
         !
         IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN
            IPANEL = 0
            JFLG = 0
            RETURN
         ENDIF

         ! IPANEL=IDATA

!yma         SUMIN = SUMIN+SUMJ

         ! Compute average radiance for completed box
         R1(JN) = SUMJ/RNJ

         ! Increment current box counters, initialize SUMJ and RNJ
         JN = JN+1
         SUMJ = 0.
         RNJ = 0.

         ! Output panel when number of boxes, NB, reaches a multiple of
         ! 2400, using then incrementing current output panel number, IPC.
         !
         IF (NB.EQ.IPC*(NHI-NLO+1)) THEN
            IPANEL = 1
            IPC = IPC+1
            NB = NB+1
            RETURN
         ENDIF

         ! Increment NB
         NB = NB+1

      ! Go back to top of loop over NB boxes
      GO TO 10
      ENDIF

      RETURN
      END  SUBROUTINE

!-----------------------------------------------------------------------
!
!     SUBROUTINE CNVRCT PERFORMS THE CONVOLUTION WITH AN ALTERNATE
!     RECTANGULAR SCANNING FUNCTION (ADJACENT BOXES OF ONE SIZE,
!     EQUAL TO 2*HWHM)
!
!     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT
!     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHI
!     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
!     BY THE SUM OF THE WEIGHTS.
!
!-----------------------------------------------------------------------
      SUBROUTINE CNVRCT( S,HWHM,R1, &!XF
                         V1I,V2I,DVI,NNI,NLO,NHI,V1,V2,DVO, &
                         IDATA,IPANEL,JFLG )
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8

      IMPLICIT NONE !REAL*8          (V)


      real     ,intent(in)    :: S(:)
      real     ,intent(in)    :: HWHM
      real     ,intent(inout) :: R1(:)
      !real     ,intent() :: XF(:)
      !
      real(r8) ,intent(in)    :: V1I,V2I
      real     ,intent(in)    :: DVI
      integer  ,intent(inout) :: NNI
      integer  ,intent(in)    :: NLO,NHI
      real(r8) ,intent(in)    :: V1,V2
      real     ,intent(in)    :: DVO
      !
      integer  ,intent(in)    :: IDATA
      integer  ,intent(inout) :: IPANEL
      integer  ,intent(inout) :: JFLG


      !---Local variables
      !
      real     ,SAVE :: SUMJ
      real     ,SAVE :: RNJ
      integer  ,SAVE :: JN
      integer  ,SAVE :: NB
      integer  ,SAVE :: IPC

      real(r8) :: VLFT,VCNT,VRGT
      real     :: WGTL,WGTR
      real     :: RATIO
      !real     :: SUMIN

      INTEGER  :: I,     IH,   IL,   ILPR
      INTEGER  :: NEP,   NNIV2
      REAL     :: RL,    RR
      REAL(r8) :: VLBLL, VLBLR



      ! LBLRTM flags
      ! JFLG = -1:  first time through; increment NB when box is full
      ! JFLG =  0:  subsequent calls:increment NB when box full
      ! JFLG =  1:  out of data, return for more; do not increment NB
      ! IDATA =-1:  first time through; or, need more data
      ! IDATA = 0:  data present
      ! IDATA = 1:  no data present
      ! IPANEL=-1:  first time through; or, after panel written
      ! IPANEL= 0:  panel not full
      ! IPANEL= 1:  panel is full

      !CALL CPUTIM (TIME0)
      RATIO = DVO/DVI

      ! During first call or if entering after writing a panel,
      ! initialize: SUMJ (radiance sum), and
      !             RNJ (accumulator for number of input points in
      !                  current box, i.e. the sum of the weights)
      !             JN (box counter from 1 at VFT)
      ! During first call only,
      ! initialize: NB (box counter from 1 at V1),
      !             IPC (output panel counter).
      !
      IF (IPANEL.EQ.-1) THEN
         SUMJ = 0.
         RNJ = 0.
         JN = NLO
         IF (JFLG.EQ.-1) THEN
            NB = 1
            IPC = 1
         ENDIF
      ENDIF

      ! Check that number of points in current panel, NNI, is correct.

      NNIV2 = (V2I-V1I)/DVI+1.0001
      IF (NNI.GT.NNIV2) NNI = NNIV2

      ! Top of loop over NB boxes

   10 IF (NLO.LE.NHI) THEN

         VCNT = V1+(NB-1)*DVO
         IF (VCNT.GT.V2) THEN
            IPANEL = 1
            RETURN
         ENDIF

         VLFT = VCNT-HWHM
         VRGT = VCNT+HWHM

         ! Find lbl panel indices for points which fall within current box.

         RL = (VLFT-V1I)/DVI+1
         RR = (VRGT-V1I)/DVI+1

         IL = INT(RL+0.5)
         IH = INT(RR+0.5)

         ! Calculate weight for each end point, inner points weighted
         ! as 1.  NEP is the number of endpoints in use.

         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = (VLBLR-VLFT)/DVI
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = (VRGT-VLBLL)/DVI
         NEP = 2

         ! Set flag if last data point on current input panel reached

         IF (IH.GT.NNI) THEN
            IH = NNI
            JFLG = 1

            ! If retrieving next panel while box sum is in progress, then
            ! use weight of 1. for temp. right endpoint at IH = NNI = 2400
            ! calculate partial sum below, then return.  If only one point
            ! is included in this sum (IL = IH = NNI), use weight of 0
            ! for right point, and add only left endpoint to sum.
            VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
            WGTR = 1.
            IF (IL.EQ.IH) THEN
               WGTR = 0.
               NEP = 1
            ENDIF
         ENDIF

         ! If returning with new panel to partially summed box, then set
         ! weight for temporary left endpoint to 1.  If it's the last
         ! point going into the box, then count it as final right endpoint
         ! and use weight of 0 for left point (since IL = IH = 1).
         !
         IF (IL.LE.1) THEN
            IL = 1
            VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
            WGTL = 1.
            IF (IL.GE.IH) THEN
               WGTL = 0.
               NEP = 1
            ENDIF
         ENDIF

         ! If retrieving next panel while box sum is not progress, then
         ! check that left edge of current output box is beyond last data
         ! panel (IL.GT.NNI), if so, go back for new panel without summing
         !
         IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
            IPANEL = 0
            JFLG = 0
            RETURN
         ENDIF

         ! If last point on current input panel is reached, and there is
         ! no more data to retrieve, then return
         !
         IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN

         ! Compute sum for current box number NB, for all points but
         ! the end points
         !
         DO 20 I = IL+1, IH-1
            SUMJ = SUMJ+S(I)
   20    CONTINUE

         ! Add weighted end points to sum
         SUMJ = SUMJ+S(IL)*WGTL
         SUMJ = SUMJ+S(IH)*WGTR

         ! Define sum of the weights, where all inner points are weighted
         ! as 1, and the end points are weighted with the fraction that
         ! occurs within the box.
         !
         RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR

         ! If out of data on current input panel, go back for more;
         ! partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV
         !
         IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN
            IPANEL = 0
            JFLG = 0
            RETURN
         ENDIF

         ! IPANEL=IDATA
!yma         SUMIN = SUMIN+SUMJ

         ! Compute average radiance for completed box
         R1(JN) = SUMJ/RNJ

         ILPR = IH+1

         ! Increment current box counters, initialize SUMJ and RNJ
         JN = JN+1
         SUMJ = 0.
         RNJ = 0.

         ! Output panel when number of boxes, NB, reaches a multiple of
         ! 2400, using then incrementing current output panel number, IPC.
         !
         IF (NB.EQ.IPC*(NHI-NLO+1)) THEN
            IPANEL = 1
            IPC = IPC+1
            NB = NB+1
            RETURN
         ENDIF

         ! Increment NB
         NB = NB+1

         ! Go back to top of loop over NB boxes
      GO TO 10
      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!      SUBROUTINE FLTRFN (IFILE)
      SUBROUTINE clblm_FLTRFN ( spect, argV1,argDV,argNPTS,argXF, &
                                filtOut )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8, &
                                     NFLTPT !maximum number of points in the incoming filter
      USE Module_Spectrum     ,ONLY: CLBLM_Spectrum
      USE Module_Config       ,ONLY: IPR
      IMPLICIT NONE

      type(CLBLM_Spectrum)   ,intent(in)  :: spect
      real(r8)               ,intent(in)  :: argV1          !!wavenumber of initial filter value
      real                   ,intent(in)  :: argDV          !!wavenumber increment between filter values
      integer                ,intent(in)  :: argNPTS        !!number of filter values
      real                   ,intent(in)  :: argXF(argNPTS) !!NPTS values of filter function
      real                   ,intent(out) :: filtOut



      integer ,PARAMETER :: i_2=2

      real ,allocatable :: XF(:)

      real      :: totRFILTR
      integer   :: processedInPoints
      real(r8)  :: processedInV2

      integer   :: I, IPRT, NREN
      real      :: S(4650)
      integer   :: NPTF, NPTS
      real      :: DVF
      real(r8)  :: V1F,V2F
      real      :: DV
      real(r8)  :: V1,V2
      real(r8)  :: VBOT,VTOP,VFT
      integer   :: NPTF_HALF
      real(r8)  :: V1F_CENTER
      integer   :: ILO,IHI
      integer   :: NNI
      real      :: DVI
      real(r8)  :: V1I,V2I
      real      :: DVC
      real      :: SUMFLT
      real      :: RFILTR
      integer   :: IDATA, IEOFSC



      NREN = 0
      IPRT = 1

      !yma   10 READ (IRD,900) V1F,DVF,NPTF,JEMIT,IUNIT,IFILST,NIFILS,junit,HEDDR
      V1F = argV1
      DVF = argDV

      NPTF = argNPTS
      allocate(XF(NPTF))
      XF(1:NPTF) = argXF(1:argNPTS)

      DVC = spect%DV


      !     Test to ensure NPTF is less than NFLTPT, the maximum number
      !     of filter points allowed
      IF (NPTF.GT.NFLTPT) THEN
         WRITE(IPR,*) 'FLTRFN: NPTS > NFLTPT limit', NFLTPT
         STOP 'FLTRFN: NPTS > NFLTPT limit'
      ENDIF

      IF (V1F.LT.0) RETURN

      WRITE (IPR,905)

      !     DVF < 0 option flags V1F value to be the center frequency
      !     Check that there NPTF is odd (to ensure a center frequency),
      !     save center frequency value, and reset V1F to endpoint value.
      if (DVF.lt.0.) then

         dvf = abs(dvf)

         if (mod((nptf-1),i_2).ne.0) then
            write(*,*) 'Use of V1F as center frequency requires odd     &
     &           number of points'
            write(ipr,*) 'Use of V1F as center frequency requires odd   &
     &           number of points, stopping in FLTFRN'
            stop 'FLTRFN'
         endif

         V1F_center = V1F
         nptf_half = (abs(nptf)-1)/2
         V1F = V1F_center - DVF* REAL(nptf_half)
         write(ipr,*) ' ``````````````````````````````'
         write(ipr,*) ' Use of V1F as center frequency:'
         write(ipr,*) '   V center = ',V1F_center
         write(ipr,*) '   V1F      = ',V1F
         write(ipr,*) " ''''''''''''''''''''''''''''''"
      endif

      IEOFSC = 0

      NPTS = NPTF
   30 V2F = V1F+DVF* REAL(NPTS-1)
      WRITE (IPR,910) V1F,V2F,DVF,NPTF,-9,-9,-9,-9,-9,'-9'  !WRITE (IPR,910) V1F,V2F,DVF,NPTF,JEMIT,JABS,IUNIT,IFILST,NIFILS,HEDDR
      V1 = V1F
      V2 = V2F
      DV = DVF
      !yma WRITE (IPR,CVAR) (XF(I),I=1,NPTS)

      !     MAKE ADJUSTMENT FOR END POINT CORRECTIONS
      XF(1) = 0.5*XF(1)
      XF(NPTS) = 0.5*XF(NPTS)

   40 SUMFLT = 0.0
      DO I = 1, NPTS
         SUMFLT = SUMFLT+XF(I)
      END DO
      SUMFLT = SUMFLT*DVF

   60 RFILTR = 0.0
      VFT = V1
      VBOT = V1
      VTOP = V2
      !TIMRDF = 0.0
      !TIMCNV = 0.0

      !     JEMIT=-1 FILTER PASSED OVER ABSORPTION
      !     JEMIT=0  FILTER PASSED OVER TRANSMISSION
      !     JEMIT=1  FILTER PASSED OVER EMISSION

      processedInPoints = 0
processedInV2 = spect%V1-spect%DV

      IDATA = -1
   80 continue
         !yma CALL CPUTIM (TIMEO)
         !yma CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)
         call panelInput( S,&
                          NREN,&
                          V1I,V2I,DVI,NNI, & !out
                          VBOT, VTOP,&
                          ILO,IHI, & !out
                          IEOFSC,IDATA,& !out
                          spect, processedInPoints,&
                          processedInV2 )
         !yma CALL CPUTIM (TIME)
         !yma TIMRDF = TIMRDF+TIME-TIMEO
         IF (IEOFSC.NE.1) GO TO 90
         CALL CNVFLT( S,RFILTR,XF, DVF,V1F,V2F,DVI,V1I,ILO,IHI,VFT )
         IF (IDATA.EQ.1) GO TO 90
      GO TO 80

   90 continue
      totRFILTR = RFILTR*DVC
      filtOut = totRFILTR/SUMFLT

      deallocate(XF)


  900 FORMAT (2F10.4,6I5,8A4,A3)
  905 FORMAT ('1',/'   ***  FILTER ***',8(' ********** '))
  910 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',   &
     &        I5,/,'0',10X,', JEMIT= ',I2,', JABS= ',I2,                &
     &        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,   &
     &        8A4,A3)
  915 FORMAT (A80)
  920 FORMAT ('0  RESULT FROM SCANNING FUNCTION INCONSISTENT WITH ',    &
     &        'FILTER REQUEST')
  925 FORMAT ('0',//,8(' -------  '),/,'0',10A8,2X,2(1X,A8,1X))
  930 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
  935 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',  &
     &        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/     &
     &        '0 V2(CM-1) = ',F12.6)
  940 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',     &
     &        1PE10.3,/(5X,A6,' = ',1PE10.3))
  945 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',   &
     &        I5,/'0 , IEMIT= ',I2,', JEMIT= ',I2,', JABS= ',I2,        &
     &        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,   &
     &        8A4,A3)
  950 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3)
  955 FORMAT ('0  INTEGRATED TRANSMISSION = ',1PE14.5,                  &
     &        '  NORMALIZATION OF  THE FILTER = ',E14.5,/               &
     &        '0 UNNORMALIZED INTEGRATED TRANSMISSION =  ',E14.5)
  960 FORMAT ('0  INTEGRATED ABSORPTION = ',1PE14.5,                    &
     &        '  NORMALIZATION OF  THE FILTER = ',E14.5,/               &
     &        '0 UNNORMALIZED INTEGRATED ABSORPTION =    ',E14.5)
!  965 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,
!     *        '  NORMALIZATION OF THE',' FILTER = ',E14.5)
  965 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,                       &
     &        '  NORMALIZATION OF THE',' FILTER = ',1PE14.5,            &
     &        ' NORM. EMISSION',E14.5)

  970 format (a4,i2.2)
  975 format (1p,e14.5,1p,e14.5,1p,e14.5,0p)
  980 format (' Filter output:')

      END SUBROUTINE

!------------------------------------------------------------------------
!------------------------------------------------------------------------
      SUBROUTINE CNVFLT( S,RFILTR,XF, DVF,V1F,V2F,DVI,V1I,ILO,IHI,VFT )
!------------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real     ,intent(in)    :: S(:)
      real     ,intent(inout) :: RFILTR
      real     ,intent(in)    :: XF(:)
      !
      real     ,intent(in)    :: DVF
      real(r8) ,intent(in)    :: V1F, V2F
      real     ,intent(in)    :: DVI
      real(r8) ,intent(in)    :: V1I
      integer  ,intent(in)    :: ILO, IHI
      real(r8) ,intent(inout) :: VFT


      INTEGER  :: I,     IFL,   IMAX, IMIN
      REAL     :: P
      REAL(r8) :: V1S,   V2S,   VXF1, VXF2
      REAL     :: XDVIF, XIF0


      !yma CALL CPUTIM (TIMEO)
      IMIN = (V1F-V1I)/DVI+1.0001
      !IMIN = (V1F-V1I)/DVI+1.5
      IMIN = MAX(IMIN,ILO)
      IMAX = (V2F+V1F-V1I)/DVI+1.0001 !YMa: why not IMAX = (V2F-V1I)/DVI+1.0001 ??
      !IMAX = (V2F+V1F-V1I)/DVI+1.5
      IMAX = MIN(IMAX,IHI)
      XIF0 = (V1I-V1F)/DVF+1.0001
      !XIF0 = (V1I-V1F)/DVF+1.5
      XDVIF = DVI/DVF
      
      !--- YMa: Out of XF bound may happen. Reduce the IMAX to avoid this. 
      if ( int(XIF0+XDVIF* REAL(IMAX))+1 >size(XF) ) then
         do i=imax-1,imin,-1
            if ( int(XIF0+XDVIF* REAL(i))+1 <=size(XF) ) then
               IMAX=i
               exit
            endif
         enddo
      endif
      if ( int(XIF0+XDVIF* REAL(IMAX))+1 >size(XF) ) then
         STOP '--- CNVFLT(): Subscript is greater than the upper bound of XF.'
      endif

      
      DO 10 I = IMIN, IMAX
         IFL = XIF0+XDVIF* REAL(I)
         !IFL = XIF0+XDVIF* REAL(I-1)
         
         ! Linearly interpolate filter function XF to avoid
         ! discontinuities in output spectrum

         v1s = v1i+(i-1)*dvi
         v2s = v1i+i*dvi
         vxf1 = v1f + (ifl-1)*dvf
         vxf2 = v1f + (ifl)*dvf
         p = (vxf2-v2s)/(vxf2-vxf1)
         RFILTR = RFILTR+S(I)*(p*XF(IFL)+(1.-p)*xf(IFL+1))
   10 END DO
      IF (IMAX.LT.IHI) VFT = VFT+(( REAL(IHI)- REAL(ILO))+1.0)*DVI
      !yma CALL CPUTIM (TIME)
      !yma TIMCNV = TIMCNV+TIME-TIMEO

      RETURN

      END  SUBROUTINE



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE clblm_INTERP(spect, V1, V2, DVO, argI4PT)
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum,&
                                    CLBLM_Spectrum_init
      USE Module_Config      ,ONLY: IPR
      IMPLICIT NONE

      type(CLBLM_Spectrum) ,intent(inout)  :: spect
      real(r8)             ,intent(in)     :: V1,V2
      real                 ,intent(in)     :: DVO !output DV,
      integer    ,OPTIONAL ,intent(in)     :: argI4PT !need explicit interface.


      !---Local variables
      !
      integer  ,PARAMETER :: I_2400 = 2400
      integer  ,PARAMETER :: NUMCOF = 201

      real     :: C1(0:202),C2(0:202),C3(0:202),C4(0:202)!,RSTAT(3)
      real(r8) :: P,PP
      integer  :: I_0,IP
      real     :: XNUMCF

      integer  :: I4PT
      integer  :: IBOUND !IBOUND = 4
      real     :: BOUND
      real(r8) :: VBOT, VTOP

      real, allocatable :: outBuffer(:)

      integer  :: processedInPoints
      integer  :: processedOutPoints
real(r8) :: processedInV2

      real     :: T(2410), R(2401)
      real(r8) :: bufV1, bufV2
      real     :: bufDV
      integer  :: bufNLIM
      integer  :: ILO, IHI
      integer  :: IEOFSC,IDATA

      real(r8) :: V1I,V2I, V1J,V2J, VI,VJ, VRATIO
      real     :: DVI, DVJ
      integer  :: NNI, NNJ
      integer  :: I,J,J1,J2,II




      !--- For default using 4 points interpolation scheme.
      I4PT = 1
      if (present(argI4PT)) I4PT = argI4PT


      ! SET UP FOUR POINT INTERPOLATION COEFICIENTS FOR P FOR 201
      ! POINTS BETWEEN 0 AND 1.0, with an extra point at each end
      I_0 = 0
      IF (I4PT.NE.0) THEN

         XNUMCF = REAL(NUMCOF)
         DO IP = I_0, NUMCOF+1
            P = ( REAL(IP)-1.0)/(XNUMCF-1.0)
            PP = P**2
            C1(IP) = -P/2.0*(1-P)**2
            C2(IP) = 1.0-PP*(3.0-2.0*P)+PP/2.0*(1.0-P)
            C3(IP) = PP*(3.0-2.0*P)+P/2.0*(1.0-P)**2
            C4(IP) = -PP/2*(1.0-P)
         ENDDO

      ENDIF



      ! V1 IS REQUESTED LOWER V, VBOT = V1-VBOUND.  VBOUND ALLOWS FOR
      ! 4 POINT INTERPOLATION OF THE FIRST DATA POINT.
      ! THE FIRST PANEL OF DATA IS STORED IN T, STARTING AT T(5)
      ! WITH A CORRESPONDING WAVENUMBER V1I.
      ! T(1-4) ARE USED TO STORE THE LAST IBOUND POINTS FROM THE
      ! PREVIOUS PANEL, BUT ARE ZEROED OUT FOR THE FIRST PANEL.
      ! THE INTERPOLATED POINTS ARE STORED IN THE ARRAY R.
      ! THE INDEX I REFERS TO THE INPUT POINTS, J TO THE OUTPUT POINTS.
      ! P VARIES FROM 0 TO 1 AND IS THE FRACTIONAL DISTANCE OF THE
      ! CURRENT OUTPUT WAVENUMBER VJ TO THE NEXT LOWEST INPUT WAVENUMBER
      ! RELATIVE TO THE INPUT DV: P = (VJ-VI(II))/DVI
      ! INRANG IS 0 IF VJ IS WITHIN THE RANGE OF THE INPUT DATA, -1
      ! IF VJ IS LESS THAN THE INPUT DATA, AND 1 IF IT IS GREATER
      !
      ! INITIALIZE THE VARIABLES

      ! NEED TO SAVE LAST IBOUND POINTS OF EACH PANEL TO ATTACH TO NEXT
      IBOUND = 4

      ! VBOT IS LOWEST NEEDED WAVENUMBER, VTOP IS HIGHEST
      BOUND =  REAL(IBOUND)*DVO
      VBOT = V1-BOUND
      VTOP = V2+BOUND

      !--- Allocate the output buffer. The size is an approximate at this time.
      allocate( outBuffer( ceiling( (VTOP-VBOT)/DVO +1.) ))
      processedInPoints = 0
      processedOutPoints = 0
processedInV2 = spect%V1-spect%DV

      ! ZERO OUT T(1 TO IBOUND)
      DO I = 1, IBOUND
         T(I) = 0.0
      END DO

      !--- load the first 2400 points into T, start from T(5)
      call panelInput( T,&
                       IBOUND,&
                       bufV1,bufV2,bufDV,bufNLIM, & !out
                       VBOT, VTOP,&
                       ILO,IHI, & !out
                       IEOFSC,IDATA,& !out
                       spect, processedInPoints,&
                       processedInV2)
      DVI = bufDV
      V1I = bufV1
      V2I = bufV2
      NNI = bufNLIM

      DVJ = DVO
      V1J = V1
      VJ  = V1J
      VRATIO = DVJ/DVI

      !RMIN = 1.0E15
      !RMAX = -1.0
      !RSUM = 0.0

      ! EXTRAPOLATE DOWN TO V1I-DVI SO THAT THE POINT I=4 IS AVAILABLE
      ! FOR THE FIRST PANEL.  THIS ALLOWS 4 POINT INTERPOLATION BETWEEN
      ! V1I AND V1I+DVI
      !T(4) = 2.0*T(5)-T(6)
      T(IBOUND) = 2.0*T(IBOUND+1)-T(IBOUND+2)

      ! LOOP OVER THE OUTPUT PANELS
      ! IF V1J .LT. V1I, THEN ZERO FILL UP TO V1I.
   20 IF (V1J.LT.V1I) THEN
         J1 = 1
         J2 = MIN(INT((V1I-V1J)/DVJ)+1,I_2400)

         ! FILL IN
         DO J = J1, J2
            R(J) = 0.0
         ENDDO

         V2J = V1J+DVJ*(J2-1)
         NNJ = J2

         !CALL OTPANL (R,JFILE,NPTS)
         outBuffer(processedOutPoints+1:&
                   processedOutPoints+NNJ) = R(1:NNJ)
         processedOutPoints = processedOutPoints + NNJ

         V1J = V2J+DVJ
         VJ = V1J

         GO TO 20
      ENDIF

      ! AT THIS POINT, VJ >= V1I
   40 CONTINUE

      ! I INDEXES THE LARGEST VI .LE. VJ
      ! AND AT THIS POINT SHOULD .GE. 1.
      I = (VJ-V1I)/DVI+1.00001
      IF (I.LT.1) THEN
         WRITE (IPR,*) ' INTERP-ERROR: I SHOULD >= 1, IS ',I
         print*, ' INTERP-ERROR: I SHOULD >= 1, IS ',I
         STOP
      ENDIF
      VI = V1I+DVI* REAL(I-1)

      ! P IS INCREMENTED BY ADDING DVJ/DVI BUT WILL BE REINITIALIZED
      ! HERE FOR EACH OUTPUT PANEL TO AVOID THE ACCUMULATION OF
      ! TRUNCATION ERRORS
      P = (VJ-VI)/DVI

      J1 = INT((VJ-V1J)/DVJ+1.001)
      J2 = MIN(INT((V2-V1J)/DVJ+1.001),INT((V2I-DVI-V1J)/DVJ+1.),I_2400)

      ! LOOP OVER A SINGLE OUTPUT PANEL
      IF (I4PT.GT.0) THEN

         ! 4 POINT INTERPOLATION
         DO J = J1, J2

            ! PERFORM INTERPOLATION
            IP = P*XNUMCF+1.00001
            R(J) = C1(IP)*T(I-1)+C2(IP)*T(I)+C3(IP)*T(I+1)+ C4(IP)*T(I+2)

            !! ACCUMULATE STATISTICS
            !RMIN = MIN(RMIN,R(J))
            !RMAX = MAX(RMAX,R(J))
            !RSUM = RSUM+R(J)

            ! INCREMENT P AND I
            P = P+VRATIO
            IF (P.GE.1.0) THEN
               I = I+P
               P = P- REAL(INT(P))
            ENDIF

         ENDDO !DO J = J1, J2

      ELSE

         ! LINEAR INTERPOLATION
         DO J = J1, J2

            ! PERFORM INTERPOLATION
            R(J) = T(I)*(1.0-P)+T(I+1)*P

            !! ACCUMULATE STATISTICS
            !RMIN = MIN(RMIN,R(J))
            !RMAX = MAX(RMAX,R(J))
            !RSUM = RSUM+R(J)

            ! INCREMENT P AND I
            P = P+VRATIO
            IF (P.GE.1.0) THEN
               I = I+P
               P = P- REAL(INT(P))
            ENDIF
         ENDDO !DO J = J1, J2

      ENDIF !!IF (I4PT.GT.0) THEN

      ! VJ IS THE FREQUENCY OF THE NEXT OUTPUT POINT (NOT THE LAST
      ! POINT IS THE CURRENT PANEL)
      VJ = V1J+DVJ*J2

      ! IF THE OUTPUT PANEL IS FULL OR IF V2 REACHED,
      ! WRITE OUT THE PANEL
      IF (J2.GE.I_2400.OR.VJ.GE.V2) THEN
         NNJ = J2
         V2J = V1J+DVJ*(J2-1)

         !CALL OTPANL (R,JFILE,NPTS)
         outBuffer(processedOutPoints+1:&
                   processedOutPoints+NNJ) = R(1:NNJ)
         processedOutPoints = processedOutPoints + NNJ

         V1J = V2J+DVJ
         J2 = 0
      ENDIF

      ! IF REACHED V2, THEN FINISH
      IF (VJ.GE.V2) GO TO 100

      ! IF THE INPUT FILE REACHED AN EOF, THEN ZERO FILL TO END
      IF (IEOFSC.LE.0) GO TO 80

      ! IF THE DATA FROM CURRENT INPUT PANEL IS EXHAUSTED, GET MORE
      IF (I.GE.NNI-2) THEN

         ! SHIFT THE LAST IBOUND POINTS DOWN TO T(1-4)
         DO II = 1, IBOUND
            T(II) = T(II+NNI-IBOUND)
         ENDDO

         ! GET THE NEXT PANEL OF DATA AND RESET I
         !CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)
         call panelInput( T,&
                          IBOUND,& !IBOUND=4
                          bufV1,bufV2,bufDV,bufNLIM, & !out
                          VBOT, VTOP,&
                          ILO,IHI, & !out
                          IEOFSC,IDATA,& !out
                          spect, processedInPoints,&
                          processedInV2)

         DVI = bufDV
         V1I = bufV1
         V2I = bufV2
         NNI = bufNLIM

         IF (IEOFSC.LE.0) THEN

            ! IF EOF ON INPUT FILE, THEN EXTRAPOLATE OUT TWO MORE
            ! POINTS BEYOND I=NNI SO THAT 4 POINT INTERPOLATION CAN
            ! BE PERFORMED UP TO VJ=V2I. (ACTUALLY, ONLY T(NNI+1) NEED
            ! BE EXTRAPOLATED, T(NNI+2) NEED ONLY BE DEFINED.)
            ! EXTEND THE INPUT PANEL BY ONE POINT  AND LOOP AROUND THE
            ! INTERPOLATION ONE LAST TIME
            T(NNI+1) = 2.0*T(NNI)-T(NNI-1)
            T(NNI+2) = 0.0
            V2I = V2I+DVI
         ENDIF

      ENDIF

      ! LOOP BACK
      GO TO 40

   80 CONTINUE
      J1 = J2+1
      J2 = MIN(INT((V2-V1J)/DVJ+1.0001),I_2400)

      DO J = J1, J2
         R(J) = 0.0
      ENDDO
      VJ = V1J+DVJ*J2

      ! IF THE OUTPUT PANEL IS FULL OR IF V2 REACHED,
      ! WRITE OUT THE PANEL
      IF (J2.GE.I_2400.OR.VJ.GE.V2) THEN
         NNJ = J2
         V2J = V1J+DVJ*(J2-1)

         !CALL OTPANL (R,JFILE,NPTS)
         outBuffer(processedOutPoints+1:&
                   processedOutPoints+NNJ) = R(1:NNJ)
         processedOutPoints = processedOutPoints+NNJ

         V1J = V2J+DVJ
         J2 = 0
      ENDIF

      ! IF REACHED V2, THEN FINISH
      IF (VJ.LT.V2) GO TO 80

  100 CONTINUE

      !RSTAT(1) = RSUM*DVJ
      !RSTAT(2) = RMIN
      !RSTAT(3) = RMAX

      !--- Fill the output spectrum object. The input spect is
      !    deallocated and reused for the output spect.
      !
      call CLBLM_Spectrum_init( spect, V1,DVO,processedOutPoints )
      spect%spect( spect%indV1:&
                   spect%indV2 ) = outBuffer(1:processedOutPoints)

      deallocate(outBuffer)


      END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE panelInput( iBuffer,NREN,&
                             VMIN,VMAX,DVI,NNI, &
                             VBOT,VTOP,ILO,IHI, &
                             IEOFSC,IDATA,&
                             spect,processedNLIM, &
                             processedV2 )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: EXPMIN, ARGMIN, r8=>kind_r8
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum

      IMPLICIT NONE

      real                  ,intent(inout) :: iBuffer(:)
      integer               ,intent(in)    :: NREN
      real(r8)              ,intent(out)   :: VMIN, VMAX
      real                  ,intent(out)   :: DVI
      integer               ,intent(out)   :: NNI
      real(r8)              ,intent(in)    :: VBOT, VTOP
      integer               ,intent(out)   :: ILO
      integer               ,intent(out)   :: IHI
      integer               ,intent(out)   :: IEOFSC
      integer               ,intent(out)   :: IDATA
      type(CLBLM_Spectrum)  ,intent(in)    :: spect
      integer               ,intent(inout) :: processedNLIM
      real(r8)              ,intent(inout) :: processedV2 !for the purpose of mimicing LBLRTM spectral increment method


      integer ,PARAMETER :: panelSize = 2400

      integer :: i,j,k
      integer :: NLIM
      integer :: NLOW



      !--- Load a panel. Only the points that higher than VBOT get loaded
      IEOFSC = 1
      DO
         NLIM = panelSize
         if ( processedNLIM+NLIM > spect%NLIM ) NLIM=spect%NLIM-processedNLIM

         if (NLIM==0) then !No more data available in the inSpect
            IEOFSC = 0
            RETURN
         endif

         DVI  = spect%DV
!         VMIN = spect%V1 + processedNLIM*DVI
VMIN = processedV2+DVI  !for mimicing LBLRTM scheme.
         VMAX = VMIN + (NLIM-1)*DVI


         processedNLIM = processedNLIM + NLIM
processedV2 = VMAX !for mimicing LBLRTM scheme.

         IF (VMAX.GE.VBOT) EXIT

      ENDDO



      !--- NREN is the number reminder points from the last panel
      NLOW = NREN+1
      IF (NREN.LE.0) NLOW = 1

      VMIN = VMIN-(NLOW-1)*DVI
      NNI = NLIM+NLOW-1

      IDATA = 0  !There are still not-processed points in the inSpect
      iBuffer(NLOW:NNI) = spect%spect( spect%indV1+processedNLIM-NLIM : &
                                       spect%indV1+processedNLIM-1 )

      ILO = 1
      IHI = NNI
      !IF (VMAX.LE.VTOP) RETURN
      !IHI = (VTOP-VMIN)/DVI+1.5
      !IDATA = 1
      !RETURN
      IF (VMAX.GT.VTOP) THEN
         IHI = (VTOP-VMIN)/DVI+1.5
         IDATA = 1 !All needed points have been done.
      ENDIF


      END SUBROUTINE


!-----------------------------------------------------------------------
!     OUTPUTS THE RESULTS OF THE SCANNING FUNCTION TO an output buffer
!-----------------------------------------------------------------------
      !SUBROUTINE PNLRCT(R1,JFILE,SUMR,NPTS)
      SUBROUTINE panelOutput( oBuffer, processedOutPoints,&
                              R1,DVO,NLO,NHI,NSHIFT,NLIMF,MAXF,&
                              VFT,V2,ISTOP )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      IMPLICIT NONE

      real      ,intent(inout) :: oBuffer(:)
      integer   ,intent(inout) :: processedOutPoints
      real      ,intent(inout) :: R1(:)
      real      ,intent(in)    :: DVO
      integer   ,intent(inout) :: NLO,NHI
      integer   ,intent(in)    :: NSHIFT,NLIMF,MAXF
      real(r8)  ,intent(inout) :: VFT
      real(r8)  ,intent(in)    :: V2
      integer   ,intent(out)   :: ISTOP
      !integer   ,intent(out)   :: IPANEL

      integer :: I,J,IJ
      integer :: NNHI,NLIM
      real    :: DV
      !real(r8):: V1P,V2P


      !CALL CPUTIM (TIME0)

      !SUMOUT = SUMR(1)
      !SMIN = SUMR(2)
      !SMAX = SUMR(3)

      DV = DVO
      ISTOP = 0
      NNHI = (V2-VFT)/DV+1.5
      IF (NHI.GE.NNHI) ISTOP = 1
      IF (ISTOP.EQ.1) NHI = NNHI
      NLIM = NHI-NLO+1

      !V1P = VFT+ REAL(NLO-1)*DV
      !V2P = VFT+ REAL(NHI-1)*DV
      ! V1P IS FIRST FREQ OF PANEL
      ! V2P IS LAST  FREQ OF PANEL

      !CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)
      !CALL BUFOUT (JFILE,R1(NLO),NLIM)
      oBuffer(processedOutPoints+1:&
              processedOutPoints+NLIM) = R1(NLO:NHI)
      processedOutPoints = processedOutPoints+NLIM

      VFT = VFT+ REAL(NLIMF-1)*DV
      !IF (NPTS.GT.0) THEN
      !   WRITE (IPR,900) V1P,V2P,DVO,NLIM
      !   WRITE (IPR,905)
      !   NNPTS = NPTS
      !   IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1
      !   IJLIM = NLIM-NNPTS+1
      !   DO 10 IJ = 1, NNPTS
      !      IK = IJ+IJLIM-1
      !      VI = V1P+ REAL(IJ-1)*DVO
      !      VK = V1P+ REAL(IK-1)*DVO
      !      JJ = NLO+IJ-1
      !      KK = NLO+IK-1
      !      WRITE (IPR,910) IJ,VI,R1(JJ),IK,VK,R1(KK)
      !10   CONTINUE
      !ENDIF

      !NLIMHI = NLIM+NLO-1
      !DO I = NLO, NLIMHI
      !   SMIN = MIN(SMIN,R1(I))
      !   SMAX = MAX(SMAX,R1(I))
      !   SUMOUT = SUMOUT+R1(I)
      !ENDDO

      IF (ISTOP.EQ.1) GO TO 50
      DO J = NLIMF, MAXF
         R1(J-NLIMF+1) = R1(J)
      ENDDO
      DO J = MAXF-NLIMF+2, MAXF
         R1(J) = 0.
      ENDDO
      NLO = NSHIFT+1
   50 continue
      !IPANEL = -1
      !SUMR(1) = SUMOUT
      !SUMR(2) = SMIN
      !SUMR(3) = SMAX
      !SUMR(4) = DVO
      !CALL CPUTIM (TIME)
      !TIMPNL = TIMPNL+TIME-TIME0
      RETURN

  900    FORMAT('0 V1P =',F12.5,' V2P =',F12.5,' DVOUT =',F12.8,        &
     &   ' NLIM =',I10)
  905    FORMAT('0')
  910    FORMAT(I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)

      END SUBROUTINE

END MODULE


