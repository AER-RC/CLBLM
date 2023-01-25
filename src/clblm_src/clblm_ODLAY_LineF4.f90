!
! CREATION HISTORY:
!       Modified from LBLRTM v12.9
!       Yingtao Ma, AER@NOAA/NESDIS/STAR
!       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_LineF4

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: LBLF4Q, &
             LINF4Q



CONTAINS !=================== Module Contains ==========================


!-----------------------------------------------------------------------
!
!     SUBROUTINE LBLF4 DOES A LINE BY LINE CALCULATION
!     USING FUNCTION F4.
!
!-----------------------------------------------------------------------
      SUBROUTINE LBLF4Q( JRAD,V1R4,V2R4, &
                         R4,RR4,NPTR4,DVR4,BOUND4,TAVE,DPTFAC,DPTMIN,NLTE,&
                         LNFIL4,NEGEPP_FLAG,ILNFLG,JCNVF4,ILIN4,ILIN4T )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, RADCN2, ONEPL
      USE Module_LineData         ,ONLY: LINE_SHRINK
      USE Module_ODLAY_CommonSub  ,ONLY: RADFN, RADFNI, AVRAT

      IMPLICIT NONE !REAL*8           (V)

      integer   ,intent(in)    :: JRAD
      real(r8)  ,intent(in)    :: V1R4
      real(r8)  ,intent(inout) :: V2R4 !V2R4 may be adjusted
      !
      real      ,intent(out)   :: R4(:)
      real      ,intent(out)   :: RR4(:)
      integer   ,intent(out)   :: NPTR4
      real      ,intent(in)    :: DVR4
      real      ,intent(in)    :: BOUND4
      real      ,intent(in)    :: TAVE
      real      ,intent(in)    :: DPTFAC
      real      ,intent(in)    :: DPTMIN
      logical   ,intent(in)    :: NLTE
      integer   ,intent(in)    :: LNFIL4
      integer   ,intent(in)    :: NEGEPP_FLAG
      integer   ,intent(in)    :: ILNFLG
      integer   ,intent(inout) :: JCNVF4
      integer   ,intent(inout) :: ILIN4  !out for information
      integer   ,intent(inout) :: ILIN4T !out for information


      !--- Local variables
      !
      type(LINE_SHRINK)  :: LINE
      integer            :: JLBLF4=0
      real               :: tmpRADVI

      real(r8)  :: VLO, VHI
      integer   :: ILO, IHI
      integer   :: LIMIN
      integer   :: LIMOUT
      real      :: DPTFC
      real      :: DPTMN

      INTEGER  :: I,      IEOF,    LIMP2
      INTEGER  :: NPTSI1, NPTSI2
      REAL     :: BETA,   RADVI,     RDEL
      REAL     :: RDLAST, TIM0,      TIM00,     TIM1
      REAL     :: TIM2,   TIM3,      TIM4,      TIM5
      REAL(r8) :: VI,     VITST
      REAL     :: XKT




      !CALL CPUTIM (TIM00)

!yma  DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)
      DPTMN = DPTMIN/RADFN(V2R4,TAVE/RADCN2)
      DPTFC = DPTFAC
      LIMIN = 1000
      LIMOUT = 2500
      JLBLF4 = 1

!     SET IEOF EQUAL TO -1 FOR FIRST READ

      IEOF = -1

!yma      V1R4 = V1
!yma      V2R4 = V2
      NPTR4 = (V2R4-V1R4)/DVR4+ONEPL
      NPTR4 = MIN(NPTR4,LIMOUT)
      V2R4 = V1R4+DVR4* REAL(NPTR4-1)

      LIMP2 = LIMOUT+2
      DO I = 1, LIMP2
         R4(I) = 0.
      ENDDO
      if ( NLTE ) then !yma,151201
         DO I = 1, LIMP2
            RR4(I) = 0.
         END DO
      endif

      BETA = RADCN2/TAVE
      VLO = V1R4-BOUND4
      VHI = V2R4+BOUND4
   20 continue
         !CALL CPUTIM (TIM0)
         CALL RDLN4Q( LINE,IEOF, &
                      LNFIL4,NEGEPP_FLAG,VLO,VHI,ILO,IHI,ILIN4T )
         !CALL CPUTIM (TIM1)


         IF (IEOF.EQ.2) THEN
            !TF4 = TF4+TIM1-TIM00
            RETURN
         ENDIF

         !TF4RDF = TF4RDF+TIM1-TIM0
         !TIM2 = TIM1
         IF (IEOF.EQ.1.AND.IHI.EQ.0) GO TO 30

         CALL CNVF4Q( LINE, NLTE, &
                      R4,RR4,V1R4,V2R4,DVR4,NPTR4,ILO,IHI,BOUND4,&
                      ILNFLG,ILIN4,DPTFC,DPTMN,JCNVF4 )

         !CALL CPUTIM (TIM3)
         !TF4CNV = TF4CNV+TIM3-TIM2

!        IF IHI EQUALS -1 THEN END OF CONVOLUTION

         IF (IHI.EQ.-1) GO TO 30
      GO TO 20

   30 continue !CALL CPUTIM (TIM4)

      IF (JRAD.EQ.1) THEN

!     RADIATION FIELD

         XKT = 1./BETA
         VITST = V1R4-DVR4
         RDLAST = -1.
         NPTSI1 = 0
         NPTSI2 = 0

   40    NPTSI1 = NPTSI2+1

         VI = V1R4+DVR4* REAL(NPTSI1-1)
         RADVI = RADFNI(VI,DVR4,XKT,VITST,RDEL,RDLAST)

!v128         NPTSI2 = (VITST-V1R4)/DVR4+1.001
         NPTSI2 = (VITST-V1R4)/DVR4+0.001
         NPTSI2 = MIN(NPTSI2,NPTR4)

         tmpRADVI = RADVI !yma,151201
         DO 50 I = NPTSI1, NPTSI2
            !yma VI = VI+DVR4
            R4(I) = R4(I)*tmpRADVI
            tmpRADVI = tmpRADVI+RDEL
   50    CONTINUE
         if ( NLTE ) then !yma,151201
            tmpRADVI = RADVI
            DO I = NPTSI1, NPTSI2
               !yma VI = VI+DVR4
               RR4(I) = RR4(I)*tmpRADVI
               tmpRADVI = tmpRADVI+RDEL
            ENDDO
         endif

         IF (NPTSI2.LT.NPTR4) GO TO 40
      ENDIF

      !CALL CPUTIM (TIM5)
      !TF4PNL = TF4PNL+TIM5-TIM4
      !TF4 = TF4+TIM5-TIM00

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE CNVF4Q( SHRUNK, NLTE, &
                         R4,RR4,V1R4,V2R4,DVR4,NPTR4,ILO,IHI,BOUND4,&
                         ILNFLG,ILIN4,DPTFC,DPTMN,JCNVF4 )
!-----------------------------------------------------------------------
      !USE omp_lib ,ONLY: omp_get_max_threads, omp_get_num_threads, omp_get_thread_num
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, ONEPL
      USE Module_LineData         ,ONLY: LINE_SHRINK
      USE Module_ODLAY_CommonSub  ,ONLY: AVRAT

      IMPLICIT NONE !REAL*8           (V)

      type(LINE_SHRINK) ,intent(inout) :: SHRUNK
      logical           ,intent(in)    :: NLTE
      !
      real              ,intent(inout) :: R4(:)
      real              ,intent(inout) :: RR4(:)
      real(r8)          ,intent(in)    :: V1R4, V2R4
      real              ,intent(in)    :: DVR4
      integer           ,intent(in)    :: NPTR4
      integer           ,intent(in)    :: ILO
      integer           ,intent(inout) :: IHI
      real              ,intent(in)    :: BOUND4
      integer           ,intent(in)    :: ILNFLG
      integer           ,intent(inout) :: ILIN4
      real              ,intent(in)    :: DPTFC
      real              ,intent(in)    :: DPTMN
      integer           ,intent(inout) :: JCNVF4


      !--- Local variables
      !
      real  ,SAVE  :: RDVCHI
      real  ,SAVE  :: CHI(0:250) !no use
      real  ,SAVE  :: A3, B3
      real  ,SAVE  :: RECPI
      real  ,SAVE  :: ZSQBND

      real, allocatable :: tCHI(:)
      integer           :: jjj
      real              :: chival
      integer           :: lastLine
      logical           :: sppEQzero, doNLTE

      real      ,PARAMETER :: ZBND = 64.
      character ,PARAMETER :: HREJ = '0'
      character ,PARAMETER :: HNOREJ  ='1'
      integer   ,PARAMETER :: I_1 = 1
      integer   ,PARAMETER :: I_250 = 250

      CHARACTER :: FREJ(1250)
      INTEGER   :: I,          ISUBL,      IZ
      INTEGER   :: JJ,         JMAX,      JMIN
      REAL      :: ALFADI,     ALFALI,     ALFAVI,    ALFLI2
      REAL      :: ALFVI2,     BNDSQ,      CUPCON,    DPTRAT
      REAL      :: DPTRAT_R,   DVCHI,      F4BND,     F4FN
      REAL      :: F4FR,       FCNTR_FN,   FCNT_FN,   FRBND
      REAL      :: RALFVI,     REC_ALFVI2,            R_BNDSQ
      REAL      :: SIL,        SIV,        SIV_A3,    SIV_B3
      REAL      :: SPEAK,      SRL,        SRV,       SRV_A3
      REAL      :: SRV_B3
      REAL(r8)  :: VNUI,       VNULST
      REAL      :: XJJ,        XM,         XMSQ,      XNUI
      REAL      :: ZETAI,      ZETDIF,    ZVSQ



      VNULST = V2R4+BOUND4

      IF (JCNVF4.NE.1234) then

         JCNVF4 = 1234

         ! Obtain CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE:

         ! CONSTANTS FOR FOURTH FUNCTION LINE SHAPE:
         RECPI = 1./(2.*ASIN(1.))
         ZSQBND = ZBND*ZBND
         A3 = (1.+2.*ZSQBND)/(1.+ZSQBND)**2
         B3 = -1./(1.+ZSQBND)**2

      endif

      BNDSQ = BOUND4*BOUND4
      r_bndsq = 1./bndsq


      !--- START OF LOOP OVER LINES
      !
      IF (ILNFLG.EQ.2) READ(16)(FREJ(I),I=ILO,IHI)


!--------------------------------------->>> yma
      !--- Find the last line
      lastLine = ILO-1
      do I = ILO, IHI
         IF (SHRUNK%VNU(I).GE.VNULST) exit
         lastLine = I
      enddo
!--------------------------------------->>><<< yma
!$OMP  PARALLEL DO &
!$OMP& schedule(dynamic) &
!$OMP& default(none) &
!$OMP& private(ALFADI,ALFALI,ZETAI,IZ,ZETDIF,ALFAVI,RALFVI,SIL,SIV,SPEAK, &
!$OMP&         JJ,VNUI,XNUI,JMIN,JMAX,ALFLI2,ALFVI2,XJJ,F4BND,dptrat, &
!$OMP&         rec_alfvi2,siv_a3,siv_b3,sppEQzero,doNLTE,SRL,SRV,FRBND, &
!$OMP&         srv_a3,srv_b3,  &
!$OMP&         XM,XMSQ,ZVSQ,fcnt_fn,F4FN,chival,fcntr_fn,F4FR ) &
!$OMP& shared(ILO,ONEPL,AVRAT,RECPI,A3,B3, SHRUNK, &
!$OMP&        V1R4,NPTR4,DVR4,I_1,ILNFLG,HNOREJ,DPTMN,DPTFC,HREJ, &
!$OMP&        FREJ,BOUND4,VNULST,BNDSQ,NLTE,&
!$OMP&        r_bndsq,tCHI,RDVCHI,ZSQBND, &
!$OMP&        lastLine) &
!$OMP& reduction(+:R4,RR4,ILIN4)

      !yma DO 60 I = ILO, IHI
      DO 60 I = ILO, lastLine

         IF (SHRUNK%SP(I).EQ.0..AND.SHRUNK%SPP(I).EQ.0.) GO TO 60
         ALFADI = SHRUNK%EPP(I)
         ALFALI = SHRUNK%ALFA(I)
         ZETAI = ALFALI/(ALFALI+ALFADI)
         IZ = 100.*ZETAI + ONEPL
         ZETDIF = 100.*ZETAI - REAL(IZ-1)

         ALFAVI = ( AVRAT(IZ) + ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)) ) *   &
         (ALFALI+ALFADI)
         RALFVI = 1./ALFAVI
         SIL = SHRUNK%SP(I) * RECPI * ALFALI
         SIV= (ALFALI*RALFVI) * SHRUNK%SP(I)*RECPI*RALFVI

         !yma ! nlte line coupling constant
         !cupcon=shrunk%spp(i)/shrunk%sp(i)

         SPEAK = A3*(ABS(SIV))

         JJ = (SHRUNK%VNU(I)-V1R4)/DVR4+1.
         JJ = MAX(JJ,I_1)
         JJ = MIN(JJ,NPTR4)

         ! SPEAK is used for line rejection
         if (ILNFLG==0 .or. ILNFLG==1) then !IF (ILNFLG.LE.1) THEN
            FREJ(I) = HNOREJ
            ! No rejection for line-coupled lines (SHRUNK%SPP ne. 0)
            IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) .and. shrunk%spp(i).eq.0.) THEN
               FREJ(I) = HREJ
               GO TO 60
            ENDIF
         ELSE if (ILNFLG==2) then
            IF (FREJ(I).EQ.HREJ) GOTO 60
         ENDIF

         !!$OMP ATOMIC
         ILIN4 = ILIN4+1

         VNUI = SHRUNK%VNU(I)

   30    CONTINUE

         XNUI = VNUI-V1R4
         JMIN = (XNUI-BOUND4)/DVR4+2.

         !yma IF (VNUI.GE.VNULST) GO TO 70  !this modification is for openMP requirment
         IF (JMIN.GT.NPTR4) GO TO 60
         JMIN = MAX(JMIN,I_1)
         JMAX = (XNUI+BOUND4)/DVR4+1.
         IF (JMAX.LT.JMIN) GO TO 50
         JMAX = MIN(JMAX,NPTR4)
         ALFLI2 = ALFALI*ALFALI
         ALFVI2 = ALFAVI*ALFAVI
         XJJ = REAL(JMIN-1)*DVR4

         F4BND = SIL/(ALFLI2+BNDSQ)
         !yma FRBND = SRL/(ALFLI2+BNDSQ)


         !--- FOURTH FUNCTION CONVOLUTION
         !
         dptrat = shrunk%spp(i)/(shrunk%sp(i)*alfavi)
         !yma dptrat_r = shrunk%spp(i)/(shrunk%srad(i)*alfavi)
         rec_alfvi2 = 1./ALFVI2
         siv_a3 = SIV*A3
         siv_b3 = SIV*B3
         !yma srv_a3 = srv*a3
         !yma srv_b3 = srv*b3


         sppEQzero = SHRUNK%SPP(I).EQ.0. !yma127 SPP(I).EQ.0.
         doNLTE = ( NLTE .and. abs(SHRUNK%SRAD(I))>0. )

         if (doNLTE) then
            SRL = SHRUNK%SRAD(I) * RECPI * ALFALI
            SRV= (ALFALI*RALFVI) * SHRUNK%SRAD(I)*RECPI*RALFVI
            FRBND = SRL/(ALFLI2+BNDSQ)
            !yma dptrat_r = shrunk%spp(i)/(shrunk%srad(i)*alfavi)
            srv_a3 = srv*a3
            srv_b3 = srv*b3
         endif




         IF (SHRUNK%MOL(I).EQ.2.) THEN ! special treatment for co2

            DO 40 JJ = JMIN, JMAX
               XM = (XJJ-XNUI)
               XMSQ = XM*XM
               ZVSQ = XMSQ * rec_alfvi2
               fcnt_fn = (2.-(xmsq*r_bndsq))*f4bnd

               if (ZVSQ.LE.ZSQBND) then
                  F4FN = (siv_A3 + ZVSQ * siv_B3) - fcnt_fn
               else
                  F4FN = SIL/(ALFLI2+XMSQ) - fcnt_fn
               endif

               IF (sppEQzero) then !(SHRUNK%SPP(I).EQ.0.) THEN
                  !test chival = tCHI( int(RDVCHI*ABS(XM)+0.5) )

                  !!$OMP ATOMIC
                  R4(JJ) = R4(JJ) + F4FN!test *chival !CHI(ISUBL)

                  if (doNLTE) then !(NLTE==.TRUE. .and. SHRUNK%SRAD(I)/=0.) then
                     fcntr_fn= (2.-(xmsq*r_bndsq))*frbnd

                     if (ZVSQ.LE.ZSQBND) then
                        F4FR = (srv_A3 + ZVSQ * srv_B3) - fcntr_fn
                     else
                        F4FR = SRL/(ALFLI2+XMSQ) - fcntr_fn
                     endif

                     !!$OMP ATOMIC
                     RR4(JJ) = RR4(JJ) + F4FR!test *chival !CHI(ISUBL)
                  endif
               ELSE
                  F4FN = f4fn + xm*dptrat*f4fn
                  !!$OMP ATOMIC
                  R4(JJ) = R4(JJ)+F4FN
               ENDIF

               XJJ = XJJ+DVR4
   40       CONTINUE


         ELSE ! all molecules other than co2:

            DO 45 JJ = JMIN, JMAX
               XM = (XJJ-XNUI)
               XMSQ = XM*XM
               ZVSQ = XMSQ * rec_alfvi2

               if (ZVSQ.LE.ZSQBND) then
                  F4FN = (siv_A3 + ZVSQ * siv_B3) - F4BND
               else
                  F4FN = SIL/(ALFLI2+XMSQ)-F4BND
               endif

               IF (sppEQzero) then !(SHRUNK%SPP(I).EQ.0.) THEN !yma; This is different from that in oprop.f90 ?
                  !test chival = tCHI( int(RDVCHI*ABS(XM)+0.5) )

                  !!$OMP ATOMIC
                  R4(JJ) = R4(JJ) + F4FN!test *chival !CHI(ISUBL)

                  if (doNLTE) then !(NLTE==.TRUE. .and. SHRUNK%SRAD(I)/=0.) then
                     if (ZVSQ.LE.ZSQBND) then
                        F4FR = (srv_A3 + ZVSQ * srv_B3) - FRBND
                     else
                        F4FR = SRL/(ALFLI2+XMSQ)-FRBND
                     endif
                     !!$OMP ATOMIC
                     RR4(JJ) = RR4(JJ)+F4FR!test *chival !CHI(ISUBL)
                  endif
               ELSE
                  F4FN = f4fn + xm*dptrat*f4fn
                  !!$OMP ATOMIC
                  R4(JJ) = R4(JJ)+F4FN
               ENDIF

               XJJ = XJJ+DVR4
   45       CONTINUE

         ENDIF

   50    IF (VNUI.GT.0..AND.VNUI.LE.25.) THEN

            ! THE CALCULATION FOR NEGATIVE VNU(I) IS FOR VAN VLECK WEISSKOPF
            VNUI = -SHRUNK%VNU(I)
            SHRUNK%SPP(I) = -SHRUNK%SPP(I)
            GO TO 30

         ENDIF

   60 END DO
!$OMP END PARALLEL DO


      IF (ILNFLG.EQ.1) WRITE(16)(FREJ(I),I=ILO,IHI)

      if (lastLine<IHI) IHI = -1

      RETURN

!      ! IF END OF CONVOLUTION, SET IHI=-1 AND RETURN
!   70 CONTINUE
!      IF (ILNFLG.EQ.1) WRITE(16)(FREJ(I),I=ILO,IHI)
!      IHI = -1
!      RETURN

      END SUBROUTINE

!-----------------------------------------------------------------------
!     SET UP CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE
!     POLYNOMIAL MATCHED TO AN EXPONENTIAL AT X0 = 10 CM-1
!-----------------------------------------------------------------------
      subroutine chi_fn( chi,dvchi )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real ,intent(out) :: chi(0:250)
      real ,intent(in)  :: dvchi

      !--- Local variables
      INTEGER :: ISUBL
      REAL    :: ASUBL, BSUBL, C2, C4
      REAL    :: C6,    F,     FI, X0
      REAL    :: Y0,    Y1,    Y2, Z0
      REAL    :: Z1,    Z2

      DATA ASUBL / 0.800 /,BSUBL / 10.0 /


!     0 - 25 cm-1 on 0.1 cm-1 grid

      X0 = 10.
      Y0 = asubl*EXP(-x0/bsubl)
      F  = 1./bsubl
      Y1 = -F*Y0
      Y2 = Y1*((BSUBL-1)/X0-F)
      Z0 = (Y0-1)/X0**2
      Z1 = Y1/(2*X0)
      Z2 = Y2/2.
      C6 = (Z0-Z1+(Z2-Z1)/4.)/X0**4
      C4 = (Z1-Z0)/X0**2-2.*X0**2*C6
      C2 = Z0-X0**2*C4-X0**4*C6

      DO 10 ISUBL = 0, 250
         FI = DVCHI* REAL(ISUBL)
         IF (FI.LT.X0) THEN
            CHI(ISUBL) = 1.+C2*FI**2+C4*FI**4+C6*FI**6
         ELSE
            CHI(ISUBL) = asubl*EXP(-fi/bsubl)
         ENDIF

!**%%$$

         chi(isubl) = 1.

   10 END DO

      return
      END SUBROUTINE


!-----------------------------------------------------------------------
!
!     SUBROUTINE RDLIN4Q INPUTS THE LINE DATA FROM LNFIL4
!
!-----------------------------------------------------------------------
      SUBROUTINE RDLN4Q( LINE,IEOF, &
                         LNFIL4,NEGEPP_FLAG,VLO,VHI,ILO,IHI,ILIN4T )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8
      USE Module_LineData     ,ONLY: LINE_SHRINK

      IMPLICIT NONE !REAL*8           (V)
     TYPE :: INPUT_HEADER
        SEQUENCE
        REAL(8)    :: VMIN, VMAX
        INTEGER(4) :: NREC, NWDS
     END TYPE INPUT_HEADER


      integer           ,intent(inout) :: IEOF
      type(LINE_SHRINK) ,intent(out)   :: LINE
      !
      integer           ,intent(in)    :: LNFIL4
      integer           ,intent(in)    :: NEGEPP_FLAG
      real(r8)          ,intent(in)    :: VLO,VHI
      integer           ,intent(out)   :: ILO,IHI
      integer           ,intent(inout) :: ILIN4T

      !--- Local variables
      INTEGER  :: LEOF,  NREC
      INTEGER  :: NWDS
      REAL     :: DUM(2)
      REAL(r8) :: VMAX,  VMIN
      type(input_header) :: rdlnpnl
      CHARACTER*8      :: HLINID(10)
      CHARACTER*1      :: CNEGEPP(8)
      real*4 :: xspace(4096)
      integer*4 ::n_resetepp(64), n_negepp(64)

      IF (IEOF.EQ.-1) THEN

!     BUFFER PAST FILE HEADER

         REWIND LNFIL4
         ILIN4T = 0
         read (LNFIL4) HLINID
         ! print *,'HLINID', HLINID
         READ (HLINID(7),'(8a1)') CNEGEPP
         ! print *,'CNEGEPP', CNEGEPP
         ! print *, 'negepp_flag', negepp_flag
         IF (CNEGEPP(8).eq.'^') THEN
            read (LNFIL4) n_negepp,n_resetepp,xspace
         ENDIF

         ! CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
         ! IF (LEOF.EQ.0) STOP 'RDLIN4; TAPE9 EMPTY'
         ! IF (LEOF.EQ.-99) THEN
         !    IEOF = 2

         !    RETURN

         ! ENDIF
         ! if (negepp_flag.eq.1) CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
      ENDIF
      IEOF = 0
      ILO = 1
      IHI = 0


   10 READ(LNFIL4,END=20,ERR=20) VMIN,VMAX,NREC,NWDS !yma "ERR=20" follows RDLN4() from oprop.f90
!   10 CALL BUFIN (LNFIL4,LEOF,LINPNL(1),NPHDRL)
!      IF (LEOF.EQ.0) GO TO 20
      ILIN4T = ILIN4T+NREC
      IF (VMAX.LT.VLO) THEN
         CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
         GO TO 10
      ELSE
         CALL BUFIN (LNFIL4,LEOF,LINE,NWDS)
      ENDIF
      IHI = NREC
      IF (LINE%VNU(NREC).GT.VHI) GO TO 20

      RETURN

   20 IEOF = 1

  950 FORMAT (8a1)

      RETURN
      END SUBROUTINE



!-----------------------------------------------------------------------
!
!     SUBROUTINE LINF4 READS THE LINES AND SHRINKS THE LINES FOR LBLF4
!
!-----------------------------------------------------------------------
      SUBROUTINE LINF4Q( V1L4,V2L4, &
                         LINFIL,LNFIL4,NEGEPP_FLAG,DVR4,V1,V2,&
                         PAVE,TAVE,molNames,W,WTOT,WKI,ISOTPL_FLAG,&
                         DPTFAC,DPTMIN,NLTE,speciesBroad,NUMSTATE,RATSTATE,&
                         NPHDRL,NLNGTH,NOPR, &
                         numNarrowLines, narrowWidth )
!-----------------------------------------------------------------------
      USE Module_ConstParam       ,ONLY: r8=>kind_r8, RADCN2, ONEPL, &
                                         MXMOL,MXISOTPL,MAX_ISO,MaxISOTPL_smass,MAXSTATE, &
                                         P0=>Press1013, TEMP0=>Temp296, &
                                         molNum, xsMolNum
      USE Module_FileIO        ,ONLY: NWDL
      USE Module_LineData         ,ONLY: MXBRDMOL, NLINEREC, LINE_DATA, LINE_SHRINK, INPUT_BLOCK
      USE Module_ODLAY_CommonSub  ,ONLY: MOLEC, RADFN, line_exception, AVRAT
      USE Module_Config           ,ONLY: IPR, ioFiles
      USE Module_XSect            ,ONLY: xsTbl

      IMPLICIT NONE !REAL*8           (V)

      real(r8)    ,intent(in)  :: V1L4
      real(r8)    ,intent(in)  :: V2L4
      !
      integer     ,intent(in)  :: LINFIL
      integer     ,intent(in)  :: LNFIL4
      integer     ,intent(out) :: NEGEPP_FLAG
      real        ,intent(in)  :: DVR4
      real(r8)    ,intent(in)  :: V1, V2
      real        ,intent(in)  :: PAVE
      real        ,intent(in)  :: TAVE
      character(*),intent(in)  :: molNames(:)
      real        ,intent(in)  :: W(:) !(60)
      real        ,intent(in)  :: WTOT
      real        ,intent(in)  :: WKI(:,:)         !(MXMOL,MXISOTPL)
      integer     ,intent(in)  :: ISOTPL_FLAG(:,:) !(MXMOL,MXISOTPL)
      real        ,intent(in)  :: DPTFAC
      real        ,intent(in)  :: DPTMIN
      logical     ,intent(in)  :: NLTE
      logical     ,intent(in)  :: speciesBroad
      integer     ,intent(in)  :: NUMSTATE(:)      !(MXMOL)
      real        ,intent(in)  :: RATSTATE(:,:)    !(MAXSTATE*MAX_ISO,MXMOL)
      integer     ,intent(in)  :: NPHDRL
      integer     ,intent(in)  :: NLNGTH
      integer     ,intent(in)  :: NOPR
      integer     ,intent(inout):: numNarrowLines
      real        ,intent(in)  :: narrowWidth




      !--- Local variables
      !
      character(*),parameter :: routineName='LINF4Q'

      integer  :: LNGTH4
      integer  :: IOUT(250)
      integer  :: ILO,IHI,ILAST
      integer  :: IST
      integer  :: LIMIN
      real     :: DPTFC, DPTMN
      real(r8) :: VHI, VLO


      character(8) ,PARAMETER :: h_linf4=' linf4  '
      real         ,PARAMETER :: TEMPLC(4)=[ 200.0,250.0,296.0,340.0 ] !TEMPERATURES FOR LINE COUPLING COEFFICIENTS
      integer      ,PARAMETER :: I_100=100
      integer      ,PARAMETER :: I_1000=1000

      TYPE(INPUT_BLOCK) :: rawLines !double convolution needs original line data.
      TYPE(INPUT_BLOCK) :: narrowLines
      TYPE(LINE_DATA)   :: BUFR
      TYPE(LINE_SHRINK) :: SHRUNK
      integer           :: MEFDP(64)
      integer           :: jrad4
      Real(r8)          :: ALFSUM
      real              :: SCOR(MXMOL,MXISOTPL)
      real              :: RHOSLF(MXMOL)
      real              :: ALFD1(MXMOL,MXISOTPL)
      real              :: A(4),B(4)
      real*4            :: amol
      real              :: ALFV,ALFL,ALFAD,ZETA,FZETA,ZETDIF
      integer           :: nrw,IZ
      integer           :: narrowLineFile
      logical           :: checkNarrowLine

      INTEGER  :: I,            IEOF,      IFLAG,          IJ
      INTEGER  :: ILC,          ILINHI,    ILINLO,         INDLOW
      INTEGER  :: INDUPP,       ISO
      INTEGER  :: J, M
      INTEGER  :: MFULL
      INTEGER  :: NLIN,         NLINS,     NLOW,           NUPP
      REAL     :: ALFA0I,    ALFA_TMP(MXBRDMOL)
      REAL     :: BETA
      REAL     :: BETA0,        BETACR,    DELTA
      REAL     :: DELTMP,       FNLTE,     FREQ,           GAMMA1
      REAL     :: GAMMA2,       GI,        HWHMSI,         PAVP0
      REAL     :: PAVP2,        RECTLC,         RHORAT
      REAL     :: RLOW,           RUPP
      REAL     :: SLOPEG
      REAL     :: SLOPEY,       SUI,       SUMS
      REAL     :: TMPCOR,   TRATIO
      REAL     :: TMPCOR_ARR(MXBRDMOL),    TMPDIF
      REAL     :: YI

      integer :: im, imol, is, ind
      real    :: v1x,v2x
      integer :: moNum( size(molNames) )
      integer :: xsNum( size(molNames) )



      !--- Local common blocks
      !
      REAL*8   :: VNULO,VNUHI
      INTEGER  :: JLIN,NLNGT4,lstdum
      real     :: RCDHDR(2)
      COMMON /TPANEL_linf4q/ VNULO,VNUHI,JLIN,NLNGT4,lstdum
      EQUIVALENCE (RCDHDR(1),VNULO)

      REAL    :: SD,AD,EPD,SPPD
      REAL*8  :: VD
      INTEGER :: MOLD,ILS2D, IWD3(2)
      COMMON /NGT4_linf4q/ VD,SD,AD,EPD,MOLD,SPPD,ILS2D
      EQUIVALENCE (IWD3(1),VD)


      real*8    :: nrwVLO,nrwVHI
      integer*4 :: nrwLINES,nrwNWDS, headLSTW, nrwHead(2), lenHead=6
      COMMON /common_narrowlineHeader/nrwVLO,nrwVHI,nrwLINES,nrwNWDS,headLSTW
      EQUIVALENCE (nrwHead(1),nrwVLO)
      integer :: lenRec=9750 !39*250 4-bytes words


      MEFDP(:) = 0 !DATA MEFDP / 64*0 /

      jrad4 = 0 !data jrad4 /0/

      !CALL CPUTIM (TIMEL0)

      ILS2D = -654321
      NLNGT4 = NWDL(IWD3,ILS2D)*1250
      LNGTH4 = NLNGT4
      PAVP0 = PAVE/P0
      PAVP2 = PAVP0*PAVP0
      DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)
      DPTFC = DPTFAC
      LIMIN = 1000


      checkNarrowLine = numNarrowLines<0
! Betsy Berry 7/21/21 disabled the narrow line calculations to compare to LBL (which doesn't have that capability)
      checkNarrowline=.FALSE.
      numNarrowLines = 0


      do im=1,size(molNames)
         moNum(im) = molNum( molNames(im) )
         xsNum(im) = xsMolNum( molNames(im) )
      enddo


      !CALL CPUTIM(TPAT0)
      CALL MOLEC( 1,SCOR,RHOSLF,ALFD1, &
                  PAVE,TAVE,W,WTOT,molNames )
      !CALL CPUTIM(TPAT1)
      !TMOLN4 = TMOLN4 + TPAT1-TPAT0

!yma       TIMR = 0.
!yma       TIMS = 0.
      SUMS = 0.
      ILAST = 0
      ILINLO = 0
      ILINHI = 0
      ILO = 1
      IST = 1
      NLINS = 0
      NLIN = 0

      VLO = V1L4
      VHI = V2L4


      !CALL CPUTIM(TPAT0)
      !--- Transfer the line file header to scratch files
      !
      call lnfilhd_4( linfil,lnfil4,v1,v2, NEGEPP_FLAG )

      if (checkNarrowLine) then
         inquire( file=trim(ioFiles%narrowLineFile), number=narrowLineFile )
         call lnfilhd_4( linfil,narrowLineFile,v1,v2, NEGEPP_FLAG )
      endif


      !CALL CPUTIM(TPAT1)
      !TBUFFR = TBUFFR + TPAT1-TPAT0

!       TEMPERATURE CORRECTION TO INTENSITY
!       TEMPERATURE AND PRESSURE CORRECTION TO HALF-WIDTH

      TRATIO = TAVE/TEMP0
      RHORAT = (PAVE/P0)*(TEMP0/TAVE)

      BETA = RADCN2/TAVE
      BETA0 = RADCN2/TEMP0
      BETACR = BETA-BETA0
      DELTMP = ABS(TAVE-TEMP0)
      !CALL CPUTIM(TPAT0)
      CALL MOLEC( 2,SCOR,RHOSLF,ALFD1, &
                  PAVE,TAVE,W,WTOT,molNames )
      !CALL CPUTIM(TPAT1)
      !TMOLN4 = TMOLN4 + TPAT1-TPAT0

!     FIND CORRECT TEMPERATURE AND INTERPOLATE FOR Y AND G

      DO 10 ILC = 1, 4
         IF (TAVE.LE.TEMPLC(ILC)) GO TO 20
   10 END DO
   20 IF (ILC.EQ.1) ILC = 2
      IF (ILC.EQ.5) ILC = 4
      RECTLC = 1.0/(TEMPLC(ILC)-TEMPLC(ILC-1))
      TMPDIF = TAVE-TEMPLC(ILC)

      nrw = 0
      IJ = 0
   30 continue
         !CALL CPUTIM (TIM0)
         CALL RDLNFL (IEOF,ILINLO,ILINHI, BUFR, LINFIL,NLNGTH,VLO,NOPR,IOUT, rawLines)
         !CALL CPUTIM (TIM1)
         !TIMR = TIMR+TIM1-TIM0

         IF (IEOF.GE.1) GO TO 60 !previous block was the last block in the file. eof reached. write out the SHRUNK and narrowlines

         DO 50 J = ILINLO, ILINHI
            YI = 0.
            GI = 0.
            GAMMA1 = 0.
            GAMMA2 = 0.
            I = IOUT(J)
            IFLAG = BUFR%IFLG(I)
            IF (I.LE.0) GO TO 50

            MFULL = BUFR%MOL(I)
            M = MOD(BUFR%MOL(I),I_100)

            ! ISO=(MOD(BUFR%MOL(I),I_1000)-M)/100   IS PROGRAMMED AS:

            ISO = MOD(BUFR%MOL(I),I_1000)/100
            if (iso.eq.0) iso=10


            !--- Check if this line is from a selected molecule and when
            ! both line data and xs data are available if this line prefer xs over line,
            ! if so, skip this line.
            !
            if ( .not.any(moNum==m) ) then !This line is not from selected molecules
               GO TO 50
            else  !line is from a selected mol, check if it is also a xs molecule.

               ind = minloc( abs(moNum-m), 1 )
               if ( xsNum(ind)>0 ) then !this molecule is also a xs molecule and has xs data in V1~V2

                  !This line is from a xs molecule, unless the lineOrXs marker ==1, the absorption
                  !is handled using x-section data. This assums xs data covers all possible intervals.
                  if( all( xsTbl%lineOrXs(:,xsNum(ind)) /=1) ) GO TO 50

                  !At least one xs spectral interval is with lineOrXs marker==1
                  do is = 1,xsTbl%numSpect( xsNum(ind) )
                     v1x = xsTbl%V1FX(is,xsNum(ind))
                     v2x = xsTbl%V2FX(is,xsNum(ind))
                     if ( v1x < BUFR%VNU(I).and.BUFR%VNU(I) < v2x ) then
                        if( xsTbl%lineOrXs(is, xsNum(ind)) /=1 ) then
                           GO TO 50 !This interval handeled using xs section
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
               !      if ( v1x < BUFR%VNU(I).and.BUFR%VNU(I) < v2x ) then
               !         GO TO 50 !this line lies in a xs sub interval and is treaded by xs
               !      endif
               !   enddo
               !endif
            endif


            ! check if lines are within allowed molecular and isotopic limits
            if (m.gt.mxmol .or. m.lt. 1) then
               call line_exception (1,ipr,h_linf4,m,mxmol,iso,MaxISOTPL_smass)
            go to 50
            else if (iso .gt. MaxISOTPL_smass(m)) then
               call line_exception (2,ipr,h_linf4,m,mxmol,iso,MaxISOTPL_smass)
            go to 50
            endif

            IF (ISOTPL_FLAG(M,ISO).EQ.0) THEN
               SUI = BUFR%SP(I)*W(M)
            ELSE
               SUI = BUFR%SP(I)*WKI(M,ISO)
            ENDIF

            IF (SUI.EQ.0.) GO TO 50
            IF (BUFR%VNU(I).LT.VLO) GO TO 50
            IJ = IJ+1

            ! Y'S AND G'S ARE STORED IN I+1 POSTION OF VNU,S,ALFA0,EPP...
            ! A(1-4),  B(1-4) CORRESPOND TO TEMPERATURES TEMPLC(1-4) ABOVE

            IF (IFLAG.EQ.1.OR.IFLAG.EQ.3) THEN
               A(1) = BUFR%VNU(I+1)
               B(1) = BUFR%SP(I+1)
               A(2) = BUFR%ALFA(I+1)
               B(2) = BUFR%EPP(I+1)
               A(3) = TRANSFER(BUFR%MOL(I+1),amol )  ! real representation of mol
               B(3) = BUFR%HWHM(I+1)
               A(4) = BUFR%TMPALF(I+1)
               B(4) = BUFR%PSHIFT(I+1)

               ! CALCULATE SLOPE AND EVALUATE

               SLOPEY = (A(ILC)-A(ILC-1))*RECTLC
               SLOPEG = (B(ILC)-B(ILC-1))*RECTLC
               IF (IFLAG.EQ.1) THEN
                  YI = A(ILC)+SLOPEY*TMPDIF
                  GI = B(ILC)+SLOPEG*TMPDIF
               ELSE
                  GAMMA1 = A(ILC)+SLOPEY*TMPDIF
                  GAMMA2 = B(ILC)+SLOPEG*TMPDIF
               ENDIF
            ENDIF

            ! IFLAG = 2 IS RESERVED FOR LINE COUPLING COEFFICIENTS ASSOCIATED
            !           WITH AN EXACT TREATMENT (NUMERICAL DIAGONALIZATION)
            ! IFLAG = 3 TREATS LINE COUPLING IN TERMS OF REDUCED WIDTHS

            SHRUNK%VNU(IJ) = BUFR%VNU(I) + RHORAT * BUFR%PSHIFT(I)
            SHRUNK%ALFA(IJ) = BUFR%ALFA(I)
            SHRUNK%EPP(IJ) = BUFR%EPP(I)
            SHRUNK%MOL(IJ) = M

            if(sum(bufr%brd_mol_flg(:,i)).gt.0.AND. speciesBroad) then
               shrunk%vnu(ij) = shrunk%vnu(ij)+sum(rhoslf(1:mxbrdmol)*bufr%brd_mol_flg(:,i) &
     &              *(bufr%brd_mol_shft(:,i)-bufr%pshift(i)))
            endif


            IF (jrad4.EQ.1) SUI = SUI*SHRUNK%VNU(IJ)

            IF (SHRUNK%VNU(IJ).EQ.0.) SUI = 2.*SUI

            !   TREAT TRANSITIONS WITH UNKNOWN EPP AS SPECIAL CASE
            !
            !>> an epp value between -0.9999 and 0.  cm-1 is taken as valid
            !
            !>> an epp value of -1. is assumed set by hitran indicating an unknown
            !   value: no temperature correction is performed
            !
            !>> for an epp value of less than -1., it is assumed that value has
            !   been provided as a reasonable value to be used for purposes of
            !   temperature correction.  epp is set positive

            if (shrunk%epp(ij).le.-1.001) shrunk%epp(ij) = abs(shrunk%epp(ij))

            if (shrunk%epp(ij).le.-0.999) MEFDP(M) = MEFDP(M)+1

            ! temperature correction:

            if (shrunk%epp(ij) .gt. -0.999) then
               SUI = SUI*SCOR(m,iso)* EXP(-SHRUNK%EPP(ij)*BETACR)*  &
               (1.+EXP(-SHRUNK%VNU(ij)* BETA))
            endif

            SUMS = SUMS+SUI

            ! TEMPERATURE CORRECTION OF THE HALFWIDTH
            ! SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN

            TMPCOR = TRATIO**BUFR%TMPALF(I)
            ALFA0I = SHRUNK%ALFA(IJ)*TMPCOR
            HWHMSI = BUFR%HWHM(I)*TMPCOR
            SHRUNK%ALFA(IJ) = ALFA0I*(RHORAT-RHOSLF(m))+HWHMSI*RHOSLF(m)

            if(sum(bufr%brd_mol_flg(:,i)).gt.0 .AND. speciesBroad) then
               tmpcor_arr = tratio**bufr%brd_mol_tmp(:,i)
               alfa_tmp = bufr%brd_mol_hw(:,i)*tmpcor_arr
               alfsum = sum(rhoslf(1:mxbrdmol)*bufr%brd_mol_flg(:,i)*alfa_tmp)
               shrunk%alfa(ij) = (rhorat-sum(rhoslf(1:mxbrdmol)* &
     &              bufr%brd_mol_flg(:,i)))*alfa0i + alfsum
               if(bufr%brd_mol_flg(m,i).eq.0)   &
     &              shrunk%alfa(ij) = shrunk%alfa(ij) + rhoslf(m)*(hwhmsi-alfa0i)
            end if

! mji - Revise code to skip broadening for incoming lines with flag = -1
!      if(ibrd.gt.0) then
      if(speciesBroad) then
         if (bufr%brd_mol_flg(m,i).eq.-1) then
            shrunk%alfa(ij) = hwhmsi
         end if
      end if

            IF (IFLAG.EQ.3) SHRUNK%ALFA(IJ) = SHRUNK%ALFA(IJ)* &
     &         (1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)

            !--- For double convolution, pick out the narrow lines and
            ! save them into a scratch file for later use.
            if (checkNarrowLine) then
               stop

               ALFL  = SHRUNK%ALFA(IJ)
               ALFAD = SHRUNK%VNU(ij)*ALFD1(m,iso)
               ZETA  = ALFL/(ALFL+ALFAD)
               !ZETAI(I) = ZETA
               FZETA = 100.*ZETA
               IZ = FZETA + ONEPL
               !IZETA(I) = IZ
               ZETDIF = FZETA - REAL(IZ-1)

               ALFV = (AVRAT(IZ)+ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)))*(ALFL+ALFAD)

               !--- Check if this is a narrow line
               if (ALFV < narrowWidth) then

                  nrw = nrw+1
                  narrowLines%VNU(             nrw) = rawLines%VNU(             I)
                  narrowLines%SP(              nrw) = rawLines%SP(              I)
                  narrowLines%ALFA(            nrw) = rawLines%ALFA(            I)
                  narrowLines%EPP(             nrw) = rawLines%EPP(             I)
                  narrowLines%MOL(             nrw) = rawLines%MOL(             I)
                  narrowLines%HWHM(            nrw) = rawLines%HWHM(            I)
                  narrowLines%TMPALF(          nrw) = rawLines%TMPALF(          I)
                  narrowLines%PSHIFT(          nrw) = rawLines%PSHIFT(          I)
                  narrowLines%IFLG(            nrw) = rawLines%IFLG(            I)
                  narrowLines%BRD_MOL_FLG_IN(:,nrw) = rawLines%BRD_MOL_FLG_IN(:,I)
                  narrowLines%BRD_MOL_DAT(   :,nrw) = rawLines%BRD_MOL_DAT(   :,I)
                  narrowLines%SPEED_DEP(       nrw) = rawLines%SPEED_DEP(       I)

                              nrwVHI = narrowLines%VNU(nrw)
                  if (nrw==1) nrwVLO = narrowLines%VNU(nrw)

                  if (rawLines%IFLG(I).EQ.1 .OR. rawLines%IFLG(I).EQ.3) then
                     nrw = nrw+1
                     narrowLines%VNU(             nrw) = rawLines%VNU(             I+1)
                     narrowLines%SP(              nrw) = rawLines%SP(              I+1)
                     narrowLines%ALFA(            nrw) = rawLines%ALFA(            I+1)
                     narrowLines%EPP(             nrw) = rawLines%EPP(             I+1)
                     narrowLines%MOL(             nrw) = rawLines%MOL(             I+1)
                     narrowLines%HWHM(            nrw) = rawLines%HWHM(            I+1)
                     narrowLines%TMPALF(          nrw) = rawLines%TMPALF(          I+1)
                     narrowLines%PSHIFT(          nrw) = rawLines%PSHIFT(          I+1)
                     narrowLines%IFLG(            nrw) = rawLines%IFLG(            I+1)
                     !narrowLines%BRD_MOL_FLG_IN(:,nrw) = rawLines%BRD_MOL_FLG_IN(:,I+1)
                     !narrowLines%BRD_MOL_DAT   (:,nrw) = rawLines%BRD_MOL_DAT   (:,I+1)
                     !narrowLines%SPEED_DEP       (nrw) = rawLines%SPEED_DEP       (I+1)
                  endif

                  if (nrw >= NLINEREC-2) then !The next line may be coupled line and has two records. So "nrw+1>="
                     nrwLINES = nrw
                     nrwNWDS = lenRec
                     CALL BUFOUT_sgl (narrowLineFile,nrwHead(1),lenHead)
                     CALL BUFOUT_sgl (narrowLineFile,narrowLines,lenRec)
                     numNarrowLines = numNarrowLines + nrw
                     nrw = 0
                  endif


                  !--- skip this narrow line. It will be treated separately in a second run
                  IJ=IJ-1
                  GO TO 50

               endif !if (ALFV < narrowWidth) then
           endif !if (checkNarrowLine)


            SHRUNK%EPP(IJ) = SHRUNK%VNU(IJ)*ALFD1(m,iso)
            NLIN = NLIN+1
            SHRUNK%SP(IJ) = SUI*(1.+GI*PAVP2)
            SHRUNK%SPP(IJ) = SUI*YI*PAVP0
            SHRUNK%SRAD(IJ) = 0.0

            ! For NLTE lines:

            FREQ=SHRUNK%VNU(IJ)
            RLOW=1.0
            RUPP=1.0


            ! PICK OUT MOLECULAR TYPES WITH VIBRATIONAL STATES

!yma             IF(NUMSTATE(M).GT.0) THEN
            IF( NLTE .and. NUMSTATE(M).GT.0) THEN
               NLOW=MOD(MFULL/1000,100)
               NUPP=MFULL/100000
               ! DELTA=EXP(-FREQ/XKT)
               ! xkt=tave/radcn2=1/beta
               DELTA=EXP(-FREQ*BETA)
               IF (NLOW.GT.NUMSTATE(M)) STOP 'NLOW GT NUMSTATE IN LINF4Q'
               INDLOW=NLOW + (ISO-1)*MAXSTATE
               IF (NLOW.GT.0) RLOW=RATSTATE(INDLOW,M)
               IF (NUPP.GT.NUMSTATE(M)) STOP 'NUPP GT NUMSTATE IN LINF4Q'
               INDUPP=NUPP + (ISO-1)*MAXSTATE
               IF (NUPP.GT.0) RUPP=RATSTATE(INDUPP,M)
               FNLTE=SHRUNK%SP(IJ)/(1.0-DELTA)
               SHRUNK%SP(IJ)=FNLTE*(RLOW-RUPP*DELTA)
               SHRUNK%SRAD(IJ)=FNLTE*(RLOW-RUPP)
            ELSE
               RLOW=0.
               RUPP=0.
            END IF

            ! RLOW AND RUPP NOW SET

            IF (SHRUNK%VNU(IJ).GT.VHI) THEN !line VNU exceed the interval limit V2L4, No further lines needed. write out the SHURNK and narrowlines.
               IEOF = 1  !IEOF=1 flag to return from this routine.
               GO TO 60
            ENDIF

   50    END DO
!??? what if IJ>LIMIN in middle of loop 50?

         IF (IJ.LT.LIMIN.AND.IEOF.EQ.0) THEN
            !CALL CPUTIM (TIM2)
            !TIMS = TIMS+TIM2-TIM1
            GO TO 30 !SHRUNK buffer not full yet, get more un-shrunk lines.
         ENDIF

   60    continue !Three ways to reach here: 1)EOF goto 60; 2)>VHI goto 60; 3)IJ==LIMIN
         !
         ! SHRUNK: sssssssss---------------------
         !         |       ||          |
         !         1       |ILO        IHI=IJ (before calling shrink)
         !                 |
         !                IHI (after calling shrink)
         !
         ! "s" means shrunk lines, "-" un-shunk lines
         !

         !CALL CPUTIM (TIM2)
         IHI = IJ
         !TIMS = TIMS+TIM2-TIM1

         !CALL CPUTIM(TPAT0)
         CALL SHRINQ( SHRUNK, ILO,IHI,DVR4 )
         !CALL CPUTIM(TPAT1)
         !TSHRNK = TSHRNK + TPAT1-TPAT0
         IJ = ILO-1
         IF (IHI.LT.LIMIN.AND.IEOF.EQ.0) GO TO 30 !SHUNK is full of shrunk lines, write the 1000 shrunk lines.

         !--- For double convolution, when checking for narrow lines at the first
         ! pass, there may be cases when all lines are narrow lines, so nothing to shrink.
         if (IHI>0) then
            VNULO = SHRUNK%VNU(1)
            VNUHI = SHRUNK%VNU(IHI)
            JLIN = IHI

            IF (JLIN.GT.0) THEN
               !CALL CPUTIM(TPAT0)
               CALL BUFOUT (LNFIL4,RCDHDR(1),NPHDRL)
               CALL BUFOUT (LNFIL4,SHRUNK,NLNGT4)
               !CALL CPUTIM(TPAT1)
               !TBUFFR = TBUFFR + TPAT1-TPAT0
            ENDIF
            NLINS = NLINS+IHI-IST+1
         endif

         IF (IEOF.EQ.1) GO TO 70
         IJ = 0
         ILO = 1
      GO TO 30 !Haven't reach the V2L4, empty the SHRUNK buffer and loop again.
   70 CONTINUE


      if ( nrw>0 ) then !There are lines remain in narrowlines structure, write them out.
         nrwLINES = nrw
         nrwNWDS = lenRec
         CALL BUFOUT_sgl (narrowLineFile,nrwHead(1),lenHead)
         CALL BUFOUT_sgl (narrowLineFile,narrowLines,lenRec)
         numNarrowLines = numNarrowLines + nrw
         nrw = 0
      endif


      DO 80 M = 1, size(molNames)
         imol = molNum( molNames(m) )
         IF (MEFDP(imol).GT.0) WRITE (IPR,905) MEFDP(imol),imol
   80 END DO
      !CALL CPUTIM (TIMEL1)
      !TIME = TIMEL1-TIMEL0
      !IF (NOPR.EQ.0) THEN
      !   WRITE (IPR,910) TIME,TIMR,TIMS,NLIN,NLINS
      !   L4TIM=TIME
      !   L4TMR=TIMR
      !   L4TMS=TIMS
      !   L4NLN=NLIN
      !   L4NLS=NLINS
      !   LOTHER = TSHRNK+TBUFFR+TMOLN4
      !ENDIF
      RETURN

  900 FORMAT ('0  *****  LINF4 - VNU LIMITS DO NOT INTERSECT WITH ',    &
     &        'LINFIL - LINF4 NOT USED *****',/,'   VNU = ',F10.3,      &
     &        ' - ',F10.3,' CM-1     LINFIL = ',F10.3,' - ',F10.3,      &
     &        ' CM-1')
  905 FORMAT ('0*************************',I5,' STRENGTHS FOR',         &
     &        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,            &
     &        ' SET TO ZERO')
  910 FORMAT ('0',20X,'TIME',11X,'READ',9X,'SHRINQ',6X,'NO. LINES',3X,  &
     &        'AFTER SHRINQ',/,2X,'LINF4 ',2X,3F15.3,2I15)

      END SUBROUTINE

!-----------------------------------------------------------------------
!     this subroutine buffers past the file header for LINF4
!-----------------------------------------------------------------------
      subroutine lnfilhd_4( linfil,linfil4,v1,v2, NEGEPP_FLAG)
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_FileIO   ,ONLY: ENDFIL_4 !yma 151225
      USE Module_Config      ,ONLY: IPR

      IMPLICIT NONE !REAL*8           (V)

      integer     ,intent(in)  :: LINFIL
      integer     ,intent(in)  :: LINFIL4
      real(r8)    ,intent(in)  :: V1,V2
      integer     ,intent(out) :: NEGEPP_FLAG


      !--- Local Variables
      !
      integer*4    :: N_NEGEPP(64)
      integer*4    :: N_RESETEPP(64)
      real*4       :: XSPACE(4096)
      integer*4    :: ILINLC,   ILINNL
      integer*4    :: IREC,     IRECTL
      integer*4    :: LINCNT
      integer*4    :: LINMOL
      integer*4    :: MCNTLC(64)
      integer*4    :: MCNTNL(64)
      integer*4    :: MOLCNT(64)
      integer*4    :: lnfil,lnfil4
      real*4       :: FLINHI,   FLINLO
      real*4       :: SUMSTR(64)
      character*8  :: BMOLID(64)
      character*8  :: HID1(2)
      character*8  :: HLINID(10)
      character*1  :: CNEGEPP(8)



      lnfil = linfil
      lnfil4= linfil4


      REWIND lnfil

      read (lnfil,end=777)    HLINID,BMOLID,MOLCNT,MCNTLC,              &
     &                MCNTNL,SUMSTR,LINMOL,FLINLO,FLINHI,               &
     &                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1

      READ (HLINID(7),950) CNEGEPP
      IF (CNEGEPP(8).eq.'^') THEN
         read (lnfil) n_negepp,n_resetepp,xspace
      endif

      go to 5

  777 STOP 'Linf4: LINFIL DOES NOT EXIST'

    5 continue


      IF (V1.GT.FLINHI.OR.V2.LT.FLINLO) THEN
         CALL ENDFIL_4 (LNFIL4)
         WRITE (IPR,900) V1,V2,FLINLO,FLINHI
         RETURN
      ENDIF

      write (lnfil4) HLINID,BMOLID,MOLCNT,MCNTLC, MCNTNL,SUMSTR,     &
      LINMOL,FLINLO,FLINHI, LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1

      IF (CNEGEPP(8).eq.'^') THEN
         write (lnfil4) n_negepp,n_resetepp,xspace
         negepp_flag = 1
      endif

      return

  900 FORMAT ('0  *****  LINF4 - VNU LIMITS DO NOT INTERSECT WITH ',    &
     &        'LINFIL - LINF4 NOT USED *****',/,'   VNU = ',F10.3,      &
     &        ' - ',F10.3,' CM-1     LINFIL = ',F10.3,' - ',F10.3,      &
     &        ' CM-1')
  950 FORMAT (8a1)

      END SUBROUTINE

!-----------------------------------------------------------------------
!
!     SUBROUTINE RDLNFL INPUTS THE LINE DATA FROM LINFIL
!
!-----------------------------------------------------------------------
      SUBROUTINE RDLNFL( IEOF,ILO,IHI,BUFR, LINFIL,NLNGTH,VLO,NOPR,IOUT, RDLNBUF )
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_LineData    ,ONLY: MXBRDMOL, &
                                    INPUT_HEADER, INPUT_BLOCK, LINE_DATA
      USE Module_Config      ,ONLY: IPR

      IMPLICIT NONE !REAL*8           (V)

      integer            ,intent(out)   :: IEOF
      integer            ,intent(inout) :: ILO
      integer            ,intent(out)   :: IHI
      TYPE(LINE_DATA)    ,INTENT(OUT)   :: BUFR
      !
      integer            ,intent(in)    :: LINFIL
      integer            ,intent(in)    :: NLNGTH
      real(r8)           ,intent(in)    :: VLO
      integer            ,intent(in)    :: NOPR
      integer            ,intent(out)   :: IOUT(:) !(250)
      TYPE(INPUT_BLOCK)  ,intent(out)   :: RDLNBUF !double convolution needs orignal line data.


      !--- Local variables
      !
      integer*4 ,PARAMETER :: i_1=1

      TYPE(INPUT_HEADER) :: RDLNPNL
      integer*4          :: lnfl,leof,npnlhd
      real*4             :: rdpnl(2),dum(2),xmol(250)

      INTEGER  :: I,         IJ,        ILNGTH
      INTEGER  :: IPASS
      INTEGER  :: J,         J1
      REAL     :: AMOLB(250)
      integer  :: M
      real     :: rvmr




      IPASS = 1
      IF (ILO.GT.0) IPASS = 2

      ILNGTH = NLNGTH*250

      IEOF = 0
      ILO = 1
      IHI = 0

      lnfl = linfil
      npnlhd = 6

!   10 CALL BUFIN_sgl(Lnfl,LEOF,rdpnl(1),npnlhd)
   10 CALL BUFIN_sgl(Lnfl,LEOF,RDLNPNL,npnlhd)


      IF (LEOF.EQ.0) GO TO 30
      IF (RDLNPNL%VMAX.LT.VLO) THEN
         CALL BUFIN_sgl(lnfl,LEOF,dum(1),i_1)
         GO TO 10
      ELSE
!         CALL BUFIN_sgl(Lnfl,LEOF,vlin(1),NWDS)
         CALL BUFIN_sgl(Lnfl,LEOF,rdlnbuf,RDLNPNL%NWDS)
      ENDIF

      IF ((IPASS.EQ.1).AND.(RDLNBUF%VNU(1).GT.VLO)) WRITE (IPR,900)

      IJ = 0

!     precision conversion occurs here:
!     incoming on right: vlin is real*8;  others are real*4 and integer*4

      do 15 i=1,RDLNPNL%NREC
         BUFR%IFLG(i)  = RDLNBUF%IFLG(i)
         BUFR%VNU(i)   = RDLNBUF%VNU(i)
         BUFR%SP(i)    = RDLNBUF%SP(i)
         BUFR%ALFA(i)   = RDLNBUF%ALFA(i)
         BUFR%EPP(i)   = RDLNBUF%EPP(i)

         !yma if (BUFR%IFLG(i) .ge.  0) then
         !yma    BUFR%MOL(i)  = RDLNBUF%MOL(i)   ! int*4 to int*8
         !yma else   ! mol read as integer, treated as real
         !yma     !amolb(i)  = xmol(i)
         !yma    xmol(i)  = transfer (RDLNBUF%mol(i), xmol(i))  !int*4 to real*4
         !yma    amolb(i) = xmol(i)     ! real*4 to real*8 (if compiled as r8)
         !yma    bufr%mol(i) = transfer (amolb(i), bufr%mol(i))  ! real*8 to int*8
         !yma endif
         BUFR%MOL(i)  = RDLNBUF%MOL(i) !The above transfer doesn't work for i4 and r8 setting, so forced BUFR%MOL to be integer*4

         BUFR%HWHM(i) = RDLNBUF%HWHM(i)
         BUFR%TMPALF(i)= RDLNBUF%TMPALF(i)
         BUFR%PSHIFT(i)= RDLNBUF%PSHIFT(i)
         do j=1,mxbrdmol
            bufr%brd_mol_flg(j,i)=rdlnbuf%brd_mol_flg_in(j,i)
         end do
         j = 1
         do j1 = 1,mxbrdmol
            bufr%brd_mol_hw(j1,i) = rdlnbuf%brd_mol_dat(j,i)
            bufr%brd_mol_tmp(j1,i) = rdlnbuf%brd_mol_dat(j+1,i)
            bufr%brd_mol_shft(j1,i) = rdlnbuf%brd_mol_dat(j+2,i)
            j = j+3
         end do
         bufr%speed_dep(i) = rdlnbuf%speed_dep(i)

!  HITRAN provides widths for broadening by air; LBLRTM and MONORTM have always treated these widths as foreign
!  This assumption is valid for most species, but not for N2 or O2. We now adjust the HITRAN widths to obtain
!  true foreign widths. Similar ajdustment is applied if self shift information is available.
         M = MOD(bufr%mol(I),100)
         if (M.eq.7 .AND. bufr%IFLG(i).ge.0) then
            rvmr = 0.21
            bufr%ALFA(i) = ( bufr%ALFA(i)-rvmr*bufr%HWHM(i))/(1.0-rvmr)
            if (bufr%brd_mol_flg(m,i) .gt. 0) then
               bufr%pshift(i) = (bufr%pshift(i)-rvmr*bufr%brd_mol_shft(m,i))/(1.0-rvmr)
            endif
         endif
         if (M.eq.22 .AND. bufr%IFLG(i).ge.0) then
             rvmr = 0.79
             bufr%ALFA(i) = ( bufr%ALFA(i)-rvmr*bufr%HWHM(i))/(1.0-rvmr)
! Currently SBS broadening is only code for the first seven HITRAN species
! When it becomes available for N2, the next three lines should become executable
!            if (bufr%brd_mol_flg(m,i) .gt. 0) then
!               bufr%pshift(i) = (bufr%pshift(i)-rvmr*bufr%brd_mol_shft(m,i))/(1.0-rvmr)
!            endif
         endif


   15 continue

      DO 20 J = 1, RDLNPNL%NREC
         IF (BUFR%IFLG(J).GE.0) THEN
            IJ = IJ+1
            IOUT(IJ) = J
         ENDIF
   20 END DO
      IHI = IJ
      RETURN
   30 IF (NOPR.EQ.0) WRITE (IPR,905)
      IEOF = 1
      RETURN

  900 FORMAT ('0 FIRST LINE USED IN RDLNFL--- CHECK THE LINEFIL  ')
  905 FORMAT ('0 EOF ON LINFIL IN RDLNFL -- CHECK THE LINFIL ')

      END SUBROUTINE

!-----------------------------------------------------------------------
!
!   SUBROUTINE SHRINK COMBINES LINES FALLING IN A WAVENUMBER INTERVAL
!   DVR4/2 INTO A SINGLE EFFECTIVE LINE TO REDUCE COMPUTATION
!
!   TYPE :: LINE_SHRINK
!      REAL(8), DIMENSION(1250)  :: VNU
!      REAL,    DIMENSION(1250)  :: SP, ALFA, EPP
!      INTEGER, DIMENSION(1250)  :: MOL
!      REAL,    DIMENSION(1250)  :: SPP, SRAD
!   END TYPE LINE_SHRINK
!
!      COMMON VNU(1250),S(1250),ALFAL(1250),ALFAD(1250),MOL(1250),       &
!     &       SPP(1250),SRAD(1250)
!
!-----------------------------------------------------------------------
      SUBROUTINE SHRINQ( SHRUNK, ILO,IHI,DVR4 )
!-----------------------------------------------------------------------
      USE Module_ConstParam   ,ONLY: r8=>kind_r8
      USE Module_LineData     ,ONLY: LINE_SHRINK

      IMPLICIT NONE !REAL*8           (V)

      TYPE(LINE_SHRINK) ,intent(inout) :: SHRUNK
      integer           ,intent(inout) :: ILO,IHI
      real              ,intent(in)    :: DVR4

      !--- Local variables
      INTEGER  ::  I,      J
      REAL     ::  DV,     SUMAD, SUMAD2,   SUMAL
      REAL     ::  SUMAL2, SUMC,  SUMC2,    SUMS
      REAL     ::  SUMS2,  SUMV,  SUMV2
      REAL(r8) ::  VLMT


      J = ILO-1
      DV = DVR4/2.
      VLMT = SHRUNK%VNU(ILO)+DV

!     INITIALIZE NON-CO2 SUMS

      SUMAL = 0.
      SUMAD = 0.
      SUMS = 0.
      SUMV = 0.
      SUMC = 0.

!     INITIALIZE CO2 SUMS

      SUMAL2 = 0.
      SUMAD2 = 0.
      SUMS2 = 0.
      SUMV2 = 0.
      SUMC2 = 0.

      DO 20 I = ILO, IHI

!     To prevent underflow issues in CONVF4 we set SP < 1.0e-35 to zero
         IF (SHRUNK%SP(I).lt.1.0e-35) THEN
            SHRUNK%SP(I)= 0.0
            SHRUNK%SPP(I) =0.0
         ENDIF


!     IF LINE COUPLING, DON'T SHRINK LINE

         IF (SHRUNK%SPP(I).NE.0.0) THEN
            J = J+1
            SHRUNK%VNU(J) = SHRUNK%VNU(I)
            SHRUNK%SP(J) = SHRUNK%SP(I)
            SHRUNK%ALFA(J) = SHRUNK%ALFA(I)
            SHRUNK%EPP(J) = SHRUNK%EPP(I)
            SHRUNK%SPP(J) = SHRUNK%SPP(I)
            SHRUNK%MOL(J) = SHRUNK%MOL(I)
            IF (SHRUNK%MOL(J).NE.2) SHRUNK%MOL(J) = 0 !yma, SHRINK()in oprop.f90 doesn't have this

            GO TO 10
         ENDIF

!     NON-CO2 LINES OF MOLECULAR INDEX IT.NE.2   ARE LOADED
!     INTO SUMS IF THE FREQUENCY WITHIN DV GROUP

         IF (SHRUNK%MOL(I).NE.2) THEN
            SUMV = SUMV+SHRUNK%VNU(I)*SHRUNK%SP(I)
            SUMS = SUMS+SHRUNK%SP(I)
            SUMAL = SUMAL+SHRUNK%SP(I)*SHRUNK%ALFA(I)
            SUMAD = SUMAD+SHRUNK%SP(I)*SHRUNK%EPP(I)
            SUMC = SUMC+SHRUNK%SPP(I)
         ELSE

!     CO2 LINES LOADED     (MOL .EQ. 2)

            SUMV2 = SUMV2+SHRUNK%VNU(I)*SHRUNK%SP(I)
            SUMS2 = SUMS2+SHRUNK%SP(I)
            SUMAL2 = SUMAL2+SHRUNK%SP(I)*SHRUNK%ALFA(I)
            SUMAD2 = SUMAD2+SHRUNK%SP(I)*SHRUNK%EPP(I)
            SUMC2 = SUMC2+SHRUNK%SPP(I)
         ENDIF

!     IF LAST LINE OR VNU GREATER THAN LIMIT THEN STORE SUMS

   10    IF (I.LT.IHI) THEN
            IF (SHRUNK%VNU(I+1).LE.VLMT) GO TO 20
         ENDIF

         VLMT = SHRUNK%VNU(I)+DV

!     ASSIGN NON-CO2 LINE AVERAGES TO 'GROUP' LINE J

         IF (SUMS.GT.0.) THEN
            J = J+1
            SHRUNK%SP(J) = SUMS
            SHRUNK%ALFA(J) = SUMAL/SUMS
            SHRUNK%EPP(J) = SUMAD/SUMS
            SHRUNK%VNU(J) = SUMV/SUMS
            SHRUNK%SPP(J) = SUMC
            SHRUNK%MOL(J) = 0
            SUMAL = 0.
            SUMAD = 0.
            SUMS = 0.
            SUMV = 0.
            SUMC = 0.
         ENDIF

!     ASSIGN CO2 LINE AVERAGES

         IF (SUMS2.GT.0.) THEN
            J = J+1
            SHRUNK%SP(J) = SUMS2
            SHRUNK%ALFA(J) = SUMAL2/SUMS2
            SHRUNK%EPP(J) = SUMAD2/SUMS2
            SHRUNK%VNU(J) = SUMV2/SUMS2
            SHRUNK%MOL(J) = 2
            SHRUNK%SPP(J) = SUMC2
            SUMAL2 = 0.
            SUMAD2 = 0.
            SUMS2 = 0.
            SUMV2 = 0.
            SUMC2 = 0.
         ENDIF

   20 END DO

      ILO = J+1
      IHI = J

      RETURN
      END SUBROUTINE

END MODULE
