!=======================================================================
!
!                             SELF

!     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
!
   if ((V2.gt.-20.0).and.(V1.lt.20000.) .and. xself.gt.0.) then
      sh2ot0 = 0.
      sh2ot1 = 0.
!
      CALL SL296 (V1C,V2C,DVC,NPTC,SH2OT0,v1ss,v2ss)
      CALL SL260 (V1C,V2C,DVC,NPTC,SH2OT1,v1ss,v2ss)
!
!           Loop calculating self continuum optical depth
!
      TFAC = (TAVE-T0)/(260.-T0)

!  MT_CKD_3.5  All previous IR corrections now included in stored coefficients
!  rather than correction functions.

      DO 20 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         SH2O = 0.
         IF (SH2OT0(J).GT.0.) THEN
            SH2O = SH2OT0(J)*(SH2OT1(J)/SH2OT0(J))**TFAC
         ENDIF
!              ---------------------------------------------------------
!
         cself(j) = WK(1)*(SH2O*Rself)
!
!********************************************
         v1h=V1C
         dvh=DVC
         npth=NPTC
!
         csh2o(j)=1.e-20 * sh2o * xself
!********************************************
!
!              ---------------------------------------------------------
!              Radiation field
!
         IF (JRAD.EQ.1) cself(j) = cself(j)*RADFN(VJ,XKT)
!              ---------------------------------------------------------

20    CONTINUE
!
!           Interpolate to total optical depth grid

      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,cself,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif
!
!=======================================================================
