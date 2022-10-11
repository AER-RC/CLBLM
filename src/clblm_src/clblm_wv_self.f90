!=======================================================================
!                                                                       
!                             SELF                                      
                                                                        
!     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1            
!                                                                       
      if ((V2.gt.-20.0).and.(V1.lt.20000.) .and. xself.gt.0.) then      

         sh2ot0 = 0.
         sh2ot1 = 0.

            CALL SL296 (V1C,V2C,DVC,NPTC,SH2OT0, V1ABS,V2ABS )
            CALL SL260 (V1C,V2C,DVC,NPTC,SH2OT1, V1ABS,V2ABS ) 
!                                                                       
!           Loop calculating self continuum optical depth               
!                                                                       
            TFAC = (TAVE-T0)/(260.-T0)                                  
!                                                                       
!-----------------------------------------------------------------------
!          CORRECTION TO SELF CONTINUUM   mt_ckd_2.4  Nov 2008    sac   
!-----------------------------------------------------------------------
!                                                                       
            f1    = 0.25                                                
            beta  = 350.                                                
! ***
!   Correction from RHUBC-II    mt_ckd_3.0   Nov 2016
! ***
            f1_rhu    = 0.08                                                
            beta_rhu  = 40.                                                
            n_s   = 6 
            
            DO 20 J = 1, NPTC                                           
               VJ = V1C+DVC* REAL(J-1)                                  
               SH2O = 0.                                                
               IF (SH2OT0(J).GT.0.) THEN                                
                  SH2O = SH2OT0(J)*(SH2OT1(J)/SH2OT0(J))**TFAC          
                  SFAC = 1.                                             
!                                                                       
                  IF (VJ .GE. 820. .AND. VJ .LE. 960.) THEN             
                     JFAC = (VJ - 820.)/10. + 0.00001                   
                     SFAC = XFACREV(JFAC)                               
                  ENDIF                                                 
!v128 !   ***                                                                 
!v128 !     Correction to the self continuum     mt_ckd_2.5   Jan 2010        
!v128 !   ***                                                                 
                  sfac = sfac * ( 1 + ( f1/(1+(VJ/beta)**n_s) ) )       
                  sfac = sfac * ( 1 + ( f1_rhu/(1+(VJ/beta_rhu)**n_s)))       
                                                                        
                  SH2O = SFAC * SH2O                                    
!                                                                       
               ENDIF                                                    
!              ---------------------------------------------------------
!                                                                       
               cself(j) = WK(indH2O)*(SH2O*Rself) !WK(1)*(SH2O*Rself)                            
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
                                                                        
   20       CONTINUE                                                    
!                                                                       
!           Interpolate to total optical depth grid                     
                                                                        
            CALL XINT (V1C,V2C,DVC,cself,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)
                                                                        
         endif                                                          
!                                                                       
!=======================================================================
