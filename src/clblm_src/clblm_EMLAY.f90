!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!
MODULE Module_EMLAY
   
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: planck, planckFnDrv, planckFunction
   PUBLIC :: layerEmis, rlin

   
CONTAINS !====================== MODULE CONTAINS =======================   

   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   real ELEMENTAL FUNCTION planck (V, XKT) 
   !--------------------------------------------------------------------
      USE Module_ConstParam, ONLY: r8=>kind_r8, RADCN1, RADCN2
      !RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07; 
      !RADCN2 = planck*clight/boltz
      IMPLICIT NONE
        
      real(r8)      ,intent(in)  :: V
      real          ,intent(in)  :: XKT  !XKT = Tave/RADCN2 = kb*Tave/(h*c) 
      
      real :: cv3,evt
      
      if (XKT<=1.e-6) then
         planck=0.
         RETURN
      endif

      !planck = RADCN1*(V**3) / ( EXP(V/XKT)-1. )       
      cv3 = RADCN1*V*V*V
      evt = EXP(V/XKT)
      planck = cv3 / ( evt - 1. ) 
            
   END FUNCTION

   !--------------------------------------------------------------------
   ! Calculates Planck function and its derivative wrt temperature
   !--------------------------------------------------------------------
   ELEMENTAL SUBROUTINE planckFnDrv (V, T, planck, dBdT) 
   !--------------------------------------------------------------------
      USE Module_ConstParam, ONLY: r8=>kind_r8, RADCN1, RADCN2
      !RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07; 
      !RADCN2 = planck*clight/boltz
      IMPLICIT NONE
        
      real(r8) ,intent(in)  :: V
      real     ,intent(in)  :: T  
      real     ,intent(out) :: planck
      real     ,intent(out) :: dBdT !Derivative w.r.t to temperature
      
      real :: XKT,cv3,evt
      
      XKT = T/RADCN2   !XKT = Tave/RADCN2 = kb*Tave/(h*c)    
      if (XKT<=1.e-6) then
         planck=0.         
         dBdT=0
         RETURN
      endif

      !planck = RADCN1*(V**3) / ( EXP(V/XKT)-1. )       
      cv3 = RADCN1*V*V*V
      evt = EXP(V/XKT)
      planck = cv3 / ( evt - 1. )       
      dBdT = planck*planck * (evt*V)/(RADCN2*cv3*XKT*XKT)
      
   END SUBROUTINE


   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   real ELEMENTAL FUNCTION planckFunction (V,T) 
   !--------------------------------------------------------------------
      USE Module_ConstParam, ONLY: r8=>kind_r8, RADCN1, RADCN2 
      !RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07; 
      !RADCN2 = planck*clight/boltz
      IMPLICIT NONE
        
      real(r8)  ,intent(in) :: V
      real      ,intent(in) :: T
      
      if (T<=0.) then
         planckFunction=0.
         RETURN
      endif
      
      !planckFunction = RADCN1*(V**3) / ( EXP(V*RADCN2/T)-1. )
      planckFunction = RADCN1*V*V*V / ( EXP(V*RADCN2/T)-1. )
      
   END FUNCTION
   
 

   !--------------------------------------------------------------------
   ! layerEmis calculate thermal emission from an atmospheric layer. It always calculates 
   ! emission from "Ta" side. For example, if Ta=Ttop and Tb=Tbot, the subroutine calculates
   ! upwelling emission. To calculate downwelling emission, set Ta=Tbot and Tb=Ttop.
   !--------------------------------------------------------------------
   SUBROUTINE layerEmis( lyrEm, Ta, Tb, Tave, lyrTx, lyrOD, linInTau, &
                         NLTE, nlteParam )
   !--------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8, RADCN2
      USE Module_Config     ,ONLY: IPR
      USE Module_Spectrum   ,ONLY: CLBLM_Spectrum, &
                                   CLBLM_Spectrum_init
      !USE Module_PlanckFunc ,ONLY: planck
      IMPLICIT NONE
      
      type(CLBLM_spectrum)            ,intent(out) :: lyrEm
      real                            ,intent(in)  :: Ta, Tb, Tave
      type(CLBLM_spectrum)            ,intent(in)  :: lyrTx
      type(CLBLM_spectrum)            ,intent(in)  :: lyrOD
      integer                         ,intent(in)  :: linInTau
      logical               ,optional ,intent(in)  :: NLTE
      type(CLBLM_spectrum)  ,optional ,intent(in)  :: nlteParam
      
      
      integer           :: ip, ind1, ind2
      logical           :: nlteFlag
      real              :: ODVI, abso, hlf, c_nlte
      real              :: XKTA, XKTB, XKT, bba, bbb, bb
      real ,allocatable :: em(:)
      real              :: DV
      real(r8)          :: vi, V1
      
      
      nlteFlag = .false.
      if (present(NLTE)) nlteFlag=NLTE

      
      XKTA = Ta   / RADCN2
      XKTB = Tb   / RADCN2
      XKT  = Tave / RADCN2
      V1   = lyrOD%V1
      DV   = lyrOD%DV
      ind1 = lyrOD%indV1
      ind2 = lyrOD%indV2
      allocate(em( ind1:ind2 )) !;EMLYRA(:) = 0. 

      do ip = ind1, ind2

         odvi   = lyrOD%spect(ip)
         abso   = 1-lyrTx%spect(ip)
         if (nlteFlag) then
            c_nlte = nlteParam%spect(ip)
            abso   = (1. - c_nlte/odvi) * abso
         endif
         
         if (linInTau/=0) call rlin( odvi, hlf )          
         vi = V1 + REAL(ip-1)*DV
         
         select case(linInTau)
         case (0) !not linear-in-tau
            bb     = planck(vi,XKT )
            em(ip) = abso * bb
         case (1) !LBLRTM formulation
            bba    = planck(vi,XKTA)
            bb     = planck(vi,XKT )
            em(ip) = abso * (bba + 2*(bb-bba)*hlf)
         case (2) !standard formulation of linear-in-tau
            bba    = planck(vi,XKTA)
            bbb    = planck(vi,XKTB)
            em(ip) = abso * (bba + (bbb-bba)*hlf)
            !em(ip) = abso*( bba*(1.-hlf) + bbb*hlf )            
         case default
            STOP '--- layerEmis(): Invalid value for linInTau.'
         end select
      enddo
      
      call CLBLM_Spectrum_init( lyrEm,  lyrOD%V1, lyrOD%DV, lyrOD%NLIM )                  
      call move_alloc( em, lyrEm%spect )

   END SUBROUTINE
   


  !----------------------------------------------------------------------------
  ! For linear-In-Tau approach, compute the Planck radiance weight 
  ! ( hlf=(1/od - T/(1-T)) ) of the far boundary of the layer, 
  ! and its derivative (dhlf)
  !----------------------------------------------------------------------------
  SUBROUTINE rlin( od, hlf,dhlf )
  !----------------------------------------------------------------------------
      IMPLICIT NONE    
      real          ,intent(in)  :: od
      real          ,intent(out) :: hlf
      real ,optional,intent(out) :: dhlf
      
      REAL, PARAMETER :: dbAvg2dtau=-1./12!Weighting factor derivative wrt optical depth, for linear-in-tau, in low-OD limit
      real :: tex,trans
      
      
      if (od.lt.1.e-03) then      
         hlf  = .5+od*dbAvg2dtau
         !dhlf = dbAvg2dtau         
      else if (od.gt.20.) then      
         hlf  = 1./od
         !dhlf = -1./od**2         
      else      
         !hlf= 1. - 2.*(trans/(trans-1.) + 1./od)
         !hlf=1./od-trans/od
         tex  = exp(od)
         hlf  = 1./(od)-1./(tex-1.)
         !dhlf = -1./od**2+tex/(tex-1)**2
      endif

      
      if (present(dhlf)) then
         if (od.lt.1.e-03) then      
            dhlf = dbAvg2dtau         
         else if (od.gt.20.) then      
            dhlf = -1./od**2         
         else      
            dhlf = -1./od**2+tex/(tex-1)**2
         endif
      endif
      
   END SUBROUTINE
      
END MODULE
