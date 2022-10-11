!
! CREATION HISTORY:
!       Modified from LBLRTM v12.9
!       Yingtao Ma, AER@NOAA/NESDIS/STAR
!       yma@aer.com; yingtao.ma@noaa.gov
!

Module Module_Solar

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: solarRadiance



CONTAINS !===================== MODULE CONTAIS =========================


!-----------------------------------------------------------------------
!     SUBROUTINE solarIrradiance read solar irradiance data from the binary
!     file SOLAR.RAD.  The following are input and output options:
!
!       ISOLVAR = -1 => uses solar source file with a single component:
!                      no temporal variability; assumes Kurucz extraterrestrial
!                      solar irradiance, which yields a solar constant of 1368.22 Wm-2,
!                      unless scaled by SCON.
!
!       ISOLVAR =  0 => uses solar source file with a single component:
!                      no temporal variability; assumes NRLSSI2 extraterrestrial
!                      solar irradiance, which yields a solar constant of 1360.85 Wm-2,
!                      (for the spectral range 100-50000 cm-1 with quiet sun, facular
!                      and sunspot contributions fixed to the mean of
!                      Solar Cycles 13-24 and averaged over the mean solar cycle),
!                      unless scaled by SCON.

!       SCON   = 0.0 => no scaling of internal solar irradiance
!              > 0.0 => ISOLVAR = -1 or 0
!                       Total solar irradiance is scaled to SCON
!              > 0.0 => ISOLVAR = 1
!                       integral of total solar irradiance averaged over solar cycle
!                       is scaled to SCON
!
!       ISOLVAR =  1 => uses solar source file with multiple components;
!                      facular brightening and sunspot blocking amplitudes
!                      are by default determined by the fraction of the
!                      way into the solar cycle (see SOLCYCFRAC) or can
!                      be scaled (see SOLVAR)
!       ISOLVAR =  2 => uses solar source file with multiple components;
!                      facular brightening and sunspot blocking amplitudes
!                      are determined by the Mg and SB indeces (see SOLVAR)
!
!       SOLCYCFRAC     Solar cycle fraction (0-1); fraction of the way through the mean 11-year
!                      cycle with 0 and 1 defined as the minimum phase of the solar cycle
!                      (ISOLVAR=1 only)

!       SOLVAR         Solar variability scaling factors or indices (ISOLVAR=1,2 only)
!                      ISOLVAR = 1 =>
!                         SOLVAR(1)    Facular (Mg) index amplitude scale factor
!                         SOLVAR(2)    Sunspot (SB) index amplitude scale factor

!                         ISOLVAR = 2 =>
!                         SOLVAR(1)    Facular (Mg) index as defined in the NRLSSI2 model;
!                                      used for modeling specific solar activity
!                         SOLVAR(2)    Sunspot (SB) index as defined in the NRLSSI2 model;
!                                      used for modeling specific solar activity
!
!     Output radiance goes to TAPE13.
!
!
!-----------------------------------------------------------------------
   SUBROUTINE solarRadiance(solRad, V1O,V2O, JulDay, solarSource, solarVarOption, &
                                    SCON, SOLCYCFRAC, faculaVar, spotVar)
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: pi, r8=>kind_r8
!      USE Module_Utility     ,ONLY: getlun
      USE Module_Config      ,ONLY: IPR,ioFiles
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum, CLBLM_Spectrum_init, XINT
      IMPLICIT NONE

      type(CLBLM_Spectrum) ,intent(out)    :: solRad  !Solar radiance at TOA in (W/cm2 sr cm-1)
      real(r8)             ,intent(in)     :: V1O,V2O
      integer              ,intent(in)     :: JulDay
      integer              ,intent(in)     :: solarSource
      integer              ,intent(in)     :: solarVarOption
      real                 ,intent(in)     :: SCON
      real                 ,intent(in)     :: SOLCYCFRAC
      real                 ,intent(in)     :: faculaVar
      real                 ,intent(in)     :: spotVar


      !integer ,PARAMETER :: NSOL=2000001
      integer           :: ISOLVAR
      integer           :: i,j,k, np
      real              :: solvar(2)
      character(256)    :: solarFile
!      integer           :: ISOLFL, LSEOF
      real              :: svar(3)
      real              :: theta, XJUL_SCALE
      real              :: conv_ster, SCAL_FAC
!      real              :: SOLAR(2400)  !SOLAR(-1:NSOL)
!      real              :: SOLAR3(2400,3)
!      real*8            :: V1P,V2P
!      real              :: DVP
!      integer           :: NLIM
      real ,allocatable :: tempSolIrr(:)
!      real              :: tempV1O,tempV2O, tempV
!      integer           :: numPts,np

      real(r8)          :: v1_sol, v2_sol
      real              :: dv_sol
      integer           :: nPts_sol    !number of spectral points readed
      integer           :: nComp_sol   !number solar irradiance components
      real ,allocatable :: solData(:,:)

!      real*4 :: xhdr_s(265)
!      real*4 :: dv_sol
!      real*8 :: v1_sol
!      real*8 :: v2_sol
!      equivalence (dv_sol,xhdr_s(218)),&
!                  (v1_sol,xhdr_s(219)),&
!                  (v2_sol,xhdr_s(221))


      solvar(:) = [faculaVar,spotVar]

      if     (solarSource==1) then;                         ISOLVAR = -1 !Use Kurucz solar irradiance
      elseif (solarSource==2 .and. solarVarOption==1) then; ISOLVAR =  0 !Use NRLSSI2 solar data without temporal variability
      elseif (solarSource==2 .and. solarVarOption==2) then; ISOLVAR =  1 !Use NRLSSI2 solar data, FB and SB terms determined by the solar cyble fraction or can be scaled by SOLVAR(:)
      elseif (solarSource==2 .and. solarVarOption==3) then; ISOLVAR =  2 !Use NRLSSI2 solar data, FB and SB terms determined by determined by the Mg and SB indeces carried in SOLVAR(:)
      endif

      if     (ISOLVAR==-1) then; solarFile = trim(ioFiles%solarPath)//trim(ioFiles%solarFile_Kurucz)
      elseif (ISOLVAR== 0) then; solarFile = trim(ioFiles%solarPath)//trim(ioFiles%solarFile_1comp_NRL)
      elseif (ISOLVAR>= 1) then; solarFile = trim(ioFiles%solarPath)//trim(ioFiles%solarFile_3comp_NRL)
      endif


!      !--- Open file SOLAR.RAD
!      ! Note that the file SOLAR.RAD is always single precision. Provision
!      ! has been made to deal with the case in which the current program i
!      ! double precision.
!      !
!      ISOLFL = getlun()
!      OPEN(UNIT=ISOLFL,FILE=trim(ioFiles%solarFile),FORM='UNFORMATTED', STATUS='OLD')


      !--- Calculate solar scaling factors
      call scale_solar (isolvar,scon,solcycfrac,solvar,svar)


      !--- Calculate Earth distance to sun given Julian Day JulDay. Used to
      ! scale solar source function. Formula taken from "Atmospheric Radia
      ! Transfer", J. Lenoble, 1993.
      !
      ! Test validity of JulDay
      if ((JulDay.lt.0).or.(JulDay.gt.366)) then
         write(*,*) 'JulDay = ',JulDay,' is out of range 0-366.'
         write(ipr,*) 'JulDay = ',JulDay,' is out of range 0-366.'
         stop 'Stopped in SOLINT'
      endif

      !If JulDay = 0 , then set XJUL_SCALE to 1
      if (JulDay.eq.0) then
         XJUL_SCALE = 1.0
         write(ipr,*) 'JulDay = 0, no scaling of solar source function'
      else
         theta = 2*pi*( REAL(JulDay)-1.)/365.
         XJUL_SCALE = 1.00011 + 0.034221*cos(theta) +  &
                      1.28E-3*sin(theta) + 7.19E-4*cos(2.*theta) + &
                      7.7E-5*sin(2.*theta)
         write(ipr,*) 'JulDay = ',JulDay, ', scale factor for solar source function = ',XJUL_SCALE
      endif


      !--- Combine Julian day scaling with solar source scaling
      ! when there is no solar variability
      if (isolvar.le.0) then
          xjul_scale = xjul_scale*solvar(1)
      endif


!      !--- Read file header of solar radiance file and determine dv ratio
!      read (isolfl) (xhdr_s(i),i=1,264)
!      !dv = dv_sol
!      !v1 = v1_sol
!      !v2 = v2_sol
!
!
!      allocate( tempSolIrr( ceiling( (V2O-V1O)/dv_sol +1. ) + 10 ) ) !get 10 more points for safety
!      numPts = 0
!      Do !till all solar data in between V1O to V2O are read
!
!         !--- Read in the solar data
!         if (isolvar.le.0) then
!            CALL SOLIN_sgl (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR,LSEOF)
!         else
!            CALL SOLIN_tri (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR3,LSEOF)
!            call comb_solar (solar3,solar,nlim,isolvar,svar)
!         end if
!
!         if ( V2P < V1O ) CYCLE !not yet reach the requested V1
!
!
!         !--- If it is the last panel needed, reset NLIM and V2P
!         if ( V2P > V2O ) then
!            !--- Find the smallest V greater then V2O
!            do i =1,NLIM
!               tempV = V1P + real(i-1)*DVP
!               if ( tempV > V2O ) exit
!            enddo
!
!            !--- Cut the panel
!            NLIM = i
!            V2P  = V1P + (NLIM-1)*DVP
!            LSEOF = 0 !set flag to exit loop
!         endif
!
!
!         if (numpts==0) then !If first panel
!
!            !--- Find the largest V less then V1O
!            do i =1,NLIM
!               tempV = V1P + real(i-1)*DVP
!               if ( tempV > V1O ) exit
!            enddo
!
!            np = NLIM-(i-1)+1
!            tempSolIrr( 1:np ) = solar(i-1:NLIM)
!            tempV1O = tempV-DVP
!            numPts = numPts+np
!         endif
!
!         tempSolIrr(numPts+1:numPts+NLIM) = solar(1:NLIM)
!         tempV2O = V2P
!         numPts = numPts + NLIM
!
!         if (LSEOF==0) EXIT
!      enddo


      !--- Read in the solar data
      CALL readSolarData( trim(solarFile), V1O,V2O, &
                          v1_sol, v2_sol, dv_sol, nPts_sol, nComp_sol, solData)

      allocate( tempSolIrr( nPts_sol ) )
      if (isolvar.le.0) then
         tempSolIrr(:) = solData(:,1)
      else
         call comb_solar( solData, tempSolIrr, nPts_sol, isolvar, svar)
      end if

      !--- Solar irradiance is input from SOLAR.RAD (W/m2 cm-1).
      ! Convert to radiance units (W/cm2 sr cm-1) by multiplying
      ! by 1/6.8e-5.
      conv_ster = 1./(1.e4*6.8e-5)

      !--- Combine XJUL_SCALE and conv_ster into one scale factor SCAL_FAC
      SCAL_FAC = conv_ster*XJUL_SCALE

      !--- Scale the solar irradiance
      DO i = 1,nPts_sol
         tempSolIrr(i) = tempSolIrr(i)*SCAL_FAC
      ENDDO

      !--- Fill the output structure
      np = ceiling((V2O-V1O)/dv_sol +1.)
      CALL CLBLM_Spectrum_init( solRad, V1O, real(dv_sol), np )
      CALL XINT( v1_sol, v2_sol, real(dv_sol), tempSolIrr, 1.,&
                 V1O, real(dv_sol), solRad%spect, 1, np)

   END SUBROUTINE


!-----------------------------------------------------------------------
! write solar irradiance data to netcdf file
!-----------------------------------------------------------------------
   subroutine writeSolarData( filename, V1,V2,DV,nPts,nComp,solData, info )
!-----------------------------------------------------------------------
    USE Module_ConstParam ,ONLY: r8=>kind_r8
    USE NETCDF
    USE Module_FileIO  ,ONLY: checkNetCDFCall, createNetCDFFile

    character(*)           ,intent(in) :: filename
    real(r8)               ,intent(in) :: V1,V2
    real                   ,intent(in) :: DV
    integer                ,intent(in) :: nPts   !number of spectral points
    integer                ,intent(in) :: nComp  !number solar irradiance components
    real                   ,intent(in) :: solData(:,:)
    character(*) ,optional ,intent(in) :: info


    integer(4) :: ncid_out, v1ID,v2ID,dvID, soldatID, nPtsID, nCompID


    if (.not. createNetCDFFile(trim(filename), ncid_out)) then
        stop "Cannot create netcdf file for solar irradiance, this should not happen"
    endif

    if (present(info)) then
      call checkNetCDFCall( nf90_put_att(ncid_out, NF90_GLOBAL, "SolarDataInfor", trim(info)))
    endif

    call checkNetCDFCall(nf90_def_dim(ncid_out, "numPoints",     int(size(solData,1),4), nPtsID))
    call checkNetCDFCall(nf90_def_dim(ncid_out, "numComponents", int(size(solData,2),4), nCompID))

    call checknetcdfcall(nf90_def_var(ncid_out, "StartWavenumber",   NF90_DOUBLE, v1ID))
    call checknetcdfcall(nf90_def_var(ncid_out, "EndWavenumber",     NF90_DOUBLE, v2ID))
    call checknetcdfcall(nf90_def_var(ncid_out, "SpectralIncrement", NF90_DOUBLE,  dvID))
    call checknetcdfcall(nf90_def_var(ncid_out, "SolarData",         NF90_DOUBLE, [nPtsID, nCompID], soldatID))

    call checknetcdfcall(nf90_enddef(ncid_out))

    call checknetcdfcall(nf90_put_var(ncid_out, v1ID, V1))
    call checknetcdfcall(nf90_put_var(ncid_out, v2ID, V2))
    call checknetcdfcall(nf90_put_var(ncid_out, dvID, DV))
    call checknetcdfcall(nf90_put_var(ncid_out, soldatID, solData))

    call checknetcdfcall(nf90_close(ncid_out))
   END SUBROUTINE


!-----------------------------------------------------------------------
! Read solar irradiance data in user provided spectral range [V1,V2]. The output
! data range is determined as  i1 = floor( (V1-V1all)/DVall + 1. ) and
! i2 = ceiling( (V2-V1all)/DVall + 1. )
!-----------------------------------------------------------------------
   SUBROUTINE readSolarData( filename, V1,V2, V1out,V2out,DVout,nPtsOut,nComp,solData, info)
!-----------------------------------------------------------------------
      USE Module_ConstParam ,ONLY: r8=>kind_r8
      USE NETCDF
      USE Module_FileIO, only : checkNetCDFCall, createNetCDFFile

      character(*)           ,intent(in)    :: filename
      real(r8)               ,intent(in)    :: V1,V2
      real(r8)               ,intent(out)   :: V1out,V2out
      real                   ,intent(out)   :: DVout
      integer                ,intent(out)   :: nPtsOut  !number of spectral points readed
      integer                ,intent(out)   :: nComp    !number solar irradiance components
      real      ,allocatable ,intent(out)   :: solData(:,:)
      character(*) ,optional ,intent(out)   :: info

      character(*) ,parameter :: routineName = 'readSolarData'

      integer        :: i1,i2
      integer(4)     :: nPtsAll_int4, nComp_int4
      real(r8)       :: V1all, V2all, DVall
      integer(4)     :: ncid, v1ID,v2ID,dvID, soldatID, nPtsID, nCompID
      character(256) :: infoStr



      call checkNetCDFCall(nf90_open(trim(filename), NF90_NOWRITE, ncid))
      !call checkNEtCDFcall(nf90_inquire(ncid, nvars, ndims, nglobalatts, nunlimdimid))

      call checkNetCDFCall(nf90_get_Att( ncid, NF90_GLOBAL, "SolarDataInfor", infoStr))
      if (present(info)) info = trim(infoStr)

      !--- Get dimensions
      call checkNEtCDFcall(nf90_inq_dimid(ncid, 'numPoints',     nPtsID))
      call checkNetcdfcall(nf90_inq_dimid(ncid, 'numComponents', nCompID))
      call checknetcdfcall(nf90_inquire_dimension(ncid, nPtsID,  len=nPtsAll_int4))
      call checknetcdfcall(nf90_inquire_dimension(ncid, nCompID, len=nComp_int4))
      nComp = nComp_int4

      !--- Get V1all,V2all,DVall of whole data set
      !
      call checknetcdfcall(nf90_inq_varid(ncid, 'StartWavenumber', v1ID))
      call checknetcdfcall(nf90_get_var(ncid, v1ID, V1all))

      call checknetcdfcall(nf90_inq_varid(ncid, 'EndWavenumber', v2ID))
      call checknetcdfcall(nf90_get_var(ncid, v2ID, V2all))

      call checknetcdfcall(nf90_inq_varid(ncid, 'SpectralIncrement', dvID))
      call checknetcdfcall(nf90_get_var(ncid, dvID, DVall))


      !--- Read solar data in between requested V1 and V2
      !
      print *,'V1', V1
      print *,'V1all', V1
      if ( V1<V1all .or. V1>V2all ) STOP '---'//routineName//'(): Requested V1 is out of data range.'
      if ( V2<V1all .or. V2>V2all ) STOP '---'//routineName//'(): Requested V2 is out of data range.'

      i1 = floor( (V1-V1all)/DVall + 1. )
      i2 = ceiling( (V2-V1all)/DVall + 1. )
      if (i1<1) i1=1
      if (i2>nPtsAll_int4) i2=nPtsAll_int4
      nPtsOut = i2-i1+1

      if (allocated(solData)) deallocate(solData)
      allocate( solData( nPtsOut, nComp ))

      call checknetcdfcall(nf90_inq_varid(ncid, 'SolarData', soldatID))
      call checknetcdfcall(nf90_get_var( ncid, soldatID, solData, &
                                         start=int([i1,1],4), count=int([nPtsOut,nComp],4) ))

      call checknetcdfcall(nf90_close(ncid))

      V1out = V1all + (i1-1)*DVall
      V2out = V1all + (i2-1)*DVall
      DVout = DVall

   END SUBROUTINE

!-----------------------------------------------------------------------
!                         Written:    January 2017
!
!                  IMPLEMENTATION:    K. E. Cady-Pereira
!
! Calculates the scaling factors  of the solar term from the options
! selected by the user
!-----------------------------------------------------------------------
   SUBROUTINE scale_solar (isolvar,scon,solcycfrac,solvar,svar)
!-----------------------------------------------------------------------
      use solar_cycle

      integer ,intent(in)    :: isolvar
      real    ,intent(in)    :: scon
      real    ,intent(in)    :: solcycfrac
      real    ,intent(inout) :: solvar(:)
      real    ,intent(out)   :: svar(:)

      !dimension solvar(2)   ! input solar variability scale factors or indices
      !real svar(3)          ! final quiet, facular and sunspot term scale factors
      real    :: indsolvar(2)     ! calculated solar index amplitude scale factors or indices
      real    :: nsfm1_inv, nsfm1_inv_2
      real    :: wgt,tmp_f_0,tmp_s_0
      integer :: isfid
      real    :: fraclo, frachi
      real    :: intfrac
      real    :: svar_f_0, svar_s_0


      indsolvar(:) = 1.0
      svar(:) = 1.0

      if (isolvar.eq.1) then
         indsolvar=solvar

! Adjust amplitude scaling to be 1.0 at solar min to be the requested indsolvar at solar max,
!  and to vary between those values at other solcycfrac.
         if (indsolvar(1).ne.1.0.or.indsolvar(2).ne.1.0) then
            if (solcycfrac.ge.0.0.and.solcycfrac.lt.solcyc_min) then
               wgt = (solcycfrac+1.0-solcyc_max)/(1.0+solcyc_min-solcyc_max)
               indsolvar(:) = indsolvar(:) + wgt * (1.0-indsolvar(:))
            endif
            if (solcycfrac.ge.solcyc_min.and.solcycfrac.le.solcyc_max) then
               wgt = (solcycfrac-solcyc_min)/(solcyc_max-solcyc_min)
               indsolvar(:) = 1.0 + wgt * (indsolvar(:)-1.0)
            endif
            if (solcycfrac.gt.solcyc_max.and.solcycfrac.le.1.0) then
               wgt = (solcycfrac-solcyc_max)/(1.0+solcyc_min-solcyc_max)
               indsolvar(:) = indsolvar(:) + wgt * (1.0-indsolvar(:))
            endif
         endif

!   Interpolate svar_f_0 and svar_s_0 from lookup tables using provided solar cycle fraction
!   Lookup tables points are located at the middle of each month, except for first and last points,
!   which correspond to the first and last days of the 11-year solar cycle.
	      if (solcycfrac .le. 0.0) then
	         tmp_f_0 = mgavgcyc(1)
	         tmp_s_0 = sbavgcyc(1)
	      elseif (solcycfrac .ge. 1.0) then
	         tmp_f_0 = mgavgcyc(nsolfrac)
	         tmp_s_0 = sbavgcyc(nsolfrac)
	      else
	         nsfm1_inv = 1.0 / (nsolfrac-2)
            nsfm1_inv_2 = nsfm1_inv/2.0
            if (solcycfrac.le.nsfm1_inv_2) then
               isfid = 0
               fraclo = 0.0
               frachi = nsfm1_inv_2
            elseif (solcycfrac.gt.(1.0-nsfm1_inv_2)) then
               isfid = nsolfrac-1
               fraclo = 1.0-nsfm1_inv_2
               frachi = 1.0
            else
	            isfid = floor((solcycfrac-nsfm1_inv_2) * (nsolfrac-3)) + 2
	            fraclo = (isfid-2) * nsfm1_inv+nsfm1_inv_2
	            frachi = (isfid-1) * nsfm1_inv+nsfm1_inv_2
            endif
	         intfrac = (solcycfrac - fraclo) / (frachi - fraclo)
	         tmp_f_0 = mgavgcyc(isfid) + intfrac * (mgavgcyc(isfid+1) - mgavgcyc(isfid))
	         tmp_s_0 = sbavgcyc(isfid) + intfrac * (sbavgcyc(isfid+1) - sbavgcyc(isfid))
	      endif
	      svar_f_0 = tmp_f_0
	      svar_s_0 = tmp_s_0
      endif

      if (isolvar.eq.2) indsolvar(:) = solvar

      if (scon <= 0.0) then
!   No solar cycle and no solar variability (Kurucz solar irradiance)
         if (isolvar .eq. -1) then
            solvar(:) = 1.0
         endif

!   Mean solar cycle with no solar variability (NRLSSI2 model solar irradiance)
!   Quiet sun, facular, and sunspot terms averaged over the mean solar cycle
!   (defined as average of Solar Cycles 13-24).
         if (isolvar .eq. 0) then
            solvar(:) = 1.0
         endif

!   Mean solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot terms interpolated from LUTs to input solar cycle
!   fraction for mean solar cycle. Scalings defined below to convert from
!   averaged Mg and SB terms to Mg and SB terms interpolated here.
!   (Includes optional facular and sunspot amplitude scale factors)
         if (isolvar .eq. 1) then
            svar(2) = indsolvar(1) * (svar_f_0 - Foffset) / (svar_f_avg - Foffset)
            svar(3) = indsolvar(2) * (svar_s_0 - Soffset) / (svar_s_avg - Soffset)
         endif

!   Specific solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot index terms input directly to model specific
!   solar cycle.  Scalings defined below to convert from averaged
!   Mg and SB terms to specified Mg and SB terms.
         if (isolvar .eq. 2) then
            svar(2)= (indsolvar(1) - Foffset) / (svar_f_avg - Foffset)
            svar(3) = (indsolvar(2) - Soffset) / (svar_s_avg - Soffset)
         endif

     endif

     if (scon.gt.0) then

!   No solar cycle and no solar variability (Kurucz or NRLSSI2 model solar irradiance)
!   Quiet sun, facular, and sunspot terms averaged over the mean solar cycle
!   (defined as average of Solar Cycles 13-24).
!   Scale from internal solar constant to requested solar constant.
         if (isolvar .eq. -1) then
            solvar(:) = scon/scon_kurucz
         else if (isolvar.eq.0) then
            solvar(:) = scon/scon_nrlssi2
         end if

!   Mean solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot terms interpolated from LUTs to input solar cycle
!   fraction for mean solar cycle. Scalings defined below to convert from
!   averaged Mg and SB terms to Mg and SB terms interpolated here.
!   Scale internal solar constant to requested solar constant.
!   (Includes optional facular and sunspot amplitude scale factors)

         if (isolvar .eq. 1) then
            svar(1)= (scon - (indsolvar(1) * Fint + indsolvar(2) * Sint)) / Iint
            svar(2)= indsolvar(1) * (svar_f_0 - Foffset) / (svar_f_avg - Foffset)
            svar(3)= indsolvar(2) * (svar_s_0 - Soffset) / (svar_s_avg - Soffset)
         endif

      endif


      return
   END SUBROUTINE

!-----------------------------------------------------------------------
!
!                         Written:    January 2017
!
!                  IMPLEMENTATION:    K. E. Cady-Pereira
!
! Combines the three components (quiet, facular and sunspot contributions) according
! to the options selected by the user
!
!-----------------------------------------------------------------------
   subroutine comb_solar (solar3,solar,nlim,isolvar,svar)
!-----------------------------------------------------------------------
      !use solar_cycle

      real    ,intent(in)    :: solar3(:,:)
      real    ,intent(inout) :: solar(:)
      integer ,intent(in)    :: nlim
      integer ,intent(in)    :: isolvar
      real    ,intent(in)    :: svar(:)

      !real :: SOLAR3(2400,3)
      !real :: SOLAR(*)
      !integer :: nlim, isolvar
      !real :: svar(3)

      !if ((isolvar.eq.0.OR.isolvar.eq.-1).AND.scon.eq.0.0) then
      !if (isolvar.eq.0) then
      !   solar(1:nlim) = solar3(1:nlim,1)+solar3(1:nlim,2)+solar3(1:nlim,3)
      !endif

      if (isolvar.gt.0.AND.isolvar.le.2) then
         solar(1:nlim) = solar3(1:nlim,1)*svar(1) &
                        +solar3(1:nlim,2)*svar(2) &
                        +solar3(1:nlim,3)*svar(3)
      end if

      return
   END SUBROUTINE























!-----------------------------------------------------------------------
!               LAST MODIFICATION:    1 November 1995
!
!                  IMPLEMENTATION:    P.D. Brown
!
!             ALGORITHM REVISIONS:    S.A. Clough
!                                     P.D. Brown
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!     SUBROUTINE SOLIN_sgl inputs files for use with solar radiation
!     calculations for interpolation in SOLINT.  Reads files with
!     one record per panel. SOLIN_sgl reads single precision files
!-----------------------------------------------------------------------
   SUBROUTINE SOLIN_sgl (V1P,V2P,DVP,NLIM,KFILE,XARRAY,KEOF)
!-----------------------------------------------------------------------
      real*8  ,intent(out)   :: V1P,V2P
      real    ,intent(out)   :: DVP
      integer ,intent(out)   :: NLIM
      integer ,intent(in)    :: KFILE
      real    ,intent(out)   :: XARRAY(:)
      integer ,intent(inout) :: KEOF

      integer   :: i
      real*4    :: xarray_s(2410)
      integer*4 :: kfil_s, keof_s, nphdr_s, nphdrf
      data nphdrf / 6 /

      real*8    :: V1PBF,V2PBF
      real*4    :: dvpbf
      integer*4 :: nlimbf
      real*4    :: PNLHDR(2)
      COMMON /BUFPNL_s/ V1PBF,V2PBF,dvpbf,nlimbf
      EQUIVALENCE (PNLHDR(1),V1PBF)

      !real*4 dvpbf,pnlhdr,xarray_s
      !integer*4 kfil_s,keof_s,nphdr_s,nphdrf,nlimbf


      kfil_s = KFILE
      keof_s = KEOF

      CALL BUFIN_sgl (kfil_s,keof_s,pnlhdr(1),nphdrf)

      KEOF = keof_s
      IF (KEOF.LE.0) RETURN

      CALL BUFIN_sgl (kfil_s,keof_s,xarray_s(1),nlimbf)

      KEOF = keof_s

      V1P = V1PBF
      V2P = V2PBF
      DVP = dvpbf
      NLIM = nlimbf

      ! The variable XARRAY (either single or double, depending on the
      ! complile option specified, is set equal to the real*4 variable xar
      do i=1,nlimbf
         XARRAY(i) = xarray_s(i)
      enddo

      RETURN
   END SUBROUTINE


!-----------------------------------------------------------------------
!     SUBROUTINE SOLIN_tri inputs files for use with solar radiation
!     calculations for interpolation in SOLINT.  Reads files with
!     three records per panel. SOLIN_tri reads single precision files
!
!                  Written  :         January 2017
!
!                  IMPLEMENTATION:    K. E. Cady-Pereira
!-----------------------------------------------------------------------
   SUBROUTINE SOLIN_tri (V1P,V2P,DVP,NLIM,KFILE,XARRAY,KEOF)
!-----------------------------------------------------------------------
      real*8  ,intent(out)   :: V1P,V2P
      real    ,intent(out)   :: DVP
      integer ,intent(out)   :: NLIM
      integer ,intent(in)    :: KFILE
      real    ,intent(out)   :: XARRAY(:,:)
      integer ,intent(inout) :: KEOF

      integer   :: ivar
      real*4    :: xarray_s(2410)
      integer*4 :: kfil_s, keof_s, nphdr_s, nphdrf
      data nphdrf / 6 /

      real*8    :: V1PBF,V2PBF
      real*4    :: dvpbf
      integer*4 :: nlimbf
      real*4    :: PNLHDR(2)
      COMMON /BUFPNL_s/ V1PBF,V2PBF,dvpbf,nlimbf
      EQUIVALENCE (PNLHDR(1),V1PBF)


      kfil_s = KFILE
      keof_s = KEOF

      CALL BUFIN_sgl (kfil_s,keof_s,pnlhdr(1),nphdrf)

      KEOF = keof_s
      IF (KEOF.LE.0) RETURN

      ! The variable XARRAY (either single or double, depending on the
      ! complile option specified, is set equal to the real*4 variable xar
      do ivar=1,3
         CALL BUFIN_sgl (kfil_s,keof_s,xarray_s(1),nlimbf)
         xarray(1:nlimbf,ivar) = xarray_s(1:nlimbf)
         !if (ivar.eq.1) then
         !   do il=1,nlimbf
         !      print *, xarray(il,ivar)
         !   end do
         !endif
      end do

      KEOF = keof_s

      V1P = V1PBF
      V2P = V2PBF
      DVP = dvpbf
      NLIM = nlimbf

      RETURN
      END SUBROUTINE




END MODULE
