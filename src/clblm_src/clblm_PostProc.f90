!
! CREATION HISTORY:
!       Written by:     Yingtao Ma, AER@NOAA/NESDIS/STAR
!                       yma@aer.com; yingtao.ma@noaa.gov
!

MODULE Module_PostProc

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: postProc_scan, &
             postProc_filter, &
             expandBound, &
             applySRF, &
             CLBLM_ChannelSRF, &
             CLBLM_FilterFunct, &
             CLBLM_FilterFunct_init, &
             readFilterFunct, &
             writeFilterFunct


   ! If spectral scan is requested, the limits of the scanned spectrum 
   ! are expanded by HWHM*expandBound(functID). 
   !    functID for Scanning Function                             
   !    =  1: Boxcar                            
   !    =  2: Triangle                                             
   !    =  3: Gauss                                                
   !    =  4: Sinc2                                                
   !    =  5: Sinc                                                 
   !    =  6: Beer                                                 
   !    =  7: Hamming                                              
   !    =  8: Hanning                                              
   !    =  9: Norton-Beer, Weak                                    
   !    = 10: Norton-Beer, Moderate                                
   !    = 11: norton-Beer, Strong                                  
   !    = 12: Brault (needs input parameter PARM)                  
   !    = 13: Kaiser-Bessel (needs input parameter PARM)           
   !    = 14: Kiruna (asymetric, c1*sinc(u)+c2*sinc(u-u1)     
   ! * Please be noted that functID is differ from the value of JFN 
   !   used internally in clblm_FFTSCN and clblm_SCANFN. JFN = functID-1
   !   
   integer ,PARAMETER :: NumFunct = 14
   real    ,PARAMETER :: expandBound(NumFunct) = [&
                         1., 2., 4., 55., 160., 20., 20., 20., &
                         40., 40., 20., 100., 10., 160.]

                         
   TYPE :: CLBLM_ChannelSRF
      integer            :: chIndx
      character(12)      :: chID       = 'NA'
      integer            :: spectUnit  = 0    !Default: 0, ‘cm-1’; others: 1,’GHz’; 2,’mu’ (microns)
      integer            :: polarStat  = -1   !If =-1, unpoloarized
      real               :: centerFreq
      real               :: startFreq
      real               :: deltaFreq
      integer            :: nSamp             !Number of sampling points
      real,  allocatable :: filtFunct(:)      ![nSampl] filter function value
   END TYPE
   TYPE :: CLBLM_FilterFunct
      character(256)                      :: info         !Any userful sensor information.
      integer                             :: WMO_satID
      integer                             :: WMO_sensorID
      integer                             :: nChan = 0    !number of sensor channels
      type(CLBLM_ChannelSRF) ,allocatable :: SRF(:)
   END TYPE

   

   INTERFACE postProc_scan
      module procedure postProc_scan_1D
      module procedure postProc_scan_2D
      module procedure postProc_scan_struct
   END INTERFACE

   INTERFACE postProc_filter
      module procedure postProc_filter_singleStruct
      module procedure postProc_filter_structArray
      module procedure postProc_filter_array
   END INTERFACE
   
CONTAINS !=======================Module Contains========================             



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE CLBLM_FilterFunct_init( this, nChan, nSampl, &
                                      chIndx, chID, spectUnit, polarStat, &
                                      info, WMO_satID, WMO_sensorID )
!-----------------------------------------------------------------------
      type(CLBLM_FilterFunct) ,intent(out) :: this
      integer                 ,intent(in)  :: nChan
      integer                 ,intent(in)  :: nSampl(:)
      integer        ,optional,intent(in)  :: chIndx(:)
      character(*)   ,optional,intent(in)  :: chID(:)
      integer        ,optional,intent(in)  :: spectUnit
      integer        ,optional,intent(in)  :: polarStat(:)
      character(*)   ,optional,intent(in)  :: info
      integer        ,optional,intent(in)  :: WMO_satID
      integer        ,optional,intent(in)  :: WMO_sensorID

      integer :: ic

      !--- Deallocate the current object, if exist.
      !if (allocated(this%SRF)) then
      !   do ic =1,this%nChan
      !      call CLBLM_ChannelSRF_final( this%SRF(ic) )
      !   enddo
      !   deallocate(this%SRF)
      !endif
      if (allocated(this%SRF)) deallocate(this%SRF)


      !-- Allocate the filter function structure.
      allocate( this%SRF(nChan) )
      do ic =1,nChan
                                 this%SRF(ic)%chIndx    = ic
         if (present(chIndx))    this%SRF(ic)%chIndx    = chIndx(ic)
         if (present(chID))      this%SRF(ic)%chID      = chID(ic)
         if (present(spectUnit)) this%SRF(ic)%spectUnit = spectUnit
         if (present(polarStat)) this%SRF(ic)%polarStat = polarStat(ic)

         allocate( this%SRF(ic)%filtFunct( nSampl(ic) ))
      enddo

      if (present(info))         this%info = info
      if (present(WMO_satID))    this%WMO_satID = WMO_satID
      if (present(WMO_sensorID)) this%WMO_sensorID = WMO_sensorID
      this%nChan = nChan

   END SUBROUTINE


   
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProc_scan_2D( nSpectr, V1in, V2in, DVin, V1out, V2out, DVout, &
                                 data_in, convolution_flags, &
                                 data_out, padding_value )                              
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Config      ,ONLY: CLBLM_Post_Ctrl
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum,&
                                    clblm_spectrum_init
      IMPLICIT NONE
     
      integer               ,intent(in)    :: nSpectr
      real(r8)              ,intent(in)    :: V1in, V2in
      real                  ,intent(in)    :: DVin
      real(r8)              ,intent(inout) :: V1out, V2out
      real                  ,intent(inout) :: DVout
      real                  ,intent(in)    :: data_in (:,:)      ![npts,nSpectr]
      type(CLBLM_Post_Ctrl) ,intent(in)    :: convolution_flags
      real                  ,intent(out)   :: data_out (:,:)       ![:,nSpectr)
      real        ,optional ,intent(in)    :: padding_value
     
      character(*), parameter :: routineName='postProc_scan'      
      integer                 :: is, npIn
      type(CLBLM_Spectrum)    :: spect


      !--- Loop over spectra
      do is = 1,nSpectr

         npIn = size(data_in,1)      
         call CLBLM_Spectrum_init( spect, V1in, DVin, npIn )            
         spect%spect(1:npIn) = data_in(1:npIn, is)
         
         call postProc_scan_struct( spect, V1out,V2out,DVout, &
                                     convolution_flags, padding_value )
                              
         !--- Copy data to the output array
         data_out(1:spect%NLIM,is) = spect%spect(1:spect%NLIM) !move_alloc( spect%spect, data_out )

      enddo

      V1out = spect%V1
      V1out = spect%V2
      DVout = spect%DV
                  
   END SUBROUTINE




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProc_scan_1D( V1in, V2in, DVin, V1out, V2out, DVout, &
                                 data_in, convolution_flags, &
                                 data_out, padding_value )                              
!-----------------------------------------------------------------------
      USE Module_ConstParam  ,ONLY: r8=>kind_r8
      USE Module_Spectrum    ,ONLY: CLBLM_Spectrum,&
                                    clblm_spectrum_init
      USE Module_Config      ,ONLY: CLBLM_Post_Ctrl
      IMPLICIT NONE
     
      real(r8)              ,intent(in)    :: V1in, V2in
      real                  ,intent(in)    :: DVin
      real(r8)              ,intent(inout) :: V1out, V2out
      real                  ,intent(inout) :: DVout
      real                  ,intent(in)    :: data_in(:)      ![npts,nSpectr], may be deallocated
      type(CLBLM_Post_Ctrl) ,intent(in)    :: convolution_flags
      real                  ,intent(out)   :: data_out(:)       ![:,nSpectr)
      real        ,optional ,intent(in)    :: padding_value
      
      character(*), parameter :: routineName='postProc_scan'      
      integer                 :: nPtsIn
      type(CLBLM_Spectrum)    :: spect
      real                    :: filtOut 

      
      nPtsIn = size(data_in)
      call CLBLM_Spectrum_init( spect, V1in, DVin, nPtsIn )            
      spect%spect(1:nPtsIn) = data_in(1:nPtsIn) !move_alloc( data_in, spect%spect )
      
      call postProc_scan_struct( spect, V1out,V2out,DVout, &
                                  convolution_flags, padding_value )
                           
      !--- Copy data to the output array
      data_out(1:spect%NLIM) = spect%spect(1:spect%NLIM) !move_alloc( spect%spect, data_out )
      V1out = spect%V1
      V1out = spect%V2
      DVout = spect%DV
      
   END SUBROUTINE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProc_scan_struct( spect, V1out,V2out,DVout, &
                                     postCtrl, paddingValue )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum
      USE Module_Config         ,ONLY: CLBLM_Post_Ctrl
      USE Module_FFT            ,ONLY: clblm_FFTSCN
      USE Module_ScanFilter     ,ONLY: clblm_SCANFN, &
                                       clblm_FLTRFN
      IMPLICIT NONE
     
      type(CLBLM_Spectrum)   ,intent(inout) :: spect
      real(r8)               ,intent(inout) :: V1out,V2out
      real                   ,intent(inout) :: DVout
      type(CLBLM_Post_Ctrl)  ,intent(in)    :: postCtrl
      real          ,OPTIONAL,intent(in)    :: paddingValue
      
      character(*), parameter :: routineName='postProc_scan'
      integer  :: ic
      real(r8) :: V1filt
      real     :: DVfilt
      real     :: filtOut

            
      IF ( postCtrl%functID /=0 ) then !not filtering, do scanning 
      
         if ( postCtrl%FFT ) then !FFT scan
            
            !--- 
            ! FFT needs the length of the working record to be power of 2,
            ! If necessary, the end of the last record with be filled with
            ! padding value 0's (for radiance) or 1's (for transmittance)
            if (.not.present(paddingValue)) STOP '--- '//routineName//'(): Padding vlaue must be specified when calling FFT scan.'
            
            CALL clblm_FFTSCN( spect, &
                               V1out, V2out, DVout, &
                               postCtrl%functID, &
                               postCtrl%HWHM, &
                               paddingValue, &
                               postCtrl%boxcarHW, &
                               postCtrl%deconvPreScan, &
                               postCtrl%functParams(:) )
                                           
         else! ( .NOT.postCtrl%FFT ) then !Scan in wavenumber space
         
            CALL clblm_SCANFN( spect, &
                               V1out, V2out, DVout, &
                               postCtrl%functID, &
                               postCtrl%HWHM )
            
         endif  !if ( .NOT.postCtrl%FFT )
      
      ENDIF
      
   END SUBROUTINE
   

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProc_filter_singleStruct( spect, filtFunct, dataOut )
!-----------------------------------------------------------------------
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum
      USE Module_ScanFilter     ,ONLY: clblm_FLTRFN
      IMPLICIT NONE
     
      type(CLBLM_Spectrum)    ,intent(in)  :: spect
      type(CLBLM_FilterFunct) ,intent(in)  :: filtFunct
      real                    ,intent(out) :: dataOut(:)  ![nChan]
      
      character(*), parameter :: routineName='postProc_filter'
      integer  :: is,ic,NPsrf
      real     :: V1srf, DVsrf, filtOut


      !TYPE :: CLBLM_ChannelSRF
      !   integer            :: chIndx
      !   character(12)      :: chID       = 'NA'
      !   integer            :: spectUnit  = 0    !Default: 0, ‘cm-1’; others: 1,’GHz’; 2,’mu’ (microns)
      !   integer            :: polarStat  = -1   !If =-1, unpoloarized
      !   real               :: centerFreq
      !   real               :: startFreq
      !   real               :: deltaFreq
      !   integer            :: nSamp
      !   real,  allocatable :: filtFunct(:)      ![nSamp] filter function value
      !END TYPE
      !TYPE :: CLBLM_FilterFunct
      !   character(256)                      :: info         !Any userful sensor information.
      !   integer                             :: WMO_satID
      !   integer                             :: WMO_sensorID
      !   integer                             :: nChan = 0    !number of sensor channels
      !   type(CLBLM_ChannelSRF) ,allocatable :: SRF(:)
      !END TYPE
      
      do ic = 1,filtFunct%nChan
      
         V1srf = filtFunct%SRF(ic)%startFreq
         DVsrf = filtFunct%SRF(ic)%deltaFreq
         NPsrf = filtFunct%SRF(ic)%nSamp
         
         CALL clblm_FLTRFN ( spect, &
                             V1srf, DVsrf, NPsrf, &
                             filtFunct%SRF(ic)%filtFunct, &
                             filtOut )
                             
         dataOut(ic) = filtOut         
      enddo      
   END SUBROUTINE
   
     
   
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProc_filter_structArray( spectArray, nSpect, filtFunct, dataOut )
!-----------------------------------------------------------------------
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum
      USE Module_ScanFilter     ,ONLY: clblm_FLTRFN
      IMPLICIT NONE
     
      type(CLBLM_Spectrum)    ,intent(in)  :: spectArray(:) ![nSpect]
      integer                 ,intent(in)  :: nSpect
      type(CLBLM_FilterFunct) ,intent(in)  :: filtFunct
      real                    ,intent(out) :: dataOut(:,:)  ![nChan,nSpect]
      
      character(*), parameter :: routineName='postProc_filter'
      integer  :: is,ic,NPsrf
      real     :: V1srf, DVsrf, filtOut


      !TYPE :: CLBLM_ChannelSRF
      !   integer            :: chIndx
      !   character(12)      :: chID       = 'NA'
      !   integer            :: spectUnit  = 0    !Default: 0, ‘cm-1’; others: 1,’GHz’; 2,’mu’ (microns)
      !   integer            :: polarStat  = -1   !If =-1, unpoloarized
      !   real               :: centerFreq
      !   real               :: startFreq
      !   real               :: deltaFreq
      !   integer            :: nSamp
      !   real,  allocatable :: filtFunct(:)      ![nSamp] filter function value
      !END TYPE
      !TYPE :: CLBLM_FilterFunct
      !   character(256)                      :: info         !Any userful sensor information.
      !   integer                             :: WMO_satID
      !   integer                             :: WMO_sensorID
      !   integer                             :: nChan = 0    !number of sensor channels
      !   type(CLBLM_ChannelSRF) ,allocatable :: SRF(:)
      !END TYPE
      
      do is = 1,nSpect         
      do ic = 1,filtFunct%nChan
      
         V1srf = filtFunct%SRF(ic)%startFreq
         DVsrf = filtFunct%SRF(ic)%deltaFreq
         NPsrf = filtFunct%SRF(ic)%nSamp
         
         CALL clblm_FLTRFN ( spectArray(is), &
                             V1srf, DVsrf, NPsrf, &
                             filtFunct%SRF(ic)%filtFunct, &
                             filtOut )
                             
         dataOut(ic,is) = filtOut
         
      enddo
      enddo

      
   END SUBROUTINE
   
   
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE postProc_filter_array( nSpectr, V1in, V2in, DVin, data_in, &
                                 nChan, V1SRF, DVSRF, nSRFPts, SRF, data_out )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum, &
                                       CLBLM_Spectrum_init
      USE Module_ScanFilter     ,ONLY: clblm_FLTRFN
      IMPLICIT NONE
      
      integer  ,intent(in)  :: nSpectr
      real(r8) ,intent(in)  :: V1in,V2in
      real     ,intent(in)  :: DVin
      real     ,intent(in)  :: data_in(:,:)         ![nPts,nSpectr]
      integer  ,intent(in)  :: nChan
      real     ,intent(in)  :: V1SRF(:)             ![nChan]
      real     ,intent(in)  :: DVSRF(:)             ![nChan]
      integer  ,intent(in)  :: nSRFPts (:)          ![nChan]
      real     ,intent(in)  :: SRF(:,:)             ![nSRFPts,nChan]
      real     ,intent(out) :: data_out(:,:)        ![nChan,nSpectr]
      
      type(CLBLM_Spectrum) :: spect
      integer :: is,ic,NPin
      real    :: filtOut
      
      
      do is = 1,nSpectr
      
         NPin = ceiling( (V2in-V1in)/DVin + 1. )
         NPin = min( NPin,size(data_in,is) )
         
         call CLBLM_Spectrum_init( spect, V1in,DVin,NPin)
         spect%spect(1:NPin) = data_in(1:NPin,is)
         
         do ic = 1,nChan
         
            call clblm_FLTRFN( spect, V1SRF(ic), DVSRF(ic), nSRFPts(ic), &
                               SRF(1:nSRFPts(ic),ic), filtOut )
                                
            data_out(ic,is) = filtOut         
         enddo
      enddo

   END SUBROUTINE

   
   
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE applySRF( nSpectr, V1in, V2in, DVin, data_in, &
                        V1_SRF, DV_SRF, nSRFpts, SRF, data_out )
!-----------------------------------------------------------------------
      USE Module_ConstParam     ,ONLY: r8=>kind_r8
      USE Module_Spectrum       ,ONLY: CLBLM_Spectrum, &
                                       CLBLM_Spectrum_init
      USE Module_ScanFilter     ,ONLY: clblm_FLTRFN
      IMPLICIT NONE
      
      integer  ,intent(in)  :: nSpectr
      real(r8) ,intent(in)  :: V1in,V2in
      real     ,intent(in)  :: DVin
      real     ,intent(in)  :: data_in(:,:)  ![nPts,nSpectr]
      real     ,intent(in)  :: V1_SRF
      real     ,intent(in)  :: DV_SRF
      integer  ,intent(in)  :: nSRFpts
      real     ,intent(in)  :: SRF(:)        ![nSRFpts]
      real     ,intent(out) :: data_out(:)   ![nSpectr]
      
      type(CLBLM_Spectrum) :: spect
      integer :: NPin, is
      real    :: filtOut
      
      
      do is = 1,nSpectr
      
         NPin = ceiling( (V2in-V1in)/DVin + 1. )
         NPin = min( NPin,size(data_in,is) )
         
         call CLBLM_Spectrum_init( spect, V1in,DVin,NPin)
         spect%spect(1:NPin) = data_in(1:NPin,is)
         
         call clblm_FLTRFN ( spect, V1_SRF, DV_SRF, nSRFpts, SRF, filtOut )

         data_out(is) = filtOut         
      enddo

   END SUBROUTINE
   
   
   
!-----------------------------------------------------------------------
! Read in instrument response functions from a netcdf file
!-----------------------------------------------------------------------
   SUBROUTINE readFilterFunct(filename, filtFunct)
!-----------------------------------------------------------------------
   USE NETCDF
   USE Module_FileIO ,ONLY: check             

   character(len=*)         ,intent(in)  :: filename
   type(CLBLM_FilterFunct)  ,intent(out) :: filtFunct
      
   integer(4)              :: ncid, grpID, sampleid, chindxid, &
                              spectunitid, polarstatid, centerfreqid, &
                              startfreqid, deltafreqid, ffuncid
   integer(4) ,allocatable :: groupIDs(:)
   integer                 :: ig
   integer(4)              :: nSamp, numgrps
   character               :: tempStr*12

   !TYPE :: CLBLM_ChannelSRF
   !   integer            :: chIndx
   !   character(12)      :: chID       = 'NA'
   !   integer            :: spectUnit  = 0    !Default: 0, ‘cm-1’; others: 1,’GHz’; 2,’mu’ (microns)
   !   integer            :: polarStat  = -1   !If =-1, unpoloarized
   !   real               :: centerFreq
   !   real               :: startFreq
   !   real               :: deltaFreq
   !   integer            :: nSamp
   !   real,  allocatable :: filtFunct(:)      ![nSamp] filter function value
   !END TYPE
   !TYPE :: CLBLM_FilterFunct
   !   character(256)                      :: info         !Any userful sensor information.
   !   integer                             :: WMO_satID
   !   integer                             :: WMO_sensorID
   !   integer                             :: nChan = 0    !number of sensor channels
   !   type(CLBLM_ChannelSRF) ,allocatable :: SRF(:)
   !END TYPE
   
   
   call check( nf90_open(trim(filename), NF90_NOWRITE, ncid ))

   call check( nf90_get_att(ncid, NF90_GLOBAL, "FilterInformation", filtFunct%info ))
   call check( nf90_get_att(ncid, NF90_GLOBAL, "WMO_satID",         filtFunct%WMO_satID ))
   call check( nf90_get_att(ncid, NF90_GLOBAL, "WMO_sensorID",      filtFunct%WMO_sensorID))
   call check( nf90_get_att(ncid, NF90_GLOBAL, "numChannels",       filtFunct%nChan) )
   
   allocate( filtFunct%SRF( filtFunct%nChan ))
      
   allocate( groupIDs( filtFunct%nChan ) )
   call check( nf90_inq_grps(ncid, numgrps, groupIDs) )
   
   do ig = 1,filtFunct%nChan
   
      grpID = groupIDs(ig)
      call check( nf90_inq_dimid( grpID, 'NSAMP', sampleid ))
      call check( nf90_inquire_dimension( grpID, sampleid, len=nSamp ))
    
      if (allocated( filtFunct%SRF(ig)%filtFunct )) deallocate( filtFunct%SRF(ig)%filtFunct )
      allocate( filtFunct%SRF(ig)%filtfunct( nSamp ))
      filtFunct%SRF(ig)%nSamp = nSamp    

      call check( nf90_get_Att( grpID, NF90_GLOBAL, "CHANNEL_ID",  tempstr)); filtFunct%SRF(ig)%chid = trim(tempstr)
      call check( nf90_get_att( grpID, NF90_GLOBAL, "CHINDX",      filtFunct%SRF(ig)%chindx))
      call check( nf90_get_att( grpID, NF90_GLOBAL, "SPECT_UNIT",  filtFunct%SRF(ig)%spectunit))
      call check( nf90_get_att( grpID, NF90_GLOBAL, "POLAR_STAT",  filtFunct%SRF(ig)%polarstat))
      call check( nf90_get_att( grpID, NF90_GLOBAL, "CENTER_FREQ", filtFunct%SRF(ig)%centerFreq))
      call check( nf90_get_att( grpID, NF90_GLOBAL, "START_FREQ",  filtFunct%SRF(ig)%startFreq))
      call check( nf90_get_att( grpID, NF90_GLOBAL, "DELTA_FREQ",  filtFunct%SRF(ig)%deltafreq))

      call check( nf90_inq_varid( grpID, "FILTER_FUNCTION", ffuncid))

      call check( nf90_get_var( grpID, ffuncid, filtFunct%SRF(ig)%filtfunct))         
   enddo
   
   call check(nf90_close(ncid))
   END SUBROUTINE

!-----------------------------------------------------------------------
! Write instrument response functions to netcdf file
!-----------------------------------------------------------------------
   SUBROUTINE writeFilterFunct(filename, filtFunct)
!-----------------------------------------------------------------------
   USE NetCDF
   USE Module_FileIO ,ONLY: check, createNetCDF4File

   character(len=*)        ,intent(in) :: filename
   type(CLBLM_FilterFunct) ,intent(in) :: filtFunct

   integer(4)              :: ncid, sampleid, chindxid, spectunitid, polarstatid, &
                              centerfreqid, startfreqid, deltafreqid, ffuncid
   integer(4) ,allocatable :: grpID(:)
   integer                 :: ig
   character               :: grpName*64
    
   !TYPE :: CLBLM_ChannelSRF
   !   integer            :: chIndx
   !   character(12)      :: chID       = 'NA'
   !   integer            :: spectUnit  = 0    !Default: 0, ‘cm-1’; others: 1,’GHz’; 2,’mu’ (microns)
   !   integer            :: polarStat  = -1   !If =-1, unpoloarized
   !   real               :: centerFreq
   !   real               :: startFreq
   !   real               :: deltaFreq
   !   integer            :: nSamp
   !   real,  allocatable :: filtFunct(:)      ![nSamp] filter function value
   !END TYPE
   !TYPE :: CLBLM_FilterFunct
   !   character(256)                      :: info         !Any userful sensor information.
   !   integer                             :: WMO_satID
   !   integer                             :: WMO_sensorID
   !   integer                             :: nChan = 0    !number of sensor channels
   !   type(CLBLM_ChannelSRF) ,allocatable :: SRF(:)
   !END TYPE
       
   if (.not. allocated( filtFunct%SRF )) then
      write(*,*) 'error, cant write blank function, filtfunct not allocated'
      stop
   endif
    
   call createNetCDF4File( trim(filename), ncid, enddef=.FALSE. )

   call check( nf90_put_att(ncid, NF90_GLOBAL, "FilterInformation", filtFunct%info ))
   call check( nf90_put_att(ncid, NF90_GLOBAL, "WMO_satID",         filtFunct%WMO_satID ))
   call check( nf90_put_att(ncid, NF90_GLOBAL, "WMO_sensorID",      filtFunct%WMO_sensorID))
   call check( nf90_put_att(ncid, NF90_GLOBAL, "numChannels",       filtFunct%nChan) )
   
   allocate( grpID( filtFunct%nChan ))
   do ig = 1, filtFunct%nChan
      
      !--- Create a group
      write( grpName, '("ChennelSRF-",I3.3)') ig
      call check( nf90_def_grp(ncid, trim( grpName ), grpID(ig) ) )
      
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "CHANNEL_ID",  trim(filtFunct%SRF(ig)%chid)))
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "CHINDX",      filtFunct%SRF(ig)%chindx))
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "SPECT_UNIT",  filtFunct%SRF(ig)%spectunit))
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "POLAR_STAT",  filtFunct%SRF(ig)%polarstat))
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "CENTER_FREQ", filtFunct%SRF(ig)%centerFreq))
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "START_FREQ",  filtFunct%SRF(ig)%startFreq))
      call check( nf90_put_att(  grpID(ig), NF90_GLOBAL, "DELTA_FREQ",  filtFunct%SRF(ig)%deltafreq))
                                                                       
      call check( nf90_def_dim(  grpID(ig), "NSAMP", int(filtFunct%SRF(ig)%nSamp,4), sampleid ))      
      call check( nf90_def_var(  grpID(ig), "FILTER_FUNCTION", NF90_DOUBLE, (/ sampleid /), ffuncid))
   enddo
     
   call check( nf90_enddef( ncid ))

   do ig = 1, filtFunct%nChan   
      call check( nf90_put_var( grpID(ig), ffuncid, filtFunct%SRF(ig)%filtfunct(1:filtFunct%SRF(ig)%nsamp)))
   enddo       

   END SUBROUTINE
    
END MODULE 

