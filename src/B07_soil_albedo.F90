! Carrer_soilalbedo_to_GISS.f
! @auth 10/14/2018 N.Y.Kiang
      
!     Takes Carrer soil albedo annual (min,mean,max, and std) already interpolated to desired spatial resolution and generates GISS GCM boundary conditions input versions.  Creates downscaled 1 km bare soil bright fraction ('bs_brightratio.nc') for ModelE grey bare soil albedo scheme.
!     Versions:
!  1) Spectral albedo TOA, GISS 6 bands.  This is the Carrer albedos from MODIS VIS (300-700 nm) and NIR (700-5000 nm) translated to GISS 6 bands, utilizing Judith Lean's top-of-the-atmosphere 1850 preindustrial solar irradiance spectrum(S0=1360.6718 W/m2).
!  2) Spectral albedo Surf, GISS 6 bands.  Like TOA, but using a surface irradiance spectrum from
!     a) NREL's Direct+circumsolar irradiance
!     b) Brian Cairns (GISS) annual averages for: i) 60S-60N, ii) arctic zonal, iii) mid-latitude zonal, and iv) tropical zonal.
!  3) Bare soil bright fraction 'bs_brightratio.nc' for GISS ModelE grey soil
!        albedo scheme.

!###  To do:  Need to pull out the spectral code from hntr4_Carrer_nc.f
!###  To do:  Doublecheck MODIS NIR 700-4000 nm vs. Carrer 700-5000 nm??      
!###  To do:  May want to correct spectral weighting for 300-400 nm to be lower.
!###  To do:  Get Brian Cairns zonal spectral irradiances.

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

! Use the filled-in version of Carrer MODIS albedo
#define USE_FILLED

!------------------------------------------------------------------------------
      
program Carrer_soilalbedo_to_GISS


    use carrer_mod
    use netcdf
    use chunker_mod
    use ent_labels_mod
    use ent_params_mod
    use assign_laimax_mod
    use hntr_mod

!    use convertnc
!    use netcdf
!    use EntGVSD_netcdf_util
    
    implicit none

    
!    integer, parameter :: NBANDS_GISS = 6

    !integer, parameter :: NSPEC = 3785 !Number of rows in Lean spectrum
    !character*118 :: GISSBANDStxt = !len = 6*19 + 4
    !Below parameters are moved to EntGVSD_util.f
!     integer, parameter :: NBANDS_GISS = 6 !GISS 6 spectral bands
!      character*101, parameter :: GISSBANDStxt =  &
!          'VIS (330-770), NIR1 (770-860), NIR2 (860-1250), '// &
!          'NIR3 (1250-1500), NIR4 (1500-2200), NIR5 (2200-4000)'


    ! BISS Band Defintions
    integer, parameter :: VIS_GISS = 1
    integer, parameter :: NBANDS_GISS = 6
    character*4, parameter :: sbands_giss(NBANDS_GISS) = &
        (/"VIS ", "NIR1", "NIR2", "NIR3", "NIR4", "NIR5"/)
    character*(*), parameter :: sbands_giss_long(NBANDS_GISS) = &
        (/"VIS (300-770nm)   ", &
          "NIR1 (770-860nm)  ", &
          "NIR2 (860-1250nm) ", &
          "NIR3 (1250-1500nm)", &
          "NIR4 (1500-2200nm)", &
          "NIR5 (2200-4000nm)"/)
    real*8, parameter :: bandpoints_giss(NBANDS_GISS+1) = &
        (/300,770,860,1250,1500,2200,4000/)

    !GISS/MODIS band sections
    integer, parameter :: nm300_400 = 1
    integer, parameter :: nm400_700 = 2
    integer, parameter :: nm700_770 = 3
    integer, parameter :: nm770_860 = 4
    integer, parameter :: nm860_1250 = 5
    integer, parameter :: nm1250_1500 = 6
    integer, parameter :: nm1500_2200 = 7
    integer, parameter :: nm2200_4000 = 8

    !GISS ModelE grey albedo fractions
    integer, parameter :: BRIGHT = 1
    integer, parameter :: DARK = 2
    character*6, parameter :: sbright_dark(2) = (/"bright", "dark  "/)
    character*11, parameter :: sbright_dark_long(2) = (/"bright soil", "dark soil  "/)

    ! Convert spectrum in one set of bands to another set of bands while conserving energy
    !From Judith Lean 0 km solar surface irradiance (from 2006 version)
    !Fraction of shortwave 300-4000 nm in MODIS & GISS band portions.
    real*4, parameter :: fracSW_MG(8) =  (/ &
        0.105992883, &        !1 300-400 nm
        0.992030205,&	!2 400-700 nm
        0.174912449,&	!3 700-770 nm
        0.203020351,&	!4 770-860 nm
        0.477761295,&	!5 860-1250 nm
        0.089450965,&	!6 1250-1500 nm
        0.17866227,&  !7 1500-2200 nm
        0.051105118 & !8 2200-4000 nm
      /)

      ! Fraction of shortwave (SW) surface irradiance (200-4000 nmm)
      !  in GISS bands.
      ! From solar.lean_TOA.0km_toNancyin2006.xlsx.
      real*4, parameter :: fracSW_Lean_0km_2006(NBANDS_GISS) = (/ &
          0.56, 0.0893, 0.2102, 0.0394, 0.0786, 0.0225  /)
      

    real*4 :: fracSW_GISS(NBANDS_GISS) !fraction of SW W m-2 nm-1 300-4000 nm

    integer :: band(NBANDS_GISS)

    integer :: k, m, b, dimlon, dimlat, dim(2), dimband
    integer :: dg
    real*4 :: s

    type(Chunker_t) :: chunker,chunkerhr,chunkere
    ! Input files
    type(ChunkIO_t) :: ioall_lc, io_lcice, io_lcwater
    real*4 :: lcice, lcwater
    type(ChunkIO_t) :: io_albmodis(NBANDS_MODIS)
    real*4 :: albmodis(NBANDS_MODIS)
    ! Output Files
    type(ChunkIO_t) :: io_albsw
    real*4 :: albsw
    type(ChunkIO_t) :: io_albgiss(NBANDS_GISS), io_albgisse(NBANDS_GISS)
    real*4 :: albgiss(NBANDS_GISS)
    integer, parameter :: BRIGHT_DARK = 2
    type(ChunkIO_t) :: ioall_fracbd(NBANDS_GISS), io_fracbd(BRIGHT_DARK,NBANDS_GISS)
    type(ChunkIO_t) :: ioall_fracgrey, io_fracgrey(BRIGHT_DARK)
    type(ChunkIO_t) :: ioall_fracgreye, io_fracgreye(BRIGHT_DARK)
    !type(ChunkIO_t) :: ioall_fracgrey_hr, io_fracgrey_hr(BRIGHT_DARK)
    type(ChunkIO_t) :: io_bs_brightratio_hr
    real*4 :: bsbr

    ! Regrid to ModelE
    type(HntrSpec_t) :: spec_e, spec_k
    type(HntrCalc_t) :: hntr_e    ! Preparation to regrid


    INTEGER :: ichunk,jchunk,ic,jc
    integer :: i1km,j1km   ! Indexing into hi-res 1km variable
    integer :: iband,istat
    real, dimension(:,:), allocatable, target :: wta,wta_fracbd,wta_fracbd_hr    ! Surmised landmask
    real, dimension(:,:), allocatable, target :: wta1
    type(EntSet_t) :: ent2
    type(FileInfo_t) :: info
    integer, parameter :: nhr(2) = (/ IM1km/IMK, JM1km/JMK /)   ! # of hi-res gridcells per lo-res
    type(ReadWrites_t) :: rw
    call rw%init(THIS_OUTPUTS_DIR, 'B07_soil_albedo', 30,30)

    call init_ent_labels

    !* Select spectral band irradiance fractions
    fracSW_GISS(:) = fracSW_Lean_0km_2006(:)


    call chunker%init(IMK, JMK, IMH*2,JMH*2, 'forplot', 100, 100, 10,(/18,15/), outputs_dir=THIS_OUTPUTS_DIR)
    call chunkere%init(IM2, JM2, IM2,JM2, 'forplot', 100, 100, 10,(/18,15/), &
        outputs_dir=THIS_OUTPUTS_DIR) ! ModelE resolution
    call chunkerhr%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 100, 100, 10,(/18,15/), outputs_dir=THIS_OUTPUTS_DIR)

    spec_e = hntr_spec(chunkere%chunk_size(1), chunkere%ngrid(2), 0d0, 180d0*60d0 / chunkere%ngrid(2))
    spec_k = hntr_spec(chunker%chunk_size(1), chunker%ngrid(2), 0d0, 180d0*60d0 / chunker%ngrid(2))
    hntr_e = hntr_calc(spec_e, spec_k, 0d0)   ! datmis=0


    allocate(wta(chunker%chunk_size(1), chunker%chunk_size(2)))
    allocate(wta1(chunker%chunk_size(1), chunker%chunk_size(2)))
    allocate(wta_fracbd(chunker%chunk_size(1), chunker%chunk_size(2)))
    allocate(wta_fracbd_hr(chunkerhr%chunk_size(1), chunkerhr%chunk_size(2)))
    wta1 = 1d0

    ! ------------- Open Input Files

    ! --------- LC
    ent2 = make_ent2()
    call chunker%file_info(info, ent2, LAI_SOURCE, 'M', 'lc', 2004, 'ent17', '1.1')
    call chunker%nc_open(ioall_lc, &
        chunker%outputs_dir, trim(info%dir), trim(info%leaf)//'.nc', trim(info%vname), 0)
    call chunker%nc_reuse_var(ioall_lc, io_lcice, (/1,1,ent2%svm(SNOW_ICE)/))
    call chunker%nc_reuse_var(ioall_lc, io_lcwater, (/1,1,ent2%svm(CV_WATER)/))

    ! ------------ albmodis
#ifdef USE_FILLED
    do iband=1,NBANDS_MODIS
        ! Read from 2D NetCDF var
        call chunker%nc_open(io_albmodis(iband), chunker%outputs_dir, &
            'tmp/carrer/', &
            'albfill_'//trim(sbands_modis(iband))//'.nc', &
            'albfill_'//trim(sbands_modis(iband))//'_MEAN', 1)
    end do
#else
    do iband=1,NBANDS_MODIS
        ! Read from 3D NetCDF var
        call chunker%nc_open(io_albmodis(iband), chunker%outputs_dir, &
            'tmp/carrer/', &
            'albmodis_'//trim(sbands_modis(iband))//'.nc', &
            'albmodis_'//trim(sbands_modis(iband)), SMEAN)
    end do
#endif

    ! ===================== Open Output Files

    ! ------------ albsw
    call clear_file_info(info)
    info%vname = 'albsw'
    info%long_name = 'Carrer soil albedo shoftwave (300-4000 nm) annual mean 2004'
    info%units = '1'
    info%file_metadata_type = 'carrer'
    call chunker%nc_create1(io_albsw, weighting(wta,1d0,0d0), &
        'soilalbedo/', 'soilalbedo_2HX2_EntGVSD_v1.1_CarrerGISS_SW_annual_2004', info)

    ! ----------- albgiss
    do iband=1,NBANDS_GISS
        call clear_file_info(info)
        info%vname = 'albgiss_'//trim(sbands_giss(iband))
        info%long_name = 'Carrer soil albedo '//trim(sbands_giss_long(iband))//' annual mean 2004'
        info%units = '1'
        info%file_metadata_type = 'carrer'
        call chunker%nc_create1(io_albgiss(iband), weighting(wta,1d0,0d0), &
            'soilalbedo/', &
            'soilalbedo_5km_EntGVSD_v1.1_CarrerGISS_'//trim(sbands_giss(iband))//'_annual_2004', info)

        call chunkere%nc_create1(io_albgisse(iband), weighting(wta,1d0,0d0), &
            'soilalbedo/', &
            'soilalbedo_2HX2_EntGVSD_v1.1_CarrerGISS_'//trim(sbands_giss(iband))//'_annual_2004', info, &
            create_lr=.false.)
    end do

    ! ----------- fracbd
    do iband=1,NBANDS_GISS
        call clear_file_info(info)
        info%vname = 'fracbd_'//trim(sbands_giss(iband))
        info%long_name = 'Bright/Dark Soil (GISS '//trim(sbands_giss(iband))//' band)'
        info%units = '1'
        info%file_metadata_type = 'carrer'

        call chunker%nc_create1(ioall_fracbd(iband), weighting(wta_fracbd,1d0,0d0), &
            'soilalbedo/', &
            'soilalbedo_fracbd_5km_EntGVSD_v1.1_CarrerGISS_'//trim(sbands_giss(iband))//'_annual_2004', info, &
            sbright_dark, sbright_dark_long)
        do k=1,BRIGHT_DARK
            call chunker%nc_reuse_var( &
                ioall_fracbd(iband), io_fracbd(k,iband), &
                (/1,1,k/), weighting(wta_fracbd,1d0,0d0))
        end do
    end do

    ! -------------- fracgrey
    call clear_file_info(info)
    info%vname = 'fracgrey'
    info%long_name = 'Bright/Dark Soil in (grey avg of GISS bands)'
    info%units = '1'
    info%file_metadata_type = 'carrer'
    call chunker%nc_create1(ioall_fracgrey, weighting(wta_fracbd,1d0,0d0), &
        'soilalbedo/', &
        'soilalbedo_5km_EntGVSD_v1.1_CarrerGISS_fracgrey_annual_2004', info, &
        sbright_dark, sbright_dark_long)
    do k=1,2
        call chunker%nc_reuse_var( &
            ioall_fracgrey, io_fracgrey(k), &
            (/1,1,k/), weighting(wta_fracbd,1d0,0d0))
    end do

    call chunkere%nc_create1(ioall_fracgreye, weighting(wta_fracbd,1d0,0d0), &
        'soilalbedo/', &
        'soilalbedo_2HX2_fracgrey', info, &
        sbright_dark, sbright_dark_long, create_lr=.false.)
    do k=1,2
        call chunkere%nc_reuse_var( &
            ioall_fracgreye, io_fracgreye(k), &
            (/1,1,k/), weighting(wta_fracbd,1d0,0d0))
    end do

    ! ---------------- bs_brightratio_hr (1km resolution) 
    call clear_file_info(info)
    info%vname = 'bs_brightratio'
    info%long_name = 'Bright Fraction of Soil'
    info%units = '1'
    info%file_metadata_type = 'carrer'
    call chunkerhr%nc_create1( &
        io_bs_brightratio_hr, weighting(wta_fracbd_hr,1d0,0d0), &
        'soilalbedo/', &
        'soilalbedo_V1km_bs_brightratio', info)

    call chunker%nc_check(rw=rw)
    call chunkerhr%nc_check(rw=rw)
    call chunkere%nc_check(rw=rw)
    call rw%write_mk

#ifdef JUST_DEPENDENCIES
    STOP 0
#endif

    ! ================== Main Loop

#ifdef ENTGVSD_DEBUG
    do jchunk = dbj0,dbj1
    do ichunk = dbi0,dbi1
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif

        call chunker%move_to(ichunk,jchunk)
        call chunkerhr%move_to(ichunk,jchunk)
        call chunkere%move_to(ichunk,jchunk)
        wta = 0

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            ! ---------- Read inputs
            lcice = io_lcice%buf(ic,jc)
            lcwater = io_lcwater%buf(ic,jc)
            do iband=1,NBANDS_MODIS
                albmodis(iband) = io_albmodis(iband)%buf(ic,jc)
            end do

            ! Compute overall NetCDF index of current cell
            !ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
            !jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

            ! Infer soil mask
            ! TODO: Get this from LC instead
            if (albmodis(1) /= FillValue) then
                wta(ic,jc) = 1
            else
                wta(ic,jc) = 0
            end if


            !* Calculate total SW albedo --------------------------------------
            if ((albmodis(VIS_MODIS).eq.FillValue).or. &
                (albmodis(NIR_MODIS).eq.FillValue)) then
                albsw = FillValue
            else
                albsw = ( &
                     albmodis(VIS_MODIS) * sum(fracSW_MG(1:2)) &
                   + albmodis(NIR_MODIS) * sum(fracSW_MG(3:8)) &
                ) / sum(fracSW_MG(1:8)) 
            endif

            !* Put spectral breakdown into GISS bands ------------------------
            do iband=1,NBANDS_GISS
                albgiss(iband) = FillValue
            end do

            if ((albmodis(VIS_MODIS).eq.FillValue).or. &
                (albmodis(NIR_MODIS).eq.FillValue)) then
                do iband=1,NBANDS_GISS
                    albgiss(iband) = FillValue
                end do
            else
                albgiss(VIS_GISS) = &
                   ( albmodis(VIS_MODIS) * &
                   (fracSW_MG(nm300_400) + fracSW_MG(nm400_700)) + &
                   albmodis(NIR_MODIS)*fracSW_MG(nm700_770)) &
                   /( fracSW_MG(nm300_400) + fracSW_MG(nm400_700) + &
                   fracSW_MG(nm700_770) )
                do k=2,NBANDS_GISS
                    albgiss(k) = albmodis(NIR_MODIS)
                end do
            endif

            !------ Old VEG bare_bright vs. bare_dark soil fractions
            !* Calculate GISS ModelE grey albedo bright vs. dark fractions
            !Same netcdf names as in Ent_pfts.
            !1) Non-grey albedos: each GISS band has its own bright|dark fractions.
            !  For partial-land cells (that have ice and/or water) use
            ! soil albedo of nearest all-land neighbor.  If none, default to
            ! to cell's mixed albedo.
            do iband=1,NBANDS_GISS
            do k=1,2
                io_fracbd(k,iband)%buf(ic,jc) = FillValue
            end do
            end do
            do k=1,NBANDS_GISS
                if (albgiss(k) /= FillValue) then
                    if (abs(1d0 - (lcice+lcwater)) < 1d-5) then
                        ! ------ All permanent ice, no soil.
                        io_fracbd(BRIGHT,k)%buf(ic,jc) = FillValue
                        io_fracbd(DARK,k)%buf(ic,jc) = FillValue
                    else if (((lcice+lcwater) > 0.).and. &
                        ((lcice+lcwater) < 1.0)) then
                        ! --------- Partial ground
                        !Search for nearest grid cell up to 5 that is all ground
                        s = FillValue
#if 0
! Not needed at 6km
                        do dg=1,5
                        do jj=j-dg,j+dg
                        do ii=i-dg,i+dg
                            if ((lcice+lcwater(ii,jj)).eq.0.) then
                                !All ground found
                                s = albgiss(ii,jj,k)
                                exit
                            end if
                        end do
                        if (s.ne.FillValue) exit
                        end do
                        if (s.ne.FillValue) exit
                        end do
#endif
                        if (s.eq.FillValue) then !No nearby ground cells found 
                            !io_fracbd(BRIGHT,k)%buf(ic,jc) = FillValue
                            !io_fracbd(DARK,k)%buf(ic,jc) = FillValue
                            !*noice - If cannot assign land albedo, cut the fractional cells.
                            !  Loses Alaska peninsula
                            !fracgrey(i,j,:) = FillValue
                            !*fringeice - Or, since there is some ground, assign i,j value there.
                            !  Leaves some coastal ice fringe.
                            io_fracgrey(BRIGHT)%buf(ic,jc) = &
                                min(0.5, albgiss(k))/0.5 !
                            io_fracgrey(DARK)%buf(ic,jc) = 1.0 - io_fracgrey(BRIGHT)%buf(ic,jc)
                        else
                            io_fracbd(BRIGHT,k)%buf(ic,jc) = min(s,0.5)/0.5
                            io_fracbd(DARK,k)%buf(ic,jc) = 1.0 - io_fracbd(BRIGHT,k)%buf(ic,jc)
                        end if
                    else
                        ! -------- All ground
                        io_fracbd(BRIGHT,k)%buf(ic,jc) = min(albgiss(k),0.5)/0.5
                        io_fracbd(DARK,k)%buf(ic,jc) = 1.0 - io_fracbd(BRIGHT,k)%buf(ic,jc)
                    end if
                else
                    io_fracbd(BRIGHT,k)%buf(ic,jc) = FillValue
                    io_fracbd(DARK,k)%buf(ic,jc) = FillValue
                end if
            end do   ! k=1,NBANDS

            if (io_fracbd(1,1)%buf(ic,jc) == FillValue) then
                wta_fracbd(ic,jc) = 0
            else
                wta_fracbd(ic,jc) = 1
            end if

            !2) Grey shortwave albedos - exclude 100% ice+water, nearest-neighbor fill
            !  For partial-land cells (that have water) use
            ! soil albedo of nearest all-land/non-ice neighbor.  If none, default to
            ! to cell's mixed albedo.
            io_fracgrey(BRIGHT)%buf(ic,jc) = 0
            io_fracgrey(DARK)%buf(ic,jc) = 0

            if (albsw.ne.FillValue) then 
               if ((lcice+lcwater).eq.1.0) then !No ground
                    io_fracgrey(BRIGHT)%buf(ic,jc) = FillValue
                    io_fracgrey(DARK)%buf(ic,jc) = FillValue
               elseif (((lcice+lcwater).gt.0.).and. &
                      ((lcice+lcwater).lt.1.0)) then !Partial
                  !Search for nearest grid cell up to 5 that is all ground
                  s = FillValue
#if 0
! Not needed at 6km
                  do dg=1,5
                  do ii=i-dg,i+dg
                     do jj=j-dg,j+dg
                        if ((lcice+lcwater(ii,jj)).eq.0.) then
                          !All ground found
                          s = albSW(ii,jj)
                          exit
                       endif
                    enddo
                    if (s.ne.FillValue) exit
                  enddo
                  if (s.ne.FillValue) exit
                  enddo
#endif
                  if (s.eq.FillValue) then !No nearby ground cells found
!                    !*If cannot assign land albedo, cut the fractional cells.
!                    !Loses Alaska peninsula
                     !fracgrey(i,j,:) = FillValue
!                    !*Or, since there is some ground, assign i,j value there.
!                    !Leaves some coastal ice fringe.
                     ! *** This is where had to cap bright/dark fraction to .5
                     io_fracgrey(BRIGHT)%buf(ic,jc) = min(0.5, albsw)/0.5 !
                     io_fracgrey(DARK)%buf(ic,jc) = 1d0 - io_fracgrey(BRIGHT)%buf(ic,jc)
                  else !Found nearby all-ground cell
                     io_fracgrey(BRIGHT)%buf(ic,jc) = min(0.5,s)/0.5
                     io_fracgrey(DARK)%buf(ic,jc) = 1d0 - io_fracgrey(BRIGHT)%buf(ic,jc)
                  endif
               else    !Cell is all ground with SW
                  io_fracgrey(BRIGHT)%buf(ic,jc) = min(0.5, albsw)/0.5
                  io_fracgrey(DARK)%buf(ic,jc) = 1d0 - io_fracgrey(BRIGHT)%buf(ic,jc)
               endif
            else
                io_fracgrey(BRIGHT)%buf(ic,jc) = FillValue
                io_fracgrey(DARK)%buf(ic,jc) = FillValue
            endif

            ! ------------ Regrid fracgrey
            bsbr = io_fracgrey(BRIGHT)%buf(ic,jc)
            do j1km = (jc-1)*nhr(2)+1, jc*nhr(2)
            do i1km = (ic-1)*nhr(1)+1, ic*nhr(1)
                wta_fracbd_hr(i1km,j1km) = wta_fracbd(ic,jc)
                io_bs_brightratio_hr%buf(i1km,j1km) = bsbr
            end do
            end do


            ! ---------- Store outputs
            io_albsw%buf(ic,jc) = albsw
            do k=1,NBANDS_GISS
                io_albgiss(k)%buf(ic,jc) = albgiss(k)
            end do


        end do
        end do

        ! ------------ Regrid to GISS resolution
print *,'AA1'
        do k=1,NBANDS_GISS
            call hntr_e%regrid4( &
                io_albgisse(k)%buf, io_albgiss(k)%buf, &
                wta1, 1d0, 0d0, &    ! weighting,
                io_albgisse(k)%startB(2), io_albgisse(k)%chunker%chunk_size(2))
        end do
print *,'AA2'
        do k=1,2
            call hntr_e%regrid4( &
                io_fracgreye(k)%buf, io_fracgrey(k)%buf, &
                wta1, 1d0, 0d0, &    ! weighting,
                io_fracgreye(k)%startB(2), io_fracgreye(k)%chunker%chunk_size(2))
        end do
print *,'AA3'


        call chunker%write_chunks
        call chunkerhr%write_chunks
        call chunkere%write_chunks
    end do
    end do


    call chunker%close_chunks
    call chunkerhr%close_chunks
    call chunkere%close_chunks

end program Carrer_soilalbedo_to_GISS
            
