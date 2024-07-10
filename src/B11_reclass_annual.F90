! Reclass annual LAIMAX from 20 land cover class scheme to 18+waterice class scheme.
!
! Takes 20-cover classes (Ent 17 PFTs + 3 non-veg) and converts to
! 18-cover + water_ice classes (Ent 16 PFTs and bright and dark bare soil fractions +
! water_ice).  In v1.0, water and snow/ice were set to undef, but now they are
! made into a layer.
! Merges C3 and C4 crops into one crop cover type for Ent 16 PFTs, combines
! water and permanent ice into one cover type, and converts barse/sparse cover 
! into equivalent veg type and bare soil bright and dark fractions, preserving total LAI of grid cell.
!
! Author: Nancy Kiang, Carlo Monte, Elizabeth Fischer
!
! See slib/cropmerge_laisparse_splitbare.f90

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module b11_mod
    use netcdf
    use chunker_mod
    use ent_labels_mod
    use ent_params_mod
    use gcm_labels_mod
    use cropmerge_laisparse_splitbare_mod

implicit none

    integer, parameter :: one = 1

CONTAINS

subroutine do_reindex(esub)
    type(GcmEntSet_t), intent(IN), target :: esub

    class(EntSet_t), pointer :: esub_p

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: io_lc(NENT20)
    type(ChunkIO_t) :: io_laiin(NENT20,one)
    type(ChunkIO_t) :: io_bs
    type(ChunkIO_t) :: io_simin(NENT20,one)
    type(ChunkIO_t) :: io_TCinave
    ! Output files
    type(ChunkIO_t) :: io_lcout(esub%ncover,one)
    type(ChunkIO_t) :: io_laiout(esub%ncover,one)
    type(ChunkIO_t) :: io_lc_checksum(one)
    type(ChunkIO_t) :: io_lchgt_checksum(one)
    type(ChunkIO_t) :: io_lclai_checksum(one)
    type(ChunkIO_t) :: io_simout(esub%ncover,one)
    type(FileInfo_t) :: info, overmeta
    integer :: k,ksub
#if (defined HGT_GEDI)
    integer, parameter :: hgt_year = 2020
    character*5, parameter :: hgt_ver = '1.1.2'
#elif (defined HGT_POTAPOV)
    integer, parameter :: hgt_year = 2021
    character*5, parameter :: hgt_ver = '1.1.3'
#endif

    esub_p => esub

    call clear_file_info(overmeta)

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', &
        100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO = 41
        120, 10, &     ! # files to write >= N_LAIMAX + 3*CHECKSUMS
        outputs_dir=THIS_OUTPUTS_DIR)

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')

    ! laimax
    call chunker%nc_open_set(ent20, io_laiin(:,1), &
        LAI_SOURCE, 'M', 'laimax', LAI_YEAR, 'ent17', '1.1')

    ! Bare Soil Brightness Ratio
    call chunker%nc_open(io_bs, chunker%outputs_dir, 'soilalbedo/', &
        'soilalbedo_1km_bs_brightratio_fill.nc', 'bs_brightratio', 1)

    call chunker%nc_open_set(ent20, io_simin(:,1), &
#if (defined HGT_GEDI) || (defined HGT_POTAPOV)
        LAI_SOURCE, 'Ha', 'hgt', hgt_year, 'ent17', hgt_ver &
#else
    ! Simard heights
        LAI_SOURCE, 'M', 'hgt', LAI_YEAR, 'ent17', '1.1' &
#endif
        )

    ! Climate statistics (we want TCinave = temperature [C])
    call chunker%nc_open_input(io_TCinave, &
        INPUTS_URL, INPUTS_DIR, &
        'climstats/CRU-TS3.22_GPCC-V6/', &
        'TCinave.nc', 'TCinave', 1)

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! PURE WITH WATER_ICE  

overmeta%global_data_source = "lc : Moderate Resolution Imaging Spectroradiometer (MODIS) "// &
  "MCD12Q1 L3 V051,, Land Cover, 500 m, annual (Friedl et "// &
  "al. 2010, doi:10.1016/j.rse.2009.08.016)"
    ! LC_pure
    call chunker%nc_create_set( &
        esub_p, io_lcout(:,1), &
        repeat_weights(esub%ncover, chunker%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1', overmeta=overmeta)

#if (defined HGT_GEDI)
overmeta%global_data_source = "hgt: P. Potapov et al. (2020) Mapping and monitoring global "// &
  "forest canopy height through integration of GEDI and Landsat data. Remote Sensing of Environment,"// &
  " 112165. https://doi.org/10.1016/j.rse.2020.112165"
overmeta%data_source = "hgt: GEDI heights (Landsat) (P. Potapov et al. 2020, https://doi.org/10.1016/j.rse.2020.112165)"
#elif (defined HGT_POTAPOV)
overmeta%global_data_source = "Potapov et  al. (2021) Remote Sensing of "// &
   "Environment, Volume 253. https://doi.org/10.1016/j.rse.2020.112165. Landsat/GEDI "// &
   "30 m global forest heights upscaled to 1 km mean and standard deviation. "// &
   "Personal communication, Peter Potapov, potapov@umd.edu."
overmeta%data_source = "hgt: GEDI heights (Landsat) (P. Potapov et al. 2020, https://doi.org/10.1016/j.rse.2020.112165)"
#else
overmeta%global_data_source = "hgt:  RH100 heights (Simard et al. 2011, doi:10.1029/2011jg001708)"
#endif
    ! ENTPFT heights in ENT16 indices
    call chunker%nc_create_set( &
        esub_p, io_simout(:,1), lc_weights(io_lcout(:,1), 1d0, 0d0), &
#if (defined HGT_GEDI) || (defined HGT_POTAPOV)
        LAI_SOURCE, 'Ha', 'hgt', hgt_year, 'pure', hgt_ver, &
#else
        LAI_SOURCE, 'M', 'hgt', LAI_YEAR, 'pure', '1.1', &
#endif
        overmeta=overmeta)

overmeta%global_data_source = "lai and laimax: Beijing Normal University LAI data product, "// &
  "1 km (Yuan et al. 2011, doi:10.1016/j.rse.2011.01.001)."
    ! laimax_pure
    call chunker%nc_create_set( &
        esub_p, io_laiout(:,1), lc_weights(io_lcout(:,1), 1d0, 0d0), &
        LAI_SOURCE, 'M', 'laimax', LAI_YEAR, 'pure', '1.1', overmeta=overmeta)

    ! PURE_NOH2O scale out water_ice

    !


    ! ------------- Checksums
    !   
    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lc_checksum(1), &
        weighting(chunker%wta1, 1d0, 0d0), &   ! Sum all cover
        info%dir, info%leaf, info%vname, &
        'SUM(lc)', info%units)

    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lchgt', LAI_YEAR, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lchgt_checksum(1), &
        weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(lc*height)', info%units)

    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lclaimax', LAI_YEAR, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lclai_checksum(1), &
        weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(LC*LAI)', info%units)

    call chunker%nc_check('B11_reclass_annual')
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call cropmerge_laisparse_splitbare(esub, chunker, one, &
#ifdef ENTGVSD_DEBUG
        dbj0,dbj1, &
        dbi0,dbi1, &
#else
        1,chunker%nchunk(2), &
        1,chunker%nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc, io_laiin, io_bs, io_TCinave, &
        io_laiout, &
        io_lclai_checksum=io_lclai_checksum, &
        io_lc_checksum=io_lc_checksum, &
        io_lcout=io_lcout, &
        io_simin=io_simin, &
        io_simout=io_simout, &
        io_lchgt_checksum=io_lchgt_checksum)

    call chunker%close_chunks

end subroutine do_reindex
end module b11_mod

! ====================================================================

program convert
    use b11_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub
    type(GcmEntSet_t), target :: esubnoh2o

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    call do_reindex(esub)


end program convert
