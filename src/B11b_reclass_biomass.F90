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

module b11b_mod
    use netcdf
    use chunker_mod
    use ent_labels_mod
    use ent_params_mod
    use gcm_labels_mod
    use cropmerge_laisparse_splitbare_mod

implicit none

    integer, parameter :: one = 1
#if (defined BIOMASS_SPAWN)
    integer, parameter :: n_biomass = 2
#elif (defined BIOMASS_GEDI) || (BIOMASS_XU)
    integer, parameter :: n_biomass = 1
#else
    integer, parameter :: n_biomass = -1
#endif

CONTAINS

subroutine do_reindex(esub)
    type(GcmEntSet_t), intent(IN), target :: esub

    class(EntSet_t), pointer :: esub_p

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: io_lc(NENT20)
    type(ChunkIO_t) :: io_biomassin(NENT20,n_biomass)
    type(ChunkIO_t) :: io_bs
    !type(ChunkIO_t) :: io_simin(NENT20,one)
    type(ChunkIO_t) :: io_TCinave
    ! Output files
    type(ChunkIO_t) :: io_lcout(esub%ncover,one)
    type(ChunkIO_t) :: io_biomassout(esub%ncover,n_biomass)
    type(ChunkIO_t) :: io_lc_checksum(one)
    type(ChunkIO_t) :: io_lchgt_checksum(one)
    type(ChunkIO_t) :: io_lclai_checksum(one)
    type(ChunkIO_t) :: io_simout(esub%ncover,one)
    type(FileInfo_t) :: info, overmeta
    integer :: k,ksub

    esub_p => esub

    call clear_file_info(overmeta)
#if (defined BIOMASS_SPAWN)
    overmeta%data_source = 'https://doi.org/10.1038/s41597-020-0444-4'
    overmeta%global_data_source = &
        'biomass: Center for Sustainability and the Global Environment, ' // &
            'Nelson Institute for Environmental Studies, University ' // &
            'of Wisconsin-Madison, Biomass Carbon Density, 300m, ' // &
            '(Spawn et al. 2010, doi:10.1038/s41597-020-0444-4)'
#elif (defined BIOMASS_GEDI)
    overmeta%data_source = 'https://doi.org/10.3334/ORNLDAAC/2017'
    overmeta%global_data_source = &
        'biomass: Dubayah, R.O., J. Armston, S.P. Healey, Z. Yang, P.L. ' // &
            'Patterson, S. Saarela, G. Stahl, L. Duncanson, and J.R. Kellner. 2022.' // &
            'GEDI L4B Gridded Aboveground Biomass Density, Version 2. ORNL DAAC, Oak ' // &
            'Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/2017. NOTE: DRYBIOMASS'
#elif (defined BIOMASS_XU)
    overmeta%data_source = 'https://doi.org/10.1126/sciadv.abe9829'
    overmeta%global_data_source = &
        'biomass: Liang Xu et al. ,Changes in global terrestrial live biomass ' // &
            'over the 21st century.Sci. Adv.7,eabe9829(2021).DOI:10.1126/sciadv.abe9829'
#endif

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'QXQ', &
        100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO = 41
        120, 10, &     ! # files to write >= N_LAIMAX + 3*CHECKSUMS
        outputs_dir=THIS_OUTPUTS_DIR)

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')

    ! biomass
#if (defined BIOMASS_SPAWN)
    call chunker%nc_open_set(ent20, io_biomassin(:,1), &
        'Spawn2020', '', 'biomass_agb', 2010, 'ent17', '1.1.2')

    call chunker%nc_open_set(ent20, io_biomassin(:,2), &
        'Spawn2020', '', 'biomass_bgb', 2010, 'ent17', '1.1.2')
#elif (defined BIOMASS_GEDI)
    call chunker%nc_open_set(ent20, io_biomassin(:,1), &
        'GEDI', '', 'biomass_agb', 2022, 'ent17', '1.1.2', &
        heightsource='H2')
#elif (defined BIOMASS_XU)
    call chunker%nc_open_set(ent20, io_biomassin(:,1), &
        'Xu', '', 'biomass_agb', 2004, 'ent17', '1.1.2', &
        heightsource='X')
#endif

    ! Bare Soil Brightness Ratio
    call chunker%nc_open(io_bs, chunker%outputs_dir, 'soilalbedo/', &
        'soilalbedo_1km_bs_brightratio_fill.nc', 'bs_brightratio', 1)

    ! Climate statistics (we want TCinave = temperature [C])
    call chunker%nc_open_input(io_TCinave, &
        INPUTS_URL, INPUTS_DIR, &
        'climstats/CRU-TS3.22_GPCC-V6/', &
        'TCinave.nc', 'TCinave', 1)

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! PURE WITH WATER_ICE  

    ! biomass_pure
#if (defined BIOMASS_SPAWN)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        esub_p, io_biomassout(:,1), lc_weights(io_lc, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_agb', 2010, 'pure', '1.1.2', &
        overmeta=overmeta)

    call chunker%nc_create_set( & ! lcweights are dummy!!
        esub_p, io_biomassout(:,2), lc_weights(io_lc, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_bgb', 2010, 'pure', '1.1.2', &
        overmeta=overmeta)
#elif (defined BIOMASS_GEDI)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        esub_p, io_biomassout(:,1), lc_weights(io_lc, 0d0, 1d0), &
        'GEDI', '', 'biomass_agb', 2022, 'pure', '1.1.2', &
        overmeta=overmeta, heightsource='H2')
#elif (defined BIOMASS_XU)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        esub_p, io_biomassout(:,1), lc_weights(io_lc, 0d0, 1d0), &
        'Xu', '', 'biomass_agb', 2004, 'pure', '1.1.2', &
        overmeta=overmeta, heightsource='X')
#endif

    ! PURE_NOH2O scale out water_ice

    !


    ! ------------- Checksums
    !   
!   call chunker%file_info(info, esub_p, &
!       LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1', &
!       varsuffix = '_checksum')
!   call chunker%nc_create(io_lc_checksum(1), &
!       weighting(chunker%wta1, 1d0, 0d0), &   ! Sum all cover
!       info%dir, info%leaf, info%vname, &
!       'SUM(lc)', info%units)

!   call chunker%file_info(info, esub_p, &
!       LAI_SOURCE, 'M', 'lchgt', LAI_YEAR, 'pure', '1.1', &
!       varsuffix = '_checksum')
!   call chunker%nc_create(io_lchgt_checksum(1), &
!       weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! Scale by _lc
!       info%dir, info%leaf, info%vname, &
!       'SUM(lc*height)', info%units)

!   call chunker%file_info(info, esub_p, &
!       LAI_SOURCE, 'M', 'lclaimax', LAI_YEAR, 'pure', '1.1', &
!       varsuffix = '_checksum')
!   call chunker%nc_create(io_lclai_checksum(1), &
!       weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! Scale by _lc
!       info%dir, info%leaf, info%vname, &
!       'SUM(LC*LAI)', info%units)

    call chunker%nc_check('B11b_reclass_biomass')
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call cropmerge_laisparse_splitbare(esub, chunker, n_biomass, &
#ifdef ENTGVSD_DEBUG
        dbj0,dbj1, &
        dbi0,dbi1, &
#else
        1,chunker%nchunk(2), &
        1,chunker%nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc, io_biomassin, io_bs, io_TCinave, &
        io_biomassout)
        !io_lclai_checksum=io_lclai_checksum, &
        !io_lc_checksum=io_lc_checksum, &
        !io_lcout=io_lcout, &
        !io_simin=io_simin, &
        !io_simout=io_simout, &
        !io_lchgt_checksum=io_lchgt_checksum)

    call chunker%close_chunks

end subroutine do_reindex
end module b11b_mod

! ====================================================================

program convert
    use b11b_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub
    type(GcmEntSet_t), target :: esubnoh2o

#if (defined BIOMASS_SPAWN) || (defined BIOMASS_GEDI) ||\
    (defined BIOMASS_XU)
#else
    write(*,*) "Missing argument for biomass: -b SPAWN, GEDI, or XU"
    stop 1
#endif

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    call do_reindex(esub)


end program convert
