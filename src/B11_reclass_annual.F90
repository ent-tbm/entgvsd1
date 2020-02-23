! Reclass of annual LAIMAX
!
! Takes 20-cover classes (Ent 17 PFTs + 3 non-veg) and converts to
! 18-cover classes (Ent 16 PFTs and bright and dark bare soil fractions).
! Merges C3 and C4 crops into one crop cover type for Ent 16 PFTs, excludes water
! and permanent ice, and converts barse/sparse cover into equivalent veg type and
! bare soil bright and dark fractions, preserving total LAI of grid cell.
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
    type(FileInfo_t) :: info
    integer :: k,ksub

    esub_p => esub

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
        'soilalbedo_1km_bs_brightratio.nc', 'bs_brightratio', 1)

    ! Simard heights
    call chunker%nc_open_set(ent20, io_simin(:,1), &
        LAI_SOURCE, 'M', 'hgt', LAI_YEAR, 'ent17', '1.1')

    ! Climate statistics (we want TCinave = temperature [C])
    call chunker%nc_open_input(io_TCinave, &
        INPUTS_URL, INPUTS_DIR, &
        'climstats/CRU-TS3.22_GPCC-V6/', &
        'TCinave.nc', 'TCinave', 1)

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! LC_pure
    call chunker%nc_create_set( &
        esub_p, io_lcout(:,1), &
        repeat_weights(esub%ncover, chunker%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1')

    ! ENTPFT heights in ENT16 indices
    call chunker%nc_create_set( &
        esub_p, io_simout(:,1), lc_weights(io_lcout(:,1), 1d0, 0d0), &
        LAI_SOURCE, 'M', 'hgt', LAI_YEAR, 'pure', '1.1')

    ! laimax_pure
    call chunker%nc_create_set( &
        esub_p, io_laiout(:,1), lc_weights(io_lcout(:,1), 1d0, 0d0), &
        LAI_SOURCE, 'M', 'laimax', LAI_YEAR, 'pure', '1.1')

    ! ------------- Checksums
    !  checksum land  laimax
    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lc_checksum(1), &
        weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(lc)', info%units)

    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lchgt', LAI_YEAR, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lchgt_checksum(1), &
        weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! TODO: Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(lc*height)', info%units)

    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lclaimax', LAI_YEAR, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lclai_checksum(1), &
        weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! TODO: Scale by _lc
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

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    call do_reindex(esub)


end program convert
