! Trims off tiny fractions and preserves total LAI for the gridcell.
! Alters the cover amount if the LAI is smaller than the main one.
! If LAI is a little biger, might increase.
!
! I don't see where the program is doing that.

! A04 only creates the pure dataset.  It does NOT trim.

module a04_mod
    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
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

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', &
        100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO = 41
        120, 10)     ! # files to write >= N_LAIMAX + 3*CHECKSUMS

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc, &
        'BNU', 'M', 'lc', 2004, 'ent17', '1.1')

    ! laimax
    call chunker%nc_open_set(ent20, io_laiin(:,1), &
        'BNU', 'M', 'laimax', 2004, 'ent17', '1.1')

    ! Bare Soil Brightness Ratio
    call chunker%nc_open(io_bs, LC_LAI_ENT_DIR, 'carrer/', &
        'V1km_bs_brightratio.nc', 'bs_brightratio', 1)

    ! Simard heights
    call chunker%nc_open_set(ent20, io_simin(:,1), &
        'BNU', 'M', 'hgt', 2004, 'ent17', '1.1')

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! LC_pure
    call chunker%nc_create_set( &
        esub_p, io_lcout(:,1), &
        repeat_weights(esub%ncover, chunker%wta1, 1d0, 0d0), &
        'BNU', 'M', 'lc', 2004, 'pure', '1.1')

    ! ENTPFT heights in ENT16 indices
    call chunker%nc_create_set( &
        esub_p, io_simout(:,1), lc_weights(io_lcout(:,1), 1d0, 0d0), &
        'BNU', 'M', 'hgt', 2004, 'pure', '1.1')

    ! laimax_pure
    call chunker%nc_create_set( &
        esub_p, io_laiout(:,1), lc_weights(io_lcout(:,1), 1d0, 0d0), &
        'BNU', 'M', 'laimax', 2004, 'pure', '1.1')

    ! ------------- Checksums
    !  checksum land  laimax
    call chunker%file_info(info, esub_p, &
        'BNU', 'M', 'lc', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lc_checksum(1), &
        weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(lc)', info%units)

    call chunker%file_info(info, esub_p, &
        'BNU', 'M', 'lchgt', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lchgt_checksum(1), &
        weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! TODO: Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(lc*height)', info%units)

    call chunker%file_info(info, esub_p, &
        'BNU', 'M', 'lclai', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    call chunker%nc_create(io_lclai_checksum(1), &
        weighting(io_lc_checksum(1)%buf, 1d0, 0d0), &   ! TODO: Scale by _lc
        info%dir, info%leaf, info%vname, &
        'SUM(lc*LAI)', info%units)

    call chunker%nc_check('A04_reclass_annual')
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
        io_lc, io_laiin, io_bs, &
        io_laiout, &
        io_lclai_checksum=io_lclai_checksum, &
        io_lc_checksum=io_lc_checksum, &
        io_lcout=io_lcout, &
        io_simin=io_simin, &
        io_simout=io_simout, &
        io_lchgt_checksum=io_lchgt_checksum)

    call chunker%close_chunks

end subroutine do_reindex
end module a04_mod

! ====================================================================

program convert
    use a04_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    call do_reindex(esub)


end program convert
