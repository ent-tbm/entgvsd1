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
    use trim_tiny_mod

implicit none

    integer, parameter :: one = 1

CONTAINS

subroutine do_reindex(esub)
    type(GcmEntSet_t), intent(IN) :: esub

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: ioall_lc, io_lc(NENT20)
    type(ChunkIO_t) :: ioall_laiin(one), io_laiin(NENT20,one)
    type(ChunkIO_t) :: io_bs
    type(ChunkIO_t) :: ioall_simin(one), io_simin(NENT20,one)
    ! Output files
    type(ChunkIO_t) :: ioall_lcout(one), io_lcout(esub%ncover,one)
    type(ChunkIO_t) :: ioall_laiout(one), io_laiout(esub%ncover,one)
    type(ChunkIO_t) :: ioall_sum, io_sum_lc(one)
    type(ChunkIO_t) :: ioall_simout(one), io_simout(esub%ncover,one)
    integer :: k,ksub

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', &
        100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO = 41
        120)     ! # files to write >= N_LAIMAX + 3*CHECKSUMS

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open(ioall_lc, LC_LAI_ENT_DIR, &
        'pure/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
    do k = 1,NENT20
        call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,k/))
    enddo

    ! laimax
    call chunker%nc_open(ioall_laiin(1), LC_LAI_ENT_DIR, &
        'pure/annual/', 'entmm29_ann_laimax.nc', 'lai', 0)
    do k = 1,NENT20
        call chunker%nc_reuse_var(ioall_laiin(1), io_laiin(k,1), (/1,1,k/))
    enddo

    ! Bare Soil Brightness Ratio
    call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
        'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)

    ! Simard heights
    call chunker%nc_open(ioall_simin(1), lC_LAI_ENT_DIR, &
        'pure/annual/', 'entmm29_ann_height.nc', 'SimardHeights', 0)
    do k=1,NENT20
        call chunker%nc_reuse_var(ioall_simin(1), io_simin(k,1), &
            (/1,1,k/))
    end do

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! LC_pure
    call chunker%nc_create(ioall_lcout(1), &
        weighting(chunker%wta1, 1d0, 0d0), &
        'pure2/annual/', 'entmm29_ann_lc', 'lc', &
        'Ent Landcover (from A04)', '1', 'Land Cover', &
        esub%mvs, esub%layer_names())
    do k=1,esub%ncover
        call chunker%nc_reuse_var(ioall_lcout(1), io_lcout(k,1), &
            (/1,1,k/), weighting(io_lc(k)%buf, 1d0,0d0))
    enddo


    ! ENTPFT heights in ENT16 indices
    call chunker%nc_create(ioall_simout(1), &
        weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
        'pure2/annual/', 'entmm29_ann_height', &
        'SimardHeights', &
        'Plant Heights', 'm', 'Plant Heights', &
        esub%mvs, esub%layer_names())
    do k = 1,esub%ncover
        call chunker%nc_reuse_var(ioall_simout(1), io_simout(k,1), &
            (/1,1,k/), weighting(io_lc(k)%buf,1d0,0d0))
    end do

    ! laimax_pure
    call chunker%nc_create(ioall_laiout(1), &
        weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
        'pure2/annual/', 'entmm29_ann_laimax', 'lai', &
        'Ent maximum LAI for year', 'm^2 m-2', 'Leaf Area Index', &
        esub%mvs, esub%layer_names())
    do k=1,esub%ncover
        call chunker%nc_reuse_var(ioall_laiout(1), io_laiout(k,1), &
            (/1,1,k/), weighting(io_lcout(k,1)%buf, 1d0,0d0))
    enddo

    !  checksum land  laimax
    call chunker%nc_create(ioall_sum, &
        weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
        '16/nc/', 'V1km_EntGVSDv1.1_LAI3g16_laimax_pure_checksum')

    call chunker%nc_reuse_file(ioall_sum, io_sum_lc(1), &
        'lc_checksum', 'Checksum of LC', '1', 'checksum - Land Cover', &
        weighting(chunker%wta1, 1d0, 0d0))

    call chunker%nc_check('A04_trim_laimax_1kmx1km')
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call trim_tiny(esub, chunker, one, &

#ifdef ENTGVSD_DEBUG
        nchunk(2)*3/4,nchunk(2)*3/4+1, &
        nchunk(1)*3/4,nchunk(1)*3/4+1, &
#else
        1,nchunk(2), &
        1,nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc, io_laiin, io_bs, &
        io_laiout, &
        io_sum_lc=io_sum_lc, &
        io_lcout=io_lcout, &
        io_simin=io_simin, &
        io_simout=io_simout)

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
    type(GcmEntSet_t) :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_reindex(esub)


end program convert
