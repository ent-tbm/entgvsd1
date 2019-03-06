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


    ! Combine C3 and C4 crops into one PFT, for ModelE
    logical, parameter :: combine_crops_c3_c4 = .true.
    ! Split the bare soil into dark and light, to get the right albedo
    logical, parameter :: split_bare_soil = .true.

    integer, parameter :: one = 1

CONTAINS

subroutine do_reindex(esub)
    type(GcmEntSet_t), intent(IN) :: esub

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: io_lcin(NENT20), io_laiin(NENT20,one), io_bs
    ! Output files
    type(ChunkIO_t) :: ioall_laiout(one), io_laiout(esub%ncover,one)
    type(ChunkIO_t) :: ioall_sum, io_sum_lc(one)

    integer :: k,ksub

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', &
        100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO = 41
        120)     ! # files to write >= N_LAIMAX + 3*CHECKSUMS

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! lcmax
    do k=1,NENT20
        ! USE the land cover we computed in a previous step!!!!!!
        ! TODO: LAI3g????
        ! PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
        call chunker%nc_open(io_lcin(k), LC_LAI_ENT_DIR, &
            'EntMM_lc_laimax_1kmx1km/', &
            trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lc.nc', &
            trim(ent20%abbrev(k)), 1)
    end do

    ! laimax
    do k=1,NENT20
        call chunker%nc_open(io_laiin(k,1), LC_LAI_ENT_DIR, &
            'EntMM_lc_laimax_1kmx1km/', &
            trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lai.nc', &
            trim(ent20%abbrev(k)), 1)
    enddo

    ! Bare Soil Brightness Ratio
    call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
        'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! laimax_pure
    call chunker%nc_create(ioall_laiout(1), &
        weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
        '16/nc/', 'V1km_EntGVSDv1.1_BNU16_laimax_pure')
    do ksub=1,esub%ncover
        call chunker%nc_reuse_file(ioall_laiout(1), io_laiout(ksub,1), &
            'lai_'//trim(esub%abbrev(ksub)), trim(esub%title(ksub)), &
            'm2 m-2', trim(esub%title(ksub)), &
            weighting(chunker%wta1,1d0,0d0))    ! TODO: Weighting???
    end do

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
        io_lcin, io_laiin, io_bs, &
        io_laiout, &
        io_sum_lc=io_sum_lc)

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
