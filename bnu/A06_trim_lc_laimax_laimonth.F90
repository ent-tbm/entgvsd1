module a06_mod

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
    use trim_tiny_mod

implicit none

CONTAINS

subroutine do_reindex(esub,m0,m1)
    type(GcmEntSet_t), intent(IN) :: esub
    integer :: m0,m1    ! First and last month to work on

    type(Chunker_t) :: chunker

    ! ------ Input Files
    type(ChunkIO_t) :: ioall_lc, io_lc(NENT20)
    type(ChunkIO_t) :: ioall_lc2, io_lc2(esub%ncover)
    type(ChunkIO_t) :: ioall_laiin(m1-m0+1), io_laiin(NENT20,m1-m0+1)
    type(ChunkIO_t) :: io_bs
    ! ------ Output files
    type(ChunkIO_t) :: ioall_laiout(m1-m0+1), io_laiout(esub%ncover,m1-m0+1)


    integer :: k,ksub
    integer :: im,imonth


    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 150, 150)

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open(ioall_lc, LC_LAI_ENT_DIR, &
        'pure/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
    do k = 1,NENT20
        call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,k/))
    enddo

    ! LC written by A04; in the esub indexing scheme
    call chunker%nc_open(ioall_lc2, LC_LAI_ENT_DIR, &
        'pure2/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
    do k = 1,esub%ncover
        call chunker%nc_reuse_var(ioall_lc2, io_lc2(k), (/1,1,k/))
    enddo

    ! laiin
    do im = m0,m1
        imonth = im - m0 + 1

        call chunker%nc_open(ioall_laiin(imonth), LC_LAI_ENT_DIR, &
            'pure/monthly/', 'entmm29_'//trim(MONTH(im))//'_lai.nc', 'lai', 0)
        do k = 1,NENT20
            call chunker%nc_reuse_var(ioall_laiin(imonth), io_laiin(k,imonth), (/1,1,k/))
        enddo
    end do

    ! bs ratio
    call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
        'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)
    ! For now at least, use BS BrightRatio from LAI3G dataset
    ! call chunker%nc_open(io_bs, LC_LAI_ENT_DIR, &
    !     '', 'bs_brightratio.nc', 'bs_brightratio', 1)

    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    do im = m0,m1
        imonth = im - m0 + 1
        call chunker%nc_create(ioall_laiout(imonth), &
            weighting(chunker%wta1,1d0,0d0), &
            'pure2/monthly/', 'entmm29_'//MONTH(im)//'_lai', 'lai', &
            'Ent LAI on the given day', 'm^2 m-2', 'Leaf Area Index', &
            esub%mvs, esub%layer_names())

        do k=1,esub%ncover
            call chunker%nc_reuse_var(ioall_laiout(imonth), io_laiout(k,imonth), &
                (/1,1,k/), weighting(io_lc2(k)%buf, 1d0,0d0))
        end do   ! k
    end do   ! imonth

    call chunker%nc_check('A06_trim_lc_laimax_laimonth')
    print *,'Done opening files: nreads',chunker%nreads,'nwrite',chunker%nwrites
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call trim_tiny(esub, chunker, m1-m0+1, &

#ifdef ENTGVSD_DEBUG
        nchunk(2)*3/4,nchunk(2)*3/4+1, &
        nchunk(1)*3/4,nchunk(1)*3/4+1, &
#else
        1,nchunk(2), &
        1,nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc, io_laiin, io_bs, &
        io_laiout)

    call chunker%close_chunks

end subroutine do_reindex
end module a06_mod
! ====================================================================

program convert
    use a06_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t) :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_reindex(esub,1,6)
    call do_reindex(esub,7,12)


end program convert
