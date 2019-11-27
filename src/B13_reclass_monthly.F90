#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module b13_mod

    use netcdf
    use chunker_mod
    use ent_labels_mod
    use gcm_labels_mod
    use ent_params_mod
    use cropmerge_laisparse_splitbare_mod

implicit none

CONTAINS

subroutine do_reindex(esub,m0,m1)
    type(GcmEntSet_t), intent(IN), target :: esub
    integer :: m0,m1    ! First and last month to work on

    class(EntSet_t), pointer :: esub_p
    type(Chunker_t) :: chunker

    ! ------ Input Files
    type(ChunkIO_t) :: io_lc_ent17(NENT20)
    type(ChunkIO_t) :: io_lc_pure(esub%ncover)
    real*4, ALLOCATABLE :: sum_lc(:,:)
    type(ChunkIO_t) :: io_laiin(NENT20,m1-m0+1)
    type(ChunkIO_t) :: io_bs
    ! ------ Output files
    type(ChunkIO_t) :: io_laiout(esub%ncover,m1-m0+1)
    type(ChunkIO_t) :: io_lclai_checksum(m1-m0+1)


    type(FileInfo_t) :: info
    integer :: k,ksub
    integer :: im,imonth


    esub_p => esub
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 300, 300, 20, outputs_dir=THIS_OUTPUTS_DIR)
    allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc_ent17, &
        LAI_SOURCE, 'M', 'lc', 2004, 'ent17', '1.1')

    ! LC written by A04; in the esub indexing scheme
    call chunker%nc_open_set(esub_p, io_lc_pure, &
        LAI_SOURCE, 'M', 'lc', 2004, 'pure', '1.1')

    ! laiin
    do im = m0,m1
        imonth = im - m0 + 1

        call chunker%nc_open_set(ent20, io_laiin(:,imonth), &
            LAI_SOURCE, 'M', 'lai', 2004, 'ent17', '1.1', &
            doytype='month', idoy=im)
    end do

    ! bs ratio
    call chunker%nc_open(io_bs, chunker%outputs_dir, 'soilalbedo/', &
        'V1km_bs_brightratio.nc', 'bs_brightratio', 1)

    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    do im = m0,m1
        imonth = im - m0 + 1

        call chunker%nc_create_set( &
            esub_p, io_laiout(:,imonth), lc_weights(io_lc_pure, 1d0, 0d0), &
            LAI_SOURCE, 'M', 'lai', 2004, 'pure', '1.1', &
            doytype='month', idoy=im)

        call chunker%file_info(info, esub_p, &
            LAI_SOURCE, 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='month', idoy=im)
        call chunker%nc_create(io_lclai_checksum(imonth), &
            weighting(sum_lc, 1d0, 0d0), &
            info%dir, info%leaf, info%vname, &
            'SUM(LC*LAI)', info%units)
    end do   ! imonth

    call chunker%nc_check('B13_reclass_monthly')
    print *,'Done opening files: nreads',chunker%nreads,'nwrite',chunker%nwrites
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call cropmerge_laisparse_splitbare(esub, chunker, m1-m0+1, &

#ifdef ENTGVSD_DEBUG
        dbj0,dbj1, &
        dbi0,dbi1, &
#else
        1,chunker%nchunk(2), &
        1,chunker%nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc_ent17, io_laiin, io_bs, &
        io_laiout, &
        sum_lc = sum_lc, &
        io_lclai_checksum=io_lclai_checksum)

    call chunker%close_chunks

end subroutine do_reindex
end module b13_mod
! ====================================================================

program convert
    use b13_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub
    integer, parameter :: monthchunk = 6
    integer :: m0,m1

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    do m0=1,nmonth,monthchunk
        m1 = min(m0+monthchunk-1, nmonth)

        call do_reindex(esub,m0,m1)
    end do


end program convert
