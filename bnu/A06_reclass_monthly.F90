module a06_mod

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
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
    type(ChunkIO_t) :: io_laiin(NENT20,m1-m0+1)
    type(ChunkIO_t) :: io_bs
    ! ------ Output files
    type(ChunkIO_t) :: io_laiout(esub%ncover,m1-m0+1)


    integer :: k,ksub
    integer :: im,imonth


    esub_p => esub
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 300, 300, 20)

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc_ent17, &
        'BNU', 'M', 'lc', 2004, 'ent17', '1.1')

    ! LC written by A04; in the esub indexing scheme
    call chunker%nc_open_set(esub_p, io_lc_pure, &
        'BNU', 'M', 'lc', 2004, 'pure', '1.1')

    ! laiin
    do im = m0,m1
        imonth = im - m0 + 1

        call chunker%nc_open_set(ent20, io_laiin(:,imonth), &
            'BNU', 'M', 'lai', 2004, 'ent17', '1.1', &
            doytype='month', idoy=im)
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

        call chunker%nc_create_set( &
            esub_p, io_laiout(:,imonth), lc_weights(io_lc_pure, 1d0, 0d0), &
            'BNU', 'M', 'lai', 2004, 'pure', '1.1', &
            doytype='month', idoy=im)
    end do   ! imonth

    call chunker%nc_check('A06_reclass_monthly')
    print *,'Done opening files: nreads',chunker%nreads,'nwrite',chunker%nwrites
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call cropmerge_laisparse_splitbare(esub, chunker, m1-m0+1, &

#ifdef ENTGVSD_DEBUG
        chunker%nchunk(2)*3/4,chunker%nchunk(2)*3/4+1, &
        chunker%nchunk(1)*3/4,chunker%nchunk(1)*3/4+1, &
#else
        1,chunker%nchunk(2), &
        1,chunker%nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc_ent17, io_laiin, io_bs, &
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
    type(GcmEntSet_t), target :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    if (nmonth == 1) then
        call do_reindex(esub,1,1)
    else
        call do_reindex(esub,1,6)
        call do_reindex(esub,7,12)
    end if


end program convert
