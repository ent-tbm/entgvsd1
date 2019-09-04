module a05_mod

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

subroutine do_reindex(esub)
    type(GcmEntSet_t), intent(IN), target :: esub
    class(EntSet_t), pointer :: esub_p

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: io_lc_ent17(NENT20)
    type(ChunkIO_t) :: io_lc_pure(esub%ncover)
    real*4 :: sum_lc(:,:)
    type(ChunkIO_t) :: io_laiin(NENT20,ndoy)
    type(ChunkIO_t) :: io_bs
    ! Output files
    type(ChunkIO_t) :: io_laiout(esub%ncover,ndoy)
    type(ChunkIO_t) :: io_lclai_checksum(ndoy)

    integer :: k,ksub
    integer :: idoy

    esub_p => esub

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120, 10)
    allocate(sum_lc(chunker%chunk_size(0), chunker%chunk_size(1)))

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc_ent17, &
        'BNU', 'M', 'lc', 2004, 'ent17', '1.1')

    ! LC written by A04; in the esub indexing scheme
    call chunker%nc_open_set(esub_p, io_lc_pure, &
        'BNU', 'M', 'lc', 2004, 'pure', '1.1')

    ! lai
    do idoy = 1,ndoy
        call chunker%nc_open_set(ent20, io_laiin(:,idoy), &
            'BNU', 'M', 'lai', 2004, 'ent17', '1.1', &
            doytype='doy', idoy=idoy)
    end do

    ! Bare Soil Brightness Ratio
    call chunker%nc_open(io_bs, LC_LAI_ENT_DIR, 'carrer/', &
        'V1km_bs_brightratio.nc', 'bs_brightratio', 1)

    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    do idoy = 1,ndoy
        call chunker%nc_create_set( &
            esub_p, io_laiout(:,idoy), lc_weights(io_lc_pure, 1d0, 0d0), &
            'BNU', 'M', 'lai', 2004, 'pure', '1.1', &
            doytype='doy', idoy=idoy)

        call chunker%file_info(info, esub_p, &
            'BNU', 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='doy', idoy=idoy)
        call chunker%nc_create(io_lclai_checksum(idoy), &
            weighting(sum_lc, 1d0, 0d0), &   ! TODO: Scale by _lc
            info%dir, info%leaf, info%vname, &
            'SUM(LC*LAI)', info%units)
    end do   ! idoy

    call chunker%nc_check('A05_reclass_doy')
#ifdef JUST_DEPENDENCIES
    stop 0
#endif
    call cropmerge_laisparse_splitbare(esub, chunker, ndoy, &

#ifdef ENTGVSD_DEBUG
        chunker%nchunk(2)*3/4,chunker%nchunk(2)*3/4+1, &
        chunker%nchunk(1)*3/4,chunker%nchunk(1)*3/4+1, &
#else
        1,chunker%nchunk(2), &
        1,chunker%nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc_ent17, io_laiin, io_bs, &
        io_laiout,
        sum_lc=sum_lc,
        io_lclai_checksum=io_lclai_checksum)

    call chunker%close_chunks

end subroutine do_reindex
end module a05_mod

! ====================================================================

program convert
    use a05_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t) :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_reindex(esub)


end program convert
