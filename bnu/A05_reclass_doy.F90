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
    type(GcmEntSet_t), intent(IN) :: esub

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: ioall_lc, io_lc(NENT20)
    type(ChunkIO_t) :: ioall_lc2, io_lc2(esub%ncover)
    type(ChunkIO_t) :: ioall_laiin(ndoy), io_laiin(NENT20,ndoy)
    type(ChunkIO_t) :: io_bs
    ! Output files
    type(ChunkIO_t) :: ioall_laiout(ndoy), io_laiout(esub%ncover,ndoy)

    integer :: k,ksub
    integer :: idoy

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

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


    ! laimax
    do idoy = 1,ndoy
        call chunker%nc_open(ioall_laiin(idoy), LC_LAI_ENT_DIR, &
            'pure/doy/', 'entmm29_'//trim(DOY(idoy))//'_lai.nc', 'lai', 0)
        do k = 1,NENT20
            call chunker%nc_reuse_var(ioall_laiin(idoy), io_laiin(k,idoy), (/1,1,k/))
        enddo
    end do

    ! Bare Soil Brightness Ratio
    call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
        'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)

    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    do idoy = 1,ndoy
        call chunker%nc_create(ioall_laiout(idoy), &
            weighting(chunker%wta1,1d0,0d0), &
            'pure2/doy/', 'entmm29_'//DOY(idoy)//'_lai', 'lai', &
            'Ent LAI on the given day', 'm^2 m-2', 'Leaf Area Index', &
            esub%mvs, esub%layer_names())
        do k=1,esub%ncover
            call chunker%nc_reuse_var(ioall_laiout(idoy), io_laiout(k,idoy), &
                (/1,1,k/), weighting(io_lc2(k)%buf, 1d0,0d0))
        end do   ! k
    end do   ! idoy


    call chunker%nc_check('A05_trim_laidoy_1kmx1km')
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
        io_lc, io_laiin, io_bs, &
        io_laiout)

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
