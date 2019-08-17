program A00b_regrid
    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
    use hntr_mod

implicit none

    type(Chunker_t) :: chunker,chunkerlr
    type(FileInfo_t) :: info
    type(ChunkIO_t) :: ioall_lc, io_lc(2)
    type(ChunkIO_t) :: ioall_lcout, io_lcout(2)
    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid
    type(EntSet_t) :: ent2      ! Master Ent categories
    integer :: ichunk,jchunk, k,ic,jc

    call init_ent_labels
    call ent2%allocate(2,NENT20)
    call ent2%sub_covertype(ent20, SNOW_ICE)
    call ent2%sub_covertype(ent20, CV_WATER)

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 10, 10, 10)
    call chunkerlr%init(IMK,JMK, IMH*2,JMH*2, 'qxq', 10, 10, 10)

    ! Hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1), chunker%ngrid(2), 0d0, 180d0*60d0 / chunker%ngrid(2))
    spec_lr = hntr_spec(chunkerlr%chunk_size(1), chunkerlr%ngrid(2), 0d0, 180d0*60d0 / chunkerlr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, 0d0)   ! datmis=0


    ! ---------------- Open land cover
    call chunker%file_info(info, ent20, 'BNU', 'M', 'lc', 2004, 'ent17', '1.1')
    call chunker%nc_open(ioall_lc, &
        LC_LAI_ENT_DIR, trim(info%dir), trim(info%leaf)//'.nc.test', trim(info%vname), 0)
    do k=1,2
        call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,ent2%mvs(k)/))
    end do


    ! --------------- Open output file
    call chunkerlr%file_info(info, ent20, 'BNU', 'M', 'lc', 2004, 'ent17', '1.1')
    call chunkerlr%nc_create(ioall_lcout, &
        weighting(chunkerlr%wta1,1d0,0d0), &
        trim(info%dir), trim(info%leaf), trim(info%vname), &
        'Land Cover Fractions', '1', &
        ent2%layer_names(), create_lr=.false.)
    do k=1,2
        call chunkerlr%nc_reuse_var( &
            ioall_lcout, io_lcout(k), (/1,1,k/), &
            weighting(chunkerlr%wta1, 1d0, 0d0))
    end do


    ! -------------- Regrid!
#ifdef ENTGVSD_DEBUG
!    do jchunk = chunker%nchunk(2)*3/4,chunker%nchunk(2)*3/4+1
!    do ichunk = chunker%nchunk(1)*3/4,chunker%nchunk(1)*3/4+1
    ! Choose a smalle area to debug
    do jchunk = 11,12
    do ichunk = 5,7
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif

        call chunker%move_to(ichunk,jchunk)
        call chunkerlr%move_to(ichunk,jchunk)

!        do k=1,2
!        do jc = 1,chunker%chunk_size(2)
!        do ic = 1,chunker%chunk_size(1)
!            if (io_lc(k)%buf(ic,jc) == FillValue) io_lc(k)%buf(ic,jc) = 0
!        end do
!        end do
!        end do        

#if 1
        do k=1,2
!            print *,io_lcout(k)%startB(2), io_lcout(k)%chunker%chunk_size(2)
            print *,io_lc(k)%buf(1,1),sum(io_lc(k)%buf)
            call hntr_lr%regrid4( &
                io_lcout(k)%buf, io_lc(k)%buf, &
                chunker%wta1, 1d0, 0d0, &   ! weighting
                io_lcout(k)%startB(2), io_lcout(k)%chunker%chunk_size(2))
        end do
#else
        do k=1,2
            io_lcout(k)%buf(:,:) = 0
        end do
#endif

        !call chunker%write_chunks
        call chunkerlr%write_chunks

    end do
    end do

end program A00b_regrid
