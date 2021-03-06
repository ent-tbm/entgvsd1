! Regrids LC to 6km grid needed for Carrer albedo processing
! JUST regrids SNOW_ICE and CV_WATER.  Uses ent2 universe to do so.
! Author: Elizabeth Fischer

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

program B03_regrid_snowice
    use netcdf
    use chunker_mod
    use ent_params_mod
    use ent_labels_mod
    use gcm_labels_mod
    use hntr_mod

implicit none

    type(Chunker_t) :: chunker,chunkerlr
    type(FileInfo_t) :: info
    type(ChunkIO_t) :: ioall_lc, io_lc(2)
    type(ChunkIO_t) :: ioall_lcout, io_lcout(2) !file handle
    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid
    type(EntSet_t) :: ent2      ! Master Ent categories
    integer :: ichunk,jchunk, k,ic,jc


    type(ReadWrites_t) :: rw

    call rw%init(THIS_OUTPUTS_DIR, 'B03_regrid_snowice', 3,3)
    call init_ent_labels
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 10, 10, 10, outputs_dir=THIS_OUTPUTS_DIR)!forplot at 0.25 degrees.
    call chunkerlr%init(IMK,JMK, IMH*2,JMH*2, 'forplot', 10, 10, 10, outputs_dir=THIS_OUTPUTS_DIR)

    ! Hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1), chunker%ngrid(2), 0d0, 180d0*60d0 / chunker%ngrid(2))
    spec_lr = hntr_spec(chunkerlr%chunk_size(1), chunkerlr%ngrid(2), 0d0, 180d0*60d0 / chunkerlr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, 0d0)   ! datmis=0


    ! ---------------- INPUT: Open land cover
    ent2 = make_ent2()
    !Make file name.
    call chunker%file_info(info, ent20, LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1') 
    !Open the file
    call chunker%nc_open(ioall_lc, &
        chunker%outputs_dir, trim(info%dir), trim(info%leaf)//'.nc', trim(info%vname), 0)
    !Get handles to snow and ice land cover, mvs translates indices from full
    !list of lc indices.  Assumes 3D array (IM, JM, #lc types); /1,1 are just
    !placeholders for 3D.
    do k=1,2
        call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,ent2%mvs(k)/))
    end do


    ! --------------- OUTPUT: Open output file, create low-res (lr) version.
    call chunkerlr%file_info(info, ent2, LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')
    call chunkerlr%nc_create(ioall_lcout, &
        weighting(chunkerlr%wta1,1d0,0d0), &
        trim(info%dir), trim(info%leaf), trim(info%vname), &
        'Land Cover Fractions', '1', &
        ent2%layer_names(), ent2%long_layer_names())!, create_lr=.false.)
    do k=1,2
        call chunkerlr%nc_reuse_var( &
            ioall_lcout, io_lcout(k), (/1,1,k/), &
            weighting(chunkerlr%wta1, 1d0, 0d0))
    end do

    !Check number of file handles opened does not excide those initialized inside chunker.
    call chunker%nc_check(rw=rw)
    call chunkerlr%nc_check(rw=rw)
    !Write .mk file containing names of input and output files.
    call rw%write_mk

#ifdef JUST_DEPENDENCIES
    !Just write the dependencies for the record, otherwise do not compute anything.
    stop 0
#endif

    ! -------------- Regrid!
#ifdef ENTGVSD_DEBUG
    !Just do a few chunks when debugging.
    do jchunk = dbj0,dbj1
    do ichunk = dbi0,dbi1
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif
        !Set up chunk location in chunker, and read chunk from input file.
        call chunker%move_to(ichunk,jchunk)
        call chunkerlr%move_to(ichunk,jchunk)

        ! Convert away from FillValue so averaging / regridding works.
        do k=1,2
        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)
            if (io_lc(k)%buf(ic,jc) == FillValue) then
                io_lc(k)%buf(ic,jc) = 0  !ONLY SET TO ZERO TO REGRID LC, not lai or hgt
!            else
!                print *,ic,jc,io_lc(k)%buf(ic,jc)
            end if
        end do
        end do
        end do        

#if 1
        !Now call hntr to regrid.
        do k=1,2
!            print *,io_lcout(k)%startB(2), io_lcout(k)%chunker%chunk_size(2)
            print *,io_lc(k)%buf(1,1),sum(io_lc(k)%buf)
            call hntr_lr%regrid4( &
                io_lcout(k)%buf, io_lc(k)%buf, &
                chunker%wta1, 1d0, 0d0, &   ! weighting = wta1*1d0 + 0d0
                io_lcout(k)%startB(2), io_lcout(k)%chunker%chunk_size(2)) !startB(2)=j of lat
        end do
#else
        do k=1,2
            io_lcout(k)%buf(:,:) = 0
        end do
#endif
        !Write the chunk.
        !call chunker%write_chunks
        call chunkerlr%write_chunks

    end do
    end do


end program B03_regrid_snowice
