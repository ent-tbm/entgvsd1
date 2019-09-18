module a10_mod
    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
    use hntr_mod

implicit none
CONTAINS


subroutine make_modele(rw, esub, lcpart0, laipart0, hgtpart0, part1, im1,jm1, mods)
    type(ReadWrites_t) :: rw
    type(EntSet_t), intent(IN) :: esub
    character(*), intent(IN) :: lcpart0
    character(*), intent(IN) :: laipart0
    character(*), intent(IN) :: hgtpart0
    character(*), intent(IN) :: part1    ! Destination section
    integer, intent(IN) :: im1,jm1    ! Destination resolution
    character*(*), intent(IN) :: mods
    ! ----------- Locals
    integer :: ichunk,jchunk
    integer :: m,k
    type(Chunker_t) :: chunker0    ! Original resolution
    type(Chunker_t) :: chunker1    ! Destination resolution

    type(ChunkIO_t) :: io0_ann_lc(esub%ncover,1)
    type(ChunkIO_t) :: io0_ann_lai(esub%ncover,1)
    type(ChunkIO_t) :: io0_ann_hgt(esub%ncover,1)
    type(ChunkIO_t) :: io0_mon_lai(esub%ncover,NMONTH)

    type(ChunkIO_t) :: io1_ann_lc(esub%ncover,1)
    type(ChunkIO_t) :: io1_ann_lai(esub%ncover,1)
    type(ChunkIO_t) :: io1_ann_hgt(esub%ncover,1)
    type(ChunkIO_t) :: io1_mon_lai(esub%ncover,NMONTH)

    type(HntrSpec_t) :: spec0, spec1
    type(HntrCalc_t) :: hntr   ! Preparation to regrid

    type(FileInfo_t) :: om   ! overmeta

    call clear_file_info(om)

    om%modifications = ''
! For each make_modele(), mods should be one of:
!nocrops = crop cover set to zero and other cover scaled up to fill grid cell;  if no other cover type in grid cell, dominant cover of adjacent grid cells is used.
!ext# = LAI of crops extended to # neighboring grid cells to provide LAI for historically changing crop cover.
!ext1 = cover and ext LAI further extended by latitude across continental boundaries for ModelE users who use slightly different continental boundaries.
!other = bug fix for TBD, etc. >;





    call chunker0%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 200, 1, 5, (/1,1/))
    call chunker1%init(im1,jm1,  im1,jm1, 'forplot', 1, 200, 5, (/1,1/))

    spec0 = hntr_spec(chunker0%chunk_size(1), chunker0%ngrid(2), 0d0, 180d0*60d0 / chunker0%ngrid(2))
    spec1 = hntr_spec(chunker1%chunk_size(1), chunker1%ngrid(2), 0d0, 180d0*60d0 / chunker1%ngrid(2))
    hntr = hntr_calc(spec0, spec1, 0d0)   ! datmis=0


    ! ---------------- Input Files
    call chunker0%nc_open_set( &
        esub, io0_ann_lc(:,1), &
        'B', 'M', 'lc', 2004, lcpart0, '1.1')

    call chunker0%nc_open_set( &
        esub, io0_ann_lai(:,1), &
        'B', 'M', 'laimax', 2004, laipart0, '1.1')

    call chunker0%nc_open_set( &
        esub, io0_ann_hgt(:,1), &
        'B', 'M', 'hgt', 2004, hgtpart0, '1.1')

    do m=1,NMONTH
        call chunker0%nc_open_set( &
            esub, io0_mon_lai(:,m), &
            'B', 'M', 'lai', 2004, laipart0, '1.1', &
            doytype='month', idoy=m)
    end do

    ! ---------------- Output Files
    call chunker1%nc_create_set( &
        esub, io1_ann_lc(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        'B', 'M', 'lc', 2004, part1, '1.1', create_lr=.false., overmeta=om)

    call chunker1%nc_create_set( &
        esub, io1_ann_lai(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        'B', 'M', 'laimax', 2004, part1, '1.1', create_lr=.false., overmeta=om)

    call chunker1%nc_create_set( &
        esub, io1_ann_hgt(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        'B', 'M', 'hgt', 2004, part1, '1.1', create_lr=.false., overmeta=om)

    do m=1,NMONTH
        call chunker1%nc_create_set( &
            esub, io1_mon_lai(:,m), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
            'B', 'M', 'lai', 2004, part1, '1.1', create_lr=.false., overmeta=om, &
            doytype='month', idoy=m
    end do

    call chunker0%nc_check(rw=rw)
    call chunker1%nc_check(rw=rw)

#ifdef ENTGVSD_DEBUG
    do jchunk = 1,1
    do ichunk = 1,1
#else
    do jchunk = 1,chunker0%nchunk(2)
    do ichunk = 1,chunker0%nchunk(1)
#endif
        call chunker0%move_to(ichunk,jchunk)
        call chunker1%move_to(ichunk,jchunk)

        do k=1,esub%ncover
            call hntr%regrid4( &
                io0_ann_lc(k,1)%buf, io1_ann_lc(k,1)%buf, &
                chunker0%wta1, 1d0, 0d0, &   ! weighting
                io0_ann_lc(k,1)%startB(2), io0_ann_lc(k,1)%chunker%chunk_size(2))

            call hntr%regrid4( &
                io0_ann_lai(k,1)%buf, io1_ann_lai(k,1)%buf, &
                chunker0%wta1, 1d0, 0d0, &   ! weighting
                io0_ann_lai(k,1)%startB(2), io0_ann_lai(k,1)%chunker%chunk_size(2))

            call hntr%regrid4( &
                io0_ann_hgt(k,1)%buf, io1_ann_hgt(k,1)%buf, &
                chunker0%wta1, 1d0, 0d0, &   ! weighting
                io0_ann_hgt(k,1)%startB(2), io0_ann_hgt(k,1)%chunker%chunk_size(2))

            do m=1,NMONTH
                call hntr%regrid4( &
                    io0_mon_lai(k,m)%buf, io1_mon_lai(k,m)%buf, &
                    chunker0%wta1, 1d0, 0d0, &   ! weighting
                    io0_mon_lai(k,m)%startB(2), io0_mon_lai(k,m)%chunker%chunk_size(2))
            end do
        end do

        call chunker1%write_chunks

    end do
    end do

    call chunker0%close_chunks
    call chunker1%close_chunks

end subroutine make_modele



end module a10_mod

! =========================================================
program regrid
    use a10_mod
implicit none
    type(GcmEntSet_t), target :: esub
    class(EntSet_t), pointer :: esub_p
    type(ReadWrites_t) :: rw

    call rw%init("A10_modele", 300,300)
    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    esub_p => esub

    ! Here we put different combinations for ModelE consumption
    call make_modele(rw, esub_p, &
        'trimmed_scaled_nocrops', &    ! LC
        'trimmed_scaled', &            ! LAIMAX, LAI
        'trimmed_scaled', &            ! HGT
        'modele1', & ! output part
        IM2,JM2)



    call rw%write_mk
end program regrid
