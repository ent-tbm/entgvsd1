!  Regrids trimmed_scaled_natveg files to grid of GISS GCM ModelE for later reformatting as GCM input.
!
! Step after this will take regridded files to format directly useable by ModelE for two reasons:
!   1. Convert arrays with an lctype dimension to separate array variables with
!   netcdf names for the cover type.
!   2. Convert from NetCDF4 format to NetCDF3 for ModelE input format.
!
! Author: Elizabeth Fischer

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module b18_mod
!Author: Elizabeth Fischer
!
! Regrids higher spatial resolution (e.g. HXH) files to GISS GCM ModelE 2HX2 grid.

    use netcdf
    use chunker_mod
    use ent_labels_mod
    use gcm_labels_mod
    use ent_params_mod
    use hntr_mod

implicit none
CONTAINS


subroutine make_modele(rw, esub, lcpart0, laipart0, hgtpart0, part1, &
    im1,jm1, N_BARE, N_WATERICE, modification)
    type(ReadWrites_t) :: rw
    type(EntSet_t), intent(IN) :: esub
    character(*), intent(IN) :: lcpart0
    character(*), intent(IN) :: laipart0
    character(*), intent(IN) :: hgtpart0
    character(*), intent(IN) :: part1    ! Destination section
    integer, intent(IN) :: im1,jm1    ! Destination resolution
    integer, intent(IN) :: N_BARE, N_WATERICE
    character*(*), intent(IN) :: modification
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

    integer :: i,j, ii, jj, shapebuf(2)
    real*8 :: sum_vfc
    
    call clear_file_info(om)

    om%modification = modification
! For each make_modele(), mods should be one of:
!natveg = crop cover set to zero and other cover scaled up to fill grid cell;  if no other cover type in grid cell, dominant cover of adjacent grid cells is used.
!ext# = LAI of crops extended to # neighboring grid cells to provide LAI for historically changing crop cover.
!ext1 = cover and ext LAI further extended by latitude across continental boundaries for ModelE users who use slightly different continental boundaries.
!other = bug fix for TBD, etc. >;





    call chunker0%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 300,   1, 17, (/1,1/), outputs_dir=THIS_OUTPUTS_DIR)
    call chunker1%init(im1,jm1,  im1,jm1, 'forplot',       1, 300, 17, (/1,1/), outputs_dir=THIS_OUTPUTS_DIR)

    spec0 = hntr_spec(chunker0%chunk_size(1), chunker0%ngrid(2), 0d0, 180d0*60d0 / chunker0%ngrid(2))
    spec1 = hntr_spec(chunker1%chunk_size(1), chunker1%ngrid(2), 0d0, 180d0*60d0 / chunker1%ngrid(2))
    hntr = hntr_calc(spec1, spec0, 0d0)   ! datmis=0


    ! ---------------- Input Files
    call chunker0%nc_open_set( &
        esub, io0_ann_lc(:,1), &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, lcpart0, '1.1')

    call chunker0%nc_open_set( &
        esub, io0_ann_lai(:,1), &
        LAI_SOURCE, 'M', 'laimax', LAI_YEAR, laipart0, '1.1')

    call chunker0%nc_open_set( &
        esub, io0_ann_hgt(:,1), &
        LAI_SOURCE, 'M', 'hgt', LAI_YEAR, hgtpart0, '1.1')

    do m=1,NMONTH
        call chunker0%nc_open_set( &
            esub, io0_mon_lai(:,m), &
            LAI_SOURCE, 'M', 'lai', LAI_YEAR, laipart0, '1.1', &
            doytype='month', idoy=m)
    end do

    ! ---------------- Output Files
    call chunker1%nc_create_set( &
        esub, io1_ann_lc(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, part1, '1.1', create_lr=.false., overmeta=om)

    call chunker1%nc_create_set( &
        !esub, io1_ann_lai(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        esub, io1_ann_lai(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'laimax', LAI_YEAR, part1, '1.1', create_lr=.false., overmeta=om)

    call chunker1%nc_create_set( &
        esub, io1_ann_hgt(:,1), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'hgt', LAI_YEAR, part1, '1.1', create_lr=.false., overmeta=om)

    do m=1,NMONTH
        call chunker1%nc_create_set( &
            esub, io1_mon_lai(:,m), repeat_weights(esub%ncover, chunker1%wta1, 1d0, 0d0), &
            LAI_SOURCE, 'M', 'lai', LAI_YEAR, part1, '1.1', create_lr=.false., overmeta=om, &
            doytype='month', idoy=m)
    end do

    call chunker0%nc_check(rw=rw)
    call chunker1%nc_check(rw=rw)

#ifdef JUST_DEPENDENCIES
    return
#endif

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
                io1_ann_lc(k,1)%buf, io0_ann_lc(k,1)%buf, &
                chunker0%wta1, 1d0, 0d0, &   ! weighting
                io1_ann_lc(k,1)%startB(2), io1_ann_lc(k,1)%chunker%chunk_size(2))

!            if ((k.eq.esub%crops_herb).or.(k.eq.esub%crops_woody)) then
              call hntr%regrid4( &
                io1_ann_lai(k,1)%buf, io0_ann_lai(k,1)%buf, &
                io0_ann_lc(k,1)%buf, 1d0, 0d0, &   ! weighting
                io1_ann_lai(k,1)%startB(2), io1_ann_lai(k,1)%chunker%chunk_size(2))

              call hntr%regrid4( &
                io1_ann_hgt(k,1)%buf, io0_ann_hgt(k,1)%buf, &
                io0_ann_lc(k,1)%buf, 1d0, 0d0, &   ! weighting
                io1_ann_hgt(k,1)%startB(2), io1_ann_hgt(k,1)%chunker%chunk_size(2))

              do m=1,NMONTH
                call hntr%regrid4( &
                    io1_mon_lai(k,m)%buf, io0_mon_lai(k,m)%buf, &
                    io0_ann_lc(k,1)%buf, 1d0, 0d0, &   ! weighting
                    io1_mon_lai(k,m)%startB(2), io1_mon_lai(k,m)%chunker%chunk_size(2))
              end do

! !           else  !non-crop cover
! !             call hntr%regrid4( &
! !               io1_ann_lai(k,1)%buf, io0_ann_lai(k,1)%buf, &
! !               io0_ann_lc(k,1)%buf, 1d0, 0d0, &   ! weighting
! !               io1_ann_lai(k,1)%startB(2), io1_ann_lai(k,1)%chunker%chunk_size(2))
! !
!              call hntr%regrid4( &
!                io1_ann_hgt(k,1)%buf, io0_ann_hgt(k,1)%buf, &
!                io0_ann_lc(k,1)%buf, 1d0, 0d0, &   ! weighting
!                io1_ann_hgt(k,1)%startB(2), io1_ann_hgt(k,1)%chunker%chunk_size(2))
!
!              do m=1,NMONTH
!                call hntr%regrid4( &
!                    io1_mon_lai(k,m)%buf, io0_mon_lai(k,m)%buf, &
!                    io0_ann_lc(k,1)%buf, 1d0, 0d0, &   ! weighting
!                    io1_mon_lai(k,m)%startB(2), io1_mon_lai(k,m)%chunker%chunk_size(2))
!              end do
!            endif
        end do

        ! Rescale lc fractions so that the land fractions sum to 1 and exclude water_ice in same pixel.
        ! Will leave water_ice that is 100% of grid cell.
        !N_BARE = esub%ncover - 1!esub%bare_dark. Crappy programming, cannot access cover types with EntSet_t
        !N_WATERICE = esub%ncover !esub%water_ice
        shapebuf = shape(io1_ann_lc(1,1)%buf)
        ii = shapebuf(1)
        jj = shapebuf(2)
        print *, 'ii, jj', ii, jj
        do i=1,ii
           do j=1,jj
              sum_vfc = 0.
              do k=1,N_BARE
                 sum_vfc = sum_vfc + io1_ann_lc(k,1)%buf(i,j)
              enddo
              !grid cells with land scale land to sum to 1 and set water_ice to 0.
             !If all land got trimmed out, then set water_ice to 1. 
              if ( sum_vfc > 0. ) then
                 do k=1,N_BARE
                    io1_ann_lc(k,1)%buf(i,j) = io1_ann_lc(k,1)%buf(i,j)/ sum_vfc
                 enddo
                 io1_ann_lc(N_WATERICE,1)%buf(i,j) = 0.
              else if ((io1_ann_lc(N_WATERICE,1)%buf(i,j).gt.0.).and.(io1_ann_lc(N_WATERICE,1)%buf(i,j).lt.1.)) then
                 print *, 'ERR with water_ice: sum_vfc, vfc', sum_vfc, io1_ann_lc(N_WATERICE,1)%buf(i,j)
                 io1_ann_lc(N_WATERICE,1)%buf(i,j) = 1.
              end if
              !sum_vfc = sum(vfc(1:N_BARE))
           enddo
        enddo

        call chunker1%write_chunks
     end do
  end do

    call chunker0%close_chunks
    call chunker1%close_chunks

end subroutine make_modele



end module b18_mod

! =========================================================
program regrid
    use b18_mod
implicit none
    type(GcmEntSet_t), target :: esub
    class(EntSet_t), pointer :: esub_p
    type(ReadWrites_t) :: rw
    integer :: N_BARE, N_WATERICE

    call rw%init(THIS_OUTPUTS_DIR, "B18_modele", 300,300)
    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    N_WATERICE = esub%water_ice
    N_BARE = N_WATERICE - 1 !esub%bare_dark, assume water_ice is last cover type and only non-land.
    esub_p => esub

    ! Here we put different combinations for ModelE consumption
    call make_modele(rw, esub_p, &
        'trimmed_scaled_natveg', &    ! LC
        'trimmed_scaled_natveg', &            ! LAIMAX, LAI
        'trimmed_scaled_natveg', &            ! HGT
        'modelE_natveg', & ! output part
        IM2,JM2, N_BARE, N_WATERICE, &
        'reformat for GISS ModelE, natveg, no crops')


    call make_modele(rw, esub_p, &
        'trimmed_scaled', &    ! LC
        'trimmed_scaled', &            ! LAIMAX, LAI
        'trimmed_scaled', &            ! HGT
        'modelE', & ! output part
        IM2,JM2, N_BARE, N_WATERICE, &
        'reformat for GISS ModelE, retain observed crop cover')


    call rw%write_mk
end program regrid
