!simard.f
!Assign tree heights from Simard et al. (2011) to tree PFTs in 
! EntMM_lc_maxlai_#x#.bin generated by modis_entpftg.f.
!
! It does the same thing as BNU/A01*, but instead of assigning LAImax,
! it assigns height to the cover fractions.  Both LAImax and Height
! are 1 km scale assignments to the 1 km subgrid cover fractions (all
! subgrid fractions in a 1 km grid cell get the same LAImax and
! Height, except for shrubs and grasses, which get assigned some
! ballpark heights).
!
!-----------------------------------------------------------------
      
program simard

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use geom_mod

implicit none
      
    type(Chunker_t) :: chunker
    integer :: jchunk, ichunk    ! Index of current chunk
    integer :: jc, ic    ! Index WITHIN current chunk
    integer :: jj, ii            ! Index in full space

    integer, parameter :: NENT19 = NENT20-1
    type(EntSet_t) :: ent19

    ! ------ Inputs
    type(ChunkIO_t) :: io_sim    ! Simard heights
    type(ChunkIO_t) :: ioall_lc    ! Simard heights
    type(ChunkIO_t), dimension(NENT20) :: io_lc
    ! ------ Outputs
    type(ChunkIO_t) :: ioall_out
    type(ChunkIO_t), dimension(NENT20) :: io_out

    real*4 :: SHEIGHT, OHEIGHT

    integer :: k

! Table assigns standard height, based on Simard Height and PFT index
! <Assigned height>(i) = SHEIGHT * std_heights(i*2) + std_heights(i*2+1)
!  ...where SHEIGHT is the Simard height from the input file.
real*8, parameter :: heights_form(2,NENT19) = RESHAPE( (/ &
    1d0, 0d0,    & ! TREE    1 - evergreen broadleaf early successional
    1d0, 0d0,    & ! TREE    2 - evergreen broadleaf late successional
    1d0, 0d0,    & ! TREE    3 - evergreen needleleaf early successional
    1d0, 0d0,    & ! TREE    4 - evergreen needleleaf late successional
    1d0, 0d0,    & ! TREE    5 - cold deciduous broadleaf early successional
    1d0, 0d0,    & ! TREE    6 - cold deciduous broadleaf late successional
    1d0, 0d0,    & ! TREE    7 - drought deciduous broadleaf
    1d0, 0d0,    & ! TREE    8 - deciduous needleleaf
    0d0, 0.365d0,& ! SHRUB   9 - cold adapted shrub
    0d0, 2.0d0,  & ! SHRUB  10 - arid adapted shrub
    0d0, 1.5d0,  & ! GRASS  11 - C3 grass perennial
    0d0, 1.5d0,  & ! GRASS  12 - C4 grass
    0d0, 0.5d0,  & ! GRASS  13 - C3 grass - annual
    0d0, 0.5d0,  & ! GRASS  14 - arctic C3 grass
    0d0, 0.5d0,  & ! HERB   15 - crops C3 herb
    0d0, 0.5d0,  & ! HERB   16 - crops C4 herb
    1d0, 0d0,    & ! TREE   17 - crops woody
    0d0, 0.0d0,  & ! BARREN 18 - Permanent snow/ice
    0d0, 0.0d0,  & ! BARREN 19 - Bare or sparsely vegetated, urban
    0d0, 0.0d0   & ! BARREN 20 - water
/), (/ 2,NENT19 /))



call init_ent_labels

! -----------------------------------------------------
! Set up ent19 = ent20 without water
! (indices will be the same, but only if water is at the end)
call ent19%allocate(NENT20-1, NENT20)
do k=1,NENT20
    if (k == CV_WATER) cycle
    call ent19%sub_covertype(ent20, k)
end do
! -----------------------------------------------------

call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

! ----------------------------------------------------------------------
!     GET NC FILES IDs

! ------------- Inputs     
!     SIMARD HEIGHTS
call chunker%nc_open_gz(io_sim, &
    DATA_DIR, DATA_INPUT, &
    'height/', &
    'simard_forest_heights.nc', 'heights', 1)

!     ENTPFTLC
do k = 1,NENT19
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', &
        itoa2(ent19%mvs(k))//'_'//trim(ent19%abbrev(k))//'_lc.nc', &
        trim(ent19%abbrev(k)), 1)
end do

! ---------------- Outputs
! ENTPFT heights
call chunker%nc_create(ioall_out, &
    weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
    'height/', 'EntGVSDmosaic17_height_1kmx1km_lai3g', &
    'SimardHeights', &
    'Plant Heights', 'm', 'Plant Heights', &
    ent19%mvs, ent19%layer_names())

do k = 1,NENT19
    call chunker%nc_reuse_var(ioall_out, io_out(k), &
        (/1,1,k/), 'w', weighting(io_lc(k)%buf,1d0,0d0))
end do

! Quit if we had any problems opening files
call chunker%nc_check('A01h_veg_heights')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

! --------------------------------- main loop
! Use these loop bounds for testing...
! it chooses a land area in Asia

#ifdef ENTGVSD_DEBUG
do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
do jchunk = 1,nchunk(2)
do ichunk = 1,nchunk(1)
#endif

    call chunker%move_to(ichunk,jchunk)

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

        ! Compute overall NetCDF index of current cell
        ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
        jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1


        do k = 1,NENT19
            OHEIGHT = undef   ! Default if nothing in this PFT
            if (io_lc(k)%buf(ic,jc) > 0) then
                ! Lookup what the height should be
                SHEIGHT = io_sim%buf(ic,jc)
                OHEIGHT = SHEIGHT * heights_form(1,k) + heights_form(2,k)

                ! TODO: Better way: inverse weighted average of
                ! gridcell height, so tree heights are taller.  So
                ! average gridcell height (averaging over non-zero
                ! LC's) is same as the Simard height.
            end if

            io_out(k)%buf(ic,jc) = OHEIGHT
        end do    ! p

    end do    ! ic
    end do    ! jc
    
    call chunker%write_chunks

end do ! ichunk
end do ! jchunk

call chunker%close_chunks

      
end program simard

