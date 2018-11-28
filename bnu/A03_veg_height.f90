!simard.f
!Assign tree heights from Simard et al. (2011) to tree PFTs in 
! EntMM_lc_maxlai_#x#.bin generated by modis_entpftg.f.

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A03_veg_height.f

! gfortran -o myExe arrayutil.o A03_veg_height.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

! ./myExe

!-----------------------------------------------------------------
      
program simard

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod

implicit none
      
    integer, parameter :: ENTPFTNUM = 19 !17 Ent PFTs + snow/ice + bare

    integer, parameter :: IMH = 720 !long at 0.5 degrees
    integer, parameter :: JMH = 360 !lat at 0.5 degrees

    integer, parameter :: X1km = 43200 !long at 1 km
    integer, parameter :: Y1km = 21600 !lat at 1 km

    integer, parameter :: IM1km = X1km !long at 1 km
    integer, parameter :: JM1km = Y1km !lat at 1 km

    type(Chunker_t) :: chunker
    integer :: jchunk, ichunk    ! Index of current chunk
    integer :: jc, ic    ! Index WITHIN current chunk
    integer :: jj, ii            ! Index in full space

    ! ------ Inputs
    type(ChunkIO_t) :: io_sim    ! Simard heights
    type(ChunkIO_t) :: ioall_lc    ! Simard heights
    type(ChunkIO_t), dimension(ENTPFTNUM) :: io_lc
    ! ------ Outputs
    type(ChunkIO_t) :: ioall_out
    type(ChunkIO_t), dimension(ENTPFTNUM) :: io_out

    real*4 :: SHEIGHT, OHEIGHT

    integer :: i,j,k,p



! Table assigns standard height, based on Simard Height and PFT index
! <Assigned height>(i) = SHEIGHT * std_heights(i*2) + std_heights(i*2+1)
!  ...where SHEIGHT is the Simard height from the input file.
real*8, parameter :: heights_form(2,ENTPFTNUM) = RESHAPE( (/ &
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
    0d0, 0.0d0   & ! BARREN 19 - Bare or sparsely vegetated, urban
/), (/ 2,ENTPFTNUM /))



!**   INPUT Files at 1km x 1km 
call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)


! ----------------------------------------------------------------------
!     GET NC FILES IDs
      
!     SIMARD HEIGHTS
call chunker%nc_open_gz(io_sim, &
    DATA_DIR, DATA_INPUT, &
    'height/', &
    'simard_forest_heights.nc', 'heights', 1)

!     ENTPFTLC
do k = 1,ENTPFTNUM
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(ent19(k)%file1)//trim(ent19(k)%file2)//'_lc.nc', &
        trim(ent19(k)%file2), 1)
end do

! ENTPFT heights
call chunker%nc_create(ioall_out, &
    weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
    'height/', 'EntGVSDmosaic17_height_1kmx1km_lai3g', &
    'SimardHeights', &
    'Plant Heights', 'm', 'Plant Heights', &
    index_array(ent19), title_array(ent19))
do p = 1,ENTPFTNUM
    call chunker%nc_reuse_var(ioall_out, io_out(p), &
        (/1,1,p/), 'w', weighting(io_lc(p)%buf,1d0,0d0))
end do

! Quit if we had any problems opening files
call chunker%nc_check('A03_veg_heights')
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


        do p = 1,ENTPFTNUM
            OHEIGHT = undef   ! Default if nothing in this PFT
            if (io_lc(p)%buf(ic,jc) > 0) then
                ! Lookup what the height should be
                SHEIGHT = io_sim%buf(ic,jc)
                OHEIGHT = SHEIGHT * heights_form(1,p) + heights_form(2,p)
            end if

            io_out(p)%buf(ic,jc) = OHEIGHT
        end do    ! p

    end do    ! ic
    end do    ! jc

    
    call chunker%write_chunks

end do ! ichunk
end do ! jchunk

call chunker%close_chunks

      
end program simard

