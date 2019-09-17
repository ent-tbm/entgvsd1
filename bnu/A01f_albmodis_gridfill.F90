module a01f_gridfill_mod

    use carrer_mod
    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use geom_mod
    use assign_laimax_mod
    use gridfill_mod
    
    implicit none
CONTAINS

subroutine do_gridfill(iband)
    integer :: iband

    ! ---------- Locals
    integer :: k, m, b, dimlon, dimlat, dim(2), dimband

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: ioall_lc, io_lcice, io_lcwater
    type(ChunkIO_t) :: ioall_albmodis, io_albmodis  ! io_albmodis is SMEAN
    ! Output Files
    type(ChunkIO_t) :: io_albfill ! sam as io_albmodis


    INTEGER :: ichunk,jchunk,ic,jc
    real, dimension(:,:), allocatable, target :: wta    ! Surmised landmask
    type(EntSet_t) :: ent2
    type(FileInfo_t) :: info

    ! ---------- Parameters for gridfill
    integer:: nlat, nlon
    real*8, dimension(:,:), allocatable :: grid
    real*8 :: FillValue8
    integer :: itermax
    real*8 :: tolerance, relaxc
    logical :: initzonal
    real*8 :: resmax   ! OUT
    integer :: numiter ! OUT

    print *,'=========== BEGIN do_gridfill',iband
    call ent2%allocate(2,NENT20)
    call ent2%sub_covertype(ent20, SNOW_ICE)
    call ent2%sub_covertype(ent20, CV_WATER)

    call chunker%init(IMK, JMK, IMH*2,JMH*2, 'forplot', 100, 100, 10, (/1,1/))
    allocate(wta(chunker%chunk_size(1), chunker%chunk_size(2)))

    ! ------------- Open Input Files

    ! --------- LC
    call chunker%file_info(info, ent20, 'BNU', 'M', 'lc', 2004, 'ent17', '1.1')
    call chunker%nc_open(ioall_lc, &
        LC_LAI_ENT_DIR, trim(info%dir), trim(info%leaf)//'.nc', trim(info%vname), 0)
    call chunker%nc_reuse_var(ioall_lc, io_lcice, (/1,1,ent2%svm(SNOW_ICE)/))
    call chunker%nc_reuse_var(ioall_lc, io_lcwater, (/1,1,ent2%svm(CV_WATER)/))

    ! ------------ albmodis
    call chunker%nc_open(ioall_albmodis, LC_LAI_ENT_DIR, &
        'carrer/', &
        'albmodis_'//trim(sbands_modis(iband))//'.nc', &
        'albmodis_'//trim(sbands_modis(iband)), 0)
    call chunker%nc_reuse_var( &
        ioall_albmodis, io_albmodis, &
        (/1,1,SMEAN/))

    ! ===================== Open Output Files

    ! ------------ albfill
    call chunker%nc_create(io_albfill, weighting(wta,1d0,0d0), &
        'carrer/', &
        'albfill_'//trim(sbands_modis(iband)), &
        'albfill_'//trim(sbands_modis(iband))//'_MEAN', &
        'albmodis MEAN with missing values filled in', '1')


    call chunker%nc_check(MAIN_PROGRAM_FILE)

    ! ================= Inputs for gridfill
    ! https://ocefpaf.github.io/python4oceanographers/blog/2014/10/20/gridfill/
    allocate(grid(chunker%chunk_size(2), chunker%chunk_size(1)))
    FillValue8 = FillValue
    itermax = 10000
    tolerance = .01
    relaxc = 0.5
    initzonal = .true.

    ! ================== Main Loop

    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)

        call chunker%move_to(ichunk,jchunk)
        wta = 1

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)
            grid(jc,ic) = io_albmodis%buf(ic,jc)   ! Transpose for gridfill
            ! We will want to fill gridcells that are not 100% land
            if (io_lcice%buf(ic,jc) + io_lcwater%buf(ic,jc) > 1.e-5) then
                grid(jc,ic) = FillValue8
            end if
        end do
        end do

        print *,'BEGIN poisson_fill()'
        call poisson_fill( &
            chunker%chunk_size(2), chunker%chunk_size(1), &
            grid, FillValue8, itermax, tolerance, &
            relaxc, initzonal, .true., &
            ! -- Output &
            resmax, numiter)
        print *,'END poisson_fill()',numiter,resmax
                
        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)
!            if (abs(1.0 - io_lcice%buf(ic,jc) + io_lcwater%buf(ic,jc)) < 1.e-5) then
!                ! Cell is all ice or land; Albedo calculator doesn't need a filled-in value
!                io_albfill%buf(ic,jc) = FillValue
!            else
                io_albfill%buf(ic,jc) = grid(jc,ic)
!            end if
        end do
        end do

        call chunker%write_chunks
    end do
    end do


    call chunker%close_chunks
    deallocate(grid)

end subroutine do_gridfill

end module a01f_gridfill_mod
         
program Carrer_soilalbedo_gridfill
    use a01f_gridfill_mod
    use paths_mod
    implicit none

    integer :: iband

    MAIN_PROGRAM_FILE = 'A01f_albmodis_gridfill'
    call init_ent_labels
    do iband=1,NBANDS_MODIS
        call do_gridfill(iband)
    end do

end program Carrer_soilalbedo_gridfill
