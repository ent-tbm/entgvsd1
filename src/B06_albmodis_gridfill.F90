! Fills in small regions of missing data in the soil albedo files.
!
! Uses Poisson relaxation, similar to NCL's poisson_grid_fill()
! function
! https://www.ncl.ucar.edu/Document/Functions/Built-in/poisson_grid_fill.shtml
!
! Author: Elizabeth Fischer
!
! See: slib/gridfill.f90

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module b06_gridfill_mod

    use carrer_mod
    use netcdf
    use chunker_mod
    use ent_labels_mod
    use ent_params_mod
    use assign_laimax_mod
    use gridfill_mod

    implicit none
CONTAINS

subroutine do_gridfill(rw, iband)
    type(ReadWrites_t) :: rw
    integer :: iband

    ! ---------- Locals
    integer :: k, m, b, dimlon, dimlat, dim(2), dimband

    type(Chunker_t) :: chunker
    ! Input files
    type(ChunkIO_t) :: ioall_lc, io_lcice, io_lcwater
    type(ChunkIO_t) :: ioall_albmodis, io_albmodis(NSTATS_FILLABLE)  ! io_albmodis is SMEAN
    ! Output Files
    type(ChunkIO_t) :: ioall_albfill, io_albfill(NSTATS_FILLABLE) ! sam as io_albmodis


    INTEGER :: ichunk,jchunk,ic,jc,id
    real, dimension(:,:,:), allocatable, target :: wta    ! Surmised landmask
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

    real*8 :: val

    ! Latitude extent of northern and southern hemisphere oceans
    integer :: max_SH_ocean, max_NH_ocean
    ! temporaries
    integer :: max_NH_land
    logical :: allocean,SH_isset

    print *,'=========== BEGIN do_gridfill',iband
    call chunker%init(IMK, JMK, IMH*2,JMH*2, 'forplot', 100, 100, 10, (/1,1/), outputs_dir=THIS_OUTPUTS_DIR)
    allocate(wta(chunker%chunk_size(1), chunker%chunk_size(2), NSTATS_FILLABLE ))

    ! ------------- Open Input Files

    ! --------- LC
    ent2 = make_ent2()
    call chunker%file_info(info, ent2, LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')
    call chunker%nc_open(ioall_lc, &
        chunker%outputs_dir, trim(info%dir), trim(info%leaf)//'.nc', trim(info%vname), 0)
    call chunker%nc_reuse_var(ioall_lc, io_lcice, (/1,1,ent2%svm(SNOW_ICE)/))
    call chunker%nc_reuse_var(ioall_lc, io_lcwater, (/1,1,ent2%svm(CV_WATER)/))

    ! ------------ albmodis
    call chunker%nc_open(ioall_albmodis, chunker%outputs_dir, &
        'tmp/carrer/', &
        'albmodis_'//trim(sbands_modis(iband))//'.nc', &
        'albmodis_'//trim(sbands_modis(iband)), 0)
    do id=1,NSTATS_FILLABLE
        call chunker%nc_reuse_var( &
            ioall_albmodis, io_albmodis(id), &
            (/1,1,id/))
    end do

    ! ===================== Open Output Files

    ! ------------ albfill
     call clear_file_info(info)
     info%vname = 'albfill_'//trim(sbands_modis(iband))
     info%long_name = 'albmodis statistics with missing values filled in'
     info%units = '1'
     info%file_metadata_type = 'carrer'
     call chunker%nc_create1(ioall_albfill, weighting(wta(:,:,1),1d0,0d0), &
         'tmp/carrer/', &
         'albfill_'//trim(sbands_modis(iband)), info, &
         sstats(1:NSTATS_FILLABLE), sstats(1:NSTATS_FILLABLE))
     do id=1,NSTATS_FILLABLE
         call chunker%nc_reuse_var(ioall_albfill, io_albfill(id), (/1,1,id/), &
             weighting(wta(:,:,id), 1d0, 0d0))
     end do

    call chunker%nc_check(rw=rw)
#ifdef JUST_DEPENDENCIES
    return
#endif

    ! ================= Inputs for gridfill
    ! https://ocefpaf.github.io/python4oceanographers/blog/2014/10/20/gridfill/
    allocate(grid(chunker%chunk_size(2), chunker%chunk_size(1)))
    FillValue8 = FillValue
    itermax = 10000
    tolerance = .01
    relaxc = 0.5
    initzonal = .true.

    ! ================== Main Loop

    ! Just one chunk, but we "loop" anyway.
    ! There is one chunk because the Poisson gridfill doesn't work
    ! across boundaries
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)

        call chunker%move_to(ichunk,jchunk)

        do id=1,NSTATS_FILLABLE
            print *,'================ Filling ',trim(sbands_modis(iband)),'-',trim(sstats(id))
            ! NOTE: temporary variable grid is transposed compared to our standard.
            max_SH_ocean = -1   ! Northernmost SH latitude that is all ocean
            SH_isset = .false.
            max_NH_land = -1  ! Southermost NH latitude that is all ocean
            do jc = 1,chunker%chunk_size(2)
                allocean = .true.
                do ic = 1,chunker%chunk_size(1)
                    grid(jc,ic) = io_albmodis(id)%buf(ic,jc)   ! Transpose for gridfill
                    ! We will want to fill gridcells that are not 100% land
                    if (io_lcice%buf(ic,jc) + io_lcwater%buf(ic,jc) > 1.e-5) then
                        grid(jc,ic) = FillValue8
                    else
                        allocean = .false.
                    end if
                end do

                ! Compute max_SH_ocean and max_NH_ocean
                if (allocean) then
                    if (.not.SH_isset) then
                        max_SH_ocean = jc
                    end if
                else
                    SH_isset = .true.
                    ! Actually northermost latitude that is all land; convert below
                    max_NH_land = jc
                end if
            end do
            max_NH_ocean = max_NH_land + 1


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
                ! poisson_fill() created albedo of 0 at the poles, for latitudes
                ! with no data at all.  Replace this with FillValue.
                val = grid(jc,ic)    ! NOTE: grid is transposed compared to standard
                if ((val == 0) .or. (jc<=max_SH_ocean) .or. (jc>=max_NH_ocean)) then
                    io_albfill(id)%buf(ic,jc) = FillValue
                    wta(ic,jc,id)=0
                else
                    io_albfill(id)%buf(ic,jc) = val
                    wta(ic,jc,id)=1
                end if
            end do
            end do

            ! This will result in some re-writing; but we want it written sooner
            ! call chunker%write_chunks

        end do    ! id=1,NSTATS_FILLABLE

        ! Write out without any extra rewriting
        call chunker%write_chunks
    end do
    end do


    call chunker%close_chunks
    deallocate(grid)

end subroutine do_gridfill

end module b06_gridfill_mod

program Carrer_soilalbedo_gridfill
    use b06_gridfill_mod
    use ent_params_mod
    implicit none

    integer :: iband
    type(ReadWrites_t) :: rw

    call rw%init(THIS_OUTPUTS_DIR, 'B06_albmodis_gridfill', 20,20)

    call init_ent_labels
    do iband=1,NBANDS_MODIS
        call do_gridfill(rw, iband)
    end do
    call rw%write_mk

end program Carrer_soilalbedo_gridfill
