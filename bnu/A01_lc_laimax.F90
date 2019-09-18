!  Program to assign 1kmx1km BNU LAImax to EntPFTs lc

!------------------------------------------------------------------------
program lc_laimax

use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod
use geom_mod
use assign_laimax_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none

integer, parameter :: one = 1

type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t) :: io_laiin(one)
type(ChunkIO_t) :: io_lc(NENT20)
real*4, allocatable :: sum_lc(:,:)
! Output files
type(ChunkIO_t) :: io_laiout(NENT20,one)
type(ChunkIO_t) :: io_err(NENT20,one)
type(ChunkIO_t) :: io_lclai_checksum(one)

type(FileInfo_t) :: info
integer :: k
integer :: ichunk,jchunk, ic,jc, ii,jj

MAIN_PROGRAM_FILE='A01_lc_laimax'
call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 100, 120,10)
allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

! ================= Input Files
!      LAI max
call chunker%nc_open_gz(io_laiin(1), DATA_DIR, DATA_INPUT, &
    'LAI/', 'global_30s_2004_max.nc', 'lai', 1)

! --- ENTPFTLC: Open outputs written by A00
call chunker%nc_open_set(ent20, io_lc, &
    'B', 'M', 'lc', 2004, 'ent17', '1.1')

! ================= Output Files

call chunker%nc_create_set( &
    ent20, io_laiout(:,1), lc_weights(io_lc, 1d0, 0d0), &
    'B', 'M', 'laimax', 2004, 'ent17', '1.1')

call chunker%nc_create_set( &
    ent20, io_err(:,1), lc_weights(io_lc, 1d0, 0d0), &
    'B', 'M', 'laimax', 2004, 'ent17', '1.1', &
    varsuffix='_err')

call chunker%file_info(info, ent20, &
    'B', 'M', 'lclaimax', 2004, 'ent17', '1.1', &
    varsuffix = '_checksum')
call chunker%nc_create( &
    io_lclai_checksum(1),  weighting(sum_lc,1d0,0d0), &
    info%dir, info%leaf, info%vname, &
    'Sum(LC*LAI)', info%units)

! ====================== Done Opening Files

! Quit if we had any problems opening files
call chunker%nc_check(MAIN_PROGRAM_FILE)
#ifdef JUST_DEPENDENCIES
stop 0
#endif

call assign_laimax(chunker, &
#ifdef ENTGVSD_DEBUG
    dbj0,dbj1, &
    dbi0,dbi1, &
#else
    1,chunker%nchunk(2), &
    1,chunker%nchunk(1), &
#endif
    io_laiin, io_lc, io_laiout, &
    sum_lc=sum_lc, io_lclai_checksum=io_lclai_checksum, &
    io_err=io_err)

call chunker%close_chunks

end program lc_laimax
