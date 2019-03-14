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
type(ChunkIO_t) :: ioall_lc,io_lc(NENT20)
! Output files
type(ChunkIO_t) :: ioall_laiout, io_laiout(NENT20,one)
type(ChunkIO_t) :: ioall_err, io_err(NENT20,one)
type(ChunkIO_t) :: io_checksum_lclai(one)

integer :: k
integer :: ichunk,jchunk, ic,jc, ii,jj

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

! ================= Input Files
!      LAI max
call chunker%nc_open_gz(io_laiin(1), DATA_DIR, DATA_INPUT, &
    'LAI/', 'global_30s_2004_max.nc', 'lai', 1)

! --- ENTPFTLC: Open outputs written by A00
call chunker%nc_open(ioall_lc, LC_LAI_ENT_DIR, &
    'pure/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
do k = 1,NENT20
    call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,k/))
enddo

! ================= Output Files


call chunker%nc_create(ioall_laiout, &
    weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
    'pure/annual/', 'entmm29_ann_laimax', 'lai', &
    'Ent maximum LAI for year', 'm^2 m-2', 'Leaf Area Index', &
    ent20%mvs, ent20%layer_names())
do k=1,NENT20
    call chunker%nc_reuse_var(ioall_laiout, io_laiout(k,1), &
        (/1,1,k/), weighting(io_lc(k)%buf, 1d0,0d0))
enddo


call chunker%nc_create(ioall_err, &
    weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
    'pure/annual/', 'entmm29_ann_laierr', 'lai', &
    'Ent maximum LAI for year', 'm^2 m-2', 'Leaf Area Index', &
    ent20%mvs, ent20%layer_names())
do k=1,NENT20
    call chunker%nc_reuse_var(ioall_err, io_err(k,1), &
        (/1,1,k/), weighting(io_lc(k)%buf, 1d0,0d0))
enddo

call chunker%nc_create(io_checksum_lclai(1),  weighting(chunker%wta1,1d0,0d0), &
    'EntMM_lc_laimax_1kmx1km/checksum_lclai/', &
    'lclai', &
    'Sum(LC*LAI) - LAI_orig == 0', 'm2 m-2', 'Sum of LC*LAI')

! ====================== Done Opening Files

! Quit if we had any problems opening files
call chunker%nc_check('A01_lc_laimax')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

call assign_laimax(chunker, &
#ifdef ENTGVSD_DEBUG
    nchunk(2)*3/4,nchunk(2)*3/4+1, &
    nchunk(1)*3/4,nchunk(1)*3/4+1, &
#else
    1,nchunk(2), &
    1,nchunk(1), &
#endif
    io_laiin, io_lc, io_laiout, io_checksum_lclai, &
    io_err=io_err)

call chunker%close_chunks

end program lc_laimax
