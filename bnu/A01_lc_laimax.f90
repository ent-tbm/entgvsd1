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
type(ChunkIO_t) :: io_lai(one)
type(ChunkIO_t) :: io_lc(NENT20)
! Output files
type(ChunkIO_t) :: io_laiout(NENT20,one)
type(ChunkIO_t) :: io_err(NENT20,one)
type(ChunkIO_t) :: io_checksum_lclai(one)

integer :: k
integer :: ichunk,jchunk, ic,jc, ii,jj

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

! ================= Input Files
!      LAI max
call chunker%nc_open_gz(io_lai(1), DATA_DIR, DATA_INPUT, &
    'LAI/', 'global_30s_2004_max.nc', 'lai', 1)

! --- ENTPFTLC: Open outputs written by A00
do k = 1,NENT20
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lc.nc', &
        trim(ent20%abbrev(k)), 1)
enddo

! ================= Output Files
do k = 1,NENT20
    call chunker%nc_create(io_laiout(k,1),  weighting(io_lc(k)%buf,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lai', &
        trim(ent20%abbrev(k)), &
        ent20%title(k), 'm2 m-2', TITLE_LAI)

    call chunker%nc_create(io_err(k,1),  weighting(io_lc(k)%buf,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_err', &
        trim(ent20%abbrev(k)), &
        ent20%title(k), 'm2 m-2', TITLE_LAI)
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
    io_lai, io_lc, io_laiout, io_checksum_lclai, io_err)

call chunker%close_chunks

end program lc_laimax
