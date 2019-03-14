!  Program to assign 1kmx1km monthly LAI to EntPFTs
!------------------------------------------------------------------------

program lc_lai_monthly

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


type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(nmonth)
type(ChunkIO_t), target :: ioall_lc,io_lc(NENT20)
! Output files
type(ChunkIO_t) :: ioall_laiout(nmonth), io_laiout(NENT20,nmonth)
type(ChunkIO_t) :: io_checksum_lclai(nmonth)
type(ChunkIO_t) :: io_checksum_lclai_allmonths

integer :: imonth,k

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 320)

!* Input file.

! ===================== Input Files
! Monthly LAI
do imonth = 1,nmonth
    call chunker%nc_open_gz(io_lai(imonth), DATA_DIR, DATA_INPUT, &
        'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', 'lai', 1)
enddo

! ENTPFTLC: Outputs written by A00
!     ENTPFTLC: Outputs written by A00
call chunker%nc_open(ioall_lc, LC_LAI_ENT_DIR, &
    'pure/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
do k = 1,NENT20
    call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,k/))
enddo

! =================== Output Files
do imonth=1,nmonth
    call chunker%nc_create(ioall_laiout(imonth), &
        weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
        'pure/monthly/', 'entmm29_'//MONTH(imonth)//'_lai', 'lai', &
        'LAI output of A03', 'm2 m-2', 'LAI', &
        ent20%mvs, ent20%layer_names())
    do k=1,NENT20
        call chunker%nc_reuse_var(ioall_laiout(imonth), io_laiout(k,imonth), &
            (/1,1,k/), weighting(io_lc(k)%buf, 1d0,0d0))
    end do

    call chunker%nc_create(io_checksum_lclai(imonth), &
        weighting(chunker%wta1,1d0,0d0), &
        'pure/monthly/checksum/', &
        'entmm29_'//MONTH(imonth)//'_lai_checksum', &
        'Sum(LC*LAI) - LAI_orig == 0', 'm2 m-2', 'Sum of LC*LAI')
enddo

call chunker%nc_create(io_checksum_lclai_allmonths, &
    weighting(chunker%wta1,1d0,0d0), &
    'pure/monthly/checksum/', &
    'entmm29_allmonths_lai_checksum', &
    'Sum(LC*LAI) - LAI_orig == 0', 'm2 m-2', 'Sum of LC*LAI')




call chunker%nc_check('A03_lc_lai_monthly')
#ifdef JUST_DEPENDENCIES
stop 0
#endif


! ====================== Done Opening Files

! Quit if we had any problems opening files
call chunker%nc_check('A03_lc_lai_monthly')
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
    io_lai, io_lc, io_laiout, io_checksum_lclai, &
    io_checksum_lclai_alldoy = io_checksum_lclai_allmonths)

call chunker%close_chunks


end program lc_lai_monthly
