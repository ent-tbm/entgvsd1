!  Program to assign 1kmx1km monthly LAI to EntPFTs
!------------------------------------------------------------------------

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

program lc_lai_monthly

use netcdf
use chunker_mod
use ent_labels_mod
use ent_params_mod
use assign_laimax_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(nmonth), io_lc(NENT20)
real*4, allocatable :: sum_lc(:,:)
! Output files
type(ChunkIO_t) :: io_laiout(NENT20,nmonth)
type(ChunkIO_t) :: io_lclai_checksum(nmonth)
type(ChunkIO_t) :: io_lclai_checksum_allmonths

type(FileInfo_t) :: info
integer :: imonth,k

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

!* Input file.

! ===================== Input Files
! Monthly LAI
do imonth = 1,nmonth
    call chunker%nc_open_input(io_lai(imonth), &
        INPUTS_URL, LAI_INPUTS_DIR, &
        'lai/BNU/monthly/'//sLAI_YEAR//'/', 'global_30s_'//sLAI_YEAR//'_'//MONTH(imonth)//'.nc', 'lai', 1)
enddo

! --- ENTPFTLC: Open outputs written by A00
call chunker%nc_open_set(ent20, io_lc, &
    LAI_SOURCE, 'M', 'lc', 2004, 'ent17', '1.1')

! =================== Output Files
do imonth=1,nmonth
    call chunker%nc_create_set( &
        ent20, io_laiout(:,imonth), lc_weights(io_lc, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'lai', 2004, 'ent17', '1.1', &
        doytype='month', idoy=imonth)

    call chunker%file_info(info, ent20, LAI_SOURCE, 'M', 'lclai', 2004, 'ent17', '1.1', &
        doytype='month', idoy=imonth, varsuffix='_checksum')
    call chunker%nc_create(io_lclai_checksum(imonth), &
        weighting(sum_lc,1d0,0d0), &
        info%dir, info%leaf, info%vname, &
        info%long_name, info%units)
enddo

call chunker%file_info(info, ent20, LAI_SOURCE, 'M', 'lclai', 2004, 'ent17', '1.1', &
    varsuffix='_allmonth_checksum')
call chunker%nc_create(io_lclai_checksum_allmonths, &
    weighting(sum_lc,1d0,0d0), &
    info%dir, info%leaf, info%vname, &
    info%long_name, info%units)


! Quit if we had any problems opening files
call chunker%nc_check('B10_lc_lai_monthly')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

! ====================== Done Opening Files

call assign_laimax(chunker, &
#ifdef ENTGVSD_DEBUG
    dbj0,dbj1, &
    dbi0,dbi1, &
#else
    1,chunker%nchunk(2), &
    1,chunker%nchunk(1), &
#endif
    io_lai, io_lc, io_laiout, &
    sum_lc=sum_lc, io_lclai_checksum=io_lclai_checksum, &
    io_lclai_checksum_alldoy = io_lclai_checksum_allmonths)

call chunker%close_chunks


end program lc_lai_monthly
