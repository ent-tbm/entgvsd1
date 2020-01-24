
!  Program to assign 1kmx1km BNU LAI of selected DOY to EntPFTs

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

!------------------------------------------------------------------------

program lc_lai_doy

use netcdf
use chunker_mod
use ent_labels_mod
use ent_params_mod
use assign_laimax_mod

implicit none

type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(ndoy)      ! Input...
type(ChunkIO_t), target :: io_lc(NENT20)
real*4, allocatable :: sum_lc(:,:)
! Output files
type(ChunkIO_t) :: io_laiout(NENT20,ndoy)
type(ChunkIO_t) :: io_lclai_checksum(ndoy)

type(FileInfo_t) :: info
integer :: idoy,k

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 100, 120, 10, outputs_dir=THIS_OUTPUTS_DIR)
allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

!* Input file.

! ================= Input Files
!     DOY LAI
do idoy=1,ndoy
    call chunker%nc_open_input(io_lai(idoy), &
        INPUTS_URL, LAI_INPUTS_DIR, &
        'lai/BNU/doy/2004/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai', 1)
enddo

!     ENTPFTLC: Outputs written by A00
call chunker%nc_open_set(ent20, io_lc, &
    LAI_SOURCE, 'M', 'lc', 2004, 'ent17', '1.1')

! ================= Output Files
do idoy = 1,ndoy
    call chunker%nc_create_set( &
        ent20, io_laiout(:,idoy), lc_weights(io_lc, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'lai', 2004, 'ent17', '1.1', &
        doytype='doy', idoy=idoy)

    call chunker%file_info(info, ent20, LAI_SOURCE, 'M', 'lclai', 2004, 'ent17', '1.1', &
        doytype='doy', idoy=idoy, varsuffix='_checksum')
    call chunker%nc_create(io_lclai_checksum(idoy), &
        weighting(sum_lc,1d0,0d0), &
        info%dir, info%leaf, info%vname, &
        info%long_name, info%units)
enddo

! ====================== Done Opening Files

! Quit if we had any problems opening files
call chunker%nc_check('B09_lc_lai_doy')
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
    io_lai, io_lc, io_laiout, &
    sum_lc=sum_lc, io_lclai_checksum=io_lclai_checksum)

call chunker%close_chunks




end program lc_lai_doy
