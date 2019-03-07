
!  Program to assign 1kmx1km BNU LAI of selected DOY to EntPFTs

!------------------------------------------------------------------------

program lc_lai_doy

use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod
use geom_mod
use assign_laimax_mod

implicit none

type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(ndoy)
type(ChunkIO_t), target :: io_lc(NENT20)
! Output files
type(ChunkIO_t) :: ioall_laiout(ndoy), io_laiout(NENT20,ndoy)
type(ChunkIO_t) :: io_checksum_lclai(ndoy)

integer :: idoy,k

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

!* Input file.

! ================= Input Files
!     DOY LAI
do idoy=1,ndoy
    call chunker%nc_open_gz(io_lai(idoy), DATA_DIR, DATA_INPUT, &
        'LAI/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai', 1)
enddo

!     ENTPFTLC: Outputs written by A00
do k = 1,NENT20
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lc.nc', &
        trim(ent20%abbrev(k)), 1)
end do

! ================= Output Files
do idoy = 1,ndoy
    call chunker%nc_create(ioall_laiout(idoy), &
        weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
        'nc/', 'EntMM_lc_lai_'//DOY(idoy)//'_1kmx1km', 'EntPFT', &
        'LAI output of A02', 'm2 m-2', 'LAI', &
        ent20%mvs, ent20%layer_names())
    do k=1,NENT20
        call chunker%nc_reuse_var(ioall_laiout(idoy), io_laiout(k,idoy), &
            (/1,1,k/), 'w', weighting(io_lc(k)%buf, 1d0,0d0))
    end do

    call chunker%nc_create(io_checksum_lclai(idoy), &
        weighting(chunker%wta1,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/checksum_lclai/', &
        'lclai_'//DOY(idoy), &
        'Sum(LC*LAI) - LAI_orig == 0', 'm2 m-2', 'Sum of LC*LAI')
enddo

call chunker%nc_check('A02_lc_lai_doy')
#ifdef JUST_DEPENDENCIES
stop 0
#endif



! ====================== Done Opening Files

! Quit if we had any problems opening files
call chunker%nc_check('A02_lc_lai_doy')
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
    io_lai, io_lc, io_laiout, io_checksum_lclai)

call chunker%close_chunks




end program lc_lai_doy
