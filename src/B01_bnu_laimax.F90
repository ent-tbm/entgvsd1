! Computes annual maximum LAI from BNU monthly LAI.
! Authors:  Elizabeth Fischer
!------------------------------------------------------------------------

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

program bnu_laimax

use netcdf
use chunker_mod
use ent_params_mod
use ent_labels_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(nallmonth)
! Output files
type(ChunkIO_t) :: io_laimax

integer :: imonth,ic,jc,ichunk,jchunk
real*4 :: lai,laimax

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 20, 5, 5, (/1,15/), outputs_dir=THIS_OUTPUTS_DIR)

!* Input file.

! ===================== Input Files
! Monthly LAI
do imonth = 1,nallmonth
    call chunker%nc_open_input(io_lai(imonth), INPUTS_URL, LAI_INPUTS_DIR, &
        'lai/BNU/monthly/'//sLAI_YEAR//'/', 'global_30s_'//sLAI_YEAR//'_'//ALLMONTH(imonth)//'.nc', 'lai', 1)
enddo

! =================== Output Files
call chunker%nc_create(io_laimax, weighting(chunker%wta1,1d0,0d0), &
    'tmp/bnu/', 'bnu_laimax', &
    'laimax', 'Maximum of montly LAI', 'm^2 m-2')

call chunker%nc_check('B01_bnu_laimax')
#ifdef JUST_DEPENDENCIES
stop 0
#endif


! ====================== Done Opening Files

do jchunk = 1,chunker%nchunk(2)
do ichunk = 1,chunker%nchunk(1)

    call chunker%move_to(ichunk,jchunk)

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

        laimax = FillValue
        do imonth=1,nallmonth
            lai = io_lai(imonth)%buf(ic,jc)
            if (lai /= FillValue) then
                if (laimax == FillValue) then
                    laimax = lai
                else
                    laimax = max(laimax, lai)
                end if
            end if
        end do
        io_laimax%buf(ic,jc) = laimax

    end do
    end do
    call chunker%write_chunks
end do
end do


call chunker%close_chunks


end program bnu_laimax
