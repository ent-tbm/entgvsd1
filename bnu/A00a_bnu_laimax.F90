! Program computes annual maximum LAI from BNU monthly LAI.
! Authors:  Elizabeth Fischer, Carlo Montes
!------------------------------------------------------------------------

program bnu_laimax

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
type(ChunkIO_t), target :: io_lai(nallmonth)
! Output files
type(ChunkIO_t) :: io_laimax

integer :: imonth,ic,jc,ichunk,jchunk
real*4 :: lai,laimax

MAIN_PROGRAM_FILE='A00a_bnu_laimax'
call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 20, 5, 5, (/1,15/))

!* Input file.

! ===================== Input Files
! Monthly LAI
do imonth = 1,nallmonth
    call chunker%nc_open_gz(io_lai(imonth), DATA_DIR, DATA_INPUT, &
        'LAI/BNUMonthly/', 'global_30s_2004_'//ALLMONTH(imonth)//'.nc', 'lai', 1)
enddo

! =================== Output Files
call chunker%nc_create(io_laimax, weighting(chunker%wta1,1d0,0d0), &
    'bnu/', &
    'bnu_laimax', 'laimax', &
    'Maximum of montly LAI', 'm^2 m-2')

call chunker%nc_check(MAIN_PROGRAM_FILE)
#ifdef JUST_DEPENDENCIES
stop 0
#endif


! ====================== Done Opening Files

! Quit if we had any problems opening files
#ifdef JUST_DEPENDENCIES
stop 0
#endif

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
