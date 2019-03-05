
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

integer, parameter :: ndoy = 2
character*3, parameter :: DOY(ndoy) = &
     (/ &
     "017","201" &
     /)

integer :: f, p, z


type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(ndoy)
type(ChunkIO_t), target :: io_lc(NENT20)
! Output files
type(ChunkIO_t) :: ioall_out(ndoy), io_out(NENT20,ndoy)
type(ChunkIO_t) :: io_checksum_lclai(2)

!integer :: startA(1),startB(2),countA(1),countB(2)
!integer :: startX(1),startY(1),countX(1),countY(1)
!integer :: start3d(4),count3d(4)
integer :: dd(4),varidx(12),varidy(12)
integer :: layer_indices(20)
integer :: ichunk,jchunk,ic,jc
real*4, dimension(:,:), pointer :: wbuf

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

!* Input file.

! ================= Input Files
!     DOY LAI
do k = 1,2
    call chunker%nc_open_gz(io_lai(k), DATA_DIR, DATA_INPUT, &
        'LAI/', 'global_30s_2004_'//DOY(k)//'.nc', 'lai', 1)
enddo

!     ENTPFTLC: Outputs written by A00
do k = 1,NENT20
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lc.nc', &
        trim(ent20%abbrev(k)), 1)
end do

! ================= Output Files
do idoy = 1,ndoy
    call chunker%nc_create(ioall_out(idoy), &
        weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
        'nc/', 'EntMM_lc_lai_'//DOY(k)//'_1kmx1km', 'EntPFT', &
        'LAI output of A02', 'm2 m-2', 'LAI', &
        ent20%index, ent20%layer_names())
    do k=1,20
        call chunker%nc_reuse_var(ioall_out(k), io_out(k,idoy), &
            (/1,1,k/), 'w', weighting(io_lc(k)%buf, 1d0,0d0))
    end do

    call chunker%nc_create(io_checksum_lclai(k), &
        weighting(chunker%wta1,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/checksum_lclai/', &
        'lclai_'//DOY(k), &
        'Sum(LC*LAI) - LAI_orig == 0', 'm2 m-2', 'Sum of LC*LAI')
enddo

call chunker%nc_check('A02_lc_lai_doy')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

!-----------------------------------------------------------------
!     Loop for every grid point

! Use these loop bounds for testing...
! it chooses a land area in Asia
#ifdef ENTGVSD_DEBUG
do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
do jchunk = 1,nchunk(2)
do ichunk = 1,nchunk(1)
#endif

    call chunker%move_to(ichunk,jchunk)

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

         
!**   LAI data ------------------------------------------------------------------------------

        do k = 1,2
            LAI = io_lai(k)%buf(ic,jc)

            CHECKSUM = 0d0        
            do p = 1,NENT20


                if (p.eq.1) then
                    LCIN_in = io_water%buf(ic,jc)
                else
                    LCIN_in = io_lc(p-1)%buf(ic,jc)
                end if

                ! If this LC type does NOT participate in this 1km grid cell,
                ! then LAI needs to be zero.
                if ((LCIN_in <= 0).or.(LCIN_in == undef)) then
                    laiout = 0d0
                else   ! non-zero LC type
                    ! Use the single LAI for all LC types in this gridcell.
                    laiout = LAI
                end if

                ! A02 and A03 will probably have same error cells at
                ! A01.  No need to do discrepancy file for them.

                ! Compute checksum
                CHECKSUM = CHECKSUM + LCIN_in * laiout
                io_out(p,k)%buf(ic,jc) = laiout
            end do ! p=1,NENT20
            CHECKSUM = CHECKSUM - LAI

            io_checksum_lclai(k)%buf(ic,jc) = CHECKSUM
        end do   ! k=1,2
    end do    ! ic
    end do    !jc

    call chunker%write_chunks

enddo    ! ichunk
enddo    ! jchunk

call chunker%close_chunks

end program lc_lai_doy
