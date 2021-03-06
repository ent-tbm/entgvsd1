module assign_laimax_mod
    use chunker_mod
    use ent_params_mod
    use ent_labels_mod

implicit none
CONTAINS

! This subroutine is the common "guts" of
!    B08_lc_laimax, B09_lc_lai_doy, B10_lc_lai_monthly
! It is called by all three, after they open their (different) files
!
! Arguments:
!
! chunker:
!    Overall I/O control
! jc0,jc1,ic0,ic1:
!    Chunk indices, loop bounds
! ---------------- INPUT
! io_lai(ndoy):
!    File handles for LAImax (B08) or LAI (B09, B10)
!    Indexed by "ndoy", the number of separate "days" (or months) to process.
!    ndoy = (B08: 1 annual instance; B09: 2 specific days of year; B10: 12 months)
! io_lc(NENT20):
!    Fractional area for each land cover type (see ent_labels)
! --------------- OUTPUT
! io_laiout(NENT20, ndoy):
!    Per-LC LAI, for each of ndoy
! io_lclai_checksum(ndoy):
!    \sum_{land cover types} LC * LAI, for each of ndoy
! -------------- OUTPUT (optional)
! sum_lc:
!    A single buffer, fill with sum_{land cover types} LC
!    To be used elsewhere by calling function.
! io_lclai_checksum_alldoy:
!    \sum{doy} io_lclai_checksum, can be used for annual checksum over all months
!    NOTE: Not used
! io_err(NENT20, ndoy)   (used only for B08)
!    Indicates if LC>0 but LAIMmax == 0
subroutine assign_laimax(chunker, &
    jc0,jc1, &
    ic0,ic1, &
    io_lai, io_lc, &                        ! INPUT
    io_laiout, io_lclai_checksum, &   ! OUTPUT
    sum_lc, io_lclai_checksum_alldoy, io_err)    ! Optional OUTPUT

    type(Chunker_t) :: chunker
    integer :: jc0,jc1, ic0,ic1
    ! Input files
    type(ChunkIO_t) :: io_lai(:)      ! io_lai(ndoy)
    type(ChunkIO_t) :: io_lc(NENT20)
    ! Output files
    type(ChunkIO_t) :: io_laiout(NENT20,size(io_lai,1))
    type(ChunkIO_t) :: io_lclai_checksum(size(io_lai,1))
    real*4, dimension(:,:), OPTIONAL :: sum_lc  ! Weighting buffer
    type(ChunkIO_t), OPTIONAL :: io_lclai_checksum_alldoy
    type(ChunkIO_t), OPTIONAL :: io_err(NENT20,size(io_lai,1))

    ! -------------- Local Vars
    integer :: ichunk,jchunk,ic,jc,ii,jj
    integer :: k
    real*8 :: CHECKSUM
    integer :: idoy,ndoy

    ndoy = size(io_lai,1)

    ! Use these loop bounds for testing...
    ! it chooses a land area in Asia
    do jchunk = jc0,jc1
    do ichunk = ic0,ic1

        call chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            ! Compute overall NetCDF index of current cell
            ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
            jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

            ! Compute weighting for checksums
            if (present(sum_lc)) then
                sum_lc(ic,jc) = 0
                do k = 1,NENT20
                    sum_lc(ic,jc) = sum_lc(ic,jc) + io_lc(k)%buf(ic,jc)
                end do
            end if

            if (present(io_lclai_checksum_alldoy)) then
                io_lclai_checksum_alldoy%buf(ic,jc) = 0.0
            end if

            do idoy=1,ndoy

                ! Compute <original lai> - sum(LC*LAI), should equal 0
                CHECKSUM = 0d0
                do k = 1,NENT20

                    ! If this LC type does NOT participate in this 1km grid cell,
                    ! then LAI needs to be zero.
                    if ((io_lc(k)%buf(ic,jc) <= 0d0).or.(io_lc(k)%buf(ic,jc) == FillValue)) then
                        ! If lc==0, fix LAImax=0 (enforce the invariant)
                        io_laiout(k,idoy)%buf(ic,jc) = 0d0
                    else    ! non-zero LC type
                        ! Use the single LAI for all LC types in this gridcell.
                        io_laiout(k,idoy)%buf(ic,jc) = io_lai(idoy)%buf(ic,jc)

                        ! Problem if lc>0 but LAImax==0
                        if (present(io_err)) then
                            if (io_laiout(k,idoy)%buf(ic,jc)==0) then
                                io_err(k,idoy)%buf(ic,jc) = io_lc(k)%buf(ic,jc)
                            end if
                        end if
                    end if

                    CHECKSUM = CHECKSUM + &
                        io_lc(k)%buf(ic,jc) * io_laiout(k,idoy)%buf(ic,jc)
                end do
!                CHECKSUM = CHECKSUM - io_lai(idoy)%buf(ic,jc)
                io_lclai_checksum(idoy)%buf(ic,jc) = CHECKSUM

                if (present(io_lclai_checksum_alldoy)) then
                    io_lclai_checksum_alldoy%buf(ic,jc) = &
                        io_lclai_checksum_alldoy%buf(ic,jc) + CHECKSUM
                end if
            end do

        end do
        end do

        call chunker%write_chunks

    end do
    end do

end subroutine assign_laimax

end module assign_laimax_mod
