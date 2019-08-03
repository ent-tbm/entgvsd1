module assign_laimax_mod
    use chunker_mod
    use chunkparams_mod
    use ent_labels_mod

implicit none
CONTAINS

subroutine assign_laimax(chunker, &
    jc0,jc1, &
    ic0,ic1, &
    io_lai, io_lc, &                        ! INPUT
    io_laiout, io_checksum_lclai, &   ! OUTPUT
    io_checksum_lclai_alldoy, io_err)    ! Optional OUTPUT

    type(Chunker_t) :: chunker
    integer :: jc0,jc1, ic0,ic1
    ! Input files
    type(ChunkIO_t) :: io_lai(:)      ! io_lai(ndoy)
    type(ChunkIO_t) :: io_lc(NENT20)
    ! Output files
    type(ChunkIO_t) :: io_laiout(NENT20,size(io_lai,1))
    type(ChunkIO_t) :: io_checksum_lclai(size(io_lai,1))
    type(ChunkIO_t), OPTIONAL :: io_checksum_lclai_alldoy
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

            if (present(io_checksum_lclai_alldoy)) then
                io_checksum_lclai_alldoy%buf(ic,jc) = 0.0
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
                CHECKSUM = CHECKSUM - io_lai(idoy)%buf(ic,jc)
                io_checksum_lclai(idoy)%buf(ic,jc) = CHECKSUM

                if (present(io_checksum_lclai_alldoy)) then
                    io_checksum_lclai_alldoy%buf(ic,jc) = &
                        io_checksum_lclai_alldoy%buf(ic,jc) + CHECKSUM
                end if
            end do

        end do
        end do

        call chunker%write_chunks

    end do
    end do

end subroutine assign_laimax

end module assign_laimax_mod
