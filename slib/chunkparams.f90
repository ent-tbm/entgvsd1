module chunkparams_mod
implicit none
    integer, parameter :: chunk_rank=2
    integer, parameter :: nchunk(chunk_rank)=(/18,15/)   ! (lon, lat) (IM, JM) for chunks
!    integer, parameter :: nchunk(chunk_rank)=(/36,30/)   ! (lon, lat) (IM, JM) for chunks

    ! Combine C3 and C4 crops into one PFT, for ModelE
    logical, parameter :: combine_crops_c3_c4 = .true.
    ! Split the bare soil into dark and light, to get the right albedo
    logical, parameter :: split_bare_soil = .true.


CONTAINS

function make_chunksizes(im,jm) result(chunksizes)
    integer, intent(IN) :: im,jm
    integer, dimension(chunk_rank) :: chunksizes
    ! ----- Locals
    integer :: i
    integer :: dims(2)

    chunksizes(1) = im / nchunk(1)
    chunksizes(2) = jm / nchunk(2)
end function make_chunksizes

end module chunkparams_mod
