module chunkparams_mod
implicit none
    integer, parameter :: chunk_rank=2
!    integer, parameter :: nchunk(chunk_rank)=(/18,15/)   ! (lon, lat) (IM, JM) for chunks
!    integer, parameter :: nchunk(chunk_rank)=(/36,30/)   ! (lon, lat) (IM, JM) for chunks

    ! Combine C3 and C4 crops into one PFT, for ModelE
    logical, parameter :: combine_crops_c3_c4 = .true.
    ! Split the bare soil into dark and light, to get the right albedo
    logical, parameter :: split_bare_soil = .true.


CONTAINS

end module chunkparams_mod
