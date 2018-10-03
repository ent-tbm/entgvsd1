program example_ioutil

use paths_mod
use chunker_mod
use chunkparams_mod
implicit none

integer :: ncid
      integer, parameter :: IMH = 720 !long at 0.5 degrees
      integer, parameter :: JMH = 360 !lat at 0.5 degrees
    integer, parameter :: X1km = 43200 !long at 1 km
    integer, parameter :: Y1km = 21600 !lat at 1 km

    integer, parameter :: IM1km = X1km !long at 1 km
    integer, parameter :: JM1km = Y1km !lat at 1 km

    type(Chunker_t) :: chunker
    type(ChunkIO_t), target :: io_var, io_out
    real*4 :: mean
    integer :: ic,jc
    integer :: ichunk,jchunk

    call chunker%init(IM1km, JM1km, IMH*2, JMH*2, 17, 17)
    call chunker%nc_open(io_var, &
        DATA_DIR, DATA_INPUT, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax')
    call chunker%nc_create(io_out, &
        'test_dir/', 'test', 'test')


    call chunker%nc_check

!    do jchunk = 1,nchunk(2)
    do jchunk = 8,8
    do ichunk = 1,nchunk(1)
        print *,'chunk',ichunk,jchunk
        call chunker%move_to(ichunk,jchunk)

        mean = 0
        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)
            mean = mean + io_var%buf(ic,jc)
            io_out%buf(ic,jc) = ichunk*100+ic+jc
        end do
        end do
        call chunker%write_chunks

    end do
    end do
!    mean = mean / (chunker%chunk_size(1) * chunker%chunk_size(2))

!    print *,io_out%buf


    call chunker%close_chunks

    print *,'mean', mean

end program example_ioutil
