program example_ioutil

use paths_mod
use chunker_mod
use ent_params_mod
implicit none

integer :: ncid
      integer, parameter :: IMH = 720 !long at 0.5 degrees
      integer, parameter :: JMH = 360 !lat at 0.5 degrees
    integer, parameter :: X1km = 43200 !long at 1 km
    integer, parameter :: Y1km = 21600 !lat at 1 km

    integer, parameter :: IM1km = X1km !long at 1 km
    integer, parameter :: JM1km = Y1km !lat at 1 km

    type(Chunker_t) :: chunker
    type(ChunkIO_t), target :: io_var, io_out0, io_outb, io_outc
    type(ChunkIO_t), dimension(:), allocatable :: io_out
    real*4 :: mean
    integer :: ic,jc
    integer :: ichunk,jchunk
    real*4, dimension(:,:), allocatable :: wta1_buf
    

character(len=46), dimension(4), parameter :: layer_names = (/ &
     'evergreen broadleaf early succ               ', &
     'evergreen broadleaf late succ                ', &
     'evergreen needleleaf early succ              ', &
     'evergreen needleleaf late succ               ' &
    /)

integer :: layer_indices(4) = (/ 1,2,3,4 /)
integer, parameter :: nlayers = 4
integer :: k

    call chunker%init(IM1km, JM1km, IMH*2, JMH*2, 17, 17)
    call chunker%nc_open_gz(io_var, &
        DATA_DIR, DATA_INPUT, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax', 1)

    ! ---- Open a single file with a single variable
    call chunker%nc_create(io_outb, &
        weighting(chunker%wta1,1d0,0d0), &
        'test_dir/', 'testb', 'testb', 'long_name', 'km', 'title')

    ! ---- Open a second variable in that file
    call chunker%nc_reuse_file(io_outb, io_outc, &
        'testc', 'testc', 'km', 'title', weighting(chunker%wta1,1d0,0d0))

    ! ---- Open a single file with many layers
    call chunker%nc_create(io_out0, &
        weighting(chunker%wta1,1d0,0d0), &
        'test_dir/', 'test', 'test', 'long_name', 'km', 'title', &
        layer_indices, layer_names)
    allocate(io_out(nlayers))
    do k=1,nlayers
        call chunker%nc_reuse_var(io_out0, io_out(k), (/1,1,k/), weighting(chunker%wta1,1d0,0d0))
    end do


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
            do k=1,nlayers
                io_out(k)%buf(ic,jc) = k*(ichunk*100+ic+jc)
                if (k==1) then
                    io_outb%buf(ic,jc) = io_out(k)%buf(ic,jc)
                    io_outc%buf(ic,jc) = io_out(k)%buf(ic,jc)
                end if
            end do
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
