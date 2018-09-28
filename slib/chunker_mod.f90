module chunker_mod

    use entgvsd_config_mod
    use, intrinsic :: iso_fortran_env
    use netcdf
    use paths_mod

implicit none
private
    public :: Chunker_t,ChunkIO_t,chunk_rank,nchunk

    integer, parameter :: chunk_rank=2
    integer, parameter :: nchunk(chunk_rank)=(/18,15/)   ! (lon, lat) (IM, JM) for chunks

! ------ One file
type ChunkIO_t
    type(Chunker_t), pointer :: chunker
    character(200) :: leaf    ! Leaf name of file, for identification
    integer :: fileid, varid
    logical :: own_fileid   ! True if we own the fileid
    real*4, dimension(:,:), allocatable :: buf

end type ChunkIO_t

type ChunkIO_ptr
    type(ChunkIO_t), pointer :: ptr
end type ChunkIO_ptr

! ------ All the files, plus chunking info
type Chunker_t
    integer :: ngrid(chunk_rank)         ! Size of fine grid. (IM,JM) for grid
    integer :: chunk_size(chunk_rank)    ! Number of fine grid cells in each chunk (im,jm)
    integer :: cur(chunk_rank)     ! Index of current chunk

    ! Same of a low-res version of the chunker...
    integer :: ngrid_lr(chunk_rank)         ! Size of fine grid. (IM,JM) for grid
    integer :: chunk_size_lr(chunk_rank)    ! Number of fine grid cells in each chunk (im,jm)
    integer :: cur_lr(chunk_rank)     ! Index of current chunk


    ! Keep track if we've globally seen an error
    integer :: nerr = 0
    integer :: max_reads, max_writes
    integer :: nreads, nwrites
    type(ChunkIO_ptr), dimension(:), allocatable :: reads,writes
contains
    procedure :: init
    procedure :: write_chunks
    procedure :: read_chunks
    procedure :: move_to
    procedure :: close_chunks
    procedure :: nc_create
    procedure :: nc_open
    procedure :: nc_check
end type Chunker_t

CONTAINS

! @param im,jm Size of overall grid
! @param max_reads, max_writes Maximum number of read/write files
subroutine init(this, im, jm, im_lr, jm_lr, max_reads, max_writes)
    class(Chunker_t) :: this
    integer, intent(IN) :: im,jm
    integer, intent(IN) :: im_lr,jm_lr
    integer :: max_reads, max_writes
    ! ------ Locals
    integer :: i

    ! Set up chunk parameters
    this%ngrid(1) = im
    this%ngrid(2) = jm
    do i=1,chunk_rank
        this%chunk_size(i) = this%ngrid(i) / nchunk(i)
    end do

    ! Set up chunk parameters (lr)
    this%ngrid_lr(1) = im_lr
    this%ngrid_lr(2) = jm_lr
    do i=1,chunk_rank
        this%chunk_size_lr(i) = this%ngrid_lr(i) / nchunk(i)
    end do

    ! Allocate space to refer to managed chunk buffers
    this%max_reads = max_reads
    allocate(this%reads(this%max_reads))
    this%nreads = 0
    this%max_writes = max_writes
    allocate(this%writes(this%max_writes))
    this%nwrites = 0

    ! Invalid chunk index; means we have nothing to write on move_to()
    this%cur = (/0,0/)

end subroutine init


subroutine write_chunks(this)
    class(Chunker_t), target, intent(IN) :: this
    character :: rw
    ! ---------- Locals
    integer :: i
    integer :: startB(2)
    integer :: nerr
    integer :: err
    real*4, dimension(:,:), allocatable :: buf_lr  ! Buffer for lo-res version of chunk

    startB = (/ &
        (this%cur(1)-1) * this%chunk_size(1) + 1, &
        (this%cur(2)-1) * this%chunk_size(2) + 1/)
    write(*, '(A I3 I3)', advance="no") 'Writing Chunk',this%cur

    nerr = 0
    do i=1,this%nwrites
        err=NF90_PUT_VAR( &
           this%writes(i)%ptr%fileid, this%writes(i)%ptr%varid, &
           this%writes(i)%ptr%buf,startB,this%chunk_size)
        write(*,'(A)',advance="no") '.'
        if (err /= NF90_NOERR) then
            write(ERROR_UNIT,*) 'Error writing ',trim(this%writes(i)%ptr%leaf)
            nerr = nerr + 1
        end if

        ! Regrid chunk to low-res


    end do
    write(*,*)

    if (nerr > 0) then
        write(ERROR_UNIT,*) 'Errors encountered writing chunks, exiting'
        stop -1
    end if
end subroutine write_chunks

subroutine read_chunks(this)
    class(Chunker_t), target, intent(IN) :: this
    ! ---------- Locals
    integer :: i
    integer :: startB(2)
    integer :: nerr
    integer :: err

    startB = (/ &
        (this%cur(1)-1) * this%chunk_size(1) + 1, &
        (this%cur(2)-1) * this%chunk_size(2) + 1/)
    write(*, '(A I3 I3)', advance="no") 'Reading Chunk',this%cur

    nerr = 0
    do i=1,this%nreads
        err=NF90_GET_VAR( &
           this%reads(i)%ptr%fileid, this%reads(i)%ptr%varid, &
           this%reads(i)%ptr%buf,startB,this%chunk_size)
        write(*,'(A)',advance="no") '.'
        if (err /= NF90_NOERR) then
            write(ERROR_UNIT,*) 'Error reading ',trim(this%reads(i)%ptr%leaf)
            nerr = nerr + 1
        end if
    end do
    write(*,*)

    if (nerr > 0) then
        write(ERROR_UNIT,*) 'Errors encountered reading chunks, exiting'
        stop -1
    end if
end subroutine read_chunks


! Sets the NetCDF start and count arrays used to read/write a single chunk
subroutine move_to(this, ci, cj)
    class(Chunker_t), target :: this
    integer, intent(IN) :: ci,cj    ! Index of chunk
    ! ---------- Locals
    integer ::i
    integer :: startB(2), countB(2)

    this%cur(1) = ci
    this%cur(2) = cj

    call this%read_chunks
end subroutine move_to

subroutine close0(this, files, nfiles)
    type(Chunker_t) :: this
    type(ChunkIO_ptr), dimension(:) :: files
    integer :: nfiles
    ! ---------------- Locals
    integer :: err
    integer :: nerr
    integer :: i

    nerr = 0
    do i=1,nfiles
        if (files(i)%ptr%own_fileid) then
            err = nf90_close(files(i)%ptr%fileid)
            if (err /= NF90_NOERR) then
                write(ERROR_UNIT,*) 'Error closing ',trim(this%writes(i)%ptr%leaf)
                nerr = nerr + 1
            end if
        end if
    end do

    if (nerr > 0) then
        write(ERROR_UNIT,*) 'Errors encountered closing files, exiting'
        stop -1
    end if
end subroutine close0

subroutine close_chunks(this)
    class(Chunker_t) :: this
    ! ---------------- Locals
    call close0(this, this%writes, this%nwrites)
    call close0(this, this%reads, this%nreads)
end subroutine close_chunks

! -------------------------------------------------------------------

! Creates a NetCDF file for writing (by chunks)
subroutine nc_create(this, cio, dir, leaf, vname)
    class(Chunker_t) :: this
    type(ChunkIO_t), target :: cio
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    character*(*), intent(in) :: vname
    ! --------- Locals
    integer :: err
    character(2048) :: path_name
    character(4192) :: cmd
    character(10) :: lon_s, lat_s,fmt
    logical :: exist

    cio%leaf = leaf
    cio%own_fileid = .false.
    cio%fileid = -1


    path_name = LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'.nc'
    print *,'Opening ',trim(path_name)
    inquire(FILE=trim(path_name), EXIST=exist)
    if (.not.exist) then
        ! ------- Create the file from a template
        write(lon_s, '(I10)') nchunk(1)
        write(lat_s, '(I10)') nchunk(2)

        ! cdl = ENTGVSD_PROJECT_SOURCE_DIR//"/templates/"//cdl_leaf
        cmd = 'entgvsd_create_nc '// &
              LC_LAI_ENT_DIR//' '// &
              LC_LAI_ENT_ORIG//' '// &
              TEMPLATE_DIR//' '//trim(dir)//' '//trim(leaf)// &
              ' nchunkspec=lon/'//trim(adjustl(lon_s))//',lat/'//trim(adjustl(lat_s))
              
        print *,'----------------------------'
        print *,trim(cmd)
        call execute_command_line(cmd, .true., err)

        if (err /= 0) then
            write(ERROR_UNIT,*) 'Error running command to create', &
                LC_LAI_ENT_DIR//'/'//LC_LAI_ENT_ORIG//'/'//TEMPLATE_DIR// &
                '/'//trim(dir)//'/'//trim(leaf)
            this%nerr = this%nerr + 1
            return
        end if
    end if

    ! ------- Open the file we just created
    err = NF90_OPEN(path_name, NF90_WRITE, cio%fileid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening after create',trim(path_name),err
        this%nerr = this%nerr + 1
        return
    end if
    cio%own_fileid = .true.

    ! Allocate write buffer
    allocate(cio%buf(this%chunk_size(1), this%chunk_size(2)))

    ! Open NetCDF array
    err = NF90_INQ_VARID(cio%fileid,vname,cio%varid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error getting varid ',trim(leaf),err
        this%nerr = this%nerr + 1
        return
    end if

    ! Store pointer to this cio
    this%nwrites = this%nwrites + 1
    if (this%nwrites > this%max_writes) then
        write(ERROR_UNIT,*) 'Exceeded maximum number of write handles', cio%leaf
        stop -1
    end if
    this%writes(this%nwrites)%ptr => cio
end subroutine nc_create

subroutine nc_open(this, cio, iroot, oroot, dir, leaf, vname)
    class(Chunker_t) :: this
    type(ChunkIO_t), target :: cio
    character*(*), intent(in) :: iroot   ! Read (maybe compressed) file from here
    character*(*), intent(in) :: oroot   ! Write or link to here, then open
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    character*(*), intent(in) :: vname
    ! --------- Local vars
    character(4192) :: cmd
    integer :: err
    logical :: exist

    cio%fileid = -1

    ! -------- Decompress / link the file if it doesn't exist
    inquire(FILE=oroot//dir//leaf, EXIST=exist)
    if (.not.exist) then
        ! cdl = ENTGVSD_PROJECT_SOURCE_DIR//"/templates/"//cdl_leaf
        ! cmd = NCGEN//' -k nc4 -o '//ofname//' '//trim(cdl)
        cmd = 'entgvsd_link_input '//trim(iroot)//' '//trim(oroot)//' '//trim(dir)//' '//trim(leaf)
        print *,cmd

        call execute_command_line(cmd, .true., err)
        if (err /= 0) then
            write(ERROR_UNIT,*) 'Error running entgvsd_link_input',leaf,err
            this%nerr = this%nerr + 1
            return
        end if
    end if

    ! ------- Now open the file
    err = nf90_open(oroot//dir//leaf,NF90_NOWRITE,cio%fileid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening',trim(leaf),err
        this%nerr = this%nerr + 1
        return
    end if

    err = NF90_INQ_VARID(cio%fileid,vname,cio%varid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error getting varid ',trim(leaf),err
        this%nerr = this%nerr + 1
        return
    end if

    ! Store pointer to this cio
    this%nreads = this%nreads + 1
    if (this%nreads > this%max_reads) then
        write(ERROR_UNIT,*) 'Exceeded maximum number of read handles', cio%leaf
        stop -1
    end if
    this%reads(this%nreads)%ptr => cio

    ! Allocate write buffer
    allocate(cio%buf(this%chunk_size(1), this%chunk_size(2)))

end subroutine nc_open


subroutine nc_check(this)
    class(Chunker_t) :: this
    ! ------ Locals

    if (this%nerr > 0) then
        write(ERROR_UNIT,*) 'Previous errors encountered creating/opening files'
        stop -1
    end if
end subroutine nc_check

end module chunker_mod
! ======================================================
