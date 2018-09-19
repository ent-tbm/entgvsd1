module chunker_mod

    use entgvsd_config_mod
    use, intrinsic :: iso_fortran_env
    use netcdf

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
    float*4, dimension(:,:), allocatable :: buf

end type ChunkIO_t

! ------ All the files, plus chunking info
type Chunker_t
    integer :: ngrid(chunk_rank)         ! Size of fine grid. (IM,JM) for grid
    integer :: chunk_size(chunk_rank)    ! Number of fine grid cells in each chunk (im,jm)
    integer :: cur(chunk_rank)     ! Index of current chunk

    ! Keep track if we've globally seen an error
    integer :: nerr = 0
    integer :: max_reads, max_writes
    integer :: nreads, nwrites
    type(ChunkerBuf_t), pointer, dimension(:), allocatable :: reads,writes
contains
    procedure :: init
    procedure :: write_chunks
    procedure :: read_chunks
    procedure :: move_to
    procedure :: close
    procedure :: nc_create
    procedure :: nc_open
    procedure :: nc_check
end type Chunker_t

CONTAINS

! @param im,jm Size of overall grid
! @param max_reads, max_writes Maximum number of read/write files
subroutine init(this, im, jm, max_reads, max_writes)
    class(Chunker_t) :: this
    integer, intent(IN) :: im,jm
    integer :: max_reads, max_writes
    ! ------ Locals
    integer :: i

    ! Set up chunk parameters
    this%ngrid(1) = im
    this%ngrid(2) = jm
    do i=1,chunk_rank
        this%chunk_size(i) = this%ngrid(i) / nchunk(i)
    end do

    ! Allocate space to refer to managed chunk buffers
    this%max_reads = max_reads
    allocate(reads(this%max_reads))
    nreads = 0
    this%max_writes = max_writes
    allocate(writes(this%max_writes))
    nwrites = 0

    ! Invalid chunk index; means we have nothing to write on move_to()
    cur_start = (/0,0/)

end subroutine init


subroutine write_chunks(this)
    class(Chunker_t), target, intent(IN) :: this
    character :: rw
    ! ---------- Locals
    integer :: i
    integer :: startB(2)
    integer :: nerr

    startB = (/ &
        (this%cur(1)-1) * this%chunk_size(1) + 1, &
        (this%cur(2)-1) * this%chunk_size(2) + 1/)

    nerr = 0
    do i=1,this%nwrites
        err=NF90_PUT_VAR( &
           this%writes(i)%fileid, this%writes(i)%varid, &
           this%writes(i)%buf,startB,this%chunk_size)
        if (err /= NF90_NOERR) then
            write(ERROR_UNIT,*) 'Error writing ',trim(this%writes(i)%leaf)
            nerr = nerr + 1
        end if
    end do

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

    startB = (/ &
        (this%cur(1)-1) * this%chunk_size(1) + 1, &
        (this%cur(2)-1) * this%chunk_size(2) + 1/)

    nerr = 0
    do i=1,this%nreads
        err=NF90_GET_VAR( &
           this%reads(i)%fileid, this%reads(i)%varid, &
           this%reads(i)%buf,startB,this%chunk_size)
        if (err /= NF90_NOERR) then
            write(ERROR_UNIT,*) 'Error reading ',trim(this%reads(i)%leaf)
            nerr = nerr + 1
        end if
    end do

    if (nerr > 0) then
        write(ERROR_UNIT,*) 'Errors encountered reading chunks, exiting'
        stop -1
    end if
end subroutine read_chunks


! Sets the NetCDF start and count arrays used to read/write a single chunk
subroutine move_to(this, ci, cj)
    class(Chunker_t), target, intent(IN) :: this
    integer, intent(IN) :: ci,cj    ! Index of chunk
    ! ---------- Locals
    integer ::i
    integer :: startB(2), countB(2)

    ! ---------- Write out old chunks
    if (cur_start(1) > 0) call this%read_chunks

    this%cur(1) = ci
    this%cur(2) = cj

    call this%write_chunks
end subroutine move_to

subroutine close0(this, files, nfiles)
    type(Chunker_t), target, intent(IN) :: this
    type(ChunkIO), pointer, dimension(:) :: files
    integer :: nfiles
    ! ---------------- Locals
    integer :: err
    integer :: nerr
    nerr = 0
    do i=1,nfiles
        if (files(i)%own_fileid) then
            err = nf90_close(files(i)%fileid)
            if (err /= NF90_NOERR) then
                write(ERROR_UNIT,*) 'Error closing ',trim(this%writes(i)%leaf)
                nerr = nerr + 1
            end if
        end if
    end do

    if (nerr > 0) then
        write(ERROR_UNIT,*) 'Errors encountered closing files, exiting'
        stop -1
    end if
end subroutine close0

subroutine close0(this)
    call close0(this, this%writes, this%nwrites)
    call close0(this this%reads, this%nreads)
end subroutine

! -------------------------------------------------------------------

! Creates a NetCDF file for writing (by chunks)
subroutine nc_create(this, cio, dir, leaf)
    class(Chunker_t) :: this
    type(ChunkIO_t) :: cio
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    ! --------- Locals
    integer :: err
    character(4192) :: cmd
    character(10) :: lon_s, lat_s,fmt

    cio%leaf = leaf
    cio%own_fileid = .false.
    cio%fileid = -1

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

    ! ------- Open the file we just created
    err = NF90_OPEN(LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'.nc', NF90_WRITE, cio%fileid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening after create',LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'.nc',err
        this%nerr = this%nerr + 1
        return
    end if
    cio%own_fileid = .true.

    ! Allocate write buffer
    allocate(cio%buf(im,jm))

    ! Store pointer to this cio
    this%nwrites = this%nwrites + 1
    if (this%nwrites > this%max_writes) then
        write(ERROR_UNIT,*) 'Exceeded maximum number of write handles', this%leaf
        stop -1
    end if
    this%writes(this%nwrites) => cio
end subroutine nc_create

subroutine nc_open(this, cio, iroot, oroot, dir, leaf)
    class(Chunker_t) :: this
    type(ChunkIO_t) :: cio
    character*(*), intent(in) :: iroot   ! Read (maybe compressed) file from here
    character*(*), intent(in) :: oroot   ! Write or link to here, then open
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    ! --------- Local vars
    character(4192) :: cmd
    integer :: err
    logical :: exist

    ncid = -1

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
            entgvsd_errs = entgvsd_errs + 1
            return
        end if
    end if

    ! ------- Now open the file
    err = nf90_open(oroot//dir//leaf,NF90_NOWRITE,ncid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening',trim(leaf),err
        entgvsd_errs = entgvsd_errs + 1
        return
    end if

    ! Store pointer to this cio
    this%nreads = this%nreads + 1
    if (this%nreads > this%max_reads) then
        write(ERROR_UNIT,*) 'Exceeded maximum number of read handles', this%leaf
        stop -1
    end if
    this%reads(this%nreads) => cio
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
module paths_mod


implicit none


! Hardcoded directories for input / output

! ---------- Original values (Carlos on discover)
!      --- Input (read-only)
!      CHARACTER (LEN=*), PARAMETER :: DATA_DIR='../../data/'
!      !CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR='../lc_lai_giss/'
!      CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR=
!     &     '../../lc_lai_for_1km1km/'
!      --- Output
!      CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR='../lc_lai_ent/'

! Values for Elizabeth's work on gibbs
!      --- Input (read-only)
CHARACTER (LEN=*), PARAMETER :: DATA_DIR= &
   '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/data/'

CHARACTER (LEN=*), PARAMETER :: DATA_INPUT= &
   '/home2/rpfische/git/entgvsd/inputs/data/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR= &
   '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/BNU/'// &
   'lc_lai_giss'

CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_INPUT= &
   '/home2/rpfische/git/entgvsd/inputs/lc_lai_giss/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR= &
    '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/'// &
    'lc_lai_for_1kmx1km/'
CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_INPUT= &
   '/home2/rpfische/git/entgvsd/inputs/lc_lai_giss/lc_lai_for_1kmx1km/'




!      --- Output
CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR= &
   '/home2/rpfische/git/entgvsd/bnu/lc_lai_ent/'

CHARACTER (LEN=*), PARAMETER :: TEMPLATE_DIR= &
   '/home2/rpfische/git/entgvsd/templates/'


! Place of output files originally written by carlos
! Used to create templates.

CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_ORIG= &
   '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/BNU/lc_lai_ent'

! ===========================================================

end module paths_mod
