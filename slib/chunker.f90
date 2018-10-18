module chunker_mod

    use iso_c_binding
    use entgvsd_config_mod
    use, intrinsic :: iso_fortran_env
    use netcdf
    use paths_mod
    use chunkparams_mod
    use hntr_mod
!    use entgvsd_netcdf_util

implicit none
private
    public :: Chunker_t,ChunkIO_t,weighting

! Describes how to weight a variable when regridding
type Weighting_t
    real*4, dimension(:,:), pointer :: buf
    real*8 :: MM,BB
end type Weighting_t

! ------ One file
! Read/Write Cases
!   1. 1 file, 1 single-layer var, 1 buffer:
!      => own_fileid=.true.
!      => allocated(buf)
!      => size(base,1) == 2
!   2. borrowed file, own single-layer var in another ChunkIO_t
!      => own_fileid=.false.
!      => size(base,1) == 2
!      => allocated(buf)
!   3. borrowed file, write in multi-layer var in another ChunkIO_t
!      => own_fileid=.false.
!      => size(base,1) == 3
!      => allocated(buf)
! 
!   4. 1 file, 1 multi-layer var, 0 buffers
!      (Never write directly from this ChunkIO_t)
!      => own_fileid=.true.
!      => /allocated(buf)
!   5. 1 file, many single-layer vars, 0 buffers
!      => own_fileid=.true.
!      => /allocated(buf)

type ChunkIO_t
    type(Chunker_t), pointer :: chunker
    character(1024) :: path
    character(200) :: leaf    ! Leaf name of file, for identification
    integer :: fileid, varid
    integer :: fileid_lr, varid_lr
    integer, dimension(:), allocatable :: base   ! Base in NetCDF var where to write our buf
    logical :: own_fileid,own_fileid_lr   ! True if we own the fileid
    real*4, dimension(:,:), allocatable :: buf    ! im,jm,nlayers (1 if just a single value)
    type(Weighting_t) :: wta
end type ChunkIO_t

type ChunkIO_ptr
    type(ChunkIO_t), pointer :: ptr
end type ChunkIO_ptr

! ------ All the files, plus chunking info
type Chunker_t
    integer :: ngrid(chunk_rank)         ! Size of fine grid. (IM,JM) for grid
    integer :: chunk_size(3)    ! Number of fine grid cells in each chunk (im,jm)
    integer :: cur(chunk_rank)     ! Index of current chunk

    ! Same of a low-res version of the chunker...
    integer :: ngrid_lr(chunk_rank)         ! Size of fine grid. (IM,JM) for grid
    integer :: chunk_size_lr(3)    ! Number of fine grid cells in each chunk (im,jm)
    integer :: cur_lr(chunk_rank)     ! Index of current chunk


    ! Keep track if we've globally seen an error
    integer :: nerr = 0
    integer :: max_reads, max_writes
    integer :: nreads, nwrites
    type(ChunkIO_ptr), dimension(:), allocatable :: reads,writes

    real*4, dimension(:,:), allocatable :: wta1  ! Default weights
contains
    procedure :: init
    procedure :: write_chunks
    procedure :: read_chunks
    procedure :: move_to
    procedure :: close_chunks
    procedure :: nc_reuse_var
    procedure :: nc_reuse_file
    procedure :: nc_create
    procedure :: nc_open_gz
    procedure :: nc_open
    procedure :: nc_check
end type Chunker_t

CONTAINS

subroutine handle_nf90_error(status, message)
    integer, intent(in) :: status
    character*(*) :: message
  
    if(status /= nf90_NoErr)  then
        write(*,*) status,message, ' ',trim(nf90_strerror(status))
        write(*,*) 'DONT FORGET TO DELETE THE OLD NETCDF FILES FIRST'
        STOP
    end if
end subroutine handle_nf90_error

subroutine calc_lon_lat_4X5(IM,JM,lon,lat)
    integer,intent(in) :: IM,JM
    real*4,intent(inout) :: lon(IM),lat(JM)
    !---Local---
    real*4 :: XX
    integer :: i,j

    !* Get and write lat, lon.
    ! Calc lats & Longs

    XX = -90.000000 + 1.00000 !  (4. * .5 = 1.0)
    lat(1) = XX
    XX = XX + 3.00000
    lat(2) = XX
    DO J=3,(JM-1)
       XX = XX + 4.00000
       lat(J) = XX
       write(0,*) 'Lat ',J,XX
    END DO
    lat(JM) = 90.000000 - 1.00000

    XX = -180.00000 + 5*.5 
    DO I=1,IM
       lon(I) = XX
       write(0,*) 'Lon ',I,XX
       XX = XX + 5.000000
    END DO
end subroutine calc_lon_lat_4X5


subroutine calc_lon_lat_2HX2(IM,JM,lon,lat)
    integer,intent(in) :: IM,JM
    real*4,intent(inout) :: lon(IM),lat(JM)
    !---Local---
    real*4 :: XX
    integer :: i,j,m

    !* Get and write lat, lon.
    ! Calc lats & Longs

    XX = -90.000000 + 1.00000   !  (2. * .5 = 1.0)
    DO J=1,JM
       lat(J) = XX
       write(0,*) 'Lat ',J,XX
       XX = XX + 2.00000
    END DO
    lat(JM) = 90.000000 - 1.00000

    XX = -178.75000
    DO I=1,IM
       lon(I) = XX
       write(0,*) 'Lon ',I,XX
       XX = XX + 2.500000
    END DO
end subroutine calc_lon_lat_2HX2

!************************************************************************
subroutine calc_lon_lat(IM,JM,lon,lat)
    integer,intent(in) :: IM,JM 
    !character(len=*) :: res
    real*4,intent(inout) :: lon(IM),lat(JM)
    !---Local---
    real*4 :: XX
    integer :: i,j
    !real*4 :: di, dj

    !*Assumes lon (-180 to 180) and lat (-90 to 90)

    if (IM.eq.72) then   !if (res.eq.'4X5') then
       call calc_lon_lat_4X5(IM,JM,lon,lat)
       !di = 5.0
       !dj = 4.0
       return
    else if (IM.eq.144) then !if (res.eq.'2HX2') then
       call calc_lon_lat_2HX2(IM,JM,lon,lat)
       !di = 2.5
       !dj = 2.0
       return
    else if (IM.eq.360) then  !if (res.eq.'1x1') then
       !call calc_lon_lat_1x1(longout,latout,lon,lat)
       !di = 1.0
       !dj = 1.0
    else if (IM.eq.720) then !if (res.eq.'HXH') then
       !di = 0.5
       !dj = 0.5
    elseif (IM.eq.1440) then !0.25 degree
       !res = 'QXQ'
    elseif (IM.eq.7200) then
      !res = '6km'
    elseif (IM.eq.43200) then
      !res = '1km'
    elseif (IM.eq.86400) then
      !res = '500m'
    else
       write(*,*) 'Fix resolution specs'
       return
    endif

    do i=1,IM
       !lon(i) = -180.0 + (i - 0.5)*di
       lon(i) = -180. + (360./IM)*i - (360./IM)*0.5
    enddo

    do j=1,JM
       !lat(j) = -90.0 + (j - 0.5)*dj
       lat(j) = -90. + (180./JM)*j - (180./JM)*0.5

    enddo
end subroutine calc_lon_lat     
!************************************************************************

! Constructor for Weighting_t
function weighting(buf,MM,BB) result(wta)
    real*4, dimension(:,:), target :: buf
    real*8 :: MM,BB
    type(Weighting_t) :: wta

    wta%buf => buf
    wta%MM = MM
    wta%BB = BB
end function weighting

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
    this%chunk_size(3)=1

    ! Allocate weight buffer
    if (allocated(this%wta1)) deallocate(this%wta1)
    allocate(this%wta1(this%chunk_size(1), this%chunk_size(2)))
    this%wta1 = 1.0

    ! Set up chunk parameters (lr)
    this%ngrid_lr(1) = im_lr
    this%ngrid_lr(2) = jm_lr
    do i=1,chunk_rank
        this%chunk_size_lr(i) = this%ngrid_lr(i) / nchunk(i)
    end do
    this%chunk_size_lr(3)=1

    ! Allocate space to refer to managed chunk buffers
    this%max_reads = max_reads
    if (allocated(this%reads)) deallocate(this%reads)
    allocate(this%reads(this%max_reads))
    this%nreads = 0
    this%max_writes = max_writes
    if (allocated(this%writes)) deallocate(this%writes)
    allocate(this%writes(this%max_writes))
    this%nwrites = 0

    ! Invalid chunk index; means we have nothing to write on move_to()
    this%cur = (/0,0/)

end subroutine init

! Helper function
subroutine setup_startB(this, cio, startB_hr, startB_lr)
    type(Chunker_t) :: this
    type(ChunkIO_t), target :: cio
    integer, dimension(:), allocatable :: startB_hr
    integer, dimension(:), allocatable, OPTIONAL :: startB_lr
    ! ------- Locals

    if (allocated(startB_hr)) deallocate(startB_hr)
    if (size(cio%base,1) == 3) then
        allocate(startB_hr(3))
        startB_hr(3) = cio%base(3)
    else
        allocate(startB_hr(2))
    end if
    startB_hr(1) = cio%base(1) + (this%cur(1)-1) * this%chunk_size(1)
    startB_hr(2) = cio%base(2) + (this%cur(2)-1) * this%chunk_size(2)


    if (present(startB_lr)) then
        if (allocated(startB_lr)) deallocate(startB_lr)
        if (size(cio%base,1) == 3) then
            allocate(startB_lr(3))
            startB_lr(3) = cio%base(3)
        else
            allocate(startB_lr(2))
        end if
        startB_lr(1) = cio%base(1) + (this%cur(1)-1) * this%chunk_size_lr(1)
        startB_lr(2) = cio%base(2) + (this%cur(2)-1) * this%chunk_size_lr(2)
    end if
end subroutine setup_startB


subroutine write_chunks(this)
    class(Chunker_t), target, intent(IN) :: this
    character :: rw
    ! ---------- Locals
    integer :: i,k
    integer, allocatable :: startB_hr(:), startB_lr(:)
    integer :: nerr
    integer :: err
    real*4, dimension(:,:), allocatable :: buf_lr  ! Buffer for lo-res version of chunk
    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr
    character(128) :: vname
    integer :: xtype,ndims
    integer :: dimids(3)
    integer :: len
    type(ChunkIO_t), pointer :: cio
    integer :: ic,jc
    real*4, parameter :: undef = -1.e30   ! Missing data in NetCDF
    real*4, dimension(:,:), pointer :: my_wta

    ! Hntr stuff
    spec_hr = hntr_spec(this%chunk_size(1), this%ngrid(2), 0d0, 180d0*60d0 / this%ngrid(2))
    spec_lr = hntr_spec(this%chunk_size_lr(1), this%ngrid_lr(2), 0d0, 180d0*60d0 / this%ngrid_lr(2))
    hntr = hntr_calc(spec_lr, spec_hr, 0d0)

    ! Buffers used to write
    allocate(buf_lr(this%chunk_size_lr(1), this%chunk_size_lr(2)))

    write(*, '(A I3 I3)', advance="no") 'Writing Chunk',this%cur

    nerr = 0
    do i=1,this%nwrites
        cio => this%writes(i)%ptr

        ! Skip ChunkIO_t that are just for keeping a file open
        if (.not.allocated(cio%buf)) cycle

        call setup_startB(this, cio, startB_hr, startB_lr)
        call setup_startB(this, cio, startB_hr, startB_lr)

        ! Store the hi-res chunk
        err=NF90_PUT_VAR( &
           cio%fileid, cio%varid, &
           cio%buf,startB_hr,this%chunk_size)
        if (err /= NF90_NOERR) then
            write(ERROR_UNIT,*) 'Error writing ',trim(cio%leaf),cio%varid,err
            nerr = nerr + 1
        end if
        err=nf90_sync(cio%fileid)

        ! Regrid chunk to low-res
        call hntr%regrid4(buf_lr, cio%buf, cio%wta%buf,cio%wta%MM,cio%wta%BB, &
            startB_lr(2), this%chunk_size_lr(2))

        ! Convert zeros to NaN in regridded data (for plotting)
        do jc=1,this%chunk_size_lr(2)
        do ic=1,this%chunk_size_lr(1)
            if (buf_lr(ic,jc) == 0) buf_lr(ic,jc) = undef
        end do
        end do
        
        ! Store the lo-res chunk
        err=NF90_PUT_VAR( &
           cio%fileid_lr, cio%varid_lr, &
           buf_lr, startB_lr, this%chunk_size_lr)
        if (err /= NF90_NOERR) then
            write(ERROR_UNIT,*) 'Error writing lo-res ',trim(cio%leaf),cio%varid_lr,err
            err= nf90_inquire_variable( &
                cio%fileid_lr, &
                cio%varid_lr, vname, xtype,ndims,dimids)
            ! More error output...
            write(ERROR_UNIT,*) trim(vname),xtype,ndims,dimids
            write(ERROR_UNIT,*) lbound(buf_lr),ubound(buf_lr),startB_lr,this%chunk_size_lr
            err = nf90_inquire_dimension(cio%fileid_lr,dimids(1),vname,len)
            write(ERROR_UNIT,*) 'dim1 ',trim(vname),len
            err = nf90_inquire_dimension(cio%fileid_lr,dimids(2),vname,len)
            write(ERROR_UNIT,*) 'dim2 ',trim(vname),len
            write(ERROR_UNIT,*) 
            nerr = nerr + 1
        end if

!        err=nf90_sync(cio%fileid_lr)

        ! Display progress
        write(*,'(A1)',advance="no") '.'

    end do

    do i=1,this%nwrites
        cio => this%writes(i)%ptr
        err=nf90_sync(cio%fileid_lr)
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
    integer, allocatable :: startB(:)
    integer :: nerr
    integer :: err
    type(ChunkIO_t), pointer :: cio

    write(*, '(A I3 I3)', advance="no") 'Reading Chunk',this%cur

    nerr = 0
    do i=1,this%nreads
        cio => this%reads(i)%ptr
        ! Skip ChunkIO_t that are just for keeping a file open
        if (.not.allocated(cio%buf)) cycle

        call setup_startB(this, this%reads(i)%ptr, startB)

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

    this%cur(1) = ci
    this%cur(2) = cj

    call this%read_chunks
end subroutine move_to

! Helper function
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
        if (files(i)%ptr%own_fileid_lr) then
            err = nf90_close(files(i)%ptr%fileid_lr)
            if (err /= NF90_NOERR) then
                write(ERROR_UNIT,*) 'Error closing lr ',trim(this%writes(i)%ptr%leaf)
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

! ===========================================================================
! Helper functions for creating

integer function my_nf90_inq_varid(ncidin,varname,varid)
  use netcdf
  integer,intent(in) :: ncidin
  character(len=*),intent(in) :: varname
  integer,intent(inout) :: varid
  !--- Local ----
  integer :: status
  
  status = nf90_inq_varid(ncidin,varname,varid)
  my_nf90_inq_varid = status
end function my_nf90_inq_varid

integer function my_nf90_inq_put_var_real32(ncidout, &
     varname,varid,varreal32)
use netcdf
integer,intent(in) :: ncidout
character(len=*),intent(in) :: varname
integer,intent(inout) :: varid
real*4, dimension(:) :: varreal32
!--- Local ---
integer :: status

status = my_nf90_inq_varid(ncidout,trim(varname),varid)
status = nf90_put_var(ncidout,varid, varreal32)
my_nf90_inq_put_var_real32 = status
end function my_nf90_inq_put_var_real32


! Opens/Creates just the file and dimensions
! Returns dim IDs
function my_nf90_create_ij(filename,IM,JM, ncid,dimids, layer_indices, layer_names) result(status)

    character(len=*), intent(in) :: filename
    integer,intent(in) :: IM,JM
    integer,intent(out) :: ncid
    integer,intent(inout) :: dimids(:)
    integer :: status    ! Return variable
    integer, dimension(:), OPTIONAL :: layer_indices
    character(len=*), dimension(:), OPTIONAL :: layer_names
    !--- Local ----
    integer :: idlon, idlat, varid
    integer :: idlayer_indices, idlayer_names
    real*4 :: lon(IM), lat(JM)
    !character*80 :: filenc
    character*10 :: res
    integer :: len1(1)
    integer :: strdimids(2)
    integer :: startS(2), countS(2)
    integer :: nlayers

    if (present(layer_indices)) then
        nlayers = size(layer_indices,1)
    else
        nlayers = 1
    end if


    if (IM.eq.72) then
      res = '72X46' !5X4
    elseif (IM.eq.144) then
      res = '2HX2' !'144X90'
    elseif (IM.eq.360) then
      res = '1X1'
    elseif (IM.eq.720) then
      res = 'HXH'
    elseif (IM.eq.1440) then
      res = 'QXQ'
    elseif (IM.eq.7200) then
      res = '6km'
    elseif (IM.eq.43200) then
      res = '1km'
    elseif (IM.eq.86400) then
      res = '500m'
    endif

!    inquire(FILE=trim(path_name), EXIST=exist)
!    if (exist) then

    status=nf90_open(filename, NF90_WRITE, ncid) !Get ncid if file exists
    if (status == NF90_NOERR) then
         dimids(1)=1
         dimids(2)=2
         if (nlayers > 1) dimids(3)=3    ! layers
         write(0,*) 'Netcdf file was previously created. ',trim(filename)
        return
    end if

    ! netcdf output file needs to be created
    write(0,*) 'Creating ',trim(filename)
    status=nf90_create(filename, NF90_NOCLOBBER+NF90_HDF5, ncid)
    if (status /= NF90_NOERR) return
    status=nf90_def_dim(ncid, 'lon', IM, dimids(1))
    if (status /= NF90_NOERR) return
    status=nf90_def_dim(ncid, 'lat', JM, dimids(2))
    if (status /= NF90_NOERR) return
    if (nlayers > 1) then
        status=nf90_def_dim(ncid, 'layers', nlayers, dimids(3))
        if (status /= NF90_NOERR) return
        strdimids(2)=dimids(3)

        status=nf90_def_dim(ncid, 'layer_name_len', len(layer_names(1)), strdimids(1))
        if (status /= NF90_NOERR) return

        status=nf90_def_var(ncid, 'layer_indices', NF90_INT, idlayer_indices)
        if (status /= NF90_NOERR) return

!        strdimids(1)=size(layer_names,1)
!        strdimids(2)=nlayers
        status=nf90_def_var(ncid, 'layer_names', NF90_CHAR, strdimids, idlayer_names)
        if (status /= NF90_NOERR) return

    end if


    ! Create lon
    status=nf90_def_var(ncid, 'lon', NF90_FLOAT, dimids(1), idlon)
    if (status /= NF90_NOERR) return
    status=nf90_def_var_deflate(ncid,idlon,1,1,4)
    if (status /= NF90_NOERR) return
    len1(1) = im
    status=nf90_def_var_chunking(ncid,idlon,NF90_CHUNKED, len1)
    if (status /= NF90_NOERR) return

    ! Create lat
    status=nf90_def_var(ncid, 'lat', NF90_FLOAT, dimids(2), idlat)
    if (status /= NF90_NOERR) return
    status=nf90_def_var_deflate(ncid,idlat,1,1,4)
    if (status /= NF90_NOERR) return
    len1(1) = jm
    status=nf90_def_var_chunking(ncid,idlat,NF90_CHUNKED, len1)
    if (status /= NF90_NOERR) return


    status=nf90_put_att(ncid, idlon, 'long_name', 'longitude')
    if (status /= NF90_NOERR) return
    status=nf90_put_att(ncid, idlat, 'long_name', 'latitude')
    if (status /= NF90_NOERR) return
    status=nf90_put_att(ncid, idlon, 'units', 'degrees_east')
    if (status /= NF90_NOERR) return
    status=nf90_put_att(ncid, idlat, 'units', 'degrees_north')
    if (status /= NF90_NOERR) return
    status=nf90_put_att(ncid, idlon, '_FillValue', -1.e30)
    if (status /= NF90_NOERR) return
    status=nf90_put_att(ncid, idlat, '_FillValue', -1.e30)

    status=nf90_enddef(ncid)
    if (status /= NF90_NOERR) return

    call calc_lon_lat(IM,JM,lon,lat)
    status=my_nf90_inq_put_var_real32(ncid,'lon',varid,lon)
    if (status /= NF90_NOERR) return
    status=my_nf90_inq_put_var_real32(ncid,'lat',varid,lat)


    if (nlayers>1) then
        status=nf90_put_var(ncid,idlayer_indices,layer_indices)
        if (status /= NF90_NOERR) return

        startS = (/ 1,1 /)
        countS = (/ len(layer_names(1)), size(layer_names,1) /)
        status=nf90_put_var(ncid,idlayer_names, layer_names,startS,countS)

        if (status /= NF90_NOERR) return
    end if

!    status=nf90_close(ncid)
!    if (status /= NF90_NOERR) return
end function my_nf90_create_ij

! Lookup dimension IDs and sizes by name
function nc_lookup_dims(ncid, dim_names, dim_ids, dim_sizes) result(err)
    integer, intent(IN) :: ncid
    character(len=*), dimension(:), intent(IN) :: dim_names
    integer, dimension(:), allocatable :: dim_ids
    integer, dimension(:), allocatable :: dim_sizes
    integer :: err

    integer :: k
    character*1024 :: name
    integer :: ndim

    ! Get dimension extents
    ndim = size(dim_names,1)
    allocate(dim_ids(ndim))
    allocate(dim_sizes(ndim))
    do k=1,ndim
        err = nf90_inq_dimid(ncid, trim(dim_names(k)), dim_ids(k))
        if (err /= NF90_NOERR) return
        err = nf90_inquire_dimension(ncid, dim_ids(k), name, dim_sizes(k))
        if (err /= NF90_NOERR) return
    end do
end function nc_lookup_dims


! Creates a variable in an already-open NetCDF file with dimensions defined
subroutine my_nf90_create_Ent_single(ncid, varid, nlayers, &
    varname,long_name,units,title)
!Creates a netcdf file for a single layer mapped Ent PFT cover variable.
    integer, intent(IN) :: ncid
    integer, intent(OUT) :: varid
    integer, intent(IN) :: nlayers

    character*(*), intent(in) :: varname
    character*(*), intent(in) :: long_name
    character*(*), intent(in) :: units,title

    !-- Local --
    integer, dimension(:), allocatable :: dimids
    integer :: k,n
    integer :: status
    integer :: ndim
    character*1024 :: text
    character(8) :: date
    integer, dimension(:), allocatable :: dim_sz

    ! Use existing variable if it exists

    status = nf90_redef(ncid)
    call handle_nf90_error(status, 'nc_redef '//trim(varname))

    ! Lookup dimensions by string name
    if (nlayers==1) then
        status=nc_lookup_dims(ncid, (/'lon', 'lat' /), dimids, dim_sz)
    else
        status=nc_lookup_dims(ncid, (/'lon   ', 'lat   ', 'layers'/), dimids, dim_sz)
    end if
    call handle_nf90_error(status, 'nc_lookup_dims '//trim(varname))


    status=nf90_inq_varid(ncid, varname, varid)
    if (status /= NF90_NOERR) then
        status=nf90_def_var(ncid, varname, NF90_FLOAT, dimids, varid)
        call handle_nf90_error(status, 'nf90_def_var '//trim(varname))
    end if

    status=nf90_def_var_deflate(ncid,varid,1,1,1)
    status=nf90_def_var_chunking(ncid,varid,NF90_CHUNKED, &
        make_chunksizes(dim_sz(1), dim_sz(2)))

    status=nf90_put_att(ncid,varid,"long_name", trim(long_name))
    call handle_nf90_error(status,  'nf90_put_att  long_name')
    status=nf90_put_att(ncid,varid,"units", units)
    status=nf90_put_att(ncid,varid,'_FillValue',-1.e30)
    call handle_nf90_error(status, 'nf90_put_att Ent vars '// &
         trim(long_name))


    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'long_name', long_name)
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'history','Sep 2018: E. Fischer,C. Montes, N.Y. Kiang')
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'title', title)
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'creator_name', 'NASA GISS')
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'creator_email', "elizabeth.fischer@columbia.edu,carlo.montes@nasa.gov,nancy.y.kiang@nasa.gov")
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'geospatial_lat_min', -90d0)
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'geospatial_lat_max', 90d0)
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'geospatial_lon_min', -180d0)
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'geospatial_lon_max', 180d0)
    status=nf90_put_att(ncid,NF90_GLOBAL, &
        'EntTBM', "Ent Terrestrial Biosphere Model")

    status=nf90_enddef(ncid)
    call handle_nf90_error(status, 'my_nf90_create_Ent')
end subroutine my_nf90_create_Ent_single

! Helper function
subroutine finish_cio_init(this, cio, rw, alloc)
    type(Chunker_t), target :: this
    type(ChunkIO_t), target :: cio
    character :: rw
    logical :: alloc

    cio%chunker => this

    ! Allocate write buffer (Hi res)
    if (alloc) then
        if (allocated(cio%buf)) deallocate(cio%buf)
        allocate(cio%buf(this%chunk_size(1), this%chunk_size(2)))
    end if

    ! Store pointer to this cio
    if (rw == 'w') then
        this%nwrites = this%nwrites + 1
        if (this%nwrites > this%max_writes) then
            write(ERROR_UNIT,*) 'Exceeded maximum number of write handles', cio%leaf
            stop -1
        end if
        this%writes(this%nwrites)%ptr => cio
    else
        ! Store pointer to this cio
        this%nreads = this%nreads + 1
        if (this%nreads > this%max_reads) then
            write(ERROR_UNIT,*) 'Exceeded maximum number of read handles', cio%leaf
            stop -1
        end if
        this%reads(this%nreads)%ptr => cio
    end if
end subroutine finish_cio_init

! ===========================================================================
! Opening for writing

! Point the ChunkIO_t to a single layer of an existing multi-layer file
subroutine nc_reuse_var(this, cio0, cio, base, rw, wta)
    class(Chunker_t) :: this
    type(ChunkIO_t), intent(IN) :: cio0
    type(ChunkIO_t), target :: cio
    integer, dimension(:), intent(IN) :: base
    character, intent(in) :: rw
    type(Weighting_t), intent(IN) :: wta

    ! Where to start other instance, in case it's a multi...
    cio%base = base
    cio%wta = wta

    ! ---- Re-use fileIDs from another instance
    cio%leaf = cio0%leaf

    cio%fileid = cio0%fileid
    cio%varid = cio0%varid
    cio%own_fileid = .false.

    cio%fileid_lr = cio0%fileid_lr
    cio%varid_lr = cio0%varid_lr
    cio%own_fileid_lr = .false.

    call finish_cio_init(this, cio, rw, .true.)
end subroutine nc_reuse_var

! Create a new single-layer variable in an existing file
subroutine nc_reuse_file(this, cio0, cio, &
    vname,long_name,units,title, wta)

    class(Chunker_t) :: this
    type(ChunkIO_t), intent(IN) :: cio0
    type(ChunkIO_t), target :: cio
    character*(*), intent(in) :: vname
    character*(*), intent(in) :: long_name,units,title
    type(Weighting_t), target :: wta

    cio%base = (/1,1/)
    cio%wta = wta
    cio%leaf = cio0%leaf

    ! ---- Re-use fileIDs from another instance
    cio%fileid = cio0%fileid
    cio%own_fileid = .false.

    cio%fileid_lr = cio0%fileid_lr
    cio%own_fileid_lr = .false.

    ! Create the netCDF variable
    call my_nf90_create_Ent_single(cio%fileid, cio%varid, 1, &
        vname, long_name, units, title)

    call my_nf90_create_Ent_single(cio%fileid_lr, cio%varid_lr, 1, &
        vname, long_name, units, title)

    call finish_cio_init(this, cio, 'w', .true.)
end subroutine nc_reuse_file

! Creates a NetCDF file with a single variable.
! The variable may be single layer or multi-layer.
subroutine nc_create(this, cio, &
wta, &   ! Weight by weights(i,j)*MM + BB
dir, leaf, &
vname,long_name,units,title, &
layer_indices, layer_names)

    class(Chunker_t) :: this
    type(ChunkIO_t), target :: cio
    type(Weighting_t), target :: wta
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    character*(*), intent(in), OPTIONAL :: vname
    character*(*), intent(in), OPTIONAL :: long_name,units,title
    integer, dimension(:), OPTIONAL :: layer_indices
    character(len=*), dimension(:), OPTIONAL :: layer_names

    ! --------- Locals
    integer :: err
    character(2048) :: path_name_lr
    character(4192) :: cmd
    character(10) :: lon_s, lat_s,fmt
    logical :: exist
    integer, dimension(:), allocatable :: dimids
    logical :: alloc_buf
    integer :: nlayer

    cio%wta = wta
    cio%leaf = leaf

    if (present(layer_indices)) then
        nlayer = size(layer_indices,1)
        allocate(dimids(3))
        cio%base = (/1,1,1/)    ! Not used
    else
        nlayer = 1
        allocate(dimids(2))
        cio%base = (/1,1/)
    end if

    ! ------ Create directory
    call execute_command_line('mkdir -p '//LC_LAI_ENT_DIR//trim(dir), &
        .true., err)

    ! ------ Open/Create hi-res file
    cio%path = LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'.nc'
    print *,'Writing ',trim(cio%path)
    if (present(layer_indices)) then
        err = my_nf90_create_ij(trim(cio%path), &
            this%ngrid(1), this%ngrid(2), &
            cio%fileid, dimids, &
            layer_indices, layer_names)
    else
        err = my_nf90_create_ij(trim(cio%path), &
            this%ngrid(1), this%ngrid(2), &
            cio%fileid, dimids)
    end if
    if (present(vname)) then
        call my_nf90_create_Ent_single(cio%fileid, cio%varid, nlayer, &
            vname, long_name, units, title)
    end if
    cio%own_fileid = .true.

    ! ---------- Open/Create lo-res file
    cio%path = LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'_lr.nc'
    print *,'Writing ',trim(cio%path)
    if (present(layer_indices)) then
        err = my_nf90_create_ij(trim(cio%path), &
            this%ngrid_lr(1), this%ngrid_lr(2), &
            cio%fileid_lr, dimids, &
            layer_indices, layer_names)
    else
        err = my_nf90_create_ij(trim(cio%path), &
            this%ngrid_lr(1), this%ngrid_lr(2), &
            cio%fileid_lr, dimids)
    end if

    if (present(vname)) then
        call my_nf90_create_Ent_single(cio%fileid_lr, cio%varid_lr, nlayer, &
            vname, long_name, units, title)
    end if
    cio%own_fileid_lr = .true.

    alloc_buf = (.not.present(layer_indices)).and.present(vname)
    call finish_cio_init(this, cio, 'w', alloc_buf)
end subroutine nc_create


! ===========================================================================
! Unzips a gzipped file in iroot, and makes it availabe in oroot
function gunzip_input_file(this, iroot, oroot, dir, leaf) result(err)
    type(Chunker_t) :: this
    character*(*), intent(in) :: iroot   ! Read (maybe compressed) file from here
    character*(*), intent(in) :: oroot   ! Write or link to here, then open
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    integer :: err
    ! -------- Locals
    character(4192) :: cmd
    logical :: exist

    ! -------- Decompress / link the file if it doesn't exist
    inquire(FILE=oroot//dir//leaf, EXIST=exist)
    if (exist) then
        err = 0   ! Success
        return
    end if

    ! cdl = ENTGVSD_PROJECT_SOURCE_DIR//"/templates/"//cdl_leaf
    ! cmd = NCGEN//' -k nc4 -o '//ofname//' '//trim(cdl)
    cmd = 'entgvsd_link_input '//trim(iroot)//' '//trim(oroot)//' '//trim(dir)//' '//trim(leaf)
    print *,trim(cmd)

    call execute_command_line(cmd, .true., err)
    if (err /= 0) then
        write(ERROR_UNIT,*) 'Error running entgvsd_link_input',leaf,err
        this%nerr = this%nerr + 1
        return 
    end if

end function gunzip_input_file

subroutine nc_open_gz(this, cio, iroot, oroot, dir, leaf, vname, k)
    class(Chunker_t) :: this
    type(ChunkIO_t), target :: cio
    character*(*), intent(in) :: iroot   ! Read (maybe compressed) file from here
    character*(*), intent(in) :: oroot   ! Write or link to here, then open
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    character*(*), intent(in) :: vname
    integer, intent(in) :: k
    ! --------- Local vars
    integer :: err

    cio%fileid = -1

    err = gunzip_input_file(this, iroot, oroot, dir, leaf)
    if (err /= 0) then
        write(ERROR_UNIT,*) 'Error unzipping ',trim(iroot)//trim(dir)//trim(leaf)
        this%nerr = this%nerr + 1
        return
    end if

    ! ------- Now open the file
    call this%nc_open(cio, oroot, dir, leaf, vname, k)

end subroutine nc_open_gz

subroutine nc_open(this, cio, oroot, dir, leaf, vname, k)
    class(Chunker_t) :: this
    type(ChunkIO_t), target :: cio
    character*(*), intent(in) :: oroot
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    character*(*), intent(in) :: vname
    integer, intent(in) :: k
    ! --------- Local vars
    integer :: err
    integer :: nlayers_dimid
    character(len=NF90_MAX_NAME) :: xname
    integer :: nlayers

    cio%fileid = -1

    print *,'Reading ', oroot//dir//leaf

    ! ------- Now open the file
    cio%path = oroot//dir//leaf
    err = nf90_open(trim(cio%path),NF90_NOWRITE,cio%fileid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening ',trim(leaf),err
        this%nerr = this%nerr + 1
        return
    end if
    cio%own_fileid = .true.

    ! Set cio%nlayers
    err = nf90_inq_dimid(cio%fileid, 'nlayers', nlayers_dimid)
    if (allocated(cio%base)) deallocate(cio%base)
    if (err == NF90_NOERR) then
        err = nf90_inquire_dimension(cio%fileid, nlayers_dimid, xname, nlayers)
        allocate(cio%base(3))
        cio%base(3) = k
    else
        nlayers = 1
        allocate(cio%base(2))
    end if
    cio%base(1) = 1
    cio%base(2) = 1

    err = NF90_INQ_VARID(cio%fileid,vname,cio%varid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error getting varid ',trim(leaf),trim(vname),err
        this%nerr = this%nerr + 1
        return
    end if

    call finish_cio_init(this, cio, 'r', (k>0))
end subroutine nc_open



subroutine nc_check(this, exename)
    class(Chunker_t) :: this
    character(len=*) :: exename
    ! ------ Locals
    integer :: i

    if (this%nerr > 0) then
        write(ERROR_UNIT,*) 'Previous errors encountered creating/opening files'
        stop -1
    end if

    open(17, FILE=trim(exename//'.mk'))
    write(17,'(AA)') trim(exename),'_INPUTS = \'
    do i=1,this%nreads
        if (this%reads(i)%ptr%own_fileid) then
            write(17,*) '   ',trim(this%reads(i)%ptr%path),' \'
        end if
    end do
    write(17,*)
    write(17,'(AA)') trim(exename),'_OUTPUTS = \'
    do i=1,this%nwrites
        if (this%writes(i)%ptr%own_fileid) then
            write(17,*) '   ',trim(this%writes(i)%ptr%path),' \'
        end if
    end do
    write(17,*)
    close(17)

end subroutine nc_check

end module chunker_mod
! ======================================================
