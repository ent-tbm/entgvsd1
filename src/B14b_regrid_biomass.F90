! Assign regrids biomass
! Author: Nancy Kiang, James Lui
!
!
!------------------------------------------------------------------------

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module regrid_biomass_mod
  use chunker_mod
  use ent_params_mod
  use ent_labels_mod
  use hntr_mod
implicit none
  real(kind=kind(1.0d0)), parameter :: FillValue8 = -1d30
  contains

subroutine regrid_biomass(chunker, chunkerlr, hntr_lr,&
    jc0,jc1,&
    ic0,ic1,&
!    io_lc,&
    io_biomass1km,&
    io_biomass_lr)

  type(Chunker_t) :: chunker, chunkerlr
  integer :: jc0, jc1, ic0, ic1
!  type(ChunkIO_t) :: io_lc(NENT20)
  type(ChunkIO_t) :: io_biomass1km(NENT20,2)
  type(ChunkIO_t) :: io_biomass_lr(NENT20,2)
  type(HntrCalc_t) :: hntr_lr

  ! local vars
  integer :: ichunk, jchunk, ic, jc, ii, jj, k, l
  real :: val, lcw ! biomass value, lc weight
  real, dimension(:,:), allocatable :: mywta1

  do jchunk = jc0, jc1
  do ichunk = ic0, ic1

    call chunker%move_to(ichunk, jchunk)
    call chunkerlr%move_to(ichunk, jchunk)

    if (allocated(mywta1)) deallocate(mywta1)
    allocate(mywta1(chunker%chunk_size(1), chunker%chunk_size(2)))
    mywta1 = chunker%wta1

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

      do l=1,2
      do k=1,NENT20

        val = io_biomass1km(k,l)%buf(ic,jc)
        if (val .eq. FillValue8) then 
          mywta1(ic,jc) = 0d0
        end if 

      enddo ! k
      enddo ! l

    enddo ! ic
    enddo ! jc

    do l=1,2
    do k=1,NENT20

      call hntr_lr%regrid4( &
        io_biomass_lr(k,l)%buf, &
        io_biomass1km(k,l)%buf, &
        mywta1, 1d0, 0d0, &
        io_biomass_lr(k,l)%startB(2), &
        io_biomass_lr(k,l)%chunker%chunk_size(2))

    enddo
    enddo

    call chunker%write_chunks
    call chunkerlr%write_chunks

  enddo ! ichunk
  enddo ! jchunk

  call chunker%close_chunks
  call chunkerlr%close_chunks

end subroutine regrid_biomass

end module regrid_biomass_mod

program regrid

use netcdf
use chunker_mod
use ent_labels_mod
use ent_params_mod
use regrid_biomass_mod
use hntr_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


type(Chunker_t) :: chunker, chunkerlr
! Input files
type(ChunkIO_t), target :: io_biomass(2), io_lc(NENT20)
real*4, allocatable :: mywta(:,:)
! Output files
type(ChunkIO_t) :: io_biomass1km(NENT20,2)
type(ChunkIO_t) :: io_biomass_lr(NENT20,2)
type(ChunkIO_t) :: io_biomass_checksum(2)
!type(ChunkIO_t) :: io_lclai_checksum(nmonth)
!type(ChunkIO_t) :: io_lclai_checksum_allmonths
type(HntrSpec_t) :: spec_hr, spec_lr
type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid

type(FileInfo_t) :: info
integer :: imonth,k

    call init_ent_labels
    call chunker%init(IM1km, JM1km, IMH,JMH, 'forplot', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    call chunkerlr%init(IMLR,JMLR,IM2,JM2, 'forplot', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    allocate(mywta(chunker%chunk_size(1), chunker%chunk_size(2)))

!hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1),chunker%ngrid(2),0d0,180d0*60d0 / chunker%ngrid(2))
    spec_lr = hntr_spec(chunkerlr%chunk_size(1),chunkerlr%ngrid(2),0d0,180d0*60d0 / chunkerlr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, FillValue8) ! datmis = FillValue


!allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

!* Input file.

! ===================== Input Files
    call chunker%nc_open_set(ent20, io_lc, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')

    call chunker%nc_open_set(ent20, io_biomass1km(:,1), &
        LAI_SOURCE, '', 'biomass', 2010, 'biomass', '1.0_aboveground')

    call chunker%nc_open_set(ent20, io_biomass1km(:,2), &
        LAI_SOURCE, '', 'biomass', 2010, 'biomass', '1.0_belowground')

! =================== Regridded Files
    call chunkerlr%nc_create_set( & ! lcweights are dummy!!
        ent20, io_biomass_lr(:,1), lc_weights(io_lc, 0d0, 1d0), &
        LAI_SOURCE, '', 'biomass', 2010, 'biomass', '1.0_aboveground', &
        create_lr=.false.)
    call chunkerlr%nc_create_set( & ! lcweights are dummy!!
        ent20, io_biomass_lr(:,2), lc_weights(io_lc, 0d0, 1d0), &
        LAI_SOURCE, '', 'biomass', 2010, 'biomass', '1.0_belowground', &
        create_lr=.false.)

! Quit if we had any problems opening files
call chunker%nc_check('B14b_regrid_biomass')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

! ====================== Done Opening Files

!call assign_laimax(chunker, &
!#ifdef ENTGVSD_DEBUG
!    dbj0,dbj1, &
!    dbi0,dbi1, &
!#else
!    1,chunker%nchunk(2), &
!    1,chunker%nchunk(1), &
!#endif
!    io_lai, io_lc, io_laiout, &
!    sum_lc=sum_lc, io_lclai_checksum=io_lclai_checksum, &
!    io_lclai_checksum_alldoy = io_lclai_checksum_allmonths)

call regrid_biomass(chunker, chunkerlr, hntr_lr,&
#ifdef ENTGVSD_DEBUG
    dbj0,dbj1, &
    dbi0,dbi0, &
#else
    1,chunker%nchunk(2), &
    1,chunker%nchunk(1), &
#endif
!    io_lc,& 
    io_biomass1km, io_biomass_lr)

call chunker%close_chunks


end program regrid
