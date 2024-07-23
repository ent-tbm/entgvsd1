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
  use gcm_labels_mod
  use ent_params_mod
  use ent_labels_mod
  use hntr_mod
implicit none
#if (defined BIOMASS_SPAWN)
  integer, parameter :: n_biomass = 2 ! 1=agb, 2=bgb
#elif (defined BIOMASS_GEDI)
  integer, parameter :: n_biomass = 1 ! 1=agb only
#else
  integer, parameter :: n_biomass = -1
#endif
real*8, parameter :: FillValue8=-1.d30
  contains

subroutine regrid_biomass(esub,&
    chunker, chunkerhxh, chunker2x2h,&
    hntr_hxh, hntr_2x2h,&
    jc0,jc1,&
    ic0,ic1,&
    io_lc_pure,&
    io_biomass1km,&
    io_biomass_hxh, io_biomass_2x2h,&
    rw)

  type(GcmEntSet_t), intent(IN), target :: esub
  type(Chunker_t) :: chunker, chunkerhxh, chunker2x2h
  integer :: jc0, jc1, ic0, ic1
  type(ChunkIO_t) :: io_lc_pure(esub%ncover)
  type(ChunkIO_t) :: io_biomass1km(esub%ncover,n_biomass)
  type(ChunkIO_t) :: io_biomass_hxh(esub%ncover,n_biomass)
  type(ChunkIO_t) :: io_biomass_2x2h(esub%ncover,n_biomass)
  type(HntrCalc_t) :: hntr_hxh, hntr_2x2h
  type(ReadWrites_t) :: rw

  ! local vars
  integer :: ichunk, jchunk, ic, jc, ii, jj, k, l
  real :: val, lcw ! biomass value, lc weight
  real, dimension(:,:), allocatable :: mywta1

  call chunker%nc_check(rw=rw)
  call chunkerhxh%nc_check(rw=rw)
  call chunker2x2h%nc_check(rw=rw)

  do jchunk = jc0, jc1
  do ichunk = ic0, ic1

    call chunker%move_to(ichunk, jchunk)
    call chunkerhxh%move_to(ichunk, jchunk)
    call chunker2x2h%move_to(ichunk, jchunk)

    do l=1,n_biomass
    do k=1,esub%ncover

      if (allocated(mywta1)) deallocate(mywta1)
      allocate(mywta1(chunker%chunk_size(1), chunker%chunk_size(2)))
      mywta1 = io_biomass1km(k,l)%buf

      do jc = 1,chunker%chunk_size(2)
      do ic = 1,chunker%chunk_size(1)

        val = io_biomass1km(k,l)%buf(ic,jc)
        if (val .eq. FillValue) then 
          mywta1(ic,jc) = 0d0
        end if 

      enddo ! ic
      enddo ! jc

      call hntr_hxh%regrid4( &
        io_biomass_hxh(k,l)%buf, &
        io_biomass1km(k,l)%buf, &
        mywta1, 1d0, 0d0, &
        io_biomass_hxh(k,l)%startB(2), &
        io_biomass_hxh(k,l)%chunker%chunk_size(2))

      where (io_biomass_hxh(k,l)%buf .lt. 0d0) io_biomass_hxh(k,l)%buf = FillValue

      call hntr_2x2h%regrid4( &
        io_biomass_2x2h(k,l)%buf, &
        io_biomass1km(k,l)%buf, &
        mywta1, 1d0, 0d0, &
        io_biomass_2x2h(k,l)%startB(2), &
        io_biomass_2x2h(k,l)%chunker%chunk_size(2))

      where (io_biomass_2x2h(k,l)%buf .lt. 0d0) io_biomass_2x2h(k,l)%buf = FillValue

    enddo
    enddo

    call chunker%write_chunks
    call chunkerhxh%write_chunks
    call chunker2x2h%write_chunks

  enddo ! ichunk
  enddo ! jchunk

  call chunker%close_chunks
  call chunkerhxh%close_chunks
  call chunker2x2h%close_chunks

end subroutine regrid_biomass

subroutine do_regrid(esub, rw)

use netcdf
use chunker_mod
use ent_labels_mod
use ent_params_mod
use hntr_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none

type(GcmEntSet_t), intent(IN), target :: esub

class(EntSet_t), pointer :: esub_p
type(Chunker_t) :: chunker, chunkerhxh, chunker2x2h
! Input files
type(ChunkIO_t), target :: io_biomass(n_biomass), io_lc_pure(esub%ncover)
real*4, allocatable :: mywta(:,:)
! Output files
type(ChunkIO_t) :: io_biomass1km(esub%ncover,n_biomass)
type(ChunkIO_t) :: io_biomass_hxh(esub%ncover,n_biomass)
type(ChunkIO_t) :: io_biomass_2x2h(esub%ncover,n_biomass)
type(ChunkIO_t) :: io_biomass_checksum(n_biomass)
!type(ChunkIO_t) :: io_lclai_checksum(nmonth)
!type(ChunkIO_t) :: io_lclai_checksum_allmonths
type(HntrSpec_t) :: spec_hr, spec_hxh, spec_2x2h
type(HntrCalc_t) :: hntr_hxh, hntr_2x2h    ! Preparation to regrid
type(ReadWrites_t) :: rw

type(FileInfo_t) :: info, overmeta
integer :: imonth,k

    esub_p => esub
    call clear_file_info(overmeta)
#if (defined BIOMASS_SPAWN)
    overmeta%data_source = 'https://doi.org/10.1038/s41597-020-0444-4'
    overmeta%global_data_source = &
        'biomass: Center for Sustainability and the Global Environment, ' // &
            'Nelson Institute for Environmental Studies, University ' // &
            'of Wisconsin-Madison, Biomass Carbon Density, 300m, ' // &
            '(Spawn et al. 2010, doi:10.1038/s41597-020-0444-4)'
#elif (defined BIOMASS_GEDI)
    overmeta%data_source = 'https://doi.org/10.3334/ORNLDAAC/2017'
    overmeta%global_data_source = &
        'biomass: Dubayah, R.O., J. Armston, S.P. Healey, Z. Yang, P.L. ' // &
            'Patterson, S. Saarela, G. Stahl, L. Duncanson, and J.R. Kellner. 2022.' // &
            'GEDI L4B Gridded Aboveground Biomass Density, Version 2. ORNL DAAC, Oak ' // &
            'Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/2017. NOTE: DRYBIOMASS'
#endif

    !call init_ent_labels
    call chunker%init(IM1km, JM1km, IMH,JMH, 'HXH', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    call chunkerhxh%init(IMLR,JMLR,IM2,JM2, '2X2H', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    call chunker2x2h%init(IM2,JM2,IM4X5,JM4X5, '4X5', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    allocate(mywta(chunker%chunk_size(1), chunker%chunk_size(2)))

!hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1),chunker%ngrid(2),0d0,180d0*60d0 / chunker%ngrid(2))
    spec_hxh = hntr_spec(chunkerhxh%chunk_size(1),chunkerhxh%ngrid(2),0d0,180d0*60d0 / chunkerhxh%ngrid(2))
    spec_2x2h = hntr_spec(chunker2x2h%chunk_size(1),chunker2x2h%ngrid(2),0d0,180d0*60d0 / chunker2x2h%ngrid(2))

    hntr_hxh = hntr_calc(spec_hxh, spec_hr, FillValue8) ! datmis = FillValue
    hntr_2x2h = hntr_calc(spec_2x2h, spec_hr, FillValue8) ! datmis = FillValue


!allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

!* Input file.

! ===================== Input Files
    call chunker%nc_open_set(esub_p, io_lc_pure, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1')

#if (defined BIOMASS_SPAWN)
    call chunker%nc_open_set(esub_p, io_biomass1km(:,1), &
        'Spawn2020', '', 'biomass_agb', 2010, 'pure', '1.1.2')

    call chunker%nc_open_set(esub_p, io_biomass1km(:,2), &
        'Spawn2020', '', 'biomass_bgb', 2010, 'pure', '1.1.2')
#elif (defined BIOMASS_GEDI)
    call chunker%nc_open_set(esub_p, io_biomass1km(:,1), &
        'GEDI', '', 'biomass_agb', 2022, 'pure', '1.1.2', heightsource='H2')
#endif

! =================== Regridded Files
#if (defined BIOMASS_SPAWN)
    call chunkerhxh%nc_create_set( & 
        esub_p, io_biomass_hxh(:,1), lc_weights(io_lc_pure, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_agb', 2010, 'pure', '1.1.2', &
        create_lr=.false., overmeta=overmeta)
    call chunkerhxh%nc_create_set( & 
        esub_p, io_biomass_hxh(:,2), lc_weights(io_lc_pure, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_bgb', 2010, 'pure', '1.1.2', &
        create_lr=.false., overmeta=overmeta)
    call chunker2x2h%nc_create_set( & 
        esub_p, io_biomass_2x2h(:,1), lc_weights(io_lc_pure, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_agb', 2010, 'pure', '1.1.2', &
        create_lr=.false., overmeta=overmeta)
    call chunker2x2h%nc_create_set( & 
        esub_p, io_biomass_2x2h(:,2), lc_weights(io_lc_pure, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_bgb', 2010, 'pure', '1.1.2', &
        create_lr=.false., overmeta=overmeta)
#elif (defined BIOMASS_GEDI)
    call chunkerhxh%nc_create_set( & 
        esub_p, io_biomass_hxh(:,1), lc_weights(io_lc_pure, 0d0, 1d0), &
        'GEDI', '', 'biomass_agb', 2022, 'pure', '1.1.2', &
        create_lr=.false., overmeta=overmeta, heightsource='H2')
    call chunker2x2h%nc_create_set( & 
        esub_p, io_biomass_2x2h(:,1), lc_weights(io_lc_pure, 0d0, 1d0), &
        'GEDI', '', 'biomass_agb', 2022, 'pure', '1.1.2', &
        create_lr=.false., overmeta=overmeta, heightsource='H2')
#endif

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

call regrid_biomass(esub,&
    chunker, chunkerhxh, chunker2x2h,&
    hntr_hxh, hntr_2x2h,&
#ifdef ENTGVSD_DEBUG
    dbj0,dbj1, &
    dbi0,dbi0, &
#else
    1,chunker%nchunk(2), &
    1,chunker%nchunk(1), &
#endif
    io_lc_pure,& 
    io_biomass1km, io_biomass_hxh, io_biomass_2x2h,&
    rw)

end subroutine do_regrid

end module regrid_biomass_mod

program regrid
    use regrid_biomass_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    type(GcmEntSet_t), target :: esub
    type(ReadWrites_t) :: rw

#if (defined BIOMASS_SPAWN) || (defined BIOMASS_GEDI)
#else
    write(*,*) "Missing argument for biomass: -b SPAWN or -b GEDI"
    stop 1
#endif
    call rw%init(THIS_OUTPUTS_DIR, "B14_regrid", 40,40)

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    call do_regrid(esub, rw)
    call rw%write_mk

end program regrid
