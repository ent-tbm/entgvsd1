! Assign 1kmx1km biomass to EntPFTs
! Author: Nancy Kiang, James Lui
!
!
!------------------------------------------------------------------------

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module assign_biomass_mod
  use chunker_mod
  use ent_params_mod
  use ent_labels_mod
  use hntr_mod
implicit none

#if (defined BIOMASS_SPAWN)
  integer, parameter :: n_biomass = 2
#elif (defined BIOMASS_GEDI) || (defined BIOMASS_XU)
  integer, parameter :: n_biomass = 1
#else
  integer, parameter :: n_biomass = -1
#endif
  contains

subroutine assign_biomass(chunker,&
   ! chunkerlr,  hntr_lr,&
    jc0,jc1,&
    ic0,ic1,&
    io_lc,&
    io_biomass,&
    io_biomassout,&
   ! io_biomassout_lr,&
    mywta,&
    checksum)

  type(Chunker_t) :: chunker!, chunkerlr
  integer :: jc0, jc1, ic0, ic1
  type(ChunkIO_t) :: io_lc(NENT20)
  type(ChunkIO_t) :: io_biomass(n_biomass) ! 1=aboveground 2=belowground
  type(ChunkIO_t) :: io_biomassout(NENT20,n_biomass)
  !type(ChunkIO_t) :: io_biomassout_lr(NENT20,2)
  real, dimension(:,:), optional :: mywta
  type(ChunkIO_t), optional :: checksum(n_biomass)
  !type(HntrCalc_t) :: hntr_lr

  ! local vars
  integer :: ichunk, jchunk, ic, jc, ii, jj, k, l
  real :: val, lcw ! biomass value, lc weight
  real, dimension(:,:), allocatable :: mywta1

  do jchunk = jc0, jc1
  do ichunk = ic0, ic1

    call chunker%move_to(ichunk, jchunk)
    !call chunkerlr%move_to(ichunk, jchunk)

    if (allocated(mywta1)) deallocate(mywta1)
    allocate(mywta1(chunker%chunk_size(1), chunker%chunk_size(2)))
    mywta1 = chunker%wta1

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

      do l=1,n_biomass
        val = io_biomass(l)%buf(ic,jc)
#if (defined BIOMASS_GEDI)
        val = val * 0.1 ! Mg Ha-1 to kg m-2
#endif

        if (present(checksum)) then
          checksum(l)%buf(ic,jc) = 0.
        endif
      do k=1,NENT20

        ! assign fillvalue if nanF
        ! assign 0 to mywta1 if nanF
        
        lcw = io_lc(k)%buf(ic,jc)

        ! filter out fillvalues from various sources
        if (val .eq. val .and. lcw .eq. lcw) then ! nanF does not equal itself, NaN = not land

          ! doesn't work the way you'd expect it to
          !if (val .le. 0d0 .or. lcw .le. 0d0 .or. lcw .eq. FillValue) then ! check if 0 or fillval
          if (lcw .le. 0. .or. lcw .eq. FillValue) then ! check if 0 or fillval
            io_biomassout(k,l)%buf(ic,jc) = 0. !FillValue
            mywta1(ic,jc) = 0.
          else
            io_biomassout(k,l)%buf(ic,jc) = max(val, 0.) ! filter negative fill values (-9999) 
          endif

          if (present(checksum)) then
            checksum(l)%buf(ic,jc) = checksum(l)%buf(ic,jc) + io_biomassout(k,l)%buf(ic,jc)
          endif

        else
          io_biomassout(k,l)%buf(ic,jc) = FillValue
          mywta1(ic,jc) = 0.
        end if 

      enddo ! k
      enddo ! l


!      write(*,*) shape(io_biomassout_lr(1,l)%buf), shape(io_biomassout(1,l)%buf), shape(mywta1)
!      write(*,*) io_biomassout_lr(1,l)%startB(2), io_biomassout_lr(1,l)%chunker%chunk_size(2)


    enddo ! ic
    enddo ! jc

    if (present(mywta)) then
      mywta = mywta1
    endif

!   do l=1,2
!   do k=1,NENT20
!
!     call hntr_lr%regrid4( &
!       io_biomassout_lr(k,l)%buf, &
!       io_biomassout(k,l)%buf, &
!       mywta1, 1d0, 0d0, &
!       io_biomassout_lr(k,l)%startB(2), &
!       io_biomassout_lr(k,l)%chunker%chunk_size(2))
!
!   enddo
!   enddo

    call chunker%write_chunks
    !call chunkerlr%write_chunks

  enddo ! ichunk
  enddo ! jchunk

  call chunker%close_chunks
  !call chunkerlr%close_chunks

end subroutine assign_biomass

end module assign_biomass_mod

program biomass

use netcdf
use chunker_mod
use ent_labels_mod
use ent_params_mod
use assign_biomass_mod
use hntr_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


type(Chunker_t) :: chunker!, chunkerlr
! Input files
type(ChunkIO_t), target :: io_biomass(n_biomass), io_lc(NENT20)
real*4, allocatable :: mywta(:,:)
! Output files
type(ChunkIO_t) :: io_biomassout(NENT20,n_biomass)
!type(ChunkIO_t) :: io_biomassout_lr(NENT20,2)
type(ChunkIO_t) :: io_biomass_checksum(n_biomass)
!type(ChunkIO_t) :: io_lclai_checksum(nmonth)
!type(ChunkIO_t) :: io_lclai_checksum_allmonths
!type(HntrSpec_t) :: spec_hr, spec_lr
!type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid

type(FileInfo_t) :: info, overmeta
integer :: imonth,k

#if (defined BIOMASS_SPAWN) || (defined BIOMASS_GEDI) ||\
    (defined BIOMASS_XU)
#else
write(*,*) "Missing argument for biomass: -b SPAWN, GEDI, or XU"
stop 1
#endif

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
#elif (defined BIOMASS_XU)
    overmeta%data_source = 'https://doi.org/10.1126/sciadv.abe9829'
    overmeta%global_data_source = &
        'biomass: Liang Xu et al. ,Changes in global terrestrial live biomass ' // &
            'over the 21st century.Sci. Adv.7,eabe9829(2021).DOI:10.1126/sciadv.abe9829'
#endif
    call init_ent_labels
    call chunker%init(IM1km, JM1km, IMH,JMH, 'forplot', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
!   call chunkerlr%init(IMLR,JMLR,IM2,JM2, 'forplot', 100, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    allocate(mywta(chunker%chunk_size(1), chunker%chunk_size(2)))

!hntr stuff
!  spec_hr = hntr_spec(chunker%chunk_size(1),chunker%ngrid(2),0d0,180d0*60d0 / chunker%ngrid(2))
!  spec_lr = hntr_spec(chunkerlr%chunk_size(1),chunkerlr%ngrid(2),0d0,180d0*60d0 / chunkerlr%ngrid(2))
!  hntr_lr = hntr_calc(spec_lr, spec_hr, FillValue8) ! datmis = FillValue


!allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

!* Input file.

! ===================== Input Files
#if (defined BIOMASS_SPAWN)
! above and belowground biomass
    call chunker%nc_open_input(io_biomass(1), &
        INPUTS_URL, INPUTS_DIR, &
        'biomass/', 'V1km_SpawnBiomass_netcdf4.nc', &
        'biomass_aboveground', 1)
    call chunker%nc_open_input(io_biomass(2), &
        INPUTS_URL, INPUTS_DIR, &
        'biomass/', 'V1km_SpawnBiomass_netcdf4.nc', &
        'biomass_belowground', 1)
#elif (defined BIOMASS_GEDI)
    call chunker%nc_open_input(io_biomass(1), &
        INPUTS_URL, INPUTS_DIR, &
        'biomass/', 'V1km_GEDI_aboveground_biomass_v2.nc', &
        'aboveground_biomass_density', 1)
#elif (defined BIOMASS_XU)
    call chunker%nc_open_input(io_biomass(1), &
        INPUTS_URL, INPUTS_DIR, &
        'biomass/', 'V1km_Xu2021_biomass_2004_v2.nc', &
        'carbon_density', 1)
#endif

! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')

! =================== Output Files
#if (defined BIOMASS_SPAWN)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        ent20, io_biomassout(:,1), lc_weights(io_lc, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_agb', 2010, 'ent17', '1.1.2', &
        create_lr=.true., overmeta=overmeta)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        ent20, io_biomassout(:,2), lc_weights(io_lc, 0d0, 1d0), &
        'Spawn2020', '', 'biomass_bgb', 2010, 'ent17', '1.1.2', &
        create_lr=.true., overmeta=overmeta)
#elif (defined BIOMASS_GEDI)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        ent20, io_biomassout(:,1), lc_weights(io_lc, 0d0, 1d0), &
        'GEDI', '', 'biomass_agb', 2022, 'ent17', '1.1.2', &
        create_lr=.true., overmeta=overmeta, heightsource='H2')
#elif (defined BIOMASS_XU)
    call chunker%nc_create_set( & ! lcweights are dummy!!
        ent20, io_biomassout(:,1), lc_weights(io_lc, 0d0, 1d0), &
        'Xu', '', 'biomass_agb', 2004, 'ent17', '1.1.2', &
        create_lr=.true., overmeta=overmeta, heightsource='X')
#endif

! =================== Regridded Files
!   call chunkerlr%nc_create_set( & ! lcweights are dummy!!
!       ent20, io_biomassout_lr(:,1), lc_weights(io_lc, 0d0, 1d0), &
!       LAI_SOURCE, '', 'biomass', 2010, 'biomass', '_aboveground', &
!       create_lr=.false.)
!   call chunkerlr%nc_create_set( & ! lcweights are dummy!!
!       ent20, io_biomassout_lr(:,2), lc_weights(io_lc, 0d0, 1d0), &
!       LAI_SOURCE, '', 'biomass', 2010, 'biomass', '_belowground', &
!       create_lr=.false.)

! =================== Checksum Files
#if (defined BIOMASS_SPAWN)
    call chunker%file_info(info, ent20, 'Spawn2020', '', 'biomass_agb', 2010, &
    'ent17', '1.1.2', varsuffix='_checksum')
    call chunker%nc_create(io_biomass_checksum(1), &
      weighting(mywta,1d0,0d0), &
      info%dir, info%leaf, info%vname, &
      info%long_name, info%units, global_data_source=overmeta%global_data_source)
      
    call chunker%file_info(info, ent20, 'Spawn2020', '', 'biomass_bgb', 2010, &
    'ent17', '1.1.2', varsuffix='_checksum')
    call chunker%nc_create(io_biomass_checksum(2), &
      weighting(mywta,1d0,0d0), &
      info%dir, info%leaf, info%vname, &
      info%long_name, info%units, global_data_source=overmeta%global_data_source)
#elif (defined BIOMASS_GEDI)
    call chunker%file_info(info, ent20, 'GEDI', '', 'biomass_agb', 2022, &
    'ent17', '1.1.2', varsuffix='_checksum', heightsource='H2')
    call chunker%nc_create(io_biomass_checksum(1), &
      weighting(mywta,1d0,0d0), &
      info%dir, info%leaf, info%vname, &
      info%long_name, info%units, global_data_source=overmeta%global_data_source)
#elif (defined BIOMASS_XU)
    call chunker%file_info(info, ent20, 'Xu', '', 'biomass_agb', 2004, &
    'ent17', '1.1.2', varsuffix='_checksum', heightsource='X')
    call chunker%nc_create(io_biomass_checksum(1), &
      weighting(mywta,1d0,0d0), &
      info%dir, info%leaf, info%vname, &
      info%long_name, info%units, global_data_source=overmeta%global_data_source)
#endif

! Quit if we had any problems opening files
call chunker%nc_check('B10b_lc_biomass_ann')
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

call assign_biomass(chunker,&
  ! chunkerlr, hntr_lr,&
#ifdef ENTGVSD_DEBUG
    dbj0,dbj1, &
    dbi0,dbi0, &
#else
    1,chunker%nchunk(2), &
    1,chunker%nchunk(1), &
#endif
    io_lc, io_biomass, io_biomassout,&
  ! io_biomassout_lr, &
    mywta=mywta, &
    checksum=io_biomass_checksum)

end program biomass
