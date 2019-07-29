module ent_labels_mod

    use, intrinsic :: iso_fortran_env

implicit none

! Output files use the following indexing schemes:
!
! entmodisraw/
!     [Master ent20 scheme]
!
! pure/, trimmed/, etc.
! make_ent_gcm_subset(combine_crops_c3_c4=.true., split_bare_soil=.true.)
! (see chunkparams.f90 for definitions of these parameters)
!
!     01  ever_br_early
!     02  ever_br_late
!     03  ever_nd_early
!     04  ever_nd_late
!     05  cold_br_early
!     06  cold_br_late
!     07  drought_br
!     08  decid_nd
!     09  cold_shrub
!     10  arid_shrub
!     11  c3_grass_per
!     12  c4_grass
!     13  c3_grass_ann
!     14  c3_grass_arct
!     15  crops_herb
!     16  crops_woody
!     17  bare_bright
!     18  bare_dark
!
! If not combine_crops_c3_c4, then this would be modified as:
!     ...
!     15  crops_c3_herb
!     16  crops_c4_herb
!     17  crops_woody
!     18  bare_bright
!     19  bare_dark
!
! If not split_bare, this would be modified as:
!     ...
!     17  bare_sparse
!
! If neither parameter is set, it would look like:
!     ...
!     15  crops_c3_herb
!     16  crops_c4_herb
!     17  crops_woody
!     18  bare_sparse
!
! NOTE: Water is NOT part of these sub-sets

integer, parameter :: ENT_ABBREV_LEN = 20

! EntSet_t encapsulates a universe of cover types used as indices.
! This is done with respect to the "master" set of 20 Ent cover types
! (see below).  An EntSet_t instance is able to convert between its
! local (sub-) universe, and the full (master) universe.
!
! EntSet_t has the following interesting members:
!    esub%mvs: Convert from local to global indexing scheme
!    esub%svm: Convert from global to local indexing scheme
!    esub%layer_names():
!        Produce an array of names for each layer.
!
! NOTE: For cover types that exist in one but not the other scheme,
!       esub%mvs and esub%svm contain -1
!
! Suppose I have a variable lc(esub%ncover): Then to index the
! DRAUGHT_BR item in lc, one uses the expression:
!     lc(esub%svm(DRAUGHT_BR))
!
! EntSet_t instances are created in two ways:

!   1. The subroutine init_ent_labels creates the master EntSet_t
!      ent20, and ent19 (same but without water).  Note that ent20%svm
!      and ent20%mvs are the identity maping.  init_ent_labels() must
!      be called at the start of any program that needs to use
!      ent_labels_mod
!
!   2. The subroutine make_ent_gcm_subset() produces cover type
!       subsets.  Indices in those subsets will change depending on
!       user options (combine_crops_c3_c4, split_bare_soil).
type EntSet_t
    ! For derived schemes...
    integer, dimension(:), allocatable :: mvs    ! Convert this indexing scheme to master (ent20)
    integer, dimension(:), allocatable :: svm    ! Convert master (ent20) indices to this indexing scheme


    character*1, dimension(:), allocatable :: cttype
    character*(ENT_ABBREV_LEN), dimension(:), allocatable :: abbrev
    character*50, dimension(:), allocatable :: title

    ! Set by EntSet setup 
    integer :: npft, nnonveg
    ! Computed
    integer :: ncover
    character*(1) :: nonveg  ! M (M for MODIS, permanent snow/ice, bare/sparse, water), G (G for GISS, bare_bright, bare_dark)

    ! Set of master indexes to be remapped
    integer, dimension(:,:), allocatable :: remap
    integer :: nremap
contains
    procedure :: allocate => EntSet_allocate
    procedure :: layer_names
!    generic, public :: allocate => EntSet_allocate
    procedure :: add_covertype
    procedure :: add_remap
    procedure :: sub_covertype
end type

integer, parameter :: NMODIS2 = 28
!type(EntSet_t) :: modis28

! List of master Ent cover types
integer, parameter :: NENT20 = 20   ! Shortcut for ent20%ncover
type(EntSet_t) :: ent20      ! Master Ent categories
integer, parameter :: EVER_BR_EARLY = 1
integer, parameter :: EVER_BR_LATE  = 2
integer, parameter :: EVER_ND_EARLY = 3
integer, parameter :: EVER_ND_LATE  = 4
integer, parameter :: COLD_BR_EARLY = 5
integer, parameter :: COLD_BR_LATE  = 6
integer, parameter :: DROUGHT_BR    = 7
integer, parameter :: DECID_ND      = 8
integer, parameter :: COLD_SHRUB    = 9
integer, parameter :: ARID_SHRUB    = 10
integer, parameter :: C3_GRASS_PER  = 11
integer, parameter :: C4_GRASS      = 12
integer, parameter :: C3_GRASS_ANN  = 13
integer, parameter :: C3_GRASS_ARCT = 14
integer, parameter :: CROPS_C3_HERB = 15
integer, parameter :: CROPS_C4_HERB = 16
integer, parameter :: CROPS_WOODY   = 17
integer, parameter :: SNOW_ICE      = 18
integer, parameter :: BARE_SPARSE   = 19
integer, parameter :: CV_WATER      = 20

integer, parameter :: NENT19 = NENT20-1
type(EntSet_t) :: ent19



character(*), parameter :: TITLE_LC = 'Ent PFT 1 km land cover fraction'
character(*), parameter :: TITLE_LAI = &
    'Maximum annual LAI (m2/m2) 2004 downscaled from 1/12 degrees'
character(*), parameter :: TITLE_CHECKSUM = 'Checksum File'

real*4, parameter :: undef = -1.e30   ! Missing data in NetCDF

#if 1
integer, parameter :: nmonth = 12
character*3, parameter :: MONTH(nmonth) = &
     (/ &
     "Jan","Feb","Mar","Apr","May","Jun", &
     "Jul","Aug","Sep","Oct","Nov","Dec" &
     /)

#else
integer, parameter :: nmonth = 1
character*3, parameter :: MONTH(nmonth) = &
     (/ &
     "Jan" &
     /)

#endif

integer, parameter :: ndoy = 2
character*3, parameter :: DOY(ndoy) = &
     (/ &
     "017","201" &
     /)
      

CONTAINS

subroutine EntSet_allocate(ents, maxcover, ncover_master)
    class(EntSet_t) :: ents
    integer :: maxcover
    integer, OPTIONAL :: ncover_master
    integer :: i

    ents%npft = 0
    ents%nnonveg = 0
    ents%ncover = 0

    allocate(ents%cttype(maxcover))
    allocate(ents%abbrev(maxcover))
    allocate(ents%title(maxcover))

    ents%nremap = 0
    if (present(ncover_master)) then
        allocate(ents%mvs(maxcover))
        allocate(ents%svm(ncover_master))
        allocate(ents%remap(2,ncover_master))

        ! Initialize an identity remap
        do i=1,maxcover
            ents%mvs(i)=i
        end do
        do i=1,ncover_master
            ents%svm(i)=i
            ents%remap(1,i)=i
            ents%remap(2,i)=i
        end do
    end if
end subroutine EntSet_allocate

function itoa2(i) result(ret)
    integer,intent(IN) :: i
    character(2) :: ret

    write(ret,'(i2.2)') i
end function

function itoa3(i) result(ret)
    integer,intent(IN) :: i
    character(3) :: ret

    write(ret,'(i3.3)') i
end function

function itoa4(i) result(ret)
    integer,intent(IN) :: i
    character(4) :: ret

    write(ret,'(i4.4)') i
end function

function itoa(i) result(ret)
    integer,intent(IN) :: i
    character(2) :: ret

    write(ret,'(i2)') i
    ret = adjustl(ret)
end function

function layer_names(ents)
    class(EntSet_t) :: ents
    character*(ENT_ABBREV_LEN+3) :: layer_names(ents%ncover)
    ! ----- Locals
    integer :: k

    do k=1,ents%ncover
        layer_names(k) = trim(ents%abbrev(k))
    end do
end function layer_names

subroutine add_covertype(ents, cttype, abbrev,title)
    class(EntSet_t) :: ents
    character, intent(IN) :: cttype    ! 'v' for vegetation, 'n' for nonveg
    character*(*), intent(IN) :: abbrev
    character*(*), intent(IN) :: title

    ents%ncover = ents%ncover + 1
    if (ents%ncover > size(ents%abbrev,1)) then
        write(ERROR_UNIT,*) 'maxcover too small',size(ents%abbrev,1)
        stop -1
    end if

    ents%cttype(ents%ncover) = cttype
    ents%abbrev(ents%ncover) = abbrev
    ents%title(ents%ncover) = title

    if (cttype == 'v') then
        ents%npft = ents%npft + 1
    else if (cttype == 'n') then
        ents%nnonveg = ents%nnonveg + 1
    else
        write(ERROR_UNIT,*) 'add_covertype() Illegal cttype: ',cttype,trim(abbrev)
        stop
    end if
    ents%ncover = ents%npft + ents%nnonveg
    if (ents%nnonveg == 3) then
        ents%nonveg = 'M'
    else if (ents%nnonveg == 2) then
        ents%nonveg = 'G'
    else
        ents%nonveg = 'X'
    end if

end subroutine add_covertype

subroutine add_remap(ents, index20)
    class(EntSet_t) :: ents
    integer:: index20

    ents%svm(index20) = ents%ncover
    ents%mvs(ents%ncover) = index20

    ! Add a remap entry
    if (ents%nremap >= size(ents%remap,2)) then
        write(ERROR_UNIT,*) 'remap array too small',size(ents%remap,2)
        stop -1
    end if
    ents%nremap = ents%nremap + 1

    ents%remap(1,ents%nremap) = ents%ncover
    ents%remap(2,ents%nremap) = index20
end subroutine add_remap

subroutine sub_covertype(ents, ent20, index20)
    class(EntSet_t) :: ents
    type(EntSet_t), intent(IN) :: ent20
    integer, intent(IN) :: index20

    call ents%add_covertype(ent20%cttype(index20), ent20%abbrev(index20), ent20%title(index20))
    call ents%add_remap(index20)
end subroutine sub_covertype


subroutine init_ent_labels
    integer :: k

    call ent20%allocate(20,20)
    call ent20%add_covertype('v', 'ever_br_early ', 'evergreen broadleaf early successional      ') !  1
    call ent20%add_covertype('v', 'ever_br_late  ', 'evergreen broadleaf late successional       ') !  2
    call ent20%add_covertype('v', 'ever_nd_early ', 'evergreen needleleaf early successional     ') !  3
    call ent20%add_covertype('v', 'ever_nd_late  ', 'evergreen needleleaf late successional      ') !  4
    call ent20%add_covertype('v', 'cold_br_early ', 'cold deciduous broadleaf early successional ') !  5
    call ent20%add_covertype('v', 'cold_br_late  ', 'cold deciduous broadleaf late successional  ') !  6
    call ent20%add_covertype('v', 'drought_br    ', 'drought deciduous broadleaf                 ') !  7
    call ent20%add_covertype('v', 'decid_nd      ', 'deciduous needleleaf                        ') !  8
    call ent20%add_covertype('v', 'cold_shrub    ', 'cold adapted shrub                          ') !  9
    call ent20%add_covertype('v', 'arid_shrub    ', 'arid adapted shrub                         ')  ! 10
    call ent20%add_covertype('v', 'c3_grass_per  ', 'C3 grass perennial                         ')  ! 11
    call ent20%add_covertype('v', 'c4_grass      ', 'C4 grass                                   ')  ! 12
    call ent20%add_covertype('v', 'c3_grass_ann  ', 'C3 grass - annual                          ')  ! 13
    call ent20%add_covertype('v', 'c3_grass_arct ', 'arctic C3 grass                            ')  ! 14
    call ent20%add_covertype('v', 'crops_c3_herb ', 'crops C3 herb                              ')  ! 15
    call ent20%add_covertype('v', 'crops_c4_herb ', 'crops C4 herb                              ')  ! 16
    call ent20%add_covertype('v', 'crops_woody   ', 'crops woody                                ')  ! 17
    call ent20%add_covertype('n', 'snow_ice      ', 'Permanent snow/ice                         ')  ! 18
    call ent20%add_covertype('n', 'bare_sparse   ', 'Bare or sparsely vegetated, urban          ')  ! 19
    call ent20%add_covertype('n', 'water         ', 'water on land                              ')  ! 20


    ! Set up ent19 = ent20 without water
    ! (indices will be the same, but only if water is at the end)
    call ent19%allocate(NENT20-1, NENT20)
    do k=1,NENT20
        if (k == CV_WATER) cycle
        call ent19%sub_covertype(ent20, k)
    end do


end subroutine init_ent_labels

end module ent_labels_mod

! ==========================================================================

module gcm_labels_mod
    use ent_labels_mod
implicit none

    type, extends(EntSet_t) :: GcmEntSet_t
        ! Indices of non-standard cover types
        integer :: crops_herb
        integer :: last_pft   ! Last PLANT functional type (as opposed to cover type in general)
        integer :: bare_bright
        integer :: bare_dark
    contains
        procedure :: allocate => GcmEntSet_allocate
    !    generic, public :: allocate => GcmEntSet_allocate
    end type GcmEntSet_t

CONTAINS

subroutine GcmEntSet_allocate(ents, maxcover, ncover_master)
    class(GcmEntSet_t) :: ents
    integer :: maxcover
    integer, OPTIONAL :: ncover_master

    call EntSet_allocate(ents, maxcover, ncover_master)

    ents%crops_herb = -1
    ents%bare_bright = -1
    ents%bare_dark = -1
    ents%NONVEG = ''
end subroutine GcmEntSet_allocate

! Produces cover a cover type subset, depending on user options.
! Parameters:
!    combine_crops_c3_c4
!         Should C3 and C4 crop types be combined into one cover type?
!    split_bare_soil
!         Should the bare soil cover type be split into bright and dark sub-types?
! Returns:
!    The Ent covertype subset
function make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil) result(esub)
    logical, intent(IN) :: combine_crops_c3_c4
    logical, intent(IN) :: split_bare_soil
    type(GcmEntSet_t) :: esub

    ! ---------- Locals
    integer :: i

    ! ------- Set up the remap
    call esub%allocate(ent20%ncover, ent20%ncover)

    call esub%sub_covertype(ent20, EVER_BR_EARLY)  ! 1
    call esub%sub_covertype(ent20, EVER_BR_LATE)
    call esub%sub_covertype(ent20, EVER_ND_EARLY)
    call esub%sub_covertype(ent20, EVER_ND_LATE)
    call esub%sub_covertype(ent20, COLD_BR_EARLY)
    call esub%sub_covertype(ent20, COLD_BR_LATE)
    call esub%sub_covertype(ent20, DROUGHT_BR)
    call esub%sub_covertype(ent20, DECID_ND)
    call esub%sub_covertype(ent20, COLD_SHRUB)
    call esub%sub_covertype(ent20, ARID_SHRUB)
    call esub%sub_covertype(ent20, C3_GRASS_PER)
    call esub%sub_covertype(ent20, C4_GRASS)
    call esub%sub_covertype(ent20, C3_GRASS_ANN)
    call esub%sub_covertype(ent20, C3_GRASS_ARCT)  ! 14

    if (combine_crops_c3_c4) then
        !  15 - crops C3+C4 
        call esub%add_covertype('v', 'crops_herb', 'crops herbacious')    ! 15
        esub%crops_herb = esub%ncover
    else
        call esub%sub_covertype(ent20, CROPS_C3_HERB)        ! 15
        call esub%sub_covertype(ent20, CROPS_C4_HERB)        ! 16
    end if
    call esub%sub_covertype(ent20, CROPS_WOODY)          ! 16/17
    esub%last_pft = esub%ncover

    if (split_bare_soil) then
        ! Change bare_spares to bare-bright
        call esub%add_covertype('n', 'bare_bright', 'Bare or sparsely vegated, urban, bright')
        esub%bare_bright = esub%ncover
        call esub%add_remap(BARE_SPARSE)

        call esub%add_covertype('n', 'bare_dark', 'Bare or sparsely vegated, urban, dark')
        esub%bare_dark = esub%ncover
        esub%NONVEG = 'G'
    else
        call esub%sub_covertype(ent20, BARE_SPARSE)
    end if

    ! Use remap to pull out EntLabels

end function make_ent_gcm_subset

end module gcm_labels_mod
! ==========================================================================

module geom_mod

implicit none

! Various standard ModelE (and related) resolutions for lon/lat grids
integer, parameter :: X1km = 43200 !long at 1 km
integer, parameter :: Y1km = 21600 !lat at 1 km

integer, parameter :: IM1km = X1km !long at 1 km
integer, parameter :: JM1km = Y1km !lat at 1 km

integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees
integer, parameter :: IM1 = 360 !long at 1 degrees
integer, parameter :: JM1 = 180 !lat at 1 degrees
integer, parameter :: IM2 = 144 !long at 2.5 degrees
integer, parameter :: JM2 = 90 !lat at 2 degrees
integer, parameter :: IM4X5 = 72 !long at 5 degrees
integer, parameter :: JM4X5 = 46 !lat at 4 degrees

! Lo-res version (what we trim on)
integer, parameter :: IMLR=IMH
integer, parameter :: JMLR=JMH



end module geom_mod