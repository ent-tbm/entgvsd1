module ent_labels_mod

implicit none

type EntLabel_t
    integer :: index
    character*3 :: file1
    character*14 :: file2
    character*50 :: title
end type EntLabel_t

type(EntLabel_t) :: ent19(19)   ! Basic Ent categories
type(EntLabel_t) :: ent20(20)   ! 0=water + Ent categories
type(EntLabel_t) :: ent18(18)   ! ent categories, crops merged, + bare bright/dark

integer, paramater :: ent20_water = 1

CONTAINS

function title_array(ents) result(titles)
    type(EntLabel_t), dimension(:) :: ents
    character*50 :: titles(size(ents,1))
    ! ------- Locals
    integer :: i

    do i=1,size(ents,1)
        titles(i) = ents(i)%title
    end do

end function title_array

function index_array(ents) result(indices)
    type(EntLabel_t), dimension(:) :: ents
    integer1 :: indices(size(ents,1))
    ! ------- Locals
    integer :: i

    do i=1,size(ents,1)
        indices(i) = ents(i)%index
    end do

end function index_array


function ent_label(index,file1,file2,title) result(label)
    integer, intent(IN) :: index
    character*(*), intent(IN) :: file1
    character*(*), intent(IN) :: file2
    character*(*), intent(IN) :: title
    type(EntLabel_t) :: label

    label%index = index
    label%file1 = file1
    label%file2 = file2
    label%title = title
end function ent_label



subroutine init_ent_labels

    integer :: i

    ent20 = (/ &
        ent_label(0, '   ', 'water         ', '0 - water                                       '), &
        ent_label(1, '01_', 'ever_br_early ', '1 - evergreen broadleaf early successional      '), &
        ent_label(2, '02_', 'ever_br_late  ', '2 - evergreen broadleaf late successional       '), &
        ent_label(3, '03_', 'ever_nd_early ', '3 - evergreen needleleaf early successional     '), &
        ent_label(4, '04_', 'ever_nd_late  ', '4 - evergreen needleleaf late successional      '), &
        ent_label(5, '05_', 'cold_br_early ', '5 - cold deciduous broadleaf early successional '), &
        ent_label(6, '06_', 'cold_br_late  ', '6 - cold deciduous broadleaf late successional  '), &
        ent_label(7, '07_', 'drought_br    ', '7 - drought deciduous broadleaf                 '), &
        ent_label(8, '08_', 'decid_nd      ', '8 - deciduous needleleaf                        '), &
        ent_label(9, '09_', 'cold_shrub    ', '9 - cold adapted shrub                          '), &
        ent_label(10,'10_', 'arid_shrub    ', '10 - arid adapted shrub                         '), &
        ent_label(11,'11_', 'c3_grass_per  ', '11 - C3 grass perennial                         '), &
        ent_label(12,'12_', 'c4_grass      ', '12 - C4 grass                                   '), &
        ent_label(13,'13_', 'c3_grass_ann  ', '13 - C3 grass - annual                          '), &
        ent_label(14,'14_', 'c3_grass_arct ', '14 - arctic C3 grass                            '), &
        ent_label(15,'15_', 'crops_c3_herb ', '15 - crops C3 herb                              '), &
        ent_label(16,'16_', 'crops_c4_herb ', '16 - crops C4 herb                              '), &
        ent_label(17,'17_', 'crops_woody   ', '17 - crops woody                                '), &
        ent_label(18,'18_', 'snow_ice      ', '18 - Permanent snow/ice                         '), &
        ent_label(19,'19_', 'bare_sparse   ', '19 - Bare or sparsely vegetated, urban          ') &
    /)

    do i=1,19
        ent19(i) = ent20(i+1)
    end do

    ! ---- Setup ent18 categories
    do i=1,14
        ent18(i) = ent19(i)
    end do
    ent18(15) = ent_label(15,'15_', 'crops_herb', '15 - Crops Herb')
    ent18(16) = ent_label(16,'16_', 'crops_woody', '17 - crops woody')
    ent18(17) = ent_label(17,'17_', 'bare_bright', '17 - Bright Bare Soil')
    ent18(18) = ent_label(18,'18_', 'bare_dark', '17 - Dark Bare Soil')

end subroutine init_ent_labels


end module ent_labels_mod





