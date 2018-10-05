program example_create_nc

use EntGVSD_netcdf_util
implicit none
    integer :: ncid,err

character(len=46), dimension(4), parameter :: layer_names = (/ &
     'evergreen broadleaf early succ               ', &
     'evergreen broadleaf late succ                ', &
     'evergreen needleleaf early succ              ', &
     'evergreen needleleaf late succ               ' &
    /)

integer :: layer_indices(4) = (/ 1,2,3,4 /)

    call my_nf90_create_Ent_single(43200,21600,4,'x.nc','myvar','My Variable','kg m-2', 'Title',ncid, &
        layer_indices, layer_names)

    err=nf90_close(ncid)
end program example_create_nc
