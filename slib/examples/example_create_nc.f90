program example_create_nc

use EntGVSD_netcdf_util
implicit none
    integer :: ncid,err

    call my_nf90_create_Ent_single(43200,21600,'x.nc','myvar','My Variable','kg m-2', ncid)

    err=nf90_close(ncid)
end program example_create_nc
