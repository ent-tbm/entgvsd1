module EntGVSD_netcdf_util
!     Routines to create netcdf files for Ent Global Vegetation Structure Data sets.

use netcdf
use convertnc
use chunkparams_mod

implicit none
save

public my_nf90_create_Ent_single, my_nf90_create_Ent
public my_nf90_create_Ent_vartime
public N_COVERTYPES17, N_COVERTYPES16
public cmon, MONTH, GISSbands
public TITLE_LC, TITLE_LAI

character(*), parameter :: TITLE_LC = 'Ent PFT 1 km land cover fraction'
character(*), parameter :: TITLE_LAI = &
    'Maximum annual LAI (m2/m2) 2004 downscaled from 1/12 degrees'
character(*), parameter :: TITLE_CHECKSUM = 'Checksum File'

!Need to move this from convertnc_util.f
character*3, parameter :: cmon(12) = &
     (/ 'JAN','FEB','MAR','APR','MAY','JUN' &
       ,'JUL','AUG','SEP','OCT','NOV','DEC' /)

!Need to move this from trim_EntMM_monthly_05x05.f
character*3, parameter :: MONTH(12) = &
     (/ &
     "Jan","Feb","Mar","Apr","May","Jun", &
     "Jul","Aug","Sep","Oct","Nov","Dec" &
     /)

!********* Ent PFTs 17 + MODIS land cover names and titles **************

integer, parameter :: N_COVERTYPES17 = 20
!character(len=13), parameter :: ent_cover_names(N_COVERTYPES) = (/
character(len=13), parameter :: ent_names17(N_COVERTYPES17) = (/ &
     "ever_br_early", &
     "ever_br_late ", &
     "ever_nd_early", &
     "ever_nd_late ", &
     "cold_br_early", &
     "cold_br_late ", &
     "drought_br   ", &
     "decid_nd     ", &
     "cold_shrub   ", &
     "arid_shrub   ", &
     "c3_grass_per ", &
     "c4_grass     ", &
     "c3_grass_ann ", &
     "c3_grass_arct", &
     "crops_c3_herb", &
     "crops_c4_herb", &
     "crops_woody  ", &
     "permanent ice", &
     "bare_sparse  ", &
     "water        " &
     /)

      character*50, parameter :: Ent_title17(N_COVERTYPES17) = &
     (/ &
     '1 - evergreen broadleaf early succ               ', &
     '2 - evergreen broadleaf late succ                ', &
     '3 - evergreen needleleaf early succ              ', &
     '4 - evergreen needleleaf late succ               ', &
     '5 - cold deciduous broadleaf early succ          ', &
     '6 - cold deciduous broadleaf late succ           ', &
     '7 - drought deciduous broadleaf                  ', &
     '8 - deciduous needleleaf                         ', &
     '9 - cold adapted shrub                           ', &
     '10 - arid adapted shrub                          ', &
     '11 - C3 grass perennial                          ', &
     '12 - C4 grass                                    ', &
     '13 - C3 grass - annual                           ', &
     '14 - arctic C3 grass                             ', &
     '15 - crops C3 herb                               ', &
     '16 - crops C4 herb                               ', &
     '17 - crops woody                                 ', &
     '18 - permanent snow/ice                          ', &
     '19 - bare or sparsely vegetated                  ', &
     '20 - water                                       ' &
     /)


!**************** end Ent PFTs 17 + MODIS cover *************************      

!********* Ent PFTs 16 names and titles ***************************
!Need to move this from trim_EntMM_monthly_05x05.f
integer, parameter :: N_COVERTYPES16 = 18
!character(len=13), parameter :: ent_cover_names(N_COVERTYPES) = (/
character(len=13), parameter :: ent_names16(N_COVERTYPES16) = (/ &
     "ever_br_early", &
     "ever_br_late ", &
     "ever_nd_early", &
     "ever_nd_late ", &
     "cold_br_early", &
     "cold_br_late ", &
     "drought_br   ", &
     "decid_nd     ", &
     "cold_shrub   ", &
     "arid_shrub   ", &
     "c3_grass_per ", &
     "c4_grass     ", &
     "c3_grass_ann ", &
     "c3_grass_arct", &
     "crops_herb   ", &
     "crops_woody  ", &
     "bare_bright  ", &
     "bare_dark    " &
     /)

character*50, parameter :: Ent_title16(N_COVERTYPES16) = (/ &
     '1 - evergreen broadleaf early succ               ', &
     '2 - evergreen broadleaf late succ                ', &
     '3 - evergreen needleleaf early succ              ', &
     '4 - evergreen needleleaf late succ               ', &
     '5 - cold deciduous broadleaf early succ          ', &
     '6 - cold deciduous broadleaf late succ           ', &
     '7 - drought deciduous broadleaf                  ', &
     '8 - deciduous needleleaf                         ', &
     '9 - cold adapted shrub                           ', &
     '10 - arid adapted shrub                          ', &
     '11 - C3 grass perennial                          ', &
     '12 - C4 grass                                    ', &
     '13 - C3 grass - annual                           ', &
     '14 - arctic C3 grass                             ', &
     '15 - crops herb                                  ', &
     '16 - crops woody                                 ', &
     '17 - bright bare soil                            ', &
     '18 - dark bare soil                              ' &
     /)

!     *****end Ent PFTs 16  **********************************************

      
!     Need to move this from convertnc_util.f
character(len=19), dimension(6), parameter :: GISSbands = (/ &
     'VIS (300-770 nm)   ', &
     'NIR1 (770-860 nm)  ', &
     'NIR2 (860-1250 nm) ', &
     'NIR3 (1250-1500 nm)', &
     'NIR4 (1500-2200 nm)', &
     'NIR5 (2200-4000 nm)' /)

contains




!     Need to move from trim_EntMM_monthly_05x05.f
subroutine my_nf90_create_Ent_single(IM,JM,nlayers,file,varname,long_name,units,title,ncid, &
layer_indices,layer_names)
!Creates a netcdf file for a single layer mapped Ent PFT cover variable.
    integer, intent(in) :: IM, JM
    integer, intent(IN) :: nlayers
    character*(*) ::file
    character*(*) :: varname
    character*(*) :: long_name
    character*(*) :: units,title
    integer, intent(out) :: ncid
    integer, dimension(:), OPTIONAL :: layer_indices
    character(len=*), dimension(:), OPTIONAL :: layer_names
    !-- Local --
    integer :: k,n
    integer :: status, varid
    integer :: dimids(3) !dim(1)=lon, dim(2)=lat
    integer :: ndim
    character*1024 :: text
    character(8) :: date

    if (nlayers == 1) then
        ndim = 2
        status = my_nf90_create_ij(trim(file), IM,JM,nlayers,ncid, dimids)
    else
        ndim = 3
        status = my_nf90_create_ij(trim(file), IM,JM,nlayers,ncid, dimids, layer_indices, layer_names)
    end if
    call handle_nf90_error(status, 'nf90_cfeate_ij '//trim(varname))

    !Define global attributes - Customize local copies of this routine.
    !call my_nf90_defglobal(file)

    status=nf90_def_var(ncid, varname, NF90_FLOAT, dimids(1:ndim), varid)
    call handle_nf90_error(status, 'nf90_def_var '//trim(varname))
    status=nf90_def_var_deflate(ncid,varid,1,1,1)
    status=nf90_def_var_chunking(ncid,varid,NF90_CHUNKED, &
        make_chunksizes(im,jm))

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


!Need to move from trim_EntMM_monthly_05x05.f
subroutine my_nf90_create_Ent(IM,JM,file,ncov, &
     ent_names,Ent_title,units,ncid)
!Creates a netcdf file of mapped Ent PFT's, each variable by
!Ent PFT cover name, including bare soil.      
    integer, intent(in) :: IM, JM
    character*(*) ::file
    integer, intent(in) :: ncov
    character(len=*) :: ent_names(ncov)
    character(len=*) :: Ent_title(ncov)
    character*(*) :: units
    integer, intent(out) :: ncid
    !-- Local --
    integer :: k,n
    integer :: status, varid
    integer :: dimids(2) !dim(1)=lon, dim(2)=lat
    character*1024 :: text
    character(8) :: date

    n=ncov

    status = my_nf90_create_ij(trim(file), IM,JM,0,ncid, dimids)
    call handle_nf90_error(status,  'nf90_create '//file)       
    !call handle_nf90_error(status,  file//' '//dimlon//' '//dimlat)

    !Define global attributes - Customize local copies of this routine.
    !call my_nf90_defglobal(file)

    !Enter define mode
    status=nf90_open(trim(file), NF90_WRITE, ncid)
    call handle_nf90_error(status,  'nf90_open '//file)       
    status=nf90_redef(ncid)

    !Add variables and their attributes
    do k=1,n
       status=nf90_def_var(ncid, trim(ent_names(k)) &
            , NF90_FLOAT, dimids, varid)
       call handle_nf90_error(status, 'nf90_def_var '// &
            trim(ent_names(k)))
       status=nf90_put_att(ncid,varid,"long_name", trim(Ent_title(k)))
       call handle_nf90_error(status,  'nf90_put_att  long_name')
       status=nf90_put_att(ncid,varid,"short_name",trim(ent_names(k)))
       status=nf90_put_att(ncid,varid,"units", units)
       status=nf90_put_att(ncid,varid,'_FillValue', -1.e30)
       call handle_nf90_error(status, 'nf90_put_att Ent vars '// trim(ent_names(k)))
    end do

    status=nf90_enddef(ncid)
    call handle_nf90_error(status, 'nf90_enddef my_nf90_create_Ent')

    status=nf90_close(ncid)
    call handle_nf90_error(status, 'nf90_close my_nf90_create_Ent')
end subroutine my_nf90_create_Ent

subroutine my_nf90_create_Ent_lc_lai_max(IM,JM,file,ncov, ent_names,Ent_title,ncid)
!THIS IS NOT THE FORMAT FOR MODELE.  
!Creates a netcdf file of mapped Ent PFT lc and laimax, each variable prefixed by var type and suffix 
!Ent PFT name.      
    integer, intent(in) :: IM, JM
    character*(*) ::file
    integer, intent(in) :: ncov
    character(len=*) :: ent_names(ncov)
    character(len=*) :: Ent_title(ncov)
    integer, intent(out) :: ncid
    !-- Local --
    integer :: k,n, v
    integer :: status, varid
    integer :: dimids(2) !dim(1)=lon, dim(2)=lat
    character*1024 :: text
    character(8) :: date
    character*7 :: pre, prefix(2) = (/ 'lc_    ', 'laimax_' /)
    character*21  :: units, unitsarr(2) = (/ &
         'cover fraction       ', &
         'm^2 leaf / m^2 ground' /)
    n=ncov

    status = my_nf90_create_ij(trim(file), IM,JM,0,ncid, dimids(1:2))
    call handle_nf90_error(status,  'nf90_open '//file)       

    !call handle_nf90_error(status,  file//' '//dimlon//' '//dimlat)

    !Define global attributes - Customize local copies of this routine.
    !call my_nf90_defglobal(file)

    !Enter define mode
    status=nf90_open(trim(file), NF90_WRITE, ncid)
    call handle_nf90_error(status,  'nf90_open '//file)       
    status=nf90_redef(ncid)

    !Add variables and their attributes
    do v = 1,2
       pre = prefix(v)
       units = unitsarr(v)
    do k=1,n
       status=nf90_def_var(ncid, trim(pre)//trim(ent_names(k)) &
            , NF90_FLOAT, dimids, varid)
       call handle_nf90_error(status, 'nf90_def_var '// &
            trim(pre)//trim(ent_names(k)))
       status=nf90_put_att(ncid,varid,"long_name", &
            trim(pre)//trim(Ent_title(k)))
       call handle_nf90_error(status,  'nf90_put_att  long_name')
       status=nf90_put_att(ncid,varid,"short_name", &
            trim(pre)//trim(ent_names(k)))
       status=nf90_put_att(ncid,varid,"units", trim(units))
       status=nf90_put_att(ncid,varid,'_FillValue', -1.e30)
       call handle_nf90_error(status, 'nf90_put_att Ent vars '// &
            trim(pre)//trim(ent_names(k)))
    end do
    end do

    status=nf90_enddef(ncid)
    call handle_nf90_error(status, 'my_nf90_create_Ent')

end subroutine my_nf90_create_Ent_lc_lai_max



!     Need to move from trim_EntMM_monthly_05x05.f
subroutine my_nf90_create_Ent_vartime(IM,JM,file,ncov, &
     ent_names, Ent_title, &
     units,ncid)
!Creates a netcdf file of mapped Ent PFT's, each variable by
!Ent PFT cover name, including bare soil, in 3D arrays for time (monthly)      
!      use netcdf
!      use convertnc
implicit none
integer, intent(in) :: IM, JM
character*(*) ::file
integer, intent(in) :: ncov
character(len=*) :: ent_names(ncov)
character(len=*) :: Ent_title(ncov)
character*(*) :: units
integer, intent(out) :: ncid
!-- Local --
integer :: k,n
integer :: status, varid
character*1024 :: text
character(8) :: date
integer :: time(12)
integer :: dimids(3) !dim(1)=lon,dim(2)=lat,dim(3)=time

    n=ncov

    status = my_nf90_create_ij(trim(file), IM,JM,0,ncid, dimids(1:2))
    call handle_nf90_error(status, 'nf90_create_ij '//trim(file))

    !Define global attributes - Customize local copies of this routine.
    !call my_nf90_defglobal(file)

    !Add time dimension
    status=nf90_open(trim(file), NF90_WRITE, ncid)
    call handle_nf90_error(status, 'nf90_open '//trim(file))
    status=nf90_redef(ncid)
    call handle_nf90_error(status, 'nf90_redef '//trim(file))
    status=nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(3))
    call handle_nf90_error(status, 'nf90_def_dim'//' time')
     status=nf90_enddef(ncid)      

    !Add variables and their attributes
    status=nf90_def_var(ncid, 'time', NF90_INT, dimids(3), varid)
    call handle_nf90_error(status, 'nf90_def_var '//'time')
    status=nf90_put_att(ncid,varid,"units", 'Month Sequence of Climatology')

    do k=1,n
       status=nf90_def_var(ncid, trim(ent_names(k)) &
            , NF90_FLOAT, dimids, varid)
       call handle_nf90_error(status, 'nf90_def_var '// &
            trim(ent_names(k)))
       status=nf90_put_att(ncid,varid,"long_name", trim(Ent_title(k)))
       call handle_nf90_error(status,  'nf90_put_att  long_name')
       status=nf90_put_att(ncid,varid,"short_name", trim(ent_names(k)))
       status=nf90_put_att(ncid,varid,"units", units)
       status=nf90_put_att(ncid,varid,'_FillValue', -1.e30)
       call handle_nf90_error(status, 'nf90_put_att Ent vars '// &
            trim(ent_names(k)))
    end do

    status=nf90_enddef(ncid)
    call handle_nf90_error(status, 'my_nf90_create_Ent_vartime')

    !Assing values to time variable.
    time(:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
    status = nf90_inq_varid(ncid, 'time',varid)
    call handle_nf90_error(status, 'nf90_inq_varid time')
    status = nf90_put_var(ncid, varid, time)
    call handle_nf90_error(status, 'nf90_put_var_int1')

end subroutine my_nf90_create_Ent_vartime


!Keep local to each program to define own values.
subroutine my_nf90_defglobal_template(file)
!Dummy template routine for assigning global attributes to data file.
!Data processing program should make a copy and customize locally.
use netcdf
use convertnc
character*(*) :: file
!---Local---
integer :: status, ncid
character*1024 :: text
character(8) :: date

     !Enter define mode
    status=nf90_open(trim(file), NF90_WRITE, ncid)
    write(0,*) 'nf90_open ', status, file, ncid
    call handle_nf90_error(status,  'nf90_open '//file)
    status=nf90_redef(ncid)
    write(0,*) 'nf90_redef ', status, ncid
    !call handle_nf90_error(status,  'nf90_redef')
    !Put metadata global attributes
    text = 'Ent Global Vegetation Structure Dataset '// &
         '(Ent GVSD) v0.2  MODIS-Monfreda '
    !write(0,*) 'len, text: ', len(trim(text)), trim(text)
    status=nf90_put_att(ncid, NF90_GLOBAL, 'Description', trim(text))
    call handle_nf90_error(status, 'nf90_global'//trim(text))
    text = 'TEMPORARY working version '// &
         'Say something about source resolution, features '// &
         'With ext1 files having crops LAI extended by 5 grid cells.'
    !write(0,*) 'len, text: ',len(trim(text)), trim(text)
    status=nf90_put_att(ncid, NF90_GLOBAL, 'Comments', trim(text))
    call handle_nf90_error(status, '')
    text = 'Institution:  NASA Goddard Institute for Space Studies'
    status=nf90_put_att(ncid, NF90_GLOBAL, 'Institution', trim(text))
    call handle_nf90_error(status, '')
    text = 'Nancy.Y.Kiang@nasa.gov'
    status=nf90_put_att(ncid, NF90_GLOBAL, 'Contact', trim(text))
    call handle_nf90_error(status, '')
    call DATE_AND_TIME(date)
    text = date//' Created'
    status=nf90_put_att(ncid, NF90_GLOBAL, 'History', trim(text))
    call handle_nf90_error(status,  'nf90_put_att NF90_GLOBAL')

    !      status=nf90_enddef(ncid)      
    !      call handle_nf90_error(status,  'nf90_enddef my_nf90_defglobal')

end subroutine my_nf90_defglobal_template


end module EntGVSD_netcdf_util
      
