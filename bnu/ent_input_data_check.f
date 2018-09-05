!  Program to assign 1kmx1km monthly LAI to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include ent_input_data_check.f

! gfortran -o myExe arrayutil.o ent_input_data_check.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

! ./myExe

      subroutine calc_lon_lat(IM,JM,lon,lat)
      implicit none
      integer,intent(in) :: IM,JM
      real*4,intent(inout) :: lon(IM),lat(JM)
      integer :: i,j
     
      do i=1,IM
         lon(i) = -180. + (360./IM)*i - (360./IM)/2.
      end do

      do j=1,JM
         lat(j) = -90. + (180./JM)*j - (180./JM)/2.
      end do
     
      end subroutine calc_lon_lat 

      
      program ent_input_data_check

      implicit none

      include 'netcdf.inc'

            !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES - short names   *
      !***************************************************
      character*50, parameter :: EntPFT_shorttitle(18) =
     &     (/
     &     "ever_br_early   ",
     &     "ever_br_late    ",
     &     "ever_nd_early   ",
     &     "ever_nd_late    ",
     &     "cold_br_early   ",
     &     "cold_br_late    ",
     &     "drought_br      ",
     &     "decid_nd        ",
     &     "cold_shrub      ",
     &     "arid_shrub      ",
     &     "c3_grass_per    ",
     &     "c4_grass        ",
     &     "c3_grass_ann    ",
     &     "c3_grass_arct   ",
     &     "crops_herb      ",
     &     "crops_woody     ",
     &     "bare_bright     ",
     &     "bare_dark       "
     &     /)


       character*40, parameter :: EntPFT_title(18) =
     &     (/
     &     '1 - Evergreen Broadleaf Early Succ      ',
     &     '2 - Evergreen Broadleaf Late Succ       ',
     &     '3 - Evergreen Needleleaf Early Succ     ',
     &     '4 - Evergreen Needleleaf Late Succ      ',
     &     '5 - Cold Deciduous Broadleaf Early Succ ',
     &     '6 - Cold Deciduous Broadleaf Late Succ  ',
     &     '7 - Drought Deciduous Broadleaf         ',
     &     '8 - Deciduous Needleleaf                ',
     &     '9 - Cold Adapted Shrub                  ',
     &     '10 - Arid Adapted Shrub                 ',
     &     '11 - C3 Grass Perennial                 ',
     &     '12 - C4 Grass                           ',
     &     '13 - C3 Grass Annual                    ',
     &     '14 - Arctic C3 Grass                    ',
     &     '15 - Crops Herb                         ',
     &     '16 - Crops Woody                        ',
     &     '17 - Bright Bare Soil                   ',
     &     '18 - Dark Bare Soil                     '
     &     /)
       
      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM1km = X1km !long at 1 km
      integer, parameter :: JM1km = Y1km !lat at 1 km
      
      real*4, parameter :: undef = -1.e30

      character*256 :: PathFilepre,PathFilepost

      character*256 :: filelaijul,filelaijan,file017,filelaimax
      character*256 :: file201,filelc
      character*20 :: inqvarin,inqvarout
      real*4 :: long(IM1km),lati(JM1km)

      integer :: fileidlaijul,fileidlaijan,fileidlaimax,fileidlc
      integer :: fileid201,fileid017
      integer :: fidlaijulnew,fidlaijannew,fidlaimaxnew
      integer :: fid201new,fid017new
      integer :: varidlaijul(18),varidlaijan(18),varidlaimax(18)
      integer :: varid017(18),varid201(18),varidx,varidy,varidlc(18)
      integer :: varidlaijulnew(18),varidlaijannew(18)
      integer :: varid017new(18),varid201new(18),varidlaimaxnew(18)

      real*4 :: laijul(IM1km,JM1km),laijan(IM1km,JM1km)
      real*4 :: lai201(IM1km,JM1km),laimax(IM1km,JM1km)
      real*4 :: lai017(IM1km,JM1km),laijulout(IM1km,JM1km)
      real*4 :: laijanout(IM1km,JM1km),lc(IM1km,JM1km)
      real*4 :: lai201out(IM1km,JM1km),lai017out(IM1km,JM1km)
      real*4 :: laimaxout(IM1km,JM1km)

      integer :: i, j, k
      integer :: err, status

      integer :: start2D(2),count2D(2),dd(2)
      integer :: dimidx,dimidy

      
      call calc_lon_lat(IM1km,JM1km,long,lati)

      start2D(1)=1
      start2D(2)=1
      count2D(1)=IM1km
      count2D(2)=JM1km

!------------------------------------------------------------
      PathFilepre= '../../LAI3g/lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_LAI3g16_lc_pure.nc'
      filelc  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(filelc,0,fileidlc)
      write(*,*) status, 'nf_open in ', filelc
      do k=1,18
         inqvarin = 'lc_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileidlc,inqvarin,varidlc(k))
         write(*,*) err,inqvarin
      enddo
      
!------------------------------------------------------------
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_Jul_pure.nc'
      filelaijul  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(filelaijul,0,fileidlaijul)
      write(*,*) status, 'nf_open in ', filelaijul
      do k=1,18
         inqvarin = 'lai_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileidlaijul,inqvarin,varidlaijul(k))
         write(*,*) err,inqvarin
      enddo

      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_Jul_pure_new2.nc'
      filelaijul  =  trim(PathFilepre)//trim(PathFilepost)
      err = NF_CREATE(filelaijul,NF_64BIT_OFFSET,fidlaijulnew)
      write(*,*) err, 'nf_create out ',trim(filelaijul)
      err=NF_DEF_DIM(fidlaijulnew,'lon',IM1km,dimidx)
      write(*,*) err
      err=NF_DEF_DIM(fidlaijulnew,'lat',JM1km,dimidy)
      write(*,*) err
      err=NF_DEF_VAR(fidlaijulnew,'lon',NF_FLOAT,1,dimidx,
     &     varidx)
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaijulnew,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaijulnew,varidx,
     &        'units',14,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(fidlaijulnew,'lat',NF_FLOAT,1,dimidy,
     &     varidy)
      err=NF_PUT_ATT_TEXT(fidlaijulnew,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaijulnew,varidy,
     &        'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fidlaijulnew,inqvarout,NF_FLOAT,2,dd,
     &       varidlaijulnew(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(fidlaijulnew,varidlaijulnew(k),
     &        'long_name',45,(EntPFT_title(k))//("- LAI"))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fidlaijulnew,varidlaijulnew(k),
     &        'units',5,"m2/m2")
         write(*,*) err
         err=NF_PUT_ATT_REAL(fidlaijulnew, varidlaijulnew(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
         write(*,*) err
      enddo
      err=NF_ENDDEF(fidlaijulnew)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fidlaijulnew,varidx,1,IM1km,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fidlaijulnew,varidy,1,JM1km,lati)
      write(*,*) err
      
!------------------------------------------------------------
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_Jan_pure.nc'
      filelaijan  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(filelaijan,0,fileidlaijan)
      write(*,*) status, 'nf_open in ', filelaijan
      do k=1,18
         inqvarin = 'lai_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileidlaijan,inqvarin,varidlaijan(k))
         write(*,*) err,inqvarin
      enddo
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_Jan_pure_new2.nc'
      filelaijan  =  trim(PathFilepre)//trim(PathFilepost)
      err = NF_CREATE(filelaijan,NF_64BIT_OFFSET,fidlaijannew)
      write(*,*) err, 'nf_create out ',trim(filelaijan)
      err=NF_DEF_DIM(fidlaijannew,'lon',IM1km,dimidx)
      write(*,*) err
      err=NF_DEF_DIM(fidlaijannew,'lat',JM1km,dimidy)
      write(*,*) err
      err=NF_DEF_VAR(fidlaijannew,'lon',NF_FLOAT,1,dimidx,
     &     varidx)
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaijannew,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaijannew,varidx,
     &        'units',14,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(fidlaijannew,'lat',NF_FLOAT,1,dimidy,
     &     varidy)
      err=NF_PUT_ATT_TEXT(fidlaijannew,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaijannew,varidy,
     &        'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fidlaijannew,inqvarout,NF_FLOAT,2,dd,
     &       varidlaijannew(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(fidlaijannew,varidlaijannew(k),
     &        'long_name',45,(EntPFT_title(k))//("- LAI"))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fidlaijannew,varidlaijannew(k),
     &        'units',5,"m2/m2")
         write(*,*) err
         err=NF_PUT_ATT_REAL(fidlaijannew, varidlaijannew(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
         write(*,*) err
      enddo
      err=NF_ENDDEF(fidlaijannew)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fidlaijannew,varidx,1,IM1km,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fidlaijannew,varidy,1,JM1km,lati)
      write(*,*) err

!------------------------------------------------------------
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_laimax_pure.nc'
      filelaimax  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(filelaimax,0,fileidlaimax)
      write(*,*) status, 'nf_open in ', filelaimax
      do k=1,18
         inqvarin = 'lai_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileidlaimax,inqvarin,varidlaimax(k))
         write(*,*) err,inqvarin
      enddo
      
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_laimax_pure_new2.nc'
      filelaimax  =  trim(PathFilepre)//trim(PathFilepost)
      err = NF_CREATE(filelaimax,NF_64BIT_OFFSET,fidlaimaxnew)
      write(*,*) err, 'nf_create out ',trim(filelaimax)
      err=NF_DEF_DIM(fidlaimaxnew,'lon',IM1km,dimidx)
      write(*,*) err
      err=NF_DEF_DIM(fidlaimaxnew,'lat',JM1km,dimidy)
      write(*,*) err
      err=NF_DEF_VAR(fidlaimaxnew,'lon',NF_FLOAT,1,dimidx,
     &     varidx)
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaimaxnew,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaimaxnew,varidx,
     &        'units',14,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(fidlaimaxnew,'lat',NF_FLOAT,1,dimidy,
     &     varidy)
      err=NF_PUT_ATT_TEXT(fidlaimaxnew,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fidlaimaxnew,varidy,
     &        'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fidlaimaxnew,inqvarout,NF_FLOAT,2,dd,
     &       varidlaimaxnew(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(fidlaimaxnew,varidlaimaxnew(k),
     &        'long_name',45,(EntPFT_title(k))//("- LAI"))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fidlaimaxnew,varidlaimaxnew(k),
     &        'units',5,"m2/m2")
         write(*,*) err
         err=NF_PUT_ATT_REAL(fidlaimaxnew, varidlaimaxnew(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
         write(*,*) err
      enddo
      err=NF_ENDDEF(fidlaimaxnew)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fidlaimaxnew,varidx,1,IM1km,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fidlaimaxnew,varidy,1,JM1km,lati)
      write(*,*) err

!------------------------------------------------------------
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_017_pure.nc'
      file017  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(file017,0,fileid017)
      write(*,*) status, 'nf_open in ', file017
      do k=1,18
         inqvarin = 'lai_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileid017,inqvarin,varid017(k))
         write(*,*) err,inqvarin
      enddo
      
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_017_pure_new2.nc'
      file017  =  trim(PathFilepre)//trim(PathFilepost)
      err = NF_CREATE(file017,NF_64BIT_OFFSET,fid017new)
      write(*,*) err, 'nf_create out ',trim(file017)
      err=NF_DEF_DIM(fid017new,'lon',IM1km,dimidx)
      write(*,*) err
      err=NF_DEF_DIM(fid017new,'lat',JM1km,dimidy)
      write(*,*) err
      err=NF_DEF_VAR(fid017new,'lon',NF_FLOAT,1,dimidx,
     &     varidx)
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fid017new,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fid017new,varidx,
     &        'units',14,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(fid017new,'lat',NF_FLOAT,1,dimidy,
     &     varidy)
      err=NF_PUT_ATT_TEXT(fid017new,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fid017new,varidy,
     &        'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fid017new,inqvarout,NF_FLOAT,2,dd,
     &       varid017new(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(fid017new,varid017new(k),
     &        'long_name',45,(EntPFT_title(k))//("- LAI"))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fid017new,varid017new(k),
     &        'units',5,"m2/m2")
         write(*,*) err
         err=NF_PUT_ATT_REAL(fid017new, varid017new(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
         write(*,*) err
      enddo
      err=NF_ENDDEF(fid017new)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fid017new,varidx,1,IM1km,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fid017new,varidy,1,JM1km,lati)
      write(*,*) err

!------------------------------------------------------------
      
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_201_pure.nc'
      file201  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(file201,0,fileid201)
      write(*,*) status, 'nf_open in ', file201
      do k=1,18
         inqvarin = 'lai_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileid201,inqvarin,varid201(k))
         write(*,*) err,inqvarin
      enddo
      
      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_BNU16_lai_201_pure_new2.nc'
      file201  =  trim(PathFilepre)//trim(PathFilepost)
      err = NF_CREATE(file201,NF_64BIT_OFFSET,fid201new)
      write(*,*) err, 'nf_create out ',trim(file201)
      err=NF_DEF_DIM(fid201new,'lon',IM1km,dimidx)
      write(*,*) err
      err=NF_DEF_DIM(fid201new,'lat',JM1km,dimidy)
      write(*,*) err
      err=NF_DEF_VAR(fid201new,'lon',NF_FLOAT,1,dimidx,
     &     varidx)
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fid201new,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fid201new,varidx,
     &        'units',14,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(fid201new,'lat',NF_FLOAT,1,dimidy,
     &     varidy)
      err=NF_PUT_ATT_TEXT(fid201new,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fid201new,varidy,
     &        'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fid201new,inqvarout,NF_FLOAT,2,dd,
     &       varid201new(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(fid201new,varid201new(k),
     &        'long_name',45,(EntPFT_title(k))//("- LAI"))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fid201new,varid201new(k),
     &        'units',5,"m2/m2")
         write(*,*) err
         err=NF_PUT_ATT_REAL(fid201new, varid201new(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
         write(*,*) err
      enddo
      err=NF_ENDDEF(fid201new)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fid201new,varidx,1,IM1km,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fid201new,varidy,1,JM1km,lati)
      write(*,*) err

!------------------------------------------------------------

      do k=1,18
         
         inqvarout = 'lc_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileidlc,varidlc(k),
     &        start2D,count2D,lc)
         write(*,*) err, 'get', k

         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileidlaijul,varidlaijul(k),
     &        start2D,count2D,laijul)
         write(*,*) err, 'get', k
         
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileidlaijan,varidlaijan(k),
     &        start2D,count2D,laijan)
         write(*,*) err, 'get', k
      
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileidlaimax,varidlaimax(k),
     &        start2D,count2D,laimax)
         write(*,*) err, 'get', k

         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileid017,varid017(k),
     &        start2D,count2D,lai017)
         write(*,*) err, 'get', k

         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileid201,varid201(k),
     &        start2D,count2D,lai201)
         write(*,*) err, 'get', k

         do i = 1,IM1km
            
            do j = 1,JM1km

               if (lc(i,j).gt.0.and.laijul(i,j).lt.0) then
                  laijulout(i,j) = 0.
               else
                  laijulout(i,j) = laijul(i,j)
               endif

               if (lc(i,j).gt.0.and.laijan(i,j).lt.0) then
                  laijanout(i,j) = 0.
               else
                  laijanout(i,j) = laijan(i,j)
               endif

               if (lc(i,j).gt.0.and.laimax(i,j).lt.0) then
                  laimaxout(i,j) = 0.
               else
                  laimaxout(i,j) = laimax(i,j)
               endif
               
               if (lc(i,j).gt.0.and.lai017(i,j).lt.0) then
                  lai017out(i,j) = 0.
               else
                  lai017out(i,j) = lai017(i,j)
               endif
               
               if (lc(i,j).gt.0.and.lai201(i,j).lt.0) then
                  lai201out(i,j) = 0.
               else
                  lai201out(i,j) = lai201(i,j)
               endif
                              
            enddo

         enddo
         
         err=NF_PUT_VARA_REAL(fidlaijulnew,varidlaijulnew(k),
     &        start2D,count2D,laijulout)
         write(*,*) err, 'put', k
         
         err=NF_PUT_VARA_REAL(fidlaijannew,varidlaijannew(k),
     &        start2D,count2D,laijanout)
         write(*,*) err, 'put', k
         
         err=NF_PUT_VARA_REAL(fidlaimaxnew,varidlaimaxnew(k),
     &        start2D,count2D,laimaxout)
         write(*,*) err, 'put', k
         
         err=NF_PUT_VARA_REAL(fid017new,varid017new(k),
     &        start2D,count2D,lai017out)
         write(*,*) err, 'put', k
         
         err=NF_PUT_VARA_REAL(fid201new,varid201new(k),
     &        start2D,count2D,lai201out)
         write(*,*) err, 'put', k

      enddo

      err = NF_CLOSE(fidlaijulnew)
      write(*,*) err, 'close'
      err = NF_CLOSE(fidlaijannew)
      write(*,*) err, 'close'
      err = NF_CLOSE(fidlaimaxnew)
      write(*,*) err, 'close'
      err = NF_CLOSE(fid017new)
      write(*,*) err, 'close'
      err = NF_CLOSE(fid201new)
      write(*,*) err, 'close'

      end program ent_input_data_check
      
      

