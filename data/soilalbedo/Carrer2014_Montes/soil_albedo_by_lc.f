!  Program to assign 1kmx1km monthly LAI to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include soil_albedo_by_lc.f

! gfortran -o myExe arrayutil.o soil_albedo_by_lc.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

! ./myExe

!------------------------------------------------------------------------

      subroutine calc_lon_lat(IM,JM,lon,lat)
      implicit none
      integer,intent(in) :: IM,JM
      real*4,intent(inout) :: lon(IM),lat(JM)
!      real*8, parameter :: deltlon = 360./IM
!      real*8, parameter ::  detlat = 180./JM
      integer :: i,j

!      deltlon = 360./IM
!      detlat = 180./JM
     
      do i=1,IM
         lon(i) = -180. + (360./IM)*i - (360./IM)/2.
      end do

      do j=1,JM
         lat(j) = -90. + (180./JM)*j - (180./JM)/2.
      end do
     
      end subroutine calc_lon_lat 

!------------------------------------------------------------------------

      program soil_albedo_by_lc

      implicit none

      include 'netcdf.inc'

!*      ENTPFTS FILES FOR LC
      character*23, parameter :: EntPFT_files(19) =
     &     (/
     &     '01_ever_br_early_lc.nc ',
     &     '02_ever_br_late_lc.nc  ',
     &     '03_ever_nd_early_lc.nc ',
     &     '04_ever_nd_late_lc.nc  ',
     &     '05_cold_br_early_lc.nc ',
     &     '06_cold_br_late_lc.nc  ',
     &     '07_drought_br_lc.nc    ',
     &     '08_decid_nd_lc.nc      ',
     &     '09_cold_shrub_lc.nc    ',
     &     '10_arid_shrub_lc.nc    ',
     &     '11_c3_grass_per_lc.nc  ',
     &     '12_c4_grass_lc.nc      ',
     &     '13_c3_grass_ann_lc.nc  ',
     &     '14_c3_grass_arct_lc.nc ',
     &     'crops_herb_lc.nc       ',
     &     '17_crops_woody_lc.nc   ',
     &     '18_snow_ice_lc.nc      ',
     &     '19_bare_sparse_lc.nc   ',
     &     'water_lc.nc            '
     &     /)

      character*14, parameter :: EntPFT_names(19) =
     &     (/
     &     'ever_br_early ',
     &     'ever_br_late  ',
     &     'ever_nd_early ',
     &     'ever_nd_late  ',
     &     'cold_br_early ',
     &     'cold_br_late  ',
     &     'drought_br    ',
     &     'decid_nd      ',
     &     'cold_shrub    ',
     &     'arid_shrub    ',
     &     'c3_grass_per  ',
     &     'c4_grass      ',
     &     'c3_grass_ann  ',
     &     'c3_grass_arct ',
     &     'crops_herb    ',
     &     'crops_woody   ',
     &     'snow_ice      ',
     &     'bare_sparse   ',
     &     'water         '
     &     /)

!     *      LONG NAMES
      character*54, parameter :: EntPFT_lnames(19) =
     &     (/
     &     "1 - Evergreen Broadleaf Early Succ - Soil Albedo      ",
     &     "2 - Evergreen Broadleaf Late Succ - Soil Albedo       ",
     &     "3 - Evergreen Needleleaf Early Succ - Soil Albedo     ",
     &     "4 - Evergreen Needleleaf Late Succ - Soil Albedo      ",
     &     "5 - Cold Deciduous Broadleaf Early Succ - Soil Albedo ",
     &     "6 - Cold Deciduous Broadleaf Late Succ - Soil Albedo  ",
     &     "7 - Drought Deciduous Broadleaf - Soil Albedo         ",
     &     "8 - Deciduous Needleleaf - Soil Albedo                ",
     &     "9 - Cold Adapted Shrub - Soil Albedo                  ",
     &     "10 - Arid Adapted Shrub - Soil Albedo                 ",
     &     "11 - C3 Grass Perennial - Soil Albedo                 ",
     &     "12 - C4 Grass - Soil Albedo                           ",
     &     "13 - C3 Grass Annual - Soil Albedo                    ",
     &     "14 - Arctic C3 Grass - Soil Albedo                    ",
     &     "15 - Crops Herb - Soil Albedo                         ",
     &     "16 - Crops Woody - Soil Albedo                        ",
     &     "17 - Bright Bare Soil - Soil Albedo                   ",
     &     "18 - Dark Bare Soil - Soil Albedo                     ",
     &     "19 - Water - Soil Albedo                              "
     &     /)
     
      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM = X1km !long at 1 km
      integer, parameter :: JM = Y1km !lat at 1 km
      
      real*4, parameter :: undef = -1e30 !AViewer undef value

      character*256 :: PathFilepre,PathFilepost

      character*50 :: PATHin, PATHout

      character*256 :: fileinvis,fileinnir,fileinsw,fileout
      character*20 :: inqvarin

      integer :: ncidinvis,ncidinnir,ncidinsw
      integer :: ncidoutvis,ncidoutnir,ncidoutsw
      integer :: varidvis,varidnir,varidsw,status

      real*4 :: lcin(IM,JM)
      real*4 :: lon(IM),lat(JM)

      integer :: i, j, k, f, p, z

      integer :: err,dimidx,dimidy,dimidz
      integer :: ycoord, xcoord

      integer :: lcfileid,visfileid,nirfileid,swfileid

      integer :: pftfileid,pftvarid
      integer :: pftvaridvis(19),pftvaridnir(19),pftvaridsw(19)
      integer :: fileidout

      character*256 :: filepft

      integer :: startA(1),startB(2),countA(1),countB(2)
      integer :: startX(1),startY(1),countX(1),countY(1)
      integer :: dd(2),varidx,varidy

      real*4 :: albvisin(IM,JM)

      real*4 :: albvisout(IM,JM)
 
!------------------------------------------------------------------------

!     define Lon and Lat
      call calc_lon_lat(IM,JM,lon,lat)
     
!------------------------------------------------------------------------
      ! albedo VIS
      fileinvis = 'VIS_Alb_soil_yearly.006.1kmx1km.nc'
      status = nf_open(fileinvis,0,ncidinvis)
      write(*,*) status, 'nf_open in ', fileinvis
      inqvarin = 'mean'
      status = nf_inq_varid(ncidinvis,inqvarin,varidvis)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albvisin(:,:) = undef
      status = nf_get_var_real(ncidinvis,varidvis,albvisin)
      write(*,*) status,'nf_get_var_real ',inqvarin
      
!-----------------------------------------------------------------
      fileout='VIS_Alb_soil_yearly.byLandCover.1kmx1km.nc'
      err = NF_CREATE(fileout,NC_NOCLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL,ncidoutvis)
      write(*,*) err, fileout, ncidoutvis
      err=NF_DEF_DIM(ncidoutvis,'lon',IM,dimidx)
      write(*,*) err, 'dimidx',dimidx
      err=NF_DEF_DIM(ncidoutvis,'lat',JM,dimidy)
      write(*,*) err, 'dimidy',dimidy
      err=NF_DEF_VAR(ncidoutvis,'lon',NF_REAL,1,dimidx,
     &     varidx)
      write(*,*) err,'varidx',varidx
      err=NF_DEF_VAR(ncidoutvis,'lat',NF_REAL,1,dimidy,
     &     varidy)
      write(*,*) err,'varidy',varidy

      dd(1)=dimidx
      dd(2)=dimidy

      do k = 1,19
         err=NF_DEF_VAR(ncidoutvis,trim(EntPFT_names(k)),
     &     NF_REAL,2,dd,pftvaridvis(k))
         write(*,*) err,'varid',pftvaridvis(k),EntPFT_names(k)
      enddo

      do k = 1,19
         err=NF_PUT_ATT_TEXT(ncidoutvis,pftvaridvis(k),
     &       'long_name',54,trim(EntPFT_lnames(k)))
         write(*,*) err, trim(EntPFT_lnames(k))
         err=NF_PUT_ATT_TEXT(ncidoutvis,pftvaridvis(k),'units',
     &       8,'fraction')
         write(*,*) err
      enddo

!      err=NF_PUT_ATT_TEXT(ncidoutvis,NF_GLOBAL,'xlabel',
!     &     'VIS soil albedo by Ent PFT Land Cover; 16 plant',
!     &     'types, bare soil, water, ice, fraction',
!     &     'of grid cell')
!      write(*,*) err
!      err=NF_PUT_ATT_TEXT(ncidoutvis,NF_GLOBAL,'history',
!     &     'MODIS soil albedo (Carrer et al 2014 RSE)', !
!     &     'version 1, July 2017')
!      write(*,*) err
!      err=NF_PUT_ATT_TEXT(ncidoutvis,NF_GLOBAL,'institution',
!     &        'NASA/GISS  C. Montes, N. Kiang')
!      write(*,*) err

      err = NF_ENDDEF(ncidoutvis)
      write(*,*) err, 'enddef'
 
      err=NF_PUT_VARA_REAL(ncidoutvis,varidx,1,IM,lon)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidoutvis,varidy,1,JM,lat)
      write(*,*) err

      do k = 1,19
         PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         PathFilepost= EntPFT_files(k)
         filepft = trim(PathFilepre)//trim(PathFilepost)
	 write(*,*) filepft
         err = NF_OPEN(filepft,0,pftfileid)
	 write(*,*) err
         inqvarin = trim(EntPFT_names(k))
         err = NF_INQ_VARID(pftfileid,inqvarin,pftvarid)
	 write(*,*) err
         lcin(:,:) = undef
         status = nf_get_var_real(pftfileid,pftvarid,lcin)
         write(*,*) status,'nf_get_var_real ',inqvarin

         albvisout(:,:) = undef
         do i = 1,IM
            do j = 1,JM
               albvisout(i,j) = albvisin(i,j)*lcin(i,j)
            enddo
         enddo
         ! VIS
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidoutvis,pftvaridvis(k),
     &      startB,countB,albvisout)
         write(*,*) err, 'put'

      enddo


      err = NF_CLOSE(ncidoutvis)
      
      end program soil_albedo_by_lc
      

