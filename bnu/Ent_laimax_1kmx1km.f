!  Program to assign 1kmx1km monthly LAI to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include Ent_laimax_1kmx1km.f

! gfortran -o myExe arrayutil.o Ent_laimax_1kmx1km.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

! ./myExe

!------------------------------------------------------------------------

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

!------------------------------------------------------------------------

      program Ent_laimax_1kmx1km

      implicit none

      include 'netcdf.inc'
      
!*      ENTPFTS FILES FOR LC
      character*23, parameter :: EntPFT_files(19) =
     &     (/
     &     '01_ever_br_early_lai.nc',
     &     '02_ever_br_late_lai.nc ',
     &     '03_ever_nd_early_lai.nc',
     &     '04_ever_nd_late_lai.nc ',
     &     '05_cold_br_early_lai.nc',
     &     '06_cold_br_late_lai.nc ',
     &     '07_drought_br_lai.nc   ',
     &     '08_decid_nd_lai.nc     ',
     &     '09_cold_shrub_lai.nc   ',
     &     '10_arid_shrub_lai.nc   ',
     &     '11_c3_grass_per_lai.nc ',
     &     '12_c4_grass_lai.nc     ',
     &     '13_c3_grass_ann_lai.nc ',
     &     '14_c3_grass_arct_lai.nc',
     &     'crops_herb_laimax.nc   ',
     &     '17_crops_woody_lai.nc  ',
     &     '18_snow_ice_lai.nc     ',
     &     '19_bare_sparse_lai.nc  ',
     &     'water_lai.nc           '
     &     /)

      character*17, parameter :: EntPFT_names(19) =
     &     (/
     &     'ever_br_early',
     &     'ever_br_late ',
     &     'ever_nd_early',
     &     'ever_nd_late ',
     &     'cold_br_early',
     &     'cold_br_late ',
     &     'drought_br   ',
     &     'decid_nd     ',
     &     'cold_shrub   ',
     &     'arid_shrub   ',
     &     'c3_grass_per ',
     &     'c4_grass     ',
     &     'c3_grass_ann ',
     &     'c3_grass_arct',
     &     'crops_herb   ',
     &     'crops_woody  ',
     &     'snow_ice     ',
     &     'bare_sparse  ',
     &     'water_lai    '
     &     /)

      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM = X1km !long at 1 km
      integer, parameter :: JM = Y1km !lat at 1 km
      
      real*4, parameter :: undef = -1.e30

      character*256 :: PathFilepre,PathFilepost

      character*256 :: filein,fileout
      character*20 :: inqvarin

      integer :: idhin(19),idhout(19)
      integer :: varidhin(19),varidhout(19)

      real*4 :: hcin(IM,JM)
      real*4 :: lon(IM),lat(JM)

      integer :: i, j, k
      integer :: err, status

      integer :: start2d(2),count2d(2)
      integer :: dd(2),varidx,varidy
 
!------------------------------------------------------------------------
!     Lon and Lat
      call calc_lon_lat(IM,JM,lon,lat)
     
!------------------------------------------------------------------------

      PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= 'Ent_laimax_1kmx1km.nc'
      fileout  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(fileout,nf_write,idhout)
      write(*,*) status, 'nf_open in ', fileout
      err = nf_inq_varid(idhout,'lon',varidx)
      write(*,*) err
      err = nf_inq_varid(idhout,'lat',varidy)
      write(*,*) err
      do k = 1,19
         if (k.eq.19) then
            inqvarin = 'lai_water'
         else
            inqvarin = trim('lai_')//trim(EntPFT_names(k))
         endif
         err=NF_INQ_VARID(idhout,inqvarin,varidhout(k))
         write(*,*) err,'varid',varidhout(k),EntPFT_names(k)
      enddo
      
      do k = 1,19
         
         PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         PathFilepost= EntPFT_files(k)
         filein  =  trim(PathFilepre)//trim(PathFilepost)
         status = nf_open(filein,0,idhin(k))
         write(*,*) status, 'nf_open in ', filein
         inqvarin = EntPFT_names(k)
         status = nf_inq_varid(idhin(k),inqvarin,varidhin(k))
         write(*,*) status,'nf_inq_varid ',inqvarin
         
         start2D(1)=1
         start2D(2)=1
         count2D(1)=IM
         count2D(2)=JM
         
         hcin(:,:) = undef
         status = NF_GET_VARA_REAL(idhin(k),varidhin(k),start2d,count2d,
     &        hcin)
         write(*,*) status,'nf_get_var_real ',varidhin(k)

	 do i = 1,IM
            do j = 1,JM
               if (hcin(i,j).lt.0.) then
                  hcin(i,j) = 0.
               else
                  hcin(i,j) = hcin(i,j)
               endif
            enddo
         enddo

         err = NF_PUT_VARA_REAL(idhout,varidhout(k),
     &        start2D,count2D,hcin)
         write(*,*) err, 'put', varidhout(k)
         
         err = NF_CLOSE(idhin(k))
         write(*,*) err

      enddo

      err=NF_PUT_VARA_REAL(idhout,varidx,1,IM,lon)
      write(*,*) err
      err=NF_PUT_VARA_REAL(idhout,varidy,1,JM,lat)
      write(*,*) err


      err = NF_CLOSE(idhout)
      write(*,*) err

      
      end program Ent_laimax_1kmx1km
      

