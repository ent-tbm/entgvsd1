!  Program to assign 1kmx1km monthly LAI to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include Ent_lai_Jan_1kmx1km.f

! gfortran -o myExe arrayutil.o Ent_lai_Jan_1kmx1km.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

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

      program Ent_lai_Jan_1kmx1km

      implicit none

      include 'netcdf.inc'

      character*17, parameter :: EntPFT_names(19) =
     &     (/
     &     'lai_ever_br_early',
     &     'lai_ever_br_late ',
     &     'lai_ever_nd_early',
     &     'lai_ever_nd_late ',
     &     'lai_cold_br_early',
     &     'lai_cold_br_late ',
     &     'lai_drought_br   ',
     &     'lai_decid_nd     ',
     &     'lai_cold_shrub   ',
     &     'lai_arid_shrub   ',
     &     'lai_c3_grass_per ',
     &     'lai_c4_grass     ',
     &     'lai_c3_grass_ann ',
     &     'lai_c3_grass_arct',
     &     'lai_crops_herb   ',
     &     'lai_crops_woody  ',
     &     'lai_snow_ice     ',
     &     'lai_bare_sparse  ',
     &     'lai_water        '
     &      /)

      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM = X1km !long at 1 km
      integer, parameter :: JM = Y1km !lat at 1 km
      
      real*4, parameter :: undef = -1.e30

      character*256 :: PathFilepre,PathFilepost

      character*256 :: filein,fileout
      character*20 :: inqvarin

      integer :: idhin,idhout,idhc
      integer :: varidhin,varidhout(18),varidhc

      real*4 :: hcin(IM,JM)
      real*4 :: lon(IM),lat(JM)

      integer :: i, j, k
      integer :: err, status

      integer :: start2d(2),count2d(2),start3d(4),count3d(4)
      integer :: dd(4),varidx,varidy
 
!------------------------------------------------------------------------
!     Lon and Lat
      call calc_lon_lat(IM,JM,lon,lat)
     
!------------------------------------------------------------------------
      
      PathFilepre= '../lc_lai_ent/nc/'
      PathFilepost= 'EntMM_lc_lai_Jan_1kmx1km.nc'
      filein  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(filein,0,idhin)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'EntPFT'
      status = nf_inq_varid(idhin,inqvarin,varidhin)
      write(*,*) status,'nf_inq_varid ',inqvarin

      PathFilepre= '../lc_lai_ent/nc/'
      PathFilepost= 'crops_herb_lc_lai_Jan_1kmx1km.nc'
      filein  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(filein,0,idhc)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'crops_herb'
      status = nf_inq_varid(idhc,inqvarin,varidhc)
      write(*,*) status,'nf_inq_varid ',inqvarin
 
      PathFilepre= '../lc_lai_ent/nc/'
      PathFilepost= 'Ent_lai_Jan_1kmx1km.nc'
      fileout  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(fileout,nf_write,idhout)
      write(*,*) status, 'nf_open in ', fileout
      err = nf_inq_varid(idhout,'lon',varidx)
      write(*,*) err
      err = nf_inq_varid(idhout,'lat',varidy)
      write(*,*) err
      do k = 1,19
         err=NF_INQ_VARID(idhout,trim(EntPFT_names(k)),
     &        varidhout(k))
         write(*,*) err,'varid',varidhout(k),EntPFT_names(k)
      enddo

!-----------------------------------------------------------------

      err=NF_PUT_VARA_REAL(idhout,varidx,1,IM,lon)
      write(*,*) err
      err=NF_PUT_VARA_REAL(idhout,varidy,1,JM,lat)
      write(*,*) err


!-----------------------------------------------------------------

      ! first 14 pfts

      do k = 1,14
         
         start3d(1) = k
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM
         count3d(3) = JM
         
         hcin(:,:) = undef
         status = NF_GET_VARA_REAL(idhin,varidhin,start3d,count3d,hcin)
         write(*,*) status,'nf_get_var_real ',varidhin

         do i = 1,IM
            do j = 1,JM
               if (hcin(i,j).lt.0.) then
                  hcin(i,j) = 0.
               else
                  hcin(i,j) = hcin(i,j)
               endif
            enddo
         enddo
         
         start2D(1)=1
         start2D(2)=1
         count2D(1)=IM
         count2D(2)=JM

         err = NF_PUT_VARA_REAL(idhout,varidhout(k),
     &        start2D,count2D,hcin)
         write(*,*) err, 'put', varidhout(k)
         
      enddo

      hcin(:,:) = undef
      status = NF_GET_VARA_REAL(idhc,varidhc,start2d,count2d,hcin)
      write(*,*) status,'nf_get_var_real ',varidhc
      do i = 1,IM
         do j = 1,JM
            if (hcin(i,j).lt.0.) then
               hcin(i,j) = 0.
            else
               hcin(i,j) = hcin(i,j)
            endif
         enddo
      enddo
      
      k = 15.
      err = NF_PUT_VARA_REAL(idhout,varidhout(k),
     &        start2D,count2D,hcin)
      write(*,*) err, 'put', varidhout(k)
      
      do k = 16,18
         
         start3d(1) = k+1
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM
         count3d(3) = JM
         
         hcin(:,:) = undef
         status = NF_GET_VARA_REAL(idhin,varidhin,start3d,count3d,hcin)
         write(*,*) status,'nf_get_var_real ',varidhin

         do i = 1,IM
            do j = 1,JM
               if (hcin(i,j).lt.0.) then
                  hcin(i,j) = 0.
               else
                  hcin(i,j) = hcin(i,j)
               endif
            enddo
         enddo
         
         start2D(1)=1
         start2D(2)=1
         count2D(1)=IM
         count2D(2)=JM
         
         err = NF_PUT_VARA_REAL(idhout,varidhout(k),
     &        start2D,count2D,hcin)
         write(*,*) err, 'put', varidhout(k)
         
      enddo
      
!-----------------------------------------------------------------

      err = NF_CLOSE(idhin)
      write(*,*) err
      err = NF_CLOSE(idhc)
      write(*,*) err
      err = NF_CLOSE(idhout)
      write(*,*) err
      
      end program Ent_lai_Jan_1kmx1km
      

