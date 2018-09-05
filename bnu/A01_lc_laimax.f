!  Program to assign 1kmx1km BNU LAImax to EntPFTs lc

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A01_lc_laimax.f

! gfortran -o myExe arrayutil.o A02_lc_laimax.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

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

!------------------------------------------------------------------------

      program lc_laimax

       ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
       ! interpolate to coarser resolutions.
      implicit none

      include 'netcdf.inc'

      !***************************************************
!*      PREFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*16, parameter :: EntPFT_files1(20) =
     &     (/
     &     'water           ',
     &     '01_ever_br_early',
     &     '02_ever_br_late ',
     &     '03_ever_nd_early',
     &     '04_ever_nd_late ',
     &     '05_cold_br_early',
     &     '06_cold_br_late ',
     &     '07_drought_br   ',
     &     '08_decid_nd     ',
     &     '09_cold_shrub   ',
     &     '10_arid_shrub   ',
     &     '11_c3_grass_per ',
     &     '12_c4_grass     ',
     &     '13_c3_grass_ann ',
     &     '14_c3_grass_arct',
     &     '15_crops_c3_herb',
     &     '16_crops_c4_herb',
     &     '17_crops_woody  ',
     &     '18_snow_ice     ',
     &     '19_bare_sparse  '
     &     /)

!***************************************************
!*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*14, parameter :: EntPFT_files2(20) =
     &     (/
     &     'water_lc      ',
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
     &     'crops_c3_herb ',
     &     'crops_c4_herb ',
     &     'crops_woody   ',
     &     'snow_ice      ',
     &     'bare_sparse   '
     &     /)

     
      integer, parameter :: XM = 43200 !long at 1 km
      integer, parameter :: YM = 21600 !lat at 1 km

      integer, parameter :: IM = XM !long at 1 km
      integer, parameter :: JM = YM !lat at 1 km
      
      integer, parameter :: longi = XM
      integer, parameter :: latin = YM

      real*4, parameter :: undef = -1e30 !AViewer undef value      !## UPDATE ME

      character*256 :: PathFilepre,PathFilepost

      integer, parameter :: NUMLAYERS = 40  !## UPDATE ME

      integer, parameter :: NUMLAYERSLC = 20  ! ONLY LC LAYERS
      integer, parameter :: NUMTITLES = 40  ! ALL LC AND LAI TITLES

      character*50 :: PATHin, PATHout
      real*4 :: LAYERIN
      real*4 :: LAYEROUT
      real*4 :: OFFIA, DLATAH, DLATA, OFFIB, DLATB, DATMIS

      character*256 :: filelai,fileout,filepft,filewater
      character*20 :: inqvarin

      real*4 :: LC(IM,JM), LAI(IM,JM), lailc(IM,JM)
      real*4 :: lon(IM),lat(JM)

      integer :: i, j, k, f, p, z

      integer :: err,dimidx,dimidy

      integer :: laifileid,laivarid

      integer :: lcfileid(20),lcvarid(20)
      integer :: varidout(20),fileidout(20)

      integer :: start2d(2),count2d(2)
      integer :: startX(1),startY(1),countX(1),countY(1)
      integer :: varidx,varidy

      call calc_lon_lat(IM,JM,lon,lat)

!      LAI max
      PathFilepre= '../../data/LAI/'
      PathFilepost = 'global_30s_2004_max.nc'
      filelai = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filelai,0,laifileid)
      write(*,*) err
      err = NF_INQ_VARID(laifileid,'lai',laivarid)
      write(*,*) err
      err = NF_INQ_VARID(laifileid,'lon',varidx)
      write(*,*) err
      err = NF_INQ_VARID(laifileid,'lat',varidy)
      write(*,*) err

!     ENTPFTLC
      do k = 1,20
         PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         filepft = trim(PathFilepre)//trim(EntPFT_files1(k))//'_lc.nc'
         err = NF_OPEN(filepft,0,lcfileid(k))
         write(*,*) err,filepft
         inqvarin = trim(EntPFT_files2(k))
         err = NF_INQ_VARID(lcfileid(k),inqvarin,lcvarid(k))
         write(*,*) err,lcfileid(k),lcvarid(k)
      enddo

!     Files out
      do k = 1,20
         PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         filepft = trim(PathFilepre)//trim(EntPFT_files1(k))//'_lai.nc'
         err = NF_OPEN(filepft,NF_WRITE,fileidout(k))
         write(*,*) err
         if (k.eq.1) then
            inqvarin = 'water_lai'
         else
            inqvarin = trim(EntPFT_files2(k))
         endif
         err = NF_INQ_VARID(fileidout(k),inqvarin,varidout(k))
         write(*,*) err,fileidout(k),varidout(k)
         err=NF_PUT_VARA_REAL(fileidout(k),varidx,1,IM,lon)
         write(*,*) err
         err=NF_PUT_VARA_REAL(fileidout(k),varidy,1,JM,lat)
      write(*,*) err

      enddo

!-----------------------------------------------------------------
!     Loop for every pft and  grid point
      
      do p = 1,NUMLAYERSLC

         write(*,*) 'processing ',EntPFT_files2(p), p
    
!**   LAI data, lon lat from LAI file                

         start2d(1)=1
         start2d(2)=1
         count2d(1)=IM
         count2d(2)=JM
         
         err = NF_GET_VARA_REAL(laifileid,laivarid,start2d,
     &              count2d,LAI)
         
!     write(*,*) err,'LAI '
               
         err = NF_GET_VARA_REAL(lcfileid(p),lcvarid(p),
     &              start2d,count2d,LC) 
!     write(*,*) err

         do i=1,IM
            do j=1,JM
               lailc(i,j) = LC(i,j)*LAI(i,j)
            enddo
         enddo
         
         err=NF_PUT_VARA_REAL(fileidout(p),varidout(p),
     &              start2d,count2d,lailc)
         
!                  write(*,*) err, 'put'
         
         
      enddo
      
      end program lc_laimax
      
      
      

