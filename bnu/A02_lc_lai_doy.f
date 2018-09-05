!  Program to assign 1kmx1km BNU LAI of selecter DOY to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A02_lc_lai_doy.f

! gfortran -o myExe arrayutil.o A02_lc_lai_doy.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

! ./myExe

!------------------------------------------------------------------------

      program lc_lai_doy

       ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
       ! interpolate to coarser resolutions.
      implicit none

      include 'netcdf.inc'

      !***************************************************
!*      PREFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*3, parameter :: EntPFT_files1(19) =
     &     (/
     &     '01_',
     &     '02_',
     &     '03_',
     &     '04_',
     &     '05_',
     &     '06_',
     &     '07_',
     &     '08_',
     &     '09_',
     &     '10_',
     &     '11_',
     &     '12_',
     &     '13_',
     &     '14_',
     &     '15_',
     &     '16_',
     &     '17_',
     &     '18_',
     &     '19_'
     &     /)

!***************************************************
!*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*14, parameter :: EntPFT_files2(19) =
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
     &     'crops_c3_herb ',
     &     'crops_c4_herb ',
     &     'crops_woody   ',
     &     'snow_ice      ',
     &     'bare_sparse   '
     &     /)

     
      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM1km = X1km !long at 1 km
      integer, parameter :: JM1km = Y1km !lat at 1 km
      
      integer, parameter :: longin = 1
      integer, parameter :: latin = 1
      integer, parameter :: longout = 1 !Should be same at longin
      integer, parameter :: latout = 1

      real*4, parameter :: undef = -1e30 !AViewer undef value      !## UPDATE ME

      character*3, parameter :: DOY(2) =
     &     (/
     &     "017","201"
     &     /)

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

      real*4 :: LCIN_in, LAI, lai_lc
      real*4 :: lon,lat

      integer :: i, j, k, f, p, z

      integer :: err,dimidx,dimidy,dimidz
      integer :: ycoord, xcoord

      integer :: laifileid(12),waterfileid
      integer :: laivarid(12),watervarid

      integer :: pftfileid(19),pftvarid(19)
      integer :: varidout(12),fileidout(12)

      integer :: startA(1),startB(2),countA(1),countB(2)
      integer :: startX(1),startY(1),countX(1),countY(1)
      integer :: start3d(4),count3d(4)
      integer :: dd(4),varidx(12),varidy(12)


      !* Input file.

!     Monthly LAI
!      PathFilepre= '../../data/LAI/LAI3gMonthly/'
!      PathFilepost = 'LAI3g_'//MONTH(7)//'_2004_1kmx1km.nc'
!      filelai = trim(PathFilepre)//trim(PathFilepost)
!      err = NF_OPEN(filelai,NF_WRITE,laifileid)
!      err = NF_INQ_VARID(laifileid,'lai',laivarid)
!      err = NF_INQ_VARID(laifileid,'lon',varidx)
!      err = NF_INQ_VARID(laifileid,'lat',varidy)

!     DOY LAI
      do k = 1,2
         PathFilepre= '../../data/LAI/'
         PathFilepost = 'global_30s_2004_'//DOY(k)//'.nc'
         filelai = trim(PathFilepre)//trim(PathFilepost)
         err = NF_OPEN(filelai,NF_WRITE,laifileid(k))
	 write(*,*) err, filelai
         err = NF_INQ_VARID(laifileid(k),'lai',laivarid(k))
	 write(*,*) err
         err = NF_INQ_VARID(laifileid(k),'lon',varidx(k))
	 write(*,*) err
         err = NF_INQ_VARID(laifileid(k),'lat',varidy(k))
	 write(*,*) err
      enddo

!     Water LC
      PathFilepre= '../../LAI3g/lc_lai_ent/'
      PathFilepost= 'EntMM_lc_laimax_1kmx1km/'
      filewater = trim(PathFilepre)//trim(PathFilepost)//
     &     'water_lc.nc'
      err = NF_OPEN(filewater,NF_WRITE,waterfileid)
      write(*,*) err
      err = NF_INQ_VARID(waterfileid,'water_lc',watervarid)
      write(*,*) err

!     ENTPFTLC
      do k = 1,19
         PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         filepft = trim(PathFilepre)//trim(EntPFT_files1(k))//
     &        trim(EntPFT_files2(k))//'_lc.nc'
         err = NF_OPEN(filepft,NF_WRITE,pftfileid(k))
      	 write(*,*) err
         inqvarin = trim(EntPFT_files2(k))
         err = NF_INQ_VARID(pftfileid(k),inqvarin,pftvarid(k))
	 write(*,*) err
      enddo

!     Fileout
      do k = 1,2
         fileout = '../lc_lai_ent/nc/EntMM_lc_lai_'
     &        //DOY(k)//'_1kmx1km.nc'
         err = NF_OPEN(fileout,NF_WRITE,fileidout(k))
	 write(*,*) err, fileout
         inqvarin = 'EntPFT'
         err = NF_INQ_VARID(fileidout(k),inqvarin,varidout(k))
	 write(*,*) err
      enddo


!-----------------------------------------------------------------
!     Loop for every grid point
      
      do ycoord = 1,JM1km

         do xcoord = 1,IM1km

            if (xcoord.eq.1) then
               startY(1) = ycoord
               countY(1) = 1
               err = NF_GET_VARA_REAL(laifileid(1),varidy,startY,
     &              countY,lat)
            endif
           
!**   INPUT Files at 1km x 1km 

            startB(1)=xcoord
            startB(2)=ycoord
            countB(1)=1
            countB(2)=1
            
               
!**   LAI data ------------------------------------------------------------------------------

            do k = 1,2

               !**   lon lat from LAI file                

               if (ycoord.eq.1) then
                  startX(1) = xcoord
                  countX(1) = 1
                  err = NF_GET_VARA_REAL(laifileid(k),varidx(k),startX,
     &                 countX,lon)
               endif
               
               if (xcoord.eq.1) then
                  write(*,*) 'ycoord', ycoord, 'lat', lat, DOY(k)
               endif
               
               
               if (ycoord.eq.1) then
                  err = NF_PUT_VARA_REAL(fileidout(k),varidx(k),startX,
     &                 countX,lon)
               endif
               if (xcoord.eq.1) then
                  err = NF_PUT_VARA_REAL(fileidout(k),varidy(k),startY,
     &                 countY,lat)
               endif
                              
               err = NF_GET_VARA_REAL(laifileid(k),laivarid(k),startB,
     &              countB,LAI)
               
              
               do p = 1,NUMLAYERSLC

                  start3d(1) = p
                  start3d(2) = xcoord
                  start3d(3) = ycoord
                  count3d(1) = 1
                  count3d(2) = 1
                  count3d(3) = 1

                  if (p.eq.1) then
                     
                     err = NF_GET_VARA_REAL(waterfileid,watervarid,
     &                    startB,countB,LCIN_in)
!                     write(*,*) err
                     
                     lai_lc = LCIN_in*LAI

                     err = NF_PUT_VARA_REAL(fileidout(k),varidout(k),
     &                    start3d,count3d,lai_lc)
!                     write(*,*) err

		  endif

                  if (p.gt.1 .and. p.le.NUMLAYERSLC) then

                     err = NF_GET_VARA_REAL(pftfileid(p-1),
     &                    pftvarid(p-1),startB,countB,LCIN_in) 
!                     write(*,*) err
                     
                     lai_lc = LCIN_in*LAI

                     err=NF_PUT_VARA_REAL(fileidout(k),varidout(k),
     &                    start3d,count3d,lai_lc)
                     
!                     write(*,*) err

                  endif

               enddo
                
            enddo

          enddo

      enddo

      end program lc_lai_doy
      
