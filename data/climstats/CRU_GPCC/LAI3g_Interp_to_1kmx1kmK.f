!* LAI3g_Interp_to_1kmx1km.f - Specifically set up to interpolate 1/12º LAI3g data to 1 km resolution

!* COMPILATION MUST BE MADE LOGGED IN "discover-sp3" IN THE FOLLOWING WAY:

!* ulimit -s unlimited
!* module purge
!* module load other/comp/gcc-4.9.2-sp3
!* gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
!* gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include LAI3g_Interp_to_1kmx1km.f
!* gfortran -o myExe arrayutil.o LAI3g_Interp_to_1kmx1km.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

!* In this case, the executable is called "myExe"

!!!!  CHECK LAND MASK FOR COASTLINES AND ASSIGN WTIN VALUES.

!* To avoid memory problems:
!  > ulimit -s 512000


!------------------------------------------------------------------------

      program LAI3g_Interp_to_1kmx1kmK
      
      implicit none
      include 'netcdf.inc'
      integer, parameter :: IM1km = 43200 !long at 1 km
      integer, parameter :: JM1km = 21600 !lat at 1 km
      integer, parameter :: IM10km = 4320 !long at 10 km
      integer, parameter :: JM10km = 2160 !lat at 10 km
      integer, parameter :: IM12 = 4320 !long at 1/12 degrees
      integer, parameter :: JM12 = 2160 !lat at 1/12 degrees
      integer, parameter :: IMH = 720 !long at 0.5 degrees
      integer, parameter :: JMH = 360 !lat at 0.5 degrees
      integer, parameter :: IM1 = 360 !long at 1 degrees
      integer, parameter :: JM1 = 180 !lat at 1 degrees
      integer, parameter :: IM2 = 144 !long at 2.5 degrees
      integer, parameter :: JM2 = 90 !lat at 2 degrees
      integer, parameter :: IM4X5 = 72 !long at 5 degrees
      integer, parameter :: JM4X5 = 46 !lat at 4 degrees

      integer, parameter :: longin = IM12
      integer, parameter :: latin = JM12
      integer, parameter :: longout = IM1km
      integer, parameter :: latout = JM1km

      real*4, parameter :: undef_in = -9999.0
      real*4, parameter :: undef_out = -1.e30
      
      character*80 :: TITLE !, TITLECHECK
      character*256 :: filein, fileout, filecheck
      character*50 :: PATHEnt, PATHfile
      !real*4 :: varin(longin,latin)
      real*4, ALLOCATABLE :: varin(:,:)
      real*4 :: varout(longout,latout)
      real*4, ALLOCATABLE :: LAYERIN(:,:)
      !real*4, dimension(longin,latin) :: LAYERIN
      real*4, dimension(longout,latout) :: LAYEROUT
      real*4, dimension(longin,latin) :: WTIN
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS
      real*4 :: lai(longout,latout)
      real*4 :: lon(longout),lat(latout)
      integer :: ncidin,ncidout,varid,status
      integer :: i, j
      character*20 :: inqvarin, inqvarout
 
      ALLOCATE(varin(longin,latin) )
      ALLOCATE( LAYERIN(longin,latin) )

      LAYEROUT(:,:) = undef_out
      DLATA = (180./latin) * 60.
      OFFIA=0.0
      DLATB = (180./latout) * 60.
      OFFIB=0.0
      DATMIS=undef_out
      
      !* Setup grids.
      write(*,*) longin,latin,OFFIA,DLATA,
     &     longout,latout,OFFIB,DLATB,DATMIS
       
      write(*,*) 'Calling HNTR40'
      write(*,*) 'Finished HNTR40'

      !* Input file.
      filein = 'cru_ts3.22_TS_means_1981-2010_HXH_bis.nc'
      status = nf_open(filein,NF_WRITE,ncidin)
      write(*,*) status, 'nf_open in ',trim(filein)

      end program LAI3g_Interp_to_1kmx1kmK

