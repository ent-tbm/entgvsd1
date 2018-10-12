! This program converts EntMM 17 PFTs to Ent 16 PFTs + bright + dark.
! To combine C3 and C4 crops: #define COMBINE_CROPS_C3_C4
! 9/12/13 Fixed bright/dark soil:  need to do partitioning after
!         each step of trim/scale/no crops to account for new
!         bare soil cover added, especially in crop grid cells.
! 1/16/13 Added subroutine replace_crops for nocrops to replace cover and LAI 
!         of crops with dominant natural veg in grid cell. If no natural veg,
!         then searches adjacent grid cells.
!         Added subroutine fill_crops for _ext version of crop LAI, to fill
!         in crop LAI in some grid cells that have no crops in MODIS cover
!         but may have crop cover in Pongratz historical cover.  Called once
!         for single grid extension.  Can be called again for further filling.
! 3/17/14 Added ext1 files at end for lai max, lai monthly, and height, by
!         replacing 15-crops in laic, laim, hm, and hsd with the ext values for
!         crops.
! Compile the program with:
!
! ifort -cpp convert_VEG5.f -convert big_endian
!

! GFORTRAN COMPILATION ON DISCOVER:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A05_trim_lai_doy_1kmx1km.f
! gfortran -o myExe arrayutil.o A05_trim_lai_doy_1kmx1km.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
! ./myExe

! BEFORE RUNNING: mkdir ../lc_lai_ent16

#define COMBINE_CROPS_C3_C4
#define SPLIT_BARE_SOIL

      

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

!      subroutine omp_set_num_threads(num_threads)
!      integer, intent(in) :: num_threads

!------------------------------------------------------------------------

      program convert
      use omp_lib
      use conversions
      
      implicit none
      
      include 'netcdf.inc'


      !***************************************************
      !*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
      !***************************************************
      character*17, parameter :: EntPFT_files(20) =
     &     (/
     &     'water            ',
     &     '01_ever_br_early ',
     &     '02_ever_br_late  ',
     &     '03_ever_nd_early ',
     &     '04_ever_nd_late  ',
     &     '05_cold_br_early ',
     &     '06_cold_br_late  ',
     &     '07_drought_br    ',
     &     '08_decid_nd      ',
     &     '09_cold_shrub    ',
     &     '10_arid_shrub    ',
     &     '11_c3_grass_per  ',
     &     '12_c4_grass      ',
     &     '13_c3_grass_ann  ',
     &     '14_c3_grass_arct ',
     &     '15_crops_c3_herb ',
     &     '16_crops_c4_herb ',
     &     '17_crops_woody   ',
     &     '18_snow_ice      ',
     &     '19_bare_sparse   '
     &     /)

            character*14, parameter :: EntPFT_lcfields(20) =
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

            character*14, parameter :: EntPFT_laifields(20) =
     &     (/
     &     'water_lai     ',
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


      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES                 *
      !***************************************************
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
      

      real*4, parameter :: undef = -1.e30

      ! define input file resolution
      !* .5x.5
      !integer, parameter :: IM=720, JM=360, KM=20
      !character*(*), parameter :: res_in="05x05"
      !character*(*), parameter :: res_in_int="720x360"
      !* 1x1
      !integer, parameter :: IM=360, JM=180, KM=20
      !character*(*), parameter :: res_in="1x1"
      !character*(*), parameter :: res_in_int="360x180"
      !* 2.5x2
      integer, parameter :: IM1km=43200, JM1km=21600
      integer, parameter :: IM=1, JM=1, KM=20
!      integer, parameter :: IM=1440, JM=720, KM=20
!      integer, parameter :: IM=144, JM=90, KM=20
!      character*(*), parameter :: res_in="2.5x2"
!      character*(*), parameter :: res_in="144x90"
!      character*(*), parameter :: res_in_int="144x90"
!      character*(*), parameter :: res_in="1440x720"
!      character*(*), parameter :: res_in_int="1440x720"
      character*(*), parameter :: res_in="1kmx1km"
      character*(*), parameter :: res_out="1kmx1km"

      integer, parameter :: divx = IM
      integer, parameter :: divy = JM

      integer :: fileid,dimidx,dimidy,dimidz,dd(4)
      integer :: varidxlai,varidxlc,varidylai,varidylc
      integer :: startA(1),startB(2),countA(1),countB(2),lenx,leny

      integer :: start2d(2),count2d(2),start3d(4),count3d(4)
      real*4 :: lcin,laiin,hin,hstd
      integer :: fidlc(KM),fileidlaimax(KM)
      integer :: varidlc(KM),varidlaimax(KM)
      integer :: fidlaiin(12),fileidlcm(12)
      integer :: varidlaiin(12),varidlcm(KM)
      integer :: varidlaiout(KM),countvid

      character*80 :: fileheight,filestd
      character*80 :: filelaiin,filelcm(12)
      character*80 :: filelc, filelaimax(KM)
      integer :: fileidheight,varidheight
      integer :: fileidstd,varidstd,err
      character*80 :: filebs
      integer :: fileidbs,varidbs
      character*60 :: PathFilepre, PathFilepost

      character*80 :: fllaimaxpure,fllcmaxpure
      integer :: fidlaimaxpure,fidlcmaxpure,fidlaiout(12)
      integer :: varidaiimaxpure(18),varidlcmaxpure(18)

      integer, parameter :: N_HT = 19 ! number of height layers input

!      real*4 vf(IM,JM,KM), lai(IM,JM,KM)
      character*80 :: title(KM), title_tmp,title12(12,KM),titlehn(2,KM)

      ! Input values
      ! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
      integer, parameter :: IMn=IM, JMn=JM

      character*(*), parameter :: file_checksum = 
     &     "../lc_lai_ent16/EntMM16_checksum_"//res_out//".ij"

      ! Input values, max, monthly
      real*4 vfn(KM), lain(KM), area
      real*4 a
      real*4 vfnm(12,KM),lainm(12,KM),aream
      real*4 am
      real*4 hmn(KM),hsdn(KM)
      ! Converted values
      real*4 vfc(KM), laic(KM), dvf, s
      character*80 :: titlec(18)
      ! Converted values - monthly
      real*4 vfm(12,KM), laim(12,KM)
      real*4 laimnc(KM), vfmnc(KM),laicropnc(12)

      character*80 :: titlem(12,18)
      ! Converted values - heights
      real*4 vfh(KM),hm(KM), hsd(KM)
      real*4 hmnc, hsdnc
      real*4 vfcnc(IM1km,JM1km,18),laicnc(IM1km,JM1km,18)
      real*4 vfcncout(IM1km,JM1km),laicncout(IM1km,JM1km)
      character*80 :: titleh(2,18) !1-h, 2-hsd
      ! Converted values - crop ext 
      real*4 laicropext, laimcropext(12)
      real*4 hmcropext,hsdcropext
      ! Vars for calculating nocrops
      integer naturalfound, flag, nonaturalcount !if no natural veg in vicinity

      real*4 vf_xx, lai_xx
      real*4 vf_yy, lai_yy
      real*4 LAYER
      character*80 :: title_xx="xx"
      character*80 :: title_yy="yy"
      character*80 :: titlefoo

      ! bs_brightratio = bare soil brightratio
      real*4 :: bs_brightratio, vfc_tmp

      integer i, j, k, io, in, jn, maxpft, kx, m
      real*8 lat
      real*4 foo
      integer N_VEG             ! number of PFTs in output
      integer N_BARE            ! index of bare soil in output
      integer count

      integer :: ncidin,ncidout,varid,status,varidn(18)
      integer :: varids(18)
      character*20 :: inqvarin, inqvarout
      real*4 :: long(IM1km),lati(JM1km)
      character*256 :: fileoutnc
      real*4, dimension(12) :: time
      time = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
            
      !define Lon and Lat
      call calc_lon_lat(IM1km,JM1km,long,lati)
      
      !------------------------------------------------------------------------
      ! OPEN INPUT FILES

      ! lcmax
      do k=1,20
         PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
	 PathFilepost=  EntPFT_files(k)
         filelc = trim(PathFilepre)//trim(PathFilepost)//'_lc.nc'
         err = NF_OPEN(filelc,0,fidlc(k))
         write(*,*) err,fidlc(k)
         inqvarin = EntPFT_lcfields(k)
         err = NF_INQ_VARID(fidlc(k),inqvarin,varidlc(k))
         write(*,*) err,fidlc(k),varidlc(k),EntPFT_lcfields(k)
      enddo

      ! bs ratio
      PathFilepre= '../lc_lai_ent/'
      PathFilepost= 'bs_brightratio.nc'
      filebs = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filebs,0,fileidbs)
      write(*,*) err,filebs
      inqvarin = 'bs_brightratio'
      err = NF_INQ_VARID(fileidbs,inqvarin,varidbs)
      write(*,*) err, 'bs brightratio'
      
      ! lc & laimax file of Ent 17 PFTs: JUST FOR TITLE(K)
      open(1,file="../lc_lai_ent/EntMM_lc_laimax_144x90.bin",
     &     form="unformatted",status="old")
      do k=1,KM 
         read(1) title(k)
         !write(*,*) title(k)
      enddo
      close(1)
      
      do m=1,12
         open(1,file="../lc_lai_ent/EntMM_lc_laimax_144x90.bin",
     &        form="unformatted",status="old")
         do k=1,KM
            read(1) title12(m,k)
            !write(*,*) title12(m,k), KM
         enddo
         close(1)
      enddo

      ! simard heights titles
      open(1,file='../../'//
     &     'data/height/'//
     &      'EntGVSDmosaic17_height_144x90.ij',
     &     form='unformatted',status="old")
      titlehn(1,1) = 'NO WATER LAYER'
      do k=1,N_HT !height
         read(1) titlehn(1,k+1)
         write(*,*) titlehn(1,k+1)
      enddo
      titlehn(2,1) = 'NO WATER LAYER'
      do k=1,N_HT !stdev
         read(1) titlehn(2,k+1)
         write(*,*) titlehn(2,k+1)
      enddo

      !------------------------------------------------------------------------

      do m = 1,2 ! doy files

         ! lai doy
         PathFilepre= '../lc_lai_ent/nc/'
         PathFilepost= 'EntMM_lc_lai_'//trim(DOY(m))//
     &        '_1kmx1km.nc'
         filelaiin = trim(PathFilepre)//trim(PathFilepost)
         err = NF_OPEN(filelaiin,0,fidlaiin(m))
	 write(*,*) err, filelaiin
         inqvarin = 'EntPFT'
         err = NF_INQ_VARID(fidlaiin(m),inqvarin,varidlaiin(m))
         write(*,*) err

        ! CREATE OUTPUT NETCDF FILES
         fileoutnc =
     &        '../lc_lai_ent16/nc/V1km_EntGVSDv1.1_BNU16_lai_'
     &        //DOY(m)//'_pure.nc'
         err = NF_CREATE(fileoutnc,NF_64BIT_OFFSET,fidlaiout(m))
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(fidlaiout(m),'lon',IM1km,dimidx)
         err=NF_DEF_DIM(fidlaiout(m),'lat',JM1km,dimidy)
         err=NF_DEF_VAR(fidlaiout(m),'lon',NF_FLOAT,1,dimidx,varidxlai)
         err=NF_PUT_ATT_TEXT(fidlaiout(m),varidxlai,
     &        'long_name',9,"longitude")
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fidlaiout(m),varidxlai,
     &        'units',12,"degrees east")
         write(*,*) err
         err=NF_DEF_VAR(fidlaiout(m),'lat',NF_FLOAT,1,dimidy,varidylai)
         err=NF_PUT_ATT_TEXT(fidlaiout(m),varidylai,
     &        'long_name',8,"latitude")
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fidlaiout(m),varidylai,
     &        'units',13,"degrees north")
         write(*,*) err
         dd(1)=dimidx
         dd(2)=dimidy
         do k=1,18
            inqvarout = 'lai_'//EntPFT_shorttitle(k)
            err=NF_DEF_VAR(fidlaiout(m),inqvarout,NF_FLOAT,2,dd,
     &           varidlaiout(k))
            write(*,*) err,inqvarout
            err=NF_PUT_ATT_TEXT(fidlaiout(m),varidlaiout(k),
     &           'long_name',45,EntPFT_title(k)//"- LAI")
            write(*,*) err
            err=NF_PUT_ATT_TEXT(fidlaiout(m),varidlaiout(k),
     &           'units',5,"m2/m2")
            write(*,*) err
            err=NF_PUT_ATT_REAL(fidlaiout(m),varidlaiout(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
            write(*,*) err
         enddo
         err=NF_ENDDEF(fidlaiout(m))
         write(*,*) err
         err=NF_PUT_VARA_REAL(fidlaiout(m),varidxlai,1,IM1km,long)
         write(*,*) err
         err=NF_PUT_VARA_REAL(fidlaiout(m),varidylai,1,JM1km,lati)
         write(*,*) err
         
         
         
         do j = 1,JM1km
            
            do i = 1,IM1km
               
               if (i.eq.1) then
                  write(*,*) 'lat', j
               endif
               
               do k = 1,20
               
                  start3d(1) = k
                  start3d(2) = i
                  start3d(3) = j
                  count3d(1) = 1
                  count3d(2) = divx
                  count3d(3) = divy
                  
                  ! get lc max
                  err = NF_GET_VARA_REAL(fidlc(k),varidlc(k),
     &                 start2d,count2d,lcin)
!                 write(*,*) err, 'lc max input',k,start2d(i),start2d(j)
                  vfn(k)=lcin
                  
                  err = NF_GET_VARA_REAL(fidlaiin(m),varidlaiin(m),
     &                 start3d,count3d,laiin)
                  lain(k) = laiin
                  
               enddo
               
               start2d(1) = i
               start2d(2) = j
               count2d(1) = 1
               count2d(2) = 1
               
            ! get bs bright ratio
               err = NF_GET_VARA_REAL(fileidbs,varidbs,
     &              start2d,count2d,bs_brightratio)
!                 write(*,*) err, 'bs bright ratio'
               
                 !Check if mismatch lc or lai values (one is zero and the other not)
!                call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)
            

                !* Convert to GISS 16 pfts format

                  ! first 14 pfts though grass are the same, ignore WATER
                  !  lc laimax
               vfc(1:14) = vfn(2:15)
               laic(1:14) = lain(2:15)
!     titlec(1:14) = title(2:15)
                  !  lc lai monthly
!                  vfm(:,1:14) = vfnm(:,2:15)
!                  laim(:,1:14) = lainm(:,2:15)
!                 do m=1,12
!                    titlem(m,1:14) = title12(m,2:15)
!                 enddo
                  !  heights
!                 vfh(1:14) = vfn(2:15) !Should be the same cover from MODIS.
!                 hm(1:14) = hmn(2:15) 
!                 hsd(1:14) = hsdn(2:15)
!                 titleh(1,1:14) = titlehn(1,2:15)
!                 titleh(2,1:14) = titlehn(2,2:15)
            
                  ! crops
#     ifdef COMBINE_CROPS_C3_C4
            
                 !lc laimax
               a = vfn(16) + vfn(17)
               
               if ( a > 0. ) then
                  laic(15) = (vfn(16)*lain(16)
     &                 + vfn(17)*lain(17)) / a
                  vfc(15) = a
               else
                  laic(15) = 0.
                  vfc(15) = 0.
               endif
                  !lc lai monthly
!                  do m=1,12
!                     am = vfnm(m,16) + vfnm(m,17)
!                     if (am > 0. ) then
!                        laim(m,15) = (vfnm(m,16)*lainm(m,16)
!     &                       + vfnm(m,17)*lainm(m,17)) / am
!                        vfm(m,15) = am
!                     else
!                        laim(m,15) = 0.
!                        vfm(m,15) = 0.
!                     endif
!                  enddo
            !heights - DO NOT AVERAGE. PRESERVE HEIGHTS. LAI will scale density
!            a = vfn(16) + vfn(17) !input cover
!            if ( a > 0. ) then
!               hm(i,j,15) = (vfn(i,j,16)*hmn(i,j,16)
!     &              + vfn(i,j,17)*hmn(i,j,17)) / a
!               if ((hmn(16)>0.).and.(hmn(17)>0.)) then
                  !average if both exist
!                  hm(15) = (vfn(16)*hmn(16)
!     &                 + vfn(17)*hmn(17)) / a
!               else
                  !don't average if only one or none exists
!                  hm(15) = max(hmn(16),hmn(17))
!               endif
               !Sum of squares for sd.  Don't weight if only one or less exists
!               if ((hmn(16)>0.).and.(hmn(17)>0.)) then
!                  hsd(15) = sqrt((vfn(16)*hsdn(16)**2
!     &                 + vfn(17)*hsdn(17)**2) / a)
!               else
!                  hsd(15) = max(hsd(16),hsd(17))
!               endif
!               vfh(15) = a
!            else
!               hm(15) = 0.
!               hsd(15) = 0.
!               vfh(15) = 0.
!            endif
            
!            write(*,*) "Re-doing crops.."
!            titlec(15) = title(16)
!            titlec(15)(1:18) = "15 - crops herb   "
!            do m=1,12
!               titlem(m,15) = title12(m,16)
!               titlem(m,15)(1:18) = "15 - crops herb   "
!            enddo
!            titleh(1,15) = titlehn(1,16)
!            titleh(2,15) = titlehn(2,16)
!           titleh(1,15)(1:18) = "15 - crops herb   "
!            titleh(2,15)(1:18) = "15 - crops herb   "
            
!            write(*,*) titlec(15)
      ! crops woody
               vfc(16) = vfn(18)
               laic(16) = lain(18)
!            titlec(16) = '16 - '//title(18)(6:80)
!            write(*,*) "titlec: "
!            write(*,*) titlec(16)
!                  do m=1,12
!                  vfm(m,16) = vfnm(m,18)
!                  laim(m,16) = lainm(m,18)
!               titlem(m,16) = '16 - '//title12(m,18)(6:80)
!               enddo
!            vfh(16) = vfn(18)
!            titleh(1,16) = '16 - '//titlehn(1,18)(6:80)
!            titleh(2,16) = '16 - '//titlehn(2,18)(6:80)
      ! bare soil
               vfc(17) = vfn(20)
               laic(17) = lain(20)
!            titlec(17) = title(20)
!            vfm(:,17) = vfnm(:,20)
!            laim(:,17) = lainm(:,20)
!            do m=1,12
!               titlem(m,17) = title12(m,20)
!            enddo
!                  vfh(17) = vfn(20)
!                  hm(17) = hmn(20)
!                 hsd(17) = hsdn(20)
!            titleh(1,17) = titlehn(1,20)
!            titleh(2,17) = titlehn(2,20)
               N_VEG = 16
               N_BARE = 17
            
!            titlefoo = 'vfn16'
!            write(92) titlefoo, vfn(16)
!            titlefoo = 'vfn17'
!            write(92) titlefoo, vfn(17)
!            titlefoo = 'Crops 15 after combining C3 and C4'
!            write(92) titlefoo, vfc(15)
               
#     else
      !crops
               vfc(15:17) = vfn(16:18)
               laic(15:17) = lain(16:18)
            
!            titlec(15:17) = title(16:18)
!                  do m=1,12
!                     vfm(m,15:17) = vfnm(m,16:18)
!                     laim(m,15:17) = lainm(m,16:18)
!               titlem(m,15:17) = title12(m,16:18)
!                  enddo
!                  vfh(15:17) = vfn(16:18)
!                  hm(15:17) = hmn(16:18)
!            titleh(1,15:17) = titlehn(1,16:18)
!            titleh(2,15:17) = titlehn(2,16:18)
! bare soil
               vfc(18) = vfn(20)
               laic(18) = lain(20)
!            titlec(18) = title(20)
!            vfm(:,18) = vfnm(:,20)
!            laim(:,18) = lainm(:,20)
!            titlem(:,18) = title12(:,20)
!            vfh(18) = vfh(20)
!            hm(18) = hmn(20)
!            hsd(18) = hsd(20)
!            titleh(1,18) = titlehn(1,20)
!            titleh(2,18) = titlehn(2,20)
               N_VEG = 17
               N_BARE = 18
#     endif
            
 !           do k=1,N_BARE
 !              write(3) titlec(k), vfc(k)
 !           enddo
 !           do k=1,N_BARE
 !              write(3) titlec(k), laic(k)
 !           enddo
            
      ! check if "bare" soil is not bare
               vf_xx = 0.
               if( vfc(N_BARE) > .01 .and. laic(N_BARE) > .5 ) then
                  vf_xx = vfc(N_BARE)
                  lai_xx = laic(N_BARE)
               endif
 !           write(4) title_xx, vf_xx
 !           write(4) title_xx, lai_xx
            
               vf_yy = 0.
               if( vfc(10) > .1 .and. laic(10) < .5 ) then
                  vf_yy = vfc(10)
                  lai_yy = laic(10)
               endif
!            write(4) title_yy, vf_yy
!            write(4) title_yy, lai_yy
            
               vf_yy = 0.
               lai_yy = 0.
               if( vfc(N_BARE) > .1 .and. laic(10) < .01 
     &              .and. laic(9) < .01 .and. laic(11) < .01 
     &              .and. laic(12) < .01 .and. laic(13) < .01 ) then
                  vf_yy = vfc(N_BARE)
                  lai_yy = laic(N_BARE)
               endif
            
!            write(4) title_yy, vf_yy
!            write(4) title_yy, lai_yy
            
            
         !!!! do conversions !!!!
            
         ! convert sparse veg to cold adapted shrub 9 if present
               s = sum(vfc(1:N_BARE))
!                  if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
!               write(*,*) 'ERROR orig:  max and monthly lc different'
!     &               ,s,sum(vfm(1,1:N_BARE))
!               write(*,*) vfc(1:N_BARE)
!               write(*,*) vfm(j,1:N_BARE)
!            endif
               if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15
     &              .and. laic(N_BARE) > .0
     &              .and. vfc(9) > .0 ) then
                  call convert_vf(vfc(N_BARE), laic(N_BARE),
     &                 vfc(9), laic(9), laic(9) )
!     lai >= lai(9)
!                     do m=1,12
!                        call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
!     &                       vfm(m,9),vfm(m,9), vfc(9))
!                     enddo
               
!                     call convert_vfh(
!     &                    vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
!     &                    vfh(9),hm(9),hsd(9), vfc(9))
               endif
            
               s = sum(vfc(1:N_BARE))
!            if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
!               write(*,*) 'ERROR sparse:  max and monthly lc different'
!     &              ,s,sum(vfm(1,1:N_BARE))
!               write(*,*) vfc(1:N_BARE)
!               write(*,*) vfm(1,1:N_BARE)
!            endif

            ! convert sparse veg to arid adapted shrub 10 if present
               if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15
     &              .and. laic(N_BARE) > .0
     &              .and. vfc(10) > .0 ) then
               
                  call convert_vf(vfc(N_BARE), laic(N_BARE),
     &                 vfc(10), laic(10), laic(10) )
                                ! lai >= lai(10)
!                        do m=1,12
!                           call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
!     &                          vfm(m,10),laim(m,10),vfc(10))
!                        enddo
               
!                        call convert_vfh(
!     &                       vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
!     &                       vfh(10),hm(10),hsd(10), vfc(10))
               endif
                     
       ! convert the rest of sparse veg to crop 15 if present
               if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0
     &              .and. vfc(15) > .0 ) then
!             print *, 'Converting spare to crop/bare',i,j,
!     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), vfc(i,j,15)
                  call convert_vf(vfc(N_BARE), laic(N_BARE),
     &                 vfc(15), laic(15), laic(15))
!             print *, 'After conversion:            ',i,j,
!     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), 
!     &            vfc(i,j,15), laic(i,j,15)
!                           do m=1,12
!                              call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
!     &                             vfm(m,15),laim(m,15),vfc(15))
!                           enddo
!                           call convert_vfh(
!     &                          vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
!     &                          vfh(15),hm(15),hsd(15), vfc(15))
               endif
                        
      ! convert the rest of sparse veg to pft with biggest fraction
      ! (if present)
               if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
                  
                  maxpft = maxloc( vfc(1:16), 1 )
!               print *, "max pft is ",maxpft
                  if ( vfc(maxpft) < .0001 ) cycle
                  
                  call convert_vf(vfc(N_BARE), laic(N_BARE),
     &                 vfc(maxpft), laic(maxpft), laic(maxpft))
!                        do m=1,12
!                           call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE)
!     &                          ,vfm(m,maxpft),laim(m,maxpft),vfc(maxpft))
!                        enddo
                        
!                        call convert_vfh(
!     &                       vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
!     &                       vfh(maxpft),hm(maxpft),hsd(maxpft), 
!     &                       vfc(maxpft))
                  
               endif
! convert the rest of sparse veg to arid adapted shrub 10
               if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
                  
                  call convert_vf(vfc(N_BARE), laic(N_BARE),
     &                 vfc(10), laic(10), .0 )
!                        do m=1,12
!                           call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
!     &                          vfm(m,10),laim(m,10),vfc(10))
!                        enddo
               
!                        if (vfc(10) > 0.) then
!                           hm(10) = 2.0 !Check simard.f Set_shrub_height for value!
!                           hsd(10) = 0.
!                        endif
!                        call convert_vfh(
!     &                       vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
!     &                       vfh(10),hm(10),hsd(10), vfc(10))
                  
               endif
               
            
#ifdef SPLIT_BARE_SOIL
               call split_bare_soil(N_VEG,KM,N_BARE
     &              ,bs_brightratio,vfc,laic,res_out)
#endif
            
      ! check titles
!            write(*,*) 'titlec:'
!               do k=1,N_BARE
!               write(*,*) trim(titlec(k))
!               enddo
!            write(*,*) 'titlem:'
!                     do m=1,12
!               write(*,*) MONTH(m)
!                        do k=1,N_BARE
!     write(*,*) trim(titlem(m,k))
!     enddo
!     enddo

               do k=1,18
                  if (laic(k).le.0.) then
                     laic(k) = undef
                  endif
               enddo
               
               ! correct height=undef when land cover>0            
               do k=1,KM
                  if (vfc(k).gt.0.and.laic(k).eq.undef) then
   		      laic(k) = 0.
	          else
	              laic(k) = laic(k)
                  endif
               enddo

               do k=1,18
                  laicnc(i,j,k) = laic(k)
               enddo
               
            enddo
            
         enddo
         
         startB(1)=1
         startB(2)=1
         countB(1)=IM1km
         countB(2)=JM1km

      ! save netcdf lai 
         do k=1,18
            laicncout(:,:) = laicnc(:,:,k)
            err=NF_PUT_VARA_REAL(fidlaiout(m),varidlaiout(k),
     &           startB,countB,laicncout)
            write(*,*) err, 'put'
         enddo
         err = NF_CLOSE(fidlaiout(m))
         err = NF_CLOSE(fidlaiin(m))

         write(*,*) err, 'nf_close in ',fidlaiin(m)
         write(*,*) err, 'nf_close out ',fidlaiin(m)
         
      enddo
      
      
      
      end program convert

