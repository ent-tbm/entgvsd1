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

!subroutine omp_set_num_threads(num_threads)
!integer, intent(in) :: num_threads

!------------------------------------------------------------------------

program convert

use conversions
use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod
    
implicit none

integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees


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
!integer, parameter :: IM=1440, JM=720, KM=20
!integer, parameter :: IM=144, JM=90, KM=20
!character*(*), parameter :: res_in="2.5x2"
!character*(*), parameter :: res_in="144x90"
!character*(*), parameter :: res_in_int="144x90"
!character*(*), parameter :: res_in="1440x720"
!character*(*), parameter :: res_in_int="1440x720"
character*(*), parameter :: res_in="1kmx1km"
character*(*), parameter :: res_out="1kmx1km"

integer, parameter :: divx = IM
integer, parameter :: divy = JM

type(Chunker_t), allocatable :: chunker
integer :: ichunk,jchunk, ic,jc
! ------ Input Files
type(ChunkIO_t) :: io_lcin(KM), io_bs
type(ChunkIO_t) :: ioall_laiin, io_laiin(KM)
! ------ Output files
type(ChunkIO_t) :: ioall_laiout, io_laiout(18)

real*4 :: lcin,laiin,hin,hstd

integer, parameter :: N_HT = 19 ! number of height layers input

!real*4 vf(IM,JM,KM), lai(IM,JM,KM)
character*80 :: title(KM), title_tmp,title12(12,KM),titlehn(2,KM)

! Input values
! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
integer, parameter :: IMn=IM, JMn=JM

character*(*), parameter :: file_checksum =  &
     "../lc_lai_ent16/EntMM16_checksum_"//res_out//".ij"

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
      

call init_ent_labels

do m = 1,12 ! monthly files
    print *,' ======================= MONTH m = ',m
    allocate(chunker)
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 100, 120)

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! lcmax
    do k=1,20
        ! TODO: LAI3g????
        ! PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
        call chunker%nc_open(io_lcin(k), LC_LAI_ENT_DIR, &
            'EntMM_lc_laimax_1kmx1km/', &
            trim(ent20(k)%file1)//trim(ent20(k)%file2)//'_lc.nc', &
            trim(ent20(k)%file2), 1)
    enddo

    ! bs ratio
    call chunker%nc_open(io_bs, LC_LAI_ENT_DIR, &
        '', 'bs_brightratio.nc', 'bs_brightratio', 1)


    titlehn(1,1) = 'NO WATER LAYER'
    titlehn(2,1) = 'NO WATER LAYER'
    do k=2,20
        titlehn(1,k) = trim(ent20(k)%title)//' height'
        titlehn(2,k) = trim(ent20(k)%title)//' stdev'
    end do

    !------------------------------------------------------------------------
    call chunker%nc_open(ioall_laiin, LC_LAI_ENT_DIR, &
        'nc/', 'EntMM_lc_lai_'//trim(MONTH(m))//'_1kmx1km.nc', 'EntPFT', 0)

    do k=1,20
        call chunker%nc_reuse_var(ioall_laiin, io_laiin(k), &
            (/1,1,k/), 'r', weighting(chunker%wta1,1d0,0d0))
    end do

    ! CREATE OUTPUT NETCDF FILES
    call chunker%nc_create(ioall_laiout, &
        weighting(chunker%wta1,1d0,0d0), &
        '16/nc/', 'V1km_EntGVSDv1.1_BNU16_lai_'//MONTH(m)//'_pure.nc')

    do k=1,18
        call chunker%nc_reuse_file(ioall_laiout, io_laiout(k), &
            'lai_'//trim(ent18(k)%file2), &
            trim(ent18(k)%title), 'm2 m-2', &
            trim(ent18(k)%title), &
            weighting(chunker%wta1,1d0,0d0))   ! TODO: What convert LC to 18-class system for weighting???

    end do


    call chunker%nc_check('A06_trim_lc_laimax_laimonth')
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    !-----------------------------------------------------------------
    !     Loop for every grid point

    ! Use these loop bounds for testing...
    ! it chooses a land area in Asia
#ifdef ENTGVSD_DEBUG
    do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
    do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
    do jchunk = 1,nchunk(2)
    do ichunk = 1,nchunk(1)
#endif

        call chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)
   
            do k = 1,20

                ! get lc max
                lcin = io_lcin(k)%buf(ic,jc)
                vfn(k)=lcin


                laiin = io_laiin(k)%buf(ic,jc)
                lain(k) = laiin
            
             end do
         
            ! get bs bright ratio
            bs_brightratio = io_bs%buf(ic,jc)
         
            ! Check if mismatch lc or lai values (one is zero and the other not)
            ! call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)
     
            !* Convert to GISS 16 pfts format

            ! first 14 pfts though grass are the same, ignore WATER
            !  lc laimax
            vfc(1:14) = vfn(2:15)
            laic(1:14) = lain(2:15)

            !     titlec(1:14) = title(2:15)
                  !  lc lai monthly
            !      vfm(:,1:14) = vfnm(:,2:15)
            !      laim(:,1:14) = lainm(:,2:15)
            !     do m=1,12
            !        titlem(m,1:14) = title12(m,2:15)
            !     enddo
                  !  heights
            !     vfh(1:14) = vfn(2:15) !Should be the same cover from MODIS.
            !     hm(1:14) = hmn(2:15) 
            !     hsd(1:14) = hsdn(2:15)
            !     titleh(1,1:14) = titlehn(1,2:15)
            !     titleh(2,1:14) = titlehn(2,2:15)
      
            ! crops
#     ifdef COMBINE_CROPS_C3_C4
      
            !lc laimax
            a = vfn(16) + vfn(17)
         
            if ( a > 0. ) then
                 laic(15) = (vfn(16)*lain(16) &
                     + vfn(17)*lain(17)) / a
                 vfc(15) = a
            else
                 laic(15) = 0.
                 vfc(15) = 0.
            endif

            !lc lai monthly
            !            do m=1,12
            !               am = vfnm(m,16) + vfnm(m,17)
            !               if (am > 0. ) then
            !                  laim(m,15) = (vfnm(m,16)*lainm(m,16)
            !     &                       + vfnm(m,17)*lainm(m,17)) / am
            !                  vfm(m,15) = am
            !               else
            !                  laim(m,15) = 0.
            !                  vfm(m,15) = 0.
            !               endif
            !            enddo
                  !heights - DO NOT AVERAGE. PRESERVE HEIGHTS. LAI will scale density
            !      a = vfn(16) + vfn(17) !input cover
            !      if ( a > 0. ) then
            !         hm(i,j,15) = (vfn(i,j,16)*hmn(i,j,16)
            !     &              + vfn(i,j,17)*hmn(i,j,17)) / a
            !         if ((hmn(16)>0.).and.(hmn(17)>0.)) then
                        !average if both exist
            !            hm(15) = (vfn(16)*hmn(16)
            !     &                 + vfn(17)*hmn(17)) / a
            !         else
                        !don't average if only one or none exists
            !            hm(15) = max(hmn(16),hmn(17))
            !         endif
                     !Sum of squares for sd.  Don't weight if only one or less exists
            !         if ((hmn(16)>0.).and.(hmn(17)>0.)) then
            !            hsd(15) = sqrt((vfn(16)*hsdn(16)**2
            !     &                 + vfn(17)*hsdn(17)**2) / a)
            !         else
            !            hsd(15) = max(hsd(16),hsd(17))
            !         endif
            !         vfh(15) = a
            !      else
            !         hm(15) = 0.
            !         hsd(15) = 0.
            !         vfh(15) = 0.
            !      endif
                  
            !      write(*,*) "Re-doing crops.."
            !      titlec(15) = title(16)
            !      titlec(15)(1:18) = "15 - crops herb   "
            !      do m=1,12
            !         titlem(m,15) = title12(m,16)
            !         titlem(m,15)(1:18) = "15 - crops herb   "
            !      enddo
            !      titleh(1,15) = titlehn(1,16)
            !      titleh(2,15) = titlehn(2,16)
            !     titleh(1,15)(1:18) = "15 - crops herb   "
            !      titleh(2,15)(1:18) = "15 - crops herb   "
                  
            !      write(*,*) titlec(15)
            ! crops woody
            vfc(16) = vfn(18)
            laic(16) = lain(18)
            !      titlec(16) = '16 - '//title(18)(6:80)
            !      write(*,*) "titlec: "
            !      write(*,*) titlec(16)
            !            do m=1,12
            !            vfm(m,16) = vfnm(m,18)
            !            laim(m,16) = lainm(m,18)
            !         titlem(m,16) = '16 - '//title12(m,18)(6:80)
            !         enddo
            !      vfh(16) = vfn(18)
            !      titleh(1,16) = '16 - '//titlehn(1,18)(6:80)
            !      titleh(2,16) = '16 - '//titlehn(2,18)(6:80)
            ! bare soil
            vfc(17) = vfn(20)
            laic(17) = lain(20)
            !      titlec(17) = title(20)
            !      vfm(:,17) = vfnm(:,20)
            !      laim(:,17) = lainm(:,20)
            !      do m=1,12
            !         titlem(m,17) = title12(m,20)
            !      enddo
            !            vfh(17) = vfn(20)
            !            hm(17) = hmn(20)
            !           hsd(17) = hsdn(20)
            !      titleh(1,17) = titlehn(1,20)
            !      titleh(2,17) = titlehn(2,20)
            N_VEG = 16
            N_BARE = 17
      
            !      titlefoo = 'vfn16'
            !      write(92) titlefoo, vfn(16)
            !      titlefoo = 'vfn17'
            !      write(92) titlefoo, vfn(17)
            !      titlefoo = 'Crops 15 after combining C3 and C4'
            !      write(92) titlefoo, vfc(15)

#     else
            !crops
            vfc(15:17) = vfn(16:18)
            laic(15:17) = lain(16:18)

            !      titlec(15:17) = title(16:18)
            !            do m=1,12
            !               vfm(m,15:17) = vfnm(m,16:18)
            !               laim(m,15:17) = lainm(m,16:18)
            !         titlem(m,15:17) = title12(m,16:18)
            !            enddo
            !            vfh(15:17) = vfn(16:18)
            !            hm(15:17) = hmn(16:18)
            !      titleh(1,15:17) = titlehn(1,16:18)
            !      titleh(2,15:17) = titlehn(2,16:18)
            ! bare soil
            vfc(18) = vfn(20)
            laic(18) = lain(20)
            !      titlec(18) = title(20)
            !      vfm(:,18) = vfnm(:,20)
            !      laim(:,18) = lainm(:,20)
            !      titlem(:,18) = title12(:,20)
            !      vfh(18) = vfh(20)
            !      hm(18) = hmn(20)
            !      hsd(18) = hsd(20)
            !      titleh(1,18) = titlehn(1,20)
            !      titleh(2,18) = titlehn(2,20)
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
            !      write(4) title_yy, vf_yy
            !      write(4) title_yy, lai_yy
          
            vf_yy = 0.
            lai_yy = 0.
            if( vfc(N_BARE) > .1 .and. laic(10) < .01  &
                  .and. laic(9) < .01 .and. laic(11) < .01  &
                  .and. laic(12) < .01 .and. laic(13) < .01 ) then
                vf_yy = vfc(N_BARE)
                lai_yy = laic(N_BARE)
            endif
          
            !      write(4) title_yy, vf_yy
            !      write(4) title_yy, lai_yy
          
          
            !!!! do conversions !!!!

            ! convert sparse veg to cold adapted shrub 9 if present
            s = sum(vfc(1:N_BARE))
            !            if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
            !         write(*,*) 'ERROR orig:  max and monthly lc different'
            !     &               ,s,sum(vfm(1,1:N_BARE))
            !         write(*,*) vfc(1:N_BARE)
            !         write(*,*) vfm(j,1:N_BARE)
            !      endif
            if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15 &
                 .and. laic(N_BARE) > .0 &
                 .and. vfc(9) > .0 ) then
               call convert_vf(vfc(N_BARE), laic(N_BARE), &
                    vfc(9), laic(9), laic(9) )
                 !     lai >= lai(9)
                 !               do m=1,12
                 !                  call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
                 !     &                       vfm(m,9),vfm(m,9), vfc(9))
                 !               enddo
                          
                 !               call convert_vfh(
                 !     &                    vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
                 !     &                    vfh(9),hm(9),hsd(9), vfc(9))
            end if

      
            s = sum(vfc(1:N_BARE))
            !      if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
            !         write(*,*) 'ERROR sparse:  max and monthly lc different'
            !     &              ,s,sum(vfm(1,1:N_BARE))
            !         write(*,*) vfc(1:N_BARE)
            !         write(*,*) vfm(1,1:N_BARE)
            !      endif

            ! convert sparse veg to arid adapted shrub 10 if present
            if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15 &
                .and. laic(N_BARE) > .0 &
                .and. vfc(10) > .0 ) &
            then
                call convert_vf(vfc(N_BARE), laic(N_BARE), &
                     vfc(10), laic(10), laic(10) )
                          ! lai >= lai(10)
                !                  do m=1,12
                !                     call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
                !     &                          vfm(m,10),laim(m,10),vfc(10))
                !                  enddo
                         
                !                  call convert_vfh(
                !     &                       vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
                !     &                       vfh(10),hm(10),hsd(10), vfc(10))
            end if
               
            ! convert the rest of sparse veg to crop 15 if present
            if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 &
                .and. vfc(15) > .0 ) &
            then
                !       print *, 'Converting spare to crop/bare',i,j,
                !     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), vfc(i,j,15)
                            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                                 vfc(15), laic(15), laic(15))
                !       print *, 'After conversion:            ',i,j,
                !     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), 
                !     &            vfc(i,j,15), laic(i,j,15)
                !                     do m=1,12
                !                        call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
                !     &                             vfm(m,15),laim(m,15),vfc(15))
                !                     enddo
                !                     call convert_vfh(
                !     &                          vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
                !     &                          vfh(15),hm(15),hsd(15), vfc(15))
            end if

            ! convert the rest of sparse veg to pft with biggest fraction
            ! (if present)
            if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
            
                maxpft = maxloc( vfc(1:16), 1 )
                ! print *, "max pft is ",maxpft
                if ( vfc(maxpft) < .0001 ) cycle

                call convert_vf(vfc(N_BARE), laic(N_BARE), &
                    vfc(maxpft), laic(maxpft), laic(maxpft))
                !    do m=1,12
                !       call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE) &
                !                  ,vfm(m,maxpft),laim(m,maxpft),vfc(maxpft))
                !    enddo
                    
                !    call convert_vfh( &
                !               vfh(N_BARE),hm(N_BARE),hsd(N_BARE), &
                !               vfh(maxpft),hm(maxpft),hsd(maxpft), &
                !               vfc(maxpft))
                                                    
            end if

            ! convert the rest of sparse veg to arid adapted shrub 10
            if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
                call convert_vf(vfc(N_BARE), laic(N_BARE), &
                    vfc(10), laic(10), .0 )
                !    do m=1,12
                !       call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
                !                  vfm(m,10),laim(m,10),vfc(10))
                !    enddo
                  
                !    if (vfc(10) > 0.) then
                !       hm(10) = 2.0 !Check simard.f Set_shrub_height for value!
                !       hsd(10) = 0.
                !    endif
                !    call convert_vfh(
                !               vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
                !               vfh(10),hm(10),hsd(10), vfc(10))
                                            
            end if

#ifdef SPLIT_BARE_SOIL
            call split_bare_soil(N_VEG,KM,N_BARE, &
                bs_brightratio,vfc,laic,res_out)
#endif
      
            ! check titles
            !      write(*,*) 'titlec:'
            !         do k=1,N_BARE
            !         write(*,*) trim(titlec(k))
            !         enddo
            !      write(*,*) 'titlem:'
            !               do m=1,12
            !         write(*,*) MONTH(m)
            !                  do k=1,N_BARE
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
               io_laiout(k)%buf(ic,jc) = laic(k)
               ! laicnc(i,j,k) = laic(k)
           enddo

        end do
        end do

        call chunker%write_chunks

    end do
    end do

    call chunker%close_chunks
    deallocate(chunker)
end do    ! m=1,2 (for MONTH)

end program convert

