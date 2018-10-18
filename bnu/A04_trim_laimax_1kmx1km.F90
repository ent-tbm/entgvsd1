
! This program converts EntMM 17 PFTs to Ent 16 PFTs + bright + dark.
! Converts lc, laimax, monthly lc and lai, and Simard heights.
! Original from trim_EntMM_monthly_noht.f, which did not do heights.
! It does NOT interpolate from fine to coarse grid -- this should be done
!  prior to using this program, so input and output are the same resolution.
! To change input/output resolution edit lines below 
! "define input file resolution" and "new (interpolated) values"
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

#define COMBINE_CROPS_C3_C4
#define SPLIT_BARE_SOIL



!------------------------------------------------------------------------
!------------------------------------------------------------------------

!      subroutine omp_set_num_threads(num_threads)
!      integer, intent(in) :: num_threads

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

type(Chunker_t) :: chunker
integer :: ichunk,jchunk, ic,jc
! Input files
type(ChunkIO_t) :: io_lcin(KM), io_laiin(KM), io_bs
! Output files
type(ChunkIO_t) :: ioall_laiout, io_laiout(18)
type(ChunkIO_t) :: ioall_sum, io_sum_lc, io_sum_lai


real*4 :: lcin,laiin,hin,hstd

character*80 :: title_tmp,title12(12,KM),titlehn(2,KM)

! Input values
! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
integer, parameter :: IMn=IM, JMn=JM

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

!integer i, j, 
integer :: k, io, in, jn, maxpft, kx, m
real*8 lat
!real*4 foolc(IM1km,JM1km),foolai(IM1km,JM1km)
integer N_VEG             ! number of PFTs in output
integer N_BARE            ! index of bare soil in output
integer count

call init_ent_labels
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
end do

! laimax
do k=1,20
    call chunker%nc_open(io_laiin(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(ent20(k)%file1)//trim(ent20(k)%file2)//'_lai.nc', &
        trim(ent20(k)%file2), 1)
enddo


!------------------------------------------------------------------------
!------------------------------------------------------------------------
! CREATE OUTPUT NETCDF FILES
! bs ratio
call chunker%nc_create(io_bs, &
    weighting(io_lcin(20)%buf, 1d0, 0d0), &    ! TODO: Scale ???
    '', 'bs_brightratio', 'bs_brightratio', &
    'Bare Soil Bright Ratio', '1', 'BrightRatio')

! laimax_pure
call chunker%nc_create(ioall_laiout, &
    weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
    '16/nc/', 'V1km_EntGVSDv1.1_BNU16_laimax_pure')
do k=1,18
    call chunker%nc_reuse_file(ioall_laiout, io_laiout(k), &
        'lai_'//trim(ent18(k)%file2), trim(ent18(k)%title), &
        'm2 m-2', trim(ent18(k)%title), &
        weighting(chunker%wta1,1d0,0d0))    ! TODO: Weighting???
end do

!  checksum land  laimax
call chunker%nc_create(ioall_sum, &
    weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
    '16/nc/', 'V1km_EntGVSDv1.1_LAI3g16_laimax_pure_checksum')

call chunker%nc_reuse_file(ioall_sum, io_sum_lc, &
    'lc_checksum', 'Checksum of LC', '1', 'checksum - Land Cover', &
    weighting(chunker%wta1, 1d0, 0d0))

call chunker%nc_reuse_file(ioall_sum, io_sum_lai, &
    'lai_checksum', 'Checksum of LAI', 'm2 m-2', 'checksum - LAI', &
    weighting(chunker%wta1, 1d0, 0d0))

call chunker%nc_check('A04_trim_laimax_1kmx1km')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

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

        do k = 1,KM
            lcin = io_lcin(k)%buf(ic,jc)
            vfn(k)=lcin

            ! get lai max
            laiin = io_laiin(k)%buf(ic,jc)
            lain(k) = laiin         
        end do

        ! height file - insert dummy WATER layer at beginning toi avoid
        ! confusion in numbering with vfc and vfm.

        ! get bs bright ratio
        bs_brightratio = io_bs%buf(ic,jc)

        ! Check if mismatch lc or lai values (one is zero and the other not)
        ! call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)
      

        !* Convert to GISS 16 pfts format

        ! first 14 pfts though grass are the same, ignore WATER
        !  lc laimax
        vfc(1:14) = vfn(2:15)
        laic(1:14) = lain(2:15)
        ! titlec(1:14) = title(2:15)
        !  lc lai monthly
        !            do m=1,12
        !               titlem(m,1:14) = title12(m,2:15)
        !            enddo
        !  heights
        !            vfh(1:14) = vfn(2:15) !Should be the same cover from MODIS.
        !            hm(1:14) = hmn(2:15) 
        !            hsd(1:14) = hsdn(2:15)
        !            titleh(1,1:14) = titlehn(1,2:15)
        !            titleh(2,1:14) = titlehn(2,2:15)
      
        ! crops

#ifdef COMBINE_CROPS_C3_C4
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
        !            vfh(16) = vfn(18)
        !            titleh(1,16) = '16 - '//titlehn(1,18)(6:80)
        !            titleh(2,16) = '16 - '//titlehn(2,18)(6:80)
        ! bare soil
        vfc(17) = vfn(20)
        laic(17) = lain(20)
        !            titlec(17) = title(20)
        !            do m=1,12
        !               titlem(m,17) = title12(m,20)
        !            enddo
        !            vfh(17) = vfn(20)
        !            hm(17) = hmn(20)
        !            hsd(17) = hsdn(20)
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
     
#else
        !crops
        vfc(15:17) = vfn(16:18)
        laic(15:17) = lain(16:18)

        !            titlec(15:17) = title(16:18)
        !            vfh(15:17) = vfn(16:18)
        !            hm(15:17) = hmn(16:18)
        !            titleh(1,15:17) = titlehn(1,16:18)
        !            titleh(2,15:17) = titlehn(2,16:18)
        ! bare soil
        vfc(18) = vfn(20)
        laic(18) = lain(20)
        !            titlec(18) = title(20)
        !            titlem(:,18) = title12(:,20)
        !            vfh(18) = vfh(20)
        !            hm(18) = hmn(20)
        !            hsd(18) = hsd(20)
        !            titleh(1,18) = titlehn(1,20)
        !            titleh(2,18) = titlehn(2,20)
        N_VEG = 17
        N_BARE = 18
#endif
      
#if 0
! TODO: Where does write(3) go to anyway?
       do k=1,N_BARE
          write(3) titlec(k), vfc(k)
       enddo
       do k=1,N_BARE
          write(3) titlec(k), laic(k)
       enddo
#endif
      
        ! check if "bare" soil is not bare
        vf_xx = 0.
        if( vfc(N_BARE) > .01 .and. laic(N_BARE) > .5 ) then
            vf_xx = vfc(N_BARE)
            lai_xx = laic(N_BARE)
        endif
#if 0
! TODO: Where does write(4) go to anyway?
        write(4) title_xx, vf_xx
        write(4) title_xx, lai_xx
#endif
      
        vf_yy = 0.
        if( vfc(10) > .1 .and. laic(10) < .5 ) then
            vf_yy = vfc(10)
            lai_yy = laic(10)
        endif
        !            write(4) title_yy, vf_yy
        !            write(4) title_yy, lai_yy
      
        vf_yy = 0.
        lai_yy = 0.
        if( vfc(N_BARE) > .1 .and. laic(10) < .01  &
            .and. laic(9) < .01 .and. laic(11) < .01  &
            .and. laic(12) < .01 .and. laic(13) < .01 ) &
        then

            vf_yy = vfc(N_BARE)
            lai_yy = laic(N_BARE)
        end if
      
        !            write(4) title_yy, vf_yy
        !            write(4) title_yy, lai_yy

        !!!! do conversions !!!!
      
        ! convert sparse veg to cold adapted shrub 9 if present
        s = sum(vfc(1:N_BARE))
        if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
        !               write(*,*) 'ERROR orig:  max and monthly lc different'
        !     &               ,s,sum(vfm(1,1:N_BARE))
        !               write(*,*) vfc(1:N_BARE)
        !               write(*,*) vfm(j,1:N_BARE)
        endif

        if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15 &
           .and. laic(N_BARE) > .0 &
           .and. vfc(9) > .0 ) &
        then
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(9), laic(9), laic(9) )
        !     lai >= lai(9)

        !              call convert_vfh(
        !    &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
        !    &              vfh(9),hm(9),hsd(9), vfc(9))
        endif
      
        s = sum(vfc(1:N_BARE))

        ! convert sparse veg to arid adapted shrub 10 if present
        if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15 &
           .and. laic(N_BARE) > .0 &
           .and. vfc(10) > .0 ) &
        then
         
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(10), laic(10), laic(10) )
                          ! lai >= lai(10)
            !               do m=1,12
            !                  call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
            !     &                 vfm(m,10),laim(m,10),vfc(10))
            !               enddo
                      
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(10),hm(10),hsd(10), vfc(10))
        end if
         
        ! convert the rest of sparse veg to crop 15 if present
        if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 &
              .and. vfc(15) > .0 ) &
        then
            !             print *, 'Converting spare to crop/bare',i,j,
            !     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), vfc(i,j,15)
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(15), laic(15), laic(15))
            !             print *, 'After conversion:            ',i,j,
            !     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), 
            !     &            vfc(i,j,15), laic(i,j,15)
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(15),hm(15),hsd(15), vfc(15))
        end if
      
        ! convert the rest of sparse veg to pft with biggest fraction
        ! (if present)
        if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
            maxpft = maxloc( vfc(1:16), 1 )
            !               print *, "max pft is ",maxpft
            if ( vfc(maxpft) < .0001 ) cycle
            
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(maxpft), laic(maxpft), laic(maxpft))
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(maxpft),hm(maxpft),hsd(maxpft), 
            !     &              vfc(maxpft))
                                 
        end if

        ! convert the rest of sparse veg to arid adapted shrub 10
        if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(10), laic(10), .0 )

            !               if (vfc(10) > 0.) then
            !                  hm(10) = 2.0  !Check simard.f Set_shrub_height for value!
            !                  hsd(10) = 0.
            !               endif
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(10),hm(10),hsd(10), vfc(10))
                                 
        end if

#ifdef SPLIT_BARE_SOIL
        call split_bare_soil(N_VEG,KM,N_BARE &
            ,bs_brightratio,vfc,laic, &
            res_out)
#endif
      
        ! check titles
        !            write(*,*) 'titlec:'
        !            do k=1,N_BARE
        !               write(*,*) trim(titlec(k))
        !            enddo
        !            write(*,*) 'titlem:'

        do k=1,18
            if (laic(k).le.0.) then
                laic(k) = undef
            end if
        end do

        ! correct height=undef when land cover>0            
        do k=1,18
            if (vfc(k).gt.0.and.laic(k).eq.undef) then
                laic(k) = 0.
            else
                laic(k) = laic(k)
            end if
        end do

        do k=1,18
            io_laiout(k)%buf(ic,jc) = laic(k)
            !laicnc(i,j,k) = laic(k)
        end do

        ! checksum lc & laimax
        io_sum_lc%buf(ic,jc) = 0.
        do k=1,18
            io_sum_lc%buf(ic,jc) = io_sum_lc%buf(ic,jc) + laic(k)
        end do
    end do  ! ic
    end do  ! jc

    call chunker%write_chunks

end do ! ichunk
end do ! jchunk

end program convert
