! Trims off tiny fractions and preserves total LAI for the gridcell.
! Alters the cover amount if the LAI is smaller than the main one.
! If LAI is a little biger, might increase.
!
! I don't see where the program is doing that.

! A04 only creates the pure dataset.  It does NOT trim.

! Combine C3 and C4 crops into one PFT, for ModelE
#define COMBINE_CROPS_C3_C4

! Split the bare soil into dark and light, to get the right albedo
#define SPLIT_BARE_SOIL

program convert

use conversions
use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod

implicit none

! -------------------------------------------------------
#ifdef COMBINE_CROPS_C3_C4
integer, parameter :: NREMAPS = 16   ! Size of the remap array for automated index remaping
#else
integer, parameter :: NREMAPS = 18
#endif


! Table remaps input PFTs to output PFTs
! ouptut_c[remap[1]] = input_n[remap[2]]
! Column 1 = index in output arrays for 16 or 17 PFTs
! Column 2 = index in input arrays (water + 17 PFTs + permanent ice + bare/sparse)
real*8, parameter :: remap(2,NREMAPS) = RESHAPE( (/ &
     1,2   & ! -> 1 - evergreen broadleaf early successional       
    ,2,3   & ! -> 2 - evergreen broadleaf late successional        
    ,3,4   & ! -> 3 - evergreen needleleaf early successional      
    ,4,5   & ! -> 4 - evergreen needleleaf late successional       
    ,5,6   & ! -> 5 - cold deciduous broadleaf early successional  
    ,6,7   & ! -> 6 - cold deciduous broadleaf late successional   
    ,7,8   & ! -> 7 - drought deciduous broadleaf                  
    ,8,9   & ! -> 8 - deciduous needleleaf                         
    ,9,10  & ! -> 9 - cold adapted shrub                           
    ,10,11  & ! -> 10 - arid adapted shrub                          
    ,11,12  & ! -> 11 - C3 grass perennial                          
    ,12,13  & ! -> 12 - C4 grass                                    
    ,13,14  & ! -> 13 - C3 grass - annual                           
    ,14,15  & ! -> 14 - arctic C3 grass                             
#ifdef COMBINE_CROPS_C3_C4
     ! 16 PFTs
              ! -> 15 - crops C3+C4 
    ,16,18  & ! -> 16 - crops woody
    ,17,20  & ! -> 17 - bare sparse (bare bright after split)
              ! -> 18 - bare dark (if SPLIT_BARE_SOIL)
#else
     ! 17 PFTs
    ,15,16  & ! -> 15 - crops C3
    ,16,17  & ! -> 16 - crops C4
    ,17,18  & ! -> 17 - crops woody
    ,18,20  & ! -> 18 - bare_sparse (bare bright after split)
              ! -> 19 - bare dark (if SPLIT_BARE_SOIL)
#endif
/), (/ 2,ENTPFTNUM /))

! ------- Specific indices in orginal indexing scheme
integer, parameter :: crops_c3_herb_n = 16
integer, parameter :: crops_c4_herb_n = 17

! ------- Specific indices in destination indexing scheme
integer, parameter :: cold_adapted_shrub_c = 9
integer, parameter :: arid_adapted_shrub_c = 10
#ifdef COMBINE_CROPS_C3_C4
integer, parameter :: crops_herb_c = 15
integer, parameter :: crops_woody_c = 16
integer, parameter:: bare_sparse_c = 17
#    ifdef SPLIT_BARE_SOIL
     integer, parameter :: bare_bright_c = 17
     integer, parameter :: bare_dark_c = 18
#    endif
#else
integer, parameter :: crops_woody_c = 17
integer, parameter:: bare_sparse_c = 18
#    ifdef SPLIT_BARE_SOIL
     integer, parameter :: bare_bright_c = 18
     integer, parameter :: bare_dark_c = 19
#    endif
#endif

integer, parameter :: last_pft_c = crops_woody_c
#ifdef SPLIT_BARE_SOIL
integer, parameter :: ncover_c = bare_dark_c
#else
integer, parameter :: ncover_c = bare_sparse_c
#endif

! Labels of our output covers
type(EntLabel_t), dimension(ncover_c) :: ent_c
! -------------------------------------------------------

integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees

integer, parameter :: IM1km=43200, JM1km=21600
integer, parameter :: IM=1, JM=1, KM=20
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
real*4 vfn(KM)    ! io_lcin
real*4 lain(KM)   ! io_laiin
real*4 hmn(KM)     ! io_simard
real*4 :: bs_brightratio   !are soil brightratio
! Renumbered input values
real*4 vfc(KM)    ! = vfn = io_lcin
real*4 laic(18)   ! = lain = io_lcaiin
real*4 hm(KM)     ! = hmn = io_simard



real*4 area
real*4 a
real*4 vfnm(12,KM),lainm(12,KM),aream
real*4 am
real*4 hsdn(KM)
! Converted values
character*80 :: titlec(18)
! Converted values - monthly
real*4 vfm(12,KM), laim(12,KM)
real*4 laimnc(KM), vfmnc(KM),laicropnc(12)

character*80 :: titlem(12,18)
! Converted values - heights
real*4 vfh(KM), hsd(KM)

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

real*4 :: vfc_tmp

!integer i, j, 
integer :: ii,jj   ! Position in overall NetCDF variable
integer :: k, io, in, jn, maxpft, kx, m
real*8 lat
!real*4 foolc(IM1km,JM1km),foolai(IM1km,JM1km)
integer N_VEG             ! number of PFTs in output
integer N_BARE            ! index of bare soil in output
integer count
real*4 :: val

call init_ent_labels

call chunker%init(IM1km, JM1km, IMH*2,JMH*2, &
    100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO + SIMARD = 42
    120)     ! # files to write >= N_LAIMAX + 3*CHECKSUMS

!------------------------------------------------------------------------
! OPEN INPUT FILES

! lcmax
do k=1,20
    ! USE the land cover we computed in a previous step!!!!!!
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

! Bare Soil Brightness Ratio
call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
    'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)

! Our processed Simard heights
call chunker%nc_open(ioall_simard, LC_LAI_ENT_DIR, &
    'height/', 'EntGVSDmosaic17_height_1kmx1km_lai3g.nc', &
    'SimardHeights', 0)
do k=1,19
    call chunker%nc_reuse_var(ioall_simard, io_simard(k), (/1,1,k/), 'r')
end do

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! CREATE OUTPUT NETCDF FILES

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

! Get our Ent labels
do ri=1,NREMAPS
    ent_c(remaps(1,ri)) = ent20(remaps(2,ri))


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

       ! Compute overall NetCDF index of current cell
       ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
       jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

        do k = 1,KM   ! KM=20
            vfn(k) = io_lcin(k)%buf(ic,jc)   ! LC
            lain(k) = io_laiin(k)%buf(ic,jc)   ! LAImax
        end do

        ! height file - insert dummy WATER layer at beginning to avoid
        ! confusion in numbering with vfc and vfm.
        hmn(1) = 0
        do k=1,19
            hmn(k+1) = io_simard(k)%buf(ic,jc)
        end do

        ! get bs bright ratio
        bs_brightratio = io_bs%buf(ic,jc)

        ! Check if mismatch lc or lai values (one is zero and the other not)
        ! call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)


        ! =============== Convert to GISS 16 pfts format
        !* Convert to GISS 16 pfts format

        do ri=1,NREMAPS
            vfc(remaps(1,ri)) = vfn(remaps(2,ri))
            laic(remaps(1,ri)) = lain(remaps(2,ri))
            hm(remaps(1,ri)) = hmn(remaps(2,ri))
        end do
! TODO: rename hm -> hmc


#ifdef COMBINE_CROPS_C3_C4
        !lc laimax
        a = vfn(crops_c3_herb_n) + vfn(crops_c4_herb_n)
        if ( a > 0. ) then
            laic(crops_herb_c) = (vfn(crops_c3_herb_n)*lain(crops_c3_herb_n) &
              + vfn(crops_c4_herb_n)*lain(crops_c4_herb_n)) / a
            vfc(crops_herb_c) = a
        else
            laic(crops_herb_c) = 0.
            vfc(crops_herb_c) = 0.
        endif
#endif


        ! =========== Partition bare/sparse LAI to actual LC types.
        ! Bare sparse has no vegetation type assigned to it; but it
        ! has non-zero LAI that must be treated.  So we take that
        ! non-zero LAI over sparse veg, assign it to the most likely
        ! vegetation type.  Then we have to assign a cover fraction
        ! for the type, and updated cover fraction for now completely
        ! bare soil.

        !!!! do conversions !!!!

        ! ----------- The conditionals below are mutually exclusive.
        ! convert_vf(), when run, sets laic(bare_sparse_c), which will
        ! cause all subsequent conditional blocks to not run.

        ! --------------
        ! Convert to shrub if there's a small fraction and the shrubs already exist
        ! convert sparse veg to cold adapted shrub (9) if present
        if( vfc(bare_sparse_c) > 0e0 .and. vfc(bare_sparse_c) < .15 &
           .and. laic(bare_sparse_c) > 0e0 &
           .and. vfc(cold_adapted_shrub_c) > 0e0 ) &
        then
            ! Preserve total LAI, but put it all in that vegetation
            ! type.
            call convert_vf( &
                vfc(bare_sparse_c), laic(bare_sparse_c), &
                vfc(cold_adapted_shrub_c), laic(cold_adapted_shrub_c), &
                laic(cold_adapted_shrub_c) )
        endif

        ! convert sparse veg to arid adapted shrub 10 if present
        if( vfc(bare_sparse_c) > 0e0 .and. vfc(bare_sparse_c) < .15 &
           .and. laic(bare_sparse_c) > 0e0 &
           .and. vfc(arid_adapted_shrub_c) > 0e0 ) &
        then
         
            call convert_vf( &
                vfc(bare_sparse_c), laic(bare_sparse_c), &
                vfc(arid_adapted_shrub_c), laic(arid_adapted_shrub_c), &
                laic(arid_adapted_shrub_c) )
        end if
        ! --------------

        ! Convert to crops if they exist (any size fraction)
        ! convert sparse veg to crop 15 if present
        if( vfc(bare_sparse_c) > 0e0 .and. laic(bare_sparse_c) > 0e0 &
              .and. vfc(crops_herb_c) > 0e0 ) &
        then
            call convert_vf( &
                vfc(bare_sparse_c), laic(bare_sparse_c), &
                vfc(crops_herb_c), laic(crops_herb_c), &
                laic(crops_herb_c))
        end if
      
        ! Else... Convert to largest existing PFT
        ! convert sparse veg to pft with biggest fraction
        ! (if present)
        if( vfc(bare_sparse_c) > 0e0 .and. laic(bare_sparse_c) > 0e0 ) then
            maxpft = maxloc( vfc(1:last_pft_c), 1 )
            ! TODO: Why is cycle here?  Check Nancy's original code
            ! if ( vfc(maxpft) < .0001 ) cycle

            call convert_vf(vfc(bare_sparse_c), laic(bare_sparse_c), &
                vfc(maxpft), laic(maxpft), laic(maxpft))
        end if

        ! If none of the above was true, then:
        !  a) It's not a smal lfraction
        !  b) There's no arid shrubs or cold shrubs
        !  c) No crops
        !  d) No clear categorization of other vegetation types in the grid cell
        ! ...so just assign it to arid shrub.
        !
        ! convert sparse veg to arid adapted shrub 10
        if( vfc(bare_sparse_c) > 0e0 .and. laic(bare_sparse_c) > 0e0 ) then
            call convert_vf(vfc(bare_sparse_c), laic(bare_sparse_c), &
                vfc(arid_adapted_shrub_c), laic(arid_adapted_shrub_c), 0e0 )
        end if
        ! -----------------------------------------------------------------

        ! Splits bare soil into bright and dark fractions, so we get
        ! the proper albedo in the GCM
#ifdef SPLIT_BARE_SOIL
        call split_bare_soil(N_VEG,KM,bare_sparse_c, &
            bs_brightratio,vfc,laic, res_out)
#endif


        do k=1,ncover_c
            io_laiout(k)%buf(ic,jc) = laic(k)
        end do

        ! checksum lc & laimax
        io_sum_lc%buf(ic,jc) = 0.
        do k=1,ncover_c
            io_sum_lc%buf(ic,jc) = io_sum_lc%buf(ic,jc) + laic(k)
        end do
    end do  ! ic
    end do  ! jc

    call chunker%write_chunks

end do ! ichunk
end do ! jchunk

call chunker%close_chunks

end program convert
