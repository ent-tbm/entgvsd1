!     modis_entpftg.f
!     Nancy Kiang, updated by Carlo Montes (August 2015)
!     Process LAI3g, Monfreda, CRU and GPCC climate files into Ent PFT
!     cover and LAI annual maximum.
!     Set longin and latin for 1km resolution and
!     edit file names.
!     Modified from modis_entpftg.f, which had used the MODIS PFT crop partition
!     of "cereal crop" and "broadleaf crop" to estimate C3 and C4 crops, 
!     but correlation is rough.  Here, sum all MODIS crop cover, and use
!     Monfreda to partition MODIS total crop cover into C3, C4, and growth form.
!     8/26/13 - Added FIX_MODIS29_SOUTHPOLE_BUG to fix error in MODIS29.
!     9/5/13 - Added FIX_MODIS29_NORTHPOLE_BUG to fix error in MODIS29.

! This produces the following checksums:
!
!# Checksom on original MODIS LAI??  LC??
!io_checksum = EntMM29lc_lai_for_1kmx1km
!   sum_{1...MODIS28 LCLASS} partit_io(k) + io_waterpart
!    partit_io = 2004/PART_SUB_1km_2004_geo.PARTITION/...
!    io_waterpart = <same>/water
!
!io_checksum2 = EntLandcover_check_sum_Jun_1kmx1km
!    sum_{1...ENTPFTNUM} ENTPFTLC(k) + WATERLC
!    NOTE: ENTPFTLC(k) == contents of file io_lc(k) = EntMM_lc_laimax_1kmx1km/...LAI
!
!io_checksum3 = EntLAI_check_sum_Jun_1kmx1km
!    sum_{1..ENTPFTNUM} ENTPFTLC(k) * ENTPFTLAIMAX(k) + WATERLC*WATERLAI
!
!# Number of non-zero PFTs in each gridcell
!io_npftgrid = EntPFTs_percell_check_sum_Jun_1kmx1km
!
!# LC of dominant PFT
!io_dompftlc = EntdominantPFT_LC_check_sum_Jun_1kmx1km
!
!# Index of dominant PFT
!io_dompft = EntdominantPFT_check_sum_Jun_1kmx1mk




module modis_ent_mod
    use ent_labels_mod
    use geom_mod

implicit none
!     private
!     save

public ENTPFTLC,ENTPFTLAIMAX !, ENTCOVSUM
public Zero_ENTPFT
public Set_pft, Set_Shrubtype, Set_Grasstype,Set_Broadleaftype
public Set_Woodysavannasshrub_miscat
public FIX_MODIS29_NORTHPOLE_BUG
public FIX_MODIS29_SOUTHPOLE_BUG

integer, parameter :: latout = 1
integer :: im, jm

integer, parameter :: LCLASS = 28 !Excludes WATER in 29 MODIS layers.
real*4 :: ENTPFTLC(NENT20)
real*4 :: ENTPFTLAIMAX(NENT20)
real*4 :: ENTPFTLC_NC(NENT20)
real*4 :: ENTPFTLAIMAX_NC(NENT20)

contains

!------------------------------------------------------------------------------------------
real*4 function latdeg(latj, JM)
!     Convert pixel j to latitude degrees (-90 to +90), center of box.
!     See GEOM_B.f
!     GISS GCM counts grid box 1 at south pole, 1/2 box at poles for 4x5
implicit none
    integer, intent(IN) :: latj, JM

    if (JM.eq.JM4X5) then     !72x46
       latdeg = -90. + (latj-1)*4.
    else                      !anything finer resolution
       latdeg = 180./JM*latj - 90. - 180./JM + 180./JM/2.
    endif
end function latdeg
!------------------------------------------------------------------------------------------

subroutine Zero_ENTPFT
ENTPFTLC = 0.  
ENTPFTLAIMAX = 0. 
end subroutine Zero_ENTPFT

subroutine Set_pft(pft,LC_IN,LAI_IN)
implicit none
    integer, intent(IN) :: pft
    real*4, intent(IN) :: LC_IN
    real*4, intent(IN) :: LAI_IN

    ENTPFTLC(pft) = ENTPFTLC(pft) + LC_IN
    ENTPFTLAIMAX(pft) = ENTPFTLAIMAX(pft)  &
         + LC_IN*LAI_IN
end subroutine Set_pft
!------------------------------------------------------------------------------
subroutine Set_Shrubtype(MATEMP,Pmave,LC_IN,LAI_IN)  
implicit none
    real*4, intent(IN) :: MATEMP,Pmave,LC_IN, LAI_IN
    !------

    real*4 :: pft_cold_shrub
    real*4 :: pft_arid_shrub

    pft_cold_shrub = 0.
    pft_arid_shrub = 0.

    if (MATEMP.lt.278.15) then !5 C cut-off
       pft_cold_shrub = LC_IN
    else
       pft_arid_shrub= LC_IN
    endif

    call Set_pft(COLD_SHRUB, pft_cold_shrub,LAI_IN)
    call Set_pft(ARID_SHRUB, pft_arid_shrub,LAI_IN)
end subroutine Set_Shrubtype
!------------------------------------------------------------------------------------------

subroutine Set_Woodysavannasshrub_miscat( &
     C4CLIMFRAC,MATEMP,LC_IN,LAI_IN,j)
!     Feng 11 Woody savannas-shrub: boreal-> 4-evergreen needleleaf, 
!     8-deciduous needleleaf
!     Feng 11 Woody savannas-shrub: rest of world -> 10-arid adapted shrub.
!     Fixes MODIS PFT product miscategorization of Woody savannas->shrub
!     in boreal zone.  The error results in shrub being then classified
!     as Ent 9 low cold-adapted tundra shrubs where there are actually 
!     short trees. Uses C4CLIMFRAC and MATEMP as masks to pick out boreal.
!     Makes boreal zone mix of dominant evergreen needleleaf then deciduous
!     needleleaf where this "shrub" cover is high, all deciduous needleleaf
!     where the cover is medium (open and wetland slopes), and all
!     evergreen needleleaf in very sparse areas.  
!     The rest of the world, assign 10-arid-adapted shrub.
implicit none
    real*4 :: C4CLIMFRAC,MATEMP,LC_IN, LAI_IN

    !------
    integer :: i,j
    real*4 pft_ever_nd_late ! evergreen needleleaf late successional
    real*4 pft_decid_nd     ! deciduous needleleaf
    real*4 pft_arid_shrub   ! arid adapted shrub

    pft_ever_nd_late = 0.
    pft_decid_nd = 0.
    pft_arid_shrub = 0.

    if ((C4CLIMFRAC.le.0.1).and.(latdeg(j,JM1km).gt.45.)) then
       if (MATEMP.ge.(-5.0+273.15)) then !Western Europe boreal
          pft_ever_nd_late = LC_IN
       else                !Russia/Siberia boreal
          if (LC_IN.ge.0.7) then !Closed forest shades out Larix
             pft_ever_nd_late = LC_IN 
          elseif ((LC_IN.lt.0.7).and.(LC_IN.ge.0.3)) &
                  then     !Mix Pinus and Larix
             pft_ever_nd_late = LC_IN - 0.15
             pft_decid_nd = 0.15
          elseif ((LC_IN.lt.0.3).and.(LC_IN.ge.0.15)) &
                  then     !Larix in open/wetlands areas
             pft_decid_nd = LC_IN
          else             !Pinus in sparse areas
             pft_ever_nd_late = LC_IN
          endif
       endif
    else                   !Rest of the world arid shrub
       pft_arid_shrub = LC_IN
    endif

    call Set_pft(EVER_ND_LATE, pft_ever_nd_late, LAI_IN)
    call Set_pft(DECID_ND,     pft_decid_nd,     LAI_IN)
    call Set_pft(ARID_SHRUB,   pft_arid_shrub,   LAI_IN)

end subroutine Set_Woodysavannasshrub_miscat
!------------------------------------------------------------------------------------------

subroutine Set_Grasstype(C4CLIMFRAC,MATEMP,Pdry,ClimMedit, &
     LC_IN,LAI_IN)
implicit none
    real*4, intent(IN) :: C4CLIMFRAC,MATEMP,Pdry
    real*4, intent(IN) :: ClimMedit,LC_IN,LAI_IN
    !--------- Local Vars
    real*4 :: C3

    real*4 :: pft_c3_grass_per     ! C3 grass perennial
    real*4 :: pft_c3_grass_ann  ! C3 grass annual
    real*4 :: pft_c3_grass_arct    ! C3 grass arctic

    pft_c3_grass_per = 0.
    pft_c3_grass_ann = 0.
    pft_c3_grass_arct = 0.

    !     * C4 grass
    call Set_pft(C4_GRASS, LC_IN*C4CLIMFRAC,LAI_IN)

    !     * C3 grass
    C3 = 1. - C4CLIMFRAC
    if (MATEMP.lt.270) then
       pft_c3_grass_arct = LC_IN*C3
    elseif (((Pdry.ge.3.).and.(MATEMP.ge.278.)).or. &
            (ClimMedit.eq.1.)) then
       pft_c3_grass_ann = LC_IN*C3
    else
       pft_c3_grass_per = LC_IN*C3
    endif

    call Set_pft(C3_GRASS_PER,  pft_c3_grass_per,  LAI_IN)
    call Set_pft(C3_GRASS_ANN,  pft_c3_grass_ann,  LAI_IN)
    call Set_pft(C3_GRASS_ARCT, pft_c3_grass_arct, LAI_IN)

end subroutine Set_Grasstype

subroutine Set_Broadleaftype(MATEMP,Pdry,C4climfrac,Tcold &
     ,Pmave,ClimMedit,LC_IN,LAI_IN)
implicit none
    real*4,intent(in) :: MATEMP,Pdry,C4climfrac
    real*4,intent(in) :: Tcold,Pmave,ClimMedit
    real*4,intent(in) :: LC_IN,LAI_IN
    !------- Local Vars
    real*4 :: pft_cold_br_late  ! cold deciduous broadleaf late successional
    real*4 :: pft_drough_br     ! drought deciduous broadleaf

    pft_cold_br_late = 0.
    pft_drough_br = 0.

    if ((C4climfrac.lt.0.8).and. &
         (MATEMP.lt.(18.0+273.15)).and. &
         (Tcold.le.8.).and. &
         (ClimMedit.eq.0.)) &
         then
       pft_cold_br_late = LC_IN
    else
       pft_drough_br = LC_IN
    endif

    call Set_pft(COLD_BR_LATE, pft_cold_br_late,LAI_IN) !cold-deciduous
    call Set_pft(DROUGHT_BR, pft_drough_br,LAI_IN) !drought-deciduous

end subroutine Set_Broadleaftype

subroutine FIX_MODIS29_NORTHPOLE_BUG(LC_IN,LAI_IN,lat)
!     8/26/13 MODIS29 bug fix:  
!     "18.Croplands-->Cereal crop Cover" has non-zero
!     cover at North Pole where it should be "27. Permanent Ice."
!     Assign to Ent17-18Permanent Ice
implicit none
    real*4,intent(inout) :: LC_IN !Input: "18.Croplands-->Cereal crop Cover"
    real*4 :: LAI_IN     !Input: "18.Croplands-->Cereal crop Cover"
    real*4 :: lat
    integer :: northpole

    northpole = 79.5

    if ((lat.gt.79.5).and.(LC_IN.gt.0.)) then
       ! Wrong crops at north pole
       ENTPFTLC(SNOW_ICE) =  ENTPFTLC(SNOW_ICE) + LC_IN
       ENTPFTLAIMAX(SNOW_ICE) = ENTPFTLAIMAX(SNOW_ICE)  &
            + 0.        !If LAI_IN not zero, then MODIS is getting algae...
       LC_IN = 0.             !Don't let later calcs usu wrong LC_IN
    endif
end subroutine FIX_MODIS29_NORTHPOLE_BUG

subroutine FIX_MODIS29_SOUTHPOLE_BUG(LC_IN,LAI_IN,lat)
!     8/26/13 MODIS29 bug fix:  
!     "18.Croplands-->Cereal crop Cover" has non-zero
!     cover at South Pole corners where it should be "27. Permanent Ice."
!     Assign to Ent17-18Permanent Ice
implicit none
    real*4,intent(inout) :: LC_IN !Should be "18.Croplands-->Cereal crop Cover"
    real*4 :: LAI_IN    !Should be "18.Croplands-->Cereal crop Cover"
    real*4 :: lat
    integer :: southpole

    southpole = -82.5

    if ((lat.lt.-82.5).and.(LC_IN>0.)) then
       ! Wrong crops at south pole
       ENTPFTLC(SNOW_ICE) =  ENTPFTLC(SNOW_ICE) + LC_IN
       ENTPFTLAIMAX(SNOW_ICE) = ENTPFTLAIMAX(SNOW_ICE)  &
            + 0.              !If it's not zero, then MODIS is getting algae...
       LC_IN = 0.             !Don't let later calcs usu wrong LC_IN
    endif
end subroutine FIX_MODIS29_SOUTHPOLE_BUG

end module modis_ent_mod
!------------------------------------------------------------------------
program modis_ent
!     Read in MODIS GISS layer and Monfreda 0.5x0.5,1x1,or 4x5 degree files, 
!     and convert veg types to Ent PFTs as GISS layers.
!     To convert to Ent data structures with mixed canopies, use program
!     in Ent repository.
!     Output is files prefixed "EntMM" for Ent-MODIS-Monfreda.
use modis_ent_mod
use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod

implicit none

character*80 :: LAICASE

character*80 :: TITLE, TITLECHECK
character*256 :: file_crops, file_C4clim, file_clim
character*256 :: file_C4clima, file_C4climb
character*256 :: file_clima, file_climb
character*256 :: file_modis_pre !Path prefix for monthly MODIS files
character*256 :: file_modis_sfx !Suffix for monthy MODIS files
character*256 :: file_modis

character*256 :: file_EntMM, file_EntMMA
character*256 :: filecheck
character*256 :: PathMonfreda, PathCropFile
character*10 :: RESOUT
character*50 :: PATHEnt, PATHfile

real*4 :: MODIS29(1+LCLASS) !WATER + 28 LCLASS

real*4 :: LAIMAX

real*4 :: LIN

!To avoid confusion in array indexing, water is read as separate layer
real*4 :: WATERLC
real*4 :: CROPSHERBNORM !Monfreda herb fraction of crop cover
real*4 :: CROPSSHRUBNORM !Monfreda shrub fraction of crop cover
real*4 :: CROPSTREENORM !Monfred tree fraction of crop cover
real*4 :: CROPSC4HERBFRAC !Monfreda C4 fraction of herb crop cover
real*4 :: CROPSTOT !Sum of MODIS crop cover
real*4 :: CROPSC3HERB != CROPSHERBNORM*(1-C4HERBFRAC)*CROPSTOT
real*4 :: CROPSC4HERB != CROPSHERBNORM*C4HERBFRAC*CROPSTOT
real*4 :: C4CLIMFRAC !Used to partition grasses
real*4 :: Tcold  
real*4 :: Pdry        
real*4 :: MAT
real*4 :: TCinave
real*4 :: Pmave !Winds up not being used
real*4 :: ClimMedit !ClimMediterranean climate
real*4 :: CHECKSUM
integer :: NPFTGRID  
real*4 :: DOMPFTLC  
integer :: DOMPFT

real*4, allocatable :: lon(:),lat(:)

integer :: i, j, k, f, m, p
real*4 :: diff
integer :: ncidin,ncidout,varid,status

integer :: err,fileid,fileidin,fileidout,dimidx,dimidy,dimidz
integer :: dd(4),varidx
integer :: varidy,varidz,myvar
integer :: startB(2),countB(2)
integer :: startX(1),startY(1),countX(1),countY(1)
integer :: lenx,leny,lenz
type(Chunker_t) :: chunker
integer :: jchunk, ichunk    ! Index of current chunk
integer :: jc, ic    ! Index WITHIN current chunk
integer :: jj, ii            ! Index in full space

!      real*4 :: inbuf(1,1)  ! Buffer reading NetCDF

type(ChunkIO_t) :: io_laiin,io_04crops, io_05crops
type(ChunkIO_t) :: io_06crops, io_04cropsm, io_C4norm
type(ChunkIO_t) :: io_Tcold, io_Pdry, io_Pmave, io_TCinave
type(ChunkIO_t) :: io_CMedit
type(ChunkIO_t) :: io_waterpart,io_checksum
type(ChunkIO_t) :: io_wateroutA,io_checksum2
type(ChunkIO_t) :: io_npftgrid,io_dompftlc
type(ChunkIO_t) :: io_dompft

real*4 :: WATERLAI

type(ChunkIO_t) :: partit_io(28)
type(ChunkIO_t) :: ioall_lc, io_lc(NENT20)
type(ChunkIO_t) :: ioall_laiout, io_laiout(NENT20)
type(ChunkIO_t) :: ioall_laicheck, io_laicheck(NENT20)
#ifdef COMPUTE_LAI
type(ChunkIO_t) :: io_checksum3
#endif

    call init_ent_labels

RESOUT = '1kmx1km'

! ----------------------------------------------------------------------
!     GET FILES AND VARS IDs

!**   INPUT Files at 1km x 1km 
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 100, 120)

!     LAI
call chunker%nc_open_gz(io_laiin, &
    DATA_DIR, DATA_INPUT, &
    'LAI/', &
    'LAI3gMax_1kmx1km.nc', 'laimax', 1)

!     Get variable IDs
err = NF90_INQ_VARID(io_laiin%fileid,'lon',varidx)
err = NF90_INQ_VARID(io_laiin%fileid,'lat',varidy)
write(*,*) err, 'variables IDs'

!     CROPS

call chunker%nc_open_gz(io_04crops, &
    DATA_DIR, DATA_INPUT,  &
    'crops/', &
    '04_Monfreda_herb_crops_1kmx1km.nc', 'crops', 1)

call chunker%nc_open_gz(io_05crops, &
    DATA_DIR, DATA_INPUT,  &
    'crops/', &
    '05_Monfreda_shrub_crops_1kmx1km.nc', 'crops', 1)

call chunker%nc_open_gz(io_06crops, &
    DATA_DIR, DATA_INPUT,  &
    'crops/', &
    '06_Monfreda_tree_crops_1kmx1km.nc', 'crops', 1)

call chunker%nc_open_gz(io_04cropsm, &
    DATA_DIR, DATA_INPUT,  &
    'crops/', &
    '08_Monfreda_c4_crops_multi1_1kmx1km.nc', 'crops', 1)

!     CLIMSTATS

call chunker%nc_open_gz(io_C4norm, &
     DATA_DIR, DATA_INPUT,  &
     'climstats/', &
     'CRU_GPCC_C4norm_1981-2010_1kmx1km.nc', 'C4climate', 1)

call chunker%nc_open_gz(io_Tcold, &
     DATA_DIR, DATA_INPUT, &
     'climstats/', 'Tcold.nc', 'Tcold', 1)

call chunker%nc_open_gz(io_Pdry, &
     DATA_DIR, DATA_INPUT,  &
     'climstats/', 'Pdry.nc', 'Pdry', 1)
     
call chunker%nc_open_gz(io_Pmave, &
     DATA_DIR, DATA_INPUT, &
     'climstats/', 'Pmave.nc', 'Pmave', 1)

call chunker%nc_open_gz(io_TCinave, &
      DATA_DIR, DATA_INPUT, &
      'climstats/', 'TCinave.nc', 'TCinave', 1)

call chunker%nc_open_gz(io_CMedit, &
     DATA_DIR, DATA_INPUT,  &
     'climstats/', 'ClimMedit.nc', 'ClimMedit', 1)

!     WATERLC MODIS PARTITION
call chunker%nc_open_gz(io_waterpart, &
     LC_LAI_FOR_1KM1KM_DIR, LC_LAI_FOR_1KM1KM_INPUT, &
     '2004/', &
     'PART_SUB_1km_2004_geo.PARTITION_00.nc', 'PARTITION_0', 1)


! ===================================================
! Create Output Files

!      ENTPFTLC -- Land Cover
call chunker%nc_create(ioall_lc, &
    weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
    'pure/annual/', 'entmm29_ann_lc', 'lc', &
    'Ent Landcover (from A00)', '1', 'Land Cover', &
    ent20%mvs, ent20%abbrev)
! Open water first because it's used to weight others
call chunker%nc_reuse_var(ioall_lc, io_lc(CV_WATER), &
    (/1,1,CV_WATER/), weighting(chunker%wta1, 1d0,0d0))
do k=1,NENT20
    if (k == CV_WATER) cycle
    call chunker%nc_reuse_var(ioall_lc, io_lc(k), &
        (/1,1,k/), weighting(io_lc(CV_WATER)%buf, -1d0,1d0))   ! Land-weighted LC
end do

!     CHECKSUM
call chunker%nc_create(io_checksum, weighting(chunker%wta1,1d0,0d0), &
    'pure/annual/checksum/', 'modis_ann_lc_checksum', &
    'EntMM29lc_lai_for_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)

!     CHECKSUM
call chunker%nc_create(io_checksum2, weighting(chunker%wta1,1d0,0d0), &
     'pure/annual/checksum/', 'entmm29_ann_lc_checksum', &
     'EntLandcover_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)




! ------------------------------------------------------------
! Low-res version computed specially for these
!     NPFTGRID  Number of PFTs in a gridcell
call chunker%nc_create(io_npftgrid, weighting(chunker%wta1,1d0,0d0), &
    'pure/annual/checksum/', 'entmm29_ann_npftgrid', &
     'EntPFTs_percell_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)
io_npftgrid%regrid_lr => accum_lr_stats

!     DOMPFTLC   Dominant PFT's LC in a gridcell
call chunker%nc_create(io_dompftlc, weighting(io_lc(CV_WATER)%buf,-1d0,1d0), &  ! LC is Land-weighted &
    'pure/annual/checksum/', 'entmm29_ann_dompftlc', &
    'EntdominantPFT_LC_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)
io_dompftlc%regrid_lr => nop_regrid_lr

!     DOMPFT     Dominant PFT index in a gridcell (int)
call chunker%nc_create(io_dompft, weighting(io_dompftlc%buf,1d0,0d0), &
    'pure/annual/checksum/', 'entmm29_ann_dompft', &
     'EntdominantPFT_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)
io_dompft%regrid_lr => nop_regrid_lr
! ------------------------------------------------------------

! MODIS PARTITION FILES
do k = 1,LCLASS
   call chunker%nc_open_gz(partit_io(k), &
        LC_LAI_FOR_1KM1KM_DIR, LC_LAI_FOR_1KM1KM_INPUT, &
        '2004/',  &
        'PART_SUB_1km_2004_geo.PARTITION_'//itoa2(k)//'.nc', &
        'PARTITION_'//trim(itoa(k)), 1)    ! var name
enddo

! ENTPFTLAIMAX
call chunker%nc_create(ioall_laiout, &
    weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
    'pure/annual/', 'entmm29_ann_laimax', 'lai', &
    'Ent maximum LAI for year', 'm^2 m-2', 'Leaf Area Index', &
    ent20%mvs, ent20%layer_names())
#ifdef COMPUTE_LAI
do k=1,NENT20
#else
do k=CV_WATER,CV_WATER    ! Compute water LAI even if not the other LAIs
#endif
    call chunker%nc_reuse_var(ioall_lc, io_laiout(k), &
        (/1,1,k/), weighting(io_lc(k)%buf, 1d0,0d0))
enddo


!! TODO: Fix these if we set COMPUTE_LAI
!! ENTPFTLAIMAX CHECKSUM
!#ifdef COMPUTE_LAI
!do k=1,NENT20
!    if (k == CV_WATER) cycle
!    call chunker%nc_create(io_laicheck(k), &
!        weighting(chunker%wta1,1d0,0d0), &
!        'checksum/', &
!        itoa2(k)//'_'//trim(ent20%abbrev(k))//'_lai', &
!        trim(ent20%abbrev(k)), &
!        ent20%title(k), '1', TITLE_CHECKSUM, 1)
!enddo
!#endif
!
!! CHECKSUM
!#ifdef COMPUTE_LAI
!call chunker%nc_create(io_checksum3, weighting(chunker%wta1,1d0,0d0), &
!    'checksum/', 'EntLAI_check_sum_Jun_1kmx1km', &
!    'EntLAI_check_sum_Jun_1kmx1km', &
!    'checksum', '1', TITLE_CHECKSUM)
!#endif


! Quit if we had any problems opening files
call chunker%nc_check('A00_LAI3g_modis_entpftrevcrop')
#ifdef JUST_DEPENDENCIES
stop 0
#endif



   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'long_name','Land Cover')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL,'history', &
        'June 2017: C. Montes, N.Y. Kiang')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'title', 'Ent PFT 1 km land cover fraction')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'creator_name', 'NASA GISS')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'creator_email', &
        'elizabeth.fischer@columbia.edu,'// &
        'carlo.montes@nasa.gov'// &
        'nancy.y.kiang@nasa.gov')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'geospatial_lat_min', '-90')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'geospatial_lat_max', '90')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'geospatial_lon_min', '-180')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'geospatial_lon_max', '180')
   err = NF90_PUT_ATT(ioall_lc%fileid,NF90_GLOBAL, &
        'EntTBM', 'Ent Terrestrial Biosphere Model')

   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'long_name','Maximum LAI over a Year')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
         'history', 'June 2017: C. Montes, N.Y. Kiang,'// &
         'downscaled from 1/12 degree to 1km resolution')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'institution', 'Original data:  LAI3g,'// &
        'Zhu Z.C. et al. 2013 RemSens 5(2):927-948.,'// &
        'Scaling: NASA Goddard Institute for Space Studies')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
         'title', 'Maximum annual LAI (m2/m2) 2004'// &
         'downscaled from 1/12 degrees')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'creator_name', 'NASA GISS')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'creator_email', &
        'elizabeth.fischer@columbia.edu,'// &
        'carlo.montes@nasa.gov,'// &
        'nancy.y.kiang@nasa.gov')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'geospatial_lat_min', '-90')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'geospatial_lat_max', '90')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'geospatial_lon_min', '-180')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'geospatial_lon_max', '180')
   err = NF90_PUT_ATT(ioall_laiout%fileid,NF90_GLOBAL, &
        'EntTBM', 'Ent Terrestrial Biosphere Model')



!-----------------------------------------------------------------
!     Read lat and lon values

startY(1)=1
countY(1)=JM1km
allocate(lat(countY(1)))
err=nf90_get_var(io_laiin%fileid, varidy, lat, startY, countY)

startX(1)=1
countX(1)=IM1km
allocate(lon(countX(1)))
err=nf90_get_var(io_laiin%fileid, varidx, lon, startX, countX)

!-----------------------------------------------------------------

!-----------------------------------------------------------------

! Use these loop bounds for testing...
! it chooses a land area in Asia

#ifdef ENTGVSD_DEBUG
do jchunk = chunker%nchunk(2)*3/4,chunker%nchunk(2)*3/4+1
do ichunk = chunker%nchunk(1)*3/4,chunker%nchunk(1)*3/4+1
#else
do jchunk = 1,chunker%nchunk(2)
do ichunk = 1,chunker%nchunk(1)
#endif

   call chunker%move_to(ichunk,jchunk)

   do jc = 1,chunker%chunk_size(2)
   do ic = 1,chunker%chunk_size(1)

       ! Compute overall NetCDF index of current cell
       ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
       jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

       !**   LAI data
       LAIMAX=io_laiin%buf(ic,jc)

       !**   Crop files
      CROPSHERBNORM=io_04crops%buf(ic,jc)
      CROPSSHRUBNORM=io_05crops%buf(ic,jc)
      CROPSTREENORM=io_06crops%buf(ic,jc)
      CROPSC4HERBFRAC=io_04cropsm%buf(ic,jc)
      if (CROPSC4HERBFRAC.eq.undef) then
          CROPSC4HERBFRAC = 0
      endif

       !** Input C4 climate file
      C4CLIMFRAC=io_C4norm%buf(ic,jc)

       !* Input climate statistics files
       !###  Now getting MAT from Climstats file.
       !###  file tas is in K, Climstats is in C.
      Tcold=io_Tcold%buf(ic,jc)
      Pdry=io_Pdry%buf(ic,jc)
      Pmave=io_Pmave%buf(ic,jc)

      TCinave=io_TCinave%buf(ic,jc)
      MAT = TCinave + 273.15 !Convert to Kelvin

      ClimMedit=io_CMedit%buf(ic,jc)

!----------------------------------------------------------------------
!**   ASSIGN MODIS COVER CONVERSION TO ENT PFT **!

!     * Loop through LCLASS for LC, convert to Ent pfts.
!     * Convert LAIMAX also to Ent pfts.
!     ENTPFTLC(:,:,:) = 0.  !Moved to subroutine Zero_ENTPFT
!     ENTPFTLAIMAX(:,:,:) = 0. !Move to subroutine Zero_ENTPFT
      call Zero_ENTPFT
!     ENTCOVSUM(:,:,:) = 0.
      LIN = 0.
      NPFTGRID = 0
      DOMPFTLC = 0.
      DOMPFT = 0
      
      CHECKSUM = 0.0

      WATERLC = io_waterpart%buf(ic,jc)
      ENTPFTLC(CV_WATER) = WATERLC

      MODIS29(1) = WATERLC
      CHECKSUM = CHECKSUM + WATERLC

      WATERLAI = LAIMAX*WATERLC
      ENTPFTLAIMAX(CV_WATER) = WATERLAI


      do k = 1,LCLASS

         LIN = 0.

!     Input file for land cover data
!     MODIS PARTITION FILES
         LIN = partit_io(k)%buf(ic,jc)
!              if (LIN.gt.0.) write(*,*) LIN,k,'kkkk'

         MODIS29(k+1) = LIN
         CHECKSUM = CHECKSUM + LIN

!     if (k.eq.0) then       !*MODIS Water*!
!     call Set_pft(0,LIN,LAIMAX(k,:,:))
         if (k.eq.1) then !*MODIS Evergreen needleleaf forest*
            p = EVER_ND_LATE         !*4. evergreen needleleaf late successional*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.2) then !*MODIS Evergreen broadleaf forest*!
            p = EVER_BR_LATE         !*2. Ent  evergreen broadleaf late succ*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.3) then !*MODIS Deciduous needleleaf forest*!
            p = DECID_ND         !*8. Ent  deciduous needleleaf*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.4) then !*MODIS Deciduous broadleaf forest*!
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.5) then !*MODIS Mixed forest-->Evergreen needleleaf forest*!
            p = EVER_ND_LATE         !*4. Ent  evergreen needleleaf late successional*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.6) then !*MODIS Mixed forest-->Deciduous broadleaf forest*!
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.7) then !*MODIS Closed shrublands*!
            call Set_Shrubtype(MAT,Pmave, LIN, LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif (k.eq.8) then !*MODIS Open shrublands*!
            call Set_Shrubtype(MAT,Pmave, LIN, LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif (k.eq.9) then !*MODIS Woody savannas-->Evergreen needleleaf forest*!
             p = EVER_ND_LATE         !*4. Ent  evergreen needleleaf late successional*!
             call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.10) then !*MODIS Woody savannas-->Deciduous broadleaf forest*!
!     This layer has some cover (<10%) in boreal and temperate zones, so necessary
!     to partition into 6-Ent cold decid trees and 7-Ent drought decid trees 
!     p = 6               !*6. Ent  cold deciduous broadleaf late successional*!
!     p = 7               !*7. Ent  drought deciduous broadleaf*!
!     call Set_pft(p,LIN,LAIMAX(k,:,:))
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.11) then !*MODIS Woody savannas-->shrub*!
            !     OLD before 5/24/2013: 
            !     call Set_Shrubtype(MAT,Pmave, LIN,LAIMAX(k,:,:)) !*9,10. Ent cold- and arid-shrub*! 
            !     REVISED 5/24/2013 to correct boreal zone
            call Set_Woodysavannasshrub_miscat( &
                 C4CLIMFRAC,MAT,LIN,LAIMAX,jj)
         elseif (k.eq.12) then !*MODIS Woody savannas-->Grass*!
            call Set_Grasstype(C4CLIMFRAC, MAT,Pdry,ClimMedit &
                 ,LIN,LAIMAX)
         elseif (k.eq.13)then !*MODIS Savannas-->Grass*!
            call Set_Grasstype(C4CLIMFRAC, MAT,Pdry,ClimMedit &
                 ,LIN,LAIMAX)
         elseif (k.eq.14)then !*MODIS Savannas-->Shrub*!
            call Set_Shrubtype(MAT,Pmave, LIN,LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif (k.eq.15)then !*MODIS Grassland*!
            call Set_Grasstype(C4CLIMFRAC, MAT,Pdry,ClimMedit &
                 ,LIN,LAIMAX)
         elseif (k.eq.16)then !*MODIS Permanent wetlands-->Evergreen needleforest*!
            ! Ent  evergreen needleleaf late successional
            call Set_pft(EVER_ND_LATE,LIN,LAIMAX)
         elseif (k.eq.17)then !*MODIS Permanent wetlands-->Shrub*!
            call Set_Shrubtype(MAT,Pmave, LIN,LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif ((k.eq.18).or.(k.eq.25)) then !*MODIS 18. Croplands-->Cereal crop *!
            ! MODIS 25. Cropland/natural vegetation-->Cereal crop *!
            if (k.eq.18) then
               call FIX_MODIS29_NORTHPOLE_BUG(LIN,LAIMAX, lat(latout))
               call FIX_MODIS29_SOUTHPOLE_BUG(LIN,LAIMAX, lat(latout))
            endif
            call Set_pft(CROPS_C3_HERB,LIN*(1-CROPSC4HERBFRAC),LAIMAX) !*15. Ent C3 herb crop
            call Set_pft(CROPS_C4_HERB,LIN*CROPSC4HERBFRAC, LAIMAX) !*16. Ent C4 herb crop
         elseif ((k.eq.19).or.(k.eq.26))then !*MODIS 19. Croplands-->Broadleaf crop=~C4 crops, wide leaves *!
            ! *MODIS 26. Cropland/natural vegetation-->Broadleaf crop *!
            if (k.eq.19) then
               call FIX_MODIS29_SOUTHPOLE_BUG(LIN,LAIMAX, lat(latout))
            endif

            call Set_pft(CROPS_C4_HERB, LIN, LAIMAX) !*16. Ent C4 herb crop
         elseif (k.eq.20)then !*MODIS Urban and built up*!
            ! * 19 - Bare or sparsely vegetated
            call Set_pft(BARE_SPARSE,LIN,LAIMAX)
         elseif (k.eq.21) then !*MODIS Cropland/natural vegetation-->Evergreen needleleaf forest*!
            p=EVER_ND_LATE           !*4. Ent  evergreen needleleaf late successional*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.22)then !*MODIS Cropland/natural vegetation-->Deciduous broadleaf forest*!
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.23) then !*MODIS Cropland/natural vegetation-->Shrub*!
            call Set_Shrubtype(MAT,Pmave, LIN,LAIMAX) !*9,10. Ent cold- and arid-adapted shrub*!
         elseif (k.eq.24)then !*MODIS Cropland/natural vegetation-->Grass*!
            call Set_Grasstype(C4CLIMFRAC, MAT,Pdry,ClimMedit, LIN,LAIMAX)
         elseif (k.eq.27)then !*MODIS Permanent snow/ice*!
            ! *Ent Permanent snow/ice
            call Set_pft(SNOW_ICE,LIN,LAIMAX)
         elseif (k.eq.28)then !*MODIS Barren or sparsely vegetated*!
            call Set_pft(BARE_SPARSE,LIN,LAIMAX)
         endif
      end do   ! k=1,lclass

        ! --------- CHECK FOR COVER SUMS TO 1.0 -------------------*!
        ! After attempting re-scaling with this section, I determined that
        ! the numerical precision is just not good enough to rescale all 
        ! all cells to produce cover sums of 1.0.
        ! The loops below cut the number of non-unitary sums from 25749 (10%)
        ! to 4123 (!%).  Iterating again is bad, as negative cover values start 
        ! occurring.
        !   The precision to 7-8 decimal places propagates through to the final
        ! VEG outputs, so I'll let Igor re-scale to keep things 0-1 when he
        ! does his smearing.

        ! This confirms that the modis_c2_gissfortran output does not sum to 1.
        io_checksum%buf(ic,jc)=CHECKSUM
       
        ! Check that ENTPFTLC cover sums to 1.0, and rescale <>1.0.
        CHECKSUM = 0.0
        do k=1,NENT20
           CHECKSUM = CHECKSUM + ENTPFTLC(k)
        enddo

!        if (CHECKSUM.ne.1.0) then
!            ! Rescale all cover to sum to 1.
!            do k=1,NENT20
!                ENTPFTLC(k) = ENTPFTLC(k)/CHECKSUM
!            enddo
!
!             CHECKSUM = 0.0
!             do k=1,NENT20
!                 CHECKSUM = CHECKSUM + ENTPFTLC(k)
!             enddo
!        endif

        CHECKSUM = 0.0
        do k=1,NENT20
           CHECKSUM = CHECKSUM + ENTPFTLC(k)
        enddo

        ! ---------------------------------------------------------------
        ! Finish by cover fraction weighted sum.
        do k=1,NENT20
            if (k == CV_WATER) cycle

            if (ENTPFTLC(k).le.0.) then
                ENTPFTLAIMAX(k) = 0.
            else
                NPFTGRID = NPFTGRID + 1
            endif
            ! Update dominant pft LC and number.
            if (ENTPFTLC(k).gt.DOMPFTLC) then
                DOMPFTLC = ENTPFTLC(k)
                DOMPFT = k
            endif
        end do   !k=1,ENTPFTNUM
        ! ---------------------------------------------------------------

      ! Checksumminmax file.
      CHECKSUM = 0.0

      
      do k=1,NENT20
         CHECKSUM = CHECKSUM + ENTPFTLC(k)
         ! Output file for land cover
         io_lc(k)%buf(ic,jc) = ENTPFTLC(k)
      end do  ! k=1,ENTPFTNUM

      io_checksum2%buf(ic,jc)=CHECKSUM

      CHECKSUM=0.0
      TITLE = "WATER (LAI)"

#ifdef COMPUTE_LAI
      do k=1,NENT20
#else
      do k=CV_WATER,CV_WATER    ! Compute water LAI even if not the other LAIs
#endif
         CHECKSUM = CHECKSUM + ENTPFTLC(k)*ENTPFTLAIMAX(k)
         io_laiout(k)%buf(ic,jc)=ENTPFTLAIMAX(k)
      end do  ! k=1,ENTPFTNUM
#ifdef COMPUTE_LAI
      io_checksum3%buf(ic,jc)=CHECKSUM
#endif

      io_npftgrid%buf(ic,jc)=NPFTGRID
      io_dompftlc%buf(ic,jc)=DOMPFTLC
      io_dompft%buf(ic,jc)=DOMPFT

#ifdef COMPUTE_LAI
      do k=1,NENT20
         if (k == CV_WATER) cycle
         io_laicheck(k)%buf(ic,jc)=ENTPFTLAIMAX(k)
      end do  ! k=1,NENT20
#endif


   end do   ! ic=1,chunker%chunk_size(1)
   end do   ! jc=1,chunker%chunk_size(2)

    call chunker%write_chunks

end do     ! ichunk=1,chunker%nchunk(1)
end do     ! jchunk=1,chunker%nchunk(2)

   
call chunker%close_chunks


! ----------------------- Functions inside the program
CONTAINS

! Accumulate stats on low-res variables
subroutine accum_lr_stats(this)
    class(ChunkIO_t) :: this
    ! -------- Locals
    integer :: jc, ic    ! Index WITHIN current chunk (LR)
    integer :: k
    real*4 :: lc
    integer :: NPFTGRID, DOMPFT
    real*4 :: DOMPFTLC

    ! Loop through low-res chunk
    do jc = 1,this%chunker%chunk_size_lr(2)
    do ic = 1,this%chunker%chunk_size_lr(1)

        NPFTGRID = 0
        DOMPFT = 0
        DOMPFTLC = 0.0
        do k=1,NENT20
            if (k == CV_WATER) cycle

            lc = io_lc(k)%buf_lr(ic,jc) 
            if ((lc == 0).or.(lc==undef)) cycle

            ! Count number of PFTs in the gridcell
            NPFTGRID = NPFTGRID + 1.

            !     Update dominant pft LC and number.
            if (lc.gt.DOMPFTLC) then
                DOMPFTLC = lc
                DOMPFT = k
            endif
        end do   !k=1,ENTPFTNUM

        io_npftgrid%buf_lr(ic,jc) = NPFTGRID
        io_dompft%buf_lr(ic,jc) = DOMPFT
        io_dompftlc%buf_lr(ic,jc) = DOMPFTLC

    end do
    end do
end subroutine accum_lr_stats

end program modis_ent
