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
!    NOTE: ENTPFTLC(k) == contents of file io_entpft(k) = EntMM_lc_laimax_1kmx1km/...LAI
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
public Debug_set_broadleaf

!     real*4, parameter :: undef = -1e30

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
integer :: latj, JM

if (JM.eq.JM4X5) then     !72x46
   latdeg = -90. + (latj-1)*4.
else                      !anything finer resolution
!     m = (90 - -90)/(1-(JM+1)) = 180./JM
!     b = -90 - 180/JM*1 = -90 - 180/JM
   latdeg = 180./JM*latj - 90. - 180./JM + 180./JM/2.
endif
!      print*, latdeg
end function latdeg
!------------------------------------------------------------------------------------------

subroutine Zero_ENTPFT
ENTPFTLC = 0.  
ENTPFTLAIMAX = 0. 
end subroutine Zero_ENTPFT

subroutine Set_pft(pft,LC_IN,LAI_IN)
implicit none
integer :: pft
real*4 :: LC_IN
real*4 :: LAI_IN
!-----
!     real*4,dimension(:,:),pointer :: LC1 !Mask where LC_IN=0. !BEWARE COMPILER BAD POINTERS
integer :: i,j

!     Quality check write statements
!     do i=1,im
!     do j=1,jm
!     if ((LC_IN(i,j).gt.0.).and.(LAI_IN(i,j).eq.0.)) then
!     !write(*,*) "Non-zero LC for LAI=0",pft, i,j
!     &              ,LC_IN(i,j),LAI_IN(i,j)
!     elseif ((LC_IN(i,j).eq.0.).and.(LAI_IN(i,j).gt.0.)) then
!     !write(*,*) "Zero LC for non-zero LAI",pft, i,j
!     &              ,LC_IN(i,j),LAI_IN(i,j)
!     endif
!     enddo
!     enddo

ENTPFTLC(pft) = ENTPFTLC(pft) + LC_IN
ENTPFTLAIMAX(pft) = ENTPFTLAIMAX(pft)  &
     + LC_IN*LAI_IN

!     ENTCOVSUM(pft,:,:) = ENTCOVSUM(pft,:,:) + LC_IN(:,:)

!     deallocate(LC1)
end subroutine Set_pft
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine Set_Shrubtype(MATEMP,Pmave,LC_IN,LAI_IN)  
!     *9,10. Ent  cold- and arid-adapted shrub*! 
implicit none
real*4 :: MATEMP,Pmave,LC_IN, LAI_IN
!------
integer :: i,j
real*4 :: PFT9,PFT10

PFT9 = 0.
PFT10 = 0.

!     if (MATEMP(i,j).lt.273.15) then !Old 0 C cut-off
!     PFT9(i,j) = LC_IN(i,j)
!     else
!     PFT10(i,j)= LC_IN(i,j)
!     endif
!     if (Pmave(i,j).lt.100.) then      !Alt
!     if ((MATEMP(i,j).lt.278.15)) then ! 5 C cut-off
!     PFT9(i,j) = LC_IN(i,j) !cold
!     else
!     PFT10(i,j)= LC_IN(i,j) !arid
!     endif
!     else if (MATEMP(i,j).lt.278.15) then ! 5 C cut-off
!     PFT9(i,j) = LC_IN(i,j) !cold
!     else
!     PFT10(i,j)= LC_IN(i,j) !arid
!     endif
if (MATEMP.lt.278.15) then !5 C cut-off
   PFT9 = LC_IN
else
   PFT10= LC_IN
endif


call Set_pft(9,PFT9,LAI_IN)
call Set_pft(10,PFT10,LAI_IN)
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
real*4 :: PFT4, PFT8, PFT10

!      allocate(PFT4)     !evergreen needleleaf
!      allocate(PFT8)     !deciduous needleleaf
!      allocate(PFT10)    !arid shrub
PFT4 = 0.
PFT8 = 0.
PFT10 = 0.
!      do j=1,JM1km
   if ((C4CLIMFRAC.le.0.1).and.(latdeg(j,JM1km).gt.45.)) then
      if (MATEMP.ge.(-5.0+273.15)) then !Western Europe boreal
         PFT4 = LC_IN
      else                !Russia/Siberia boreal
         if (LC_IN.ge.0.7) then !Closed forest shades out Larix
            PFT4 = LC_IN 
         elseif ((LC_IN.lt.0.7).and.(LC_IN.ge.0.3)) &
                 then     !Mix Pinus and Larix
            PFT4 = LC_IN - 0.15
            PFT8 = 0.15
         elseif ((LC_IN.lt.0.3).and.(LC_IN.ge.0.15)) &
                 then     !Larix in open/wetlands areas
            PFT8 = LC_IN
         else             !Pinus in sparse areas
            PFT4 = LC_IN
         endif
      endif
   else                   !Rest of the world arid shrub
      PFT10 = LC_IN
   endif
!      enddo
call Set_pft(4,PFT4,LAI_IN)
call Set_pft(8,PFT8,LAI_IN)
call Set_pft(10,PFT10,LAI_IN)
!      deallocate(PFT4)
!      deallocate(PFT8)
!      deallocate(PFT10)

end subroutine Set_Woodysavannasshrub_miscat
!------------------------------------------------------------------------------------------

subroutine Set_Grasstype(C4CLIMFRAC,MATEMP,Pdry,ClimMedit, &
     LC_IN,LAI_IN)
!     * Ent PFTs for grass:
!     * 11- C3 grass perennial
!     * 12 - C4 grass
!     * 13 - C3 grass - annual
!     * 14- arctic C3 grass
implicit none
real*4,intent(in) :: C4CLIMFRAC,MATEMP,Pdry, &
     ClimMedit,LC_IN,LAI_IN
!---------
!     real*4,dimension(:,:),pointer :: C3,PFT11,PFT13,PFT14
real*4 :: C3, PFT11, PFT13, PFT14
integer :: i,j

!     allocate(C3(im,jm))
!     allocate(PFT11(im,jm))
!     allocate(PFT13(im,jm))
!     allocate(PFT14(im,jm))
C3 = 0.
PFT11 = 0.
PFT13 = 0.
PFT14 = 0.

!     * C4 grass
call Set_pft(12,LC_IN*C4CLIMFRAC,LAI_IN)

!     * C3 grass
C3 = 1. - C4CLIMFRAC
if (MATEMP.lt.270) then
   PFT14 = LC_IN*C3
elseif (((Pdry.ge.3.).and.(MATEMP.ge.278.)).or. &
        (ClimMedit.eq.1.)) then
   PFT13 = LC_IN*C3
else
   PFT11 = LC_IN*C3
endif

call Set_pft(11,PFT11,LAI_IN)
call Set_pft(13,PFT13,LAI_IN)
call Set_pft(14,PFT14,LAI_IN)
!     call Set_pft(11, LC_IN*(1.-C4CLIMFRAC),LAI_IN)

!     deallocate(C3)
!     deallocate(PFT11)
!     deallocate(PFT13)
!     deallocate(PFT14)

end subroutine Set_Grasstype

subroutine Set_Broadleaftype(MATEMP,Pdry,C4climfrac,Tcold &
     ,Pmave,ClimMedit,LC_IN,LAI_IN)
implicit none
real*4,intent(in) :: MATEMP,Pdry,C4climfrac &
     ,Tcold,Pmave,ClimMedit &
     ,LC_IN,LAI_IN
!-------
integer :: i,j
!     real*4,dimension(:,:),pointer :: PFT6,PFT7 !BEWARE BAD COMPILER POINTER MEMORY!
real*4 :: PFT6, PFT7
real*4 :: PFT6clim !## DEBUG
character*80 :: TITLE     !## DEBUG

!     allocate(PFT6(im,jm))
!     allocate(PFT7(im,jm))
PFT6 = 0.
PFT7 = 0.
PFT6clim = 0.        !## DEBUG

if ((C4climfrac.lt.0.8).and. &
     (MATEMP.lt.(18.0+273.15)).and. &
     (Tcold.le.8.).and. &
     (ClimMedit.eq.0.)) &
     then
   PFT6 = LC_IN
   PFT6clim = 1.0    !## DEBUG
else
   PFT7 = LC_IN
endif

TITLE = 'PFT6clim in Set_Broad' !##DEBUG
write(40) TITLE, PFT6clim

call Set_pft(6,PFT6,LAI_IN) !cold-deciduous
call Set_pft(7,PFT7,LAI_IN) !drought-deciduous
!     deallocate(PFT6)
!     deallocate(PFT7)

end subroutine Set_Broadleaftype

subroutine Debug_set_broadleaf(modispftclass,LC_IN,LAI_IN)
!     Debug print messages only
implicit none
character*71, intent(in) :: modispftclass
real*4, intent(in) :: LC_IN, LAI_IN
!-----
character*80 :: TITLE

TITLE = 'LC_IN k='//modispftclass
write(40) TITLE,LC_IN
TITLE = 'LAI_IN k='//modispftclass
write(40) TITLE, LAI_IN
TITLE = 'k='//modispftclass(1:3)//' ENTPFTLC 6 colddecidbr'
write(40) TITLE, ENTPFTLC(6)
TITLE = 'k='//modispftclass(1:3)//' ENTPFTLC 7 drought decid'
write(40) TITLE, ENTPFTLC(7)
end subroutine Debug_set_broadleaf

subroutine FIX_MODIS29_NORTHPOLE_BUG(LC_IN,LAI_IN,lat)
!     8/26/13 MODIS29 bug fix:  
!     "18.Croplands-->Cereal crop Cover" has non-zero
!     cover at North Pole where it should be "27. Permanent Ice."
!     Assign to Ent17-18Permanent Ice

implicit none
integer :: pft18ice=18
real*4,intent(inout) :: LC_IN !Input: "18.Croplands-->Cereal crop Cover"
real*4 :: LAI_IN     !Input: "18.Croplands-->Cereal crop Cover"
real*4 :: lat
integer :: northpole

northpole = 79.5

if ((lat.gt.79.5).and.(LC_IN.gt.0.)) then
!         print *,'Wrong crops at north pole:',lat,LC_IN
   ENTPFTLC(pft18ice) =  &
        ENTPFTLC(pft18ice) + LC_IN
   ENTPFTLAIMAX(pft18ice) = ENTPFTLAIMAX(pft18ice)  &
        + 0.              !If LAI_IN not zero, then MODIS is getting algae...
   LC_IN = 0.             !Don't let later calcs usu wrong LC_IN
endif
end subroutine FIX_MODIS29_NORTHPOLE_BUG

subroutine FIX_MODIS29_SOUTHPOLE_BUG(LC_IN,LAI_IN,lat)
!     8/26/13 MODIS29 bug fix:  
!     "18.Croplands-->Cereal crop Cover" has non-zero
!     cover at South Pole corners where it should be "27. Permanent Ice."
!     Assign to Ent17-18Permanent Ice

implicit none
integer :: pft18ice=18
real*4,intent(inout) :: LC_IN !Should be "18.Croplands-->Cereal crop Cover"
real*4 :: LAI_IN    !Should be "18.Croplands-->Cereal crop Cover"
real*4 :: lat
integer :: southpole

southpole = -82.5

if ((lat.lt.-82.5).and.(LC_IN>0.)) then
!         print *,'Wrong crops at south pole:',lat,LC_IN
   ENTPFTLC(pft18ice) =  &
        ENTPFTLC(pft18ice) + LC_IN
   ENTPFTLAIMAX(pft18ice) = ENTPFTLAIMAX(pft18ice)  &
        + 0.              !If it's not zero, then MODIS is getting algae...
   LC_IN = 0.             !Don't let later calcs usu wrong LC_IN
endif
end subroutine FIX_MODIS29_SOUTHPOLE_BUG

end module modis_ent_mod


!------------------------------------------------------------------------
subroutine calc_lon_lat_HxH(IMH,JMH,lon,lat)
implicit none
integer,intent(in) :: IMH,JMH
real*4,intent(inout) :: lon(IMH),lat(JMH)
integer :: i,j

do i=1,IMH
   lon(i) = -180.0 + 0.5*i - 0.5
end do

do j=1,JMH
   lat(j) = -90.0 + 0.5*j - 0.5
end do

end subroutine calc_lon_lat_HxH
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
real*4 :: NPFTGRID  
real*4 :: DOMPFTLC  
real*4 :: DOMPFT

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

type(ChunkIO_t) :: io_lai, io_04crops, io_05crops
type(ChunkIO_t) :: io_06crops, io_04cropsm, io_C4norm
type(ChunkIO_t) :: io_Tcold, io_Pdry, io_Pmave, io_TCinave
type(ChunkIO_t) :: io_CMedit
type(ChunkIO_t) :: io_waterpart,io_checksum
type(ChunkIO_t) :: io_wateroutA,io_checksum2
type(ChunkIO_t) :: io_npftgrid,io_dompftlc
type(ChunkIO_t) :: io_dompft

real*4 :: WATERLAI

type(ChunkIO_t) :: partit_io(28)
type(ChunkIO_t) :: io_entpft(NENT20)
type(ChunkIO_t) :: io_entpftlaimax(NENT20)
type(ChunkIO_t) :: io_entpftlaimaxcheck(NENT20)
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
call chunker%nc_open_gz(io_lai, &
    DATA_DIR, DATA_INPUT, &
    'LAI/', &
    'LAI3gMax_1kmx1km.nc', 'laimax', 1)

!     Get variable IDs
err = NF90_INQ_VARID(io_lai%fileid,'lon',varidx)
err = NF90_INQ_VARID(io_lai%fileid,'lat',varidy)
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

!     CHECKSUM
call chunker%nc_create(io_checksum, weighting(chunker%wta1,1d0,0d0), &
    'checksum/', 'EntMM29lc_lai_for_1kmx1km', &
    'EntMM29lc_lai_for_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)

!     CHECKSUM
call chunker%nc_create(io_checksum2, weighting(chunker%wta1,1d0,0d0), &
     'checksum/', 'EntLandcover_check_sum_Jun_1kmx1km', &
     'EntLandcover_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)


!      ENTPFTLC

! Open water first because we use it later
call chunker%nc_create(io_entpft(CV_WATER), weighting(chunker%wta1,1d0,0d0), &    ! LC land-weighted &
    'EntMM_lc_laimax_1kmx1km/', &
    itoa2(CV_WATER)//'_'//trim(ent20%abbrev(CV_WATER))//'_lc', &
    trim(ent20%abbrev(CV_WATER)), &
    ent20%title(CV_WATER), '1', TITLE_LC)
do k = 1,ent20%ncover
    if (k == CV_WATER) cycle

    call chunker%nc_create(io_entpft(k), weighting(io_entpft(CV_WATER)%buf,-1d0,1d0), &    ! LC land-weighted &
        'EntMM_lc_laimax_1kmx1km/', &
        itoa2(k)//'_'//trim(ent20%abbrev(k))//'_lc', &
        trim(ent20%abbrev(k)), &
        ent20%title(k), '1', TITLE_LC)
enddo


! ------------------------------------------------------------
! Low-res version computed specially for these
!     NPFTGRID
call chunker%nc_create(io_npftgrid, weighting(chunker%wta1,1d0,0d0), &
    'checksum/', 'EntPFTs_percell_check_sum_Jun_1kmx1km', &
     'EntPFTs_percell_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)
io_npftgrid%regrid_lr => accum_lr_stats

!     DOMPFTLC
call chunker%nc_create(io_dompftlc, weighting(io_entpft(CV_WATER)%buf,-1d0,1d0), &  ! LC is Land-weighted &
    'checksum/', 'EntdominantPFT_LC_check_sum_Jun_1kmx1km', &
     'EntdominantPFT_LC_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)
io_dompftlc%regrid_lr => nop_regrid_lr

!     DOMPFT
call chunker%nc_create(io_dompft, weighting(io_dompftlc%buf,1d0,0d0), &
    'checksum/', 'EntdominantPFT_check_sum_Jun_1kmx1km', &
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
#ifdef COMPUTE_LAI
do k=1,NENT20
#else
do k=CV_WATER,CV_WATER    ! Compute water LAI even if not the other LAIs
#endif
    call chunker%nc_create(io_entpftlaimax(k), &
        weighting(io_entpft(k)%buf,1d0,0d0), &   ! Weight by entpft LC
        'EntMM_lc_laimax_1kmx1km/', &
        itoa2(k)//'_'//trim(ent20%abbrev(k))//'_lai', &
        trim(ent20%abbrev(k)), &
        ent20%title(k), 'm2 m-2', TITLE_LAI)
enddo


! ENTPFTLAIMAX CHECKSUM
#ifdef COMPUTE_LAI
do k=1,NENT20
    if (k == CV_WATER) cycle
    call chunker%nc_create(io_entpftlaimaxcheck(k), &
        weighting(chunker%wta1,1d0,0d0), &
        'checksum/', &
        itoa2(k)//'_'//trim(ent20%abbrev(k))//'_lai', &
        trim(ent20%abbrev(k)), &
        ent20%title(k), '1', TITLE_CHECKSUM, 1)
enddo
#endif

! CHECKSUM
#ifdef COMPUTE_LAI
call chunker%nc_create(io_checksum3, weighting(chunker%wta1,1d0,0d0), &
    'checksum/', 'EntLAI_check_sum_Jun_1kmx1km', &
    'EntLAI_check_sum_Jun_1kmx1km', &
    'checksum', '1', TITLE_CHECKSUM)
#endif


! Quit if we had any problems opening files
call chunker%nc_check('A00_LAI3g_modis_entpftrevcrop')
#ifdef JUST_DEPENDENCIES
stop 0
#endif


!-----------------------------------------------------------------
!     Read lat and lon values

startY(1)=1
countY(1)=JM1km
allocate(lat(countY(1)))
err=nf90_get_var(io_lai%fileid, varidy, lat, startY, countY)

startX(1)=1
countX(1)=IM1km
allocate(lon(countX(1)))
err=nf90_get_var(io_lai%fileid, varidx, lon, startX, countX)

!-----------------------------------------------------------------

!-----------------------------------------------------------------

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

       !**   LAI data
       LAIMAX=io_lai%buf(ic,jc)

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
      NPFTGRID = 0.
      DOMPFTLC = 0.
      DOMPFT = 0.
      
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
            p = 4         !*4. evergreen needleleaf late successional*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.2) then !*MODIS Evergreen broadleaf forest*!
            p = 2         !*2. Ent  evergreen broadleaf late succ*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.3) then !*MODIS Deciduous needleleaf forest*!
            p = 8         !*8. Ent  deciduous needleleaf*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.4) then !*MODIS Deciduous broadleaf forest*!
!     p = 6               !*6. Ent  cold deciduous broadleaf late succ*!
!     call Set_pft(6,LIN,LAIMAX(k,:,:))
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.5) then !*MODIS Mixed forest-->Evergreen needleleaf forest*!
            p = 4         !*4. Ent  evergreen needleleaf late successional*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.6) then !*MODIS Mixed forest-->Deciduous broadleaf forest*!
!     p = 6               !*6. Ent  cold deciduous broadleaf late successional*!
!     p = 5               !*5. Ent cold deciduous broadleaf early successional*!
!     call Set_pft(6,LIN,LAIMAX(k,:,:))
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.7) then !*MODIS Closed shrublands*!
            call Set_Shrubtype(MAT,Pmave, LIN, LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif (k.eq.8) then !*MODIS Open shrublands*!
            call Set_Shrubtype(MAT,Pmave, LIN, LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif (k.eq.9) then !*MODIS Woody savannas-->Evergreen needleleaf forest*!
!     p = 3               !*3. Ent  evergreen needleleaf early successional*!
         p = 4         !*4. Ent  evergreen needleleaf late successional*!
!               write(*,*) 'Partitioned MODIS Woody savannas'//
!     &           '-->Evergreen needleleaf forest'
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
!     *4. Ent  evergreen needleleaf late successional*!
            call Set_pft(4,LIN,LAIMAX)
         elseif (k.eq.17)then !*MODIS Permanent wetlands-->Shrub*!
            call Set_Shrubtype(MAT,Pmave, LIN,LAIMAX) !*9,10. Ent  cold- and arid-adapted shrub*!
         elseif ((k.eq.18).or.(k.eq.25)) then !*MODIS 18. Croplands-->Cereal crop *!
!     *MODIS 25. Cropland/natural vegetation-->Cereal crop *!
            if (k.eq.18) then
               call FIX_MODIS29_NORTHPOLE_BUG(LIN,LAIMAX, &
                  lat(latout))
               call FIX_MODIS29_SOUTHPOLE_BUG(LIN,LAIMAX, &
                  lat(latout))
            endif
            call Set_pft(15,LIN*(1-CROPSC4HERBFRAC),LAIMAX) !*15. Ent C3 herb crop
            call Set_pft(16,LIN*CROPSC4HERBFRAC, LAIMAX) !*16. Ent C4 herb crop
         elseif ((k.eq.19).or.(k.eq.26))then !*MODIS 19. Croplands-->Broadleaf crop=~C4 crops, wide leaves *!
!     *MODIS 26. Cropland/natural vegetation-->Broadleaf crop *!
            if (k.eq.19) then
               call FIX_MODIS29_SOUTHPOLE_BUG(LIN,LAIMAX, &
                  lat(latout))
            endif

            call Set_pft(16,LIN, LAIMAX) !*16. Ent C4 herb crop
         elseif (k.eq.20)then !*MODIS Urban and built up*!
!     * 19 - Bare or sparsely vegetated
            call Set_pft(19,LIN,LAIMAX)
         elseif (k.eq.21) then !*MODIS Cropland/natural vegetation-->Evergreen needleleaf forest*!
            p=4           !*4. Ent  evergreen needleleaf late successional*!
            call Set_pft(p,LIN,LAIMAX)
         elseif (k.eq.22)then !*MODIS Cropland/natural vegetation-->Deciduous broadleaf forest*!
!     p = 6               !*6. Ent  cold deciduous broadleaf late successional*!
!     call Set_pft(p,LIN,LAIMAX(k,:,:))
            call Set_Broadleaftype(MAT,Pdry,C4CLIMFRAC,TCOLD &
                 ,Pmave,ClimMedit,LIN,LAIMAX)
         elseif (k.eq.23) then !*MODIS Cropland/natural vegetation-->Shrub*!
            call Set_Shrubtype(MAT,Pmave, LIN,LAIMAX) !*9,10. Ent cold- and arid-adapted shrub*!
         elseif (k.eq.24)then !*MODIS Cropland/natural vegetation-->Grass*!
            call Set_Grasstype(C4CLIMFRAC, MAT,Pdry,ClimMedit &
                 ,LIN,LAIMAX)
         elseif (k.eq.27)then !*MODIS Permanent snow/ice*!
!     *Ent Permanent snow/ice
            call Set_pft(18,LIN,LAIMAX)
         elseif (k.eq.28)then !*MODIS Barren or sparsely vegetated*!
!     * 19 - Bare or sparsely vegetated
            call Set_pft(19,LIN,LAIMAX)
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

        if (CHECKSUM.ne.1.0) then
            ! Rescale all cover to sum to 1.
            do k=1,NENT20
                ENTPFTLC(k) = ENTPFTLC(k)/CHECKSUM
            enddo

             CHECKSUM = 0.0
             do k=1,NENT20
                 CHECKSUM = CHECKSUM + ENTPFTLC(k)
             enddo
        endif

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
                NPFTGRID = NPFTGRID + 1.
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
         io_entpft(k)%buf(ic,jc) = ENTPFTLC(k)
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
         io_entpftlaimax(k)%buf(ic,jc)=ENTPFTLAIMAX(k)
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
         io_entpftlaimaxcheck(k)%buf(ic,jc)=ENTPFTLAIMAX(k)
      end do  ! k=1,NENT20
#endif


   end do   ! ic=1,chunker%chunk_size(1)
   end do   ! jc=1,chunker%chunk_size(2)

    call chunker%write_chunks

end do     ! ichunk=1,nchunk(1)
end do     ! jchunk=1,nchunk(2)

do k=1,NENT20
   
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'long_name',ent20%title(k))
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL,'history', &
        'June 2017: C. Montes, N.Y. Kiang')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'title', 'Ent PFT 1 km land cover fraction')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'creator_name', 'NASA GISS')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'creator_email', &
        'elizabeth.fischer@columbia.edu,'// &
        'carlo.montes@nasa.gov'// &
        'nancy.y.kiang@nasa.gov')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'geospatial_lat_min', '-90')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'geospatial_lat_max', '90')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'geospatial_lon_min', '-180')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'geospatial_lon_max', '180')
   err = NF90_PUT_ATT(io_entpft(k)%fileid,NF90_GLOBAL, &
        'EntTBM', 'Ent Terrestrial Biosphere Model')

#ifndef COMPUTE_LAI
   if (k==CV_WATER) cycle
#endif
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'long_name',ent20%title(k))
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
         'history', 'June 2017: C. Montes, N.Y. Kiang,'// &
         'downscaled from 1/12 degree to 1km resolution')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'institution', 'Original data:  LAI3g,'// &
        'Zhu Z.C. et al. 2013 RemSens 5(2):927-948.,'// &
        'Scaling: NASA Goddard Institute for Space Studies')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
         'title', 'Maximum annual LAI (m2/m2) 2004'// &
         'downscaled from 1/12 degrees')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'creator_name', 'NASA GISS')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'creator_email', &
        'elizabeth.fischer@columbia.edu,'// &
        'carlo.montes@nasa.gov,'// &
        'nancy.y.kiang@nasa.gov')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'geospatial_lat_min', '-90')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'geospatial_lat_max', '90')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'geospatial_lon_min', '-180')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'geospatial_lon_max', '180')
   err = NF90_PUT_ATT(io_entpftlaimax(k)%fileid,NF90_GLOBAL, &
        'EntTBM', 'Ent Terrestrial Biosphere Model')
end do

call chunker%close_chunks


! ----------------------- Functions inside the program
CONTAINS

! Accumulate stats on low-res variables
subroutine accum_lr_stats(this, startB, startB_lr)
    class(ChunkIO_t), intent(INOUT) :: this
    integer, dimension(:), intent(IN) :: startB, startB_lr
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

            lc = io_entpft(k)%buf_lr(ic,jc) 
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
