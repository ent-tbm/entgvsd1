!  Program to assign 1kmx1km monthly LAI to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include soil_albedo_by_lc.f

! gfortran -o myExe arrayutil.o soil_albedo_by_lc.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

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

      program soil_albedo_by_lc_bis

      implicit none

      include 'netcdf.inc'

!*      ENTPFTS FILES FOR LC
      character*22, parameter :: EntPFT_files(19) =
     &     (/
     &     '01_ever_br_early_lc.nc',
     &     '02_ever_br_late_lc.nc ',
     &     '03_ever_nd_early_lc.nc',
     &     '04_ever_nd_late_lc.nc ',
     &     '05_cold_br_early_lc.nc',
     &     '06_cold_br_late_lc.nc ',
     &     '07_drought_br_lc.nc   ',
     &     '08_decid_nd_lc.nc     ',
     &     '09_cold_shrub_lc.nc   ',
     &     '10_arid_shrub_lc.nc   ',
     &     '11_c3_grass_per_lc.nc ',
     &     '12_c4_grass_lc.nc     ',
     &     '13_c3_grass_ann_lc.nc ',
     &     '14_c3_grass_arct_lc.nc',
     &     'crops_herb_lc.nc      ',
     &     '17_crops_woody_lc.nc  ',
     &     '18_snow_ice_lc.nc     ',
     &     '19_bare_sparse_lc.nc  ',
     &     'water_lc.nc           '
     &     /)

      character*13, parameter :: EntPFT_names(19) =
     &     (/
     &     'ever_br_early',
     &     'ever_br_late ',
     &     'ever_nd_early',
     &     'ever_nd_late ',
     &     'cold_br_early',
     &     'cold_br_late ',
     &     'drought_br   ',
     &     'decid_nd     ',
     &     'cold_shrub   ',
     &     'arid_shrub   ',
     &     'c3_grass_per ',
     &     'c4_grass     ',
     &     'c3_grass_ann ',
     &     'c3_grass_arct',
     &     'crops_herb   ',
     &     'crops_woody  ',
     &     'snow_ice     ',
     &     'bare_sparse  ',
     &     'water        '
     &     /)

      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM = X1km !long at 1 km
      integer, parameter :: JM = Y1km !lat at 1 km
      
      real*4, parameter :: undef = -1.e30

      character*256 :: PathFilepre,PathFilepost

      character*256 :: filein,fileout
      character*20 :: inqvarin

      integer :: ncidinvis,ncidinnir,ncidinsw
      integer :: ncidoutvis,ncidoutnir,ncidoutsw
      integer :: varidvis,varidnir,varidsw

      real*4 :: lcin(IM,JM)
      real*4 :: lon(IM),lat(JM)

      integer :: i, j, k
      integer :: err, status

      integer :: pftfileid(19),pftvarid(19)
      integer :: pftvaridvis(19),pftvaridnir(19),pftvaridsw(19)

      character*256 :: filepft

      integer :: start2D(2),count2D(2)
      integer :: dd(2),varidx,varidy

      real*4 :: albvisin(IM,JM),albnirin(IM,JM),albswin(IM,JM)

      real*4 :: albvisout(IM,JM),albnirout(IM,JM),albswout(IM,JM)
 
!------------------------------------------------------------------------

!     define Lon and Lat
      call calc_lon_lat(IM,JM,lon,lat)
     
!------------------------------------------------------------------------
      ! albedo VIS input
      filein = 'VIS_Alb_soil_yearly.006.1kmx1km.nc'
      status = nf_open(filein,0,ncidinvis)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'mean'
      status = nf_inq_varid(ncidinvis,inqvarin,varidvis)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albvisin(:,:) = undef
      status = nf_get_var_real(ncidinvis,varidvis,albvisin)
      write(*,*) status,'nf_get_var_real ',inqvarin

      ! albedo NIR input
      filein = 'NIR_Alb_soil_yearly.006.1kmx1km.nc'
      status = nf_open(filein,0,ncidinnir)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'mean'
      status = nf_inq_varid(ncidinnir,inqvarin,varidnir)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albnirin(:,:) = undef
      status = nf_get_var_real(ncidinnir,varidnir,albnirin)
      write(*,*) status,'nf_get_var_real ',inqvarin

      ! albedo SW input
      filein = 'SW_Alb_soil_yearly.006.1kmx1km.nc'
      status = nf_open(filein,0,ncidinsw)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'mean'
      status = nf_inq_varid(ncidinsw,inqvarin,varidsw)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albswin(:,:) = undef
      status = nf_get_var_real(ncidinsw,varidsw,albswin)
      write(*,*) status,'nf_get_var_real ',inqvarin
      
!-----------------------------------------------------------------
      ! albedo VIS output
      fileout='VIS_Alb_soil_yearly.byLandCover.1kmx1km.nc'
      err = nf_open(fileout,NF_WRITE,ncidoutvis)
      write(*,*) err, fileout, ncidoutvis
      err = NF_INQ_VARID(ncidoutvis,'lon',varidx)
      write(*,*) err
      err = NF_INQ_VARID(ncidoutvis,'lat',varidy)
      write(*,*) err
      do k = 1,19
         err=NF_INQ_VARID(ncidoutvis,trim(EntPFT_names(k)),
     &        pftvaridvis(k))
         write(*,*) err,'varid',pftvaridvis(k),EntPFT_names(k)
      enddo

      ! albedo NIR output
      fileout='NIR_Alb_soil_yearly.byLandCover.1kmx1km.nc'
      err = nf_open(fileout,NF_WRITE,ncidoutnir)
      write(*,*) err, fileout, ncidoutnir
      err = NF_INQ_VARID(ncidoutnir,'lon',varidx)
      write(*,*) err
      err = NF_INQ_VARID(ncidoutnir,'lat',varidy)
      write(*,*) err
      do k = 1,19
         err=NF_INQ_VARID(ncidoutnir,trim(EntPFT_names(k)),
     &        pftvaridnir(k))
         write(*,*) err,'varid',pftvaridnir(k),EntPFT_names(k)
      enddo

      ! albedo SW output
      fileout='SW_Alb_soil_yearly.byLandCover.1kmx1km.nc'
      err = nf_open(fileout,NF_WRITE,ncidoutsw)
      write(*,*) err, fileout, ncidoutsw
      err = NF_INQ_VARID(ncidoutsw,'lon',varidx)
      write(*,*) err
      err = NF_INQ_VARID(ncidoutsw,'lat',varidy)
      write(*,*) err
      do k = 1,19
         err=NF_INQ_VARID(ncidoutsw,trim(EntPFT_names(k)),
     &        pftvaridsw(k))
         write(*,*) err,'varid',pftvaridsw(k),EntPFT_names(k)
      enddo

!-----------------------------------------------------------------

      err=NF_PUT_VARA_REAL(ncidoutvis,varidx,1,IM,lon)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidoutvis,varidy,1,JM,lat)
      write(*,*) err

      err=NF_PUT_VARA_REAL(ncidoutnir,varidx,1,IM,lon)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidoutnir,varidy,1,JM,lat)
      write(*,*) err

      err=NF_PUT_VARA_REAL(ncidoutsw,varidx,1,IM,lon)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidoutsw,varidy,1,JM,lat)
      write(*,*) err

!-----------------------------------------------------------------

      do k = 1,19
         
         PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         PathFilepost= EntPFT_files(k)
         filepft = trim(PathFilepre)//trim(PathFilepost)
	 write(*,*) filepft
         err = NF_OPEN(filepft,0,pftfileid(k))
	 write(*,*) err
	 if (k.eq.19) then
	    inqvarin = trim(EntPFT_names(k))//'_lc'
	 else
	    inqvarin = trim(EntPFT_names(k))
	 endif
         err = NF_INQ_VARID(pftfileid(k),inqvarin,pftvarid(k))
	 write(*,*) err
         lcin(:,:) = undef
         status = nf_get_var_real(pftfileid(k),pftvarid(k),lcin)
         write(*,*) status,'nf_get_var_real ',inqvarin

         do i = 1,IM
            do j = 1,JM
               albvisout(i,j) = albvisin(i,j)*lcin(i,j)
               albnirout(i,j) = albnirin(i,j)*lcin(i,j)
               albswout(i,j) = albswin(i,j)*lcin(i,j)
            enddo
         enddo

         start2D(1)=1
         start2D(2)=1
         count2D(1)=IM
         count2D(2)=JM

         ! VIS
         err = NF_PUT_VARA_REAL(ncidoutvis,pftvaridvis(k),
     &      start2D,count2D,albvisout)
         write(*,*) err, 'put', pftvaridvis(k)

         ! NIR
         err = NF_PUT_VARA_REAL(ncidoutnir,pftvaridnir(k),
     &      start2D,count2D,albnirout)
         write(*,*) err, 'put', pftvaridnir(k)

         ! SW
         err = NF_PUT_VARA_REAL(ncidoutsw,pftvaridsw(k),
     &      start2D,count2D,albswout)
         write(*,*) err, 'put', pftvaridsw(k)

         err = NF_CLOSE(pftfileid(k))
         write(*,*) err

      enddo

      err = NF_CLOSE(ncidoutvis)
      write(*,*) err
      err = NF_CLOSE(ncidoutnir)
      write(*,*) err
      err = NF_CLOSE(ncidoutsw)
      write(*,*) err
      
      end program soil_albedo_by_lc_bis
      

