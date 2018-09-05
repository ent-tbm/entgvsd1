! Combine C3 and C3 crops to create crops_herb array

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include C3C4_crops_2_crops_herb.f

! gfortran -o myExe arrayutil.o C3C4_crops_2_crops_herb.f -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

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

      program C3C4_crops_2_crops_herb

      implicit none
      
      include 'netcdf.inc'

      character*3, parameter :: MONTH(12) =
     &     (/
     &     "Jan","Feb","Mar","Apr","May","Jun",
     &     "Jul","Aug","Sep","Oct","Nov","Dec"
     &     /)
      
      real*4, parameter :: undef = -1.e30

      integer, parameter :: IM=43200., JM=21600.

      integer :: fileid,dimidx,dimidy,dimidz,dd(2),varidx,varidy,varidz
      integer startB(2), countB(2)
      character*256 :: PathFilepre,PathFilepost
      character*256 :: filein,fileinc3,fileinc4,fileout
      character*20 :: inqvarin, inqvarout
      character*80 :: title
      character*37 :: long_title

      integer :: err

      real*4 :: C3in(IM,JM),LCC3(IM,JM)
      real*4 :: C4in(IM,JM),LCC4(IM,JM)
      real*4 :: crops_herb(IM,JM)

      integer i, j, k, m

      integer :: ncidin,ncidout,varid,status
      real*4 :: long(IM),lati(JM), a
      integer :: start3d(3),count3d(3)
      
      !------------------------------------------------------------------------

      ! define Lon and Lat
      call calc_lon_lat(IM,JM,long,lati)

      !------------------------------------------------------------------------
      ! LC
      PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= '15_crops_c3_herb_lc.nc'
      fileinc3 = trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(fileinc3,0,ncidin)
      write(*,*) status, 'nf_open in ', fileinc3
      inqvarin = 'crops_c3_herb'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      LCC3(:,:) = undef
      status = nf_get_var_real(ncidin,varid,LCC3)
      write(*,*) status,'nf_get_var_real ',inqvarin
      err = NF_CLOSE(ncidin)

      PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= '16_crops_c4_herb_lc.nc'
      fileinc4 = trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(fileinc4,0,ncidin)
      write(*,*) status, 'nf_open in ', fileinc4
      inqvarin = 'crops_c4_herb'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      LCC4(:,:) = undef
      status = nf_get_var_real(ncidin,varid,LCC4)
      write(*,*) status,'nf_get_vara_real ',inqvarin
      err = NF_CLOSE(ncidin)

      crops_herb(:,:) = undef
      do i = 1,IM
         do j = 1,JM
            crops_herb(i,j) = LCC3(i,j) + LCC4(i,j)
         enddo
      enddo

      PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= 'crops_herb_lc.nc'
      fileout = trim(PathFilepre)//trim(PathFilepost)

      err = NF_CREATE(fileout,NF_CLOBBER,ncidout)
      write(*,*) err, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      write(*,*) err, 'dimidx',dimidx
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      write(*,*) err, 'dimidy',dimidy
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      write(*,*) err,'varidx',varidx
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      write(*,*) err,'varidy',varidy
      dd(1)=dimidx
      dd(2)=dimidy
      err=NF_DEF_VAR(ncidout,'crops_herb',NF_REAL,2,dd,varid)
      write(*,*) err,'varid',varid
      long_title = 'Herbaceous C3/C4 crops'
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     trim(long_title))
      write(*,*) err
!      err=NF_PUT_ATT_REAL(ncidout,varid,'_FillValue',nf_real,1,undef)
!      write(*,*) err
      err=NF_ENDDEF(ncidout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      write(*,*) err
      startB(1)=1
      startB(2)=1
      countB(1)=IM
      countB(2)=JM
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,crops_herb)
      write(*,*) err, 'put'
      err = NF_CLOSE(ncidout)

      !------------------------------------------------------------------------
      ! LAImax
      PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= '15_crops_c3_herb_lai.nc'
      fileinc3 = trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(fileinc3,0,ncidin)
      write(*,*) status, 'nf_open in ', fileinc3
      inqvarin = 'crops_c3_herb'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      C3in(:,:) = undef
      status = nf_get_var_real(ncidin,varid,C3in)
      write(*,*) status,'nf_get_var_real ',inqvarin
      err = NF_CLOSE(ncidin)

      PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= '16_crops_c4_herb_lai.nc'
      fileinc4 = trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(fileinc4,0,ncidin)
      write(*,*) status, 'nf_open in ', fileinc4
      inqvarin = 'crops_c4_herb'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      C4in(:,:) = undef
      status = nf_get_var_real(ncidin,varid,C4in)
      write(*,*) status,'nf_get_var_real ',inqvarin
      err = NF_CLOSE(ncidin)

      crops_herb(:,:) = undef
      do i = 1,IM
         do j = 1,JM
	    a = LCC3(i,j) + LCC4(i,j)
               if ( a > 0. ) then
                  crops_herb(i,j) = (LCC3(i,j)*C3in(i,j)
     &              + LCC4(i,j)*C4in(i,j)) / a
               else
                  crops_herb(i,j) = 0.
               endif
         enddo
      enddo

      PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= 'crops_herb_laimax.nc'
      fileout = trim(PathFilepre)//trim(PathFilepost)

      err = NF_CREATE(fileout,NF_CLOBBER,ncidout)
      write(*,*) err, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      write(*,*) err, 'dimidx',dimidx
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      write(*,*) err, 'dimidy',dimidy
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      write(*,*) err,'varidx',varidx
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      write(*,*) err,'varidy',varidy
      dd(1)=dimidx
      dd(2)=dimidy
      err=NF_DEF_VAR(ncidout,'crops_herb',NF_REAL,2,dd,varid)
      write(*,*) err,'varid',varid
      long_title = 'Herbaceous C3/C4 crops'
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     trim(long_title))
      write(*,*) err
!      err=NF_PUT_ATT_REAL(ncidout,varid,'_FillValue',nf_real,1,undef)
!      write(*,*) err
      err=NF_ENDDEF(ncidout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      write(*,*) err
      startB(1)=1
      startB(2)=1
      countB(1)=IM
      countB(2)=JM
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,crops_herb)
      write(*,*) err, 'put'
      err = NF_CLOSE(ncidout)

      !------------------------------------------------------------------------
      ! Monthly LAI
      do k = 1,12
         PathFilepre= '../lc_lai_ent/nc/'
         PathFilepost = 'EntMM_lc_lai_'//MONTH(k)//'_1kmx1km.nc'
         filein = trim(PathFilepre)//trim(PathFilepost)
         status = nf_open(filein,0,ncidin)
         write(*,*) status, 'nf_open in ', filein
         inqvarin = 'EntPFT'
         status = nf_inq_varid(ncidin,inqvarin,varid)
         write(*,*) status,'nf_inq_varid ',inqvarin
         start3d(1) = 16
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM
         count3d(3) = JM
         C3in(:,:) = undef
         status = NF_GET_VARA_REAL(ncidin,varid,start3d,count3d,C3in)
         write(*,*) status,'nf_get_var_real C3 ',inqvarin
         start3d(1) = 17
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM
         count3d(3) = JM
         C4in(:,:) = undef
         status = NF_GET_VARA_REAL(ncidin,varid,start3d,count3d,C4in)
         write(*,*) status,'nf_get_var_real C4 ',inqvarin
         err = NF_CLOSE(ncidin)

         crops_herb(:,:) = undef
         do i = 1,IM
            do j = 1,JM
               a = LCC3(i,j) + LCC4(i,j)
                  if ( a > 0. ) then
                     crops_herb(i,j) = (LCC3(i,j)*C3in(i,j)
     &                 + LCC4(i,j)*C4in(i,j)) / a
                  else
                     crops_herb(i,j) = 0.
                  endif
            enddo
         enddo
         
         PathFilepre= '../lc_lai_ent/nc/'
         PathFilepost = 'crops_herb_lc_lai_'//MONTH(k)//'_1kmx1km.nc'
         fileout = trim(PathFilepre)//trim(PathFilepost)

         err = NF_CREATE(fileout,NF_CLOBBER,ncidout)
         write(*,*) err, ncidout
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         write(*,*) err, 'dimidx',dimidx
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         write(*,*) err, 'dimidy',dimidy
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         write(*,*) err,'varidx',varidx
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         write(*,*) err,'varidy',varidy
         dd(1)=dimidx
         dd(2)=dimidy
         err=NF_DEF_VAR(ncidout,'crops_herb',NF_REAL,2,dd,varid)
         write(*,*) err,'varid',varid
         long_title = 'Herbaceous C3/C4 crops'
         err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &        trim(long_title))
         write(*,*) err
!     err=NF_PUT_ATT_REAL(ncidout,varid,'_FillValue',nf_real,1,undef)
!     write(*,*) err
         err=NF_ENDDEF(ncidout)
         write(*,*) err
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         write(*,*) err
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         write(*,*) err
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,crops_herb)
         write(*,*) err, 'put'
         err = NF_CLOSE(ncidout)
      enddo

      end program C3C4_crops_2_crops_herb

