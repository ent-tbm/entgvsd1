!  Program to assign 1kmx1km monthly LAI to EntPFTs

! GFORTRAN COMPILATION:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include sparse_lai_2_shrubs.f

! gfortran -o myExe arrayutil.o sparse_lai_2_shrubs.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

! ./myExe

      program sparse_lai_2_shrubs

      implicit none

      include 'netcdf.inc'
      
      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM = X1km !long at 1 km
      integer, parameter :: JM = Y1km !lat at 1 km
      
      real*4, parameter :: undef = -1.e30

      character*256 :: PathFilepre,PathFilepost

      character*256 :: file
      character*20 :: inqvarin

      integer :: fileid
      integer :: varidsparse,varidshrubs

      real*4 :: sparsein(IM,JM),shrubsin(IM,JM)
      real*4 :: sparseout(IM,JM),shrubsout(IM,JM)

      integer :: i, j
      integer :: err, status

      integer :: start2d(2),count2d(2),dd(2)

      start2D(1)=1
      start2D(2)=1
      count2D(1)=IM
      count2D(2)=JM

      PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
      PathFilepost= 'Ent_laimax_1kmx1km.nc'
      file  =  trim(PathFilepre)//trim(PathFilepost)
      status = nf_open(file,nf_write,fileid)
      write(*,*) status, 'nf_open in ', file

      inqvarin = 'lai_bare_sparse'
      err=NF_INQ_VARID(fileid,inqvarin,varidsparse)
      write(*,*) err,'varid',varidsparse
      sparsein(:,:) = undef
      status = NF_GET_VARA_REAL(fileid,varidsparse,start2d,count2d,
     &     sparsein)
      write(*,*) status,'nf_get_var_real ',varidsparse

      inqvarin = 'lai_arid_shrub'
      err=NF_INQ_VARID(fileid,inqvarin,varidshrubs)
      write(*,*) err,'varid',varidshrubs
      shrubsin(:,:) = undef
      status = NF_GET_VARA_REAL(fileid,varidshrubs,start2d,count2d,
     &     shrubsin)
      write(*,*) status,'nf_get_var_real ',varidshrubs

      shrubsout(:,:) = 0.
      do i = 1,IM
         do j = 1,JM
            if (sparsein(i,j).gt.0.5) then
               shrubsout(i,j) = sparsein(i,j) + shrubsin(i,j)
            else
               shrubsout(i,j) = shrubsin(i,j)
            endif
         enddo
      enddo

      sparseout(:,:) = 0.
      do i = 1,IM
         do j = 1,JM
            if (sparsein(i,j).le.0.5) then
               sparseout(i,j) = sparsein(i,j)
            else
               sparseout(i,j) = 0.
            endif
         enddo
      enddo
       
      err = NF_PUT_VARA_REAL(fileid,varidsparse,
     &     start2D,count2D,sparseout)
      write(*,*) err, 'put', varidsparse

      err = NF_PUT_VARA_REAL(fileid,varidshrubs,
     &     start2D,count2D,shrubsout)
      write(*,*) err, 'put', varidsparse
      
      err = NF_CLOSE(fileid)
      write(*,*) err

      
      end program sparse_lai_2_shrubs
      

