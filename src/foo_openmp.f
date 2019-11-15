      module bar
! I. Aleinov.  Example program to parallelize with OpenMP.
!
! To compile:
!    ifort -openmp foo_openmp.f
!
!  To run with 28 openmp threads:
!     OMP_NUM_THREADS=28 ./a.out
!  Compare that to using only 1 openmp thread. The commands below will run your program with
!     1 threads vs. 28, and show the time to run each.
!     time  .a.out
!     time OMP_NUM_THREADS=28 ./a.out  
!
!  For submitting jobs in the background, make sure that you ask for Haswell node. 
!     I.e. to submit the job use a command like
!
!     sbatch  -A s1001 -n 28 -C hasw -t 1:00:00 your_script


      implicit none
      contains
      
      subroutine do_it(i,j,y)
      !Some processing that can be called inside an OMP thread.
      integer :: i, j
      real*8 :: y(:,:)
      real*8 :: a
      
      integer k
      
      a = 1
      do k=1,10000000
         a = a * cos(i+.1d0)*sin(i+.1d0)**9
      enddo
      
      
      a = a + i * cos(2*3.14)
      y(i,j) = a + i 
      if( mod(i,22)==0 )print *, y(i,j),i,j !Skip every 22 i, just dummy output
      
      end subroutine do_it
      
      end module bar

!---------------------------------------------------------      
      program foo
      use bar
      !use netcdf
      implicit none

      include 'netcdf.inc'
      
      integer i, j, fid, varid
      integer start2d(2), count2d(2)
      integer :: err
      real*8, allocatable :: x(:,:)  !A global shared array. E.g. data file, output array.
      character(len=50) :: path="/discover/nobackup/projects/"
      character(len=128) :: fname
      allocate(x(72,46))

      fname = trim(path)//"giss/prod_input_files/"//
     &     "V72X46.1.cor2_no_crops.ext.nc"

      err=nf_open( fname, 0, fid)
      err=nf_inq_varid(fid,'evegreen',varid)
      print *, "err", err
      
!$OMP parallel do private(i,j,x,start2d,count2d)
      !Processing you want to parallelize on shared array x.
      !print *, "err", err

      do j=1,46
         do i=1,72
            start2d(1) = j
            start2d(2) = i
            count2d(1) = 1
            count2d(2) = 1
            err=nf_get_vara_real(fid,varid,start2d, count2d,x)
            call do_it(i,j,x)
         enddo
      enddo
!$OMP end parallel do

      deallocate(x)
      end

