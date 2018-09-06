!arrayutil.f
!Author:  Nancy Y. Kiang
!Routines for doing matrix arithmetic.
!*Specialized to exclude UNDEF, NaN.

#define GCM_DRV

      module arrayutil
      
      implicit none
      private
      save

      public IJMULT,IJDIV,IJADD,IJSUB  !real*8 functions
      public IJMULT4,IJADD4,IJSUB4,IJDIV4 !real*4 functions

      !---- FROM OTHER MODULES ---
      !use CONSTANT, only : undef
      real*4, parameter :: undef = -1.e30
  
      contains

  !---- Utility functions -----------------------------------

      logical function NONUM(num)
      real*8 :: num
      if ((num.le.undef).or.(isNaN(num))) then
         NONUM = .true.
      else
         NONUM = .false.
      endif
      end function NONUM

      function NONUM4(num) Result(NONUMres)
      real*4 :: num
      logical :: NONUMres
      if ((num.le.undef).or.(isNaN(num))) then
         NONUMres = .true.
      else
         NONUMres = .false.
      endif
      end function NONUM4


      function CHECKNUM(num) Result(CHECKres)
      !Return 0.0 if num is undef or NaN
      real*8 :: num
      real*8 :: CHECKres
      if (NONUM(num)) then
         CHECKres = 0.d0
      else
         CHECKres = num
      endif
      end function CHECKNUM


      function CHECKNUM4(num) Result(CHECKres)
      !Return 0.0 if num is undef or NaN
      real*4 :: num
      real*4 :: CHECKres
      if (NONUM4(num)) then
         CHECKres = 0.0
      else
         CHECKres = num
      endif
      end function CHECKNUM4

      function DIV0(num, denom) Result(DIVres)
      !Returns 0 if divide by zero or divide by undef.
      real*8 :: num, denom, DIVres
      real*8 :: n,d

      n = CHECKNUM(num)
      d = CHECKNUM(denom)
      if (d.eq.0.0) then
         DIVres = 0.0
      else
         DIVres = n/d
      endif

      end function DIV0

      function DIV04(num, denom) Result(DIVres)
      !Returns 0 if divide by zero or divide by undef.
      real*4 :: num, denom, DIVres
      real*4 :: n,d

      n = CHECKNUM4(num)
      d = CHECKNUM4(denom)
      if (d.eq.0.0) then
         DIVres = 0.0
      else
         DIVres = n/d
      endif

      end function DIV04


      function IJDIV(I0,I1,J0,J1,num, denom) Result(DIVres)
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: num(I0:I1,J0:J1),denom(I0:I1,J0:J1)
      real*8 :: DIVres(I0:I1,J0:J1)
      integer :: i,j

      do i=I0,I1
         do j=J0,J1
            DIVres(i,j) = DIV0(num(i,j),denom(i,j))
         enddo
      enddo
      end function IJDIV

      function IJDIV4(I0,I1,J0,J1,num, denom) Result(DIVres)
      integer,intent(in) :: I0,I1,J0,J1
      real*4 :: num(I0:I1,J0:J1),denom(I0:I1,J0:J1)
      real*4 :: DIVres(I0:I1,J0:J1)
      integer :: i,j

      do i=I0,I1
         do j=J0,J1
            DIVres(i,j) = DIV04(num(i,j),denom(i,j))
         enddo
      enddo
      end function IJDIV4

      function IJMULT(I0,I1,J0,J1,fac1, fac2) Result(MULTres)
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*8 :: MULTres(I0:I1,J0:J1)
      integer :: i,j
      real*8 :: m1, m2

      do i=I0,I1
         do j=J0,J1
            m1 = CHECKNUM(fac1(i,j))
            m2 = CHECKNUM(fac2(i,j))
            MULTres(i,j) = m1 * m2
         enddo
      enddo
      end function IJMULT

      function IJMULT4(I0,I1,J0,J1,fac1, fac2) Result(MULTres)
      integer,intent(in) :: I0,I1,J0,J1
      real*4 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*4 :: MULTres(I0:I1,J0:J1)
      integer :: i,j
      real*4 :: m1, m2

      do i=I0,I1
         do j=J0,J1
            m1 = CHECKNUM4(fac1(i,j))
            m2 = CHECKNUM4(fac2(i,j))
            MULTres(i,j) = m1 * m2
         enddo
      enddo
      end function IJMULT4

      function IJSUB(I0,I1,J0,J1,fac1, fac2) Result(SUBres)
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*8 :: SUBres(I0:I1,J0:J1)
      integer :: i,j
      real*8 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM(fac1(i,j))
            a2 = CHECKNUM(fac2(i,j))
            SUBres(i,j) = a1 - a2
         enddo
      enddo
      end function IJSUB

      function IJSUB4(I0,I1,J0,J1,fac1, fac2) Result(SUBres)
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*4 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*4 :: SUBres(I0:I1,J0:J1)
      integer :: i,j
      real*4 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM4(fac1(i,j))
            a2 = CHECKNUM4(fac2(i,j))
            SUBres(i,j) = a1 - a2
         enddo
      enddo
      end function IJSUB4

      function IJADD(I0,I1,J0,J1,fac1, fac2) Result(ADDres)
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*8 :: ADDres(I0:I1,J0:J1)
      !---------
      integer :: i,j
      real*8 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM(fac1(i,j))
            a2 = CHECKNUM(fac2(i,j))
            ADDres(i,j) = a1 + a2
         enddo
      enddo
      end function IJADD


      function IJADD4(I0,I1,J0,J1,fac1, fac2) Result(ADDres)
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*4 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*4 :: ADDres(I0:I1,J0:J1)
      !---------
      integer :: i,j
      real*8 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM4(fac1(i,j))
            a2 = CHECKNUM4(fac2(i,j))
            ADDres(i,j) = a1 + a2
         enddo
      enddo
      end function IJADD4

!      real*8 DIMENSION(I0:I1,J0:J1) function IJzero(datij)
!      real*8 :: datij(I0:I1,J0:J1)
!      !---------
!      integer :: i,j
!
!      do i=I0,I1
!         do j=J0,J1
!            IJzero(i,j) = CHECKNUM(dat(i,j))
!         enddo
!      enddo
!
!      end function IJzero

      end module arrayutil
