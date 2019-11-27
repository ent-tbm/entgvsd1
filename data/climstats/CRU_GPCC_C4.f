C**** From sheffield.f
C**** Nancy Kiang, August 2009 | Modified by Carlo Montes, July 2015
C**** Read in CRU and GPCC 1951-1980 and 1981-2010 HxH monthly climatology and calculate
C**** several C4 climatologies, e.g.:
C****   Lloyd & Farquhar (1994): 
C****     Tave of warmest month>22 C and linear proportion over 5.94-17.80 C.
C****   Collatz et al. (1998):
C****     Precip in warmest month > 25 mm.


!** 1 km x 1 km version **!

! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include CRU_GPCC_C4.f
! gfortran -o myExe arrayutil.o CRU_GPCC_C4.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
! ./myExe

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

      subroutine calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      implicit none
      integer,intent(in) :: IM1km,JM1km
      real*4,intent(inout) :: lon(IM1km),lat(JM1km)
!      real*8, parameter :: delta= 0.08333333333333  ! 10 km
      real*8, parameter :: delta= 0.008333333333333 ! 1 km
      integer :: i,j
      
      do i=1,IM1km
         lon(i) = -180-delta/2 + delta*i
      end do

      do j=1,JM1km
         lat(j) = -90-delta/2 + delta*j
!         lat(j) = 90+delta/2 - delta*j
      end do
      
      end subroutine calc_lon_lat_1kmx1km

!------------------------------------------------------------------------

      program CRU_GPCC_C4

      implicit none
      include 'netcdf.inc'
      integer, parameter :: IM1km = 43200 !long at 1 km
      integer, parameter :: JM1km = 21600 !lat at 1 km
      integer, parameter :: IMH = 720 !long at 0.5 degrees
      integer, parameter :: JMH = 360 !lat at 0.5 degrees
      integer, parameter :: IM1 = 360 !long at 1 degrees
      integer, parameter :: JM1 = 180 !lat at 1 degrees
      integer, parameter :: IM2 = 144 !long at 2.5 degrees
      integer, parameter :: JM2 = 90 !lat at 2 degrees
      integer, parameter :: IM4X5 = 72 !long at 5 degrees
      integer, parameter :: JM4X5 = 46 !lat at 4 degrees
      
      integer, parameter :: longin = IM1km
      integer, parameter :: latin = JM1km
      integer, parameter :: longout = IM1km
      integer, parameter :: latout = JM1km

      character*3, parameter :: MONTH(12) =
     &     (/
     &     "Jan","Feb","Mar","Apr","May","Jun",
     &     "Jul","Aug","Sep","Oct","Nov","Dec"
     &     /)
      
      character*80 :: TITLE
      character*256 :: filein, fileout, fileout2,filetemp,fileprecip
      character*256 :: fileout_pre, fileout_sfx
      character*256 :: fileoutnc, fileoutnc2
      character*256 :: filetempave, fileprecave
      real*4 :: TCin(longin,latin) !Monthly mean temperature
      real*4 :: TCinave(longin,latin) !Annual mean temperature
      real*4 :: Pmave(longin,latin) !Annual mean precipitation
      real*4 :: Pmmin(longin,latin) !Monthly mean precip
      integer :: TCmaxmonth(longout,latout) !Warmest month (1-12 Jan-Dec)
      real*4 :: TCmaxmonthr4(longout,latout) !Warmest month (1-12 Jan-Dec)
      real*4 :: TCmax(longout,latout) !Mean temperature of warmest month
      real*4 :: Pmmhot(longout,latout) !Mean precip of warmest month
      real*4 :: Pmmmax(longout,latout) !Mean precip of wettest month
      real*4 :: C4FRAC(longout,latout) !Fraction of grasses that are C4
      real*4 :: TC12(longin,latin)  !1=all months have mean T>22 C, linear to 0.
      real*4 :: Pmm12(longin,latin) !1=all months have mean P>25 mm
      real*4 :: Pdry(longin,latin) !#months with P<5 mm
      real*4 :: TC22(longout,latout) !#months where monthly mean TC>22 C.
      real*4 :: TC22P25(longout,latout) !#months where TC>22 C AND Pmm>25 mm
      real*4 :: lon(longout),lat(latout)
      real*4 :: varout(longout,latout)
      real*4 :: varout2(13,longout,latout)
      integer :: i, j,  m, c4scheme, ncidin, ncidout, varid, status
      character*20 :: inqvarin, inqvarout
      character*50 :: PATHCRU, PATHGPCC, PATHtemp, PATHprec
      real*4 :: Tfrac
      real*4 :: Tcold(longout,latout) !Temperature of coldest month
      real*4 :: TC10C(longout,latout) !#months with temperature > 10C
      real*4 :: ClimMedit(longout,latout) !Mediterranean climate, 0 or 1

      TCmaxmonth(:,:) = -1.e30
      TCmax(:,:) = -1.e30
      Pmmin(:,:) = -1.e30
      Pmmhot(:,:) = -1.e30
      Pmmmax(:,:) = -1.e30
      TC12(:,:) = 12.   !Init counter
      Pmm12(:,:) = 12.  !Init counter
      Pdry(:,:) = 0.  !Init counter
      TC22(:,:) = 0.    !Init
      TC22P25(:,:) = 0.    !Init
      Tcold(:,:) = 100.
      TC10C(:,:) = 0.

      !* Get TCmax and TCmaxmonth *!
      !* Input file.
!      filetemp = './CRU_TS3.22/cru_ts3.22_TS_means_1951-1980_HXH.nc'
!      write(*,*) filetemp
!      open(10,file=filetemp,form='unformatted',status='old')

      do m=1,12

        write(*,*) 'Month ', m
C     temperature data
	fileout_pre = './CRU_GPCC/CRU_temp_1kmx1km_'
	fileout_sfx = '.nc'
	filetemp = 
     &        trim(fileout_pre)// MONTH(m) //
     &        trim(fileout_sfx)
        status = nf_open(trim(filetemp), NF_WRITE, ncidin)
        write(*,*) status, 'nf_open in ',trim(filetemp)

        inqvarin = 'tmp'
        status = nf_inq_varid(ncidin, trim(inqvarin), varid)
        write(*,*) status,'nf_inq_varid ',trim(inqvarin)

        status = nf_get_var_real(ncidin, varid, TCin)
        write(*,*) status,'nf_get_var_real ',trim(inqvarin)
        write(*,*) shape(TCin)
        status = nf_close(ncidin)

C     precipitation data
	fileout_pre = './CRU_GPCC/GPCC_prec_1kmx1km_'
	fileout_sfx = '.nc'
	fileprecip = 
     &        trim(fileout_pre)// MONTH(m) //
     &        trim(fileout_sfx)
        status = nf_open(trim(fileprecip), NF_WRITE, ncidin)
        write(*,*) status, 'nf_open in ',trim(fileprecip)

        inqvarin = 'prec'
        status = nf_inq_varid(ncidin, trim(inqvarin), varid)
        write(*,*) status,'nf_inq_varid ',trim(inqvarin)

        status = nf_get_var_real(ncidin, varid, Pmmin)
        write(*,*) status,'nf_get_var_real ',trim(inqvarin)

        write(*,*) shape(Pmmin)
        status = nf_close(ncidin)

!         TCin(m,:,:) = -9999.
!         read(10) TITLE, TCin
!         write(*,*) m, TITLE
!         read(40) TITLE, Pmmin(m,:,:)

         do i=1,longin
            do j=1,latin
               !write(*,*) i,j, TCin(i,j)
               if ((TCin(i,j))>TCmax(i,j)) then
                  TCmax(i,j) = TCin(i,j)
                  TCmaxmonth(i,j) = m
               endif
               if (((TCin(i,j)).gt.22.0).and.
     &              (Pmmin(i,j).gt.25.))
     &                 TC22P25(i,j) = TC22P25(i,j) + 1.0
               if ((TCin(i,j)).gt.22.0)
     &                 TC22(i,j) = TC22(i,j) + 1.0
               if ((TCin(i,j)).lt.22.0) TC12(i,j)=TC12(i,j)-1.0
               if (TCin(i,j).lt.Tcold(i,j)) Tcold(i,j) = TCin(i,j)
               if (TCin(i,j).gt.10.) TC10C(i,j) = TC10C(i,j) + 1.
            enddo
         enddo
      enddo

!      read(40) TITLE, Pmmin(13,:,:) !Annual mean
!      write(*,*) 13, TITLE
!      Pmmin(13,:,:) = Pmmin(13,:,:)*86400.*365*.1

      TC12(:,:) = TC12(:,:)/12.
!      TC22(:,:) = TC22(:,:)
!      do i=1,longin
!         do j=1,latin
!            if (TC12(i,j).ge.10.) then
!               TC12(i,j) = 1.0
!            elseif ((TC12(i,j).gt.1.0).and.(TC12(i,j).lt.10.)) then
!               TC12(i,j) = (1./9.)*TC12(i,j)  - 1./9.
!            else !TC12(i,j).le.1.0
!               TC12(i,j) = 0.
!            endif
!         enddo
!      enddo

!     Get annual mean temperature
      filetempave = './CRU_GPCC/CRU_temp_1kmx1km_Ave.nc'
      status = nf_open(trim(filetempave), NF_WRITE, ncidin)
      write(*,*) status, 'nf_open in ',trim(filetempave)

      inqvarin = 'tmp'
      status = nf_inq_varid(ncidin, trim(inqvarin), varid)
      write(*,*) status,'nf_inq_varid ',trim(inqvarin)

      status = nf_get_var_real(ncidin, varid, TCinave)
      write(*,*) status,'nf_get_var_real ',trim(inqvarin)

!     Get annual mean precipitation
      fileprecave = './CRU_GPCC/GPCC_prec_1kmx1km_Ave.nc'
      status = nf_open(trim(fileprecave), NF_WRITE, ncidin)
      write(*,*) status, 'nf_open in ',trim(fileprecave)

      inqvarin = 'prec'
      status = nf_inq_varid(ncidin, trim(inqvarin), varid)
      write(*,*) status,'nf_inq_varid',trim(inqvarin)

      status = nf_get_var_real(ncidin, varid, Pmave)
      write(*,*) status,'nf_get_var_real ',trim(inqvarin)


!     read(10) TITLE, TCin !Get annual mean temperature
!     close(10)
!      close(40)

      !* Get Pmmhot and Pmmmax *!
      !* Input file.
!      filein = '../../Sheffield/Sheff1951-2006_precip.bi'
!      write(*,*) filein
!      open(40,file=filein,form='unformatted',status='old')
!      fileprecip = '../../Sheffield/Sheff1951-2006_pmm.bi'
      !open(30,file=fileprecip,form='unformatted')

      do m=1,12

C     temperature data
	fileout_pre = './CRU_GPCC/CRU_temp_1kmx1km_'
	fileout_sfx = '.nc'
	filetemp = 
     &        trim(fileout_pre)// MONTH(m) //
     &        trim(fileout_sfx)
        status = nf_open(trim(filetemp), NF_WRITE, ncidin)
        write(*,*) status, 'nf_open in ',trim(filetemp)

        inqvarin = 'tmp'
        status = nf_inq_varid(ncidin, trim(inqvarin), varid)
        write(*,*) status,'nf_inq_varid ',trim(inqvarin)

        status = nf_get_var_real(ncidin, varid, TCin)
        write(*,*) status,'nf_get_var_real',trim(inqvarin)
        write(*,*) shape(TCin)
        status = nf_close(ncidin)

C     precipitation data
	fileout_pre = './CRU_GPCC/GPCC_prec_1kmx1km_'
	fileout_sfx = '.nc'
	fileprecip = 
     &        trim(fileout_pre)// MONTH(m) //
     &        trim(fileout_sfx)
        status = nf_open(trim(fileprecip), NF_WRITE, ncidin)
        write(*,*) status, 'nf_open in ',trim(fileprecip)

        inqvarin = 'prec'
        status = nf_inq_varid(ncidin, trim(inqvarin), varid)
        write(*,*) status,'nf_inq_varid ',trim(inqvarin)

        status = nf_get_var_real(ncidin, varid, Pmmin)
        write(*,*) status,'nf_get_var_real ',trim(inqvarin)

        write(*,*) shape(Pmmin)
        status = nf_close(ncidin)

         !Pmmin(:,:) = -9999.
         !read(40) TITLE, Pmmin
         !write(*,*) m, TITLE 
         !write(30) TITLE,Pmmin(:,:)*86400.*days(m)  !Fix title to mm/month
         
         do i=1,longin
            do j=1,latin
               !* Get precip in warmest month.
               !* Convert from kg/m2/s to mm/month.
               !* rho.h2o = 998.2071 kg/m3 @ 20 C
               !* rho.h2o = 997.0479 kg/m3 @ 25 C
               !* seconds per month = 60*60*24*days(m) = 86400*days(m)
               if (m.eq.TCmaxmonth(i,j)) 
!     &              Pmmhot(i,j) = Pmmin(i,j) !*86400.*days(m)
!               if ((Pmmin(i,j)*86400.*days(m)).lt.25.0)
!     &              Pmm12(i,j) = Pmm12(i,j) - 1.0
     &              Pmmhot(i,j) = Pmmin(i,j)
               if (Pmmin(i,j).gt.Pmmmax(i,j))
     &              Pmmmax(i,j) = Pmmin(i,j)
               if (Pmmin(i,j).lt.25.0)
     &              Pmm12(i,j) = Pmm12(i,j) - 1.0
               if (Pmmin(i,j).lt.7.5)
     &              Pdry(i,j) = Pdry(i,j) + 1.0
          enddo
        enddo
      enddo
      Pmm12(:,:) = Pmm12(:,:)/12.0
!      close(40)

      !* Calculate C4FRAC *!
      !* Output file.
      fileout = 'CRU_GPCC_C4frac_1981-2010_1kmx1km.bin'
      write(*,*) fileout
      open(20, file=fileout, form='unformatted')

      C4FRAC(:,:) = 0.
      c4scheme = 1
        do i=1,longin
          do j=1,latin

            if (c4scheme.eq.1) then
               !* Lloyd & Farquhar 1994:
               !* mean annual temperature linear 5.94 to 17.80 Celsius
                 if (TCinave(i,j) < 5.94) then
                 C4FRAC(i,j) = 0.0
                 elseif ((TCinave(i,j).ge.5.94).and.
     &                 (TCinave(i,j).le.17.80)) then
                 Tfrac = 0.0848*TCinave(i,j) - 0.504
!                  if (TC22(i,j).gt.0.) C4FRAC(i,j) = Tfrac
            if (((Pmmhot(i,j)>25.0*Tfrac).and.
     &                 (TCmax(i,j)>22.0*Tfrac)).or.
     &                 (TC22P25(i,j)>0.0)) C4FRAC(i,j) = Tfrac
            elseif ((TCinave(i,j).gt.17.8).and.(TC22P25(i,j)>0.0)) then
                C4FRAC(i,j) = 1.0
            endif
            endif

            if (c4scheme.eq.2) then
            !* Collatz 1998: GDD>1000, TCmax>22, Pmmax(warm season)>25 mm
            !* (Note:  Collatz is unclear, in one instance saying
            !*    TCmax>22 and Pmmhot>25 mm in the SAME MONTH, but then
            !*    elsewhere saying TCmax>22 while Pmmhot>25 mm in ANY MONTH.??
            !* (Note:  Bonan 2002 for CLM cites Collatz wrong, saying
            !*   Pmm of the *driest* month is > 25 mm (??)).
!               if ((TCmax(i,j).gt.22.0).and.(Pmmhot(i,j).gt.25.0)) !wrong
               if (TC22P25(i,j).gt.0.)
     &              C4FRAC(i,j) = 1.0
            endif

            if (c4scheme.eq.3) then
               !* NYK scheme: linear 5.94 to 22 Celsius
               if ((TCmax(i,j).ge.5.94).and.(TCmax(i,j).le.22.0)) then
                  Tfrac = 0.0498*TCmax(i,j) - 0.2959
                  if (Pmmhot(i,j).gt.25.0) C4FRAC(i,j) = Tfrac
               elseif ((TCmax(i,j).gt.22.0).and.(Pmmhot(i,j)>25.0)) then
                  C4FRAC(i,j) = 1.0
               endif
            endif

            if (c4scheme.eq.4) then
               !* NYK scheme: linear 17 to 22 Celsius
               if ((TCmax(i,j).ge.17).and.(TCmax(i,j).le.22.0)) then
                  Tfrac = 0.16*TCmax(i,j) - 2.72
                  if (Pmmhot(i,j).gt.25.0) C4FRAC(i,j) = Tfrac
               elseif ((TCmax(i,j).gt.22.0).and.(Pmmhot(i,j)>25.0)) then
                  C4FRAC(i,j) = 1.0
               endif
            endif

            if (c4scheme.eq.5) then
               !* NYK scheme: linear 22 to 27 Celsius
               if ((TCmax(i,j).ge.22).and.(TCmax(i,j).le.27.0)) then
                  Tfrac = 0.16*TCmax(i,j) - 3.52
                  if (Pmmhot(i,j).gt.25.0) C4FRAC(i,j) = Tfrac
               elseif ((TCmax(i,j).gt.22.0).and.(Pmmhot(i,j)>25.0)) then
                  C4FRAC(i,j) = 1.0
               endif
            endif

            if (c4scheme.eq.6) then
               !* NYK scheme: linear 50%@22C to 100%@25 Celsius
               !*             linear 25%@12C to 50%@22 C
               !* Approximate Cabido et al. 2008 elevational gradients
               if ((TCmax(i,j).ge.12.0).and.(TCmax(i,j).lt.22)) then
                  Tfrac = 0.03125*TCmax(i,j) - 0.1875
                  if (Pmmhot(i,j).gt.25.0) C4FRAC(i,j) = Tfrac
               elseif ((TCmax(i,j).ge.22).and.(TCmax(i,j).lt.25.0)) then
                  Tfrac = 0.1667*TCmax(i,j) - 3.1667
                  if (Pmmhot(i,j).gt.25.0) C4FRAC(i,j) = Tfrac
               elseif ((TCmax(i,j).ge.25.0).and.(Pmmhot(i,j)>25.0)) then
                  C4FRAC(i,j) = 1.0
               endif
            endif

            if (c4scheme.eq.7) then
               !* NYK scheme: linear 18 to 25 Celsius
               if ((TCmax(i,j).ge.18).and.(TCmax(i,j).le.25.0)) then
                  Tfrac = 0.1071429*TCmax(i,j) - 1.678571
                  if (Pmmhot(i,j).gt.25.0) C4FRAC(i,j) = Tfrac
               elseif ((TCmax(i,j).gt.25.0).and.(Pmmhot(i,j)>25.0)) then
                  C4FRAC(i,j) = 1.0
               endif
            endif

            if (c4scheme.eq.8) then
            !* Collatz 1998: TCmax>22, Pmmax(warm season)>25 mm
            !* (Note:  Collatz is unclear, in one instance saying
            !*    TCmax>22 and Pmmhot>25 mm in the SAME MONTH, but then
            !*    elsewhere saying TCmax>22 while Pmmhot>25 mm in ANY MONTH.??
            !*    Then for mixed C3/C4, description is very unclear, and 
            !*    they have 100% C4 if ALL months have Tmean>22 AND Pmm>25 mm.
            !* (Note:  Bonan 2002 for CLM cites Collatz wrong, saying
            !*   Pmm of the *driest* month is > 25 mm (??)).
               if ((TCmax(i,j).gt.22.0).and.(Pmmhot(i,j).gt.25.0))
     &              C4FRAC(i,j) = 1.0 * TC12(i,j)! * Pmm12(i,j)
            endif

            enddo
         enddo

      !* Mediterranean climate - based on Köppen Csa classification
      !*     with Trewartha modifications to exclude PNW and Argentina.
      write(*,*) "Calculating Mediterranan climate"
      ClimMedit(:,:) = 0.
      do i=1,longin
         do j=1,latin
            if ((TC10C(i,j).ge.8.).and. !7 Italy, 4 Köppen, 8 Trewartha.
     &           ((Tcold(i,j).lt.18.).and.(Tcold(i,j).gt.-3.)).and. !-3 Köppen, -1 Sheffield
     &           (Pmmhot(i,j).lt.55.).and. !Assume warmest month is driest month. 30 mm Köppen, 55 mm Italy
     &           (Pmmhot(i,j).lt.(0.3333*Pmmmax(i,j))).and.
! !    &           (TCmax(i,j).gt.10.).and.  !20 Italy & S.Africa, 22 Köppen
     &           (Pmave(i,j).lt.90.).and. !Annual mean precip<900 mm, Trewartha (Pmmin is cm)
     &           (TCinave(i,j).le.18.))! !Nancy added, mean ann temp, <18 gets rid of Africa,
     &           then
               ClimMedit(i,j) = 1.
            endif
         enddo
      enddo

      fileout2 = 'CRU_GPCC_climstats_1981-2010_1kmx1km.bin'
      write(*,*) fileout2
      open(50,file=fileout2,form='unformatted')

      TITLE = "TC_max_month (month)"
      TCmaxmonthr4(:,:) = TCmaxmonth(:,:)
      write(50) TITLE, TCmaxmonthr4
      TITLE = "TC_max (Celsius)"
      write(50) TITLE, TCmax
      TITLE = "Pmm_hot (mm/month)"  !Was Pmmmax, check old output files
      write(50) TITLE, Pmmhot
      TITLE = "C4 climate (#months with T>22 C and P>25 mm) 1x1"
      write(50) TITLE, TC22P25
      TITLE = "C4 climate (#months with T>22 C) 1x1"
      write(50) TITLE, TC22
      TITLE = "Temperature of coldest month (C) 1x1" 
      write(50) TITLE, Tcold
      TITLE = "Precip (fraction of months > 25 mm)  1x1"
      write(50) TITLE, Pmm12
      TITLE = "Precip (# of months < 7.5 mm)  1x1"
      write(50) TITLE, Pdry

      do m=1,12

	fileout_pre = './CRU_GPCC/GPCC_prec_1kmx1km_'
	fileout_sfx = '.nc'
	fileprecip = 
     &        trim(fileout_pre)// MONTH(m) //
     &        trim(fileout_sfx)
        status = nf_open(trim(fileprecip), NF_WRITE, ncidin)
        write(*,*) status, 'nf_open in ',trim(fileprecip)

        inqvarin = 'prec'
        status = nf_inq_varid(ncidin, trim(inqvarin), varid)
        write(*,*) status,'nf_inq_varid',trim(inqvarin)

        status = nf_get_var_real(ncidin, varid, Pmmin)
        write(*,*) status,'nf_get_var_real',trim(inqvarin)

        write(*,*) shape(Pmmin)
        status = nf_close(ncidin)

         TITLE = "Precip (mm/month) "//MONTH(m)//" GPCC 1km x 1km"
         write(50) TITLE, Pmmin(:,:)

      enddo 

      TITLE = "Mean Annual Precip (cm/yr)  GPCC 1km x 1km"
      write(50) TITLE, Pmave(:,:)

      TITLE = "Mean Annual Temperature (C)  CRU 1km x 1km"
      write(50) TITLE, TCinave(:,:)
      TITLE = "Precip of wettest month (mm/month)        "
      write(50) TITLE, Pmmmax
      TITLE = "TC>10 C (#months of year)                 "
      write(50) TITLE, TC10C
      TITLE = "Mediterranean climate                     "
      write(50) TITLE, ClimMedit

      if (c4scheme.eq.1) then
         TITLE = 
     &        "C4 climate (fraction of herbs) 
     &         1kmx1km Lloyd&Farquhar 1994"
      elseif(c4scheme.eq.2) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Collatz et al 1998"
      elseif(c4scheme.eq.3) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Kiang 5.94-22 C"
      elseif(c4scheme.eq.4) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Kiang 17-22 C"
      elseif(c4scheme.eq.5) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Kiang 22-27 C"
      elseif(c4scheme.eq.6) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Kiang/Cabido"
      elseif(c4scheme.eq.7) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Kiang 18-25"
      elseif(c4scheme.eq.8) then
         TITLE = "C4 climate (fraction of herbs) 1x1 Collatz mixed"
      endif

      write(20) TITLE, C4FRAC
      write(*,*) TITLE
      close(20)

      close(50)
      !close(30)


      !* Netcdf files
      fileoutnc = 'CRU_GPCC_C4frac_1981-2010_1kmx1km.nc'
      status = nf_open(trim(fileoutnc),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc)
      varout(:,:) = C4FRAC(:,:)
      inqvarout = 'C4FRAC'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)

      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)

      !Make and put lon and lat
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)

      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc)

!----------------------------------------------------------
      !* Climstats

      TCmaxmonthr4(:,:) = TCmaxmonth(:,:)
      varout(:,:) = TCmaxmonthr4(:,:)
      fileoutnc2 = 'TCmaxmonth.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'TCmaxmonth'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = TCmax(:,:)
      fileoutnc2 = 'TCmax.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'TCmax'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = Pmmhot(:,:)
      fileoutnc2 = 'Pmmhot.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'Pmmhot'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = TC22P25(:,:)
      fileoutnc2 = 'TC22P25.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'TC22P25'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = TC22(:,:)
      fileoutnc2 = 'TC22.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'TC22'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = Tcold(:,:)
      fileoutnc2 = 'Tcold.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'Tcold'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)

      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)
      varout(:,:) = Pmm12(:,:)
      fileoutnc2 = 'Pmm12.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'Pmm12'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = Pdry(:,:)
      fileoutnc2 = 'Pdry.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'Pdry'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = TCinave(:,:)
      fileoutnc2 = 'TCinave.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'TCinave'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = Pmave(:,:)
      fileoutnc2 = 'Pmave.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'Pmave'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = Pmmmax(:,:)
      fileoutnc2 = 'Pmmmax.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'Pmmmax'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = ClimMedit(:,:)
      fileoutnc2 = 'ClimMedit.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'ClimMedit'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      varout(:,:) = TC10C(:,:)
      fileoutnc2 = 'TC10C.nc'
      status = nf_open(trim(fileoutnc2),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileoutnc2)
      inqvarout = 'TC10C'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      inqvarout = 'lon'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lon)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      inqvarout = 'lat'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,lat)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileoutnc2)

      end program CRU_GPCC_C4

