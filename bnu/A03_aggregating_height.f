!* PROGRAM TO AGGREGATE VEGETATION HEIGHT  HEIGHT from 1kmx1km to 0.25x0.25
!* Author: Carlo Montes: carlo.montes@nasa.gov
     
!* GFORTRAN COMPILATION ON DISCOVER-SP3
      
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A04_aggregating_height.f
!  gfortran -o myExe_height arrayutil.o A04_aggregating_height.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
! ./myExe_height  

      
C**** HNTR4.FOR   Horizontal Interpolation Program Real*4   2008/06/18
C****
      Subroutine HNTR40 (IMA,JMA,OFFIA,DLATA,
     *     IMB,JMB,OFFIB,DLATB, DATMIS)
C****
C**** HNTR40 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR4.
C**** The 5 Real input values are expected to be Real*4.
C****
C**** Input: IMA = number of cells in east-west direction of grid A
C****        JMA = number of cells in north-south direction of grid A
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit Real*8 (A-H,O-Z)
      Parameter (TWOPI=6.283185307179586477d0)
      Real*4 OFFIA,DLATA, OFFIB,DLATB, DATMIS,DATMCB
      Common /HNTRCB/ SINA(0:540100),SINB(0:540100),
     *     FMIN(1080000),FMAX(1080000),GMIN(540100),GMAX(540100),
     *     IMIN(1080000),IMAX(1080000),JMIN(540100),JMAX(540100),
     *     DATMCB, INA,JNA, INB,JNB
C****
      Write (0,*) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,
     &           DLATB, DATMIS  !DEBUG NK
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>1080000 .or. JMA<1 .or. JMA>540100 .or.
     *     IMB<1 .or. IMB>1080000 .or. JMB<1 .or. JMB>540100)  GoTo 400
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB                 !  width of single A grid cell in scaled domain
      DIB = IMA                 !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
        RIB = (IBp1-1+OFFIB)*IMA !  scaled longitude of eastern edge
 110    If (RIA-RIB)  120,130,140
 120    IA  = IA  + 1
        RIA = RIA + DIA
        GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
 130    IMAX(IB) = IA
        FMAX(IB) = 0
        IA  = IA  + 1
        RIA = RIA + DIA
        IMIN(IBp1) = IA
        FMIN(IBp1) = 0
        GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
 140    IMAX(IB) = IA
        FMAX(IB) = (RIA-RIB)/DIA
        IMIN(IBp1) = IA
        FMIN(IBp1) = 1-FMAX(IB)
 150    IB = IBp1
        IMAX(IMB) = IMAX(IMB) + IMA
c      WRITE (0,915) 'IMIN=',IMIN(1:IMB)
c      WRITE (0,915) 'IMAX=',IMAX(1:IMB)
c      WRITE (0,916) 'FMIN=',FMIN(1:IMB)
c      WRITE (0,916) 'FMAX=',FMAX(1:IMB)
C**** 
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
        FJEQA = .5*(1+JMA)
        Do 210 JA=1,JMA-1
          RJA = (JA+.5-FJEQA)*DLATA !  latitude in minutes of northern edge
 210      SINA(JA) = Sin (RJA*TWOPI/(360*60))
          SINA(0)  = -1
          SINA(JMA)=  1
C****
          FJEQB = .5*(1+JMB)
          Do 220 JB=1,JMB-1
             RJB = (JB+.5-FJEQB)*DLATB !  latitude in minutes of northern edge      
 220        SINB(JB) = Sin (RJB*TWOPI/(360*60))
            SINB(0)  = -1
            SINB(JMB)=  1
C**** 
            JMIN(1) = 1
            GMIN(1) = 0
            JA = 1
            Do 350 JB=1,JMB-1
 310          If (SINA(JA)-SINB(JB))  320,330,340
 320          JA = JA + 1
              GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
 330          JMAX(JB) = JA
              GMAX(JB) = 0
              JA = JA + 1
              JMIN(JB+1) = JA
              GMIN(JB+1) = 0
              GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
 340          JMAX(JB) = JA
              GMAX(JB) = SINA(JA) - SINB(JB)
              JMIN(JB+1) = JA
              GMIN(JB+1) = SINB(JB) - SINA(JA-1)
 350        Continue
            JMAX(JMB) = JMA
            GMAX(JMB) = 0
c      WRITE (0,915) 'JMIN=',JMIN(1:JMB)
c      WRITE (0,915) 'JMAX=',JMAX(1:JMB)
c      WRITE (0,916) 'GMIN=',GMIN(1:JMB)
c      WRITE (0,916) 'GMAX=',GMAX(1:JMB)
            Return
C**** 
C**** Invalid parameters or dimensions out of range
C**** 
 400        Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,
     &           DLATB, DATMIS
            Stop 400
C**** 
  915 Format (/ 1X,A5 / (20I6))
  916 Format (/ 1X,A5 / (20F6.2))
  940 Format ('0Arguments received by HNTRP0 in order:'/
     *  2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *  ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *  2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     * ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     * ' when integrated WTA = 0'/
     * '0These arguments are invalid or out of range.')
            End subroutine HNTR40

        Subroutine HNTR4 (WTA,A,B)
C****
C**** HNTR4 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*4.
C**** 
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      Implicit Real*8 (A-H,O-Z)
      Real*4 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:540100),SINB(0:540100),
     *     FMIN(1080000),FMAX(1080000),GMIN(540100),GMAX(540100),
     *     IMIN(1080000),IMAX(1080000),JMIN(540100),JMAX(540100),
     *     DATMIS, IMA,JMA, IMB,JMB
C****
C**** Interpolate the A grid onto the B grid
C****

      Do 20 JB=1,JMB
        JAMIN = JMIN(JB)
        JAMAX = JMAX(JB)
        Do 20 IB=1,IMB
          IJB  = IB + IMB*(JB-1)
          WEIGHT= 0
          VALUE = 0
          IAMIN = IMIN(IB)
          IAMAX = IMAX(IB)
          Do 10 JA=JAMIN,JAMAX
            G = SINA(JA)-SINA(JA-1)
            If (JA==JAMIN)  G = G - GMIN(JB)
            If (JA==JAMAX)  G = G - GMAX(JB)
            Do 10 IAREV=IAMIN,IAMAX
              IA  = 1 + Mod(IAREV-1,IMA)
              IJA = IA + IMA*(JA-1)
              F   = 1
              If (IAREV==IAMIN)  F = F - FMIN(IB)
              If (IAREV==IAMAX)  F = F - FMAX(IB)
              WEIGHT = WEIGHT + F*G*WTA(IJA)
 10         VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
            B(IJB) = DATMIS
            If (WEIGHT.ne.0)  B(IJB) = VALUE/WEIGHT
 20         Continue
            Return
        End subroutine HNTR4

      Subroutine HNTR4P (WTA,A,B)
C**** 
C**** HNTR4P is similar to HNTR4 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*4.
C**** 
      Implicit Real*8 (A-H,O-Z)
      Real*4 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:540100),SINB(0:540100),
     *     FMIN(1080000),FMAX(1080000),GMIN(540100),GMAX(540100),
     *     IMIN(1080000),IMAX(1080000),JMIN(540100),JMAX(540100),
     *     DATMIS, IMA,JMA, IMB,JMB
C****
      Call HNTR4 (WTA,A,B)
C**** 
C**** Replace individual values near the poles by longitudinal mean
C****
      Do 40 JB=1,JMB,JMB-1
        BMEAN  = DATMIS
        WEIGHT = 0
        VALUE  = 0
        Do 10 IB=1,IMB
          IJB  = IB + IMB*(JB-1)
          If (B(IJB) == DATMIS)  GoTo 20
          WEIGHT = WEIGHT + 1
          VALUE  = VALUE  + B(IJB)
 10     Continue
        If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
 20     Do 30 IB=1,IMB
          IJB  = IB + IMB*(JB-1)
 30       B(IJB) = BMEAN
 40     Continue
        Return
        End subroutine HNTR4P
      
!---------------------------------------------------------------------------------------------

      subroutine calc_lon_lat_025x025(IM,JM,lon,lat)
      implicit none
      integer,intent(in) :: IM,JM
      real*4,intent(inout) :: lon(IM),lat(JM)
      integer :: i,j

      do i=1,IM
         lon(i) = -180.25 + 0.25/2 + .25*i
      enddo

      do j=1,JM
!         lat(j) = 90.25 - 0.25/2 - .25*j
         lat(j) = -90.25 + 0.25/2 + .25*j
      enddo

      end subroutine calc_lon_lat_025x025
      
!---------------------------------------------------------------------------------------------
      
      
      program aggregating_height

      implicit none

      include 'netcdf.inc'

!***************************************************
!*      PREFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*3, parameter :: EntPFT_files1(19) =
     &     (/
     &     '01_',
     &     '02_',
     &     '03_',
     &     '04_',
     &     '05_',
     &     '06_',
     &     '07_',
     &     '08_',
     &     '09_',
     &     '10_',
     &     '11_',
     &     '12_',
     &     '13_',
     &     '14_',
     &     '15_',
     &     '16_',
     &     '17_',
     &     '18_',
     &     '19_'
     &     /)

!***************************************************
!*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*14, parameter :: EntPFT_files2(19) =
     &     (/
     &     'ever_br_early ',
     &     'ever_br_late  ',
     &     'ever_nd_early ',
     &     'ever_nd_late  ',
     &     'cold_br_early ',
     &     'cold_br_late  ',
     &     'drought_br    ',
     &     'decid_nd      ',
     &     'cold_shrub    ',
     &     'arid_shrub    ',
     &     'c3_grass_per  ',
     &     'c4_grass      ',
     &     'c3_grass_ann  ',
     &     'c3_grass_arct ',
     &     'crops_c3_herb ',
     &     'crops_c4_herb ',
     &     'crops_woody   ',
     &     'snow_ice      ',
     &     'bare_sparse   '
     &     /)
      
      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km

      integer, parameter :: IM25 = 1440 !long at 0.25 degrees
      integer, parameter :: JM25 = 720 !lat at 0.25 degrees

      integer, parameter :: NUMLAYERS = 40
      integer, parameter :: NUMHEIGHTS = 19

      integer, parameter :: longin = X1km
      integer, parameter :: latin = Y1km
      integer, parameter :: longout = IM25
      integer, parameter :: latout = JM25
      integer :: err

      character*80 :: filehin,filelaiin
      character*80 :: filewlc,filewlai
      character*80 :: filehout
      integer :: fileidhin
      integer :: fileidwlc, fileidwlai
      integer :: varidwlai,varidwlc
      integer :: fileidhout
      integer :: varidhin, varidhout
      integer :: varidlonh, varidlath
      integer :: varidlatm, varidlonm

      integer :: varidlaimaxin(19)
      integer :: fileidlaimaxin(19)
      integer :: varidlcin(19)
      integer :: fileidlcin(19)
      integer :: fileidlai(19)
      integer :: varidlai(19)

      real*4, parameter :: undef = -1e30

      real*4 :: lon(IM25)
      real*4 :: lat(JM25)

      integer :: k,p,i,j,m
      real*4 :: laimaxin(X1km,Y1km)
      real*4 :: hin(X1km,Y1km),hinpre(X1km,Y1km)
      real*4 :: wlaiin(X1km,Y1km)
      
      real*4 :: hlaiprod(X1km,Y1km)
      real*4 :: hlaiprodout(IM25,JM25),laiout(IM25,JM25)
      real*4 :: heightout(IM25,JM25)
      
      character*60 :: PathFilepre, PathFilepost
      integer :: start3d(4),count3d(4),start(2),count(2)
      integer :: start4d(5),count4d(5)
      character*20 :: inqvarin, inqvarout
      real*4 :: h

      real*4 :: WTIN(longin,latin)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS

      
!     GET FILES AND VARS IDs

!     INPUT FILES

!     HEIGHT IN
      PathFilepre= '../../data/height/'
      PathFilepost= 'EntGVSDmosaic17_height_1kmx1km.nc'
      filehin = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filehin,NF_WRITE,fileidhin)
      err = NF_INQ_VARID(fileidhin,'SimardHeights',varidhin)

!     ENTPFTLAIMAX
      do k=1,19
         PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
         filelaiin = trim(PathFilepre)//trim(EntPFT_files1(k))//
     &        trim(EntPFT_files2(k))//'_lai.nc'
         err = NF_OPEN(filelaiin,NF_WRITE,fileidlaimaxin(k))
         write(*,*) err
         inqvarin = trim(EntPFT_files2(k))
	 err = NF_INQ_VARID(fileidlaimaxin(k),inqvarin,
     &        varidlaimaxin(k))
         write(*,*) err
      enddo



! ------------------------------------------------------------------
! OUTPUT FILES
      
      ! HEIGHT
      PathFilepre= '../../data/height/'
      PathFilepost= 'EntGVSDmosaic17_height_1440x720_lai3g.nc'
      filehout = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filehout,NF_WRITE,fileidhout)
      write(*,*) err
      err = NF_INQ_VARID(fileidhout,'SimardHeights',varidhout)
      write(*,*) err
      err = NF_INQ_VARID(fileidhout,'lon',varidlonh)
      write(*,*) err
      err = NF_INQ_VARID(fileidhout,'lat',varidlath)
      write(*,*) err

! ------------------------------------------------------------------
      
      call calc_lon_lat_025x025(IM25,JM25,lon,lat)

      err = NF_PUT_VARA_REAL(fileidhout,varidlonh,1,IM25,lon)
      write(*,*) err, ' lon'
      err = NF_PUT_VARA_REAL(fileidhout,varidlath,1,JM25,lat)
      write(*,*) err, ' lat'

      hlaiprod(:,:) = 0.0
      laiout(:,:) = 0.0
      heightout(:,:) = 0.0
      
      DLATA = (180./latin) * 60.
      OFFIA = 0.0
      DLATB = (180./latout) * 60.
      OFFIB = 0.0
      DATMIS = undef

      
!     MAIN LOOP

      do k = 1,NUMHEIGHTS+1

         if (k.ge.1 .and. k.le.NUMHEIGHTS) then

         write(*,*) 'main loop, layer = ', k,' of',NUMHEIGHTS
         start(1) = 1
         start(2) = 1
         count(1) = X1km
         count(2) = Y1km
         
         start3d(1) = k
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = X1km
         count3d(3) = Y1km
         
         ! read height
         err = NF_GET_VARA_REAL(fileidhin,varidhin,start3d,
     &        count3d,hinpre)
	 do i = 1,longin
	     do j = 1,latin
	        hin(i,j) = hinpre(i,j+Y1km+1-2*j)
	     enddo
	 enddo	
         write(*,*) err, ' height in', shape(hin)

         do i=1,longin
            do j=1,latin
               if (hin(i,j).eq.0) then
                  hin(i,j) = undef
               endif
            end do
         end do

         ! read laimax

         !* Setup grids
         write(*,*) longin,latin,OFFIA,DLATA,
     &        longout,latout,OFFIB,DLATB,DATMIS
         
         write(*,*) 'Calling HNTR40'
         call HNTR40(longin,latin,OFFIA,DLATA,
     *        longout,latout,OFFIB,DLATB, DATMIS)
         write(*,*) 'Finished HNTR40'

         do i=1,longin
            do j=1,latin
               !* Weights, if any
               WTIN(i,j) = 1.
               if (hin(i,j).eq.undef) WTIN(i,j) = 0.
            end do
         end do

         write(*,*) 'Calling HNTR40'
         call HNTR4P(WTIN, hin, heightout)
         write(*,*) 'Finished HNTR40'

         do i=1,longout
            do j=1,latout
               if (heightout(i,j).eq.undef) then
                  heightout(i,j) = 0
               endif
            end do
         end do

         start3d(1) = k
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM25
         count3d(3) = JM25
         err = NF_PUT_VARA_REAL(fileidhout,varidhout,start3d,count3d,
     &        heightout)
         write(*,*) 'heightout ', k, err
         
!         err = NF_CLOSE(fileidlaimaxin(k))
         
      else

         start3d(1) = k-1
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM25
         count3d(3) = JM25
         err = NF_PUT_VARA_REAL(fileidhout,varidhout,start3d,count3d,
     &        heightout)
         write(*,*) 'heightout ', k, err

      endif
      
         
      enddo
      
!      err = NF_CLOSE(fileidhin)
     
      
      end program aggregating_height
      


