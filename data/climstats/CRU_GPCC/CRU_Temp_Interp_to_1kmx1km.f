!* CRU_Temp_Interp_to_1kmx1km.f - Specifically set up to interpolate 0.5º CRU temp files

!* COMPILATION MUST BE MADE LOGGED IN "discover-sp3" IN THE FOLLOWING WAY:
!* ulimit -s unlimited
!* module purge
!* module load other/comp/gcc-4.9.2-sp3
!* gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
!* gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include CRU_Temp_Interp_to_1kmx1km.f
!* gfortran -o myExe arrayutil.o CRU_Temp_Interp_to_1kmx1km.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

!* In this case, the executable is called "myExe"

!!!!  CHECK LAND MASK FOR COASTLINES AND ASSIGN WTIN VALUES.

!* To avoid memory problems:
!  > ulimit -s 512000


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

!------------------------------------------------------------------------
      subroutine read_layer(filein, lineskip,
     &     latin,longin,latout,longout,undef_in,undef_out,LAYEROUT)
      implicit none
      character*256 :: filein
      integer :: lineskip
      integer :: latin,longin,latout,longout !#zones in lat and long
      real*4 :: undef_in, undef_out
      real*4 :: LAYEROUT(longout,latout)
      !----Local----
      real*4 :: LAYERIN(longin,latin)
      real*4 :: WTIN(longin,latin)
      integer :: i,j,skip

      write(*,*) filein
      open(10,file=filein)

      do skip=1,lineskip
         read(10,*) 
      end do
      do j=1,latin
         read(10,*) LAYERIN(:,latin-j+1)
      enddo
      do i=1,longin
         do j=1,latin
            if (LAYERIN(i,j).eq.undef_in) LAYERIN(i,j)=undef_out !AViewer undef
            !* Weights, if any
            WTIN(i,j) = 1.      !FGRNDIN
            if (LAYERIN(i,j).le.undef_out) WTIN(i,j) = 0.
         end do
      end do
!      write(*,*) 'Got to before HNTR4P'
      call HNTR4P(WTIN, LAYERIN, LAYEROUT)

      close(10)
      end subroutine read_layer

!------------------------------------------------------------------------
      subroutine calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)
      implicit none
      integer,intent(in) :: IM1km,JM1km
      real*4,intent(inout) :: lon(IM1km),lat(JM1km)
!      real*8, parameter :: delta= 0.08333333333333  ! for 10 km
      real*8, parameter :: delta= 0.008333333333333 ! for 1 km
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

      program CRU_Temp_Interp_to_1kmx1km
      
      implicit none
      include 'netcdf.inc'
      integer, parameter :: IM1km = 43200 !long at 1 km
      integer, parameter :: JM1km = 21600 !lat at 1 km
      integer, parameter :: IM10km = 4320 !long at 10 km
      integer, parameter :: JM10km = 2160 !lat at 10 km
      integer, parameter :: IM12 = 4320 !long at 1/12 degrees
      integer, parameter :: JM12 = 2160 !lat at 1/12 degrees
      integer, parameter :: IMH = 720 !long at 0.5 degrees
      integer, parameter :: JMH = 360 !lat at 0.5 degrees
      integer, parameter :: IM1 = 360 !long at 1 degrees
      integer, parameter :: JM1 = 180 !lat at 1 degrees
      integer, parameter :: IM2 = 144 !long at 2.5 degrees
      integer, parameter :: JM2 = 90 !lat at 2 degrees
      integer, parameter :: IM4X5 = 72 !long at 5 degrees
      integer, parameter :: JM4X5 = 46 !lat at 4 degrees
      character*256 :: file_CRU !30
      integer, parameter :: longin = IMH
      integer, parameter :: latin = JMH
      integer, parameter :: longout = IM1km
      integer, parameter :: latout = JM1km
      integer, parameter :: time = 13

      real*4, parameter :: undef_in = -1e30
      real*4, parameter :: undef_out = -1e30
      
      character*80 :: TITLE !, TITLECHECK
      character*256 :: filein, fileout, filecheck
      character*50 :: PATHEnt, PATHfile
      !real*4 :: varin(longin,latin)
      real*4, ALLOCATABLE :: varin(:,:,:)
      real*4 :: varout(time,longout,latout)
      real*4, ALLOCATABLE :: LAYERIN(:,:)
      !real*4, dimension(longin,latin) :: LAYERIN
      real*4, dimension(longout,latout) :: LAYEROUT
      real*4, dimension(longin,latin) :: WTIN
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS
      real*4 :: tmp(time,longout,latout)
      real*4 :: lon(longout),lat(latout)
      integer :: ncidin,ncidout,varid,status
      integer :: i, j, k, f
      character*20 :: inqvarin, inqvarout
 
      ALLOCATE( varin(time,longin,latin) )
      ALLOCATE( LAYERIN(longin,latin) )

      LAYEROUT(:,:) = undef_out
!     DIVJ = 360. / NINT(360./latout)
!     DLATA = NINT(360./latin) * 60.
      DLATA = (180./latin) * 60.
      OFFIA=0.0
      DLATB = (180./latout) * 60.
      OFFIB=0.0
      DATMIS=undef_out
      
      !* Setup grids.
      write(*,*) longin,latin,OFFIA,DLATA,
     &     longout,latout,OFFIB,DLATB,DATMIS
       
      write(*,*) 'Calling HNTR40'
      call HNTR40(longin,latin,OFFIA,DLATA,
     *     longout,latout,OFFIB,DLATB, DATMIS)
      write(*,*) 'Finished HNTR40'

      !* Input file.
      filein =  'cru_ts3.22_TS_means_1981-2010_HXH_bis.nc'
!      filein =  'cru_ts3.22_TS_means_1981-2010_10kmx10km.nc'
      status = nf_open(trim(filein),NF_WRITE,ncidin)
      write(*,*) status, 'nf_open in ',trim(filein)

      !* Output file.
      fileout = 'cru_ts3.22_TS_means_1981-2010_1kmx1km.nc'
      status = nf_open(trim(fileout),NF_WRITE,ncidout)
      write(*,*) status, 'nf_open out ',trim(fileout)

      inqvarin = 'tmp'
      status = nf_inq_varid(ncidin,trim(inqvarin),varid)
      write(*,*) status,'nf_inq_varid ',trim(inqvarin)

      varin(:,:,:) = undef_in
      status = nf_get_var_real(ncidin,varid,varin)
      write(*,*) status,'nf_get_var_real ',trim(inqvarin)

      write(*,*) 'Looping through layers'
      varout(:,:,:) = undef_out
      
      do k=1,13

        LAYERIN(:,:) = varin(k,:,:)

        LAYEROUT(:,:) = undef_out
        !* Loop through layers
        do i=1,longin
            do j=1,latin
                !write(*,*) 'i,j',i,j
                !AViewer undef
                !* Weights, if any
                WTIN(i,j) = 1.   !FGRNDIN
                if (LAYERIN(i,j).eq.undef_out) WTIN(i,j) = 0.
            end do
        end do


        write(*,*) shape(WTIN), shape(LAYERIN), shape(LAYEROUT)
                
        !write(*,*) i,j, 'Calling HNTR4P...'
        call HNTR4P(WTIN, LAYERIN, LAYEROUT)
        !write(*,*) 'Finished calling HNTR4P'

        varout(k,:,:) = LAYEROUT(:,:)
	
      end do

      call calc_lon_lat_1kmx1km(IM1km,JM1km,lon,lat)

!     Make and put lon and lat
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

      inqvarout = 'time'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,time)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)

      inqvarout = 'tmp'
      status = nf_inq_varid(ncidout,trim(inqvarout),varid)
      write(*,*) status,'nf_inq_varid out ',trim(inqvarout)
      status = nf_put_var_real(ncidout,varid,varout)
      write(*,*) status,'nf_put_var_real out ',trim(inqvarout)
      
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileout)

      DEALLOCATE(varin)
      DEALLOCATE(LAYERIN)

      end program CRU_Temp_Interp_to_1kmx1km

