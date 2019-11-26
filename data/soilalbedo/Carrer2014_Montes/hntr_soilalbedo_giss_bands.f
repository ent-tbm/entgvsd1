! GFORTRAN COMPILATION ON DISCOVER:

! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include hntr_soilalbedo_giss_bands.f
! gfortran -o myExeTrim arrayutil.o hntr_soilalbedo_giss_bands.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
! ./myExeTrim

     

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

      program hntr_soilalbedo

      implicit none
      
      include 'netcdf.inc'

!*      ENTPFTS FILES FOR LC
      character*30, parameter :: files(6) =
     &     (/
     &     'VIS_soil_albedo_giss_bands.nc ',
     &     'NIR1_soil_albedo_giss_bands.nc',
     &     'NIR2_soil_albedo_giss_bands.nc',
     &     'NIR3_soil_albedo_giss_bands.nc',
     &     'NIR4_soil_albedo_giss_bands.nc',
     &     'NIR5_soil_albedo_giss_bands.nc'
     &     /)

      character*4, parameter :: namesin(6) =
     &     (/
     &     'VIS ',
     &     'NIR1',
     &     'NIR2',
     &     'NIR3',
     &     'NIR4',
     &     'NIR5'
     &     /)

      character*12, parameter :: namesout(6) =
     &     (/
     &     'soilalb_VIS0',
     &     'soilalb_NIR1',
     &     'soilalb_NIR2',
     &     'soilalb_NIR3',
     &     'soilalb_NIR4',
     &     'soilalb_NIR5'
     &     /)

      character*14, parameter :: spbands(6) =
     &     (/
     &     '300 - 770 nm  ',
     &     '770 - 860 nm  ',
     &     '860 - 1250 nm ',
     &     '1250 - 1500 nm',
     &     '1500 - 2200 nm',
     &     '2200 - 4000 nm'
     &     /)

      real*4, parameter :: undef = -1.e30

      integer, parameter :: IMin=7200., JMin=3600.
      integer, parameter :: IMout=43200., JMout=21600.

      integer :: fileid,dimidx,dimidy,dd(2),varidx,varidy
      integer startB(2), countB(2)
      character*256 :: filein,fileout
      character*20 :: inqvarin, inqvarout
      character*80 :: title
      character*37 :: long_title

      integer :: err

      real*4 :: WTIN(IMin,JMin)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS

      real*4 :: albin(IMin,JMin)
      real*4 :: albout(IMout,JMout)

      integer i, j, k
      real*8 lat
      integer count

      integer :: start2D(2),count2D(2)

      integer :: ncidin(6),ncidout,varidin(6),varidout(6),status
      real*4 :: long(IMout),lati(JMout)
 
!------------------------------------------------------------------------

!     define Lon and Lat
      call calc_lon_lat(IMout,JMout,long,lati)

!------------------------------------------------------------------------

      DLATA = (180./JMin) * 60.
      OFFIA=0.0
      DLATB = (180./JMout) * 60.
      OFFIB=0.0
      DATMIS=undef
      
      !* Setup grids.
      write(*,*) IMin,JMin,OFFIA,DLATA,
     &     IMout,JMout,OFFIB,DLATB,DATMIS
       
      write(*,*) 'Calling HNTR40'
      call HNTR40(IMin,JMin,OFFIA,DLATA,
     *     IMout,JMout,OFFIB,DLATB, DATMIS)
      write(*,*) 'Finished HNTR40'

      !------------------------------------------------------------------------

      ! albedo output
      fileout='soilalbedo_giss_bands_1kmx1km.nc'
      err = NF_CREATE(fileout,NF_64BIT_OFFSET,ncidout)
      write(*,*) err, fileout, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IMout,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JMout,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_PUT_ATT_TEXT(ncidout,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(ncidout,varidx,
     &     'units',12,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      err=NF_PUT_ATT_TEXT(ncidout,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(ncidout,varidy,
     &     'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,6
         inqvarout = namesout(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,
     &           varidout(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(ncidout,varidout(k),
     &        'long_name',12,namesout(k))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(ncidout,varidout(k),
     &        'units',8,"fraction")
         write(*,*) err
         err=NF_PUT_ATT_TEXT(ncidout,varidout(k),
     &        'spectral_band',14,spbands(k))
         write(*,*) err
      enddo

      err = NF_PUT_ATT_TEXT(ncidout,NF_GLOBAL,
     &     'xlabel',26,'GISS GCM bands soil albedo')
      write(*,*) err
      err = NF_PUT_ATT_TEXT(ncidout,NF_GLOBAL,
     &     'history',39,'MODIS-derived VIS and NIR soil albedo')
      write(*,*) err
      err = NF_PUT_ATT_TEXT(ncidout,NF_GLOBAL,
     &     'original_data',21,'Carrer et al 2014 RSE')
      write(*,*) err
      err = NF_PUT_ATT_TEXT(ncidout,NF_GLOBAL,
     &     'institution',39,'NASA GISS, C.Montes, N.Kiang, I.Aleinov')
      write(*,*) err
         
      err=NF_ENDDEF(ncidout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IMout,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JMout,lati)
      write(*,*) err
      
      
      do k = 1,6

         filein = trim(files(k))
         status = nf_open(filein,0,ncidin(k))
         write(*,*) status, 'nf_open in ', filein
         inqvarin = namesin(k)
         status = nf_inq_varid(ncidin,inqvarin,varidin(k))
         write(*,*) status,'nf_inq_varid ',inqvarin
         albin(:,:) = undef
         status = nf_get_var_real(ncidin,varidin(k),albin)
         write(*,*) status,'nf_get_var_real ',inqvarin
         albout(:,:) = undef
         do i=1,IMin
            do j=1,JMin
               WTIN(i,j) = 1.
               if (albin(i,j).eq.undef) WTIN(i,j) = 0.
            end do
         end do
         write(*,*) shape(WTIN), shape(albin), shape(albout)
         call HNTR4P(WTIN, albin, albout)
      
         do i = 1,IMout
            do j = 1,JMout
               if (albout(i,j).lt.0) albout(i,j) = 0
            enddo
         enddo

         start2D(1)=1
         start2D(2)=1
         count2D(1)=IMout
         count2D(2)=JMout

         err = NF_PUT_VARA_REAL(ncidout,varidout(k),
     &      start2D,count2D,albout)
         write(*,*) err, 'put', varidout(k)

         err = NF_CLOSE(ncidin(k))
         write(*,*) err

      enddo

      err = NF_CLOSE(ncidout)
      write(*,*) err

      end program hntr_soilalbedo

