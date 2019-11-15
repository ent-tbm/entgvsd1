!* PROGRAM TO AGGREGATE VEGETATION DATA from 1kmx1km to coarser resolution
!* Author: Carlo Montes: carlo.montes@nasa.gov
     
!* GFORTRAN COMPILATION ON DISCOVER-SP3
      
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include aggregating_laimax.f
!  gfortran -o myExe arrayutil.o aggregating_laimax.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
! ./myExe

      
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
      
!---------------------------------------------------------------------------------------------
      
      
      program aggregating_laimax

      implicit none

      include 'netcdf.inc'

      character*50, parameter :: EntPFT_shorttitle(18) =
     &     (/
     &     "ever_br_early   ",
     &     "ever_br_late    ",
     &     "ever_nd_early   ",
     &     "ever_nd_late    ",
     &     "cold_br_early   ",
     &     "cold_br_late    ",
     &     "drought_br      ",
     &     "decid_nd        ",
     &     "cold_shrub      ",
     &     "arid_shrub      ",
     &     "c3_grass_per    ",
     &     "c4_grass        ",
     &     "c3_grass_ann    ",
     &     "c3_grass_arct   ",
     &     "crops_herb      ",
     &     "crops_woody     ",
     &     "bare_bright     ",
     &     "bare_dark       "
     &     /)
      
       character*40, parameter :: EntPFT_title(18) =
     &     (/
     &     '1 - Evergreen Broadleaf Early Succ      ',
     &     '2 - Evergreen Broadleaf Late Succ       ',
     &     '3 - Evergreen Needleleaf Early Succ     ',
     &     '4 - Evergreen Needleleaf Late Succ      ',
     &     '5 - Cold Deciduous Broadleaf Early Succ ',
     &     '6 - Cold Deciduous Broadleaf Late Succ  ',
     &     '7 - Drought Deciduous Broadleaf         ',
     &     '8 - Deciduous Needleleaf                ',
     &     '9 - Cold Adapted Shrub                  ',
     &     '10 - Arid Adapted Shrub                 ',
     &     '11 - C3 Grass Perennial                 ',
     &     '12 - C4 Grass                           ',
     &     '13 - C3 Grass Annual                    ',
     &     '14 - Arctic C3 Grass                    ',
     &     '15 - Crops Herb                         ',
     &     '16 - Crops Woody                        ',
     &     '17 - Bright Bare Soil                   ',
     &     '18 - Dark Bare Soil                     '
     &     /)

      integer, parameter :: X1km = 43200 !long at 1 km
      integer, parameter :: Y1km = 21600 !lat at 1 km
      integer, parameter :: IM1km = X1km !long at 1 km
      integer, parameter :: JM1km = Y1km !lat at 1 km

!      integer, parameter :: IMout = 1440 !long at 0.25 degrees
!      integer, parameter :: JMout = 720 !lat at 0.25 degrees
      integer, parameter :: IMout = 144 !long at 2.5x2
      integer, parameter :: JMout = 90 !lat at 2.5x2

      integer, parameter :: longin = X1km
      integer, parameter :: latin = Y1km
      integer, parameter :: longout = IMout
      integer, parameter :: latout = JMout
      integer :: err

      character*80 :: filehin,filehout
      integer :: fileidhin
      integer :: fileidhout
      integer :: varidhin(18), varidhout(18)
      integer :: varidx, varidy,dimidx,dimidy

      real*4, parameter :: undef = -1e30

      real*4 :: longi(IMout)
      real*4 :: lati(JMout)

      integer :: k,p,i,j,m
      real*4 :: hin(X1km,Y1km)
      real*4 :: hout(IMout,JMout)
      
      character*60 :: PathFilepre, PathFilepost
      integer :: start2D(2),count2D(2),dd(2)
      character*20 :: inqvarin, inqvarout
      real*4 :: h

      real*4 :: WTIN(longin,latin)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS

      call calc_lon_lat(IMout,JMout,longi,lati)

!     FILE IN

      PathFilepre= '../lc_lai_ent16/nc/'
      PathFilepost=  'V1km_EntGVSDv1.1_LAI3g16_laimax_pure.nc'
      filehin  =  trim(PathFilepre)//trim(PathFilepost)
      err = nf_open(filehin,0,fileidhin)
      write(*,*) err, 'nf_open in ', filehin
      do k=1,18
         inqvarin = 'lai_'//EntPFT_shorttitle(k)
         err = NF_INQ_VARID(fileidhin,inqvarin,varidhin(k))
         write(*,*) err,inqvarin
      enddo
      
!     FILE OUT

      PathFilepre= '../lc_lai_ent16/nc/'
      if (IMout.eq.1440) then
         PathFilepost=  'V1440x720_EntGVSDv1.1_LAI3g16_laimax_pure.nc'
      else
         PathFilepost=  'V144x90_EntGVSDv1.1_LAI3g16_laimax_pure.nc'
      endif
      filehout  =  trim(PathFilepre)//trim(PathFilepost)
      err = NF_CREATE(filehout,NF_FLOAT,fileidhout)
      write(*,*) err, 'nf_create out ',trim(filehout)
      err=NF_DEF_DIM(fileidhout,'lon',IMout,dimidx)
      write(*,*) err
      err=NF_DEF_DIM(fileidhout,'lat',JMout,dimidy)
      write(*,*) err
      err=NF_DEF_VAR(fileidhout,'lon',NF_FLOAT,1,dimidx,
     &     varidx)
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fileidhout,varidx,
     &     'long_name',9,"longitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fileidhout,varidx,
     &        'units',14,"degrees east")
      write(*,*) err
      err=NF_DEF_VAR(fileidhout,'lat',NF_FLOAT,1,dimidy,
     &     varidy)
      err=NF_PUT_ATT_TEXT(fileidhout,varidy,
     &     'long_name',8,"latitude")
      write(*,*) err
      err=NF_PUT_ATT_TEXT(fileidhout,varidy,
     &        'units',13,"degrees north")
      write(*,*) err

      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fileidhout,inqvarout,NF_FLOAT,2,dd,
     &       varidhout(k))
         write(*,*) err,inqvarout
         err=NF_PUT_ATT_TEXT(fileidhout,varidhout(k),
     &        'long_name',45,(EntPFT_title(k))//("- LAI"))
         write(*,*) err
         err=NF_PUT_ATT_TEXT(fileidhout,varidhout(k),
     &        'units',5,"m2/m2")
         write(*,*) err
         err=NF_PUT_ATT_REAL(fileidhout, varidhout(k), 
     &          '_FillValue',NF_FLOAT, 1, undef)
         write(*,*) err
      enddo
      err=NF_ENDDEF(fileidhout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fileidhout,varidx,1,IMout,longi)
      write(*,*) err
      err=NF_PUT_VARA_REAL(fileidhout,varidy,1,JMout,lati)
      write(*,*) err

      
      DLATA = (180./latin) * 60.
      OFFIA = 0.0
      DLATB = (180./latout) * 60.
      OFFIB = 0.0
      DATMIS = undef

      do k=1,18

         start2D(1)=1
         start2D(2)=1
         count2D(1)=IM1km
         count2D(2)=JM1km
 
         inqvarout = 'lai_'//EntPFT_shorttitle(k)
         err = NF_GET_VARA_REAL(fileidhin,varidhin(k),
     &        start2D,count2D,hin)
         write(*,*) err, 'get', k

         hout(:,:) = undef

         do i=1,longin
            do j=1,latin
               if (hin(i,j).eq.0) then
                  hin(i,j) = undef
               endif
            end do
         end do

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
         call HNTR4P(WTIN, hin, hout)
         write(*,*) 'Finished HNTR40'
         
         start2D(1)=1
         start2D(2)=1
         count2D(1)=IMout
         count2D(2)=JMout
         
         err=NF_PUT_VARA_REAL(fileidhout,varidhout(k),
     &        start2D,count2D,hout)
         write(*,*) err, 'put', k
         
      enddo

      err = NF_CLOSE(fileidhin)
      write(*,*) err, 'close'

      err = NF_CLOSE(fileidhout)
      write(*,*) err, 'close'

     
      
      end program aggregating_laimax
      

