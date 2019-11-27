! GFORTRAN COMPILATION ON DISCOVER:

! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include hntr_soilalbedo.f
! gfortran -o myExeTrim arrayutil.o hntr_soilalbedo.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
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
!      real*8, parameter :: deltlon = 360./IM
!      real*8, parameter ::  detlat = 180./JM
      integer :: i,j

!      deltlon = 360./IM
!      detlat = 180./JM
     
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

!      character*38, parameter :: titlsoilalb(6) =
!     &     (/
!     &     'SOIL ALBEDO CARRER interim 2004 - VIS ',
!     &     'SOIL ALBEDO CARRER interim 2004 - NIR1',
!     &     'SOIL ALBEDO CARRER interim 2004 - NIR2',
!     &     'SOIL ALBEDO CARRER interim 2004 - NIR3',
!     &     'SOIL ALBEDO CARRER interim 2004 - NIR4',
!     &     'SOIL ALBEDO CARRER interim 2004 - NIR5'
!     &     /)

      real*4, parameter :: undef = -1.e30

      integer, parameter :: IMin=7200., JMin=3600.
      integer, parameter :: IMout=43200., JMout=21600.

      integer :: fileid,dimidx,dimidy,dimidz,dd(2),varidx,varidy,varidz
      integer startB(2), countB(2)
      character*256 :: filein,fileout
      character*20 :: inqvarin, inqvarout
      character*80 :: title
      character*37 :: long_title

      integer :: err

      real*4 :: WTIN(IMin,JMin)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS

!      real*4, ALLOCATABLE :: albvisin(:,:)
!      real*4, ALLOCATABLE :: albnirin(:,:)
!      real*4, ALLOCATABLE :: albvisout(:,:)
!      real*4, ALLOCATABLE :: albnirout(:,:)

      real*4 :: albvisin(IMin,JMin)
      real*4 :: albnirin(IMin,JMin)
      real*4 :: albswin(IMin,JMin)
      real*4 :: albvisout(IMout,JMout)
      real*4 :: albnirout(IMout,JMout)
      real*4 :: albswout(IMout,JMout)

      integer i, j, k, io, in, jn, maxpft, kx, m
      real*8 lat
      integer count

      integer :: ncidin,ncidout,varid,status
      real*4 :: long(IMout),lati(JMout)
            
!      ALLOCATE( albvisin(IMin,JMin) )
!      ALLOCATE( albvisout(IMout,JMout) )
!      ALLOCATE( albnirin(IMin,JMin) )
!      ALLOCATE( albnirout(IMout,JMout) )
 
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
      ! albedo VIS
      filein = 'VIS_Alb_soil_yearly.006.forEnt.nc'
      status = nf_open(filein,0,ncidin)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'mean'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albvisin(:,:) = undef
      status = nf_get_var_real(ncidin,varid,albvisin)
      write(*,*) status,'nf_get_var_real ',inqvarin
      albvisout(:,:) = undef
      do i=1,IMin
         do j=1,JMin
            WTIN(i,j) = 1.
            if (albvisin(i,j).eq.undef) WTIN(i,j) = 0.
         end do
      end do
      write(*,*) shape(WTIN), shape(albvisin), shape(albvisout)
      call HNTR4P(WTIN, albvisin, albvisout)
      do i = 1,IMout
         do j = 1,JMout
	    if (albvisout(i,j).lt.0) albvisout(i,j) = 0
         enddo
      enddo

      err = NF_CREATE('VIS_Alb_soil_yearly.006.1kmx1km.nc',NF_CLOBBER,
     &     ncidout)
      write(*,*) err, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IMout,dimidx)
      write(*,*) err, 'dimidx',dimidx
      err=NF_DEF_DIM(ncidout,'lat',JMout,dimidy)
      write(*,*) err, 'dimidy',dimidy
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      write(*,*) err,'varidx',varidx
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      write(*,*) err,'varidy',varidy
      dd(1)=dimidx
      dd(2)=dimidy
      err=NF_DEF_VAR(ncidout,'mean',NF_REAL,2,dd,varid)
      write(*,*) err,'varid',varid
      long_title = 'SOIL ALBEDO CARRER interim 2004 - VIS'
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     trim(long_title))
      write(*,*) err
!      err=NF_PUT_ATT_REAL(ncidout,varid,'_FillValue',nf_real,1,undef)
!      write(*,*) err
      err=NF_ENDDEF(ncidout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IMout,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JMout,lati)
      write(*,*) err
      startB(1)=1
      startB(2)=1
      countB(1)=IMout
      countB(2)=JMout
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,albvisout)
      write(*,*) err, 'put'
      err = NF_CLOSE(ncidin)
      err = NF_CLOSE(ncidout)
      
      !------------------------------------------------------------------------
      ! albedo NIR
      filein = 'NIR_Alb_soil_yearly.006.forEnt.nc'
      status = nf_open(filein,0,ncidin)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'mean'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albnirin(:,:) = undef
      status = nf_get_var_real(ncidin,varid,albnirin)
      write(*,*) status,'nf_get_var_real ',inqvarin
      albnirout(:,:) = undef
      do i=1,IMin
         do j=1,JMin
            WTIN(i,j) = 1.
            if (albnirin(i,j).eq.undef) WTIN(i,j) = 0.
         end do
      end do
      write(*,*) shape(WTIN), shape(albnirin), shape(albnirout)
      call HNTR4P(WTIN, albnirin, albnirout)
      do i = 1,IMout
         do j = 1,JMout
	    if (albnirout(i,j).lt.0) albnirout(i,j) = 0
         enddo
      enddo

      err = NF_CREATE('NIR_Alb_soil_yearly.006.1kmx1km.nc',NF_CLOBBER,
     &     ncidout)
      write(*,*) err, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IMout,dimidx)
      write(*,*) err, 'dimidx',dimidx
      err=NF_DEF_DIM(ncidout,'lat',JMout,dimidy)
      write(*,*) err, 'dimidy',dimidy
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      write(*,*) err,'varidx',varidx
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      write(*,*) err,'varidy',varidy
      dd(1)=dimidx
      dd(2)=dimidy
      err=NF_DEF_VAR(ncidout,'mean',NF_REAL,2,dd,varid)
      write(*,*) err,'varid',varid
      long_title = 'SOIL ALBEDO CARRER interim 2004 - NIR'
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     trim(long_title))
      write(*,*) err
!      err=NF_PUT_ATT_REAL(ncidout,varid,'_FillValue',nf_real,1,undef)
!      write(*,*) err
      err=NF_ENDDEF(ncidout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IMout,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JMout,lati)
      write(*,*) err
      startB(1)=1
      startB(2)=1
      countB(1)=IMout
      countB(2)=JMout
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,albnirout)
      write(*,*) err, 'put'
      err = NF_CLOSE(ncidin)
      err = NF_CLOSE(ncidout)

      !------------------------------------------------------------------------
      ! albedo SW
      filein = 'SW_Alb_soil_yearly.006.nc'
      status = nf_open(filein,0,ncidin)
      write(*,*) status, 'nf_open in ', filein
      inqvarin = 'albedo'
      status = nf_inq_varid(ncidin,inqvarin,varid)
      write(*,*) status,'nf_inq_varid ',inqvarin
      albswin(:,:) = undef
      status = nf_get_var_real(ncidin,varid,albswin)
      write(*,*) status,'nf_get_var_real ',inqvarin
      albswout(:,:) = undef
      do i=1,IMin
         do j=1,JMin
            WTIN(i,j) = 1.
            if (albswin(i,j).eq.undef) WTIN(i,j) = 0.
         end do
      end do
      write(*,*) shape(WTIN), shape(albswin), shape(albswout)
      call HNTR4P(WTIN, albswin, albswout)
      do i = 1,IMout
         do j = 1,JMout
	    if (albswout(i,j).lt.0) albswout(i,j) = 0
         enddo
      enddo

      err = NF_CREATE('SW_Alb_soil_yearly.006.1kmx1km.nc',NF_CLOBBER,
     &     ncidout)
      write(*,*) err, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IMout,dimidx)
      write(*,*) err, 'dimidx',dimidx
      err=NF_DEF_DIM(ncidout,'lat',JMout,dimidy)
      write(*,*) err, 'dimidy',dimidy
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      write(*,*) err,'varidx',varidx
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      write(*,*) err,'varidy',varidy
      dd(1)=dimidx
      dd(2)=dimidy
      err=NF_DEF_VAR(ncidout,'mean',NF_REAL,2,dd,varid)
      write(*,*) err,'varid',varid
      long_title = 'SOIL ALBEDO CARRER interim 2004 - SW'
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     trim(long_title))
      write(*,*) err
!      err=NF_PUT_ATT_REAL(ncidout,varid,'_FillValue',nf_real,1,undef)
!      write(*,*) err
      err=NF_ENDDEF(ncidout)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IMout,long)
      write(*,*) err
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JMout,lati)
      write(*,*) err
      startB(1)=1
      startB(2)=1
      countB(1)=IMout
      countB(2)=JMout
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,albswout)
      write(*,*) err, 'put'
      err = NF_CLOSE(ncidin)
      err = NF_CLOSE(ncidout)
      !------------------------------------------------------------------------
      ! ij file

!      fileout = 'Carrer2014_soil_albedo_means_2004_1kmX1km.ij'

!      open(90,file=trim(fileout),
!     &   form='unformatted',status='unknown')
!      title = trim(titlsoilalb(1))
!      write(90) trim(title), albvisout(:,:)
!      print *, trim(title)
!      do k=1,5
!         title = trim(titlsoilalb(k+1))
!         write(90) trim(title), albnirout(:,:)
!         print *, trim(title)
!      enddo
!      close(90)

      !------------------------------------------------------------------------

!      DEALLOCATE(albnirin)
!      DEALLOCATE(albnirout)


      end program hntr_soilalbedo

