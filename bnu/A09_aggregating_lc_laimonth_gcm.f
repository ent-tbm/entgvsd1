!* PROGRAM TO AGGREGATE LAND COVER AND LAIMAX from 1kmx1km to 0.25x0.25
!* Author: Carlo Montes: carlo.montes@nasa.gov
     
!* GFORTRAN COMPILATION ON DISCOVER-SP3
      
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A05_aggregating_lc_laimax.f
!  gfortran -o myExe_height arrayutil.o A05_aggregating_lc_laimax.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
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
      subroutine calc_lon_lat_144x90(IM2,JM2,lon,lat)
      implicit none
      integer,intent(in) :: IM2,JM2
      real*4,intent(inout) :: lon(IM2),lat(JM2)
      integer :: i,j

      do i=1,IM2
         lon(i) = -180.0 + 2.5*i - 1.25
      end do

      do j=1,JM2
         lat(j) = -90.0 + 2*j - 1
      end do
      end subroutine calc_lon_lat_144x90
      
!---------------------------------------------------------------------------------------------
      
      
      program aggregating_lc_laimonth_gcm

      implicit none

      include 'netcdf.inc'


!***************************************************
!*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
      character*18, parameter :: EntPFT_layers(40) =
     &     (/
     &     'water_lc         ',
     &     'ever_br_early_lc ',
     &     'ever_br_late_lc  ',
     &     'ever_nd_early_lc ',
     &     'ever_nd_late_lc  ',
     &     'cold_br_early_lc ',
     &     'cold_br_late_lc  ',
     &     'drought_br_lc    ',
     &     'decid_nd_lc      ',
     &     'cold_shrub_lc    ',
     &     'arid_shrub_lc    ',
     &     'c3_grass_per_lc  ',
     &     'c4_grass_lc      ',
     &     'c3_grass_ann_lc  ',
     &     'c3_grass_arct_lc ',
     &     'crops_c3_herb_lc ',
     &     'crops_c4_herb_lc ',
     &     'crops_woody_lc   ',
     &     'snow_ice_lc      ',
     &     'bare_sparse_lc   ',
     &     'water_lai        ',
     &     'ever_br_early_lai',
     &     'ever_br_late_lai ',
     &     'ever_nd_early_lai',
     &     'ever_nd_late_lai ',
     &     'cold_br_early_lai',
     &     'cold_br_late_lai ',
     &     'drought_br_lai   ',
     &     'decid_nd_lai     ',
     &     'cold_shrub_lai   ',
     &     'arid_shrub_lai   ',
     &     'c3_grass_per_lai ',
     &     'c4_grass_lai     ',
     &     'c3_grass_ann_lai ',
     &     'c3_grass_arct_lai',
     &     'crops_c3_herb_lai',
     &     'crops_c4_herb_lai',
     &     'crops_woody_lai  ',
     &     'snow_ice_lai     ',
     &     'bare_sparse_lai  '
     &     /)

      character*3, parameter :: MONTH(12) =
     &     (/
     &     "Jan","Feb","Mar","Apr","May","Jun",
     &     "Jul","Aug","Sep","Oct","Nov","Dec"
     &     /)
      
      integer, parameter :: IM25 = 1440 !long at 0.25 degrees
      integer, parameter :: JM25 = 720 !lat at 0.25 degrees
      integer, parameter :: IM2 = 144 !long at 2.0 degrees
      integer, parameter :: JM2 = 90 !lat at 2.5 degrees

      character*(*), parameter :: res_in="1440x720"
      character*(*), parameter :: res_out="144x90"

      integer, parameter :: NUMLAYERS = 40

      integer, parameter :: longin = IM25
      integer, parameter :: latin = JM25
      integer, parameter :: longout = IM2
      integer, parameter :: latout = JM2
      integer :: err

      character*80 :: filein, fileout
      integer :: varidin
      
      real*4, parameter :: undef = -1e30

      real*4 :: lon(longout)
      real*4 :: lat(latout)

      integer :: k,p,i,j,m
      real*4 :: layerin(longin,latin)
      real*4 :: layerout(NUMLAYERS,longout,latout)
      real*4 :: layeroutpre(longout,latout)
      
      character*60 :: PathFilepre, PathFilepost
      integer :: start2d(2),count2d(2)
      character*40 :: inqvar
      
      real*4 :: WTIN(longin,latin)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS

      integer :: dd(2),dimidx,dimidy,fileidin,fileidout
      integer :: varidx,varidy,varidout(40)
      

! --------------------------------------------------------

      call calc_lon_lat_144x90(longout,latout,lon,lat)
      
      layerout(:,:,:) = 0.0
      
      DLATA = (180./latin) * 60.
      OFFIA = 0.0
      DLATB = (180./latout) * 60.
      OFFIB = 0.0
      DATMIS = undef

      
      !     MAIN LOOP

      do m = 1,12
         
!     INPUT FILE
      
         PathFilepre= '../lc_lai_ent/'
         PathFilepost= 'EntMM_lc_lai_'//MONTH(m)//'_'//res_in//'.nc '
         filein = trim(PathFilepre)//trim(PathFilepost)
         err = NF_OPEN(filein,NF_WRITE,fileidin)

!     OUTPUT FILE

         PathFilepre= '../lc_lai_ent/'
         PathFilepost= 'EntMM_lc_lai_'//MONTH(m)//'_'//res_out//'.nc '
         fileout = trim(PathFilepre)//trim(PathFilepost)
         err = NF_CREATE(fileout,NF_CLOBBER,fileidout)
         !write(*,*) err, 'nf_create out ',trim(fileout)
         err=NF_DEF_DIM(fileidout,'lon',longout,dimidx)
         err=NF_DEF_DIM(fileidout,'lat',latout,dimidy)
	 
         err=NF_DEF_VAR(fileidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(fileidout,'lat',NF_REAL,1,dimidy,varidy)

         do k = 1,NUMLAYERS

            !write(*,*) 'main loop, layer = ', k,' of',NUMLAYERS

            inqvar = trim(EntPFT_layers(k))
            err = NF_INQ_VARID(fileidin,inqvar,varidin)

            start2d(1) = 1
            start2d(2) = 1
            count2d(1) = longin
            count2d(2) = latin
            
            err = NF_GET_VARA_REAL(fileidin,varidin,start2d,
     &           count2d,layerin)
            
            do i=1,longin
               do j=1,latin
                  if (layerin(i,j).eq.0) then
                     layerin(i,j) = undef
                  endif
               end do
            end do
            
            call HNTR40(longin,latin,OFFIA,DLATA,
     *           longout,latout,OFFIB,DLATB, DATMIS)
           
            do i=1,longin
               do j=1,latin
                  WTIN(i,j) = 1.
                  if (layerin(i,j).eq.undef) WTIN(i,j) = 0.
               end do
            end do
            
            call HNTR4P(WTIN, layerin, layeroutpre)
            
            do i=1,longout
               do j=1,latout
                  if (layeroutpre(i,j).eq.undef) then
                     layeroutpre(i,j) = 0
                  endif
               end do
            end do
            
	    layerout(k,:,:) = layeroutpre(:,:)


            dd(1)=dimidx
            dd(2)=dimidy
            err=NF_DEF_VAR(fileidout,inqvar,NF_REAL,2,dd,varidout(k))
            write(*,*) err,varidout(k)

         enddo

         err=NF_ENDDEF(fileidout)
              
         start2d(1) = 1
         start2d(2) = 1 
         count2d(1) = longout
         count2d(2) = latout
         
         do k = 1,NUMLAYERS   
            err = NF_PUT_VARA_REAL(fileidout,varidout(k),start2d,
     &           count2d,layerout(k,:,:))
            !write(*,*) err
	 enddo

         err = NF_PUT_VARA_REAL(fileidout,varidx,1,longout,lon)
         err = NF_PUT_VARA_REAL(fileidout,varidy,1,latout,lat)

         err = NF_CLOSE(fileidin)
         err = NF_CLOSE(fileidout)

         !write(*,*) err, 'nf_close out ',trim(fileout)
         
      enddo
      
      end program aggregating_lc_laimonth_gcm
      
