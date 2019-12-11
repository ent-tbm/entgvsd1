C**** hntr4_monfreda2008.f - 
C**** Author:  Nancy Y. Kiang
C****  Converts Monfreda et al. (2008) crops and pasture cover fractions
C****  at 5 min x 5 min to a normalized product giving crop type as fraction 
C****  of total crop cover at other desired spatial resolutions.
C****
C****  Original Monfreda et al (2008) crops & pasture fraction are at
C****  5 min x 5 min.  HNTR4 can convert from as high res as 2 minutes.
C****  Here, convert 0.5 degree x 0.5 degree or to 1 degree x 1 degree.
C****  Can modify this program to convert that file to lower res.
C****  Because Monfreda et al. counted multiple harvest crops mutliple 
C****  times in a year, the grid fractions add up to > 1.  For the Ent GVSD,
C****  this program outputs a normalized version of the Monfreda et al.
C****  values, estimating the fraction of crop cover for a crop type, rather
C****  that fraction of grid cell area.  These ratios then can be used to
C****  to apportion satellite-observed crop cover.

!* To compile on MacBook Pro, follow config file in Ent code.
!  g95 -fno-second-underscore -O2 -fendian=big -cpp hntr4_monfreda2008.f
!  gfortran -cpp -fconvert=big-endian -O0 -fno-range-check hntr4_monfreda2008.f
!* To avoid memory problems:
!  > ulimit -s 64000

C**** Uses:
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
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *     FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *     IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *     DATMCB, INA,JNA, INB,JNB
C**** 
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or.
     *     IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
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
C     WRITE (0,915) 'IMIN=',IMIN(1:IMB)
C     WRITE (0,915) 'IMAX=',IMAX(1:IMB)
C     WRITE (0,916) 'FMIN=',FMIN(1:IMB)
C     WRITE (0,916) 'FMAX=',FMAX(1:IMB)
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
C     WRITE (0,915) 'JMIN=',JMIN(1:JMB)
C     WRITE (0,915) 'JMAX=',JMAX(1:JMB)
C     WRITE (0,916) 'GMIN=',GMIN(1:JMB)
C     WRITE (0,916) 'GMAX=',GMAX(1:JMB)
            Return
C**** 
C**** Invalid parameters or dimensions out of range
C**** 
 400        Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,
     &           DLATB, DATMIS
            Stop 400
C**** 
C     915 Format (/ 1X,A5 / (20I6))
C     916 Format (/ 1X,A5 / (20F6.2))
 940        Format ('0Arguments received by HNTRP0 in order:'/
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
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *     FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *     IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
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
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *     FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *     IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
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
      call HNTR4P(WTIN, LAYERIN, LAYEROUT)

      close(10)
      end subroutine read_layer

      subroutine Set_val(A,im,jm,valin,valout)
      !@ Set some unknown value in matrix A to valout (e.g. zero, undef).
      implicit none
      integer :: im, jm
      real*4 :: A(im,jm)
      real*4 :: valin,valout
      !---------
      integer :: i,j

      do i=1,im
         do j=1,jm
            if (A(i,j).eq.valin) A(i,j) = valout
         enddo
      enddo
      
      end subroutine Set_val

      function div_layers(A,B,undef) 
      implicit none
      real*4,dimension(:),intent(in) :: A
      real*4,dimension(:),intent(in) :: B
      real*4,intent(in) :: undef
      real*4,dimension(size(A)) :: div_layers
      !--------
      integer :: i
      real*4 :: LAYEROUT(size(A))

      do i=1,size(A)
         if ((B(i).eq.0.).or.(B(i).eq.undef)
     &        .or.(A(i).eq.undef)) then
            LAYEROUT(i) = 0.
         else
            LAYEROUT(i) = A(i)/B(i)
         endif
      enddo
      
      div_layers = (LAYEROUT)
      end function div_layers
!------------------------------------------------------------------------

      
      program hntr4_monfreda
!     Read in Monfreda et al. (2008) crop type cover for year 2000
!     5 minute x 5 minute, and use HNTRP* to 
       ! interpolate to coarser resolutions and make normalized version.
      implicit none

      !character*25, parameter :: pathin = '../../CROPS/Monfreda2008/'
      character*17, parameter :: pathin = '../5min_origdata/'
      character*0, parameter :: pathout = ''

      integer, parameter :: IM5M = 4320 !long at 5 minutes
      integer, parameter :: JM5M = 2160 !lat at 5 minutes
      integer, parameter :: IMH = 720 !long at 0.5 degrees
      integer, parameter :: JMH = 360 !lat at 0.5 degrees
      integer, parameter :: IM1 = 360 !long at 1 degrees
      integer, parameter :: JM1 = 180 !lat at 1 degrees
      integer, parameter :: IM2 = 144 !long at 2.5 degrees
      integer, parameter :: JM2 = 90 !lat at 2 degrees
      integer, parameter :: IM4X5 = 72 !long at 5 degrees
      integer, parameter :: JM4X5 = 46 !lat at 4 degrees
      
      integer, parameter :: longin = IM5M
      integer, parameter :: latin = JM5M

      !* EDIT OUTPUT RESOLUTION *!
      !integer, parameter :: longout = IM1
      !integer, parameter :: latout = JM1
      integer, parameter :: longout = IMH
      integer, parameter :: latout = JMH

      real*4, parameter :: undef_in = -9999.0
      real*4, parameter :: undef_A = -1.e30 !AViewer undef value
      
      character*80 :: TITLE
      character*256 :: filein, fileout, fileoutnorm,fileoutAviewer
      character*5 :: RESOUT
      real*4 :: LAYERIN(longin,latin)
      real*4 :: LAYEROUT(longout,latout)
      real*4 :: C3CROP(longout,latout),C4CROP(longout,latout)
      real*4 :: HERBCROP(longout,latout),SHRUBCROP(longout,latout)
      real*4 :: TREECROP(longout,latout)
      real*4 :: CROPTOT(longout,latout)
      real*4 :: HERBNORM(longout,latout), SHRUBNORM(longout,latout)
      real*4 :: TREENORM(longout,latout)
      real*4 :: C3NORM(longout,latout),C4NORM(longout,latout)
      real*4 :: C3CAP1(longout,latout),C4CAP1(longout,latout)
     &     ,HERBCAP1(longout,latout)
      interface
      function div_layers(A,B,undef)
      implicit none
      real*4,dimension(:),intent(in)::A,B
      real*4,intent(in) :: undef
      real*4,dimension(size(A)) :: div_layers
      end function div_layers
      end interface

!      real*4 :: WTIN(longin,latin)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS
      character :: lineskip
      integer :: i, j, skip
      real*4 :: temp

      LAYEROUT(:,:) = 1.E-30      
!DIVJ = 360. / NINT(360./latout)
      !DLATA = NINT(360./latin) * 60.
      DLATA = (180./latin) * 60. !THIS WORKS FOR 1x1
      OFFIA=0.0
      !DLATB = NINT(180./latout) * 60. !THIS WORKS FOR 1x1
      DLATB = NINT(180./latout *60.)  !THIS WORKS FOR 0.5x0.5
      OFFIB=0.0
!      DATMIS=-999999.
!      DATMIS=-9999.0  !PROVIDE MISSING DATA VALUE
      DATMIS = undef_A  !AViewer undef value

      !* Setup grids.
      write(*,*) longin,latin,OFFIA,DLATA,
     &     longout,latout,OFFIB,DLATB,DATMIS
      call HNTR40(longin,latin,OFFIA,DLATA,
     *    longout,latout,OFFIB,DLATB, DATMIS)

      LAYEROUT(:,:) = 0.

      !* Output files.
      if (longout.eq.IMH) then
         RESOUT = '05x05'
      elseif (longout.eq.IM1) then
         RESOUT= '1x1'
      else
         write(*,*) 'Not set up for this resolution,',longout,latout
         return
      endif

      fileout = 
     &     pathout//'Monfreda_crops_'//
     &     trim(RESOUT)//'.bin'
      write(*,*) fileout
      open(20, file=fileout, form='unformatted')
      fileoutnorm = 
     &     pathout//'Monfreda_crops_'//
     &     trim(RESOUT)//'_norm.bin'
      write(*,*) fileoutnorm
      open(30, file=fileoutnorm, form='unformatted')
      fileoutAviewer =
     &     pathout//'Monfreda_crops_'//
     &     trim(RESOUT)// '_normA.bin'
      write(*,*) fileoutAviewer
      open(50, file=fileoutAviewer, form='unformatted')

      !* Weights, if any
!      do i=1,longin
!         do j=1,latin
!            WTIN(i,j) = 1.      !FGRNDIN
!            !if (LAYERIN(i,j).le.-9999.) WTIN(i,j) = 0.
!         end do
!      end do


      !* Input file. - Format is ascii, no title, space delimited
      !* C3
      filein = pathin//'C3C4/L1C3C4.asc'
      call read_layer(filein, 6,
     &     latin,longin,latout,longout,undef_in,undef_A,LAYEROUT)
      call Set_val(LAYEROUT,longout,latout,undef_A,0.)
      TITLE = 'C3 CROPS 2000 (cover fraction) (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      C3CROP(:,:) = LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      !* C4
      filein = pathin//'C3C4/L2C3C4.asc'
      call read_layer(filein, 6,
     &     latin,longin,latout,longout,undef_in,undef_A,LAYEROUT)
      call Set_val(LAYEROUT,longout,latout,undef_A,0.)
      TITLE = 'C4 CROPS 2000 (cover fraction) (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      C4CROP(:,:) = LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      !* Herbaceous
      filein = pathin//'form/L1form.asc'
      call read_layer(filein, 6,
     &     latin,longin,latout,longout,undef_in,undef_A,LAYEROUT)
      call Set_val(LAYEROUT,longout,latout,undef_A,0.)
      TITLE = 'HERB CROPS 2000 (cover fraction) (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      HERBCROP(:,:) = LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT
      
      !* Shrub crops
      filein = pathin//'form/L2form.asc'
      call read_layer(filein, 6,
     &     latin,longin,latout,longout,undef_in,undef_A,LAYEROUT)
      call Set_val(LAYEROUT,longout,latout,undef_A,0.)
      TITLE = 'SHRUB CROPS 2000 (cover fraction) (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      SHRUBCROP(:,:) = LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      !* Tree crops
      filein = pathin//'form/L3form.asc'
      call read_layer(filein, 6,
     &     latin,longin,latout,longout,undef_in,undef_A,LAYEROUT)
      call Set_val(LAYEROUT,longout,latout,undef_A,0.)
      TITLE = 'TREE CROPS 2000 (cover fraction) (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      TREECROP(:,:) = LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT
      

      !* Sum up C3 + C4.
      LAYEROUT(:,:) = 0.d0
      LAYEROUT(:,:) = C3CROP(:,:) + C4CROP(:,:) !Will sum>1
      TITLE = 
     & 'C3+C4 CROPS 2000 (cover fraction multi) (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      CROPTOT(:,:) = LAYEROUT(:,:)
      do i=1,longout
         do j=1,latout
            if (CROPTOT(i,j).eq.0.0) then
               C3NORM(i,j) = 0.0
               C4NORM(i,j) = 0.0
            else if (CROPTOT(i,j).le.undef_A) then
               CROPTOT(i,j)= undef_A
               C3NORM(i,j) = undef_A
               C4NORM(i,j) = undef_A
            else
               C3NORM(i,j) = C3CROP(i,j)/CROPTOT(i,j)
               C4NORM(i,j) = C4CROP(i,j)/CROPTOT(i,j)
            endif
         enddo
      enddo
      TITLE = 
     &   'C3 CROPS (crop cover fraction) '//
     &' (normalized Monfreda et al. 2008)'
      write(30) TITLE, C3NORM
      LAYEROUT = C3NORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 
     &   'C4 CROPS (crop cover fraction)'//
     &' (normalized Monfreda et al. 2008)'
      write(30) TITLE, C4NORM
      LAYEROUT = C4NORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 
     & 'C3+C4 CROPS (crop cover fraction)'//
     &' (normalized Monfreda et al. 2008)'
      write(30) TITLE, C3NORM + C4NORM
      LAYEROUT = C3NORM + C4NORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT



      !* Sum up herb, shrub, tree
      LAYEROUT(:,:) = 0.d0
      LAYEROUT(:,:) = HERBCROP(:,:) + SHRUBCROP(:,:)+ TREECROP(:,:)
      TITLE = 
     & 'HERB+SHRUB+TREE CROPS 2000 (crop cover fraction)'//
     & ' (Monfreda et al. 2008)'
      write(20) TITLE, LAYEROUT
      CROPTOT(:,:) = LAYEROUT(:,:)
      do i=1,longout
         do j=1,latout
            if (CROPTOT(i,j).eq.0.0) then
               HERBNORM(i,j) = 0.0
               SHRUBNORM(i,j) = 0.0
               TREENORM(i,j) = 0.0
            elseif (CROPTOT(i,j).le.undef_A) then
               CROPTOT(i,j)=undef_A
               HERBNORM(i,j) = undef_A
               SHRUBNORM(i,j) = undef_A
               TREENORM(i,j) = undef_A
            else
               HERBNORM(i,j) = HERBCROP(i,j)/CROPTOT(i,j)
               SHRUBNORM(i,j) = SHRUBCROP(i,j)/CROPTOT(i,j)
               TREENORM(i,j) = TREECROP(i,j)/CROPTOT(i,j)
            endif
         enddo
      enddo
      TITLE = 
     & 'HERB CROPS (crop cover fraction)'//
     &' (normalized Monfreda et al. 2008)'
      write(30) TITLE, HERBNORM
      LAYEROUT = HERBNORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 
     & 'SHRUB CROPS (crop cover fraction)'//
     &'(normalized Monfreda et al. 2008)'
      write(30) TITLE, SHRUBNORM
      LAYEROUT = SHRUBNORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 
     & 'TREE CROPS (crop cover fraction)'//
     &'  (normalized Monfreda et al. 2008)'
      write(30) TITLE, TREENORM
      LAYEROUT = TREENORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 
     & 'HERB+SHRUB+TREE CROPS (crop cover fraction)'//
     &'  (normalized Monfreda et al. 2008)'
      write(30) TITLE, HERBNORM + SHRUBNORM + TREENORM
      LAYEROUT = HERBNORM + SHRUBNORM + TREENORM
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      LAYEROUT(:,:) = 0.
      do i=1,longout
         do j=1,latout
            if ((HERBNORM(i,j).eq.0.).or.(HERBNORM(i,j).eq.undef_A))
     &           then
               LAYEROUT(i,j) = HERBNORM(i,j)
            elseif (C4NORM(i,j).eq.undef_A) then
               LAYEROUT(i,j) = undef_A
            else
!               temp = (C4NORM(i,j)+C3NORM(i,j)
!     &              -SHRUBNORM(i,j)-TREENORM(i,j))
!               if (temp.le.0.) then 
!                  LAYEROUT(i,j) = 0.
!               else
!                  LAYEROUT(i,j) = C4NORM(i,j)/temp
!               endif
!!               LAYEROUT(i,j) = C4NORM(i,j)*HERBNORM(i,j)
               if (HERBCROP(i,j).le.0) then
                  LAYEROUT(i,j) = 0.
               else 
                  LAYEROUT(i,j) = C4CROP(i,j)/HERBCROP(i,j)
               endif
            endif
         enddo
      enddo
      TITLE = 'C4 CROPS (multi-harvest herb cover fraction) '// 
     &'(Monfreda et al. 2008)'
      write(30) TITLE, LAYEROUT

      do i=1,longout
         do j=1,latout
            if (LAYEROUT(i,j).eq.undef_A) LAYEROUT(i,j) = 0.
         enddo
      enddo
      write(30) TITLE, LAYEROUT

      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 'C3 HERB CROPS (cover fraction,multi-harvest)'//
     &'(Monfreda et al. 2008)'
      LAYEROUT = HERBCROP - C4CROP
      write(30) TITLE, LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 'C3 HERB CROPS (fraction of herb crop cover,'//
     &     'multi-harvest)'
      do i=1,longout
         do j=1,latout
            if (HERBCROP(i,j).le.0.0) then
               LAYEROUT(i,j) = 0.0
            elseif (C4CROP(i,j).le.0.) then
               LAYEROUT(i,j) = HERBCROP(i,j)
            else
               LAYEROUT(i,j) = (HERBCROP(i,j) - C4CROP(i,j))
     &              /HERBCROP(i,j)
            endif
         enddo
      enddo
      write(30) TITLE, LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      TITLE = 'SHRUB + TREE CROPS (cover fraction)'
      do i=1,longout
         do j=1,latout
            if (SHRUBCROP(i,j).eq.undef_A) then
               LAYEROUT(i,j) = TREECROP(i,j)
            elseif (TREECROP(i,j).eq.undef_A) then
               LAYEROUT(i,j) = SHRUBCROP(i,j)
            else
               LAYEROUT(i,j) = SHRUBCROP(i,j) + TREECROP(i,j)
            endif
         enddo
      enddo
      write(30) TITLE, LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      do i=1,longout
         do j=1,latout
            if (C3CROP(i,j).ge.1.0) then
               C3CAP1(i,j) = 1.0
            else
               C3CAP1(i,j) = C3CROP(i,j)
            endif
            if (C4CROP(i,j).ge.1.0) then
               C4CAP1(i,j) = 1.0
            else
               C4CAP1(i,j) = C4CROP(i,j)
            endif
         enddo
      enddo
      TITLE = 'C3 + C4 crops (cover fraction capped at 1)'//
     &     ' CANNOT DO THIS, >1'
      LAYEROUT(:,:) = C3CAP1(:,:) + C4CAP1(:,:)
      write(30) TITLE, LAYEROUT
      call Set_val(LAYEROUT,longout,latout,0.,undef_A)
      write(50) TITLE, LAYEROUT

      close(20)
      close(30)
      close(50)
      end program hntr4_monfreda
      
