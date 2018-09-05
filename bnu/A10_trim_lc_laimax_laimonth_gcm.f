! This program converts EntMM 17 PFTs to Ent 16 PFTs + bright + dark.
! Converts lc, laimax, monthly lc and lai, and Simard heights.
! Original from trim_EntMM_monthly_noht.f, which did not do heights.
! It does NOT interpolate from fine to coarse grid -- this should be done
!  prior to using this program, so input and output are the same resolution.
! To change input/output resolution edit lines below 
! "define input file resolution" and "new (interpolated) values"
! To combine C3 and C4 crops: #define COMBINE_CROPS_C3_C4
! 9/12/13 Fixed bright/dark soil:  need to do partitioning after
!         each step of trim/scale/no crops to account for new
!         bare soil cover added, especially in crop grid cells.
! 1/16/13 Added subroutine replace_crops for nocrops to replace cover and LAI 
!         of crops with dominant natural veg in grid cell. If no natural veg,
!         then searches adjacent grid cells.
!         Added subroutine fill_crops for _ext version of crop LAI, to fill
!         in crop LAI in some grid cells that have no crops in MODIS cover
!         but may have crop cover in Pongratz historical cover.  Called once
!         for single grid extension.  Can be called again for further filling.
! 3/17/14 Added ext1 files at end for lai max, lai monthly, and height, by
!         replacing 15-crops in laic, laim, hm, and hsd with the ext values for
!         crops.
! Compile the program with:
!
! ifort -cpp convert_VEG5.f -convert big_endian
!

! GFORTRAN COMPILATION ON DISCOVER:
! ulimit -s unlimited
! module purge
! module load other/comp/gcc-4.9.2-sp3
! module load other/ncl-6.3.0
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f
! gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include A07_trim_lc_laimax_laimonth.f
! gfortran -o myExeTrim arrayutil.o A07_trim_lc_laimax_laimonth.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf
! ./myExeTrim

! BEFORE RUNNING: mkdir ../lc_lai_ent16

#define COMBINE_CROPS_C3_C4
#define SPLIT_BARE_SOIL

     

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

      module conversions
      implicit none
      
      character*3, parameter :: MONTH(12) =
     &     (/
     &     "Jan","Feb","Mar","Apr","May","Jun",
     &     "Jul","Aug","Sep","Oct","Nov","Dec"
     &     /)
      
      contains

      subroutine convert_vf(vf1, lai1, vf2, lai2, laimin)
      !Converts LAI that is on BARE to a vegetated fraction.
      !Reduces vf1, increases vf2 by increment that maintains lai2,
      ! or only increases lai2 if lai1 is so big that the cover-weighted
      ! lai is bigger than lai2.
      !vf1, lai1 - BARE cover fraction and LAI
      !vf2, lai2 - vegetated cover fraction and max LAI
      real*4 vf1, lai1, vf2, lai2, laimin
      !---
      real*4 tot_la, tot_vf, new_lai, new_vf

      tot_la = vf1*lai1 + vf2*lai2
      tot_vf = vf1 + vf2
      new_lai = tot_la/tot_vf
      !Fix for if 100% conversion to vf2.
      if (new_lai.eq.0.0) then
         new_vf = 0.0
      else
         new_lai = max( new_lai, laimin)
         new_vf = tot_la/new_lai
      endif
      ! get rid of round-off errors
      new_vf = min( new_vf, tot_vf )

      vf2 = new_vf
      lai2 = new_lai
      vf1 = tot_vf - vf2
      lai1 = 0.

      end subroutine convert_vf

!               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
!     &              vfm(m,i,j,9),vfm(m,i,j,9), vfc(i,j,9))

      subroutine convert_vfm(vf1, lai1, vf2, lai2, vfc)
      !convert monthly vf using vfc for new vf2 derived from laimax trim. 
      !vf1, lai1 - BARE cover fraction and LAI
      !vf2, lai2 - vegetated cover fraction and monthly LAI
      !vfc = vf2 + dvf1 from laimax trim (dvf1 is positive)
      !vf1 + vf2 = vfc + (vf1-dvf1)
      !vf1*lai1 = 0*(vf1-dvf1) + (dvf1*dlai1)
      !new_lai2 = ((dvf1*dlai1) + (vf2*lai2))/(dvf1+vf2)
      !         = ((vf1*lai1) + (vf2*lai2))/vfc
      real*4, intent(inout) :: vf1, lai1, vf2, lai2
      real*4, intent(in) :: vfc
      !---
      real*4 :: tot_lai, new_lai, new_vf

      new_vf = vfc
      if (new_vf.eq.0.0) then  
         print *,'new_vf is zero',vf1,vf2,lai1,lai2,vfc
         STOP
      endif
      new_lai = (vf1*lai1 + vf2*lai2)/new_vf
!      vf2 = new_vf
      lai2 = new_lai
      vf1 = vf1 - (vfc - vf2)
      vf2 = new_vf
      lai1 = 0.

      end subroutine convert_vfm


      subroutine convert_vfh(vf1, h1, hsd1, vf2, h2, hsd2, vfc)
      !convert height data, h and hsd, based on the converted vfc.
      !Since conversion is to shrub and crops, or, if none, to the
      ! next dominant PFT cover type, there is no cover-weighted height
      ! averaging, but the h and hsd values are simply preserved (this
      ! routine effectively does nothing, since there is no height
      ! over BARE of ICE in EntGVSDmosaic_05x05.ij.
      !However, for the last conversion of any remaining BARE laimax to
      ! arid-adapted shrub cover that does not pre-exist, a new height
      ! has to be assigned.  ##DOUBLEBCHECK simard.f for height assigned
      ! to arid shrubs!
      !vf1, h1, hsd1 - BARE cover fraction, height, and height stdev
      !vf2, h2, hsd2  - vegetated cover fraction, heigth, and height stdev
      !vfc - new vegetated cover fraction that was derived from laimax.
      !vfc = vf2 + dvf1 from laimax trim (dvf1 is positive)
      !vf1 + vf2 = vfc + (vf1-dvf1)
      !vf1*h1 = 0*(vf1-dvf1) + (dvf1*dlai1)
      !new_lai2 = ((dvf1*dlai1) + (vf2*lai2))/(dvf1+vf2)
      !         = ((vf1*lai1) + (vf2*lai2))/vfc
      real*4, intent(inout) :: vf1, h1, hsd1, vf2, h2, hsd2
      real*4, intent(in) :: vfc
      !---
      real*4 tot_lai, new_h, new_hsd, new_vf2

      new_vf2 = vfc
      if (new_vf2>0.) then !vf2 is present, keep h2
         new_h = h2
         new_hsd = hsd2
      else
         new_h = h1
         new_hsd = hsd1
      endif
      vf2 = new_vf2
      h2 = new_h
      vf1 = vf1 - (vfc - vf2)
      h1 = 0.
      hsd1 = 0.

      end subroutine convert_vfh


      subroutine write_output_lai(titlec, laic, n, fileprefix
     &     , MISC, MON, resoutt)
      !Same as write_output, but only lai, not lc.
      character*80 :: titlec(:)
      real*4 :: laic(:,:,:)
      integer :: n
      character*(*) :: fileprefix
      character*(*) :: MISC, MON
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      character*4  :: MONstr
      integer :: k

      if (trim(MON)=="") then
         MONstr=""
      else
         MONstr="_"//MON
      endif

      open(90,file=fileprefix//"_lai_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!        read( titlec(k), '(a4, a46)' ) dum, pft_name
!        write( title, '(i2, " - ", a46, " (LAI) 4x5")' )
!     &       k, adjustl(pft_name)
!        write(90) title, laic(:,:,k)
         title = titlec(k)(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
         write(90) title, laic(:,:,k)
      enddo

      close(90)

      end subroutine write_output_lai


      subroutine write_output_single(titlefoo, vf, lai, fileprefix
     &     , MISC, MON, resoutt)
      !Same as write-output, but arrays are for only one veg type.
      character*80 :: titlefoo
      real*4 :: vf(:,:), lai(:,:)
      character*(*) :: fileprefix
      character*(*) :: MISC, MON
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title
      character*4  :: MONstr
      integer :: k

      if (trim(MON)=="") then
         MONstr=""
      else
         MONstr="_"//MON
      endif

      open(80,file=fileprefix//"_lc_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")


      title = titlefoo(1:48)//"  "//trim(MON)//" (cover)  "//resoutt
      write(80) title, vf(:,:)
      close(80)

      open(90,file=fileprefix//"_lai_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")



      title = titlefoo(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
      write(90) title, lai(:,:)

      close(90)

      end subroutine write_output_single


      subroutine write_output(titlec, vfc, laic, n, fileprefix
     &     , MISC, MON, resoutt)
      character*80 :: titlec(:)
      real*4 :: vfc(:,:,:), laic(:,:,:)
      integer :: n
      character*(*) :: fileprefix
      character*(*) :: MISC, MON
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      character*4  :: MONstr
      integer :: k

      if (trim(MON)=="") then
         MONstr=""
      else
         MONstr="_"//MON
      endif

      open(80,file=fileprefix//"_lc_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!        read( titlec(k), '(a4, a46)' ) dum, pft_name
!        write( title, "(i2, ' - ', a46, ' (cover) 4x5')" )
!     &       k, adjustl(pft_name)
!        titlevar = trim('cover) '//resout)
!        write( title, "(i2, ' - ', a46, titlevar)" )
!     &       k, adjustl(pft_name)
!        write(80) title, vfc(:,:,k)
         title = titlec(k)(1:48)//"  "//trim(MON)//" (cover)  "//resoutt
         write(80) title, vfc(:,:,k)
      enddo
      close(80)

      open(90,file=fileprefix//"_lai_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!        read( titlec(k), '(a4, a46)' ) dum, pft_name
!        write( title, '(i2, " - ", a46, " (LAI) 4x5")' )
!     &       k, adjustl(pft_name)
!        write(90) title, laic(:,:,k)
         title = titlec(k)(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
         write(90) title, laic(:,:,k)
      enddo

      close(90)


      end subroutine write_output


      subroutine write_output_h_single(titleh, h,hsd, filename
     &     , MISC, resoutt)
      !Same as write_output_h, but single arrays
      character*80 :: titleh(:)
      real*4 :: h(:,:), hsd(:,:)
      character*(*) :: filename
      character*(*) :: MISC
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      integer :: k

      open(90,file=filename,
     &     form="unformatted",status="unknown")

      title = titleh(1)(1:63)//"  "//trim(resoutt)
      write(90) title, h(:,:)

      title = titleh(2)(1:63)//"  "//trim(resoutt)
      write(90) title, hsd(:,:)

      close(90)
      end subroutine write_output_h_single


      subroutine write_output_h(titleh, h,hsd, n, filename
     &     , MISC, resoutt)
      character*80 :: titleh(:,:)
      real*4 :: h(:,:,:), hsd(:,:,:)
      integer :: n
      character*(*) :: filename
      character*(*) :: MISC
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      integer :: k

      open(90,file=filename,
     &     form="unformatted",status="unknown")

      do k=1,n
!         title = titleh(1,k)(1:48)//" "//trim(MISC)//
!     &        "   height (m) "//resout
         title = titleh(1,k)(1:63)//"  "//trim(resoutt)
         write(90) title, h(:,:,k)
      enddo
      do k=1,n
!         title = titleh(2,k)(1:48)//" "//trim(MISC)//
!     &        "   stdev (m) "//resout
         title = titleh(2,k)(1:63)//"  "//trim(resoutt)
         write(90) title, hsd(:,:,k)
      enddo

      close(90)
      end subroutine write_output_h



      subroutine get_bare_soil_brightratio(IMn,JMn, filename
     &     ,bs_brightratio)
      ! read bare soil from the "old"  dataset, 
      ! compute ratio of bright/total bare soil.
      ! this will be used to compute bright and dark cover fractions
      !   so that their sum preserves albedo of the bare soil.
      ! extend it to all cells with no data
      implicit none
      integer, intent(in) :: IMn, JMn
      character*(*) :: filename
      real*4, intent(out) :: bs_brightratio(IMn,JMn) !Fraction of bare that is bright.
      !---
      ! bare soil data from "old" dataset
      real*4 bsf(IMn,JMn), bsf_1(IMn,JMn), bsf_0(IMn,JMn) !1-bright cov, 0-dark cov
      character*80 :: title_bs
      integer count, i,j, ii,jj, k
      real*4 :: a, s, s1


      open(1,file=filename,
     &     form="unformatted",status="old")

      read(1) title_bs, bsf_1(:,:) !BRIGHT cover
      do k=2,9
        read(1)
      enddo
      read(1) title_bs, bsf_0(:,:) !DARK cover
      close(1)

      do j=1,JMn
        do i=1,IMn
          bsf(i,j) = bsf_1(i,j) + bsf_0(i,j)
          if ( bsf(i,j) > 0. ) then
            bs_brightratio(i,j) = bsf_1(i,j) /  bsf(i,j) 
          endif
        enddo
      enddo

      do
        count = 0
        
        do j=1,JMn
          do i=1,IMn
            if ( bsf(i,j) <= 0. ) then
              a = 0.
              s = 0.
              s1 = 0.
              do jj=max(1,j-1),min(JMn,j+1)
                do ii=max(1,i-1),min(IMn,i+1)
                  a = a + bs_brightratio(ii,jj)*bsf(ii,jj)
                  s = s + bsf(ii,jj)
                  s1 = s1 + 1.
                enddo
              enddo
              if ( s > 0 ) then
                bs_brightratio(i,j) = a/s
                bsf(i,j) = s/s1
              else
                count = count + 1
              endif
            endif
          enddo
        enddo

        print *, "count=", count
        if ( count== 0 ) exit
      enddo
      title_bs = "bare soil brightratio"
      write(91) title_bs, bs_brightratio
      title_bs = "bare soil fraction"
      write(91) title_bs, bsf


      end subroutine get_bare_soil_brightratio


      subroutine split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
      !Split BARE soil into BRIGHT and DARK cover to preserve albedo from
      !  "old" ModelE cover.  Should be called after each trim, scale, nocrops.
      !Any LAI on BARE soil should already have been moved to vegetated cover,
      !  so laic(:,:,N_BARE) should be zero.
      !This checks for cases if BARE is original total or was previously split.

      integer, intent(in) :: N_VEG, IMn,JMn,KM
      integer, intent(inout) :: N_BARE
      real*4, intent(in) :: bs_brightratio(:,:) !(IMn,JMn) !Fraction of bare that is bright.
      real*4 vfc(:,:,:), laic(:,:,:) !(IMn,JMn,KM)
      real*4 vfm(:,:,:,:), laim(:,:,:,:) !(12,IMn,JMn,KM)
      real*4 vfh(:,:,:) !(IMn,JMn,KM),
      real*4 hm(:,:,:), hsd(:,:,:) !(IMn,JMn,KM)
      character*80 :: titlec(:) !(18)
      character*80 :: titlem(:,:) !(12,18)
      character*80 :: titleh(:,:) !(2,18) !1-h, 2-hsd
      character*(*) :: res_out
      !-----Local----
      real*4 :: vfc_bare(IMn,JMn), vf_tot(IMn,JMn)
      integer :: m

      if ((N_BARE-N_VEG).eq.1) then !First time splitting
      ! increase number of fractions by 1 so that
      ! vfc(:,:,N_BARE-1) = fraction of bright soil
      ! vfc(:,:,N_BARE) = fraction of dark soil
         vf_tot(:,:) = vfc(:,:,N_BARE)  !All bare soil orig in N_BARE layer
         vfc(:,:,N_BARE) = vf_tot(:,:)*bs_brightratio(:,:)
         vfc(:,:,N_BARE+1) = vf_tot(:,:) - vfc(:,:,N_BARE)
         laic(:,:,N_BARE:N_BARE+1) = 0.
         do m=1,12
            vf_tot(:,:) = vfm(m,:,:,N_BARE)
            vfm(m,:,:,N_BARE) = vf_tot(:,:)*bs_brightratio(:,:)
            vfm(m,:,:,N_BARE+1) = vf_tot(:,:) - vfm(m,:,:,N_BARE)
            laim(m,:,:,N_BARE:N_BARE+1) = 0.
         enddo
         N_BARE = N_BARE+1
      else  !N_BARE-NVEG.eq.2, previously split and incremented N_BARE
         vf_tot(:,:) = vfc(:,:,N_BARE-1) + vfc(:,:,N_BARE)
         vfc(:,:,N_BARE-1) = vf_tot(:,:)*bs_brightratio(:,:)
         vfc(:,:,N_BARE) = vf_tot(:,:) - vfc(:,:,N_BARE-1)
         laic(:,:,N_BARE-1:N_BARE) = 0.
         do m=1,12
            vf_tot(:,:) = vfm(m,:,:,N_BARE-1) + vfm(m,:,:,N_BARE)
            vfm(m,:,:,N_BARE-1) = vf_tot(:,:)*bs_brightratio(:,:)
            vfm(m,:,:,N_BARE) = vf_tot(:,:) - vfm(m,:,:,N_BARE-1)
            laim(m,:,:,N_BARE-1:N_BARE) = 0.
         enddo
      endif

      vfh(:,:,N_BARE-1) = vfc(:,:,N_BARE-1)
      vfh(:,:,N_BARE) = vfc(:,:,N_BARE)
      hm(:,:,N_BARE) = 0.
      hsd(:,:,N_BARE) = 0.

      titlec(N_BARE-1) =   "17 - bright bare soil  "//
     &     "                           (COVER) "//res_out
      titlec(N_BARE) =     "18 - dark bare soil    "//
     &     "                           (COVER) "//res_out
      do m=1,12
            titlem(m,N_BARE-1) = "17 - bright bare soil  "//
     &           "                          "//MONTH(m)//" (cover) "
     &           //res_out
            titlem(m,N_BARE) =   "18 - dark bare soil    "//
     &           "                          "//MONTH(m)//" (cover) "
     &           //res_out
         enddo
      titleh(1,N_BARE-1) =   "17 - bright bare soil  "//
     &     "height (m)                                "//res_out
      titleh(1,N_BARE) =     "18 - dark bare soil    "//
     &     "height (m)                                "//res_out
      titleh(2,N_BARE-1) =   "17 - bright bare soil  "//
     &     "height stdev(m)                           "//res_out
      titleh(2,N_BARE) =     "18 - dark bare soil    "//
     &     "height stdev(m)                           "//res_out

      end subroutine split_bare_soil


      subroutine sum_lai_single(iu, IMn,JMn,KM,Ln, tag, vfs, lais)
      !Sum LAI over all cover for a single point in time. Output to iu.
      use arrayutil, only : IJAdd4, IJMult4, IJDiv4
      implicit none
      integer, intent(in) :: iu
      integer, intent(in) :: IMn,JMn,KM
      integer, intent(in) :: Ln !# layers to loop through
      character*(*) :: tag
      real*4 :: vfs(IMn,JMn,KM), lais(IMn,JMn,KM)
      !---Local---

      integer :: i,j,k,m
      character*80 :: TITLE
      real*4 :: L1(IMn,JMn)

      L1(:,:) = 0.0
      do k=1,Ln
            L1(:,:) = IJADD4(1,IMn,1,JMn,L1(:,:),
     &           IJMULT4(1,IMn,1,JMn,vfs(:,:,k),lais(:,:,k)))
      enddo

      TITLE = trim(tag)//'   LAI'
      write(iu) TITLE, L1(:,:)

      end subroutine sum_lai_single

      subroutine sum_lai(iu, IMn,JMn,KM,Ln, tag, vfmon, laimon)
      !Sum LAI over all cover by month.  Output to iu.
      use arrayutil, only : IJAdd4, IJMult4, IJDiv4
      implicit none
      integer, intent(in) :: iu
      integer, intent(in) :: IMn,JMn,KM,Ln
      character*(*) :: tag
      real*4 :: vfmon(12,IMn,JMn,KM), laimon(12,IMn,JMn,KM)
      !---Local---
      character*3, parameter :: MON(12) =
     &     (/
     &     "Jan","Feb","Mar","Apr","May","Jun",
     &     "Jul","Aug","Sep","Oct","Nov","Dec"
     &     /)

      integer :: i,j,k,m
      character*80 :: TITLE
      real*4 :: L1(12,IMn,JMn)

      L1(:,:,:) = 0.0
      do m=1,12
         do k=1,Ln
            L1(m,:,:) = IJADD4(1,IMn,1,JMn,L1(m,:,:),
     &           IJMULT4(1,IMn,1,JMn,vfmon(m,:,:,k),laimon(m,:,:,k)))
         enddo
      enddo

      do m=1,12
         TITLE = trim(tag)//' '//MON(m)//'   LAI'
         write(iu) TITLE, L1(m,:,:)
      enddo

      
      end subroutine sum_lai

      subroutine check_lc_lai_mismatch(KX,IX,JX,vf,lai,label,titlek)
      implicit none
      integer, intent(in) :: KX,IX,JX
      real*4, intent(in) :: vf(KX,IX,JX),lai(KX,IX,JX)
      character*3,intent(in) :: label
      character*80 :: titlek(:)
      !------
      integer :: i,j,k
      real*4 :: LAYER(IX,JX)
      character*80 :: TITLE

      do k=1,KX
         LAYER(:,:) = 0.
         do i=1,IX
            do j=1,JX
               if (((vf(i,j,k).eq.0.0).and.(lai(i,j,k).gt.0.0))
     &              .or.(vf(i,j,k).gt.0.0).and.(lai(i,j,k).eq.0.0))
     &              then
                  print *,label,' mismatch,k,i,j,cov,lai'
     &                 ,k,i,j,vf(i,j,k),lai(i,j,k)
                  LAYER(i,j) = vf(i,j,k)
               endif
            enddo
         enddo
         TITLE = label//' lc_lai mis '//titlek(k)(1:20)
         write(100) TITLE, LAYER
      enddo
      end subroutine check_lc_lai_mismatch


      subroutine replace_crops(IMn,JMn,KM,i,j,s,vfc,vfm,laic,laim,flag)
      !Replace crops in (i,j) with main natural cover in cell or adjacent.
      !The check for existence of crops in (i,j) is done before subroutine call
      !Fix 6 cells that are all crops due to MODIS error or two islands.
      !Assign LAI (cover-avg) of natural cover if from adjacent cells.
      !Do not zero crop LAI but keep for historical cover change.
      !NOTES:  Checks of adjacent cells will use replaced cover from previously
      ! processed cells (previous i,j).  I don't see this as a problem,
      ! since there won't be a consistent dominance either side of the cell.
      ! it helps reduce iteration for adjacent cells.

      integer :: IMn,JMn,KM,i, j
      real*4 :: s
      !real*4 :: vfct(:,:,:) !vfc trim before nocrops
      !real*4 :: vfmt(:,:,:,:)
      real*4 :: vfc(:,:,:) , laic(:,:,:)
      real*4 :: vfm(:,:,:,:), laim(12,IMn,JMn,KM)
      integer :: flag
      !----Local----
      integer :: m, k, ii, jj !cell to search for natural (i,j or adjacent)
      real*4 :: covmax
      integer :: covmaxk, covmaxii,covmaxjj
      integer :: naturalvegfound
      real*4 :: covsum(14), covavglai(14)
      real*4 :: covsumm(12,14),covavglaim(12,14)

      !First check just (i,j) cell for natural veg
      naturalvegfound = 0
      covmax = 0
      covmaxk = 0
      covmaxii = i
      covmaxjj = j
      do k=1,14                 !Find max non-crop, non-bare natural cover type
         if (vfc(i,j,k) > covmax)  then
            covmax = vfc(i,j,k)
            covmaxk = k
         endif
      enddo
      if (covmax.gt.0) then !Assign dominant natural veg to crop
         naturalvegfound = 1
         vfc(i,j,covmaxk) = vfc(i,j,covmaxk) 
     &        + vfc(i,j,15) + vfc(i,j,16)
         vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk)
     &        + vfm(:,i,j,15) + vfm(:,i,j,16)
         !Assign LAI in (i,j) from covmaxii,covmaxjj,covmaxk 
         laic(i,j,covmaxk) = laic(covmaxii,covmaxjj,covmaxk)
         laim(:,i,j,covmaxk) = laim(:,covmaxii,covmaxjj,covmaxk)
!         write(*,*) 'Cell has crops + natural',i,j,s
      else
!         write(*,*) "No natural veg in (i,j). Checking adjacent",i,j,s
         covmax = 0
         covmaxk = 0
         covmaxii = 0
         covmaxjj = 0
         covsum(:) = 0
         covavglai(:) = 0
         do ii=i-1,i+1
            do jj=j-1,j+1
               if ( (ii.ge.1).and.(ii.le.IMn)
     &              .and.(jj.ge.1).and.(jj.le.JMn) !in grid range
     &              .and.((ii.ne.i).or.(jj.ne.j)) ) !not the i,j center cell
     &              then
                  do k=1,14     
!                     !Find max non-crop, non-bare natural cover type
!                     if (vfc(ii,jj,k) > covmax)  then
!                        covmax = vfc(ii,jj,k)
!                        covmaxk = k
!                        covmaxii = ii
!                        covmaxjj = jj
!                     endif
                     !Sum adjacent natural veg cover by type.
                     covsum(k) = covsum(k) + vfc(ii,jj,k)
                     covavglai(k) = covavglai(k) + 
     &                    laic(ii,jj,k)*vfc(ii,jj,k)
                     covsumm(:,k) = covsumm(:,k) + vfm(:,ii,jj,k)
                     covavglaim(:,k) = covavglaim(:,k) +
     &                       laim(:,ii,jj,k)*vfm(:,ii,jj,k)
                  enddo
               endif
            enddo
         enddo
         covmax = 0
         covmaxk = 0
         do k=1,14 !Find largest adjacent natural cover
            if (covsum(k)>covmax) then
               covmax = covsum(k)
               covmaxk = k
            endif
         enddo
         if (covmax>0) then !Assign adjacent natural cover type and LAI to crop
            naturalvegfound = 1
            covavglai(k) = covavglai(k)/covsum(k)
            covavglaim(:,k) = covavglaim(:,k)/covsumm(:,k)
            vfc(i,j,covmaxk) = vfc(i,j,covmaxk) 
     &           + vfc(i,j,15) + vfc(i,j,16)
            vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk)
     &           + vfm(:,i,j,15) + vfm(:,i,j,16)
            !Assign LAI in (i,j) from covmaxii,covmaxjj,covmaxk 
            laic(i,j,covmaxk) = covavglai(covmaxk)
            laim(:,i,j,covmaxk) = covavglaim(:,covmaxk)
            !write(*,*) 'Cell or adjacent has crops + natural',i,j,s
!            write(*,*) 'Found natural veg in adjacent cells',i,j
         else
!            write(*,*) 'No natural veg in adjacent cells,',i,j
         endif
      endif

      !For IMn,JMn = 144x90
      !Fix remaining all-crop cells.  These are:  
      !- MODIS error in Antarctic, 
      !     (i,j) = (34,9),(41,9),(48,9),(34,10) -> C3 arctic grass
      !- Islands:  
      !  Mauritius (98,36) - sugar cane,tea,pasture,forest,savanna -> C4 grass 
      !  Nauru (139,45) - grassland bordered by tropicalforest -> C4 grass
      !Antarctic MODIS error
      if ((IMn.eq.144).and.(JMn.eq.90)) then
!      if ((IMn.eq.1440).and.(JMn.eq.720)) then
         if ( (((i.eq.34).or.(i.eq.41).or.(i.eq.48)).and.(j.eq.9))
     &        .or.( (i.eq.34).and.(j.eq.10) ) ) then !Antartic
            vfc(i,j,14) = vfc(i,j,14) + vfc(i,j,15) + vfc(i,j,16)
            vfm(:,i,j,14) = vfm(:,i,j,14) 
     &           + vfm(:,i,j,15) + vfm(:,i,j,16)
            laic(i,j,14) = (laic(i,j,15)*vfc(i,j,15) 
     &           + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
            do m=1,12
               laim(m,i,j,14) = ( laim(m,i,j,15)*vfm(m,i,j,15) +
     &              laim(m,i,j,16)*vfm(m,i,j,16) ) / 
     &              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
            enddo
            write(*,*) 'Replaced Antarctic crops.'
            naturalvegfound=1
         elseif ( ((i.eq.98).and.(j.eq.36))
     &           .or.((i.eq.139).and.(j.eq.45)) ) then !tropical islands
            vfc(i,j,12) = vfc(i,j,12) + vfc(i,j,15) + vfc(i,j,16)
            laic(i,j,12) = (laic(i,j,15)*vfc(i,j,15) 
     &           + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
            do m=1,12
               laim(m,i,j,12) = ( laim(m,i,j,15)*vfm(m,i,j,15) +
     &              laim(m,i,j,16)*vfm(m,i,j,16)) / 
     &              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
            enddo
            write(*,*) 'Replaced island crops.'
            naturalvegfound=1
         endif
      else
         write(*,*) 'STOP: Check array indices for grid res for nocrops'
         STOP
      endif

      if (naturalvegfound.eq.1) then
         vfc(i,j,15) = 0.
         vfc(i,j,16) = 0.
         vfm(:,i,j,15) = 0.
         vfm(:,i,j,16) = 0.
      endif

      end subroutine replace_crops


      subroutine fill_crops(IMn,JMn,vfc15,laic15
     i     ,vfm15,laim15,hm15,hsd15
     o     ,laiccrop, laimcrop, hmcrop, hsdcrop)
      !This performs in-fill once for herb crop (PFT15) LAI to create
      !an extended crop LAI data set for use with historical crop cover.
      integer :: IMn, JMn
      real*4 :: vfc15(:,:), laic15(:,:)  !(i,j)
      real*4 :: vfm15(:,:,:), laim15(:,:,:) !(m,i,j)
      real*4 :: hm15(:,:), hsd15(:,:) !(i,j)
      real*4 :: laiccrop(:,:), laimcrop(:,:,:) !OUTPUT
      real*4 :: hmcrop(:,:), hsdcrop(:,:)
      !--Local----
      integer :: i, j, m, ii, jj, nocropcells
      real*4 :: covsum15, laiavg15, covmsum15(12), laimavg15(12)
      real*4 :: hmavg15, hsdavg15
      character*80 :: titlefoo

      titlefoo = 'vfc15 in - crop fill'
      write(93) titlefoo, vfc15
      titlefoo = 'laic15 in - crop fill'
      write(93) titlefoo, laic15

      nocropcells=0
      do i=1,IMn
         do j=1,JMn
            if (vfc15(i,j).gt.0.d0) then !crops in cell, replicate
               laiccrop(i,j) = laic15(i,j)
               laimcrop(:,i,j) = laim15(:,i,j)
               hmcrop(i,j) = hm15(i,j)
               hsdcrop(i,j) = hsd15(i,j)
            else !(vfc15(i,j).eq.0.d0) then !no crops in cell, fill in LAI15
               covsum15 = 0.d0
               laiavg15 = 0.d0
               hmavg15 = 0.d0
               hsdavg15 = 0.d0
               laimavg15(:) = 0.d0
               covmsum15(:) = 0.d0
               do ii=i-1,i+1
                  do jj=j-1,j+1
                     if ( (ii.ge.1).and.(ii.le.IMn)
     &                    .and.(jj.ge.1).and.(jj.le.JMn) !in grid range
     &                    .and.((ii.ne.i).or.(jj.ne.j)) ) !not in i,j 
     &                    then
                        if (vfc15(ii,jj).gt.0.d0) then
                           covsum15 = covsum15 + vfc15(ii,jj)
                           laiavg15 = laiavg15 
     &                          + laic15(ii,jj)*vfc15(ii,jj)
                           hmavg15 = hmavg15 + hm15(ii,jj)*vfc15(ii,jj)
                           hsdavg15 = hsdavg15 
     &                          + hsd15(ii,jj)*vfc15(ii,jj)
                           covmsum15(:) = covmsum15(:) + vfm15(:,ii,jj)
                           laimavg15(:) = laimavg15(:) +
     &                          laim15(:,ii,jj)*vfm15(:,ii,jj)
                        endif
                     endif
                  enddo
               enddo
               if (covsum15>0.d0) then !adjacent crops found, assign laiavg
                  laiavg15 = laiavg15/covsum15
                  hmavg15 = hmavg15/covsum15
                  hsdavg15 = hsdavg15/covsum15
                  laimavg15(:) = laimavg15(:)/covmsum15(:)
                  laiccrop(i,j) = laiavg15
                  hmcrop(i,j) = hmavg15
                  hsdcrop(i,j) = hsdavg15
                  laimcrop(:,i,j) = laimavg15(:)
               else
                  nocropcells = nocropcells + 1
               endif
            endif
         enddo
      enddo

      write(*,*) 'Crop fill-in: no-crop cells:',nocropcells

      end subroutine fill_crops

      end module conversions

!------------------------------------------------------------------------
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
      
!------------------------------------------------------------------------

      program convert
      use conversions
      implicit none
      
      include 'netcdf.inc'


      !***************************************************
      !*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
      !***************************************************
      character*14, parameter :: EntPFT_fields(20) =
     &     (/
     &     'water         ',
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

      
      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES - short names   *
      !***************************************************
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


      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES                 *
      !***************************************************
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
      

      real*4, parameter :: undef = -1.e30

      ! define input file resolution
      !* .5x.5
      !integer, parameter :: IM=720, JM=360, KM=20
      !character*(*), parameter :: res_in="05x05"
      !character*(*), parameter :: res_in_int="720x360"
      !* 1x1
      !integer, parameter :: IM=360, JM=180, KM=20
      !character*(*), parameter :: res_in="1x1"
      !character*(*), parameter :: res_in_int="360x180"
      !* 2.5x2
!      integer, parameter :: IM=1440, JM=720, KM=20
      integer, parameter :: IM=144, JM=90, KM=20
!      character*(*), parameter :: res_in="2.5x2"
!      character*(*), parameter :: res_in="144x90"
!      character*(*), parameter :: res_in_int="144x90"
      character*(*), parameter :: res_in="144x90"
      character*(*), parameter :: res_in_int="144x90"

      integer :: fileid,dimidx,dimidy,dimidz,dd(4),varidx,varidy,varidz
      integer :: startA(1),startB(2),countA(1),countB(2),lenx,leny

      integer :: start2d(2),count2d(2),start3d(4),count3d(4)
      real*4 :: lcin(IM,JM),laiin(IM,JM),hin(IM,JM),hstd(IM,JM)
      integer :: fileidlclaim,fileidlclaimax
      integer :: varidlclaim,varidlclaimax
      character*80 :: fileheight,filestd,filelclaim,filelclaimax
      character*80 :: filetest
      integer :: fileidheight,varidheight
      integer :: fileidstd,varidstd,err
      character*60 :: PathFilepre, PathFilepost

      real*4 :: WTIN(144,90)
      real*4 :: OFFIA, DLATA, OFFIB, DLATB, DATMIS

      character*(*), parameter :: old_veg_file=
     &     "../../data/"//
     &     "V144X90_no_crops.ext"
      !* 4x5
!      integer, parameter :: IM=72, JM=46, KM=20
!      character*(*), parameter :: res_in="4x5"
!      character*(*), parameter :: res_in_int="72x46"
!      character*(*), parameter :: old_veg_file=
!     &     "/Users/nkiang/NancyResearch/GISS/Models/Ent/Code/cmrun/"//
!     &     "V72x46.1.cor2_no_crops.ext"

      !Igor's run:   
      !character*(*), parameter :: old_veg_file="V72x46_no_crops.ext"
      !character*(*), parameter :: old_veg_file="V144X90_no_crops.ext"


      integer, parameter :: N_HT = 19 ! number of height layers input

!      real*4 vf(IM,JM,KM), lai(IM,JM,KM)
      character*80 :: title(KM), title_tmp,title12(12,KM),titlehn(2,KM)

      ! Input values
      ! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
      integer, parameter :: IMn=IM, JMn=JM
      character*(*), parameter :: res_out=res_in
      character*(*), parameter :: res_out_int=res_in_int   

      character*(*), parameter :: file_checksum = 
     &     "../lc_lai_ent16/EntMM16_checksum_"//res_out_int//".ij"

      ! Input values, max, monthly
      real*4 vfn(IMn,JMn,KM), lain(IMn,JMn,KM), area(IMn,JMn), a
      real*4 vfnm(12,IMn,JMn,KM),lainm(12,IMn,JMn,KM),aream(IMn,JMn),am
      real*4 hmn(IMn,JMn,KM),hsdn(IMn,JMn,KM)
      ! Converted values
      real*4 vfc(IMn,JMn,KM), laic(IMn,JMn,KM), dvf, s
      character*80 :: titlec(18)
      ! Converted values - monthly
      real*4 vfm(12,IMn,JMn,KM), laim(12,IMn,JMn,KM)
      real*4 laimnc(IMn,JMn,KM), vfmnc(IMn,JMn,KM),laicropnc(IMn,JMn,12)

      character*80 :: titlem(12,18)
      ! Converted values - heights
      real*4 vfh(IMn,JMn,KM),hm(IMn,JMn,KM), hsd(IMn,JMn,KM)
      real*4 hmnc(IMn,JMn), hsdnc(IMn,JMn)
      real*4 vfcnc(IMn,JMn),laicnc(IMn,JMn)
      character*80 :: titleh(2,18) !1-h, 2-hsd
      ! Converted values - crop ext 
      real*4 laicropext(IMn,JMn), laimcropext(12,IMn,JMn)
      real*4 hmcropext(IMn,JMn),hsdcropext(IMn,JMn)
      ! Vars for calculating nocrops
      integer naturalfound, flag, nonaturalcount !if no natural veg in vicinity

      real*4 vf_xx(IMn,JMn), lai_xx(IMn,JMn)
      real*4 vf_yy(IMn,JMn), lai_yy(IMn,JMn)
      real*4 LAYER(IMn,JMn)
      character*80 :: title_xx="xx"
      character*80 :: title_yy="yy"
      character*80 :: titlefoo

      ! bs_brightratio = bare soil brightratio
      real*4 :: bs_brightratio(IMn,JMn), vfc_tmp(IMn,JMn)
      real*4 :: bs_brightratiopre(IMn,JMn)
      real*4 :: bs_brightratio_old(144,90)

      integer i, j, k, io, in, jn, maxpft, kx, m
      real*8 lat
      real*4 foo(IMn,JMn)
      integer, parameter :: delta_j = JM/JMn
      integer, parameter :: delta_i2 = (IM*2)/IMn
      real*8,  parameter :: delta_lat = 180/JM
      integer N_VEG             ! number of PFTs in output
      integer N_BARE            ! index of bare soil in output
      integer count

      integer :: ncidin,ncidout,varid,status,varidn(18)
      integer :: varids(18)
      character*20 :: inqvarin, inqvarout
      real*4 :: long(IM),lati(JM)
      character*256 :: fileoutnc
      real*4, dimension(12) :: time
      time = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
            
      !define Lon and Lat
      call calc_lon_lat_144x90(IM,JM,long,lati)
!      call calc_lon_lat_025x025(IM,JM,long,lati)

      ! lc & laimax file of Ent 17 PFTs
!     file lc_laimax
      PathFilepre= '../lc_lai_ent/'
      PathFilepost= 'EntMM_lc_laimax_'//res_in//'.nc'
      filelclaimax = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filelclaimax,NF_WRITE,fileidlclaimax)
      write(*,*) err

      ! get lc
      do k=1,20
         inqvarin = trim(EntPFT_fields(k))//'_lc'
         err = NF_INQ_VARID(fileidlclaimax,inqvarin,varidlclaimax)
         write(*,*) err
	 start2d(1) = 1
         start2d(2) = 1
         count2d(1) = IM
         count2d(2) = JM
         err = NF_GET_VARA_REAL(fileidlclaimax,varidlclaimax,start2d,
     &           count2d,lcin)
	 write(*,*) err, 'lc input'
         vfn(:,:,k)=lcin
      enddo

      ! get laimax
      do k=21,40
         inqvarin = trim(EntPFT_fields(k-20))//'_lai'
         err = NF_INQ_VARID(fileidlclaimax,inqvarin,varidlclaimax)
         write(*,*) err
         start2d(1) = 1
         start2d(2) = 1
         count2d(1) = IM
         count2d(2) = JM
         err = NF_GET_VARA_REAL(fileidlclaimax,varidlclaimax,start2d,
     &           count2d,laiin)
         write(*,*) err, 'laimax input'
         lain(:,:,k-20)=laiin
      enddo

      !Check if mismatch lc or lai values (one is zero and the other not)
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)

      ! monthly lc and lai files
      do m=1,12
         PathFilepre= '../lc_lai_ent/'
         PathFilepost= 'EntMM_lc_lai_'//trim(MONTH(m))//'_'
     &       // res_in//'.nc'
         filelclaim = trim(PathFilepre)//trim(PathFilepost)
         err = NF_OPEN(filelclaim,NF_WRITE,fileidlclaim)
	 write(*,*) err, 'lc monthly files'
         do k=1,20
            inqvarin = trim(EntPFT_fields(k))//'_lc'
            err = NF_INQ_VARID(fileidlclaim,inqvarin,varidlclaim)
            write(*,*) err
            start2d(1) = 1
            start2d(2) = 1
            count2d(1) = IM
            count2d(2) = JM
            err = NF_GET_VARA_REAL(fileidlclaim,varidlclaim,start2d,
     &           count2d,lcin)
	    write(*,*) err, 'lc monthly fields'
            vfnm(m,:,:,k)=lcin
         enddo

         do k=21,40
            inqvarin = trim(EntPFT_fields(k-20))//'_lai'
            err = NF_INQ_VARID(fileidlclaim,inqvarin,varidlclaim)
            write(*,*) err, 'lai monthly files'
            start2d(1) = 1
            start2d(2) = 1
            count2d(1) = IM
            count2d(2) = JM
            err = NF_GET_VARA_REAL(fileidlclaim,varidlclaim,start2d,
     &           count2d,laiin)
            write(*,*) err
            lainm(m,:,:,k-20)=laiin
         enddo
      enddo
      

      ! height file - insert dummy WATER layer at beginning toi avoid
      !               confusion in numbering with vfc and vfm.
      PathFilepre= '../../data/height/'
      PathFilepost= 'EntGVSDmosaic17_height_'//res_in//'_lai3g.nc'
      fileheight = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(fileheight,NF_WRITE,fileidheight)
      inqvarin = 'SimardHeights'
      err = NF_INQ_VARID(fileidheight,inqvarin,varidheight)
      hmn(:,:,:) = 0.
      do k=1,N_HT !height
         start3d(1) = k
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM
         count3d(3) = JM
         err = NF_GET_VARA_REAL(fileidheight,varidheight,start3d,
     &        count3d,hin)
         write(*,*) err, 'height'
         hmn(:,:,k+1)=hin
      enddo

      PathFilepre= '../../data/height/'
      PathFilepost= 'EntGVSDmosaic17_heightstd_'//res_in//'.nc'
      filestd = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filestd,NF_WRITE,fileidstd)
      write(*,*) err
      inqvarin = 'SimardHeights'
      err = NF_INQ_VARID(fileidstd,inqvarin,varidstd)
      write(*,*) err
      hsdn(:,:,:) = 0.
      do k=1,N_HT !height std
         start3d(1) = k
         start3d(2) = 1
         start3d(3) = 1
         count3d(1) = 1
         count3d(2) = IM
         count3d(3) = JM
         err = NF_GET_VARA_REAL(fileidstd,varidstd,start3d,
     &        count3d,hstd)
         write(*,*) err, 'height std'
         hsdn(:,:,k+1)=hstd
      enddo

      open(1,file='../../'//
     &     'data/height/'//
     &      'EntGVSDmosaic17_height_144x90.ij',
     &     form='unformatted',status="old")
      titlehn(1,1) = 'NO WATER LAYER'
      do k=1,N_HT !height
         read(1) titlehn(1,k+1)
         write(*,*) titlehn(1,k+1)
      enddo
      titlehn(2,1) = 'NO WATER LAYER'
      do k=1,N_HT !stdev
         read(1) titlehn(2,k+1)
         write(*,*) titlehn(2,k+1)
      enddo

!      call get_bare_soil_brightratio(IMn,JMn, old_veg_file
!     &     ,bs_brightratio)
      call get_bare_soil_brightratio(144,90, old_veg_file
     &     ,bs_brightratio_old)

      DLATA = (180./90) * 60.
      OFFIA = 0.0
      DLATB = (180./JMn) * 60.
      OFFIB = 0.0
      DATMIS = undef
      
      !* Setup grids
      write(*,*) 144,90,OFFIA,DLATA,
     &     IM,JM,OFFIB,DLATB,DATMIS
      
      write(*,*) 'Calling HNTR40'
      call HNTR40(144,90,OFFIA,DLATA,
     *     IM,JM,OFFIB,DLATB, DATMIS)
      write(*,*) 'Finished HNTR40'
      
      do i=1,144
         do j=1,90
            !* Weights, if any
            WTIN(i,j) = 1.
            if (bs_brightratio_old(i,j).eq.undef) WTIN(i,j) = 0.
         enddo
      enddo
      
      write(*,*) 'Calling HNTR40'
      call HNTR4P(WTIN, bs_brightratio_old,bs_brightratio)
      write(*,*) 'Finished HNTR40'

!      do i=1,IMn
!           do j=1,JMn
!	      bs_brightratio(i,j) = bs_brightratiopre(i,JMn+1-j) ! flip left to right
!	   enddo
!      enddo

      err = NF_CREATE('../lc_lai_ent/bs_brightratio.nc',
     &     NF_CLOBBER,ncidout)
      write(*,*) err, ncidout
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      write(*,*) err, 'dimidx'
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      write(*,*) err, 'dimidy'
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      err=NF_DEF_VAR(ncidout,'bs_brightratio',NF_REAL,2,dd,varid)
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      startB(1)=1
      startB(2)=1
      countB(1)=IM
      countB(2)=JM
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,bs_brightratio)
      write(*,*) err, 'put'
      err = NF_CLOSE(ncidout)
      
      open(100,file=file_checksum,form="unformatted",status="unknown")


      ! lc & laimax file of Ent 17 PFTs: JUST FOR TITLE(K)
      open(1,file="../lc_lai_ent/EntMM_lc_laimax_144x90.bin",
     &     form="unformatted",status="old")
      do k=1,KM 
        read(1) title(k)
        write(*,*) title(k)
      enddo
      close(1)

      do m=1,12
         open(1,file="../lc_lai_ent/EntMM_lc_laimax_144x90.bin",
     &     form="unformatted",status="old")
         do k=1,KM
            read(1) title12(m,k)
            write(*,*) title12(m,k), KM
         enddo
         close(1)
      enddo

      
      !* Convert to GISS 16 pfts format

      ! first 14 pfts though grass are the same, ignore WATER
      !  lc laimax
      vfc(:,:,1:14) = vfn(:,:,2:15)
      laic(:,:,1:14) = lain(:,:,2:15)
      titlec(1:14) = title(2:15)
      !  lc lai monthly
      vfm(:,:,:,1:14) = vfnm(:,:,:,2:15)
      laim(:,:,:,1:14) = lainm(:,:,:,2:15)
      do m=1,12
         titlem(m,1:14) = title12(m,2:15)
      enddo
      !  heights
      vfh(:,:,1:14) = vfn(:,:,2:15) !Should be the same cover from MODIS.
      hm(:,:,1:14) = hmn(:,:,2:15) 
      hsd(:,:,1:14) = hsdn(:,:,2:15)
      titleh(1,1:14) = titlehn(1,2:15)
      titleh(2,1:14) = titlehn(2,2:15)

      ! crops
#     ifdef COMBINE_CROPS_C3_C4
      do j=1,JMn
         do i=1,IMn
            !lc laimax
            a = vfn(i,j,16) + vfn(i,j,17)
            if ( a > 0. ) then
               laic(i,j,15) = (vfn(i,j,16)*lain(i,j,16)
     &              + vfn(i,j,17)*lain(i,j,17)) / a
               vfc(i,j,15) = a
            else
               laic(i,j,15) = 0.
               vfc(i,j,15) = 0.
            endif
            !lc lai monthly
            do m=1,12
               am = vfnm(m,i,j,16) + vfnm(m,i,j,17)
               if (am > 0. ) then
                  laim(m,i,j,15) = (vfnm(m,i,j,16)*lainm(m,i,j,16)
     &                 + vfnm(m,i,j,17)*lainm(m,i,j,17)) / am
                  vfm(m,i,j,15) = am
               else
                  laim(m,i,j,15) = 0.
                  vfm(m,i,j,15) = 0.
               endif
            enddo
            !heights - DO NOT AVERAGE. PRESERVE HEIGHTS. LAI will scale density
            a = vfn(i,j,16) + vfn(i,j,17)  !input cover
            if ( a > 0. ) then
!               hm(i,j,15) = (vfn(i,j,16)*hmn(i,j,16)
!     &              + vfn(i,j,17)*hmn(i,j,17)) / a
               if ((hmn(i,j,16)>0.).and.(hmn(i,j,17)>0.)) then
                  !average if both exist
                  hm(i,j,15) = (vfn(i,j,16)*hmn(i,j,16)
     &                 + vfn(i,j,17)*hmn(i,j,17)) / a
               else
                  !don't average if only one or none exists
                  hm(i,j,15) = max(hmn(i,j,16),hmn(i,j,17))
               endif
               !Sum of squares for sd.  Don't weight if only one or less exists
               if ((hmn(i,j,16)>0.).and.(hmn(i,j,17)>0.)) then
                  hsd(i,j,15) = sqrt((vfn(i,j,16)*hsdn(i,j,16)**2
     &                 + vfn(i,j,17)*hsdn(i,j,17)**2) / a)
               else
                  hsd(i,j,15) = max(hsd(i,j,16),hsd(i,j,17))
               endif
               vfh(i,j,15) = a
            else
               hm(i,j,15) = 0.
               hsd(i,j,15) = 0.
               vfh(i,j,15) = 0.
            endif
         enddo                  !j=
      enddo                     !i=
      write(*,*) "Re-doing crops.."
      titlec(15) = title(16)
      titlec(15)(1:18) = "15 - crops herb   "
      do m=1,12
         titlem(m,15) = title12(m,16)
         titlem(m,15)(1:18) = "15 - crops herb   "
      enddo
      titleh(1,15) = titlehn(1,16)
      titleh(2,15) = titlehn(2,16)
      titleh(1,15)(1:18) = "15 - crops herb   "
      titleh(2,15)(1:18) = "15 - crops herb   "

      write(*,*) titlec(15)
      ! crops woody
      vfc(:,:,16) = vfn(:,:,18)
      laic(:,:,16) = lain(:,:,18)
      titlec(16) = '16 - '//title(18)(6:80)
      write(*,*) "titlec: "
      write(*,*) titlec(16)
      do m=1,12
         vfm(m,:,:,16) = vfnm(m,:,:,18)
         laim(m,:,:,16) = lainm(m,:,:,18)
         titlem(m,16) = '16 - '//title12(m,18)(6:80)
      enddo
      vfh(:,:,16) = vfn(:,:,18)
      titleh(1,16) = '16 - '//titlehn(1,18)(6:80)
      titleh(2,16) = '16 - '//titlehn(2,18)(6:80)
      ! bare soil
      vfc(:,:,17) = vfn(:,:,20)
      laic(:,:,17) = lain(:,:,20)
      titlec(17) = title(20)
      vfm(:,:,:,17) = vfnm(:,:,:,20)
      laim(:,:,:,17) = lainm(:,:,:,20)
      do m=1,12
         titlem(m,17) = title12(m,20)
      enddo
      vfh(:,:,17) = vfn(:,:,20)
      hm(:,:,17) = hmn(:,:,20)
      hsd(:,:,17) = hsdn(:,:,20)
      titleh(1,17) = titlehn(1,20)
      titleh(2,17) = titlehn(2,20)
      N_VEG = 16
      N_BARE = 17

      titlefoo = 'vfn16'
      write(92) titlefoo, vfn(:,:,16)
      titlefoo = 'vfn17'
      write(92) titlefoo, vfn(:,:,17)
      titlefoo = 'Crops 15 after combining C3 and C4'
      write(92) titlefoo, vfc(:,:,15)

#     else
      !crops
      vfc(:,:,15:17) = vfn(:,:,16:18)
      laic(:,:,15:17) = lain(:,:,16:18)

      titlec(15:17) = title(16:18)
      do m=1,12
         vfm(m,:,:,15:17) = vfnm(m,:,:,16:18)
         laim(m,:,:,15:17) = lainm(m,:,:,16:18)
         titlem(m,15:17) = title12(m,16:18)
      enddo
      vfh(:,:,15:17) = vfn(:,:,16:18)
      hm(:,:,15:17) = hmn(:,:,16:18)
      titleh(1,15:17) = titlehn(1,16:18)
      titleh(2,15:17) = titlehn(2,16:18)
      ! bare soil
      vfc(:,:,18) = vfn(:,:,20)
      laic(:,:,18) = lain(:,:,20)
      titlec(18) = title(20)
      vfm(:,:,:,18) = vfnm(:,:,:,20)
      laim(:,:,:,18) = lainm(:,:,:,20)
      titlem(:,18) = title12(:,20)
      vfh(:,:,18) = vfh(:,:,20)
      hm(:,:,18) = hmn(:,:,20)
      hsd(:,:,18) = hsd(:,:,20)
      titleh(1,18) = titlehn(1,20)
      titleh(2,18) = titlehn(2,20)
      N_VEG = 17
      N_BARE = 18
#     endif

      do k=1,N_BARE
        write(3) titlec(k), vfc(:,:,k)
      enddo
      do k=1,N_BARE
        write(3) titlec(k), laic(:,:,k)
      enddo

      ! check if "bare" soil is not bare
      vf_xx(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .01 .and. laic(i,j,N_BARE) > .5 ) then
            vf_xx(i,j) = vfc(i,j,N_BARE)
            lai_xx(i,j) = laic(i,j,N_BARE)
          endif
        enddo
      enddo
      write(4) title_xx, vf_xx(:,:)
      write(4) title_xx, lai_xx(:,:)

      vf_yy(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,10) > .1 .and. laic(i,j,10) < .5 ) then
            vf_yy(i,j) = vfc(i,j,10)
            lai_yy(i,j) = laic(i,j,10)
          endif
        enddo
      enddo
      write(4) title_yy, vf_yy(:,:)
      write(4) title_yy, lai_yy(:,:)

      vf_yy(:,:) = 0.
      lai_yy(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .1 .and. laic(i,j,10) < .01 
     &      .and. laic(i,j,9) < .01 .and. laic(i,j,11) < .01 
     &      .and. laic(i,j,12) < .01 .and. laic(i,j,13) < .01 ) then
            vf_yy(i,j) = vfc(i,j,N_BARE)
            lai_yy(i,j) = laic(i,j,N_BARE)
          endif
        enddo
      enddo

      write(4) title_yy, vf_yy(:,:)
      write(4) title_yy, lai_yy(:,:)


      !!!! do conversions !!!!

      ! convert sparse veg to cold adapted shrub 9 if present
      do j=1,JMn
        do i=1,IMn
          s = sum(vfc(i,j,1:N_BARE))
          if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
             write(*,*) 'ERROR orig:  max and monthly lc different'
     &            ,i,j ,s,sum(vfm(1,i,j,1:N_BARE))
             write(*,*) vfc(i,j,1:N_BARE)
             write(*,*) vfm(1,i,j,1:N_BARE)
          endif
          if( vfc(i,j,N_BARE) > .0 .and. vfc(i,j,N_BARE) < .15
     &         .and. laic(i,j,N_BARE) > .0
     &         .and. vfc(i,j,9) > .0 ) then
            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,9), laic(i,j,9), laic(i,j,9) )
                                ! lai >= lai(9)
            do m=1,12
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &              vfm(m,i,j,9),vfm(m,i,j,9), vfc(i,j,9))
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,9),hm(i,j,9),hsd(i,j,9), vfc(i,j,9))
          endif

          s = sum(vfc(i,j,1:N_BARE))
          if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
             write(*,*) 'ERROR sparse:  max and monthly lc different'
     &            ,i,j, s,sum(vfm(1,i,j,1:N_BARE))
             write(*,*) vfc(i,j,1:N_BARE)
             write(*,*) vfm(1,i,j,1:N_BARE)
          endif
        enddo
      enddo

      ! convert sparse veg to arid adapted shrub 10 if present
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .0 .and. vfc(i,j,N_BARE) < .15
     &         .and. laic(i,j,N_BARE) > .0
     &         .and. vfc(i,j,10) > .0 ) then

            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,10), laic(i,j,10), laic(i,j,10) )
                                ! lai >= lai(10)
            do m=1,12
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &              vfm(m,i,j,10),laim(m,i,j,10),vfc(i,j,10))
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,10),hm(i,j,10),hsd(i,j,10), vfc(i,j,10))
          endif
        enddo
      enddo
      ! convert the rest of sparse veg to crop 15 if present
      do j=1,JMn
        do i=1,IMn
           if( vfc(i,j,N_BARE) > .0 .and. laic(i,j,N_BARE) > .0
     &         .and. vfc(i,j,15) > .0 ) then
!             print *, 'Converting spare to crop/bare',i,j,
!     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), vfc(i,j,15)
             call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &            vfc(i,j,15), laic(i,j,15), laic(i,j,15))
!             print *, 'After conversion:            ',i,j,
!     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), 
!     &            vfc(i,j,15), laic(i,j,15)
             do m=1,12
                call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &               vfm(m,i,j,15),laim(m,i,j,15),vfc(i,j,15))
             enddo
             call convert_vfh(
     &            vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &            vfh(i,j,15),hm(i,j,15),hsd(i,j,15), vfc(i,j,15))
          endif
        enddo
      enddo

      ! convert the rest of sparse veg to pft with biggest fraction
      ! (if present)
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .0 .and. laic(i,j,N_BARE) > .0 ) then

            maxpft = maxloc( vfc(i,j,1:16), 1 )
            print *, "max pft is ",maxpft
            if ( vfc(i,j,maxpft) < .0001 ) cycle
            
            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,maxpft), laic(i,j,maxpft), laic(i,j,maxpft))
            do m=1,12
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE)
     &           ,vfm(m,i,j,maxpft),laim(m,i,j,maxpft),vfc(i,j,maxpft))
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,maxpft),hm(i,j,maxpft),hsd(i,j,maxpft), 
     &           vfc(i,j,maxpft))

          endif
        enddo
      enddo
      ! convert the rest of sparse veg to arid adapted shrub 10
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .0 .and. laic(i,j,N_BARE) > .0 ) then

            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,10), laic(i,j,10), .0 )
            do m=1,12
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &              vfm(m,i,j,10),laim(m,i,j,10),vfc(i,j,10))
            enddo

            if (vfc(i,j,10) > 0.) then
               hm(i,j,10) = 2.0 !Check simard.f Set_shrub_height for value!
               hsd(i,j,10) = 0.
            endif
            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,10),hm(i,j,10),hsd(i,j,10), vfc(i,j,10))

          endif
        enddo
      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif

      ! check titles
      write(*,*) 'titlec:'
      do k=1,N_BARE
         write(*,*) trim(titlec(k))
      enddo
      write(*,*) 'titlem:'
      do m=1,12
         write(*,*) MONTH(m)
         do k=1,N_BARE
            write(*,*) trim(titlem(m,k))
         enddo
      enddo

!*************************************************************************
      call write_output(titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_pure_"//res_out//".ij"
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &     ,"max_pure","  ",res_out)
      write(*,*) "pure"

      ! save netcdf lai_max_pure
      fileoutnc =
     & '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_max_pure.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               laicnc(i,j) = laic(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,laicnc)
      enddo
      err = NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
      ! save netcdf lc_max_pure
      fileoutnc =
     & '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lc_max_pure.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               vfcnc(i,j) = vfc(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,vfcnc)
      enddo
      err = NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
!*************************************************************************

      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)
      do m=1,12
         call write_output(titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE,
     &        "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &        ,"pure",MONTH(m), res_out)
      enddo
      
      ! save netcdf lai
      do m=1,12
         
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lai_pure_'//MONTH(m)//'.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         do k=1,18
            inqvarout = EntPFT_shorttitle(k)
            err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         enddo
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do k=1,18
            do i=1,IMn
               do j=1,JMn
                  laicnc(i,j) = laim(m,i,j,k)
               enddo
            enddo
            err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &           EntPFT_title(k))
            startB(1)=1
            startB(2)=1
            countB(1)=IM
            countB(2)=JM
            err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &           laicnc)
            write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         enddo
         err = NF_CLOSE(ncidout)
         write(*,*) err, 'nf_close out ',trim(fileoutnc)
      enddo
      
      ! save netcdf lc
      do m=1,12
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lc_pure_'//MONTH(m)//'.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         do k=1,18
            inqvarout = EntPFT_shorttitle(k)
            err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         enddo
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do k=1,18
            do i=1,IMn
               do j=1,JMn
                  vfcnc(i,j) = vfm(m,i,j,k)
               enddo
            enddo
            err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &           EntPFT_title(k))
            startB(1)=1
            startB(2)=1
            countB(1)=IM
            countB(2)=JM
            err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &           vfcnc)
         enddo
         err = NF_CLOSE(ncidout)
         write(*,*) err, 'nf_close out ',trim(fileoutnc)
      enddo
      
!*************************************************************************
      call write_output_h(titleh, hm, hsd, N_BARE,
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16_height_pure"
     &     //".ij","   ",res_out)
      
      ! save netcdf height_pure
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_height_pure.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'hgt_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         inqvarout = 'stdev_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varids(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               hmnc(i,j) = hm(i,j,k)
               hsdnc(i,j) = hsd(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         err=NF_PUT_ATT_TEXT(ncidout,varids(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,hmnc)
         err=NF_PUT_VARA_REAL(ncidout,varids(k),startB,countB,hsdnc)
      enddo
      err = NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)

!*************************************************************************
      
      foo(:,:) = 0.
      do j=1,JMn
         do i=1,IMn
            do k=1,KM
               foo(i,j) = foo(i,j) + vfc(i,j,k)
            enddo
         enddo
      enddo
      titlefoo = "LC pure checksum"
      write(100) titlefoo, foo
      
      call sum_lai_single(100,IMn,JMn,KM,N_BARE,"LAI MAX", vfc,laic)
      call sum_lai(100,IMn,JMn,KM,N_BARE,"LAI PURE", vfm,laim)
      
      ! convert arid adapted shrub with lai < .15 to bare soil
      ! and restrict lai >= .15 
      do j=1,JMn
         do i=1,IMn
            if( vfc(i,j,10) > .0 .and. laic(i,j,10) < .15 ) then
               
               if ((i.eq.70).and.(j.eq.11)) then
                  print *,'vfc10 before',vfc(i,j,10),vfc(i,j,N_BARE),
     &                 laic(i,j,10),laic(i,j,N_BARE),i,j
!               STOP
               endif
               call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &              vfc(i,j,10), laic(i,j,10), .15 )
                                ! lai >= .15
               if (vfc(i,j,10).le.0.0) then
                  print *,'vfc10=',vfc(i,j,10),vfc(i,j,N_BARE),
     &                 laic(i,j,10),laic(i,j,N_BARE),i,j
!               STOP
               endif
               do m=1,12
                  if (vfc(i,j,10).gt.0.0) then  
                     call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,
     &                    N_BARE),
     &                    vfm(m,i,j,10),laim(m,i,j,10), vfc(i,j,10))
                  endif
               enddo
               call convert_vfh(
     &              vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &              vfh(i,j,10),hm(i,j,10),hsd(i,j,10), vfc(i,j,10))
               
            endif
            s = sum(vfc(i,j,1:N_BARE))
            if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
               write(*,*) 'ERROR trim:  max and monthly lc different'
     &              ,i,j, s,sum(vfm(1,i,j,1:N_BARE))
               write(*,*) vfc(i,j,1:N_BARE)
               write(*,*) vfm(1,i,j,1:N_BARE)
            endif
            
         enddo
      enddo
      
#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif
      
!*************************************************************************
      call write_output(titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_"//res_out//".ij"
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &     ,"max_trimmed","  ",res_out)
      write(*,*) "trimmed"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)
      
      ! save netcdf lai_max_trimmed
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_max_trimmed.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               laicnc(i,j) = laic(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,laicnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err = NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)

      
      ! save netcdf lc_max_trimmed
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lc_max_trimmed.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               vfcnc(i,j) = vfc(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,vfcnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err = NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)

!*************************************************************************

      do m=1,12
         call write_output(titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE,
     &        "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &        ,"trimmed",MONTH(m), res_out)
      enddo
      
      ! save netcdf lai_trimmed
      do m=1,12
          fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_trimmed_'//MONTH(m)//'.nc'
          err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
          write(*,*) err, 'nf_create out ',trim(fileoutnc)
          err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
          err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
          err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
          err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
          dd(1)=dimidx
          dd(2)=dimidy
          do k=1,18
             inqvarout = EntPFT_shorttitle(k)
             err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
          enddo
          err=NF_ENDDEF(ncidout)
          err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
          err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
          do k=1,18
              do i=1,IMn
                  do j=1,JMn
                      laicnc(i,j) = laim(m,i,j,k)
                  enddo
              enddo
             err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &            EntPFT_title(k))
             startB(1)=1
             startB(2)=1
             countB(1)=IM
             countB(2)=JM
             err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &            laicnc)
             write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
          enddo
          err=NF_CLOSE(ncidout)
          write(*,*) err, 'nf_close out ',trim(fileoutnc)
       enddo
       
      ! save netcdf lc_trimmed
       do m=1,12
          fileoutnc =
     &         '../lc_lai_ent16/nc/V'//res_out_int//
     &         '_EntMM16_lc_trimmed_'//MONTH(m)//'.nc'
          err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
          write(*,*) err, 'nf_create out ',trim(fileoutnc)
          err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
          err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
          err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
          err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
          dd(1)=dimidx
          dd(2)=dimidy
          do k=1,18
             inqvarout = EntPFT_shorttitle(k)
             err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
          enddo
          err=NF_ENDDEF(ncidout)
          err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
          err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
          do k=1,18
             do i=1,IMn
                do j=1,JMn
                   vfcnc(i,j) = vfm(m,i,j,k)
                enddo
             enddo
             err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &            EntPFT_title(k))
             startB(1)=1
             startB(2)=1
             countB(1)=IM
             countB(2)=JM
             err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &            vfcnc)
          enddo
          err=NF_CLOSE(ncidout)
          write(*,*) err, 'nf_close out ',trim(fileoutnc)
       enddo
!*************************************************************************

       foo(:,:) = 0.
       do j=1,JMn
          do i=1,IMn
             do k=1,KM
                foo(i,j) = foo(i,j) + vfc(i,j,k)
             enddo
        enddo
      enddo
      titlefoo = "LC trimmed checksum"
      write(100) titlefoo, foo
      call sum_lai(100,IMn,JMn,KM,N_BARE,"LAI trimmed", vfm, laim)
      
      
      ! rescale fractions so that they sum to 1
      write(*,*) 'Rescaling...'
      do j=1,JMn
         do i=1,IMn
            s = sum( vfc(i,j,1:N_BARE) )
            if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
               write(*,*) 'ERROR scale0:  max and monthly lc different'
     &              ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
               write(*,*) 'vfc',i,j,vfc(i,j,1:N_BARE)
               write(*,*) 'vfm',i,j,vfm(1,i,j,1:N_BARE)
            endif
            if ( s > 1.00001 ) then
!            if ( abs(s-1.0) > 0.00001 ) then
               vfc(i,j,1:N_BARE) = vfc(i,j,1:N_BARE) / s
               vfm(:,i,j,1:N_BARE) = vfm(:,i,j,1:N_BARE) / s
              !vfh(i,j,1:N_BARE) = vfh(i,j,1:N_BARE) / s !Heights are vertical, no rescale
            endif
            s = sum( vfc(i,j,1:N_BARE) ) 
            if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
               write(*,*) 'ERROR scale1:  max and monthly lc different'
     &              ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
               write(*,*) 'vfc',i,j,vfc(i,j,1:N_BARE)
               write(*,*) 'vfm',i,j,vfm(1,i,j,1:N_BARE)
            endif
         enddo
      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif

!*************************************************************************
      call write_output(titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_scaled_"
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &     ,"max_trimmed_scaled","  ",res_out)
      write(*,*) "trimmed scaled"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)

	  ! save netcdf lai_max_trimmed_scaled
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_max_trimmed_scaled.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               laicnc(i,j) = laic(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,laicnc)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
      ! save netcdf lc_max_trimmed_scaled
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lc_max_trimmed_scaled.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               vfcnc(i,j) = vfc(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,vfcnc)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
!*************************************************************************
       do m=1,12
          call write_output(titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_lai_trimmed_scaled_"
!     &        //MONTH(m)//"_"//res_out//".ij"
     &         "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &         ,"trimmed_scaled",MONTH(m), res_out)
       enddo
       
      ! save netcdf 12 months x 18 PFTs lai trimmed and crops
       do m = 1,12
          fileoutnc =
     &         '../lc_lai_ent16/nc/V'//res_out_int//
     &         '_EntMM16_lai_trimmed_scaled_'//MONTH(m)//'.nc'
          err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
          write(*,*) err, 'nf_create out ',trim(fileoutnc)
          err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
          err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
          err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
          err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
          dd(1)=dimidx
          dd(2)=dimidy
          do k=1,18
             inqvarout = EntPFT_shorttitle(k)
             err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
          enddo
          err=NF_ENDDEF(ncidout)
          err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
          err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
          do k=1,18
             do i=1,IMn
                do j=1,JMn
                   laicnc(i,j) = laim(m,i,j,k)
                enddo
             enddo
             err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &            EntPFT_title(k))
             startB(1)=1
             startB(2)=1
             countB(1)=IM
             countB(2)=JM
             err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &            laicnc)
             write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
          enddo
          err=NF_CLOSE(ncidout)
          write(*,*) err, 'nf_close out ',trim(fileoutnc)
       enddo

       ! lc
       do m = 1,12
          fileoutnc =
     &         '../lc_lai_ent16/nc/V'//res_out_int//
     &         '_EntMM16_lc_trimmed_scaled_'//MONTH(m)//'.nc'
          err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
          write(*,*) err, 'nf_create out ',trim(fileoutnc)
          err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
          err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
          err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
          err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
          dd(1)=dimidx
          dd(2)=dimidy
          do k=1,18
             inqvarout = EntPFT_shorttitle(k)
             err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
          enddo
          err=NF_ENDDEF(ncidout)
          err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
          err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
          do k=1,18
             do i=1,IMn
                do j=1,JMn
                   vfcnc(i,j) = vfm(m,i,j,k)
                enddo
             enddo
             err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &            EntPFT_title(k))
             startB(1)=1
             startB(2)=1
             countB(1)=IM
             countB(2)=JM
             err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &            vfcnc)
             write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
          enddo
          err=NF_CLOSE(ncidout)
          write(*,*) err, 'nf_close out ',trim(fileoutnc)
       enddo
       
!*************************************************************************

      call write_output_h(titleh, hm, hsd, N_BARE,
     &     "../lc_lai_ent16/V"//res_out_int//
     &     "_EntMM16_height_trimmed_scaled.ij","   ",res_out)

      foo(:,:) = 0.
      do j=1,JMn
         do i=1,IMn
            do k=1,KM
               foo(i,j) = foo(i,j) + vfc(i,j,k)
            enddo
         enddo
      enddo
      titlefoo = "LC trimmed scaled checksum"
      write(100) titlefoo, foo
      call sum_lai(100,IMn,JMn,KM,N_BARE,"LAI trimmed scaled"
     &     , vfm, laim)
      
      
    ! save netcdf simard heights trimmed all
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_height_trimmed_scaled.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'hgt_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         inqvarout = 'stdev_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varids(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               hmnc(i,j) = hm(i,j,k)
               hsdnc(i,j) = hsd(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,hmnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         err=NF_PUT_ATT_TEXT(ncidout,varids(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varids(k),startB,countB,hsdnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)


      !Generate fill-in crop cover from trimmed_scaled before doing nocrops
      !Herb crop only, since right now zero woody crops.
      laicropext(:,:) = 0.d0    !(i,j) only
      laimcropext(:,:,:) = 0.d0 !(m,i,j)
      hmcropext(:,:) = 0.d0     !(i,j)
      hsdcropext(:,:) = 0.d0    !(i,j)
      
!*************************************************************************
      call fill_crops(IMn,JMn,vfc(:,:,15),laic(:,:,15)
     i     ,vfm(:,:,:,15),laim(:,:,:,15),hm(:,:,15),hsd(:,:,15)
     o     ,laicropext(:,:),laimcropext(:,:,:)
     o     ,hmcropext,hsdcropext)

      titlefoo = "15 - max crops herb ext"
      call write_output_single(titlefoo, vfc(:,:,15), laicropext,
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &     ,"max_trimmed_scaled_crops.ext1","  ",res_out)
      write(*,*) "laimax trimmed scaled crops ext1"
      
	  ! save netcdf lai_max_trimmed_scaled
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_max_trimmed_scaled_crops.ext1.nc'    
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      k=15
      inqvarout = EntPFT_shorttitle(k)
      err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varid)
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do i=1,IMn
         do j=1,JMn
            laicnc(i,j) = laic(i,j,15)
         enddo
      enddo
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     EntPFT_title(k))
      startB(1)=1
      startB(2)=1
      countB(1)=IM
      countB(2)=JM
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,laicnc)
      write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      err=NF_CLOSE(ncidout)

      
      ! save netcdf lc_max_trimmed_scaled
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lc_max_trimmed_scaled_crops.ext1.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      k=15
      inqvarout = EntPFT_shorttitle(k)
      err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varid)
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do i=1,IMn
         do j=1,JMn
            vfcnc(i,j) = vfc(i,j,15)
         enddo
      enddo
      err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &     EntPFT_title(k))
      startB(1)=1
      startB(2)=1
      countB(1)=IM
      countB(2)=JM
      err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,vfcnc)
      write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      err=NF_CLOSE(ncidout)

!*************************************************************************
      do m=1,12
         titlefoo = "15 - crops herb ext"
         call write_output_single(titlefoo, vfm(m,:,:,15)
     &        , laimcropext(m,:,:)
     &        ,"../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &        ,"trimmed_scaled_crops.ext1",MONTH(m), res_out)
      enddo
      write(*,*) "lai month trimmed scaled crops ext1"

      ! save netcdf lai
      do m=1,12
         
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lai_trimmed_scaled_crops.ext1_'//MONTH(m)//
     &        '.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         k=15
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varid)
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do i=1,IMn
            do j=1,JMn
               laicnc(i,j) = laimcropext(m,i,j)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,laicnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         err=NF_CLOSE(ncidout)
      enddo

      ! save netcdf lc
      do m=1,12
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lc_trimmed_scaled_crops.ext1_'//MONTH(m)//'.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         k=15
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varid)
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do i=1,IMn
            do j=1,JMn
               vfcnc(i,j) = vfm(m,i,j,15)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varid,'long_name',40,
     &        EntPFT_title(k))
                  startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varid,startB,countB,vfcnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         err=NF_CLOSE(ncidout)
      enddo
      
!*************************************************************************
      titlefoo = "15 - crops herb ext"
      call write_output_h_single(titleh(:,15), hmcropext, hsdcropext,
     &     "../lc_lai_ent16/V"//res_out_int//
     &     "_EntMM16_height_trimmed_scaled_crops.ext1.ij",
     &     "   ",res_out)
      write(*,*) "height trimmed scaled crops ext1"
      
      ! save netcdf lai_max_trimmed_scaled
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_height_trimmed_scaled_crops.ext1.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      k=15
      inqvarout = 'hgt_'//EntPFT_shorttitle(k)
      err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(1))
      inqvarout = 'stdev_'//EntPFT_shorttitle(k)
      err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(2))
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do i=1,IMn
         do j=1,JMn
            hmnc(i,j) = hmcropext(i,j)
            hsdnc(i,j) = hsdcropext(i,j)
         enddo
      enddo
      err=NF_PUT_ATT_TEXT(ncidout,varidn(1),'long_name',40,
     &     EntPFT_title(k))
      err=NF_PUT_ATT_TEXT(ncidout,varidn(2),'long_name',40,
     &     EntPFT_title(k))
      startB(1)=1
      startB(2)=1
      countB(1)=IM
      countB(2)=JM
      err=NF_PUT_VARA_REAL(ncidout,varidn(1),startB,countB,hmnc)
      err=NF_PUT_VARA_REAL(ncidout,varidn(2),startB,countB,hsdnc)
      err=NF_CLOSE(ncidout)

!*************************************************************************
      !Generate nocrops
      !THIS ONLY WORKS FOR N_VEG=16 !!!
      ! remove combo crops 15 and rescale fractions so that they sum to 1
      vfc(:,:,15) = 0.
      laic(:,:,15) = 0.
!      vfm(:,:,:,15) = 0.
      laim(:,:,:,15) = 0.
!!      if ( N_VEG == 17 ) then   ! remove separate C4 crops 
!        vfc(:,:,16) = 0.       
      laic(:,:,16) = 0.        
!        vfm(:,:,:,16) = 0.
      laim(:,:,:,16) = 0.        
!!      endif
      flag = 0
      do while (flag.le.1)
         nonaturalcount=0
         LAYER(:,:) = 0
         do j=1,JMn
            do i=1,IMn
               s = sum( vfc(i,j,1:N_BARE) ) 
               if (sum(vfc(i,j,15:16))>0) then !Cell has crops
              !write(*,*) 'Cell has crops',i,j,s
!           if (sum(vfc(i,j,15:18)).eq.s) then !Cell only crops and bare
!              ! OLD HACK: if cell was all crop and bare, replace herb with C4 grass, and trees with broadleaf decid. 
!              vfc(i,j,12) = vfc(i,j,12) + vfc(i,j,15) !herb
!              vfc(i,j,6) = vfc(i,j,6) + vfc(i,j,16) !tree 
!               vfc(i,j,15) = 0.
!              vfc(i,j,16) = 0.
!              vfm(:,i,j,12) = vfm(:,i,j,12) + vfm(:,i,j,15)  !herb
!              vfm(:,i,j,6) = vfm(:,i,j,6) + vfm(:,i,j,16) !tree
!              vfm(:,i,j,15) = 0.
!              vfm(:,i,j,16) = 0.

                  if (sum(vfc(i,j,15:18)).eq.s) then
                     write(*,*) 'Cell is all crops + bare',i,j,s
                  endif
                  
!              if ((vfc(i,j,15).gt.0.).or.(vfc(i,j,16).gt.0.)) then 
              !Replace with closest non-zero natural cover
                  call replace_crops(IMn,JMn,KM,i,j,s
!     &             ,vfc,vfm,laic,laim,naturalfound)
     &                 ,vfc,vfm,vfc,vfm,naturalfound)
                  if (naturalfound.eq.0) then
                     nonaturalcount = nonaturalcount+1
                     LAYER(i,j) = 1.d0
                  endif
               endif
            enddo
         enddo
         
         do i=1,IMn
            do j=1,JMn
           ! Rescale always to account for removal of crops 
           ! and to make coasts sum to 1.
               s = sum( vfc(i,j,1:N_BARE) ) 
               if ( (s.gt.0.).and.(s.ne.1.) ) then
                  vfc(i,j,1:N_BARE) = vfc(i,j,1:N_BARE) / s
                  vfm(:,i,j,1:N_BARE) = vfm(:,i,j,1:N_BARE) / s
               endif
               
               s = sum( vfc(i,j,1:N_BARE) ) !#DEBUG
               if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then 
!              write(*,*) 'ERROR nocrops:  max and monthly lc differ'
!     &             ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
              !write(*,*) vfc(i,j,1:N_BARE)
              !write(*,*) vfm(1,i,j,1:N_BARE)
                  do m=1,12
                     vfm(m,i,j,1:N_BARE) = vfc(i,j,N_BARE)
                  enddo
               endif
            enddo
         enddo
         
         if (nonaturalcount.eq.0) then
            write(*,*) 'All crops cells successfully fixed.'
            flag = 2
         else
            flag = flag + 1
            if (flag.eq.1) then
               write(*,*) 'Some all-crop cells. Iterating once...'
     &           ,nonaturalcount
            elseif (flag.gt.1) then
               write(*,*) 'Remaining all-crop cells,',nonaturalcount
            endif
         endif
         titlefoo = 'All-crop cells with no near natural veg.'
         write(93) titlefoo, LAYER
      enddo                     !do while flag ------
      
      !## HACK -nk
      do m=1,12
         vfm(m,:,:,:) = vfc(:,:,:)
      enddo
      
#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif

!*************************************************************************
      call write_output(titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_scaled_nocrops_"
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &     ,"max_trimmed_scaled_nocrops","  ",res_out)
      write(*,*) "trimmed, scaled, no crops"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)
      
	  ! save netcdf lai_max_trimmed_scaled_nocrops
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_max_trimmed_scaled_nocrops.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               laicnc(i,j) = laic(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,laicnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
      
      ! save netcdf lc_max_trimmed_scaled_nocrops
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lc_max_trimmed_scaled_nocrops.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               vfcnc(i,j) = vfc(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,vfcnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
!*************************************************************************
      
      do m=1,12
         call write_output(titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE, 
     &        "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &        ,"trimmed_scaled_nocrops",MONTH(m),res_out)
      enddo
      
      ! save netcdf lai_trimmed_scaled_nocrops_MONTH
      do m=1,12
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lai_trimmed_scaled_nocrops_'//MONTH(m)//'.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         do k=1,18
            inqvarout = EntPFT_shorttitle(k)
            err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         enddo
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do k=1,18
            do i=1,IMn
               do j=1,JMn
                  laicnc(i,j) = laim(m,i,j,k)
               enddo
            enddo
            err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &           EntPFT_title(k))
            startB(1)=1
            startB(2)=1
            countB(1)=IM
            countB(2)=JM
            err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &           laicnc)
            write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         enddo
         err=NF_CLOSE(ncidout)
         write(*,*) err, 'nf_close out ',trim(fileoutnc)
      enddo
      
      ! save netcdf lc_trimmed_scaled_nocrops_MONTH
      do m=1,12
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lc_trimmed_scaled_nocrops_'//MONTH(m)//'.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         do k=1,18
            inqvarout = EntPFT_shorttitle(k)
            err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         enddo
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do k=1,18
            do i=1,IMn
               do j=1,JMn
                  vfcnc(i,j) = vfm(m,i,j,k)
               enddo
            enddo
            err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &           EntPFT_title(k))
            startB(1)=1
            startB(2)=1
            countB(1)=IM
            countB(2)=JM
            err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &           vfcnc)
            write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         enddo
         err=NF_CLOSE(ncidout)
         write(*,*) err, 'nf_close out ',trim(fileoutnc)
      enddo

!*************************************************************************
      call write_output_h(titleh, hm, hsd, N_BARE,
     &     "../lc_lai_ent16/V"//res_out_int//
     &     "_EntMM16_height_trimmed_scaled_nocrops.ij","   ",res_out)
      
	  ! save netcdf simard heights trimmed all
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_height_trimmed_scaled_nocrops.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'hgt_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         inqvarout = 'stdev_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varids(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               hmnc(i,j) = hm(i,j,k)
               hsdnc(i,j) = hsd(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,hmnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         err=NF_PUT_ATT_TEXT(ncidout,varids(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varids(k),startB,countB,hsdnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
      
      foo(:,:) = 0.
      do j=1,JMn
         do i=1,IMn
            do k=1,KM
               foo(i,j) = foo(i,j) + vfc(i,j,k)
            enddo
         enddo
      enddo
      titlefoo = "LC trimmed scaled nocrops checksum"
      write(100) titlefoo, foo
      
      ! write the final output
      do k=1,N_BARE
         write(7) titlec(k), vfc(:,:,k)
      enddo
      do k=1,N_BARE
         write(7) titlec(k), laic(:,:,k)
      enddo
      
      
      lai_yy(:,:) = 0.
      do j=1,JMn
         do i=1,IMn
            lai_yy(i,j) = laic(i,j,N_BARE)*vfc(i,j,N_BARE)
         enddo
      enddo
      
      write(8) title_yy, lai_yy(:,:)
      
      lai_yy(:,:) = 0.
      do k=1,KM
         if (k==1 .or. k==19) cycle
         lai_yy(:,:) = lai_yy(:,:) + vfn(:,:,k)*lain(:,:,k)
      enddo
      
      write(9) title_yy, lai_yy(:,:)
      
      lai_yy(:,:) = 0.
      do k=1,N_BARE
         lai_yy(:,:) = lai_yy(:,:) + vfc(:,:,k)*laic(:,:,k)
      enddo
      
      write(9) title_yy, lai_yy(:,:)
      
      
      
      ! find the biggest vf for cells which still have sparse veg
      vf_yy(:,:) = -1e30
      lai_yy(:,:) = 0.
      do j=1,JMn
         do i=1,IMn
            if( vfc(i,j,N_BARE) * laic(i,j,N_BARE) > .0 ) then
               
               maxpft = maxloc( vfc(i,j,1:16), 1 )
               
               vf_yy(i,j) = vfc(i,j,maxpft)
               lai_yy(i,j) = maxpft
               
            endif
         enddo
      enddo
      title_yy = "max fraction"
      write(10) title_yy, vf_yy(:,:)
      title_yy = "max pft"
      write(10) title_yy, lai_yy(:,:)
      
!*************************************************************************
      !------------------------------------------------------------
      ! convert monthly LAI to nocrops vfc.

      ! Merge crops.ext1 laimax, laimonthly, and height into nocrops
      ! to create laimax_ext1, laimonthly_ext1, and height_ext1 files.

      laic(:,:,15) = laicropext(:,:)
      laim(:,:,:,15) = laimcropext(:,:,:)
      hm(:,:,15) =  hmcropext(:,:)
      hsd(:,:,15) = hsdcropext(:,:)
      
      call write_output_lai(titlec, laic, N_BARE,
     &     "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &     ,"max_trimmed_scaled_ext1","  ",res_out)
      write(*,*) "trimmed, scaled, ext1"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)

	  ! save netcdf lai_max_trimmed_scaled_ext1 and crops
      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_lai_max_trimmed_scaled_ext1.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               laicnc(i,j) = laic(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,laicnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)
      
      
!*************************************************************************
      do m=1,12
         call write_output_lai(titlem(m,:)
     &        , laim(m,:,:,:), N_BARE, 
     &        "../lc_lai_ent16/V"//res_out_int//"_EntMM16"
     &        ,"trimmed_scaled_ext1",MONTH(m),res_out)
      enddo
      
      ! save netcdf
      do m=1,12
         fileoutnc =
     &        '../lc_lai_ent16/nc/V'//res_out_int//
     &        '_EntMM16_lai_trimmed_scaled_ext1_'//MONTH(m)//'.nc'
         err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
         write(*,*) err, 'nf_create out ',trim(fileoutnc)
         err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
         err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
         err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
         err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
         dd(1)=dimidx
         dd(2)=dimidy
         do k=1,18
            inqvarout = EntPFT_shorttitle(k)
            err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         enddo
         err=NF_ENDDEF(ncidout)
         err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
         err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
         do k=1,18
            do i=1,IMn
               do j=1,JMn
                  laicnc(i,j) = laim(m,i,j,k)
               enddo
            enddo
            err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &           EntPFT_title(k))
            startB(1)=1
            startB(2)=1
            countB(1)=IM
            countB(2)=JM
            err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,
     &           laicnc)
            write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         enddo
         err=NF_CLOSE(ncidout)
         write(*,*) err, 'nf_close out ',trim(fileoutnc)
      enddo
      
!*************************************************************************
      call write_output_h(titleh, hm, hsd, N_BARE,
     &     "../lc_lai_ent16/V"//res_out_int//
     &     "_EntMM16_height_trimmed_scaled_ext1.ij","   ",res_out)
      
	  ! save netcdf simard heights trimmed all

      fileoutnc =
     &     '../lc_lai_ent16/nc/V'//res_out_int//
     &     '_EntMM16_height_trimmed_scaled_ext1.nc'
      err = NF_CREATE(fileoutnc,NF_CLOBBER,ncidout)
      write(*,*) err, 'nf_create out ',trim(fileoutnc)
      err=NF_DEF_DIM(ncidout,'lon',IM,dimidx)
      err=NF_DEF_DIM(ncidout,'lat',JM,dimidy)
      err=NF_DEF_VAR(ncidout,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(ncidout,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,18
         inqvarout = 'hgt_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varidn(k))
         inqvarout = 'stdev_'//EntPFT_shorttitle(k)
         err=NF_DEF_VAR(ncidout,inqvarout,NF_REAL,2,dd,varids(k))
      enddo
      err=NF_ENDDEF(ncidout)
      err=NF_PUT_VARA_REAL(ncidout,varidx,1,IM,long)
      err=NF_PUT_VARA_REAL(ncidout,varidy,1,JM,lati)
      do k=1,18
         do i=1,IMn
            do j=1,JMn
               hmnc(i,j) = hm(i,j,k)
               hsdnc(i,j) = hsd(i,j,k)
            enddo
         enddo
         err=NF_PUT_ATT_TEXT(ncidout,varidn(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varidn(k),startB,countB,hmnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
         err=NF_PUT_ATT_TEXT(ncidout,varids(k),'long_name',40,
     &        EntPFT_title(k))
         startB(1)=1
         startB(2)=1
         countB(1)=IM
         countB(2)=JM
         err=NF_PUT_VARA_REAL(ncidout,varids(k),startB,countB,hsdnc)
         write(*,*) err,'nf_put_var_real out ',trim(inqvarout)
      enddo
      err=NF_CLOSE(ncidout)
      write(*,*) err, 'nf_close out ',trim(fileoutnc)

      close(100) !Checksum foo

      end program convert

