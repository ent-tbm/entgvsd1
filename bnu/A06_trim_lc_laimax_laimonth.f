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
      covmax = 0.d0
      covmaxk = 0
      covmaxii = i
      covmaxjj = j
      do k=1,14                 !Find max non-crop, non-bare natural cover type
         if (vfc(i,j,k) > covmax)  then
            covmax = vfc(i,j,k)
            covmaxk = k
         endif
      enddo
      if (covmax.gt.0.d0) then !Assign dominant natural veg to crop
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
         covsum(:) = 0.d0
         covavglai(:) = 0.d0
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
         covmax = 0.d0
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
!      if ((IMn.eq.144).and.(JMn.eq.90)) then
      if ((IMn.eq.1440).and.(JMn.eq.720)) then
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
         lat(j) = -90.25 + 0.25/2 + .25*j
      enddo

      end subroutine calc_lon_lat_025x025

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

!      subroutine omp_set_num_threads(num_threads)
!      integer, intent(in) :: num_threads

!------------------------------------------------------------------------

      program convert
      use omp_lib
      use conversions
      
      implicit none
      
      include 'netcdf.inc'


      !***************************************************
      !*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
      !***************************************************
      character*17, parameter :: EntPFT_files(20) =
     &     (/
     &     'water            ',
     &     '01_ever_br_early ',
     &     '02_ever_br_late  ',
     &     '03_ever_nd_early ',
     &     '04_ever_nd_late  ',
     &     '05_cold_br_early ',
     &     '06_cold_br_late  ',
     &     '07_drought_br    ',
     &     '08_decid_nd      ',
     &     '09_cold_shrub    ',
     &     '10_arid_shrub    ',
     &     '11_c3_grass_per  ',
     &     '12_c4_grass      ',
     &     '13_c3_grass_ann  ',
     &     '14_c3_grass_arct ',
     &     '15_crops_c3_herb ',
     &     '16_crops_c4_herb ',
     &     '17_crops_woody   ',
     &     '18_snow_ice      ',
     &     '19_bare_sparse   '
     &     /)

            character*14, parameter :: EntPFT_lcfields(20) =
     &     (/
     &     'water_lc      ',
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

            character*14, parameter :: EntPFT_laifields(20) =
     &     (/
     &     'water_lai     ',
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
      integer, parameter :: IM1km=43200, JM1km=21600
      integer, parameter :: IM=432, JM=216, KM=20
!      integer, parameter :: IM=1440, JM=720, KM=20
!      integer, parameter :: IM=144, JM=90, KM=20
!      character*(*), parameter :: res_in="2.5x2"
!      character*(*), parameter :: res_in="144x90"
!      character*(*), parameter :: res_in_int="144x90"
!      character*(*), parameter :: res_in="1440x720"
!      character*(*), parameter :: res_in_int="1440x720"
      character*(*), parameter :: res_in="1kmx1km"
      character*(*), parameter :: res_out="1kmx1km"

      integer, parameter :: divx = IM
      integer, parameter :: divy = JM
      integer, parameter :: ndiv = 100

      integer :: nthreads, tid
!      integer :: omp_get_num_threads,omp_get_threads_num

      integer :: startx(ndiv),starty(ndiv)

      integer :: fileid,dimidx,dimidy,dimidz,dd(4),varidx,varidy,varidz
      integer :: startA(1),startB(2),countA(1),countB(2),lenx,leny

      integer :: start2d(2),count2d(2),start3d(4),count3d(4)
      real*4 :: lcin(IM,JM),laiin(IM,JM),hin(IM,JM),hstd(IM,JM)
      integer :: fileidlcmax(KM),fileidlaimax(KM)
      integer :: varidlcmax(KM),varidlaimax(KM)
      integer :: fileidlaim(12),fileidlcm(12)
      integer :: varidlaim(12),varidlcm(12)
      character*80 :: fileheight,filestd
      character*80 :: filelaim(12),filelcm(12)
      character*80 :: filelcmax(KM), filelaimax(KM)
      integer :: fileidheight,varidheight
      integer :: fileidstd,varidstd,err
      character*80 :: filebs
      integer :: fileidbs,varidbs
      character*60 :: PathFilepre, PathFilepost

      character*80 :: fllaimaxpure
      integer :: fidlaimaxpure
      integer :: varidlaimaxpure(18)

      integer, parameter :: N_HT = 19 ! number of height layers input

!      real*4 vf(IM,JM,KM), lai(IM,JM,KM)
      character*80 :: title(KM), title_tmp,title12(12,KM),titlehn(2,KM)

      ! Input values
      ! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
      integer, parameter :: IMn=IM, JMn=JM

      character*(*), parameter :: file_checksum = 
     &     "../lc_lai_ent16/EntMM16_checksum_"//res_out//".ij"

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

      integer i, j, k, io, in, jn, maxpft, kx, m
      real*8 lat
      real*4 foo(IMn,JMn)
      integer N_VEG             ! number of PFTs in output
      integer N_BARE            ! index of bare soil in output
      integer count

      integer :: ncidin,ncidout,varid,status,varidn(18)
      integer :: varids(18)
      character*20 :: inqvarin, inqvarout
      real*4 :: long(IM1km),lati(JM1km)
      character*256 :: fileoutnc
      real*4, dimension(12) :: time
      time = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
            
      !define Lon and Lat
      call calc_lon_lat(IM1km,JM1km,long,lati)
      
      do i=1,ndiv
         startx(i) = -divx + 1 + divx*i
      enddo

      do i=1,ndiv
         starty(i) = -divy + 1 + divy*i
      enddo

      !------------------------------------------------------------------------
      ! OPEN INPUT FILES
      
      ! lcmax
      do k=1,20
         PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
	 PathFilepost=  EntPFT_files(k)
         filelcmax(k) = trim(PathFilepre)//trim(PathFilepost)//'_lc.nc'
         err = NF_OPEN(filelcmax(k),0,fileidlcmax(k))
         write(*,*) err,fileidlcmax(k)
         inqvarin = EntPFT_lcfields(k)
         err = NF_INQ_VARID(fileidlcmax(k),inqvarin,varidlcmax(k))
         write(*,*) err,fileidlcmax(k),varidlcmax(k),EntPFT_lcfields(k)
      enddo

      ! laimax
      do k=1,20
         PathFilepre= '../lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
	 PathFilepost=  EntPFT_files(k)
         filelaimax(k) = trim(PathFilepre)//trim(PathFilepost)//
     &        '_lai.nc'
         err = NF_OPEN(filelaimax(k),0,fileidlaimax(k))
         write(*,*) err,filelaimax(k)
         inqvarin = EntPFT_laifields(k)
         err = NF_INQ_VARID(fileidlaimax(k),inqvarin,varidlaimax(k))
         write(*,*) err
      enddo

      ! lai monthly
      do m=1,12
         PathFilepre= '../lc_lai_ent/nc/'
         PathFilepost= 'EntMM_lc_lai_'//trim(MONTH(m))//
     &        '_1kmx1km.nc'
         filelaim(m) = trim(PathFilepre)//trim(PathFilepost)
         err = NF_OPEN(filelaim(m),0,fileidlaim(m))
	 write(*,*) err, 'lai monthly files'
         inqvarin = 'EntPFT'
         err = NF_INQ_VARID(fileidlaim(m),inqvarin,varidlaim(m))
         write(*,*) err
      enddo

      ! height file - insert dummy WATER layer at beginning to avoid
      !               confusion in numbering with vfc and vfm.
      PathFilepre= '../../data/height/'
      PathFilepost= 'EntGVSDmosaic17_height_1kmx1km_lai3g.nc'
      fileheight = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(fileheight,0,fileidheight)
      inqvarin = 'SimardHeights'
      err = NF_INQ_VARID(fileidheight,inqvarin,varidheight)
      write(*,*) err, 'height'
      
      PathFilepre= '../../data/height/'
      PathFilepost= 'EntGVSDmosaic17_heightstd_1kmx1km.nc'
      filestd = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filestd,0,fileidstd)
      !write(*,*) err
      inqvarin = 'SimardHeights'
      err = NF_INQ_VARID(fileidstd,inqvarin,varidstd)
      write(*,*) err

      ! bs ratio
      PathFilepre= '../lc_lai_ent/'
      PathFilepost= 'bs_brightratio.nc'
      filebs = trim(PathFilepre)//trim(PathFilepost)
      err = NF_OPEN(filebs,0,fileidbs)
      write(*,*) err,filebs
      inqvarin = 'bs_brightratio'
      err = NF_INQ_VARID(fileidbs,inqvarin,varidbs)
      write(*,*) err, 'bs brightratio'

      ! lc & laimax file of Ent 17 PFTs: JUST FOR TITLE(K)
      open(1,file="../lc_lai_ent/EntMM_lc_laimax_144x90.bin",
     &     form="unformatted",status="old")
      do k=1,KM 
         read(1) title(k)
         !write(*,*) title(k)
      enddo
      close(1)
      
      do m=1,12
         open(1,file="../lc_lai_ent/EntMM_lc_laimax_144x90.bin",
     &        form="unformatted",status="old")
         do k=1,KM
            read(1) title12(m,k)
            !write(*,*) title12(m,k), KM
         enddo
         close(1)
      enddo

      ! simard heights titles
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

      nthreads = 28
      call omp_set_num_threads(nthreads)
      write(*,*) omp_get_num_procs(),'number of processors'
      write(*,*) omp_get_max_threads(),'max size of thread team'
      write(*,*) omp_get_num_threads(),'size of thread team'

      !------------------------------------------------------------------------
      ! CREATE OUTPUT NETCDF FILES
      ! save netcdf lai_max_pure
      fllaimaxpure =
     &     '../lc_lai_ent16/nc/V'//res_out//
     &     '_EntMM16_lai_max_pure.nc'
      err = NF_CREATE(fllaimaxpure,NF_CLOBBER,fidlaimaxpure)
      write(*,*) err, 'nf_create out ',trim(fllaimaxpure)
      err=NF_DEF_DIM(fidlaimaxpure,'lon',IM1km,dimidx)
      err=NF_DEF_DIM(fidlaimaxpure,'lat',JM1km,dimidy)
      err=NF_DEF_VAR(fidlaimaxpure,'lon',NF_REAL,1,dimidx,varidx)
      err=NF_DEF_VAR(fidlaimaxpure,'lat',NF_REAL,1,dimidy,varidy)
      dd(1)=dimidx
      dd(2)=dimidy
      do k=1,1
         inqvarout = EntPFT_shorttitle(k)
         err=NF_DEF_VAR(fidlaimaxpure,inqvarout,NF_REAL,2,dd,
     &       varidlaimaxpure(k))
         write(*,*) err,inqvarout
      enddo
      err=NF_ENDDEF(fidlaimaxpure)


      !------------------------------------------------------------------------
      ! OMP parallelization 

!$OMP& parallel do private(i,j,k,m,ndiv,err,N_HT), &
!$OMP& shared(startx,starty,divx,divy,fileidlcmax,varidlcmax), &
!$OMP& shared(start2d,count2d,start3d,count3d), &
!$OMP& shared(vfn,lcin,lain,laiin,vfn,vfnm,laiinm), &
!$OMP& shared(fileidlaim,varidlaim,hmn,fileidheight,varidheight,hin,hmn), &
!$OMP& shared(fileidstd,varidstd,hstd,hsdn,hstd), &
!$OMP& shared(fileidbs,varidbs,bs_brightratio), &
!$OMP& shared(fileidbs,varidbs,bs_brightratio)


      tid = omp_get_thread_num()
      print *, 'thread = ', tid

 !    Only master thread does this
      if (tid .eq. 0) then
        nthreads = omp_get_num_threads()
        print *, 'Number of threads = ', nthreads
      endif

      do j=1,1!ndiv

         do i=1,1!ndiv

            start2d(1) = startx(i)
            start2d(2) = starty(j)
            count2d(1) = divx
            count2d(2) = divy

            do k=1,KM

	       ! get lc max
               err = NF_GET_VARA_REAL(fileidlcmax(k),varidlcmax(k),
     &              start2d,count2d,lcin)
               write(*,*) err, 'lc max input',k,startx(i),starty(j)
               vfn(:,:,k)=lcin

               ! get lai max
               err = NF_GET_VARA_REAL(fileidlaimax(k),varidlaimax(k),
     &              start2d,count2d,laiin)
               write(*,*) err, 'lai max input'
               lain(:,:,k)=laiin

            enddo

            ! get monthly lc
	    do m=1,12
               vfnm(m,:,:,:)=vfn(:,:,:)
            enddo

            ! get monthly lai
            do m=1,12
               do k=1,20
                  start3d(1) = k
                  start3d(2) = startx(i)
                  start3d(3) = starty(j)
                  count3d(1) = 1
                  count3d(2) = divx
                  count3d(3) = divy
                  err = NF_GET_VARA_REAL(fileidlaim(m),varidlaim(m),
     &                 start3d,count3d,laiin)
                  lainm(m,:,:,k)=laiin
               enddo
               write(*,*) err, 'lai monthly fields ',MONTH(m)
            enddo

            ! height file - insert dummy WATER layer at beginning toi avoid
            ! confusion in numbering with vfc and vfm.
            hmn(:,:,:) = 0.
            do k=1,N_HT
               start3d(1) = k
               start3d(2) = startx(i)
               start3d(3) = starty(j)
               count3d(1) = 1
               count3d(2) = divx
               count3d(3) = divy
               err = NF_GET_VARA_REAL(fileidheight,varidheight,start3d,
     &              count3d,hin)
               write(*,*) err, 'height ',k
               hmn(:,:,k+1)=hin
            enddo

            hsdn(:,:,:) = 0.
            do k=1,N_HT
               start3d(1) = k
               start3d(2) = startx(i)
               start3d(3) = starty(j)
               count3d(1) = 1
               count3d(2) = divx
               count3d(3) = divy
               err = NF_GET_VARA_REAL(fileidstd,varidstd,start3d,
     &              count3d,hstd)
               write(*,*) err, 'height ',k
               hsdn(:,:,k+1)=hstd
            enddo

	    ! get bs bright ratio
            err = NF_GET_VARA_REAL(fileidbs,varidbs,
     &            start2d,count2d,bs_brightratio)
            write(*,*) err, 'bs bright ratio'

         !Check if mismatch lc or lai values (one is zero and the other not)
          call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)


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



         enddo
      enddo

!$OMP& end parallel do




      
      end program convert

