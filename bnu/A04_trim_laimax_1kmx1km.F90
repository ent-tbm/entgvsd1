
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

#define COMBINE_CROPS_C3_C4
#define SPLIT_BARE_SOIL


module conversions


implicit none

CONTAINS

! -----------------------------------------------------------------
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
! -----------------------------------------------------------------
!    call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE), &
!       vfm(m,i,j,9),vfm(m,i,j,9), vfc(i,j,9))
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
! -----------------------------------------------------------------
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
! -----------------------------------------------------------------



subroutine write_output_single(titlefoo, vf, lai, fileprefix &
     , MISC, MON, resoutt)
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

open(80,file=fileprefix//"_lc_"// &
     trim(MISC)//trim(MONstr)//".ij", &
     form="unformatted",status="unknown")


title = titlefoo(1:48)//"  "//trim(MON)//" (cover)  "//resoutt
write(80) title, vf(:,:)
close(80)

open(90,file=fileprefix//"_lai_"// &
     trim(MISC)//trim(MONstr)//".ij", &
     form="unformatted",status="unknown")



title = titlefoo(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
write(90) title, lai(:,:)

close(90)

end subroutine write_output_single


subroutine write_output(titlec, vfc, laic, n, fileprefix &
     , MISC, MON, resoutt)
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

open(80,file=fileprefix//"_lc_"// &
     trim(MISC)//trim(MONstr)//".ij", &
     form="unformatted",status="unknown")

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

open(90,file=fileprefix//"_lai_"// &
     trim(MISC)//trim(MONstr)//".ij", &
     form="unformatted",status="unknown")

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


subroutine write_output_h_single(titleh, h,hsd, filename &
     , MISC, resoutt)
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

open(90,file=filename, &
     form="unformatted",status="unknown")

title = titleh(1)(1:63)//"  "//trim(resoutt)
write(90) title, h(:,:)

title = titleh(2)(1:63)//"  "//trim(resoutt)
write(90) title, hsd(:,:)

close(90)
end subroutine write_output_h_single


subroutine write_output_h(titleh, h,hsd, n, filename &
     , MISC, resoutt)
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

open(90,file=filename, &
     form="unformatted",status="unknown")

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



subroutine get_bare_soil_brightratio(IMn,JMn, filename &
     ,bs_brightratio)
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


open(1,file=filename, &
     form="unformatted",status="old")

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

!        print *, "count=", count
  if ( count== 0 ) exit
enddo
title_bs = "bare soil brightratio"
write(91) title_bs, bs_brightratio
title_bs = "bare soil fraction"
write(91) title_bs, bsf


end subroutine get_bare_soil_brightratio


subroutine split_bare_soil(N_VEG,KM,N_BARE &
     ,bs_brightratio,vfc,laic, &
     res_out)
!Split BARE soil into BRIGHT and DARK cover to preserve albedo from
!  "old" ModelE cover.  Should be called after each trim, scale, nocrops.
!Any LAI on BARE soil should already have been moved to vegetated cover,
!  so laic(:,:,N_BARE) should be zero.
!This checks for cases if BARE is original total or was previously split.

integer, intent(in) :: N_VEG, KM
integer, intent(inout) :: N_BARE
real*4, intent(in) :: bs_brightratio !Fraction of bare that is bright.
real*4 vfc(:), laic(:) !(KM)
!      real*4 vfm(:,:), laim(:,:) !(12,KM)
!      real*4 vfh(:) !(KM),
!      real*4 hm(:), hsd(:) !(KM)
!      character*80 :: titlec(:) !(18)
!      character*80 :: titlem(:,:) !(12,18)
!      character*80 :: titleh(:,:) !(2,18) !1-h, 2-hsd
character*(*) :: res_out
!-----Local----
real*4 :: vfc_bare, vf_tot
integer :: m

if ((N_BARE-N_VEG).eq.1) then !First time splitting
! increase number of fractions by 1 so that
! vfc(:,:,N_BARE-1) = fraction of bright soil
! vfc(:,:,N_BARE) = fraction of dark soil
   vf_tot = vfc(N_BARE)  !All bare soil orig in N_BARE layer
   vfc(N_BARE) = vf_tot*bs_brightratio
   vfc(N_BARE+1) = vf_tot - vfc(N_BARE)
   laic(N_BARE:N_BARE+1) = 0.
!         do m=1,12
!            vf_tot = vfm(m,N_BARE)
!            vfm(m,N_BARE) = vf_tot*bs_brightratio
!            vfm(m,N_BARE+1) = vf_tot - vfm(m,N_BARE)
!            laim(m,N_BARE:N_BARE+1) = 0.
!         enddo
   N_BARE = N_BARE+1
else  !N_BARE-NVEG.eq.2, previously split and incremented N_BARE
   vf_tot = vfc(N_BARE-1) + vfc(N_BARE)
   vfc(N_BARE-1) = vf_tot*bs_brightratio
   vfc(N_BARE) = vf_tot - vfc(N_BARE-1)
   laic(N_BARE-1:N_BARE) = 0.
!         do m=1,12
!            vf_tot = vfm(m,N_BARE-1) + vfm(m,N_BARE)
!            vfm(m,N_BARE-1) = vf_tot*bs_brightratio
!            vfm(m,N_BARE) = vf_tot - vfm(m,N_BARE-1)
!            laim(m,N_BARE-1:N_BARE) = 0.
!         enddo
endif

!      vfh(N_BARE-1) = vfc(N_BARE-1)
!      vfh(N_BARE) = vfc(N_BARE)
!      hm(N_BARE) = 0.
!      hsd(N_BARE) = 0.

!      titlec(N_BARE-1) =   "17 - bright bare soil  "//
!     &     "                           (COVER) "//res_out
!      titlec(N_BARE) =     "18 - dark bare soil    "//
!     &     "                           (COVER) "//res_out
!      do m=1,12
!            titlem(m,N_BARE-1) = "17 - bright bare soil  "//
!     &           "                          "//MONTH(m)//" (cover) "
!     &           //res_out
!            titlem(m,N_BARE) =   "18 - dark bare soil    "//
!     &           "                          "//MONTH(m)//" (cover) "
!     &           //res_out
!         enddo
!      titleh(1,N_BARE-1) =   "17 - bright bare soil  "//
!     &     "height (m)                                "//res_out
!      titleh(1,N_BARE) =     "18 - dark bare soil    "//
!     &     "height (m)                                "//res_out
!      titleh(2,N_BARE-1) =   "17 - bright bare soil  "//
!     &     "height stdev(m)                           "//res_out
!      titleh(2,N_BARE) =     "18 - dark bare soil    "//
!     &     "height stdev(m)                           "//res_out

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
      L1(:,:) = IJADD4(1,IMn,1,JMn,L1(:,:), &
           IJMULT4(1,IMn,1,JMn,vfs(:,:,k),lais(:,:,k)))
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
character*3, parameter :: MON(12) = &
     (/ &
     "Jan","Feb","Mar","Apr","May","Jun", &
     "Jul","Aug","Sep","Oct","Nov","Dec" &
     /)

integer :: i,j,k,m
character*80 :: TITLE
real*4 :: L1(12,IMn,JMn)

L1(:,:,:) = 0.0
do m=1,12
   do k=1,Ln
      L1(m,:,:) = IJADD4(1,IMn,1,JMn,L1(m,:,:), &
           IJMULT4(1,IMn,1,JMn,vfmon(m,:,:,k),laimon(m,:,:,k)))
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
         if (((vf(i,j,k).eq.0.0).and.(lai(i,j,k).gt.0.0)) &
              .or.(vf(i,j,k).gt.0.0).and.(lai(i,j,k).eq.0.0)) &
              then
!                  print *,label,' mismatch,k,i,j,cov,lai'
!     &                 ,k,i,j,vf(i,j,k),lai(i,j,k)
!                  LAYER(i,j) = vf(i,j,k)
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
   vfc(i,j,covmaxk) = vfc(i,j,covmaxk)  &
        + vfc(i,j,15) + vfc(i,j,16)
   vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk) &
        + vfm(:,i,j,15) + vfm(:,i,j,16)
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
         if ( (ii.ge.1).and.(ii.le.IMn) &
              .and.(jj.ge.1).and.(jj.le.JMn) !in grid range &
              .and.((ii.ne.i).or.(jj.ne.j)) ) !not the i,j center cell &
              then
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
               covavglai(k) = covavglai(k) +  &
                    laic(ii,jj,k)*vfc(ii,jj,k)
               covsumm(:,k) = covsumm(:,k) + vfm(:,ii,jj,k)
               covavglaim(:,k) = covavglaim(:,k) + &
                       laim(:,ii,jj,k)*vfm(:,ii,jj,k)
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
      vfc(i,j,covmaxk) = vfc(i,j,covmaxk)  &
           + vfc(i,j,15) + vfc(i,j,16)
      vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk) &
           + vfm(:,i,j,15) + vfm(:,i,j,16)
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
   if ( (((i.eq.34).or.(i.eq.41).or.(i.eq.48)).and.(j.eq.9)) &
        .or.( (i.eq.34).and.(j.eq.10) ) ) then !Antartic
      vfc(i,j,14) = vfc(i,j,14) + vfc(i,j,15) + vfc(i,j,16)
      vfm(:,i,j,14) = vfm(:,i,j,14)  &
           + vfm(:,i,j,15) + vfm(:,i,j,16)
      laic(i,j,14) = (laic(i,j,15)*vfc(i,j,15)  &
           + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
      do m=1,12
         laim(m,i,j,14) = ( laim(m,i,j,15)*vfm(m,i,j,15) + &
              laim(m,i,j,16)*vfm(m,i,j,16) ) /  &
              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
      enddo
      write(*,*) 'Replaced Antarctic crops.'
      naturalvegfound=1
   elseif ( ((i.eq.98).and.(j.eq.36)) &
           .or.((i.eq.139).and.(j.eq.45)) ) then !tropical islands
      vfc(i,j,12) = vfc(i,j,12) + vfc(i,j,15) + vfc(i,j,16)
      laic(i,j,12) = (laic(i,j,15)*vfc(i,j,15)  &
           + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
      do m=1,12
         laim(m,i,j,12) = ( laim(m,i,j,15)*vfm(m,i,j,15) + &
              laim(m,i,j,16)*vfm(m,i,j,16)) /  &
              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
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
     ,vfm15,laim15,hm15,hsd15
     ,laiccrop, laimcrop, hmcrop, hsdcrop)
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
               if ( (ii.ge.1).and.(ii.le.IMn) &
                    .and.(jj.ge.1).and.(jj.le.JMn) !in grid range &
                    .and.((ii.ne.i).or.(jj.ne.j)) ) !not in i,j  &
                    then
                  if (vfc15(ii,jj).gt.0.d0) then
                     covsum15 = covsum15 + vfc15(ii,jj)
                     laiavg15 = laiavg15  &
                          + laic15(ii,jj)*vfc15(ii,jj)
                     hmavg15 = hmavg15 + hm15(ii,jj)*vfc15(ii,jj)
                     hsdavg15 = hsdavg15  &
                          + hsd15(ii,jj)*vfc15(ii,jj)
                     covmsum15(:) = covmsum15(:) + vfm15(:,ii,jj)
                     laimavg15(:) = laimavg15(:) + &
                          laim15(:,ii,jj)*vfm15(:,ii,jj)
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
!------------------------------------------------------------------------

!      subroutine omp_set_num_threads(num_threads)
!      integer, intent(in) :: num_threads

!------------------------------------------------------------------------

program convert

use conversions
use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use entgvsd_netcdf_util
use ent_labels_mod

implicit none



!***************************************************
!*      ENT PLANT FUNCTIONAL TYPES - short names   *
!***************************************************
character*50, parameter :: EntPFT_shorttitle(18) = &
     (/ &
     "ever_br_early   ", &
     "ever_br_late    ", &
     "ever_nd_early   ", &
     "ever_nd_late    ", &
     "cold_br_early   ", &
     "cold_br_late    ", &
     "drought_br      ", &
     "decid_nd        ", &
     "cold_shrub      ", &
     "arid_shrub      ", &
     "c3_grass_per    ", &
     "c4_grass        ", &
     "c3_grass_ann    ", &
     "c3_grass_arct   ", &
     "crops_herb      ", &
     "crops_woody     ", &
     "bare_bright     ", &
     "bare_dark       " &
     /)


!***************************************************
!*      ENT PLANT FUNCTIONAL TYPES                 *
!***************************************************
character*40, parameter :: EntPFT_title(18) = &
     (/ &
     '1 - Evergreen Broadleaf Early Succ      ', &
     '2 - Evergreen Broadleaf Late Succ       ', &
     '3 - Evergreen Needleleaf Early Succ     ', &
     '4 - Evergreen Needleleaf Late Succ      ', &
     '5 - Cold Deciduous Broadleaf Early Succ ', &
     '6 - Cold Deciduous Broadleaf Late Succ  ', &
     '7 - Drought Deciduous Broadleaf         ', &
     '8 - Deciduous Needleleaf                ', &
     '9 - Cold Adapted Shrub                  ', &
     '10 - Arid Adapted Shrub                 ', &
     '11 - C3 Grass Perennial                 ', &
     '12 - C4 Grass                           ', &
     '13 - C3 Grass Annual                    ', &
     '14 - Arctic C3 Grass                    ', &
     '15 - Crops Herb                         ', &
     '16 - Crops Woody                        ', &
     '17 - Bright Bare Soil                   ', &
     '18 - Dark Bare Soil                     ' &
     /)


real*4, parameter :: undef = -1.e30

integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees

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
integer, parameter :: IM=1, JM=1, KM=20
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

type(Chunker_t) :: chunker
integer :: ichunk,jchunk, ic,jc
! Input files
type(ChunkIO_t) :: io_lcin(KM), io_laiin(KM), io_bs
! Output files
type(ChunkIO_t) :: ioall_laiout, io_laiout(18)
type(ChunkIO_t) :: ioall_sum, io_sum_lc, io_sum_lai


real*4 :: lcin,laiin,hin,hstd

character*80 :: title_tmp,title12(12,KM),titlehn(2,KM)

! Input values
! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
integer, parameter :: IMn=IM, JMn=JM

! Input values, max, monthly
real*4 vfn(KM), lain(KM), area
real*4 a
real*4 vfnm(12,KM),lainm(12,KM),aream
real*4 am
real*4 hmn(KM),hsdn(KM)
! Converted values
real*4 vfc(KM), laic(KM), dvf, s
character*80 :: titlec(18)
! Converted values - monthly
real*4 vfm(12,KM), laim(12,KM)
real*4 laimnc(KM), vfmnc(KM),laicropnc(12)

character*80 :: titlem(12,18)
! Converted values - heights
real*4 vfh(KM),hm(KM), hsd(KM)

character*80 :: titleh(2,18) !1-h, 2-hsd
! Converted values - crop ext 
real*4 laicropext, laimcropext(12)
real*4 hmcropext,hsdcropext
! Vars for calculating nocrops
integer naturalfound, flag, nonaturalcount !if no natural veg in vicinity

real*4 vf_xx, lai_xx
real*4 vf_yy, lai_yy
real*4 LAYER
character*80 :: title_xx="xx"
character*80 :: title_yy="yy"
character*80 :: titlefoo

! bs_brightratio = bare soil brightratio
real*4 :: bs_brightratio, vfc_tmp

!integer i, j, 
integer :: k, io, in, jn, maxpft, kx, m
real*8 lat
!real*4 foolc(IM1km,JM1km),foolai(IM1km,JM1km)
integer N_VEG             ! number of PFTs in output
integer N_BARE            ! index of bare soil in output
integer count

call init_ent_labels
call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 100, 120)

!------------------------------------------------------------------------
! OPEN INPUT FILES

! lcmax
do k=1,20
    ! TODO: LAI3g????
    ! PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
    call chunker%nc_open(io_lcin(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(ent20(k)%file1)//trim(ent20(k)%file2)//'_lc.nc', &
        trim(ent20(k)%file2), 1)
end do

! laimax
do k=1,20
    call chunker%nc_open(io_laiin(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(ent20(k)%file1)//trim(ent20(k)%file2)//'_lai.nc', &
        trim(ent20(k)%file2), 1)
enddo


! bs ratio
call chunker%nc_create(io_bs, &
    weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale ???
    '', 'bs_brightratio', 'bs_brightratio', &
    'Bare Soil Bright Ratio', '1', 'BrightRatio')

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! CREATE OUTPUT NETCDF FILES
! laimax_pure
call chunker%nc_create(ioall_laiout, &
    weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
    '16/nc/', 'V1km_EntGVSDv1.1_BNU16_laimax_pure')
do k=1,18
    call chunker%nc_reuse_file(ioall_laiout, io_laiout(k), &
        'lai_'//trim(ent18(k)%file2), trim(ent18(k)%title), &
        'm2 m-2', trim(ent18(k)%title), &
        weighting(chunker%wta1,1d0,0d0))
end do

!  checksum land  laimax
call chunker%nc_create(ioall_sum, &
    weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
    '16/nc/', 'V1km_EntGVSDv1.1_LAI3g16_laimax_pure_checksum')

call chunker%nc_reuse_file(ioall_sum, io_sum_lc, &
    'lc_checksum', 'Checksum of LC', '1', 'checksum - Land Cover', &
    weighting(chunker%wta1, 1d0, 0d0))

call chunker%nc_reuse_file(ioall_sum, io_sum_lai, &
    'lai_checksum', 'Checksum of LAI', 'm2 m-2', 'checksum - LAI', &
    weighting(chunker%wta1, 1d0, 0d0))

call chunker%nc_check

! Use these loop bounds for testing...
! it chooses a land area in Asia
!do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
!do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1

do jchunk = 1,nchunk(2)
do ichunk = 1,nchunk(1)

    call chunker%move_to(ichunk,jchunk)

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

        do k = 1,KM
            lcin = io_lcin(k)%buf(ic,jc)
            vfn(k)=lcin

            ! get lai max
            laiin = io_laiin(k)%buf(ic,jc)
            lain(k) = laiin         
        end do

        ! height file - insert dummy WATER layer at beginning toi avoid
        ! confusion in numbering with vfc and vfm.

        ! get bs bright ratio
        bs_brightratio = io_bs%buf(ic,jc)

        ! Check if mismatch lc or lai values (one is zero and the other not)
        ! call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)
      

        !* Convert to GISS 16 pfts format

        ! first 14 pfts though grass are the same, ignore WATER
        !  lc laimax
        vfc(1:14) = vfn(2:15)
        laic(1:14) = lain(2:15)
        ! titlec(1:14) = title(2:15)
        !  lc lai monthly
        !            do m=1,12
        !               titlem(m,1:14) = title12(m,2:15)
        !            enddo
        !  heights
        !            vfh(1:14) = vfn(2:15) !Should be the same cover from MODIS.
        !            hm(1:14) = hmn(2:15) 
        !            hsd(1:14) = hsdn(2:15)
        !            titleh(1,1:14) = titlehn(1,2:15)
        !            titleh(2,1:14) = titlehn(2,2:15)
      
        ! crops

#ifdef COMBINE_CROPS_C3_C4
        !lc laimax
        a = vfn(16) + vfn(17)
      
        if ( a > 0. ) then
            laic(15) = (vfn(16)*lain(16) &
              + vfn(17)*lain(17)) / a
            vfc(15) = a
        else
            laic(15) = 0.
            vfc(15) = 0.
        endif
        !heights - DO NOT AVERAGE. PRESERVE HEIGHTS. LAI will scale density
        !            a = vfn(16) + vfn(17) !input cover
        !            if ( a > 0. ) then
        !               hm(i,j,15) = (vfn(i,j,16)*hmn(i,j,16)
        !     &              + vfn(i,j,17)*hmn(i,j,17)) / a
        !               if ((hmn(16)>0.).and.(hmn(17)>0.)) then
                    !average if both exist
        !                  hm(15) = (vfn(16)*hmn(16)
        !     &                 + vfn(17)*hmn(17)) / a
        !               else
                    !don't average if only one or none exists
        !                  hm(15) = max(hmn(16),hmn(17))
        !               endif
                 !Sum of squares for sd.  Don't weight if only one or less exists
        !               if ((hmn(16)>0.).and.(hmn(17)>0.)) then
        !                  hsd(15) = sqrt((vfn(16)*hsdn(16)**2
        !     &                 + vfn(17)*hsdn(17)**2) / a)
        !               else
        !                  hsd(15) = max(hsd(16),hsd(17))
        !               endif
        !               vfh(15) = a
        !            else
        !               hm(15) = 0.
        !               hsd(15) = 0.
        !               vfh(15) = 0.
        !            endif
              
        !            write(*,*) "Re-doing crops.."
        !            titlec(15) = title(16)
        !            titlec(15)(1:18) = "15 - crops herb   "
        !            do m=1,12
        !               titlem(m,15) = title12(m,16)
        !               titlem(m,15)(1:18) = "15 - crops herb   "
        !            enddo
        !            titleh(1,15) = titlehn(1,16)
        !            titleh(2,15) = titlehn(2,16)
        !           titleh(1,15)(1:18) = "15 - crops herb   "
        !            titleh(2,15)(1:18) = "15 - crops herb   "
              
        !            write(*,*) titlec(15)
        ! crops woody
        vfc(16) = vfn(18)
        laic(16) = lain(18)
        !            titlec(16) = '16 - '//title(18)(6:80)
        !            write(*,*) "titlec: "
        !            write(*,*) titlec(16)
        !            vfh(16) = vfn(18)
        !            titleh(1,16) = '16 - '//titlehn(1,18)(6:80)
        !            titleh(2,16) = '16 - '//titlehn(2,18)(6:80)
        ! bare soil
        vfc(17) = vfn(20)
        laic(17) = lain(20)
        !            titlec(17) = title(20)
        !            do m=1,12
        !               titlem(m,17) = title12(m,20)
        !            enddo
        !            vfh(17) = vfn(20)
        !            hm(17) = hmn(20)
        !            hsd(17) = hsdn(20)
        !            titleh(1,17) = titlehn(1,20)
        !            titleh(2,17) = titlehn(2,20)
        N_VEG = 16
        N_BARE = 17
      
        !            titlefoo = 'vfn16'
        !            write(92) titlefoo, vfn(16)
        !            titlefoo = 'vfn17'
        !            write(92) titlefoo, vfn(17)
        !            titlefoo = 'Crops 15 after combining C3 and C4'
        !            write(92) titlefoo, vfc(15)
     
#else
        !crops
        vfc(15:17) = vfn(16:18)
        laic(15:17) = lain(16:18)

        !            titlec(15:17) = title(16:18)
        !            vfh(15:17) = vfn(16:18)
        !            hm(15:17) = hmn(16:18)
        !            titleh(1,15:17) = titlehn(1,16:18)
        !            titleh(2,15:17) = titlehn(2,16:18)
        ! bare soil
        vfc(18) = vfn(20)
        laic(18) = lain(20)
        !            titlec(18) = title(20)
        !            titlem(:,18) = title12(:,20)
        !            vfh(18) = vfh(20)
        !            hm(18) = hmn(20)
        !            hsd(18) = hsd(20)
        !            titleh(1,18) = titlehn(1,20)
        !            titleh(2,18) = titlehn(2,20)
        N_VEG = 17
        N_BARE = 18
#endif
      
#if 0
! TODO: Where does write(3) go to anyway?
       do k=1,N_BARE
          write(3) titlec(k), vfc(k)
       enddo
       do k=1,N_BARE
          write(3) titlec(k), laic(k)
       enddo
#endif
      
        ! check if "bare" soil is not bare
        vf_xx = 0.
        if( vfc(N_BARE) > .01 .and. laic(N_BARE) > .5 ) then
            vf_xx = vfc(N_BARE)
            lai_xx = laic(N_BARE)
        endif
#if 0
! TODO: Where does write(4) go to anyway?
        write(4) title_xx, vf_xx
        write(4) title_xx, lai_xx
#endif
      
        vf_yy = 0.
        if( vfc(10) > .1 .and. laic(10) < .5 ) then
            vf_yy = vfc(10)
            lai_yy = laic(10)
        endif
        !            write(4) title_yy, vf_yy
        !            write(4) title_yy, lai_yy
      
        vf_yy = 0.
        lai_yy = 0.
        if( vfc(N_BARE) > .1 .and. laic(10) < .01  &
            .and. laic(9) < .01 .and. laic(11) < .01  &
            .and. laic(12) < .01 .and. laic(13) < .01 ) &
        then

            vf_yy = vfc(N_BARE)
            lai_yy = laic(N_BARE)
        end if
      
        !            write(4) title_yy, vf_yy
        !            write(4) title_yy, lai_yy

        !!!! do conversions !!!!
      
        ! convert sparse veg to cold adapted shrub 9 if present
        s = sum(vfc(1:N_BARE))
        if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
        !               write(*,*) 'ERROR orig:  max and monthly lc different'
        !     &               ,s,sum(vfm(1,1:N_BARE))
        !               write(*,*) vfc(1:N_BARE)
        !               write(*,*) vfm(j,1:N_BARE)
        endif

        if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15 &
           .and. laic(N_BARE) > .0 &
           .and. vfc(9) > .0 ) &
        then
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(9), laic(9), laic(9) )
        !     lai >= lai(9)

        !              call convert_vfh(
        !    &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
        !    &              vfh(9),hm(9),hsd(9), vfc(9))
        endif
      
        s = sum(vfc(1:N_BARE))

        ! convert sparse veg to arid adapted shrub 10 if present
        if( vfc(N_BARE) > .0 .and. vfc(N_BARE) < .15 &
           .and. laic(N_BARE) > .0 &
           .and. vfc(10) > .0 ) &
        then
         
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(10), laic(10), laic(10) )
                          ! lai >= lai(10)
            !               do m=1,12
            !                  call convert_vfm(vfm(m,N_BARE),laim(m,N_BARE),
            !     &                 vfm(m,10),laim(m,10),vfc(10))
            !               enddo
                      
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(10),hm(10),hsd(10), vfc(10))
        end if
         
        ! convert the rest of sparse veg to crop 15 if present
        if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 &
              .and. vfc(15) > .0 ) &
        then
            !             print *, 'Converting spare to crop/bare',i,j,
            !     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), vfc(i,j,15)
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(15), laic(15), laic(15))
            !             print *, 'After conversion:            ',i,j,
            !     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), 
            !     &            vfc(i,j,15), laic(i,j,15)
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(15),hm(15),hsd(15), vfc(15))
        end if
      
        ! convert the rest of sparse veg to pft with biggest fraction
        ! (if present)
        if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
            maxpft = maxloc( vfc(1:16), 1 )
            !               print *, "max pft is ",maxpft
            if ( vfc(maxpft) < .0001 ) cycle
            
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(maxpft), laic(maxpft), laic(maxpft))
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(maxpft),hm(maxpft),hsd(maxpft), 
            !     &              vfc(maxpft))
                                 
        end if

        ! convert the rest of sparse veg to arid adapted shrub 10
        if( vfc(N_BARE) > .0 .and. laic(N_BARE) > .0 ) then
            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                vfc(10), laic(10), .0 )

            !               if (vfc(10) > 0.) then
            !                  hm(10) = 2.0  !Check simard.f Set_shrub_height for value!
            !                  hsd(10) = 0.
            !               endif
            !               call convert_vfh(
            !     &              vfh(N_BARE),hm(N_BARE),hsd(N_BARE),
            !     &              vfh(10),hm(10),hsd(10), vfc(10))
                                 
        end if

#ifdef SPLIT_BARE_SOIL
        call split_bare_soil(N_VEG,KM,N_BARE &
            ,bs_brightratio,vfc,laic, &
            res_out)
#endif
      
        ! check titles
        !            write(*,*) 'titlec:'
        !            do k=1,N_BARE
        !               write(*,*) trim(titlec(k))
        !            enddo
        !            write(*,*) 'titlem:'

        do k=1,18
            if (laic(k).le.0.) then
                laic(k) = undef
            end if
        end do

        ! correct height=undef when land cover>0            
        do k=1,18
            if (vfc(k).gt.0.and.laic(k).eq.undef) then
                laic(k) = 0.
            else
                laic(k) = laic(k)
            end if
        end do

        do k=1,18
            io_laiout(k)%buf(ic,jc) = laic(k)
            !laicnc(i,j,k) = laic(k)
        end do

        ! checksum lc & laimax
        io_sum_lc%buf(ic,jc) = 0.
        do k=1,18
            io_sum_lc%buf(ic,jc) = io_sum_lc%buf(ic,jc) + laic(k)
        end do
    end do  ! ic
    end do  ! jc

    call chunker%write_chunks

end do ! ichunk
end do ! jchunk

end program convert
