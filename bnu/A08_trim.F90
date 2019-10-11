module a08_mod
!  Performs several steps toward producing files suitable for input the GISS GCM
!  ModelE:
!   trimmed:  trim out small subgrid fractions, and conserve total grid LAI.
!   trimmed_scaled:  recale subgrid cover fractions to ensure they sum to 1 over
!   land.
!   trimmed_scaled_nocrops:  scale out crop cover for version with all natural
!   vegetation cover, used for rescaling with historical crop cover change.
!   trimmed_scaled_crops_ext:  crops LAI and height extended a few grid cells
!   out to ensure availability of non-zero values for crops when cover changes.

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
    use hntr_mod

implicit none

type OutputSegment_t
    character*(64) :: step
    type(Chunker_t) :: chunker

    type(ChunkIO_t), allocatable :: io_ann_lc(:,:)          ! (esub%ncover,1)
    type(ChunkIO_t), allocatable :: io_ann_lai(:,:)        ! (esub%ncover,1)
    type(ChunkIO_t), allocatable :: io_ann_hgt(:,:)        ! (esub%ncover,1)
    type(ChunkIO_t), allocatable :: io_mon_lai(:,:)   ! (esub%ncover,NMONTH)

    type(ChunkIO_t) :: io_ann_lc_checksum
    type(ChunkIO_t) :: io_ann_lclaimax_checksum
    type(ChunkIO_t) :: io_ann_lchgt_checksum
    type(ChunkIO_t), allocatable :: io_mon_lclai_checksum(:)   ! (NMONTH)
contains
    procedure :: open => outputsegment_open
    procedure :: checksum => outputsegment_checksum
end type OutputSegment_t

CONTAINS

! Opens one of trimmed, trimmed_scaled or nocrops outputs
subroutine outputsegment_open(this, step, esub)
    class(OutputSegment_t) :: this    ! Set of file handles to open
    character*(*) :: step    ! 'trimmed', 'trimmed_scaled' or 'nocrops'
    type(EntSet_t), intent(IN) :: esub
    ! --------------------------- Locals
    integer :: m,k
    type(FileInfo_t) :: info

    this%step = step
    call this%chunker%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 1, 500, 30, (/1,1/))

    ! ------- Allocate file handles
    allocate(this%io_ann_lc(esub%ncover,1))
    allocate(this%io_ann_lai(esub%ncover,1))
    allocate(this%io_ann_hgt(esub%ncover,1))
    allocate(this%io_mon_lai(esub%ncover,NMONTH))

    ! Open the files
    call this%chunker%nc_create_set( &
        esub, this%io_ann_lc(:,1), repeat_weights(esub%ncover, this%chunker%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'lc', 2004, step, '1.1')

    call this%chunker%nc_create_set( &
        esub, this%io_ann_lai(:,1), repeat_weights(esub%ncover, this%chunker%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'laimax', 2004, step, '1.1')

    call this%chunker%nc_create_set( &
        esub, this%io_ann_hgt(:,1), repeat_weights(esub%ncover, this%chunker%wta1, 1d0, 0d0), &
        LAI_SOURCE, 'M', 'hgt', 2004, step, '1.1')

    do m=1,NMONTH
        call this%chunker%nc_create_set( &
            esub, this%io_mon_lai(:,m), repeat_weights(esub%ncover, this%chunker%wta1, 1d0, 0d0), &
            LAI_SOURCE, 'M', 'lai', 2004, step, '1.1', &
            doytype='month', idoy=m)
    end do

    ! ------------- Open checksum files
    allocate(this%io_mon_lclai_checksum(NMONTH))

    ! Open checksum files
    call this%chunker%file_info(info, esub, LAI_SOURCE, 'M', 'lc', 2004, step, '1.1', &
        varsuffix = '_checksum')
    call this%chunker%nc_create(this%io_ann_lc_checksum, &
        weighting(this%chunker%wta1,1d0,0d0), &
        info%dir, info%leaf, info%vname, &
        info%long_name, info%units)

    call this%chunker%file_info(info, esub, LAI_SOURCE, 'M', 'lclaimax', 2004, step, '1.1', &
        varsuffix = '_checksum')
    call this%chunker%nc_create(this%io_ann_lclaimax_checksum, &
        weighting(this%io_ann_lc_checksum%buf,1d0,0d0), &
        info%dir, info%leaf, info%vname, &
        info%long_name, info%units)


    call this%chunker%file_info(info, esub, LAI_SOURCE, 'M', 'lchgt', 2004, step, '1.1', &
        varsuffix = '_checksum')
    call this%chunker%nc_create(this%io_ann_lchgt_checksum, &
        weighting(this%io_ann_lc_checksum%buf,1d0,0d0), &
        info%dir, info%leaf, info%vname, &
        info%long_name, info%units)

    do m=1,NMONTH
        call this%chunker%file_info(info, esub, LAI_SOURCE, 'M', 'lclai', 2004, step, '1.1', &
            doytype='month', idoy=m, varsuffix='_checksum')
        call this%chunker%nc_create(this%io_mon_lclai_checksum(m), &
            weighting(this%io_ann_lc_checksum%buf,1d0,0d0), &
            info%dir, info%leaf, info%vname, &
            info%long_name, info%units)
    end do


end subroutine outputsegment_open

subroutine outputsegment_checksum(this, esub)
    class(OutputSegment_t) :: this    ! Set of file handles to open
    type(EntSet_t), intent(IN) :: esub
    ! ---------- Local Vars
    integer :: ic,jc,k,m
    real*4 :: sum

    do jc = 1,this%chunker%chunk_size(2)
    do ic = 1,this%chunker%chunk_size(1)
        ! lc_checksum = sum(lc)
        sum = 0
        do k=1,esub%ncover
            sum = sum + this%io_ann_lc(k,1)%buf(ic,jc)
        end do
        this%io_ann_lc_checksum%buf(ic,jc) = sum

        ! lclaimax_checksum = sum(LC*LAIMAX)
        sum = 0
        do k=1,esub%ncover
            sum = sum + this%io_ann_lc(k,1)%buf(ic,jc) * this%io_ann_lai(k,1)%buf(ic,jc)
        end do
        this%io_ann_lclaimax_checksum%buf(ic,jc) = sum

        ! lchgt_checksum = sum(LC*HGT)
        sum = 0
        do k=1,esub%ncover
            sum = sum + this%io_ann_lc(k,1)%buf(ic,jc) * this%io_ann_hgt(k,1)%buf(ic,jc)
        end do
        this%io_ann_lchgt_checksum%buf(ic,jc) = sum

        ! MONTHLY: lclai_checksum = sum(LC*LAI)
        do m=1,NMONTH
            sum = 0
            do k=1,esub%ncover
                sum = sum + this%io_ann_lc(k,1)%buf(ic,jc) * this%io_mon_lai(k,m)%buf(ic,jc)
            end do
            this%io_mon_lclai_checksum(m)%buf(ic,jc) = sum
        end do
    end do
    end do


end subroutine

! =======================================================================


!Converts LAI that is on BARE to a vegetated fraction.
!Reduces vf1, increases vf2 by increment that maintains lai2,
! or only increases lai2 if lai1 is so big that the cover-weighted
! lai is bigger than lai2.
!vf1, lai1 - BARE cover fraction and LAI
!vf2, lai2 - vegetated cover fraction and max LAI
subroutine convert_vf(vf1, lai1, vf2, lai2, laimin)
    real*4, intent(INOUT) :: vf1, lai1, vf2, lai2
    real*4, intent(IN) :: laimin
    ! --- Locals
    real*4 tot_la, tot_vf, new_lai, new_vf

    tot_la = vf1*lai1 + vf2*lai2
    tot_vf = vf1 + vf2
    new_lai = tot_la/tot_vf
    !Fix for if 100% conversion to vf2.
    if (new_lai.eq.0.0) then
       new_vf = vf2 !0.0 !Keep original cover
    else
       new_lai = max( new_lai, laimin)
       new_vf = tot_la/new_lai
    endif
    ! get rid of round-off errors
    !!new_vf = min( new_vf, tot_vf )

    vf2 = new_vf
    lai2 = new_lai
    vf1 = max(0.,tot_vf - vf2) !Argh, round-off errors
    lai1 = 0.
end subroutine convert_vf

!convert monthly vf using vfc for new vf2 derived from laimax trim. 
!vf1, lai1 - BARE cover fraction and LAI
!vf2, lai2 - vegetated cover fraction and monthly LAI
!vfc = vf2 + dvf1 from laimax trim (dvf1 is positive)
!vf1 + vf2 = vfc + (vf1-dvf1)
!vf1*lai1 = 0*(vf1-dvf1) + (dvf1*dlai1)
!new_lai2 = ((dvf1*dlai1) + (vf2*lai2))/(dvf1+vf2)
!         = ((vf1*lai1) + (vf2*lai2))/vfc
subroutine convert_vfm(vf1, lai1, vf2, lai2, vfc)
    real*4, intent(inout) :: vf1, lai1, vf2, lai2
    real*4, intent(in) :: vfc
    !---
    real*4 :: tot_lai, new_lai, new_vf

    new_vf = vfc
    if (new_vf.le.0.) then  
       write(ERROR_UNIT,*) 'new_vf is .leq. zero',vf1,vf2,lai1,lai2,vfc
       STOP
    endif
    new_lai = (vf1*lai1 + vf2*lai2)/new_vf
    lai2 = new_lai
    vf1 = max(0., vf1 - (vfc - vf2))   !bare, max for neg. round-off error
    vf2 = new_vf
    lai1 = 0.
end subroutine convert_vfm


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
subroutine convert_vfh(vf1, h1, hsd1, vf2, h2, hsd2, vfc)
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
    vf1 = max(0., vf1 - (vfc - vf2)) !bare, max for neg. round-off error
    h1 = 0.
    hsd1 = 0.

end subroutine convert_vfh


!Split BARE soil into BRIGHT and DARK cover to preserve albedo from
!  "old" ModelE cover.  Should be called after each trim, scale, nocrops.
!Any LAI on BARE soil should already have been moved to vegetated cover,
!  so laic(:,:,N_BARE) should be zero.
!!!!This checks for cases if BARE is original total or was previously split.
! This assumes it's NOT the first time splitting bare soil (that was done in a previous step)
subroutine do_split_bare_soil( &
    esub, bs_brightratio, &
    vfc,laic,hm,hsd, &
    vfm,laim)
implicit none
    type(GcmEntSet_t), intent(IN) :: esub
    real*4, intent(IN) :: bs_brightratio   !Fraction of bare that is bright.
    real*4, intent(INOUT), dimension(esub%ncover) :: vfc,laic,hm,hsd
    real*4, intent(INOUT), dimension(esub%ncover,NMONTH) :: vfm, laim
    !-----Local----
    real*4 :: vft
    integer :: m

    vft = vfc(esub%bare_bright) + vfc(esub%bare_dark)
    if ((vft.lt.0.and.abs(vft)<1e20).or. &
          ((vft.gt.0.and.abs(vft)<1e20).and.(bs_brightratio.lt.0.))) then !Bad data
        write(ERROR_UNIT,*) 'ERR2: bs_brightratio',vft,bs_brightratio
        !Keep existing bright and dark
    else             !Good data
        vfc(esub%bare_bright) = vft*bs_brightratio
        vfc(esub%bare_dark) = vft*(1.-bs_brightratio)
        laic(esub%bare_bright) = 0.
        laic(esub%bare_dark) = 0.
    end if

    do m=1,NMONTH
        if ((vft<0.and.abs(vft)<1e20).or. &
           ((vft>0.and.abs(vft)<1e20).and.(bs_brightratio<0.))) then !Bad data
            write(ERROR_UNIT,*) 'ERR2m: bs_brightratio',vft,abs(vft) &
              ,vfc(esub%bare_bright),vfc(esub%bare_dark),bs_brightratio
            !Keep existing bright and dark, or check why vft<0.
        else          !Good data
            vfm(esub%bare_bright,m) = vft*bs_brightratio
            vfm(esub%bare_dark,m) = vft*(1d0-bs_brightratio)
        end if
        laim(esub%bare_bright,m) = 0.
        laim(esub%bare_dark,m) = 0.
    end do

    hm(esub%bare_dark) = 0.
    hsd(esub%bare_dark) = 0.

end subroutine do_split_bare_soil

subroutine replace_crops(esub, IMn,JMn,i,j, &
         bs_brightratio, &
         vfc,vfm,laic,laim, &
         naturalvegfound)
    !Replace crops in (i,j) with main natural cover in cell or adjacent.
    !The check for existence of crops in (i,j) is done before subroutine call
    !Fix 6 cells that are all crops due to MODIS error or two islands.
    !Assign LAI (cover-avg) of natural cover if from adjacent cells.
    !Do not zero crop LAI but keep for historical cover change.
    !NOTES:  Checks of adjacent cells will use replaced cover from previously
    ! processed cells (previous i,j).  I don't see this as a problem,
    ! since there won't be a consistent dominance either side of the cell.
    ! it helps reduce iteration for adjacent cells.
    implicit none
    type(GcmEntSet_t), intent(IN) :: esub
    integer, intent(in) :: IMn,JMn,i, j
    real*4, intent(in) :: bs_brightratio(IMn,JMn) ! Fraction of bare that is bright.
    real*4 :: vfc(IMn,JMn,esub%ncover) , laic(IMn,JMn,esub%ncover)
    real*4 :: vfm(NMONTH,IMn,JMn,esub%ncover), laim(NMONTH,IMn,JMn,esub%ncover)
    !integer :: flag
    logical,intent(out) :: naturalvegfound
    !----Local----
    integer :: KM
    integer :: m, k, ii, jj !cell to search for natural (i,j or adjacent)
    real*4 :: covmax
    integer :: covmaxk, covmaxii,covmaxjj
    real*4, allocatable :: covsum(:), covavglai(:)   ! (NATVEG=14)
    real*4, allocatable :: covsumm(:,:),covavglaim(:,:)  ! (NMONTH,NATVEG)
    integer :: dg
    real*4 :: br, brcov   !bright soil ratio, br, to total bare soil, brcov

    integer :: NATVEG    ! Index of highest natural vegetation type
    integer :: c3herb_s, c4herb_s, c3_grass_arct_s, c3_grass_s

    KM = esub%ncover
    c3herb_s = esub%svm(CROPS_C3_HERB)   ! shortcut
    c4herb_s = esub%svm(CROPS_C4_HERB)   ! shortcut
    c3_grass_arct_s = esub%svm(C3_GRASS_ARCT)

    NATVEG = esub%svm(CROPS_C3_HERB)-1
    allocate(covsum(NATVEG))
    allocate(covavglai(NATVEG))
    allocate(covsumm(NMONTH, NATVEG))
    allocate(covavglaim(NMONTH, NATVEG))

    if (sum(vfc(i,j,c3herb_s:c4herb_s)).le.0.0) then
       return !No crops, no need to do
    endif
    
    !--First check just (i,j) cell for natural veg---
    naturalvegfound = .false.
    covmax = 0.d0
    covmaxk = 0
    covmaxii = i
    covmaxjj = j
    do k=1,NATVEG         !Find max non-crop, non-bare natural cover type
       if (vfc(i,j,k) > covmax)  then
          covmax = vfc(i,j,k)
          covmaxk = k
       endif
    enddo
    if (covmax.gt.0.d0) then  !Assign dominant natural veg to crop
       ! -------------- Replace crops with dominant natural vegetation type
       naturalvegfound = .true.
       vfc(i,j,covmaxk) = vfc(i,j,covmaxk)  &
            + vfc(i,j,c3herb_s)  + vfc(i,j,c4herb_s)
       vfc(i,j,c3herb_s:c4herb_s) = 0.0 !zero out crop cover - done below
       vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk) &
            + vfm(:,i,j,c3herb_s)  + vfm(:,i,j,c4herb_s)
       vfm(:,i,j,c3herb_s:c4herb_s) = 0.0 !zero out crop cover - done below
       !If in same cell, keep original LAI in (i,j) of natural veg
       !!DON'T DO ASSIGNMENT BELOW
       !!laic(i,j,covmaxk) = laic(covmaxii,covmaxjj,covmaxk)
       !!laim(:,i,j,covmaxk) = laim(:,covmaxii,covmaxjj,covmaxk)
       !write(*,*) 'Cell has crops + natural',i,j,s
    else
       ! --------- If no natural vegetation, look at adjacent cells for veg type
       write(*,*) "No natural veg in (i,j). Checking adjacent",i,j !,s
       dg = 1
       do while (.not.naturalvegfound) !Check up to 9 grid 

          if (dg > 5) then
              ! ---------------- NO NATURAL VEG FOUND IN NEARBY GRIDCELLS!
              ! Look again at immediately adjacent gridcells, for bare soil vs. crops

              ! ===> Assign crop LC and LAI to C3 grass (which started out with lc=0)
              c3_grass_s = esub%svm(C3_GRASS_PER)

              write(*,*) "Assigning to C3_GRASS_PER at",i,j

              ! (this code block copied from below)
              vfc(i,j,c3_grass_s) = vfc(i,j,c3_grass_s)  &
                   + vfc(i,j,c3herb_s) + vfc(i,j,c4herb_s)
              vfc(i,j,c3herb_s:c4herb_s) = 0.0 !zero out crop cover - done below
              vfm(:,i,j,c3_grass_s) = vfm(:,i,j,c3_grass_s) &
                   + vfm(:,i,j,c3herb_s) + vfm(:,i,j,c4herb_s)
              vfm(:,i,j,c3herb_s:c4herb_s) = 0.0 !zero out crop cover  - done below

              ! Move the crop LAI to the new C3 grass
              laic(i,j,c3_grass_s) = &
                  (laic(i,j,c3herb_s)*vfc(i,j,c3herb_s) + laic(i,j,c4herb_s)*vfc(i,j,c4herb_s)) &
                  / (vfc(i,j,c3herb_s) + vfc(i,j,c4herb_s))
              laim(:,i,j,c3_grass_s) = &
                  (laim(:,i,j,c3herb_s)*vfm(:,i,j,c3herb_s) + laim(:,i,j,c4herb_s)*vfm(:,i,j,c4herb_s)) &
                  / (vfm(:,i,j,c3herb_s) + vfm(:,i,j,c4herb_s))

              EXIT   ! do loop; done increasing dg
          end if


          covsum(:) = 0.d0
          covavglai(:) = 0.d0
          covsumm(:,:) = 0.
          covavglaim(:,:) = 0.
          do ii=max(1,i-dg),min(i+dg,IMn)
             do jj=max(1,j-dg),min(j+dg,JMn)
                if ( (ii.ge.1).and.(ii.le.IMn) &
                     .and.(jj.ge.1).and.(jj.le.JMn) & !in grid range
                     .and.((ii.ne.i).or.(jj.ne.j)) ) & !not the i,j center cell
                     then
                   do k=1,NATVEG
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
          covmax = 0.
          covmaxk = 0
          covmaxii = 0
          covmaxjj = 0
          do k=1,NATVEG           !Find largest adjacent natural cover
             if (covsum(k)>covmax) then
                covmax = covsum(k)
                covmaxk = k
             endif
          enddo
          if (covmax>0.) then
             do m=1,NMONTH  !Error check
                if (covsumm(m,covmaxk)<0.) then
                   write(*,*) 'ERRc ',i,j,ii,jj,m &
                        ,covmax,covsum(covmaxk), covsumm(m,covmaxk) &
                        ,covavglai(covmaxk),covavglaim(m,covmaxk)
                   STOP
                endif
             enddo
             ! Assign adjacent natural cover type and LAI to crop
             naturalvegfound = .true.
             !covavglai(k) = covavglai(k)/covsum(k)  !BUG FOUND
             covavglai(covmaxk) = covavglai(covmaxk)/covsum(covmaxk)
             covavglaim(:,covmaxk) = covavglaim(:,covmaxk)/ &
                  covsumm(:,covmaxk)
             vfc(i,j,covmaxk) = vfc(i,j,covmaxk)  &
                  + vfc(i,j,c3herb_s) + vfc(i,j,c4herb_s)
             vfc(i,j,c3herb_s:c4herb_s) = 0.0 !zero out crop cover - done below
             vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk) &
                  + vfm(:,i,j,c3herb_s) + vfm(:,i,j,c4herb_s)
             vfm(:,i,j,c3herb_s:c4herb_s) = 0.0 !zero out crop cover  - done below
             !Assign LAI in (i,j) from adjacent cell
             laic(i,j,covmaxk) = covavglai(covmaxk)
             laim(:,i,j,covmaxk) = covavglaim(:,covmaxk)
             !write(*,*) 'Cell or adjacent has crops + natural',i,j,s
             write(*,*) 'Found natural veg in adjacent cells' &
                  ,i,j,covmaxk,dg
             do m=1,NMONTH !Error check
                if (covavglaim(m,covmaxk)>10.) then
                   write(*,*) 'ERRc2 bad avg lai',i,j,ii,jj,m,covmaxk &
                        ,covsumm(m,covmaxk),covavglaim(m,covmaxk)
                   STOP
                endif
             enddo
          else
             write(*,*) 'No natural veg in adjacent cells'
             dg=dg+1
             write(*,*) 'Increasing dg,',i,j, dg
          endif
       end do                 !while
    endif !Checking for adjacent natural veg
    
    !!!** PRESCRIBED FIXED BY GRID CELL FOR 144X90 ONLY ***!!!
    !For IMn,JMn = 144x90
    !Fix remaining all-crop cells.  These are:  
    !- MODIS error in Antarctic, 
    !     (i,j) = (34,9),(41,9),(48,9),(34,10) -> C3 arctic grass
    !- Islands:  
    !  Mauritius (98,36) - sugar cane,tea,pasture,forest,savanna -> C4 grass 
    !  Nauru (139,45) - grassland bordered by tropicalforest -> C4 grass
    !Antarctic MODIS error
!      if ((IMn.eq.144).and.(JMn.eq.90)) then
!         if ( (((i.eq.34).or.(i.eq.41).or.(i.eq.48)).and.(j.eq.9))
!     &        .or.( (i.eq.34).and.(j.eq.10) ) ) then !Antartic
!            vfc(i,j,c3_grass_arct) = vfc(i,j,c3_grass_arct) + vfc(i,j,c3herb_s) + vfc(i,j,c4herb_s)
!            vfm(:,i,j,c3_grass_arct) = vfm(:,i,j,c3_grass_arct) 
!     &           + vfm(:,i,j,c3herb_s) + vfm(:,i,j,c4herb_s)
!            laic(i,j,c3_grass_arct) = (laic(i,j,c3herb_s)*vfc(i,j,c3herb_s) 
!     &           + laic(i,j,c4herb_s)*vfc(i,j,c4herb_s)) / (vfc(i,j,c3herb_s)+vfc(i,j,c4herb_s))
!            do m=1,12
!               laim(m,i,j,c3_grass_arct) = ( laim(m,i,j,c3herb_s)*vfm(m,i,j,c3herb_s) +
!     &              laim(m,i,j,c4herb_s)*vfm(m,i,j,c4herb_s) ) / 
!     &              ( vfm(m,i,j,c3herb_s)+vfm(m,i,j,c4herb_s) )
!            enddo
!            write(*,*) 'Replaced Antarctic crops.'
!            naturalvegfound=.true.
!         elseif ( ((i.eq.98).and.(j.eq.36))
!     &           .or.((i.eq.139).and.(j.eq.45)) ) then !tropical islands
!            vfc(i,j,12) = vfc(i,j,12) + vfc(i,j,c3herb_s) + vfc(i,j,c4herb_s)
!            laic(i,j,12) = (laic(i,j,c3herb_s)*vfc(i,j,c3herb_s) 
!     &        + laic(i,j,c4herb_s)*vfc(i,j,c4herb_s)) / (vfc(i,j,c3herb_s)+vfc(i,j,c4herb_s))
!            do m=1,12
!               laim(m,i,j,12) = ( laim(m,i,j,c3herb_s)*vfm(m,i,j,c3herb_s) +
!     &              laim(m,i,j,c4herb_s)*vfm(m,i,j,c4herb_s)) / 
!     &              ( vfm(m,i,j,c3herb_s)+vfm(m,i,j,c4herb_s) )
!            enddo
!            write(*,*) 'Replaced island crops.'
!            naturalvegfound=.true.
!         endif
!      endif

    !*** PRESCRIBED GRID CELL FIXES FOR 720X360 *********
    if (.not.naturalvegfound) then !Assign bare
       if (bs_brightratio(i,j).eq.FillValue) then
          !Find bs_brightratio in adjacent cells
          br = 0.
          brcov = 0.
          dg=0
          do while ((brcov.eq.0.).and.(dg.le.5)) !arbitrary neighbor range
             dg=dg+1
             do ii=max(1,i-dg),min(i+dg,IMn)
                do jj=max(1,j-dg),min(j-dg,JMn)
                   if (bs_brightratio(ii,jj).gt.0.) then
                      br = br + bs_brightratio(ii,jj)
                      brcov = brcov + sum(vfc(ii,jj,esub%bare_bright:esub%bare_dark))
                   endif
                enddo
             enddo
          enddo
          if (brcov.gt.0.) then
             br = br/brcov    !cover-weighted average br
             write(*,*) "Adjacent bare soil found",i,j,dg, br
          else
             write(*,*) "Adjacent bare soil NOT found",i,j,dg, br
             if (.FALSE.) then !Done only to diagnose cells w/o bratio
                write(*,*) "assigning gray 0.2"
                br=0.4        !This gives 0.5*.4 + 0.*0. = 0.2 albedo
             else
                write(*,*) "Assign specific grid cells"
                br=FillValue     !Assign FillValue to make sure no wrong values.
             endif
          endif
       else !bs_brightratio was defined
          br = bs_brightratio(i,j)
       endif

       !Assign bare soil fractions
       if (br.ne.FillValue) then  !Adjacent bs_brightratio found
          vfc(i,j,esub%bare_bright) = vfc(i,j,esub%bare_bright) + &
               sum(vfc(i,j,c3herb_s:c4herb_s))*br
          vfc(i,j,esub%bare_dark) = vfc(i,j,esub%bare_dark) + &
               sum(vfc(i,j,c3herb_s:c4herb_s))*(1. - br)
          vfc(i,j,c3herb_s:c4herb_s) = 0. !zero out crops
          vfm(:,i,j,esub%bare_bright) = vfm(:,i,j,esub%bare_bright) + &
               sum(vfm(:,i,j,c3herb_s:c4herb_s))*br
          vfm(:,i,j,esub%bare_dark) = vfm(:,i,j,esub%bare_dark) + &
               sum(vfm(:,i,j,c3herb_s:c4herb_s))*(1. - br)
          vfm(:,i,j,c3herb_s:c4herb_s) = 0. !zero out crops
          write(*,*) 'Assigned nearby bare' &
               ,i,j, dg, bs_brightratio(i,j),br
       elseif ((IMn.eq.720).and.(JMn.eq.360)) then
    !Adjacent bs_brightratio not found for only these cells
!     Replace specific grid cells ONLY VALID FOR IM=720,JM=360!!!
!     NOTE:  "Mauritius" and "Fiji/Amer.Samoa" might be off-grid and land mask.
!     assigning bare Found for IM=720,JM=360 = (lon, lat)
!     Below:  i,j,bs_brightratio(i,j), br = (lon, lat)
!     Antarctica 1: (pureesub%bare_bright has C3 crop, partial coastal cell/ice)
!     c4herb_s9  36  -1.00000002E+30  0.400000006 = (-95.75,-72.25)
!     c4herb_s8  37  -1.00000002E+30  0.400000006 = (-96.25,-71.75)
!     Antarctica 2: (pureesub%bare_bright has C3 crop, partial coastal cell/ice)
!     239  36  -1.00000002E+30  0.400000006 = (-60.75,-72.75)
!     204  36  -1.00000002E+30  0.400000006 = (-78.25,-72.25)
!     Mauritius (too east of):   (pure17 was everbroad, C4 crop, sparse, coast)
!     487 141 0.461078912 0.461078912  = (63.25, -19.75)
!     Fiji?  American Samoa? (pure17 was everbroad, C4 crop, sparse, coast)
!     694 179 0.353166729 0.353166729  = (163.25,-0.75)
          if (((i.eq.169).and.(j.eq.36)) & !Antarctica1
               .or.((i.eq.168).and.(j.eq.37)) & !Antarctica1
               .or.((i.eq.239).and.(j.eq.36)) & !Antarctica2
               .or.((i.eq.204).and.(j.eq.36))) & !Antarctica2
               then
                !* Assign Antarctica cells to arctic grass
             vfc(i,j,c3_grass_arct) = vfc(i,j,c3_grass_arct) + sum(vfc(i,j,c3herb_s:c4herb_s))
             vfc(i,j,c3herb_s:c4herb_s) = 0. !zero out crops
             laic(i,j,c3_grass_arct) = max(laic(i,j,c3_grass_arct), &
                  max(laic(i,j,c3herb_s),laic(i,j,c4herb_s)))
             do m=1,NMONTH
                vfm(m,i,j,c3_grass_arct_s)=vfm(m,i,j,c3_grass_arct) + sum(vfm(m,i,j,c3herb_s:c4herb_s))
                vfm(m,i,j,c3herb_s:c4herb_s) = 0. !zero out crops
                laim(m,i,j,c3_grass_arct)=max(laim(m,i,j,c3_grass_arct), &
                     max(laim(m,i,j,c3herb_s),laim(m,i,j,c4herb_s)))
             enddo
             write(*,*) 'Assigned c3 crops cells in Antartica',i,j
          endif
          if (((i.eq.487).and.(j.eq.141)) & !Mauritius
               .or.((i.eq.694).and.(j.eq.179))) & !Fiji or American Samoa
               then
             !* Assign tropical rainforest
             vfc(i,j,2) = vfc(i,j,2) + sum(vfc(i,j,c3herb_s:c4herb_s))
             vfc(i,j,c3herb_s:c4herb_s) = 0. !zero out crops
             laic(i,j,2) = max(laic(i,j,2), &
                  max(laic(i,j,c3herb_s), laic(i,j,c4herb_s)))
             do m=1,NMONTH
                vfm(m,i,j,2)=vfm(m,i,j,2) + sum(vfm(m,i,j,c3herb_s:c4herb_s))
                vfm(m,i,j,c3herb_s:c4herb_s) = 0. !zero out crops
                laim(m,i,j,2)=max(laim(m,i,j,2), &
                     max(laim(m,i,j,c3herb_s),laim(m,i,j,c4herb_s)))
             enddo
             write(*,*) 'Assigned c4 crops cells in islands',i,j
          endif
       else
          !Don't do any grid cell fixes
          write(*,*) 'WRONG GRID RES:  PRESCRIBED FIX NOT DONE: '
          write(*,*) '**Check array indices for grid res for nocrops'
       endif                  !Adjacent bs_brightratio not found
    endif                     !.not.(naturalvegfound)
    
    !Zero out crop cover - done above for each vfc assignment.
!      vfc(i,j,c3herb_s:c4herb_s) = 0.
!      vfm(:,i,j,c3herb_s:c4herb_s) = 0.
end subroutine replace_crops


subroutine fill_crops(IMn,JMn,io_bs, vfc15,laic15, &
    vfm15,laim15,hm15,hsd15,&
    laiccrop, laimcrop, hmcrop,hsdcrop)

    !This performs in-fill once for herb crop (PFT15) LAI to create
    !an extended crop LAI data set for use with historical crop cover.
    integer, intent(IN) :: IMn, JMn
    type(ChunkIO_t), target, intent(IN) :: io_bs
    real*4, intent(IN) :: vfc15(:,:), laic15(:,:)  !(i,j)
    real*4, intent(IN) :: vfm15(:,:,:)   ! (m,i,j)
    real*4, intent(IN) ::  laim15(:,:,:) !(m,i,j)
    real*4, intent(IN) :: hm15(:,:) !(i,j)
    real*4, intent(IN) :: hsd15(:,:) !(i,j)
    ! ------ OUTPUTS
    real*4 :: laiccrop(:,:), laimcrop(:,:,:) !OUTPUT, max monthly
    real*4 :: hmcrop(:,:)
    real*4 :: hsdcrop(:,:)
    !--Local----
    integer :: i, j, m, ii, jj, nocropcells, dg
    real*4 :: covsum15, laiavg15, covmsum15(NMONTH), laimavg15(NMONTH)
    real*4 :: hmavg15, hsdavg15
    character*80 :: titlefoo

    titlefoo = 'vfc15 in - crop fill'
    write(93) titlefoo, vfc15
    titlefoo = 'laic15 in - crop fill'
    write(93) titlefoo, laic15

    nocropcells=0
    do i=1,IMn
       do j=1,JMn

        ! Avoid cells not processed by previous steps
        if (io_bs%buf(i,j) == FillValue) cycle

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
             covmsum15(:) = 0.d0
             laimavg15(:) = 0.d0
              !do ii=i-1,i+1  !ext 1 grid cell
              !   do jj=j-1,j+1
             dg=5
! This is causing problems in trimmed_scaled_ext (but not trimmed_scaled)
! Blocky stuff in crops_herb
             do ii=i-dg,i+dg    !ext 5 grid cells = 2.5 degrees at HXH
                do jj=j-dg,j+dg
                   if ( (ii.ge.1).and.(ii.le.IMn) &
                        .and.(jj.ge.1).and.(jj.le.JMn) & !in grid range
                        .and.((ii.ne.i).or.(jj.ne.j)) ) & !not in i,j
                        then
                      if (vfc15(ii,jj).gt.0.d0) then
                         covsum15 = covsum15 + vfc15(ii,jj)
                         laiavg15 = laiavg15  &
                              + laic15(ii,jj)*vfc15(ii,jj)
                         hmavg15 = hmavg15 + hm15(ii,jj)*vfc15(ii,jj)
                         hsdavg15 = hsdavg15 + hsd15(ii,jj)*vfc15(ii,jj)
                         covmsum15(:) = covmsum15(:) + vfm15(:,ii,jj)
                         laimavg15(:) = laimavg15(:) + laim15(:,ii,jj)*vfm15(:,ii,jj)
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

subroutine check_laim(laim, msg)
    real*4 :: laim(:,:)
    character*(*), intent(IN) :: msg
    ! --------- Locals
    integer :: i,j

    do j=1,size(laim,2)
    do i=1,size(laim,1)
        if (laim(i,j) /= laim(i,j)) print *,trim(msg),i,j
    end do
    end do

end subroutine check_laim

subroutine do_part1_2_trimmed(esub, IM,JM, io_ann_lc, io_bs, io_ann_hgt, io_ann_lai, io_mon_lai,    tr, ts)
    type(GcmEntSet_t), intent(IN) :: esub
    integer, intent(IN) :: IM,JM
    type(ChunkIO_t), intent(IN) :: io_ann_lc(:), io_bs, io_ann_hgt(:,:), io_ann_lai(:,:), io_mon_lai(:,:)
    ! ------------- OUTPUT
    type(OutputSegment_t) :: tr,ts

    ! ------------- Locals
    integer :: i,j,k,m
    integer :: n_bare
    ! Annual LC and LAI (LC doesn't change so annual LC is used everywhere)
    real*4, dimension(esub%ncover) :: vfc,laic,hm,hsd,vfh
    ! Monthly LC and LAI (NOTE: Monthly LC is set to same as annual)
    real*4, dimension(esub%ncover,NMONTH) :: vfm,laim
    real*4 :: bs_brightratio
    integer :: arid_shrub_s   ! Shortcut indices

    logical :: isgood
    real*4 :: sum_vfc,sum_vfm

    ! real*4, parameter :: lc_trim_threshold = .05
    ! Wow, that <0.05 threshold sure cleans out a lot of that dark
    ! blue.  I think we can go to a threshold of < 0.1.
    real*4, parameter :: lc_trim_threshold = .1

    print *,'========================= Part 1&2: trimmed, trimmed_scaled'

    do j = 1,JM
    do i = 1,IM
        ! Avoid cells not processed by previous steps
        if (abs(io_bs%buf(i,j)) > 1.e10) cycle

        ! -------------------------- Read Inputs
        do k=1,esub%ncover
            ! Annual LC and LAI file
            vfc(k) = io_ann_lc(k)%buf(i,j)    ! vfn in 1km
            laic(k) = io_ann_lai(k,1)%buf(i,j)  ! laic in 1km

            ! Monthly LC and LAI files
            do m=1,NMONTH
                vfm(k,m) = vfc(k)   ! vfnm in 1km; Re-use annual LC
                laim(k,m) = io_mon_lai(k,m)%buf(i,j)  ! lainm in 1km
call check_laim(laim, 'check1')
            end do ! m

            ! Height file
            hm(k) = io_ann_hgt(k,1)%buf(i,j)
            ! hsd = stdev
            ! vfh = vfn (LC annual), in GISS 16 pfts format (but C3 and C4 crop are summed)
            !      ===> "LC for heights", i.e. just use the global LC
            vfh(k) = vfc(k)
            hsd(k) = 0

        end do    ! k=1,esub%ncover
        ! Bare soil brightness ratiaafo
        bs_brightratio = io_bs%buf(i,j)

        ! Make sure this gridcell has been defined in inputs
        isgood = (bs_brightratio /= FillValue)
        ! ----------------------------------------------------

        ! By definition, N_BARE indexes the last bare covertype
        !    If not yet split: = BARE_SPARSE
        !    If split already: = BARE_DARK
        ! 1:N_BARE == everything but water
        ! We assume we've already been split in A04...A06
        n_bare = esub%bare_dark

        ! ========== Part 1: trimmed

        ! ------------------------------------------------------
        ! convert arid adapted shrub with lai < .15 to bare soil
        ! and restrict lai >= .15 
        arid_shrub_s = esub%svm(ARID_SHRUB)   ! shortcut
        if( vfc(arid_shrub_s) > .0 .and. laic(arid_shrub_s) < .15 ) then

            call convert_vf(vfc(N_BARE), laic(N_BARE), &
                 vfc(arid_shrub_s), laic(arid_shrub_s), .15 )
                                ! lai >= .15
            if (isgood.and.vfc(arid_shrub_s).le.0.0) then !ERROR CHECK
               write(ERROR_UNIT,*) &
                    'b vfc10=',vfc(arid_shrub_s),vfc(N_BARE), &
                    laic(arid_shrub_s),laic(N_BARE),i,j
               STOP
            endif
            do m=1,NMONTH
               if (vfc(arid_shrub_s).gt.0.0) then
                  call convert_vfm(vfm(N_BARE,m),laim(N_BARE,m), &
                       vfm(arid_shrub_s,m),laim(arid_shrub_s,m), vfc(arid_shrub_s))
               else
                  if (isgood) then
                    write(ERROR_UNIT,*) 'vfc(arid_shrub_s) le 0.:', N_BARE,m
                    write(ERROR_UNIT,*) vfc(:)
                    write(ERROR_UNIT,*) vfm(:,m)
                    write(ERROR_UNIT,*) laic(:)
                    write(ERROR_UNIT,*) laim(:,m)
                  end if
               endif
               if ((vfm(arid_shrub_s,m).lt.0.).or.(vfm(N_BARE,m).lt.0.)) &
                    then !CHECK ERROR
                 if (isgood) then
                   write(ERROR_UNIT,*) 'vfm<0:',N_BARE,m
                   write(ERROR_UNIT,*) vfc(:)
                   write(ERROR_UNIT,*) vfm(:,m)
                   write(ERROR_UNIT,*) laic(:)
                   write(ERROR_UNIT,*) laim(:,m)
                   STOP
                end if
              endif
            enddo
call check_laim(laim, 'check1')


            ! **hm = Simard hegihts
            ! **hsd = stdev of heights
            ! **vfh = "LC for heights"
            call convert_vfh( &
                 vfh(N_BARE),hm(N_BARE),hsd(N_BARE), &
                 vfh(arid_shrub_s),hm(arid_shrub_s),hsd(arid_shrub_s), vfc(arid_shrub_s))
        end if   ! (vfc(arid_shrub_s) > .0 .and. laic(arid_shrub_s) < .15 ) then

        if (isgood.and.split_bare_soil) then
            call do_split_bare_soil(esub, bs_brightratio, &
                vfc,laic,hm,hsd,vfm,laim)
        end if


        ! ------------------ Trim small fractions, just zero them out
        ! Trimming stage is not getting rid of some tiny fractions
        ! Things over Antarctica, also over the ocean
        ! Should be getting rid of those; maybe not because there's nothing else to replace them with.
        !
        !   At 1/2 degree: if lc<.1 and laimax==0 or undef  ==> zero it out
        ! if laimax is nonzero, then... must make sure some lai of that gridcell is preserved when zero out that cover
        ! if there is vegetation in some other PFT in the same cell...
        !      ==>  zero it out
        ! If no other vegetation in that cell and that's the only point with LAI and it's that small...
        !     ==> just zero it out
        !     Then checksum to see loss/gain in LAI

        do k=1,esub%ncover
            ! Note: Do not trim bare_bright and bare_dark.  They still have to
            ! add up to total bare soil fraction and give the right
            ! weighted-sum albedo.  Or you can trim them ONLY if bare_bright +
            ! bare_dark < 0.05.
            if ((k == esub%bare_dark).or.(k == esub%bare_bright)) cycle

            if (vfc(k) < lc_trim_threshold) vfc(k) = 0
            do m=1,NMONTH
                if (vfm(k,m) < lc_trim_threshold) vfm(k,m) = 0
            end do   ! m
        end do    ! k=1,esub%ncover

        do k=1,esub%ncover
            ! Note: Do not trim bare_bright and bare_dark.  They still have to
            ! add up to total bare soil fraction and give the right
            ! weighted-sum albedo.
            if ((k == esub%bare_bright).or.(k == esub%bare_dark)) cycle

            if (vfc(k) < lc_trim_threshold) vfc(k) = 0
            do m=1,NMONTH
                if (vfm(k,m) < lc_trim_threshold) vfm(k,m) = 0
            end do   ! m
        end do    ! k=1,esub%ncover

        ! Or you can trim them ONLY if bare_bright +
        ! bare_dark < 0.05.
        if (vfc(esub%bare_bright) + vfc(esub%bare_dark) < lc_trim_threshold) then
            vfc(esub%bare_bright) = 0
            vfc(esub%bare_dark) = 0
        end if
        do m=1,NMONTH
            if (vfm(esub%bare_bright,m) + vfm(esub%bare_dark,m) < lc_trim_threshold) then
                vfm(esub%bare_bright,m) = 0
                vfm(esub%bare_dark,m) = 0
            end if
        end do

        ! -------------------------- Write Outputs (trimmed)
        do k=1,esub%ncover
            ! Annual LC and LAI file
            tr%io_ann_lc(k,1)%buf(i,j) = vfc(k)    ! vfn in 1km
            tr%io_ann_lai(k,1)%buf(i,j) = laic(k)  ! laic in 1km
            tr%io_ann_hgt(k,1)%buf(i,j) = hm(k)

            ! Monthly LC and LAI files
            do m=1,NMONTH
                !tr%io_mon_lc(k,m)%buf(i,j) = vfm(k,m)
                tr%io_mon_lai(k,m)%buf(i,j) = laim(k,m)  ! lainm in 1km
            end do ! m
        end do    ! k=1,esub%ncover
        ! ----------------------------------------------------


        ! ============= Part 2: trimmed_scaled

        ! rescale fractions so that they sum to 1 (except water)
        sum_vfc = sum(vfc(1:N_BARE))
        sum_vfm = sum(vfm(1:N_BARE, 1))
        if (abs(sum_vfc - sum_vfm) > 1.e-5) then !#DEBUG
            if (isgood) then
               write(ERROR_UNIT,*) 'ERROR trim:  max and monthly lc different' &
                    , sum_vfc,sum(vfm(1:N_BARE,1))
               write(ERROR_UNIT,*) vfc(1:N_BARE)
               write(ERROR_UNIT,*) vfm(1:N_BARE,1)
            end if
        endif

        if ( sum_vfc > 0.00001 ) then
            vfc(1:N_BARE) = vfc(1:N_BARE) / sum_vfc
            vfm(1:N_BARE,:) = vfm(1:N_BARE,:) / sum_vfc
        end if
        sum_vfc = sum(vfc(1:N_BARE))
        sum_vfm = sum(vfm(1:N_BARE, 1))
        if (abs(sum_vfc - sum_vfm) > 1.e-5) then !#DEBUG
           if (isgood) then
              write(ERROR_UNIT,*) 'ERROR scale1:  max and monthly lc different' &
                 ,i,j, sum_vfc, sum_vfm
              write(ERROR_UNIT,*) 'vfc',i,j,vfc(1:N_BARE)
              write(ERROR_UNIT,*) 'vfm',i,j,vfm(1:N_BARE,1)
           end if
        endif

        if (isgood.and.split_bare_soil) then
            call do_split_bare_soil(esub, bs_brightratio, &
                vfc,laic,hm,hsd,vfm,laim)
        end if


        ! -------------------------- Write Outputs (trimmed_scaled)
        do k=1,esub%ncover
            ! Annual LC and LAI file
            ts%io_ann_lc(k,1)%buf(i,j) = vfc(k)    ! vfn in 1km
            ts%io_ann_lai(k,1)%buf(i,j) = laic(k)  ! laic in 1km
            ts%io_ann_hgt(k,1)%buf(i,j) = hm(k)

            ! Monthly LC and LAI files
            do m=1,NMONTH
                !ts%io_mon_lc(k,m)%buf(i,j) = vfm(k,m)
                ts%io_mon_lai(k,m)%buf(i,j) = laim(k,m)  ! lainm in 1km
            end do ! m

            ts%io_ann_hgt(k,1)%buf(i,j) = hm(k)
        end do    ! k=1,esub%ncover
        ! ----------------------------------------------------


    end do   ! ic
    end do   ! jc

end subroutine do_part1_2_trimmed

subroutine do_part3_maxcrops(esub, IM,JM, io_bs, ts,   mc)
    type(GcmEntSet_t), intent(IN) :: esub
    type(ChunkIO_t), target, intent(IN) :: io_bs
    integer, intent(IN) :: IM,JM
    type(OutputSegment_t), intent(IN), target :: ts
    ! ------- Output
    type(OutputSegment_t), target :: mc
    ! ----------- Locals
    integer :: k,m
    integer :: c3herb_s,c4herb_s

    real*4, pointer :: vfc15(:,:), laic15(:,:)  !(i,j)
    real*4, allocatable :: vfm15(:,:,:), laim15(:,:,:) !(m,i,j)
    real*4, pointer :: hm15(:,:)
    real*4, allocatable :: hsd15(:,:) !(i,j)
    real*4, pointer :: laiccrop(:,:) !OUTPUT, max monthly
    real*4, allocatable ::  laimcrop(:,:,:) !OUTPUT, max monthly
    real*4, pointer :: hmcrop(:,:)
    real*4, allocatable :: hsdcrop(:,:)

real*4 :: xsum

    print *,'========================= Part 3: maxcrops'

    c3herb_s = esub%svm(CROPS_C3_HERB)   ! shortcut
    c4herb_s = esub%svm(CROPS_C4_HERB)   ! shortcut
print *,'c3herb',c3herb_s,c4herb_s
    ! Copy input to 3D array
    allocate(laim15(NMONTH, IM,JM))
    do m=1,NMONTH
        laim15(m,:,:) = ts%io_mon_lai(c3herb_s,m)%buf(:,:)
    end do

    ! Initialize dummy inputs
    allocate(vfm15(NMONTH,IM,JM))
    do m=1,NMONTH
        vfm15(m,:,:) = ts%io_ann_lc(c3herb_s,1)%buf(:,:)
    end do
    allocate(hsd15(IM,JM))
    hsd15 = 0d0

    ! Initialze outputs
    allocate(laimcrop(NMONTH, IM, JM))
    laimcrop = 0d0
    allocate(hsdcrop(IM,JM))
    hsdcrop = 0d0
    ! ----------------------------------------------------
    ! Initialze maxcrops output to be same as trimmed_scaled output
    do k=1,esub%ncover
        ! Annual LC and LAI file
        mc%io_ann_lc(k,1)%buf(:,:) = ts%io_ann_lc(k,1)%buf(:,:)    ! vfn in 1km
        mc%io_ann_lai(k,1)%buf(:,:) = ts%io_ann_lai(k,1)%buf(:,:)  ! laic in 1km
        mc%io_ann_hgt(k,1)%buf(:,:) = ts%io_ann_hgt(k,1)%buf(:,:)  ! laic in 1km

        ! Monthly LC and LAI files
        do m=1,NMONTH
            !mc%io_mon_lc(k,m)%buf(:,:) = ts%io_mon_lc(k,m)%buf(:,:)
            mc%io_mon_lai(k,m)%buf(:,:) = ts%io_mon_lai(k,m)%buf(:,:)  ! lainm in 1km
        end do ! m
    end do    ! k=1,esub%ncover

    ! Zero stuff out...
    mc%io_ann_lc(c3herb_s,1)%buf = 0d0
    mc%io_ann_hgt(c3herb_s,1)%buf = 0d0
    ! ----------------------------------------------------

    ! Aliases
    vfc15 => ts%io_ann_lai(c3herb_s,1)%buf
    laic15 => ts%io_ann_lc(c3herb_s,1)%buf
    hm15 => ts%io_ann_hgt(c3herb_s,1)%buf
    laiccrop => mc%io_ann_lc(c3herb_s,1)%buf
    hmcrop => mc%io_ann_hgt(c3herb_s,1)%buf

!do k=1,esub%ncover
!    print *,'sumA',k,sum(mc%io_ann_lai(k,1)%buf)
!end do

    !Generate fill-in crop cover from trimmed_scaled before doing nocrops
    !Herb crop only, since right now zero woody crops.
    call fill_crops( &
        IM,JM, io_bs, vfc15, laic15, &
        vfm15, laim15, hm15, hsd15, &
        laiccrop, laimcrop, hmcrop, hsdcrop)    ! outputs
    ! Copy output from 3D array
    do m=1,NMONTH
        mc%io_mon_lai(c3herb_s,m)%buf(:,:) = laimcrop(m,:,:)
    end do

!do k=1,esub%ncover
!    print *,'sumB',k,sum(mc%io_ann_lai(k,1)%buf)
!end do
end subroutine do_part3_maxcrops

subroutine do_part4_nocrops(esub, IM,JM, io_bs, ts,   nc)
    type(GcmEntSet_t), intent(IN) :: esub
    integer, intent(IN) :: IM,JM
    type(ChunkIO_t), target, intent(IN) :: io_bs
    type(OutputSegment_t), intent(IN), target :: ts
    ! ------- Output
    type(OutputSegment_t), target :: nc
    ! ----------- Locals
    integer :: i,j,k,m
    integer :: c3herb_s,c4herb_s,N_BARE,NOTBARE

    real*4, dimension(:,:), pointer :: bs_brightratio
    real*4, dimension(IM,JM,esub%ncover) :: vfc, laic
    real*4, dimension(NMONTH,IM,JM,esub%ncover) :: vfm, laim
    real*4, dimension(esub%ncover) :: hm,hsd
    real*4, dimension(esub%ncover,NMONTH) :: vfmx, laimx

    ! Temporary variables
    real*4 LAYER(IM,JM),s
    logical :: naturalfound, nonaturalfound
    integer :: naturalcount, nonaturalcount

    print *,'========================= Part 4: nocrops'

    c3herb_s = esub%svm(CROPS_C3_HERB)   ! shortcut
    c4herb_s = esub%svm(CROPS_C4_HERB)   ! shortcut

    N_BARE = esub%bare_dark
    NOTBARE = esub%svm(BARE_SPARSE) - 1

    bs_brightratio => io_bs%buf

    ! ==================== Part 4:
    ! Initialze nocrops output to be same as trimmed_scaled output
    do k=1,esub%ncover
        ! Annual LC and LAI file
        nc%io_ann_lc(k,1)%buf(:,:) = ts%io_ann_lc(k,1)%buf(:,:)
        nc%io_ann_lai(k,1)%buf(:,:) = ts%io_ann_lai(k,1)%buf(:,:)
        nc%io_ann_hgt(k,1)%buf(:,:) = ts%io_ann_hgt(k,1)%buf(:,:)

        ! Monthly LC and LAI files
        do m=1,NMONTH
            !nc%io_mon_lc(k,m)%buf(:,:) = ts%io_mon_lc(k,m)%buf(:,:)
            nc%io_mon_lai(k,m)%buf(:,:) = ts%io_mon_lai(k,m)%buf(:,:)  ! lainm in 1km
        end do ! m
    end do    ! k=1,esub%ncover
    ! ----------------------------------------------------
    nc%io_ann_lai(c3herb_s,1)%buf(:,:) = 0d0
    nc%io_ann_lai(c4herb_s,1)%buf(:,:) = 0d0
    do m=1,NMONTH
        nc%io_mon_lai(c3herb_s,m)%buf(:,:) = 0d0
        nc%io_mon_lai(c4herb_s,m)%buf(:,:) = 0d0
    end do  ! m

    ! Copy to input/output arrays
    do k=1,esub%ncover
        vfc(:,:,k) = nc%io_ann_lc(k,1)%buf(:,:)
        laic(:,:,k) = nc%io_ann_lai(k,1)%buf(:,:)
        do m=1,NMONTH
            vfm(m,:,:,k) = nc%io_ann_lc(k,1)%buf(:,:)
            laim(m,:,:,k) = nc%io_mon_lai(k,m)%buf(:,:)
        end do   ! m
    end do

    ! ------------------ Run the algorithm

    nonaturalcount=0
    LAYER(:,:) = 0.d0 !Map the cells with no natural veg
    do j=1,JM
    do i=1,IM
        s = sum( vfc(i,j,1:NOTBARE) ) 
        if (sum(vfc(i,j,c3herb_s:c4herb_s))>0.d0) then !Cell has crops
            !write(*,*) 'Cell has crops',i,j,s
            !Replace with closest non-zero natural cover
            call replace_crops(esub, IM,JM,i,j, &
                bs_brightratio, &
                vfc,vfm,laic,laim,naturalfound)
            if (.not.naturalfound) then
                nonaturalcount = nonaturalcount+1
                LAYER(i,j) = 1.
            end if
        end if !else do nothing
    end do
    end do
    write(*,*) 'Replaced crops'

    !Rescale after removing crops
    do i=1,IM
    do j=1,JM
        ! Rescale always to account for removal of crops 
        ! and to make coasts sum to 1.
        s = sum( vfc(i,j,1:N_BARE) ) 
        if ( (s.gt.0.).and.(s.ne.1.) ) then
            vfc(i,j,1:N_BARE) = vfc(i,j,1:N_BARE) / s
        endif

        do m=1,NMONTH
            s = sum( vfm(m,i,j,1:N_BARE) ) 
            if ( (s.gt.0.).and.(s.ne.1.) ) then
                vfm(m,i,j,1:N_BARE) = vfm(m,i,j,1:N_BARE) / s
            endif
        end do
       
        s = sum( vfc(i,j,1:N_BARE) ) !#DEBUG
        if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then 
            write(*,*) 'ERROR nocrops:  max and monthly lc differ' &
               ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
            !write(*,*) vfc(i,j,1:N_BARE)
            !write(*,*) vfm(1,i,j,1:N_BARE)
            do m=1,NMONTH
                vfm(m,i,j,1:N_BARE) = vfc(i,j,N_BARE)
            enddo
        endif
    end do
    end do

    !write(93) 'All-crop cells with no near natural veg.', LAYER

    if (split_bare_soil) then
        do i=1,IM
        do j=1,JM
            ! Loade data from arrays
            do k=1,esub%ncover
                hm(k) = nc%io_ann_hgt(k,1)%buf(i,j)
                hsd(k) = 0d0
                do m=1,NMONTH
                    vfmx(k,m) = vfm(m,i,j,k)
                    laimx(k,m) = laim(m,i,j,k)
                end do
            end do


            call do_split_bare_soil( &
                esub, bs_brightratio(i,j), &
                vfc(i,j,:),laic(i,j,:), &
                hm,hsd,vfmx,laimx)

            ! Store answers
            do k=1,esub%ncover
                nc%io_ann_hgt(k,1)%buf(i,j) = hm(k)
                do m=1,NMONTH
                    vfm(m,i,j,k) = vfmx(k,m)
                    laim(m,i,j,k) = laimx(k,m)
                end do
            end do

        end do
        end do

    end if   ! split_bare_soil


    ! ----------------- Copy back to buffers
    ! Copy to input/output arrays
    do k=1,esub%ncover
        nc%io_ann_lc(k,1)%buf(:,:) = vfc(:,:,k)
        nc%io_ann_lai(k,1)%buf(:,:) = laic(:,:,k)
        do m=1,NMONTH
            nc%io_ann_lc(k,1)%buf(:,:) = vfm(m,:,:,k)
            nc%io_mon_lai(k,m)%buf(:,:) = laim(m,:,:,k)
        end do   ! m
    end do

end subroutine do_part4_nocrops

subroutine do_trim(rw, esub)
    type(ReadWrites_t) :: rw
    type(GcmEntSet_t), target, intent(IN) :: esub
    ! ----------- Locals
    class(EntSet_T), pointer :: esub_p
    type(Chunker_t) :: chunker_pu    ! pure

    ! -------- Inputs: pure2
    type(ChunkIO_t) :: ioall_ann_lc,io_ann_lc(esub%ncover)
    type(ChunkIO_t) :: ioall_ann_hgt(1),io_ann_hgt(esub%ncover,1)
    type(ChunkIO_t) :: ioall_ann_lai(1),io_ann_lai(esub%ncover,1)
    type(ChunkIO_t) :: ioall_mon_lai(NMONTH),io_mon_lai(esub%ncover,NMONTH)
    type(ChunkIO_t) :: io_bs

    ! -------- Outputs
    type(OutputSegment_t) :: tr    ! trimmed
    type(OutputSegment_t) :: ts    ! trimmed_scaled
    type(OutputSegment_t) :: mc    ! maxcrops
    type(OutputSegment_t) :: nc    ! nocrops


    integer :: IM,JM
    integer :: k,m, n_bare
    real*4 :: sm
    logical :: isgood    ! .true. if this gridcell was processed by earlier stages (A00,A01,etc)
    type(FileInfo_t) :: info

    esub_p => esub

    call chunker_pu%init(IMLR,JMLR,  0,0,'', 300, 1, 30, (/1,1/))

    ! --- Inputs: Same for annual vs. monthly
    call chunker_pu%nc_open_set( &
        esub_p, io_ann_lc, &
        LAI_SOURCE, 'M', 'lc', 2004, 'purelr', '1.1')

    ! Bare Soil Brightness Ratio
    call chunker_pu%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'bs_brightratio', 2004, 'purelr', '1.1')
    call chunker_pu%nc_open(io_bs, LC_LAI_ENT_DIR, &
        info%dir, 'bs_brightratio.nc', info%vname, 1)

    ! Simard Heights
    call chunker_pu%nc_open_set( &
        esub_p, io_ann_hgt, &
        LAI_SOURCE, 'M', 'hgt', 2004, 'purelr', '1.1')

    ! laimax
    call chunker_pu%nc_open_set( &
        esub_p, io_ann_lai(:,1), &
        LAI_SOURCE, 'M', 'laimax', 2004, 'purelr', '1.1')


    do m=1,NMONTH
        call chunker_pu%nc_open_set( &
            esub_p, io_mon_lai(:,m), &
            LAI_SOURCE, 'M', 'lai', 2004, 'purelr', '1.1', &
            doytype='month', idoy=m)
    end do

    ! --------------------- Outputs: trimmed
    call tr%open('trimmed', esub_p)
    call ts%open('trimmed_scaled', esub_p)
!    call mc%open('maxcrops', esub_p)
!    call nc%open('nocrops', esub_p)
    call mc%open('trimmed_scaled_crops_ext', esub_p)
    call nc%open('trimmed_scaled_nocrops', esub_p)

    call chunker_pu%nc_check(rw=rw)
    call tr%chunker%nc_check(rw=rw)
    call ts%chunker%nc_check(rw=rw)
    call mc%chunker%nc_check(rw=rw)
    call nc%chunker%nc_check(rw=rw)

    ! Check for only one chunk
    if ((chunker_pu%nchunk(1)/=1).or.(chunker_pu%nchunk(2)/=1)) then
        write(ERROR_UNIT,*) 'nchunk must be (/1,1/) for A08, due to Part 3'
        STOP
    end if

    ! --------------------- Process things: only one chunk
    call chunker_pu%move_to(1,1)
    call tr%chunker%move_to(1,1)
    call ts%chunker%move_to(1,1)
    call mc%chunker%move_to(1,1)
    call nc%chunker%move_to(1,1)

    IM = chunker_pu%chunk_size(1)
    JM = chunker_pu%chunk_size(2)

    call do_part1_2_trimmed(esub, IM,JM, &
        io_ann_lc, io_bs, io_ann_hgt, io_ann_lai, io_mon_lai,    tr, ts)
    call do_part3_maxcrops(esub, IM,JM, io_bs, ts,    mc)
    call do_part4_nocrops(esub, IM,JM, io_bs, ts,   nc)

    ! Compute checksums for each output segment
    call tr%checksum(esub_p)
    call ts%checksum(esub_p)
    call mc%checksum(esub_p)
    call nc%checksum(esub_p)

    !call chunker_pu%write_chunks
    call tr%chunker%write_chunks
    call ts%chunker%write_chunks
    call mc%chunker%write_chunks
    call nc%chunker%write_chunks


    call chunker_pu%close_chunks
    call tr%chunker%close_chunks
    call ts%chunker%close_chunks
    call mc%chunker%close_chunks
    call nc%chunker%close_chunks

end subroutine do_trim

end module a08_mod

! =========================================================
program regrid
    use a08_mod
implicit none
    type(GcmEntSet_t), target :: esub
    type(ReadWrites_t) :: rw
    call rw%init("A08_trim", 1000,1000)

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_trim(rw, esub)
    call rw%write_mk
end program regrid
