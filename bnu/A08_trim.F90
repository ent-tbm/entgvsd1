module a08_mod
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
    type(Chunker_t) :: chunker
    type(ChunkIO_t) :: ioall_ann_lc(1), io_ann_lc(:,:)          ! (esub%ncover,1)
    type(ChunkIO_t) :: ioall_ann_lai(1), io_ann_lai(:,:)        ! (esub%ncover,1)
    type(ChunkIO_t) :: ioall_mon_lai(NMONTH), io_mon_lai(:,:)   ! (esub%ncover,NMONTH)

    procedure :: open => outputsegment_open
end type OutputSegment_t

CONTAINS

! Opens one of trimmed, trimmed_scaled or nocrops outputs
subroutine outputsegment_open(out, label, ncover)
    class(OutputSegment_t) :: out    ! Set of file handles to open
    character*(*) :: label    ! 'trimmed', 'trimmed_scaled' or 'nocrops'
    integer, intent(IN) :: ncover   ! # covertypes to allocate for; comes from esub

    call this%chunker%init(IMLR,JMLR,  0,0,'', 1, 500, (/1,1/))

    ! ------- Allocate file handles
    allocate(this%io_ann_lc(ncover,1))
    allocate(this%io_ann_lai(ncover,1))
    allocate(this%io_mon_lai(ncover,NMONTH))

    ! Open the files
    call this%chunker%nc_create(this%ioall_ann_lc(1), &
        weighting(this%chunker%wta1, 1d0, 0d0), &
        label//'/annual/', 'ent29_ann_lc', 'lc', &
        'Ent Landcover from A04', '1', 'Land Cover', &
        esub%mvs, esub%layer_names(), create_lr=.false.)
    do k=1,esub%ncover
        call this%chunker%nc_reuse_var(this%ioall_ann_lc(1), this%io_ann_lc(k,1), &
            (/1,1,k/), weighting(this%chunker%wta1, 1d0,0d0))
    enddo

    call this%chunker%nc_create(this%ioall_ann_lai(1), &
        weighting(this%chunker%wta1, 1d0, 0d0), &
        label//'/annual/', 'ent29_ann_lai', 'lai', &
        'Ent maximum LAI for year', 'm^2 m-2', 'Leaf Area Index', &
        esub%mvs, esub%layer_names(), create_lr=.false.)
    do k=1,esub%ncover
        call this%chunker%nc_reuse_var(this%ioall_ann_lai(1), this%io_ann_lai(k,1), &
            (/1,1,k/), weighting(this%ioall_ann_lc(k)%buf, 1d0,0d0))
    enddo


    do imonth=1,NMONTH
        call this%chunker%nc_create(this%ioall_mon_lai(imonth), &
            weighting(this%chunker%wta1, 1d0, 0d0), &
            label//'/annual/', 'ent29_mon_lai', 'lai', &
            'Ent monthly LAI', 'm^2 m-2', 'Leaf Area Index', &
            esub%mvs, esub%layer_names(), create_lr=.false.)
        do k=1,esub%ncover
            call this%chunker%nc_reuse_var(this%ioall_mon_lai(imonth), this%io_mon_lai(k,imonth), &
                (/1,1,k/), weighting(this%ioall_ann_lc(k)%buf, 1d0,0d0))
        enddo
    end do


end subroutine outputsegment_open
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
    if ((vft.lt.0.).or. &
          ((vft.gt.0).and.(bs_brightratio.lt.0.))) then !Bad data
        write(ERROR_UNIT,*) 'ERR2: bs_brightratio',vft,bs_brightratio
        !Keep existing bright and dark
    else             !Good data
        vfc(esub%bare_bright) = vft*bs_brightratio
        vfc(esub%bare_dark) = vft*(1.-bs_brightratio)
        laic(esub%bare_bright) = 0.
        laic(esub%bare_dark) = 0.
    end if

    do m=1,12
        if ((vft.lt.0.).or. &
           ((vft.gt.0.).and.(bs_brightratio.lt.0.))) then !Bad data
            write(ERROR_UNIT,*) 'ERR2m: bs_brightratio',vft &
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

subroutine replace_crops(IMn,JMn,KM,i,j !,s &
           ,bs_brightratio &
           ,vfc,laic,laim &
           ,naturalvegfound)  !,flag)
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
    integer, intent(in) :: IMn,JMn,KM,i, j
    real*4, intent(in) :: bs_brightratio(:,:) !(IMn,JMn) !Fraction of bare that is bright.
    real*4 :: vfc(:,:,:) , laic(:,:,:)
    !real*4 :: vfm(:,:,:,:)
    real*4 :: laim(12,IMn,JMn,KM)
    logical,intent(out) :: naturalvegfound
    !----Local----
    integer :: m, k, ii, jj !cell to search for natural (i,j or adjacent)
    real*4 :: covmax
    integer :: covmaxk, covmaxii,covmaxjj
    real*4 :: covsum(14), covavglai(14)
    !real*4 :: covsumm(12,14)
    real*4 :: covavglaim(12,14)
    integer :: dg
    real*4 :: br, brcov   !bright soil ratio, br, to total bare soil, brcov
    
    if (sum(vfc(i,j,15:16)).le.0.0) then
       return !No crops, no need to do
    endif
    
    !--First check just (i,j) cell for natural veg---
    naturalvegfound = .false.
    covmax = 0.d0
    covmaxk = 0
    covmaxii = i
    covmaxjj = j
    do k=1,14         !Find max non-crop, non-bare natural cover type
       if (vfc(i,j,k) > covmax)  then
          covmax = vfc(i,j,k)
          covmaxk = k
       endif
    enddo
    if (covmax.gt.0.d0) then  !Assign dominant natural veg to crop
       naturalvegfound = .true.
       vfc(i,j,covmaxk) = vfc(i,j,covmaxk)  &
            + vfc(i,j,15)  + vfc(i,j,16)
       vfc(i,j,15:16) = 0.0 !zero out crop cover - done below
       !vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk) &
       !     + vfm(:,i,j,15)  + vfm(:,i,j,16)
       !vfm(:,i,j,15:16) = 0.0 !zero out crop cover - done below

       !If in same cell, keep original LAI in (i,j) of natural veg
       !!DON'T DO ASSIGNMENT BELOW
       !!laic(i,j,covmaxk) = laic(covmaxii,covmaxjj,covmaxk)
       !!laim(:,i,j,covmaxk) = laim(:,covmaxii,covmaxjj,covmaxk)
       !write(*,*) 'Cell has crops + natural',i,j,s
    else
       write(*,*) "No natural veg in (i,j). Checking adjacent",i,j !,s
       dg = 1
       do while ((.not.naturalvegfound).and.(dg.le.9)) !Check up to 9 grid 
       
          covsum(:) = 0.d0
          covavglai(:) = 0.d0
          !covsumm(:,:) = 0.
          covavglaim(:,:) = 0.
          do ii=max(1,i-dg),min(i+dg,IMn)
             do jj=max(1,j-dg),min(j+dg,JMn)
                if ( (ii.ge.1).and.(ii.le.IMn) &
                     .and.(jj.ge.1).and.(jj.le.JMn) !in grid range &
                     .and.((ii.ne.i).or.(jj.ne.j)) ) !not the i,j center cell &
                     then
                   do k=1,14     
                   !Sum adjacent natural veg cover by type.
                      covsum(k) = covsum(k) + vfc(ii,jj,k)
                      covavglai(k) = covavglai(k) +  &
                           laic(ii,jj,k)*vfc(ii,jj,k)
                      !covsumm(:,k) = covsumm(:,k) + vfm(:,ii,jj,k)
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
          do k=1,14           !Find largest adjacent natural cover
             if (covsum(k)>covmax) then
                covmax = covsum(k)
                covmaxk = k
             endif
          enddo
          if (covmax>0.) then
             !do m=1,12  !Error check
             !   if (covsumm(m,covmaxk)<0.) then
             !      write(*,*) 'ERRc ',i,j,ii,jj,m &
             !           ,covmax,covsum(covmaxk), covsumm(m,covmaxk) &
             !           ,covavglai(covmaxk),covavglaim(m,covmaxk)
             !      STOP
             !   endif
             !enddo
             ! Assign adjacent natural cover type and LAI to crop
             naturalvegfound = .true.
             !covavglai(k) = covavglai(k)/covsum(k)  !BUG FOUND
             covavglai(covmaxk) = covavglai(covmaxk)/covsum(covmaxk)
             covavglaim(:,covmaxk) = covavglaim(:,covmaxk)/ covsum(covmaxk) ! /covsumm(:,covmaxk)
             vfc(i,j,covmaxk) = vfc(i,j,covmaxk)  &
                  + vfc(i,j,15) + vfc(i,j,16)
             vfc(i,j,15:16) = 0.0 !zero out crop cover - done below
             !vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk) &
             !     + vfm(:,i,j,15) + vfm(:,i,j,16)
             !vfm(:,i,j,15:16) = 0.0 !zero out crop cover  - done below
             !Assign LAI in (i,j) from adjacent cell
             laic(i,j,covmaxk) = covavglai(covmaxk)
             laim(:,i,j,covmaxk) = covavglaim(:,covmaxk)
             !write(*,*) 'Cell or adjacent has crops + natural',i,j,s
             write(*,*) 'Found natural veg in adjacent cells' &
                  ,i,j,covmaxk,dg
             do m=1,12 !Error check
                if (covavglaim(m,covmaxk)>10.) then
                   write(*,*) 'ERRc2 bad avg lai',i,j,ii,jj,m,covmaxk &
                        ,covsum(covmaxk),covavglaim(m,covmaxk)
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
!            vfc(i,j,14) = vfc(i,j,14) + vfc(i,j,15) + vfc(i,j,16)
!            vfm(:,i,j,14) = vfm(:,i,j,14) 
!     &           + vfm(:,i,j,15) + vfm(:,i,j,16)
!            laic(i,j,14) = (laic(i,j,15)*vfc(i,j,15) 
!     &           + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
!            do m=1,12
!               laim(m,i,j,14) = ( laim(m,i,j,15)*vfm(m,i,j,15) +
!     &              laim(m,i,j,16)*vfm(m,i,j,16) ) / 
!     &              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
!            enddo
!            write(*,*) 'Replaced Antarctic crops.'
!            naturalvegfound=.true.
!         elseif ( ((i.eq.98).and.(j.eq.36))
!     &           .or.((i.eq.139).and.(j.eq.45)) ) then !tropical islands
!            vfc(i,j,12) = vfc(i,j,12) + vfc(i,j,15) + vfc(i,j,16)
!            laic(i,j,12) = (laic(i,j,15)*vfc(i,j,15) 
!     &        + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
!            do m=1,12
!               laim(m,i,j,12) = ( laim(m,i,j,15)*vfm(m,i,j,15) +
!     &              laim(m,i,j,16)*vfm(m,i,j,16)) / 
!     &              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
!            enddo
!            write(*,*) 'Replaced island crops.'
!            naturalvegfound=.true.
!         endif
!       endif

    !*** PRESCRIBED GRID CELL FIXES FOR 720X360 *********
    if (.not.naturalvegfound) then !Assign bare
       if (bs_brightratio(i,j).eq.-1.e30) then
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
                      brcov = brcov + sum(vfc(ii,jj,17:18))
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
                br=-1.e30     !Assign undef to make sure no wrong values.
             endif
          endif
       else !bs_brightratio was defined
          br = bs_brightratio(i,j)
       endif

       !Assign bare soil fractions
       if (br.ne.-1.e30) then  !Adjacent bs_brightratio found
          vfc(i,j,17) = vfc(i,j,17) + &
               sum(vfc(i,j,15:16))*br
          vfc(i,j,18) = vfc(i,j,18) + &
               sum(vfc(i,j,15:16))*(1. - br)
          vfc(i,j,15:16) = 0. !zero out crops
          !vfm(:,i,j,17) = vfm(:,i,j,17) + &
          !     sum(vfm(:,i,j,15:16))*br
          !vfm(:,i,j,18) = vfm(:,i,j,18) + &
          !     sum(vfm(:,i,j,15:16))*(1. - br)
          !vfm(:,i,j,15:16) = 0. !zero out crops
          write(*,*) 'Assigned nearby bare' &
               ,i,j, dg, bs_brightratio(i,j),br
       elseif ((IMn.eq.720).and.(JMn.eq.360)) then
    !Adjacent bs_brightratio not found for only these cells
!     Replace specific grid cells ONLY VALID FOR IM=720,JM=360!!!
!     NOTE:  "Mauritius" and "Fiji/Amer.Samoa" might be off-grid and land mask.
!     assigning bare Found for IM=720,JM=360 = (lon, lat)
!     Below:  i,j,bs_brightratio(i,j), br = (lon, lat)
!     Antarctica 1: (pure17 has C3 crop, partial coastal cell/ice)
!     169  36  -1.00000002E+30  0.400000006 = (-95.75,-72.25)
!     168  37  -1.00000002E+30  0.400000006 = (-96.25,-71.75)
!     Antarctica 2: (pure17 has C3 crop, partial coastal cell/ice)
!     239  36  -1.00000002E+30  0.400000006 = (-60.75,-72.75)
!     204  36  -1.00000002E+30  0.400000006 = (-78.25,-72.25)
!     Mauritius (too east of):   (pure17 was everbroad, C4 crop, sparse, coast)
!     487 141 0.461078912 0.461078912  = (63.25, -19.75)
!     Fiji?  American Samoa? (pure17 was everbroad, C4 crop, sparse, coast)
!     694 179 0.353166729 0.353166729  = (163.25,-0.75)
          if (((i.eq.169).and.(j.eq.36)) !Antarctica1 &
               .or.((i.eq.168).and.(j.eq.37)) !Antarctica1 &
               .or.((i.eq.239).and.(j.eq.36)) !Antarctica2 &
               .or.((i.eq.204).and.(j.eq.36))) !Antarctica2 &
               then
                !* Assign Antarctica cells to arctic grass
             vfc(i,j,14) = vfc(i,j,14) + sum(vfc(i,j,15:16))
             vfc(i,j,15:16) = 0. !zero out crops
             laic(i,j,14) = max(laic(i,j,14), &
                  max(laic(i,j,15),laic(i,j,16)))
             do m=1,12
                !vfm(m,i,j,14)=vfm(m,i,j,14) + sum(vfm(m,i,j,15:16))
                !vfm(m,i,j,15:16) = 0. !zero out crops
                laim(m,i,j,14)=max(laim(m,i,j,14), &
                     max(laim(m,i,j,15),laim(m,i,j,16)))
             enddo
             write(*,*) 'Assigned c3 crops cells in Antartica',i,j
          endif
          if (((i.eq.487).and.(j.eq.141)) !Mauritius &
               .or.((i.eq.694).and.(j.eq.179))) !Fiji or American Samoa &
               then
             !* Assign tropical rainforest
             vfc(i,j,2) = vfc(i,j,2) + sum(vfc(i,j,15:16))
             vfc(i,j,15:16) = 0. !zero out crops
             laic(i,j,2) = max(laic(i,j,2), &
                  max(laic(i,j,15), laic(i,j,16)))
             do m=1,12
                !vfm(m,i,j,2)=vfm(m,i,j,2) + sum(vfm(m,i,j,15:16))
                !vfm(m,i,j,15:16) = 0. !zero out crops
                laim(m,i,j,2)=max(laim(m,i,j,2), &
                     max(laim(m,i,j,15),laim(m,i,j,16)))
             enddo
             write(*,*) 'Assigned c4 crops cells in islands',i,j
          endif
       else
          !Don't do any grid cell fixes
          write(*,*) 'WRONG GRID RES:  PRESCRIBED FIX NOT DONE: '
          write(*,*) '**Check array indices for grid res for nocrops'            endif !Correct grid resolution
       endif                  !Adjacent bs_brightratio not found
    endif                     !.not.(naturalvegfound)
    
    !Zero out crop cover - done above for each vfc assignment.
!      vfc(i,j,15:16) = 0.
!      vfm(:,i,j,15:16) = 0.
         
end subroutine replace_crops


subroutine fill_crops(IMn,JMn,vfc15,laic15m, &
    laim15,hm15,&
    laiccrop, laimcrop, hmcrop)

    !This performs in-fill once for herb crop (PFT15) LAI to create
    !an extended crop LAI data set for use with historical crop cover.
    integer, intent(IN) :: IMn, JMn
    real*4, intent(IN) :: vfc15(:,:), laic15(:,:)  !(i,j)
    !real*4, intent(IN) :: vfm15(:,:,:)   ! (m,i,j)
    real*4, intent(IN) ::  laim15(:,:,:) !(m,i,j)
    real*4, intent(IN) :: hm15(:,:) !(i,j)
    !real*4, intent(IN) :: hsd15(:,:) !(i,j)
    ! ------ OUTPUTS
    real*4 :: laiccrop(:,:), laimcrop(:,:,:) !OUTPUT, max monthly
    real*4 :: hmcrop(:,:)
    !real*4 :: hsdcrop(:,:)
    !--Local----
    integer :: i, j, m, ii, jj, nocropcells, dg
    real*4 :: covsum15, laiavg15, covmsum15(NMONTH), laimavg15(NMONTH)
    real*4 :: hmavg15, hsdavg15
    character*80 :: titlefoo
    real*4 :: hsd15   ! stand-in for hsd15 above

    titlefoo = 'vfc15 in - crop fill'
    write(93) titlefoo, vfc15
    titlefoo = 'laic15 in - crop fill'
    write(93) titlefoo, laic15

    hsd15 = 0d0        ! No stdev for now

    nocropcells=0
    do i=1,IMn
       do j=1,JMn
          if (vfc15(i,j).gt.0.d0) then !crops in cell, replicate
             laiccrop(i,j) = laic15(i,j)
             laimcrop(:,i,j) = laim15(:,i,j)
             hmcrop(i,j) = hm15(i,j)
             hsdcrop(i,j) = hsd15
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
             do ii=i-dg,i+dg    !ext 5 grid cells = 2.5 degrees at HXH
                do jj=j-dg,j+dg
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
                              + hsd15*vfc15(ii,jj)
                         !covmsum15(:) = covmsum15(:) + vfm15(:,ii,jj)
                         covmsum15(:) = covmsum15(:) + vfc15(ii,jj)
                         !laimavg15(:) = laimavg15(:) + laim15(:,ii,jj)*vfm15(:,ii,jj)
                         laimavg15(:) = laimavg15(:) + laim15(:,ii,jj)*vfc15(ii,jj)
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
                !hsdcrop(i,j) = hsdavg15
                laimcrop(:,i,j) = laimavg15(:)
             else
                nocropcells = nocropcells + 1
             endif
          endif
       enddo
    enddo

    write(*,*) 'Crop fill-in: no-crop cells:',nocropcells

end subroutine fill_crops

subroutine do_part1_2(...)
    ! ------------- Locals
    ! Annual LC and LAI (LC doesn't change so annual LC is used everywhere)
    real*4, dimension(esub%ncover) :: vfc,laic,hm,hsd,vfh
    ! Monthly LC and LAI (NOTE: Monthly LC is set to same as annual)
    real*4, dimension(esub%ncover,NMONTH) :: vfm,laim
    real*4 :: bs_brightratio
    integer :: arid_shrub_s   ! Shortcut indices

    ichunk = 1
    jchunk = 1
    ! do ichunk,jchunk...3

        call chunker_pu%move_to(ichunk,jchunk)
        call tr%chunker%move_to(ichunk,jchunk)
        call ts%chunker%move_to(ichunk,jchunk)
        call mc%chunker%move_to(ichunk,jchunk)
        call nc%chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker_pu%chunk_size(2)
        do ic = 1,chunker_pu%chunk_size(1)

            ! Compute overall NetCDF index of current cell
            ii = (ichunk-1)*chunker_pu%chunk_size(1)+(ic-1)+1
            jj = (jchunk-1)*chunker_pu%chunk_size(2)+(jc-1)+1

            ! -------------------------- Read Inputs
            do k=1,esub%ncover
                ! Annual LC and LAI file
                vfc(k) = io_ann_lc(k)%buf(ic,jc)    ! vfn in 1km
                laic(k) = io_ann_lai(k,1)%buf(ic,jc)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
                    vfm(k,imonth) = vfc(k)   ! vfnm in 1km; Re-use annual LC
                    laim(k,imonth) = io_mon_lai(k,imonth)%buf(ic,jc)  ! lainm in 1km
                end do ! imonth

                ! Height file
                hm(k) = io_ann_height(k,1)%buf(ic,jc)
                ! hsd = stdev
                ! vfh = vfn (LC annual), in GISS 16 pfts format (but C3 and C4 crop are summed)
                !      ===> "LC for heights", i.e. just use the global LC
                vfh(k) = vfc(k)
                hsd(k) = 0

                ! Bare soil brightness ratiaafo
                bs_brightratio = io_bs%buf(ic,jc)

                ! Make sure this gridcell has been defined in inputs
                isgood = (bs_brightratio /= undef)
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------

            ! By definition, N_BARE indexes the last bare covertype
            !    If not yet split: = BARE_SPARSE
            !    If split already: = BARE_DARK
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
                        laic(arid_shrub_s),laic(N_BARE),ii,jj
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

                ! **hm = Simard hegihts
                ! **hsd = stdev of heights
                ! **vfh = "LC for heights"
                call convert_vfh( &
                     vfh(N_BARE),hm(N_BARE),hsd(N_BARE), &
                     vfh(arid_shrub_s),hm(arid_shrub_s),hsd(arid_shrub_s), vfc(arid_shrub_s))
            end if   ! (vfc(arid_shrub_s) > .0 .and. laic(arid_shrub_s) < .15 ) then

            sm = sum(vfc(1:N_BARE))
            if (sm.ne.sum(vfm(1:N_BARE,1))) then !#DEBUG
                if (isgood) then
                   write(ERROR_UNIT,*) 'ERROR trim:  max and monthly lc different' &
                        , sm,sum(vfm(1:N_BARE,1))
                   write(ERROR_UNIT,*) vfc(1:N_BARE)
                   write(ERROR_UNIT,*) vfm(1:N_BARE,1)
                end if
            endif

            if (isgood.and.split_bare_soil) then
                call do_split_bare_soil(esub, bs_brightratio, &
                    vfc,laic,hm,hsd,vfm,laim)
            end if

            ! -------------------------- Write Outputs (trimmed)
            do k=1,esub%ncover
                ! Annual LC and LAI file
                tr%io_ann_lc(k,1)%buf(ic,jc) = vfc(k)    ! vfn in 1km
                tr%io_ann_lai(k,1)%buf(ic,jc) = laic(k)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
                    !tr%io_mon_lc(k,imonth)%buf(ic,jc) = vfm(k,imonth)
                    tr%io_mon_lai(k,imonth)%buf(ic,jc) = laim(k,imonth)  ! lainm in 1km
                end do ! imonth
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------


            ! ============= Part 2: trimmed_scaled

            ! rescale fractions so that they sum to 1
            write(*,*) 'Rescaling...'
            do j=1,JMn
               do i=1,IMn
                  s = sum( vfc(i,j,1:N_BARE) )
                  if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
                    if (isgood) then
                        write(ERROR_UNIT,*) 'ERROR scale0:  max and monthly lc different' &
                          ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
                        write(ERROR_UNIT,*) 'vfc',i,j,vfc(i,j,1:N_BARE)
                        write(ERROR_UNIT,*) 'vfm',i,j,vfm(1,i,j,1:N_BARE)
                    end if
                  endif
                  if ( s > 1.00001 ) then
      !            if ( abs(s-1.0) > 0.00001 ) then
                     vfc(i,j,1:N_BARE) = vfc(i,j,1:N_BARE) / s
                     vfm(:,i,j,1:N_BARE) = vfm(:,i,j,1:N_BARE) / s
                    !vfh(i,j,1:N_BARE) = vfh(i,j,1:N_BARE) / s !Heights are vertical, no rescale
                  endif
                  s = sum( vfc(i,j,1:N_BARE) ) 
                  if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
                     if (isgood) then
                        write(ERROR_UNIT,*) 'ERROR scale1:  max and monthly lc different' &
                           ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
                        write(ERROR_UNIT,*) 'vfc',i,j,vfc(i,j,1:N_BARE)
                        write(ERROR_UNIT,*) 'vfm',i,j,vfm(1,i,j,1:N_BARE)
                     end if
                  endif
               enddo
            enddo

            if (split_bare_soil) then
                call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE &
                    ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd &
                    ,titlec, titlem, titleh,res_out)
            end if


            ! -------------------------- Write Outputs (trimmed_scaled)
            do k=1,esub%ncover
                ! Annual LC and LAI file
                ts%io_ann_lc(k,1)%buf(ic,jc) = vfc(k)    ! vfn in 1km
                ts%io_ann_lai(k,1)%buf(ic,jc) = laic(k)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
                    !ts%io_mon_lc(k,imonth)%buf(ic,jc) = vfm(k,imonth)
                    ts%io_mon_lai(k,imonth)%buf(ic,jc) = laim(k,imonth)  ! lainm in 1km
                end do ! imonth
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------


        end do   ! ic
        end do   ! jc

end subroutine do_part1_2

subroutine do_part3
        ! =============== Part 3: maxcrops
        c3herb_s = esub%svn(CROPS_C3_HERB)   ! shortcut
        c4herb_s = esub%svn(CROPS_C4_HERB)   ! shortcut

        ! Copy input to 3D array
        allocate(laim15(NMONTH, size(chunker_pu%buf,1), size(chunker_pu%buf,2)))   ! INPUT
        do m=1,NMONTH
            laim15(m,:,:) = ts%io_mon_lai(c3herb_s,m)%buf(:,:)
        end do

        ! ----------------------------------------------------
        ! Initialze maxcrops output to be same as trimmed_scaled output
        do k=1,esub%ncover
            ! Annual LC and LAI file
            mc%io_ann_lc(k,1)%buf(:,:) = ts%io_ann_lc(k,1)%buf(:,:)    ! vfn in 1km
            mc%io_ann_lai(k,1)%buf(:,:) = ts%io_ann_lai(k,1)%buf(:,:)  ! laic in 1km

            ! Monthly LC and LAI files
            do imonth=1,12
                !mc%io_mon_lc(k,imonth)%buf(:,:) = ts%io_mon_lc(k,imonth)%buf(:,:)
                mc%io_mon_lai(k,imonth)%buf(:,:) = ts%io_mon_lai(k,imonth)%buf(:,:)  ! lainm in 1km
            end do ! imonth
        end do    ! k=1,esub%ncover
        ! ----------------------------------------------------

        ! Initialize output variables
        allocate(laimcrop(NMONTH, size(chunker_pu%buf,1), size(chunker_pu%buf,2)))   ! OUTPUT
        mc%io_ann_lc(c3herb_s,1)%buf
        laimcrop = 0d0
        mc%io_ann_height(c3herb_s,1)%buf = 0d0

        !Generate fill-in crop cover from trimmed_scaled before doing nocrops
        !Herb crop only, since right now zero woody crops.
        call fill_crops( &
            size(chunker_pu%buf,1), size(chunker_pu%buf,2), &
            ts%io_ann_lc(c3herb_s,1)%buf, ts%io_ann_lai(c3herb_s,1)%buf, &   ! vfc,laic
            laim15, &
            io_ann_height(c3herb_s,1)%buf, &   ! hm15
            & ! Outputs....
            mc%io_ann_lc(c3herb_s,1)%buf, laimcrop, & !laicropext(:,:),laimcropext(:,:,:)
            mc%io_ann_height(c3herb_s,1)%buf)

        ! Copy output from 3D array
        do m=1,NMONTH
            mc%io_ann_lai(c3herb_s,m)%buf(:,:) = laimcrop(m,:,:)
        end do

end subroutine do_part3_maxcrops

subroutine do_part4_nocrops
        ! ==================== Part 4:
        ! Initialze nocrops output to be same as trimmed_scaled output
        do k=1,esub%ncover
            ! Annual LC and LAI file
            nc%io_ann_lc(k,1)%buf(:,:) = ts%io_ann_lc(k,1)%buf(:,:)    ! vfn in 1km
            nc%io_ann_lai(k,1)%buf(:,:) = ts%io_ann_lai(k,1)%buf(:,:)  ! laic in 1km

            ! Monthly LC and LAI files
            do imonth=1,NMONTH
                !nc%io_mon_lc(k,imonth)%buf(:,:) = ts%io_mon_lc(k,imonth)%buf(:,:)
                nc%io_mon_lai(k,imonth)%buf(:,:) = ts%io_mon_lai(k,imonth)%buf(:,:)  ! lainm in 1km
            end do ! imonth
        end do    ! k=1,esub%ncover
        ! ----------------------------------------------------

        nc%io_ann_lai(c3herb_s,1)%buf)(:,:) = 0d0
        nc%io_ann_lai(c4herb_s,1)%buf)(:,:) = 0d0
        do imonth=1,NMONTH
            nc%io_mon_lai(c3herb_s,imonth)%buf(:,:) = 0d0
            nc%io_mon_lai(c4herb_s,imonth)%buf(:,:) = 0d0
        end do  ! imonth


        do m=1,NMONTH
            laim4(m,:,:) = nc%io_mon_lai(c3herb_s,m)%buf(:,:)
        end do



        flag = 0
        do while (flag.le.1)
            nonaturalcount=0
            LAYER(:,:) = 0.d0 !Map the cells with no natural veg

        do jc = 1,chunker_pu%chunk_size(2)
        do ic = 1,chunker_pu%chunk_size(1)
           s = sum(nc%io_ann_lc%buf(ic,jc,1:esub%ncover))
           if (sum(nc%io_ann_lc%buf(ic,jc,c3herb_s:c4herb_s))>0.d0) then !Cell has crops
              !write(*,*) 'Cell has crops',ic,jc,s
              !Replace with closest non-zero natural cover
              call replace_crops(IMn,JMn,KM,ic,jc !,s &
                   ,bs_brightratio &
                   ,vfc,vfm,laic,laim4,naturalfound)
!??     &             ,vfc,vfm,vfc,vfm,naturalfound)
              if (.not.naturalfound) then
                 nonaturalcount = nonaturalcount+1
                 LAYER(ic,jc) = 1.
              endif
           endif !else do nothing
        enddo
      enddo
      write(*,*) 'Replaced crops'

      !* DEBUG *
      write(*,*) 'Grid cells with no natural veg: ',nonaturalcount
      titlefoo = 'Grid cells with no natural veg: '
      write(996) titlefoo, LAYER
      do k=1,18
         write(*,*) k, titlec(k)
         titlefoo = 'vfc repcrops '//titlec(k)
         write(998) titlefoo, vfc(:,:,k)
         titlefoo = 'laic repcrops '//titlec(k)
         write(998) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm repcrops'//titlem(m,4)
         write(998) titlefoo, vfm(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'laim repcrops'//titlem(m,4)
         write(998) titlefoo, laim(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'vfm repcrops'//titlem(m,c3herb_s)
         write(998) titlefoo, vfm(m,:,:,c3herb_s) !Just check crops
      enddo
#endif
      !Rescale after removing crops
      do i=1,IMn
         do j=1,JMn
           ! Rescale always to account for removal of crops 
           ! and to make coasts sum to 1.
           s = sum( vfc(ic,jc,1:N_BARE) ) 
           if ( (s.gt.0.).and.(s.ne.1.) ) then
              vfc(ic,jc,1:N_BARE) = vfc(ic,jc,1:N_BARE) / s
           endif

           do m=1,12
              s = sum( vfm(m,ic,jc,1:N_BARE) ) 
              if ( (s.gt.0.).and.(s.ne.1.) ) then
                 vfm(m,ic,jc,1:N_BARE) = vfm(m,ic,jc,1:N_BARE) / s
              endif
           enddo
           
           s = sum( vfc(ic,jc,1:N_BARE) ) !#DEBUG
           if (s.ne.sum(vfm(1,ic,jc,1:N_BARE))) then 
              write(*,*) 'ERROR nocrops:  max and monthly lc differ' &
                   ,ic,jc, s,sum( vfm(1,ic,jc,1:N_BARE) ) 
              !write(*,*) vfc(ic,jc,1:N_BARE)
              !write(*,*) vfm(1,ic,jc,1:N_BARE)
              do m=1,12
                 vfm(m,ic,jc,1:N_BARE) = vfc(ic,jc,N_BARE)
              enddo
           endif
        enddo
      enddo

      !** DEBUG after rescaling crops
      do k=1,18
         write(*,*) k, titlec(k)
         titlefoo = 'vfc rescalecrops '//titlec(k)
         write(999) titlefoo, vfc(:,:,k)
         titlefoo = 'laic rescalecrops '//titlec(k)
         write(999) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm rescalecrops'//titlem(m,4)
         write(999) titlefoo, vfm(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'laim repcrops'//titlem(m,4)
         write(999) titlefoo, laim(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'vfm rescalecrops'//titlem(m,c3herb_s)
         write(999) titlefoo, vfm(m,:,:,c3herb_s) !Just check crops
      enddo

      if (.FALSE.) then !DON'T NEED TO DO THIS ANY MORE, DONE IN replace_crops
         !CAUSES SOME WRONG RESCALING
         if (nonaturalcount.eq.0) then
            write(*,*) 'All crops cells successfully fixed.'
            flag = 2
         else
            flag = flag + 1
            if (flag.eq.1) then
               write(*,*) 'Some all-crop cells. Iterating once...' &
                    ,nonaturalcount
            elseif (flag.gt.1) then
               write(*,*) 'Remaining all-crop cells,',nonaturalcount
            endif
         endif
      else
         flag = 2 !Finish the while loop
      endif

      titlefoo = 'All-crop cells with no near natural veg.'
      write(93) titlefoo, LAYER

      enddo                     !do while flag ------

      !## HACK -nk
!      do m=1,12
!         vfm(m,:,:,:) = vfc(:,:,:)
!      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE &
           ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd &
           ,titlec, titlem, titleh,res_out)
#endif

      !** DEBUG after rescaling crops and split_bare
      do k=1,18
         write(*,*) k, titlec(k)
         titlefoo = 'vfc rescalecrops '//titlec(k)
         write(1000) titlefoo, vfc(:,:,k)
         titlefoo = 'laic rescalecrops '//titlec(k)
         write(1000) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm rescalecrops'//titlem(m,c3herb_s)
         write(1000) titlefoo, vfm(m,:,:,c3herb_s) !Just check crops
      enddo
      
      call write_output(IM,JM,titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_scaled_nocrops_" &
           "lc_lai_ent16/V"//res_out_int//"_"//filepreout &
           ,"max_trimmed_scaled_nocrops_"//trim(fversion),"   ",res_out)
      write(*,*) "trimmed, scaled, no crops"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)


      do m=1,12
         call write_output(IM,JM,titlem(m,:), vfm(m,:,:,:) &
              , laim(m,:,:,:), N_BARE,  &
           "lc_lai_ent16/V"//res_out_int//"_"//filepreout &
          ,"trimmed_scaled_nocrops_"//trim(fversion),MONTH(m),res_out)
      enddo

      call write_output_h(IM,JM,N_BARE &
           ,titleh, hm, hsd &
           ,"lc_lai_ent16/V"//res_out_int//"_"//filepreout// &
           "_height_trimmed_scaled_nocrops_"//trim(fversion) &
           ,"   ",res_out)

      foo(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
           do k=1,KM
              foo(ic,jc) = foo(ic,jc) + vfc(ic,jc,k)
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
          lai_yy(ic,jc) = laic(ic,jc,N_BARE)*vfc(ic,jc,N_BARE)
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
          if( vfc(ic,jc,N_BARE) * laic(ic,jc,N_BARE) > .0 ) then

            maxpft = maxloc( vfc(ic,jc,1:c4herb_s), 1 )

            vf_yy(ic,jc) = vfc(ic,jc,maxpft)
            lai_yy(ic,jc) = maxpft

          endif
        enddo
      enddo
      title_yy = "max fraction"
      write(10) title_yy, vf_yy(:,:)
      title_yy = "max pft"
      write(10) title_yy, lai_yy(:,:)

end subroutine do_part4_nocrops

subroutine do_trim(esub)
    type(GcmEntSet_t), intent(IN) :: esub
    ! ----------- Locals
    type(Chunker_t) :: chunker_pu    ! pure

    ! -------- Inputs: pure2
    type(ChunkIO_t) :: ioall_ann_lc,io_ann_lc(esub%ncover)
    type(ChunkIO_t) :: ioall_ann_height(1),io_ann_height(esub%ncover,1)
    type(ChunkIO_t) :: ioall_ann_lai(1),io_ann_lai(esub%ncover,1)
    type(ChunkIO_t) :: ioall_mon_lai(NMONTH),io_mon_lai(esub%ncover,NMONTH)
    type(ChunkIO_t) :: io_bs

    ! -------- Outputs
    type(OutputSegment_t) :: tr    ! trimmed
    type(OutputSegment_t) :: ts    ! trimmed_scaled
    type(OutputSegment_t) :: mc    ! maxcrops
    type(OutputSegment_t) :: nc    ! nocrops


    integer :: ichunk,jchunk,ic,jc,ii,jj
    integer :: k,m,imonth, n_bare
    real*4 :: sm
    logical :: isgood    ! .true. if this gridcell was processed by earlier stages (A00,A01,etc)

    ! Part 3...
    integer :: c3herb_s,c4herb_s


    call chunker_pu%init(IMLR,JMLR,  0,0,'', 300, 1, (/1,1/))

    ! --- Inputs: Same for annual vs. monthly
    call chunker_pu%nc_open(ioall_ann_lc, LC_LAI_ENT_DIR, &
        'purelr/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
    do k = 1,esub%ncover
        call chunker_pu%nc_reuse_var(ioall_ann_lc, io_ann_lc(k), (/1,1,k/))
    enddo

    ! Bare Soil Brightness Ratio
    call chunker_pu%nc_open(io_bs, LC_LAI_ENT_DIR, &
        'purelr/annual/', 'bs_brightratio.nc', 'bs_brightratio', 1)

    ! Simard Heights
    call chunker_pu%nc_open(ioall_ann_height(1), LC_LAI_ENT_DIR, &
        'purelr/annual/', 'entmm29_ann_height.nc', 'SimardHeights', 0)
    do k = 1,esub%ncover
        call chunker_pu%nc_reuse_var(ioall_ann_height(1), io_ann_height(k,1), (/1,1,k/))
    enddo


    ! laimax
    call chunker_pu%nc_open(ioall_ann_lai(1), LC_LAI_ENT_DIR, &
        'purelr/annual/', 'entmm29_ann_laimax.nc', 'lai', 0)
    do k = 1,esub%ncover
        call chunker_pu%nc_reuse_var(ioall_ann_lai(1), io_ann_lai(k,1), (/1,1,k/))
    enddo


    do imonth=1,NMONTH
        call chunker_pu%nc_open(ioall_mon_lai(imonth), LC_LAI_ENT_DIR, &
            'purelr/monthly/', 'entmm29_'//trim(MONTH(imonth))//'_lai.nc', 'lai', 0)
        do k = 1,esub%ncover
            call chunker_pu%nc_reuse_var(ioall_mon_lai(imonth), io_mon_lai(k,imonth), (/1,1,k/))
        enddo
    end do


    ! --------------------- Outputs: trimmed
    call tr%open(tr, 'trimmed', esub%ncover)
    call ts%open(ts, 'trimmed_scaled', esub%ncover)
    call ts%open(mc, 'maxcrops', esub%ncover)
    call nc%open(nc, 'nocrops', esub%ncover)

    ! Check for only one chunk
    if ((chunker_pu%nchunk(1)/=1).or.(chunker_pu%nchunk(2)/=1)) then
        write(ERROR_UNIT,*) 'nchunk must be (/1,1/) for A08, due to Part 3'
        STOP
    end if


    ! --------------------- Process things: only one chunk
    call chunker_pu%move_to(ichunk,jchunk)
    call tr%chunker%move_to(ichunk,jchunk)
    call ts%chunker%move_to(ichunk,jchunk)
    call mc%chunker%move_to(ichunk,jchunk)
    call nc%chunker%move_to(ichunk,jchunk)

    call do_part1_2_trimmed(esub, io_ann_lc, io_bs, io_ann_height, io_ann_lai, io_mon_lai,    tr, ts)
    call do_part3_maxcrops(esub, ts,    mc)
    call do_part4_nocrops(esub, ts, mc,    nc)

    call chunker_pu%write_chunks
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
    type(GcmEntSet_t) :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_trim(esub)
end program regrid
