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

CONTAINS


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
!This checks for cases if BARE is original total or was previously split.
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
            write(*,*) 'ERR2m: bs_brightratio',vft &
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



subroutine do_trim(esub)
    type(GcmEntSet_t), intent(IN) :: esub
    ! ----------- Locals
    type(Chunker_t) :: chunker_pu    ! pure
    type(Chunker_t) :: chunker_tr    ! trimmed
    type(Chunker_t) :: chunker_ts    ! trimmed scaled
    type(Chunker_t) :: chunker_nc    ! no crops

    ! -------- Inputs: pure2
    type(ChunkIO_t) :: ioall_ann_lc,io_ann_lc(esub%ncover)
    type(ChunkIO_t) :: ioall_ann_height(1),io_ann_height(esub%ncover,1)
    type(ChunkIO_t) :: ioall_ann_lai(1),io_ann_lai(esub%ncover,1)
    type(ChunkIO_t) :: ioall_mon_lai(NMONTH),io_mon_lai(esub%ncover,NMONTH)
    type(ChunkIO_t) :: io_bs

    ! -------- Outputs: trimmed
    type(ChunkIO_t) :: ioall_ann_lc_tr(1), io_ann_lc_tr(esub%ncover,1)
    type(ChunkIO_t) :: ioall_ann_lai_tr(1), io_ann_lai_tr(esub%ncover,1)
    type(ChunkIO_t) :: ioall_mon_lai_tr(NMONTH), io_mon_lai_tr(esub%ncover,NMONTH)


    integer :: ichunk,jchunk,ic,jc,ii,jj
    integer :: k,m,imonth, n_bare
    real*4 :: sm

    ! Annual LC and LAI (LC doesn't change so annual LC is used everywhere)
    real*4, dimension(esub%ncover) :: vfc,laic,hm,hsd,vfh
    ! Monthly LC and LAI (NOTE: Monthly LC is set to same as annual)
    real*4, dimension(esub%ncover,NMONTH) :: vfm,laim
    real*4 :: bs_brightratio
    integer :: arid_shrub_s   ! Shortcut indices


    call chunker_pu%init(IMLR,JMLR,  0,0,'', 100, 1)
    call chunker_tr%init(IMLR,JMLR,  0,0,'', 1, 500)


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


    ! --------------------- Outputs: trimmed
    call chunker_tr%nc_create(ioall_ann_lc_tr(1), &
        weighting(chunker_tr%wta1, 1d0, 0d0), &
        'trimmed/annual/', 'ent29_ann_lc', 'lc', &
        'Ent Landcover from A04', '1', 'Land Cover', &
        esub%mvs, esub%layer_names(), create_lr=.false.)
    do k=1,esub%ncover
        call chunker_tr%nc_reuse_var(ioall_ann_lc_tr(1), io_ann_lc_tr(k,1), &
            (/1,1,k/), weighting(chunker_tr%wta1, 1d0,0d0))
    enddo

    call chunker_tr%nc_create(ioall_ann_lai_tr(1), &
        weighting(chunker_tr%wta1, 1d0, 0d0), &
        'trimmed/annual/', 'ent29_ann_lai', 'lai', &
        'Ent maximum LAI for year', 'm^2 m-2', 'Leaf Area Index', &
        esub%mvs, esub%layer_names(), create_lr=.false.)
    do k=1,esub%ncover
        call chunker_tr%nc_reuse_var(ioall_ann_lai_tr(1), io_ann_lai_tr(k,1), &
            (/1,1,k/), weighting(ioall_ann_lc_tr(k)%buf, 1d0,0d0))
    enddo


    do imonth=1,NMONTH
        call chunker_tr%nc_create(ioall_mon_lai_tr(imonth), &
            weighting(chunker_tr%wta1, 1d0, 0d0), &
            'trimmed/annual/', 'ent29_mon_lai', 'lai', &
            'Ent monthly LAI', 'm^2 m-2', 'Leaf Area Index', &
            esub%mvs, esub%layer_names(), create_lr=.false.)
        do k=1,esub%ncover
            call chunker_tr%nc_reuse_var(ioall_mon_lai_tr(imonth), io_mon_lai_tr(k,imonth), &
                (/1,1,k/), weighting(ioall_ann_lc_tr(k)%buf, 1d0,0d0))
        enddo
    end do


    ! --------------------- Process things

    ! Use these loop bounds for testing...
    ! it chooses a land area in Asia
#ifdef ENTGVSD_DEBUG
    do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
    do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
    do jchunk = 1,nchunk(2)
    do ichunk = 1,nchunk(1)
#endif

        call chunker_pu%move_to(ichunk,jchunk)
        call chunker_tr%move_to(ichunk,jchunk)

        do jc = 1,chunker_pu%chunk_size(2)
        do ic = 1,chunker_pu%chunk_size(1)

            ! Compute overall NetCDF index of current cell
            ii = (ichunk-1)*chunker_pu%chunk_size(1)+(ic-1)+1
            jj = (jchunk-1)*chunker_pu%chunk_size(2)+(jc-1)+1


            ! -------------------------- Read Inputs
            do k=1,esub%ncover
print *,'AA0',k
                ! Annual LC and LAI file
                vfc(k) = io_ann_lc(k)%buf(ic,jc)    ! vfn in 1km
                laic(k) = io_ann_lai(k,1)%buf(ic,jc)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
print *,'BB1',imonth
                    vfm(k,imonth) = vfc(k)   ! vfnm in 1km; Re-use annual LC
print *,'BB1.1',imonth
                    laim(k,imonth) = io_mon_lai(k,imonth)%buf(ic,jc)  ! lainm in 1km
print *,'BB1.2',imonth
                end do ! imonth

print *,'BB2'
                ! Height file
                hm(k) = io_ann_height(k,1)%buf(ic,jc)
                ! hsd = stdev
                ! vfh = vfn (LC annual), in GISS 16 pfts format (but C3 and C4 crop are summed)
                !      ===> "LC for heights", i.e. just use the global LC
                vfh(k) = vfc(k)
                hsd(k) = 0

print *,'BB3'
                ! Bare soil brightness ratiaafo
                bs_brightratio = io_bs%buf(ic,jc)
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------
print *,'AA1'

            ! By definition, N_BARE indexes the last bare covertype
            !    If not yet split: = BARE_SPARSE
            !    If split already: = BARE_DARK
            ! We assume we've already been split in A04...A06
            n_bare = esub%bare_dark

            ! ------------------------------------------------------
            ! convert arid adapted shrub with lai < .15 to bare soil
            ! and restrict lai >= .15 
            arid_shrub_s = esub%svm(ARID_SHRUB)   ! shortcut

print *,'AA2'
            if( vfc(arid_shrub_s) > .0 .and. laic(arid_shrub_s) < .15 ) then

                call convert_vf(vfc(N_BARE), laic(N_BARE), &
                     vfc(arid_shrub_s), laic(arid_shrub_s), .15 )
                                    ! lai >= .15
                if (vfc(arid_shrub_s).le.0.0) then !ERROR CHECK
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
                      write(ERROR_UNIT,*) 'vfc(arid_shrub_s) le 0.:', N_BARE,m
                      write(ERROR_UNIT,*) vfc(:)
                      write(ERROR_UNIT,*) vfm(:,m)
                      write(ERROR_UNIT,*) laic(:)
                      write(ERROR_UNIT,*) laim(:,m)
                   endif
                   if ((vfm(arid_shrub_s,m).lt.0.).or.(vfm(N_BARE,m).lt.0.)) &
                        then !CHECK ERROR
                      write(ERROR_UNIT,*) 'vfm<0:',N_BARE,m
                      write(ERROR_UNIT,*) vfc(:)
                      write(ERROR_UNIT,*) vfm(:,m)
                      write(ERROR_UNIT,*) laic(:)
                      write(ERROR_UNIT,*) laim(:,m)
                      STOP
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
            if (sm.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
               write(ERROR_UNIT,*) 'ERROR trim:  max and monthly lc different' &
                    , sm,sum(vfm(1,1:N_BARE))
               write(ERROR_UNIT,*) vfc(1:N_BARE)
               write(ERROR_UNIT,*) vfm(1,1:N_BARE)
            endif

            if (split_bare_soil) then
                call do_split_bare_soil(esub, bs_brightratio, &
                    vfc,laic,hm,hsd,vfm,laim)
            end if

print *,'AA3'
            ! -------------------------- Write Outpus (trimmed)
            do k=1,esub%ncover
                ! Annual LC and LAI file
                io_ann_lc_tr(k,1)%buf(ic,jc) = vfc(k)    ! vfn in 1km
                io_ann_lai_tr(k,1)%buf(ic,jc) = laic(k)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
                    !io_mon_lc_tr(k,imonth)%buf(ic,jc) = vfm(k,imonth)
                    io_mon_lai_tr(k,imonth)%buf(ic,jc) = laim(k,imonth)  ! lainm in 1km
                end do ! imonth
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------
        end do   ! ic
        end do   ! jc
print *,'AA4'

        call chunker_pu%write_chunks
        call chunker_tr%write_chunks

    end do    ! ichunk
    end do    ! jchunk
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
