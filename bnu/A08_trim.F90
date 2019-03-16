module a07_mod
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


subroutine do_trim(esub)
    type(GcmEntSet_t), intent(IN) :: esub
    ! ----------- Locals
    type(Chunker_t) :: chunker_pu    ! pure
    type(Chunker_t) :: chunker_tr    ! trimmed
    type(Chunker_t) :: chunker_ts    ! trimmed scaled
    type(Chunker_t) :: chunker_nc    ! no crops

    ! -------- Inputs: pure2
    type(ChunkIO_t) :: ioall_lc,io_lc(esub%ncover)
    type(ChunkIO_t) :: ioall_ann_lai(1),io_ann_lai(esub%ncover,1)
    type(ChunkIO_t) :: ioall_mon_lai(NMONTH),io_mon_lai(esub%ncover,NMONTH)


    ! --- Inputs: Same for annual vs. monthly
    call chunker%nc_open(ioall_lc, LC_LAI_ENT_DIR, &
        'purelr/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
    do k = 1,NENT20
        call chunker%nc_reuse_var(ioall_lc, io_lc(k), (/1,1,k/))
    enddo

    ! Bare Soil Brightness Ratio
    call chunker%nc_open(io_bs, LC_LAI_ENT_DIR, &
        'purelr/', 'bs_brightratio.nc', 'bs_brightratio', 1)

! Do we need this?
!    ! Simard Heights
!    call chunker%nc_open(io_height, LC_LAI_ENT_DIR, &
!        'purelr/', 'simard_height.nc', 'SimardHeights', 0)


    ! laimax
    call chunker%nc_open(ioall_laiin(1), LC_LAI_ENT_DIR, &
        'purelr/annual/', 'entmm29_ann_laimax.nc', 'lai', 0)
    do k = 1,NENT20
        call chunker%nc_reuse_var(ioall_laiin(1), io_laiin(k,1), (/1,1,k/))
    enddo

    ! --------------------- Process things

    ! Use these loop bounds for testing...
    ! it chooses a land area in Asia
    do jchunk = jc0,jc1
    do ichunk = ic0,ic1

        call chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            ! Compute overall NetCDF index of current cell
            ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
            jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1


            ! -------------------------- Read Inputs
            do k=1,esub%ncover
                ! Annual LC and LAI file
                vfc(k) = io_ann_lc(k)%buf(ic,jc)    ! vfn in 1km
                laic(k) = io_ann_lai(k)%buf(ic,jc)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
                    vfm(k,imonth) = vfn(k)   ! vfnm in 1km; Re-use annual LC
                    laim(k,month) = io_mon_lai(k,imonth)%buf(ic,jc)  ! lainm in 1km
                end imonth

                ! Height file
                hmn(k) = io_height(k)%buf(ic,jj)

                ! Bare soil brightness ratiaafo
                bs_brightratio = io_bs_brightratio%buf(ic,jc)
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------

            ! By definition, N_BARE indexes the last bare covertype
            !    If not yet split: = BARE_SPARSE
            !    If split already: = BARE_DARK
            ! We assume we've already been split in A04...A06
            n_bare = esub%bare_dark

            ! ------------------------------------------------------
            ! convert arid adapted shrub with lai < .15 to bare soil
            ! and restrict lai >= .15 
            arid_shrub_s = esub%svm(ARID_SHRUB)   ! shortcut

            if( vfc(arid_shrub_s) > .0 .and. laic(arid_shrub_s) < .15 ) then

                call convert_vf(vfc(N_BARE), laic(N_BARE), &
                     vfc(arid_shrub_s), laic(arid_shrub_s), .15 )
                                    ! lai >= .15
                if (vfc(arid_shrub_s).le.0.0) then !ERROR CHECK
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

                call convert_vfh( &
                     vfh(N_BARE),hm(N_BARE),hsd(N_BARE), &
                     vfh(arid_shrub_s),hm(arid_shrub_s),hsd(arid_shrub_s), vfc(arid_shrub_s))
            
            end if   ! (vfc(arid_shrub_s) > .0 .and. laic(arid_shrub_s) < .15 ) then

            s = sum(vfc(1:N_BARE))
            if (s.ne.sum(vfm(1,1:N_BARE))) then !#DEBUG
               write(ERROR_UNIT,*) 'ERROR trim:  max and monthly lc different' &
                    , s,sum(vfm(1,1:N_BARE))
               write(ERROR_UNIT,*) vfc(1:N_BARE)
               write(ERROR_UNIT,*) vfm(1,1:N_BARE)
            endif

#ifdef SPLIT_BARE_SOIL
TODO: split_bare_soil needs to do just one gridcell, remove the loop
            call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE &
               ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd &
               ,titlec, titlem, titleh,res_out)
#endif



            ! -------------------------- Write Outpus (trimmed)
            do k=1,esub%ncover
                ! Annual LC and LAI file
                io_ann_lc_t(k)%buf(ic,jc) = vfn(k)    ! vfn in 1km
                io_ann_lai_t(k)%buf(ic,jc) = lain(k)  ! laic in 1km

                ! Monthly LC and LAI files
                do imonth=1,12
                    io_mon_lc_t(k,imonth)%buf(ic,jc) = vfm(k,month)
                    io_mon_lai_t(k,imonth)%buf(ic,jc) = laim(k,month)  ! lainm in 1km
                end imonth
            end do    ! k=1,esub%ncover
            ! ----------------------------------------------------
        end do   ! ic
        end do   ! jc

        call chunker%write_chunks

    end do    ! ichunk
    end do    ! jchunk
end subroutine do_trim



end module a07_mod

! =========================================================
program regrid
    use a07_mod
implicit none
    call do_regrid_all_lais
end program regrid
