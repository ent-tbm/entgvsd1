! Trims off tiny fractions and preserves total LAI for the gridcell.
! Alters the cover amount if the LAI is smaller than the main one.
! If LAI is a little biger, might increase.
!
! I don't see where the program is doing that.

! A04 only creates the pure dataset.  It does NOT trim.

module a04_mod
    use conversions, only : convert_vf
    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod

implicit none


    ! Combine C3 and C4 crops into one PFT, for ModelE
    logical, parameter :: combine_crops_c3_c4 = .true.
    ! Split the bare soil into dark and light, to get the right albedo
    logical, parameter :: split_bare_soil = .true.

CONTAINS

subroutine do_reindex(esub)
    type(GcmEntSet_t), intent(IN) :: esub

    type(Chunker_t) :: chunker
    integer :: ichunk,jchunk, ic,jc
    ! Input files
    type(ChunkIO_t) :: io_lcin(NENT20), io_laiin(NENT20), io_bs
    type(ChunkIO_t) :: ioall_simard, io_simard(NENT19)
    ! Output files
    type(ChunkIO_t) :: ioall_laiout, io_laiout(esub%ncover)
    type(ChunkIO_t) :: ioall_sum, io_sum_lc, io_sum_lai


    ! Input values, max, monthly
    real*4 :: bs_brightratio   !are soil brightratio
    ! Renumbered input values
    real*4 vfc(esub%ncover)    ! = vfn = io_lcin
    real*4 laic(esub%ncover)   ! = lain = io_lcaiin
    real*4 hmc(esub%ncover)     ! = hmn = io_simard



    real*4 area
    real*4 c3c4
    real*4 am
    ! Converted values

    character*80 :: titlem(12,18)
    ! Converted values - heights
    real*4 vfh(esub%ncover), hsd(esub%ncover)

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

    real*4 :: vfc_tmp

    integer i
    integer :: ii,jj   ! Position in overall NetCDF variable
    integer :: k,k19,ksub, ri
    !integer :: io, in, jn, maxpft, kx, m
    integer :: maxpft
    real*8 lat
    !real*4 foolc(IM1km,JM1km),foolai(IM1km,JM1km)
    integer count
    real*4 :: val
    real*4 :: vf_bare_sparse

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', &
        100, &   ! # files to >= (N_VEG + N_BARE)*(LC + LAI) + BARE_BRIGHTRATIO + SIMARD = 42
        120)     ! # files to write >= N_LAIMAX + 3*CHECKSUMS

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! lcmax
    do k=1,NENT20
        ! USE the land cover we computed in a previous step!!!!!!
        ! TODO: LAI3g????
        ! PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
        call chunker%nc_open(io_lcin(k), LC_LAI_ENT_DIR, &
            'EntMM_lc_laimax_1kmx1km/', &
            trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lc.nc', &
            trim(ent20%abbrev(k)), 1)
    end do

    ! laimax
    do k=1,NENT20
        call chunker%nc_open(io_laiin(k), LC_LAI_ENT_DIR, &
            'EntMM_lc_laimax_1kmx1km/', &
            trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lai.nc', &
            trim(ent20%abbrev(k)), 1)
    enddo

    ! Bare Soil Brightness Ratio
    call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
        'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)

    ! Our processed Simard heights
    call chunker%nc_open(ioall_simard, LC_LAI_ENT_DIR, &
        'height/', 'EntGVSDmosaic17_height_1kmx1km_lai3g.nc', &
        'SimardHeights', 0)
    do k19=1,NENT19
        call chunker%nc_reuse_var(ioall_simard, io_simard(k19), (/1,1,k19/), 'r')
    end do

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    ! laimax_pure
    call chunker%nc_create(ioall_laiout, &
        weighting(chunker%wta1, 1d0, 0d0), &   ! TODO: Scale by _lc
        '16/nc/', 'V1km_EntGVSDv1.1_BNU16_laimax_pure')
    do ksub=1,esub%ncover
        call chunker%nc_reuse_file(ioall_laiout, io_laiout(ksub), &
            'lai_'//trim(esub%abbrev(ksub)), trim(esub%title(ksub)), &
            'm2 m-2', trim(esub%title(ksub)), &
            weighting(chunker%wta1,1d0,0d0))    ! TODO: Weighting???
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

    call chunker%nc_check('A04_trim_laimax_1kmx1km')
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    ! Use these loop bounds for testing...
    ! it chooses a land area in Asia
#ifdef ENTGVSD_DEBUG
    do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
    do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
    do jchunk = 1,nchunk(2)
    do ichunk = 1,nchunk(1)
#endif

        call chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

           ! Compute overall NetCDF index of current cell
           ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
           jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

            ! get bs bright ratio
            bs_brightratio = io_bs%buf(ic,jc)

            ! Check if mismatch lc or lai values (one is zero and the other not)
            ! call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)


            ! =============== Convert to GISS 16 pfts format
            ! First simple transfers
            do ri=1,esub%nremap
                ksub = esub%remap(1,ri)
                k = esub%remap(2,ri)
                k19 = ent19%svm(k)

                vfc(ksub) = io_lcin(k)%buf(ic,jc)
                hmc(ksub) = io_simard(k19)%buf(ic,jc)
                laic(ksub) = io_laiin(k)%buf(ic,jc)
            end do

            ! Then more complex stuff
            if (combine_crops_c3_c4) then
                !lc laimax
                c3c4 = io_lcin(CROPS_C3_HERB)%buf(ic,jc) + io_lcin(CROPS_C4_HERB)%buf(ic,jc)
                if (c3c4 > 0.) then
                    laic(esub%crops_herb) = ( &
                          io_lcin(CROPS_C3_HERB)%buf(ic,jc) *io_laiin(CROPS_C3_HERB)%buf(ic,jc) &
                        + io_lcin(CROPS_C4_HERB)%buf(ic,jc) *io_laiin(CROPS_C4_HERB)%buf(ic,jc) &
                        ) / c3c4
                    vfc(esub%crops_herb) = c3c4
                else
                    laic(esub%crops_herb) = 0.
                    vfc(esub%crops_herb) = 0.
                endif
            end if


            ! =========== Partition bare/sparse LAI to actual LC types.
            ! Bare sparse has no vegetation type assigned to it; but it
            ! has non-zero LAI that must be treated.  So we take that
            ! non-zero LAI over sparse veg, assign it to the most likely
            ! vegetation type.  Then we have to assign a cover fraction
            ! for the type, and updated cover fraction for now completely
            ! bare soil.

            !!!! do conversions !!!!

            ! ----------- The conditionals below are mutually exclusive.
            ! convert_vf(), when run, sets laic(esub%svm(BARE_SPARSE)), which will
            ! cause all subsequent conditional blocks to not run.

            ! --------------
            ! Convert to shrub if there's a small fraction and the shrubs already exist
            ! convert sparse veg to cold adapted shrub (9) if present
            if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. vfc(esub%svm(BARE_SPARSE)) < .15 &
               .and. laic(esub%svm(BARE_SPARSE)) > 0e0 &
               .and. vfc(esub%svm(COLD_SHRUB)) > 0e0 ) &
            then
                ! Preserve total LAI, but put it all in that vegetation
                ! type.
                call convert_vf( &
                    vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                    vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB)), &
                    laic(esub%svm(COLD_SHRUB)) )
            endif

            ! convert sparse veg to arid adapted shrub 10 if present
            if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. vfc(esub%svm(BARE_SPARSE)) < .15 &
               .and. laic(esub%svm(BARE_SPARSE)) > 0e0 &
               .and. vfc(esub%svm(ARID_SHRUB)) > 0e0 ) &
            then
             
                call convert_vf( &
                    vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                    vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB)), &
                    laic(esub%svm(ARID_SHRUB)) )
            end if
            ! --------------

            ! Convert to crops if they exist (any size fraction)
            ! convert sparse veg to crop 15 if present
            if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. laic(esub%svm(BARE_SPARSE)) > 0e0 &
                  .and. vfc(esub%crops_herb) > 0e0 ) &
            then
                call convert_vf( &
                    vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                    vfc(esub%crops_herb), laic(esub%crops_herb), &
                    laic(esub%crops_herb))
            end if
          
            ! Else... Convert to largest existing PFT
            ! convert sparse veg to pft with biggest fraction
            ! (if present)
            if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. laic(esub%svm(BARE_SPARSE)) > 0e0 ) then
                maxpft = maxloc( vfc(1:esub%last_pft), 1 )
                ! TODO: Why is cycle here?  Check Nancy's original code
                ! if ( vfc(maxpft) < .0001 ) cycle

                call convert_vf(vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                    vfc(maxpft), laic(maxpft), laic(maxpft))
            end if

            ! If none of the above was true, then:
            !  a) It's not a smal lfraction
            !  b) There's no arid shrubs or cold shrubs
            !  c) No crops
            !  d) No clear categorization of other vegetation types in the grid cell
            ! ...so just assign it to arid shrub.
            !
            ! convert sparse veg to arid adapted shrub 10
            if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. laic(esub%svm(BARE_SPARSE)) > 0e0 ) then
                call convert_vf(vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                    vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB)), 0e0 )
            end if
            ! -----------------------------------------------------------------

            ! Splits bare soil into bright and dark fractions, so we get
            ! the proper albedo in the GCM
            if (split_bare_soil) then
                vf_bare_sparse = vfc(esub%svm(BARE_SPARSE))
                vfc(esub%bare_bright) = vf_bare_sparse * bs_brightratio
                vfc(esub%bare_dark) = vf_bare_sparse - vfc(esub%bare_bright)
                laic(esub%bare_bright) = 0.
                laic(esub%bare_dark) = 0.
            end if


            do k=1,esub%ncover
                io_laiout(k)%buf(ic,jc) = laic(k)
            end do

            ! checksum lc & laimax
            io_sum_lc%buf(ic,jc) = 0.
            do k=1,esub%ncover
                io_sum_lc%buf(ic,jc) = io_sum_lc%buf(ic,jc) + laic(k)
            end do
        end do  ! ic
        end do  ! jc

        call chunker%write_chunks

    end do ! ichunk
    end do ! jchunk

    call chunker%close_chunks

end subroutine do_reindex
end module a04_mod

! ====================================================================

program convert
    use a04_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t) :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_reindex(esub)


end program convert
