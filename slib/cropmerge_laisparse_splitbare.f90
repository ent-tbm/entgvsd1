module cropmerge_laisparse_splitbare_mod
    use chunker_mod
    use chunkparams_mod
    use ent_labels_mod
    use gcm_labels_mod
    use conversions, only : convert_vf

implicit none
CONTAINS

subroutine cropmerge_laisparse_splitbare(esub, chunker, ndoy, &
    jc0,jc1, &
    ic0,ic1, &
    combine_crops_c3_c4, split_bare_soil, &
    io_lcin, io_laiin, io_bs, &               ! INPUT
    io_laiout, &                              ! OUTPUT
    io_sum_lc, io_lcout, &                      ! OUTPUT (optional)
    io_simin, io_simout)

    type(GcmEntSet_t) :: esub
    type(Chunker_t) :: chunker
    integer :: ndoy
    integer :: jc0,jc1, ic0,ic1
    logical :: combine_crops_c3_c4, split_bare_soil
    ! Input files
    type(ChunkIO_t) :: io_lcin(NENT20)
    type(ChunkIO_t) :: io_laiin(NENT20,ndoy)
    type(ChunkIO_t) :: io_bs
    ! Output files
    type(ChunkIO_t) :: io_laiout(esub%ncover, ndoy)
    type(ChunkIO_t), OPTIONAL :: io_sum_lc(ndoy)
    type(ChunkIO_t), OPTIONAL :: io_lcout(esub%ncover, ndoy)
    type(ChunkIO_t), OPTIONAL :: io_simin(NENT20,ndoy)   ! Input
    type(ChunkIO_t), OPTIONAL :: io_simout(esub%ncover,ndoy)   ! Output

    ! -------------- Local Vars
    integer :: ichunk,jchunk,ic,jc,ii,jj
    integer :: i
    integer :: k,ksub, ri
    real*8 :: CHECKSUM
    integer :: idoy

    ! Input values, max, monthly
    real*4 :: bs_brightratio   !are soil brightratio
    ! Renumbered input values
    real*4 vfc(esub%ncover)    ! = vfn = io_lcin
    real*4 laic(esub%ncover)   ! = lain = io_lcaiin

    real*4 area
    real*4 c3c4_crops
    real*4 am
    ! Converted values

    integer :: maxpft
    real*4 :: vf_bare_sparse

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

            do idoy=1,ndoy

                ! get bs bright ratio
                bs_brightratio = io_bs%buf(ic,jc)

                ! Check if mismatch lc or lai values (one is zero and the other not)
                ! call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)


                ! Clear per-grid-cell buffers
                vfc = 0
                laic = 0

                ! =============== Convert to GISS 16 pfts format
                ! First simple transfers
                do ri=1,esub%nremap
                    ksub = esub%remap(1,ri)    ! Index in esub
                    k = esub%remap(2,ri)       ! Index in ent20

                    vfc(ksub) = io_lcin(k)%buf(ic,jc)
                    laic(ksub) = io_laiin(k,idoy)%buf(ic,jc)

                    if (present(io_simout)) then
                        io_simout(ksub,1)%buf(ic,jc) = io_simin(k,1)%buf(ic,jc)
                    end if
                end do

                ! Then more complex stuff
                if (combine_crops_c3_c4) then
                    !lc laimax
                    c3c4_crops = io_lcin(CROPS_C3_HERB)%buf(ic,jc) + io_lcin(CROPS_C4_HERB)%buf(ic,jc)
                    if (c3c4_crops > 0.) then
                        laic(esub%crops_herb) = ( &
                              io_lcin(CROPS_C3_HERB)%buf(ic,jc) *io_laiin(CROPS_C3_HERB,idoy)%buf(ic,jc) &
                            + io_lcin(CROPS_C4_HERB)%buf(ic,jc) *io_laiin(CROPS_C4_HERB,idoy)%buf(ic,jc) &
                            ) / c3c4_crops
                        vfc(esub%crops_herb) = c3c4_crops
                    else
                        laic(esub%crops_herb) = 0.
                        vfc(esub%crops_herb) = 0.
                    end if

                    if (present(io_simout)) then
                        io_simout(esub%crops_herb,idoy)%buf(ic,jc) = ( &
                            io_lcin(CROPS_C3_HERB)%buf(ic,jc) * io_simin(CROPS_C3_HERB,idoy)%buf(ic,jc) + &
                            io_lcin(CROPS_C4_HERB)%buf(ic,jc) * io_simin(CROPS_C4_HERB,idoy)%buf(ic,jc) &
                            ) / c3c4_crops
                    end if
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
                    ! Original code has cycle, which I don't think was correct: because it would
                    !    have skipped putting the variables back into arrays to be written,
                    !    as well as split_bare_soil (see original A07.f)
                    ! if ( vfc(maxpft) < .0001 ) cycle
                    if ( vfc(maxpft) >= .0001 ) then
                        call convert_vf(vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                            vfc(maxpft), laic(maxpft), laic(maxpft))
                    end if
                end if

                ! If none of the above was true, then:
                !  a) It's not a small fraction
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
                    io_laiout(k,idoy)%buf(ic,jc) = laic(k)
                end do
                if (present(io_lcout)) then
                    do k=1,esub%ncover
                        io_lcout(k,idoy)%buf(ic,jc) = vfc(k)
                    end do
                end if

                ! checksum lc & laimax
                if (present(io_sum_lc)) then
                    io_sum_lc(idoy)%buf(ic,jc) = 0.
                    do k=1,esub%ncover
                        io_sum_lc(idoy)%buf(ic,jc) = io_sum_lc(idoy)%buf(ic,jc) + laic(k)
                    end do
                end if

            end do    ! idoy

        end do
        end do

        call chunker%write_chunks

    end do
    end do

end subroutine cropmerge_laisparse_splitbare

end module cropmerge_laisparse_splitbare_mod
