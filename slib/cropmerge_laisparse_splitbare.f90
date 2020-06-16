module cropmerge_laisparse_splitbare_mod
    use chunker_mod
    use ent_params_mod
    use ent_labels_mod
    use gcm_labels_mod
    use conversions, only : convert_vf

implicit none
private

public cropmerge_laisparse_splitbare, Shrubtype
public split_bare_single

CONTAINS

! Allocates some LC and LAI to either COLD_SHRUB or ARID_SHRUB
! NOTE: This must be the same logic as set_shrubtype() in
!       B02_lc_modis_entpftrevcrop.F90
!
! MATEMP:
!     Mean Annaul Temperature [K]
! LC_IN:
!     Some LC fraction to be allocated
! LAI_IN:
!     LAI corresponding to LC_IN
! lc_cs, lai_cs:
!     LC/LAI variables to update for cold_shrub
! lc_as, lai_as:
!     LC/LAI variables to update for arid_shrub
!subroutine Set_Shrubtype(MATEMP,  LC_IN,LAI_IN, lc_cs,lai_cs,  lc_as,lai_as)
!implicit none
!    real*4, intent(IN) :: MATEMP,LC_IN, LAI_IN
!    real*4, intent(INOUT) :: lc_cs,lai_cs,  lc_as,lai_as
!    !------
!
!    if (MATEMP.lt.278.15) then !5 C cut-off
!        lai_cs = (lc_cs*lai_cs + LC_IN*LAI_IN) / (lc_cs + LC_IN)
!        lc_cs = lc_cs + LC_IN
!    else
!        lai_as = (lc_as*lai_as + LC_IN*LAI_IN) / (lc_as + LC_IN)
!        lc_as = lc_as + LC_IN
!    endif
!end subroutine Set_Shrubtype


integer function Shrubtype(MATEMP)
    real*4, intent(IN) :: MATEMP !,Pmave

    if (MATEMP.lt.278.15) then !5 C cut-off                                                                                                                 
       Shrubtype = COLD_SHRUB
    else
       Shrubtype = ARID_SHRUB
    endif
end function Shrubtype

! split_bare_single splits bare soil into bright and dark fractions for only
! ModelE legacy soil albedo scheme.  Also, due to the original soil albedo
! source file lacking some data in grid cells where MODIS land cover is
! non-zero, a check for undefined bs_brightratio is done:  if it is <0
! (undefined, should be FillValue, then it uses the most recent
! bs_brightratio_last if defined, otherwise assigns typical fractons of 0.3
! bright and 0.7 dark; if bs_brightratio is defined, then it bs_brightratio_last
! is updated.
subroutine split_bare_single( &
    esub, bs_brightratio, &
    vf_bare_sparse, &
    vfc,laic, &
    bs_brightratio_last)
implicit none
    type(GcmEntSet_t), intent(IN) :: esub
    real*4, intent(IN) :: bs_brightratio   !Fraction of bare that is bright.
    real*4, intent(IN) :: vf_bare_sparse
    real*4, intent(INOUT), dimension(esub%ncover) :: vfc,laic
    real*4, intent(INOUT) :: bs_brightratio_last
    !-----Local----
    integer :: m

    if ((split_bare_soil).and.(vf_bare_sparse>0.)) then
        if (bs_brightratio.lt.0.) then  !Unknown source of ~1742 points with bs_brightratio==FillValue
            if (bs_brightratio_last.eq.FillValue) then !No nearby previous okay bs_brightratio
                !print *, ic, jc, "bs_brightratio<0 ", bs_brightratio, vfc(esub%svm(BARE_SPARSE)), &
                print *, "bs_brightratio<0 ", bs_brightratio, vfc(esub%svm(BARE_SPARSE)), &
                    'Assigning 0.3 bright, 0.7 dark x BARE_SPARSE'
                    vfc(esub%bare_bright) = vf_bare_sparse * 0.3
                    vfc(esub%bare_dark) = vf_bare_sparse*(1. - 0.3) !vf_bare_sparse - vfc(esub%bare_bright)
            else
                !print *, ic, jc, "Using bs_brightratio_last",bs_brightratio_last, vf_bare_sparse
                print *, "Using bs_brightratio_last",bs_brightratio_last, vf_bare_sparse
                    vfc(esub%bare_bright) = vf_bare_sparse * bs_brightratio_last
                    vfc(esub%bare_dark) = vf_bare_sparse*max(0., 1. - bs_brightratio_last)
            endif
        else
                    vfc(esub%bare_bright) = vf_bare_sparse * bs_brightratio
                    vfc(esub%bare_dark) = vf_bare_sparse*(1. - bs_brightratio)  !vf_bare_sparse - vfc(esub%bare_bright)
                    bs_brightratio_last = bs_brightratio
        endif
            laic(esub%bare_bright) = 0.
            laic(esub%bare_dark) = 0.
     endif

end subroutine split_bare_single

! This is the common "guts" of B11_reclass_annual, B12_reclass_doy and B13_reclass_doy
! Arguments:
!
! esub:
!    Sub-universe, created with
!    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
! chunker:
!    Overall I/O control
! jc0,jc1,ic0,ic1:
!    Chunk indices, loop bounds
! ndoy:
!    Number of separate "days" or months to process
!    ndoy = (B11: 1 annual instance; B12: 2 specific days of year; B13: 12 months)
! jc0,jc1,ic0,ic1:
!    Chunk indices, loop bounds
! combine_crops_c3_c4:
!    Combine C3 and C4 crops into one PFT, for ModelE?
! split_bare_soil:
!    Split the bare soil into dark and light, to get the right albedo?
! ---------------- INPUT
! io_lcin(NENT20):
!    Land cover franctions (on master universe)
! io_laiin(NENT20, ndoy):
!    LAI (on master universe)
! io_bs:
!    Bare soil brightness ratio (V1km_bs_brightratio.nc)
! ---------------- OUTPUT
! io_laiout(esub%ncover, ndoy):
!    LAI (on sub universe, esub)
! -------------- INPUT (optional)
! io_simin(NENT20, ndoy)
!    (Simard) heights, produced from earlier step.
! -------------- OUTPUT (optional)
! sum_lc:
!    A single buffer, fill with sum_{land cover types} LC
!    To be used elsewhere by calling function.
! io_lclai_checksum(ndoy):
!    \sum_{land cover types} LC * LAI, for each of ndoy
! io_lc_checksum(ndoy):
!    \sum_{land cover types} LC, for each of ndoy; should be 1
! io_lcout(esub%ncover, ndoy)
!    Land cover fractions (on esub universe)
! io_simout(esub%ncover, ncdoy)
!    Heights, remapped from io_simin to esub universe
! io_lchgt_checksum(ndoy)
!    \sum_{land cover types} LC * height
subroutine cropmerge_laisparse_splitbare(esub, chunker, ndoy, &
    jc0,jc1, &
    ic0,ic1, &
    combine_crops_c3_c4, split_bare_soil, &
    io_lcin, io_laiin, io_bs, io_TCinave, &  ! INPUT
    io_laiout, &                              ! OUTPUT
    sum_lc, io_lclai_checksum, &
    io_lc_checksum, io_lcout, &                      ! OUTPUT (optional)
    io_simin, io_simout, io_lchgt_checksum)

    type(GcmEntSet_t) :: esub
    type(Chunker_t) :: chunker
    integer :: ndoy
    integer :: jc0,jc1, ic0,ic1
    logical :: combine_crops_c3_c4, split_bare_soil
    ! Input files
    type(ChunkIO_t) :: io_lcin(NENT20)
    type(ChunkIO_t) :: io_laiin(NENT20,ndoy)
    type(ChunkIO_t) :: io_bs
    type(ChunkIO_t) :: io_TCinave
    ! Output files
    type(ChunkIO_t) :: io_laiout(esub%ncover, ndoy)
    real*4, dimension(:,:), OPTIONAL :: sum_lc
    type(ChunkIO_t), OPTIONAL :: io_lclai_checksum(ndoy)
    type(ChunkIO_t), OPTIONAL :: io_lc_checksum(ndoy)
    type(ChunkIO_t), OPTIONAL :: io_lcout(esub%ncover, ndoy)
    type(ChunkIO_t), OPTIONAL :: io_simin(NENT20,ndoy)   ! Input
    type(ChunkIO_t), OPTIONAL :: io_simout(esub%ncover,ndoy)   ! Output
    type(ChunkIO_t), OPTIONAL :: io_lchgt_checksum(ndoy)

    ! -------------- Local Vars
    !These logical parameters could be passed in as parameters like split_bare_soil
    logical,parameter :: split_coasts=.false.   !Allocate all LAI to land portion of grid cell that has a water fraction<1. .false. in official release.
    logical,parameter :: land_ice_lai=.true.   !Assign any land ice with lai>0. to be arctic c3 grass.
    !---------------
    integer :: ichunk,jchunk,ic,jc,ii,jj
    integer :: i
    integer :: p
    integer :: k,ksub, ri
    real*8 :: CHECKSUM
    integer :: idoy

    ! Input values, max, monthly
    real*4 :: bs_brightratio   !soil albedo brightratio
    real*4 :: bs_brightratio_last !Save last good value of bs_brightratio, hack for some that are FillValue
    ! Renumbered input values
    real*4 vfc(esub%ncover)    ! = vfn = io_lcin
    real*4 laic(esub%ncover)   ! = lain = io_lcaiin

    real*4 area
    real*4 c3c4_crops
    real*4 waterice
    !real*4 am
    ! Converted values

    integer :: maxpft
    real*4 :: vf_bare_sparse
    real*4 :: laitot, lclailand, vf_land, f

    bs_brightratio_last = FillValue

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

                ! This operation leads to rare high LAI pixels. NOT USED IN OFFICIAL RELEASE. Only need to zero out water lai.
                ! Allocate coastal water lai to io_laiin land subgrid fractions (incl. SNOW_ICE and BARE_SPARSE).
                ! No change in cover fractions, just scale land lai to make
                ! sum(lc*lai) = laitot, and set water io_laiin to 0. 
                ! Later can convert any lai on SNOW_ICE oR BARE_SPARSE.
                if (split_coasts) then
                    vf_land = 0.
                    do  i=1,BARE_SPARSE
                        vf_land = vf_land + (io_lcin(i)%buf(ic,jc))
                    enddo
                    if (((vf_land.gt.0.).and.(io_lcin(CV_WATER)%buf(ic,jc).gt.0.)) &
                        .and.(io_laiin(CV_WATER,idoy)%buf(ic,jc).gt.0.)) then
                       laitot = 0.
                       do i=1,CV_WATER
                           laitot = laitot + (io_lcin(i)%buf(ic,jc))*io_laiin(i,idoy)%buf(ic,jc) !/sum(lc all)=1.
                       enddo
                       lclailand = 0.
                       do i=1,BARE_SPARSE
                           lclailand = lclailand + (io_lcin(i)%buf(ic,jc))*io_laiin(i,idoy)%buf(ic,jc) 
                       enddo
                       if (lclailand.gt.0) then
                           ! increase lai of all land fractions proportional to their existing lai, preserving laitot.
                           f = laitot / lclailand
                           do i=1,BARE_SPARSE
                               !print *, 'Scaling up lai on land from CV_WATER', &
                               !         ic,jc,f,vf_land,i, io_lcin(i)%buf(ic,jc), io_lcin(CV_WATER)%buf(ic,jc),&
                               !         io_laiin(CV_WATER,idoy)%buf(ic,jc), laitot, io_laiin(i,idoy)%buf(ic,jc)
                               io_laiin(i,idoy)%buf(ic,jc) = f * io_laiin(i,idoy)%buf(ic,jc)
                           enddo
                       else
                           do i=1,BARE_SPARSE
                               !This operation can be wrong if a green water
                               !body is surrounded by bare ground.
                               !print *,'Moving lai to land from CV_WATER', &
                               !        ic,jc,vf_land,i, io_lcin(i)%buf(ic,jc), io_lcin(CV_WATER)%buf(ic,jc),&
                               !        io_laiin(CV_WATER,idoy)%buf(ic,jc), laitot, io_laiin(i,idoy)%buf(ic,jc) 
                               io_laiin(i,idoy)%buf(ic,jc) = laitot / vf_land
                           enddo 
                       endif
                       io_laiin(CV_WATER,idoy)%buf(ic,jc) = 0.
                    endif
                endif

                ! Convert lai>0 on permanent ice to arctic c3 grass
                if (land_ice_lai) then
                   if ( (io_lcin(SNOW_ICE)%buf(ic,jc).gt.0.).and.(io_laiin(SNOW_ICE,idoy)%buf(ic,jc).gt.0.) ) then
                        print *, 'Converting lai on snow_ice to arctic c3 grass', &
                                ic,jc,io_lcin(SNOW_ICE)%buf(ic,jc), io_laiin(SNOW_ICE,idoy)%buf(ic,jc), &
                                io_lcin(C3_GRASS_ARCT)%buf(ic,jc), io_laiin(C3_GRASS_ARCT,idoy)%buf(ic,jc)
                        laitot = ( io_lcin(C3_GRASS_ARCT)%buf(ic,jc)*io_laiin(C3_GRASS_ARCT,idoy)%buf(ic,jc) + &
                                io_lcin(SNOW_ICE)%buf(ic,jc)*io_laiin(SNOW_ICE,idoy)%buf(ic,jc) ) / &
                                ( io_lcin(C3_GRASS_ARCT)%buf(ic,jc) + io_lcin(SNOW_ICE)%buf(ic,jc))
                        print *, '...laitot', laitot
                        io_laiin(C3_GRASS_ARCT,idoy)%buf(ic,jc) = laitot
                        io_laiin(SNOW_ICE,idoy)%buf(ic,jc) = 0.
                        if (present(io_simin)) then
                            print *, 'io_simin',io_simin(C3_GRASS_ARCT,idoy)%buf(ic,jc), io_simin(SNOW_ICE,idoy)%buf(ic,jc)
                            io_simin(C3_GRASS_ARCT,idoy)%buf(ic,jc) =  heights_form(2,C3_GRASS_ARCT)
                            io_simin(SNOW_ICE,idoy)%buf(ic,jc) =  0.
                        endif
                        io_lcin(C3_GRASS_ARCT)%buf(ic,jc) = io_lcin(C3_GRASS_ARCT)%buf(ic,jc) + io_lcin(SNOW_ICE)%buf(ic,jc)
                        io_lcin(SNOW_ICE)%buf(ic,jc) = 0.
                    endif
                endif 

                ! =============== Convert to GISS 16 pfts format with water_ice
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

                ! Combine water and permanent snow/ice. 
                waterice = 0.
!                waterice = io_lcin(CV_WATER)%buf(ic,jc) + io_lcin(SNOW_ICE)%buf(ic,jc)
                if (io_lcin(CV_WATER)%buf(ic,jc) > 0.) waterice = waterice + io_lcin(CV_WATER)%buf(ic,jc)
                if (io_lcin(SNOW_ICE)%buf(ic,jc) > 0.) waterice = waterice + io_lcin(SNOW_ICE)%buf(ic,jc)
                if (waterice > 0.) then
                     laic(esub%water_ice) = ( & 
                        io_lcin(CV_WATER)%buf(ic,jc) *io_laiin(CV_WATER,idoy)%buf(ic,jc) &
                      + io_lcin(SNOW_ICE)%buf(ic,jc) *io_laiin(SNOW_ICE,idoy)%buf(ic,jc) &
                            ) / waterice
                     vfc(esub%water_ice) = waterice

                     if (present(io_simout)) then
                        io_simout(esub%water_ice,idoy)%buf(ic,jc) = ( & 
                            io_lcin(CV_WATER)%buf(ic,jc) * io_simin(CV_WATER,idoy)%buf(ic,jc) + &
                            io_lcin(SNOW_ICE)%buf(ic,jc) * io_simin(SNOW_ICE,idoy)%buf(ic,jc) &
                         ) / waterice
                     end if
                     if (present(io_simout)) then
                        io_simout(esub%water_ice,idoy)%buf(ic,jc) = 0.
                     end if

                else   ! waterice <= 0
                        laic(esub%water_ice) = 0.
                        vfc(esub%water_ice) = 0.
                        if (present(io_simout)) then
                            io_simout(esub%water_ice,idoy)%buf(ic,jc) = FillValue
                        end if
                end if


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

                        if (present(io_simout)) then
                            io_simout(esub%crops_herb,idoy)%buf(ic,jc) = ( &
                                io_lcin(CROPS_C3_HERB)%buf(ic,jc) * io_simin(CROPS_C3_HERB,idoy)%buf(ic,jc) + &
                                io_lcin(CROPS_C4_HERB)%buf(ic,jc) * io_simin(CROPS_C4_HERB,idoy)%buf(ic,jc) &
                                ) / c3c4_crops

                        end if
                    else   ! c3c4_crops <= 0
                        laic(esub%crops_herb) = 0.
                        vfc(esub%crops_herb) = 0.
                        if (present(io_simout)) then
                            io_simout(esub%crops_herb,idoy)%buf(ic,jc) = FillValue
                        end if
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

                 maxpft = maxloc( vfc(1:esub%last_pft), 1 )

                ! --------------
                ! Convert to shrub if there's a small fraction and the shrubs already exist
                ! convert sparse veg to cold adapted shrub (9) if present
                ! At coarser spatial resolutions, use the upper limit of BARE_SPARSE fraction to convert.
                if( vfc(esub%svm(BARE_SPARSE)) > 0. .and. vfc(esub%svm(BARE_SPARSE)) < .85 &
                   .and. laic(esub%svm(BARE_SPARSE)) > 0.) &
                then
                   if ( vfc(esub%svm(COLD_SHRUB)) > 0. )  then
                    ! Preserve total LAI, but put it all in that vegetation
                    ! type.
                      if (laic(esub%svm(COLD_SHRUB)).le.0.) then
                         print *, 'ERROR: no cold_shrub LAI for cold_shrub cover', ii,jj
                      else
                         !print *, ii,jj,'Converting to cold shrub', &
                         !       vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                         !       vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB))
                          call convert_vf( &
                             vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                             vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB)), &
                             laic(esub%svm(COLD_SHRUB)) )
                      endif
!                endif

                ! convert sparse veg to arid adapted shrub 10 if present
!                if( vfc(esub%svm(BARE_SPARSE)) > 0e0 & ! .and. vfc(esub%svm(BARE_SPARSE)) < .15 &
!                   .and. laic(esub%svm(BARE_SPARSE)) > 0e0 &
!                   .and. vfc(esub%svm(ARID_SHRUB)) > 0e0 ) &
                   elseif (vfc(esub%svm(ARID_SHRUB)) > 0e0 )  &
                   then
                        if (laic(esub%svm(ARID_SHRUB)).le.0.) then
                           print *, 'ERROR: no arid_shrub LAI for arid_shrub cover',ii,jj 
                        else 
                           !print *, ii,jj,'Converting to arid shrub', &
                           !     vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                           !     vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB))
                           call convert_vf( &
                             vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                             vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB)), &
                             laic(esub%svm(ARID_SHRUB)) )
                        endif

!                end if

                ! --------------
                ! Convert to crops if they exist (any size fraction)
                ! convert the rest of sparse veg to crop 15 if present
!                if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. laic(esub%svm(BARE_SPARSE)) > 0e0 &
!                      .and. vfc(esub%crops_herb) > 0e0 ) &
                   elseif (vfc(esub%crops_herb) > 0e0 ) &
                   then
                        if (laic(esub%crops_herb).le.0.) then
                            print *, 'ERROR: no crops_herb LAI for crops_herb cover',ii,jj
                        else
                           !print *, 'Converting to crops', &
                           !  vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                           !  vfc(esub%crops_herb), laic(esub%crops_herb)

                           call convert_vf( &
                             vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                             vfc(esub%crops_herb), laic(esub%crops_herb), &
                             laic(esub%crops_herb))
                        endif
!                end if

                ! Else... Convert to largest existing PFT
                ! convert sparse veg to pft with biggest fraction
                ! (if present)
!                if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. laic(esub%svm(BARE_SPARSE)) > 0e0 ) then

                   elseif (( vfc(esub%svm(maxpft)) >= 0. ).and.(laic(esub%svm(maxpft)).gt.0.)) &
                   then
                    !maxpft = maxloc( vfc(1:esub%last_pft), 1 )  !Moved to top
                    !of conditionals
                    ! TODO: Why is cycle here?  Check Nancy's original code
                    ! Original code has cycle, which I don't think was correct: because it would
                    !    have skipped putting the variables back into arrays to be written,
                    !    as well as split_bare_soil (see original A07.f)
                    ! if ( vfc(maxpft) < .0001 ) cycle
                         
                     !    if ((vfc(esub%svm(maxpft)).gt.0.0).and.(laic(esub%svm(maxpft)).gt.0.0)) then
                             call convert_vf( &
                                vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                                vfc(esub%svm(maxpft)), laic(esub%svm(maxpft)), laic(esub%svm(maxpft)))

                ! If none of the above was true, then:
                !  a) It's not a small fraction
                !  b) There's no arid shrubs or cold shrubs
                !  c) No crops
                !  d) No clear categorization of other vegetation types in the grid cell
                ! Assign to cold_shrub or arid_shrub based on temperature. Set minimum LAI of veg fraction to 0.5 to keep some bare fraction.
!                if( vfc(esub%svm(BARE_SPARSE)) > 0e0 .and. laic(esub%svm(BARE_SPARSE)) > 0e0 ) then
                   elseif (Shrubtype(io_TCinave%buf(ic,jc)+273.15).eq.COLD_SHRUB) then
                              call convert_vf(vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                                 vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB)), max(0.5,laic(esub%svm(BARE_SPARSE))))
                   else
                              call convert_vf(vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                                 vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB)), max(0.5,laic(esub%svm(BARE_SPARSE))))
                   endif

                   !call Set_Shrubtype( &
                   !     io_TCinave%buf(ic,jc)+273.15, &
                   !     vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                   !     vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB)), &
                   !     vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB)))
                     !else
                     !   print *,ii,jj, 'maxpft',maxpft, 'vfc',  vfc(esub%svm(maxpft)), 'laic',laic(esub%svm(maxpft))
                     !end if
                else ! vfc BARE_SPARSE>=0.85
                     if  ((vfc(esub%svm(BARE_SPARSE)).ge.0.85) .and. laic(esub%svm(BARE_SPARSE)) > 0.) then
                        !print *, 'vfc BARE >= 0.85',vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                        !        vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB)), &
                        !        vfc(esub%svm(ARID_SHRUB)), laic(esub%svm(ARID_SHRUB))
                        if (Shrubtype(io_TCinave%buf(ic,jc)+273.15).eq.COLD_SHRUB) then
                            call convert_vf(vfc(esub%svm(BARE_SPARSE)), laic(esub%svm(BARE_SPARSE)), &
                                vfc(esub%svm(COLD_SHRUB)), laic(esub%svm(COLD_SHRUB)), & 
                                max(0.5, laic(esub%svm(BARE_SPARSE))))
                        else
                            call convert_vf(vfc(esub%svm(BARE_SPARSE)), & 
                               laic(esub%svm(BARE_SPARSE)), vfc(esub%svm(ARID_SHRUB)), &
                                laic(esub%svm(ARID_SHRUB)), max(0.5, laic(esub%svm(BARE_SPARSE))))
                        endif
                     else
                        !print *, ii,jj, 'No lai for BARE_SPARSE>=0.85'
                     endif
                end if !if vfc and lai BARE_SPARSE are >0. AND vfc BARE_SPASE < 0.85

                ! -----------------------------------------------------------------
                if (present(io_simout)) then
                    if (vfc(esub%svm(ARID_SHRUB)) > 0 .and. &
                        io_simout(esub%svm(ARID_SHRUB),idoy)%buf(ic,jc) == FillValue) then

                        io_simout(esub%svm(ARID_SHRUB),idoy)%buf(ic,jc) = heights_form(2,ARID_SHRUB)
                    end if

                    if (vfc(esub%svm(COLD_SHRUB)) > 0 .and. &
                        io_simout(esub%svm(COLD_SHRUB),idoy)%buf(ic,jc) == FillValue) then

                        io_simout(esub%svm(COLD_SHRUB),idoy)%buf(ic,jc) = heights_form(2,COLD_SHRUB)
                    end if
                end if

                ! Splits bare soil into bright and dark fractions, only for
                ! legacy GISS ModelE scheme for calculating albedo.
                ! (We need a soil albedo wherever there is veg or bare soil type)
                ! The soil albedo grid fill interpolates along latitudes.
                ! Carrer uses same dataset, so things should line up; however,
                ! MODIS lc still has pixels at high latitude not in the Carrer 6
                ! km data set; therefore, checked for in split_bare_single.
                vf_bare_sparse = vfc(esub%svm(BARE_SPARSE))
                call split_bare_single( esub, bs_brightratio, &
                    vf_bare_sparse, &
                    vfc,laic, &
                    bs_brightratio_last)
                !if ((split_bare_soil).and.(vfc(esub%svm(BARE_SPARSE))>0.)) then
                !    vf_bare_sparse = vfc(esub%svm(BARE_SPARSE))
                !    if (bs_brightratio.lt.0.) then  !Unknown source of ~1742 points with bs_brightratio==FillValue
                !        if (bs_brightratio_last.eq.FillValue) then !No nearby previous okay bs_brightratio
                !           print *, ic, jc, "bs_brightratio<0 ", bs_brightratio, vfc(esub%svm(BARE_SPARSE)), &
                !                'Assigning 0.3 bright, 0.7 dark x BARE_SPARSE'
                !           vfc(esub%bare_bright) = vf_bare_sparse * 0.3
                !           vfc(esub%bare_dark) = vf_bare_sparse*(1. - 0.3) !vf_bare_sparse - vfc(esub%bare_bright)
                !        else
                !           print *, ic, jc, "Using bs_brightratio_last",bs_brightratio_last, vf_bare_sparse 
                !           vfc(esub%bare_bright) = vf_bare_sparse * bs_brightratio_last
                !           vfc(esub%bare_dark) = vf_bare_sparse*max(0., 1. - bs_brightratio_last) 
                !        endif
                !    else
                !        vfc(esub%bare_bright) = vf_bare_sparse * bs_brightratio
                !        vfc(esub%bare_dark) = vf_bare_sparse*(1. - bs_brightratio)  !vf_bare_sparse - vfc(esub%bare_bright)
                !        bs_brightratio_last = bs_brightratio
                !    endif
                !    laic(esub%bare_bright) = 0.
                !    laic(esub%bare_dark) = 0.
                !end if


                do k=1,esub%ncover
                    io_laiout(k,idoy)%buf(ic,jc) = laic(k)
                end do
                if (present(io_lcout)) then
                    do k=1,esub%ncover
                        io_lcout(k,idoy)%buf(ic,jc) = vfc(k)
                    end do
                end if

                ! checksum lc & laimax

                ! ============= Zero stuff out
                if (present(io_lc_checksum)) io_lc_checksum(idoy)%buf(ic,jc) = 0
                if (present(io_lclai_checksum)) io_lclai_checksum(idoy)%buf(ic,jc) = 0
                if (present(io_lchgt_checksum)) io_lchgt_checksum(idoy)%buf(ic,jc) = 0

                if (present(io_lc_checksum)) then
                    io_lc_checksum(idoy)%buf(ic,jc) = 0.
                    do k=1,esub%ncover
                        io_lc_checksum(idoy)%buf(ic,jc) = io_lc_checksum(idoy)%buf(ic,jc) + vfc(k)
                    end do
!## NK DEBUG ##
!                    if (io_lc_checksum(idoy)%buf(ic,jc).gt.1.) then
!                       print *, 'io',ic,jc, io_lc_checksum(idoy)%buf(ic,jc), vfc(:)
!                    endif
                end if

                ! Compute weighting for checksums
                if (present(sum_lc)) then
                    sum_lc(ic,jc) = 0
                    do k = 1,esub%ncover
                        sum_lc(ic,jc) = sum_lc(ic,jc) + vfc(k)
                    end do
!## NK DEBUG ##
!                    if (sum_lc(ic,jc).gt.1.) then
!                       print *, 'sum_lc',ic,jc, sum_lc(ic,jc), vfc(:)
!                    endif
                end if

                if (present(io_lclai_checksum)) then
                    io_lclai_checksum(idoy)%buf(ic,jc) = 0.
                    do k=1,esub%ncover
                        io_lclai_checksum(idoy)%buf(ic,jc) = io_lclai_checksum(idoy)%buf(ic,jc) + &
                            vfc(k) * laic(k)
                    end do
                end if

                if (present(io_lchgt_checksum)) then
                  if (present(io_simout)) then
                    io_lchgt_checksum(idoy)%buf(ic,jc) = 0.
                    do k=1,esub%ncover
                        io_lchgt_checksum(idoy)%buf(ic,jc) = io_lchgt_checksum(idoy)%buf(ic,jc) + &
                            vfc(k) * io_simout(k,idoy)%buf(ic,jc) 
                    end do
                  endif
                end if

            end do    ! idoy

        end do
        end do

        call chunker%write_chunks

    end do
    end do

end subroutine cropmerge_laisparse_splitbare

end module cropmerge_laisparse_splitbare_mod
