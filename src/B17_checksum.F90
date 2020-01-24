#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module B17_mod
!Author:  Elizabeth Fischer
!
! Calculates miscellaneous checksums: error checks, differences from
! observations, sums of lc, sums over PFTs or lc for grid totals, etc.

use chunker_mod
use ent_labels_mod
use gcm_labels_mod
use ent_params_mod
use hntr_mod

implicit none
CONTAINS

subroutine ent_diff(rw, chunker, &
    info1, &
    iroot2, idir2, ileaf2, vname2, &
    create_lr)

    type(ReadWrites_t) :: rw
    type(Chunker_t) :: chunker
    integer :: im,jm
    type(FileInfo_t) :: info1
    character*(*), intent(IN) :: iroot2,idir2,ileaf2,vname2
    logical, intent(IN), OPTIONAL :: create_lr
    real*4 :: v1,v2

    ! -------------- Locals
    type(ChunkIO_t) :: xin1, xin2, xout
    integer :: ichunk,jchunk, ic,jc

    call chunker%nc_open(xin1, &
        chunker%outputs_dir, info1%dir, trim(info1%leaf)//'.nc', info1%vname, 1)

    call chunker%nc_open(xin2, &
        iroot2, idir2, ileaf2, vname2, 1)

    call chunker%nc_create(xout, weighting(chunker%wta1,1d0,0d0), &
        'checksum/'//trim(info1%dir), trim(info1%leaf)//'_diff', trim(info1%vname)//'_diff', &
        'checksum', '1', create_lr=create_lr)
    
    call chunker%nc_check(rw=rw)

#ifdef JUST_DEPENDENCIES
    return
#endif

#ifdef ENTGVSD_DEBUG
    do jchunk = 1,1
    do ichunk = 1,1
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif
        call chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            v1 = xin1%buf(ic,jc)
            v2 = xin2%buf(ic,jc)

            if ((v1/=v1).and.(v2/=v2)) then
                xout%buf(ic,jc) = 0.   ! They match!
            else
                xout%buf(ic,jc) = v1-v2
            end if

        end do
        end do

        call chunker%write_chunks
        
    end do
    end do

end subroutine ent_diff




subroutine init_hr(chunker)
    type(Chunker_t) :: chunker
    call chunker%init(im1km,jm1km,imh,jmh,'forplot',2,1,1, nchunk=(/1,15/), outputs_dir=THIS_OUTPUTS_DIR)
end subroutine


subroutine do_B17_checksums(rw, esub_p)
    type(ReadWrites_t) :: rw
    type(EntSet_t), pointer, intent(IN) :: esub_p
    ! --------- Locals
    type(FileInfo_t) :: info
    type(Chunker_t) :: chunker
    integer :: idoy, imonth

#if 1
    ! ----------- B08
    print *,'========================= B08_lc_laimax'
    call init_hr(chunker)
    call chunker%file_info(info, ent20, &
        LAI_SOURCE, 'M', 'lclaimax', 2004, 'ent17', '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'L') then
        call ent_diff(rw, chunker, info, &
            INPUTS_DIR, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax')
    else if (LAI_SOURCE == 'B') then
        call ent_diff(rw, chunker, info, &
            chunker%outputs_dir, 'tmp/bnu/', 'bnu_laimax.nc', 'laimax')
    end if

    ! ------------ B04
    print *,'========================= B04_veg_height'
    call init_hr(chunker)
    ! NOTE: Only makes sense for forest types
    call chunker%file_info(info, ent20, LAI_SOURCE, 'M', 'lchgt', 2004, 'ent17', '1.1', &
        varsuffix='_checksum')
    call ent_diff(rw, chunker, info, &
        INPUTS_DIR, 'height/', 'simard_forest_heights.nc', 'heights')


    ! ------------ B09
    print *,'========================= B09_lc_lai_doy'
    do idoy = 1,ndoy
        call init_hr(chunker)
        call chunker%file_info(info, ent20, &
            LAI_SOURCE, 'M', 'lclai', 2004, 'ent17', '1.1', &
            doytype='doy', idoy=idoy, &
            varsuffix = '_checksum')
        call ent_diff(rw, chunker, info, &
            LAI_INPUTS_DIR, 'lai/BNU/doy/2004/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai')
    end do

    ! ----------- B10
    print *,'========================= B10_lc_lai_monthly'
    do imonth = 1,nmonth
        call init_hr(chunker)
        call chunker%file_info(info, ent20, LAI_SOURCE, 'M', 'lclai', 2004, 'ent17', '1.1', &
            doytype='month', idoy=imonth, varsuffix='_checksum')
        call ent_diff(rw, chunker, info, &
            LAI_INPUTS_DIR, 'lai/BNU/monthly/2004/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
    enddo

    ! ----------- B11
    print *,'========================= B11_reclass_annual'
    call init_hr(chunker)
    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lchgt', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    call ent_diff(rw, chunker, info, &
        INPUTS_DIR, 'height/', 'simard_forest_heights.nc', 'heights')
#endif

! Un-comment after we rerun
#if 1
    call init_hr(chunker)
    call chunker%file_info(info, esub_p, &
        LAI_SOURCE, 'M', 'lclaimax', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'L') then
        call ent_diff(rw, chunker, info, &
            INPUTS_DIR, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax')
    else if (LAI_SOURCE == LAI_SOURCE) then
        call ent_diff(rw, chunker, info, &
            chunker%outputs_dir, 'tmp/bnu/', 'bnu_laimax.nc', 'laimax')
    end if
#endif

    ! --------------- B12
    print *,'========================= B12_reclass_doy'
    do idoy = 1,ndoy
        call init_hr(chunker)
        call chunker%file_info(info, esub_p, &
            LAI_SOURCE, 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='doy', idoy=idoy)
        call ent_diff(rw, chunker, info, &
            LAI_INPUTS_DIR, 'lai/BNU/doy/2004/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai')
    end do

    ! --------------- B13
    print *,'========================= B13_reclass_monthly'
    do imonth = 1,nmonth
        call init_hr(chunker)
        call chunker%file_info(info, esub_p, &
            LAI_SOURCE, 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='month', idoy=imonth)

        call ent_diff(rw, chunker, info, &
            LAI_INPUTS_DIR, 'lai/BNU/monthly/2004/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
    end do
end subroutine

subroutine init_lr(chunker)
    type(Chunker_t) :: chunker
    call chunker%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 2, 1, 1, nchunk=(/1,1/), outputs_dir=THIS_OUTPUTS_DIR)
end subroutine

subroutine do_B16_checksums(rw, esub_p, step)
    type(ReadWrites_t) :: rw
    class(EntSet_T), pointer :: esub_p
    character*(*), intent(IN) :: step
    ! --------- Locals
    type(FileInfo_t) :: info
    type(Chunker_t) :: chunker
    integer :: idoy, imonth


    ! --------------- B16_trim
    print *,'========================= B16 ',trim(step)
    call init_lr(chunker)
    call chunker%file_info(info, esub_p, LAI_SOURCE, 'M', 'lclaimax', 2004, step, '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'L') then
        call ent_diff(rw, chunker, info, &
            chunker%outputs_dir, 'tmp/regrids/', 'LAI3gMax_1kmx1km_hxh.nc', 'laimax', &
            create_lr=.false.)
    else if (LAI_SOURCE == 'B') then
        call ent_diff(rw, chunker, info, &
            chunker%outputs_dir, 'tmp/regrids/', 'bnu_laimax_hxh.nc', 'laimax', create_lr=.false.)
    end if


    call init_lr(chunker)
    call chunker%file_info(info, esub_p, LAI_SOURCE, 'M', 'lchgt', 2004, step, '1.1', &
        varsuffix = '_checksum')
    call ent_diff(rw, chunker, info, &
        chunker%outputs_dir, 'tmp/regrids/', 'simard_forest_heights_hxh.nc', 'heights', &
        create_lr=.false.)



    do imonth=1,NMONTH
        call init_lr(chunker)
        call chunker%file_info(info, esub_p, LAI_SOURCE, 'M', 'lclai', 2004, step, '1.1', &
            doytype='month', idoy=imonth, varsuffix='_checksum')
        call ent_diff(rw, chunker, info, &
            chunker%outputs_dir, 'tmp/regrids/', 'global_30s_2004_'//MONTH(imonth)//'_hxh.nc', &
            'lai', create_lr=.false.)
    end do

end subroutine do_B16_checksums

end module B17_mod


program B17_checksums
    use B17_mod
    use ent_params_mod

implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub
    class(EntSet_t), pointer :: esub_p
    type(ReadWrites_t) :: rw
    call rw%init(THIS_OUTPUTS_DIR, 'B17_checksum', 300,300)

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    esub_p => esub
    call do_B17_checksums(rw, esub_p)
    call do_B16_checksums(rw, esub_p, 'trimmed')
    call do_B16_checksums(rw, esub_p, 'trimmed_scaled')
    call do_B16_checksums(rw, esub_p, 'trimmed_scaled_crops_ext')
    call do_B16_checksums(rw, esub_p, 'trimmed_scaled_nocrops')

    call rw%write_mk
end program B17_checksums
