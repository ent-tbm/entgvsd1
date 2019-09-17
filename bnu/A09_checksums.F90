module A09_mod

use chunker_mod
use ent_labels_mod
use geom_mod
use paths_mod
use gcm_labels_mod

implicit none
CONTAINS

subroutine ent_diff(chunker, &
    info1, &
    iroot2, idir2, ileaf2, vname2, &
    create_lr)

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
        LC_LAI_ENT_DIR, info1%dir, trim(info1%leaf)//'.nc', info1%vname, 1)

    call chunker%nc_open(xin2, &
        iroot2, idir2, ileaf2, vname2, 1)

    call chunker%nc_create(xout, weighting(chunker%wta1,1d0,0d0), &
        'checksum/'//trim(info1%dir), trim(info1%leaf)//'_diff', trim(info1%vname)//'_diff', &
        'checksum', '1', create_lr=create_lr)
    
    call chunker%nc_check('A09_checksums')  ! will get overwritten

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

            xout%buf(ic,jc) = v1-v2

        end do
        end do

        call chunker%write_chunks
        
    end do
    end do

end subroutine ent_diff

subroutine init_hr(chunker)
    type(Chunker_t) :: chunker
    call chunker%init(im1km,jm1km,imh,jmh,'forplot',2,1,1, nchunk=(/6,5/))
end subroutine


subroutine do_A09_checksums(esub_p)
    type(EntSet_t), pointer, intent(IN) :: esub_p
    ! --------- Locals
    type(FileInfo_t) :: info
    type(Chunker_t) :: chunker
    integer :: idoy, imonth

#if 0
    ! ----------- A01
    print *,'========================= A01'
    call init_hr(chunker)
    call chunker%file_info(info, ent20, &
        'BNU', 'M', 'lclaimax', 2004, 'ent17', '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'LAI3g') then
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax')
    else if (LAI_SOURCE == 'BNU') then
        call ent_diff(chunker, info, &
            LC_LAI_ENT_DIR, 'bnu/', 'bnu_laimax.nc', 'laimax')
    end if

    ! ------------ A01h
    print *,'========================= A01h'
    call init_hr(chunker)
    ! NOTE: Only makes sense for forest types
    call chunker%file_info(info, ent20, 'BNU', 'M', 'lchgt', 2004, 'ent17', '1.1', &
        varsuffix='_checksum')
    call ent_diff(chunker, info, &
        DATA_INPUT, 'height/', 'simard_forest_heights.nc', 'heights')


    ! ------------ A02
    print *,'========================= A02'
    do idoy = 1,ndoy
        call init_hr(chunker)
        call chunker%file_info(info, ent20, &
            'BNU', 'M', 'lclai', 2004, 'ent17', '1.1', &
            doytype='doy', idoy=idoy, &
            varsuffix = '_checksum')
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai')
    end do

    ! ----------- A03
    print *,'========================= A03'
    do imonth = 1,nmonth
        call init_hr(chunker)
        call chunker%file_info(info, ent20, 'BNU', 'M', 'lclai', 2004, 'ent17', '1.1', &
            doytype='month', idoy=imonth, varsuffix='_checksum')
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
    enddo

    ! ----------- A04
    print *,'========================= A04'
    call init_hr(chunker)
    call chunker%file_info(info, esub_p, &
        'BNU', 'M', 'lchgt', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    call ent_diff(chunker, info, &
        DATA_INPUT, 'height/', 'simard_forest_heights.nc', 'heights')
#endif

! Un-comment after we rerun
#if 0
    call init_hr(chunker)
    call chunker%file_info(info, esub_p, &
        'BNU', 'M', 'lclaimax', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'LAI3g') then
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax')
    else if (LAI_SOURCE == 'BNU') then
        call ent_diff(chunker, info, &
            LC_LAI_ENT_DIR, 'bnu/', 'bnu_laimax.nc', 'laimax')
    end if
#endif

    ! --------------- A05
    print *,'========================= A05'
    do idoy = 1,ndoy
        call init_hr(chunker)
        call chunker%file_info(info, esub_p, &
            'BNU', 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='doy', idoy=idoy)
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai')
    end do

    ! --------------- A06
    print *,'========================= A06'
    do imonth = 1,nmonth
        call init_hr(chunker)
        call chunker%file_info(info, esub_p, &
            'BNU', 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='month', idoy=imonth)

        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
    end do
end subroutine

subroutine init_lr(chunker)
    type(Chunker_t) :: chunker
    call chunker%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 2, 1, 1, (/1,1/))
end subroutine

subroutine do_A08_checksums(esub_p, step)
    class(EntSet_T), pointer :: esub_p
    character*(*), intent(IN) :: step
    ! --------- Locals
    type(FileInfo_t) :: info
    type(Chunker_t) :: chunker
    integer :: idoy, imonth


    ! --------------- A08_trim
    print *,'========================= A08 ',trim(step)
    call init_lr(chunker)
    call chunker%file_info(info, esub_p, 'BNU', 'M', 'lclaimax', 2004, step, '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'LAI3g') then
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax', create_lr=.false.)
    else if (LAI_SOURCE == 'BNU') then
        call ent_diff(chunker, info, &
            LC_LAI_ENT_DIR, 'bnu/', 'bnu_laimax.nc', 'laimax', create_lr=.false.)
    end if

    call init_lr(chunker)
    call chunker%file_info(info, esub_p, 'BNU', 'M', 'lchgt', 2004, step, '1.1', &
        varsuffix = '_checksum')
    call ent_diff(chunker, info, &
        DATA_INPUT, 'height/', 'simard_forest_heights.nc', 'heights', create_lr=.false.)

    do imonth=1,NMONTH
        call init_lr(chunker)
        call chunker%file_info(info, esub_p, 'BNU', 'M', 'lclai', 2004, step, '1.1', &
            doytype='month', idoy=imonth, varsuffix='_checksum')
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai', create_lr=.false.)
    end do

end subroutine do_A08_checksums

end module A09_mod


program A09_checksums
    use A09_mod
    use chunkparams_mod

implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub
    class(EntSet_t), pointer :: esub_p

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    esub_p => esub
!    call do_A09_checksums(esub_p)
    call do_A08_checksums(esub_p, 'trimmed')
!    call do_A08_checksums(esub_p, 'trimmed_scaled')
!    call do_A08_checksums(esub_p, 'trimmed_scaled_crops_ext')
!    call do_A08_checksums(esub_p, 'trimmed_scaled_nocrops')

end program A09_checksums
