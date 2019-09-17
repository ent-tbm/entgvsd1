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

    ! -------------- Locals
    type(ChunkIO_t) :: xin1, xin2, xout
    integer :: ichunk,jchunk, ic,jc

    call chunker%nc_open(xin1, &
        LC_LAI_ENT_DIR, info1%dir, trim(info1%leaf)//'.nc', info1%vname, 1)

    call chunker%nc_open(xin2, &
        iroot2, idir2, ileaf2, vname2, 1)

    call chunker%nc_create(xout, weighting(chunker%wta1,1d0,0d0), &
        'checksum/', trim(info1%leaf)//'_diff', trim(info1%vname)//'_diff', &
        'checksum', '1', create_lr=create_lr)
    
    call chunker%nc_check('A09_checksums')  ! will get overwritten


    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)

        call chunker%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            xout%buf(ic,jc) = xin1%buf(ic,jc) - xin2%buf(ic,jc)

        end do
        end do

        call chunker%write_chunks
        
    end do
    end do

end subroutine ent_diff


subroutine do_A09_checksums(esub_p)
    type(EntSet_t), pointer, intent(IN) :: esub_p
    ! --------- Locals
    type(FileInfo_t) :: info
    type(Chunker_t) :: chunker
    integer :: idoy, imonth
    call chunker%init(im1km,jm1km,imh,jmh,'forplot',2,1,1, nchunk=(/6,5/))


    ! ----------- A01
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
    ! NOTE: Only makes sense for forest types
    call chunker%file_info(info, ent20, 'BNU', 'M', 'lchgt', 2004, 'ent17', '1.1', &
        varsuffix='_checksum')
    call ent_diff(chunker, info, &
        DATA_INPUT, 'height/', 'simard_forest_heights.nc', 'heights')


    ! ------------ A02
    do idoy = 1,ndoy
        call chunker%file_info(info, ent20, &
            'BNU', 'M', 'lclai', 2004, 'ent17', '1.1', &
            doytype='doy', idoy=idoy, &
            varsuffix = '_checksum')
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai')
    end do

    ! ----------- A03
    do imonth = 1,nmonth
        call chunker%file_info(info, ent20, 'BNU', 'M', 'lclai', 2004, 'ent17', '1.1', &
            doytype='month', idoy=imonth, varsuffix='_checksum')
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
    enddo

    ! ----------- A04
    call chunker%file_info(info, esub_p, &
        'BNU', 'M', 'lchgt', 2004, 'pure', '1.1', &
        varsuffix = '_checksum')
    call ent_diff(chunker, info, &
        DATA_INPUT, 'height/', 'simard_forest_heights.nc', 'heights')

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

    ! --------------- A05
    do idoy = 1,ndoy
        call chunker%file_info(info, esub_p, &
            'BNU', 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='doy', idoy=idoy)
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'global_30s_2004_'//DOY(idoy)//'.nc', 'lai')
    end do

    ! --------------- A06
    do imonth = 1,nmonth
        call chunker%file_info(info, esub_p, &
            'BNU', 'M', 'lclai', 2004, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='month', idoy=imonth)

        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
    end do
end subroutine

subroutine do_A08_checksums(esub_p, step)
    class(EntSet_T), pointer :: esub_p
    character*(*), intent(IN) :: step
    ! --------- Locals
    type(FileInfo_t) :: info
    type(Chunker_t) :: chunker
    integer :: idoy, imonth

    call chunker%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 2, 1, 1, (/1,1/))

    ! --------------- A08_trim
    call chunker%file_info(info, esub_p, 'BNU', 'M', 'lclaimax', 2004, step, '1.1', &
        varsuffix = '_checksum')
    if (LAI_SOURCE == 'LAI3g') then
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/', 'LAI3gMax_1kmx1km.nc', 'laimax', create_lr=.false.)
    else if (LAI_SOURCE == 'BNU') then
        call ent_diff(chunker, info, &
            LC_LAI_ENT_DIR, 'bnu/', 'bnu_laimax.nc', 'laimax', create_lr=.false.)
    end if

    call chunker%file_info(info, esub_p, 'BNU', 'M', 'lchgt', 2004, step, '1.1', &
        varsuffix = '_checksum')
    call ent_diff(chunker, info, &
        DATA_INPUT, 'height/', 'simard_forest_heights.nc', 'heights', create_lr=.false.)

    do imonth=1,NMONTH
        call chunker%file_info(info, esub_p, 'BNU', 'M', 'lclai', 2004, step, '1.1', &
            doytype='month', idoy=imonth, varsuffix='_checksum')
        call ent_diff(chunker, info, &
            DATA_INPUT, 'LAI/BNUMonthly/', 'global_30s_2004_'//MONTH(imonth)//'.nc', &
            'lai')
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
    call do_A09_checksums(esub_p)
    call do_A08_checksums(esub_p, 'trimmed')
    call do_A08_checksums(esub_p, 'trimmed_scaled')
    call do_A08_checksums(esub_p, 'trimmed_scaled_crops_ext')
    call do_A08_checksums(esub_p, 'trimmed_scaled_nocrops')

end program A09_checksums
