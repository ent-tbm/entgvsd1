! Regrids original data files from 1km to 1/2 degree resolution
!    for data check comparison to 1/2 degree computed files.
!
! Includes LAIMAX, Height and monthly LAI files
!
! Author: Elizabeth Fischer


#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module B15_mod

use chunker_mod
use ent_labels_mod
use gcm_labels_mod
use ent_params_mod
use hntr_mod

implicit none
CONTAINS

subroutine regrid_control(rw, root, dir, leaf, vname)
    type(ReadWrites_t), intent(INOUT) :: rw
    character*(*), intent(in) :: root
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    character*(*), intent(in) :: vname

    type(Chunker_t) :: chunker_hr
    type(Chunker_t) :: chunker_lr
    type(ChunkIO_t) :: io_hr, io_lr

    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid

    integer :: ichunk,jchunk

    call chunker_hr%init(im1km,jm1km,imh,jmh,'forplot',4,1,1, nchunk=(/1,5/), outputs_dir=THIS_OUTPUTS_DIR)
    call chunker_hr%nc_open(io_hr, root, dir, trim(leaf)//'.nc', vname, 1)

    call chunker_lr%init(IMLR,JMLR,  IMLR,JMLR, 'forplot', 1, 4, 1, (/1,5/), outputs_dir=THIS_OUTPUTS_DIR)
    call chunker_lr%nc_create(io_lr, weighting(chunker_lr%wta1,1d0,0d0), &
        'tmp/regrids/', trim(leaf)//'_hxh', vname, vname, '1', create_lr=.false.)


    spec_lr = hntr_spec(chunker_lr%chunk_size(1), chunker_lr%ngrid(2), 0d0, 180d0*60d0 / chunker_lr%ngrid(2))
    spec_hr = hntr_spec(chunker_hr%chunk_size(1), chunker_hr%ngrid(2), 0d0, 180d0*60d0 / chunker_hr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, 0d0)   ! datmis=0

    call chunker_hr%nc_check(rw=rw)
    call chunker_lr%nc_check(rw=rw)

#ifdef JUST_DEPENDENCIES
    return
#endif


#ifdef ENTGVSD_DEBUG
    do jchunk = 1,1
    do ichunk = 1,1
#else
    do jchunk = 1,chunker_lr%nchunk(2)
    do ichunk = 1,chunker_lr%nchunk(1)
#endif
        call chunker_hr%move_to(ichunk,jchunk)
        call chunker_lr%move_to(ichunk,jchunk)

        call hntr_lr%regrid4( &
            io_lr%buf, io_hr%buf, &
            chunker_hr%wta1, 1d0, 0d0, &   ! weighting
            io_lr%startB(2), io_lr%chunker%chunk_size(2))

        call chunker_lr%write_chunks
    end do
    end do

    call chunker_lr%close_chunks
    call chunker_hr%close_chunks
end subroutine regrid_control

subroutine regrid_controls(rw, root)
    type(ReadWrites_t) :: rw
    character*(*) :: root

    integer :: imonth

    if (LAI_SOURCE == 'L') then
        call regrid_control(rw, INPUTS_DIR, 'LAI/', 'LAI3gMax_1kmx1km', 'laimax')
    else if (LAI_SOURCE == 'B') then
        call regrid_control(rw, root, 'tmp/bnu/', 'bnu_laimax', 'laimax')
    end if

#if 1
    call regrid_control(rw, &
        INPUTS_DIR, 'height/', 'simard_forest_heights', 'heights')


    do imonth=1,NMONTH
        call regrid_control(rw, &
            LAI_INPUTS_DIR, 'lai/BNU/monthly/'//sLAI_YEAR//'/', 'global_30s_'//sLAI_YEAR//'_'//MONTH(imonth), 'lai')
    end do
#endif

end subroutine regrid_controls


end module B15_mod


program A07a_regrid_controls
    use B15_mod
    use ent_params_mod

implicit none
    type(ReadWrites_t) :: rw
    call rw%init(THIS_OUTPUTS_DIR, 'B15_regrid_controls', 100,100)
#if JUST_DEPENDENCIES
    call regrid_controls(rw, MKFILES_DIR)
#else
    call regrid_controls(rw, DEFAULT_OUTPUTS_DIR)
#endif
    call rw%write_mk

end program A07a_regrid_controls
