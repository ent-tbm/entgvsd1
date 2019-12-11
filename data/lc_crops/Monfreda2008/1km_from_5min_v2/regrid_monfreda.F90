#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module regrid_monfreda_mod

use chunker_mod
use ent_labels_mod
use gcm_labels_mod
use ent_params_mod
use hntr_mod

implicit none
CONTAINS

subroutine regrid_control(rw, root, dir, ileaf, oleaf, vnames)
    type(ReadWrites_t), intent(INOUT) :: rw
    character*(*), intent(in) :: root
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: ileaf
    character*(*), intent(in) :: oleaf
    character(len=*), dimension(:) :: vnames

    type(Chunker_t) :: chunker_in
    type(Chunker_t) :: chunker_out
    type(ChunkIO_t), dimension(:), allocatable :: io_in, io_out
    type(ChunkIO_t) :: ioall_out

    type(HntrSpec_t) :: spec_in, spec_out
    type(HntrCalc_t) :: hntr_out    ! Preparation to regrid

    integer :: ichunk,jchunk
    integer :: i,nvnames
    type(FileInfo_t) :: info


    call chunker_out%init(im1km,jm1km,imh,jmh,'forplot',4,30,30, nchunk=(/1,15/), outputs_dir=THIS_OUTPUTS_DIR)
    call clear_file_info(info)
    info%use_outputs_dir = .false.

    call chunker_out%nc_create1(ioall_out, weighting(chunker_out%wta1,1d0,0d0), &
        './', trim(oleaf), info)!, vname, vname, '1', create_lr=.true.)

    nvnames = size(vnames,1)
    allocate(io_in(nvnames))
    allocate(io_out(nvnames))

    call clear_file_info(info)
    do i=1,nvnames
        info%vname = vnames(i)
        info%long_name = vnames(i)
        info%data_source = 'Monfreda'
        info%units = '1'
        call chunker_out%nc_reuse_file(ioall_out, io_out(i), info, weighting(chunker_out%wta1,1d0,0d0))
    end do

    call chunker_in%init(IM5m,JM5m,  IMh,JMh, 'forplot', 30, 30, 30, (/1,15/), outputs_dir=THIS_OUTPUTS_DIR)
    do i=1,nvnames
        call chunker_in%nc_open(io_in(i), root, dir, trim(ileaf)//'.nc', vnames(i), 1)
    end do


    spec_out = hntr_spec(chunker_out%chunk_size(1), chunker_out%ngrid(2), 0d0, 180d0*60d0 / chunker_out%ngrid(2))
    spec_in = hntr_spec(chunker_in%chunk_size(1), chunker_in%ngrid(2), 0d0, 180d0*60d0 / chunker_in%ngrid(2))
    hntr_out = hntr_calc(spec_out, spec_in, 0d0)   ! datmis=0

    call chunker_in%nc_check(rw=rw)
    call chunker_out%nc_check(rw=rw)

#ifdef JUST_DEPENDENCIES
    return
#endif


#ifdef ENTGVSD_DEBUG
    do jchunk = 1,1
    do ichunk = 1,1
#else
    do jchunk = 1,chunker_out%nchunk(2)
    do ichunk = 1,chunker_out%nchunk(1)
#endif
        call chunker_in%move_to(ichunk,jchunk)
        call chunker_out%move_to(ichunk,jchunk)

        do i=1,nvnames
            call hntr_out%regrid4( &
                io_out(i)%buf, io_in(i)%buf, &
                chunker_in%wta1, 1d0, 0d0, &   ! weighting
                io_out(i)%startB(2), io_out(i)%chunker%chunk_size(2))
        end do

        call chunker_out%write_chunks
    end do
    end do

    call chunker_out%close_chunks
    call chunker_in%close_chunks
end subroutine regrid_control



end module regrid_monfreda_mod





program A07a_regrid_controls
    use regrid_monfreda_mod
    use ent_params_mod

implicit none
    type(ReadWrites_t) :: rw
    call rw%init(THIS_OUTPUTS_DIR, 'regrid_monfreda', 100,100)

    call regrid_control(rw, './', '', &
        'Monfreda_crops_5min', 'Monfreda_crops_1km', (/ &
        'c3crop_cf            ', &
        'c4crop_cf            ', &
        'herbcrop_cf          ', &
        'shrubcrop_cf         ', &
        'treecrop_cf          ', &
        'c3c4tot_cf           ', &
        'herbshrubtreecrop_ccf' /) )

    call regrid_control(rw, './', '', &
        'Monfreda_crops_5min_norm', 'Monfreda_crops_1km_norm', (/ &
        'c3norm_ccf           ', &
        'c4norm_ccf           ', &
        'c3c4normtot_ccf      ', &
        'herbnorm_ccf         ', &
        'shrubnorm_ccf        ', &
        'treenorm_ccf         ', &
        'herbshrubtreenorm_ccf', &
        'c4multi_hcf          ', &
        'c3multi_cf           ', &
        'c3crop_hcf           ', &
        'shrubtreecrop_cf     ', &
        'c3c4crop_capped      ' /) )

end program A07a_regrid_controls



