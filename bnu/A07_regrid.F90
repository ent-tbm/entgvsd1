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

type IOFname_t
    character*(512) :: root
    character*(100) :: idir,odir
    character*(100) :: leaf

end type IOFname_t

CONTAINS

function make_fname(root, idir, odir, leaf) result(fn)
    character*(*), intent(IN) :: root, idir, odir, leaf
    type(IOFname_t) :: fn

    fn%root = root
    fn%idir = idir
    fn%odir = odir
    fn%leaf = leaf

end function make_fname

subroutine regrid_lais(esub, fname)
    type(GcmEntSet_t), intent(IN) :: esub
    type(IOFname_t),intent(IN) :: fname(:)
    ! ------------ Local Vars
    type(Chunker_t) :: chunker, chunkerlr
    type(ChunkIO_t) :: ioall_lc2, io_lc2(esub%ncover)
    type(ChunkIO_t) :: ioall_lai(size(fname,1)), io_lai(esub%ncover,size(fname,1))
    type(ChunkIO_t) :: ioall_laiout(size(fname,1)), io_laiout(esub%ncover,size(fname,1))

    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid

    integer :: jchunk, ichunk    ! Index of current chunk
    integer :: jc, ic    ! Index WITHIN current chunk
    integer :: jj, ii            ! Index in full space
    integer :: k
    integer :: ndoy,idoy

    integer, parameter :: IMLR=IMH
    integer, parameter :: JMLR=JMH

    ndoy = size(fname,1)
    do idoy=1,ndoy
        print *,'****************** BEGIN Regrid ',trim(fname(idoy)%leaf)
    end do

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 300, 320)
    call chunkerlr%init(IMLR,JMLR, IMH*2,JMH*2, 'qxq', 300, 320)

    ! Hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1), chunker%ngrid(2), 0d0, 180d0*60d0 / chunker%ngrid(2))
    spec_lr = hntr_spec(chunkerlr%chunk_size(1), chunkerlr%ngrid(2), 0d0, 180d0*60d0 / chunkerlr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, 0d0)   ! datmis=0


    ! ------------- Input Files
    ! LC written by A04; in the esub indexing scheme
    call chunker%nc_open(ioall_lc2, LC_LAI_ENT_DIR, &
        'pure2/annual/', 'entmm29_ann_lc.nc', 'lc', 0)
    do k = 1,esub%ncover
        call chunker%nc_reuse_var(ioall_lc2, io_lc2(k), (/1,1,k/))
    enddo

    ! LAI
    do idoy=1,ndoy
        call chunker%nc_open(ioall_lai(idoy), trim(fname(idoy)%root), &
            trim(fname(idoy)%idir), trim(fname(idoy)%leaf)//'.nc', 'lai', 0)
        do k = 1,esub%ncover
            call chunker%nc_reuse_var(ioall_lai(idoy), io_lai(k,idoy), (/1,1,k/))
        enddo

        ! ----------- Output Files
        call chunkerlr%nc_create(ioall_laiout(idoy), &
            weighting(chunker%wta1,1d0,0d0), &
            trim(fname(idoy)%odir), trim(fname(idoy)%leaf), 'lai', &
            'Ent LAI (annual, DOY or monthly)', 'm^2 m-2', 'Leaf Area Index', &
            esub%mvs, esub%layer_names(), &
            create_lr=.false.)
        do k=1,esub%ncover
            call chunkerlr%nc_reuse_var(ioall_laiout(idoy), io_laiout(k,idoy), &
                (/1,1,k/), weighting(io_lc2(k)%buf, 1d0,0d0))
        end do   ! k

    end do    ! idoy


#ifdef ENTGVSD_DEBUG
    do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
    do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
    do jchunk = 1,nchunk(2)
    do ichunk = 1,nchunk(1)
#endif

       call chunker%move_to(ichunk,jchunk)
       call chunkerlr%move_to(ichunk,jchunk)

        do idoy=1,ndoy
        do k=1,esub%ncover
            call hntr_lr%regrid4( &
                io_laiout(k,idoy)%buf, io_lai(k,idoy)%buf, &
                io_lc2(k)%buf, 1d0, 0d0, &        ! weighting
                io_laiout(k,idoy)%startB(2), io_laiout(k,idoy)%chunker%chunk_size(2))
        end do    ! k=1,esub%ncover
        end do    ! idoy

        call chunker%write_chunks
        call chunkerlr%write_chunks

    end do
    end do

    call chunker%close_chunks
    call chunkerlr%close_chunks
end subroutine regrid_lais


subroutine do_regrid_all

    use ent_labels_mod
    use gcm_labels_mod
    use chunkparams_mod

    type(GcmEntSet_t) :: esub
    integer, parameter :: MAXFILES = (1 + (1 + 2 + 12))   ! LC + LAI(annual + doy + monthly)
    type(IOFname_t) :: fname(MAXFILES)
    integer :: nf,i0,i1,idoy,imonth

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)


    nf = 0

    ! ----------- Annual
!    nf = nf + 1
!    fname(nf) = make_fname(LC_LAI_ENT_DIR, &
!        'pure/annual/', 'purelr/annual/', 'entmm29_ann_lc')

    nf = nf + 1
    fname(nf) = make_fname(LC_LAI_ENT_DIR, &
        'pure2/annual/', 'purelr/annual/', 'entmm29_ann_laimax')

    ! ---------- DOY
    do idoy=1,NDOY
        nf = nf + 1
        fname(nf) = make_fname(LC_LAI_ENT_DIR, &
            'pure2/doy/', 'purelr/doy/', 'entmm29_'//DOY(idoy)//'_lai')
    end do

    ! ---------- MONTH
    do imonth=1,NMONTH
        nf = nf + 1
        fname(nf) = make_fname(LC_LAI_ENT_DIR, &
            'pure2/monthly/', 'purelr/monthly/', 'entmm29_'//MONTH(imonth)//'_lai')
    end do


    do i0=1,nf,10
        i1 = min(nf,i0 + 10 - 1)

!        print *,i0,i1
        call regrid_lais(esub, fname(i0:i1))
    end do

end subroutine do_regrid_all



end module a07_mod

! =========================================================
program regrid
    use a07_mod
implicit none
    call do_regrid_all
end program regrid
