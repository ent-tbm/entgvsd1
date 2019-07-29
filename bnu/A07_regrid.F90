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
    character*(100) :: ileaf,oleaf
    logical :: lc_weighting
    character*(40) :: vname
    character*(100) :: long_name
    character*(30) :: units
    character*(60) :: title
end type IOFname_t

CONTAINS

function make_fname( &
    root, idir, ileaf, odir, oleaf, lc_weighting, &
    vname, long_name, units) result(fn)
    character*(*), intent(IN) :: root, idir, ileaf, odir, oleaf
    logical, intent(IN) :: lc_weighting
    character*(*), intent(IN) :: vname, long_name, units
    type(IOFname_t) :: fn

    fn%root = root
    fn%idir = idir
    fn%ileaf = ileaf
    fn%odir = odir
    fn%oleaf = oleaf
    fn%lc_weighting = lc_weighting
    fn%vname = vname
    fn%long_name = long_name
    fn%units = units

end function make_fname


!type IOFname_t
!    type(FileInfo_t) :: iinfo, oinfo
!    logical :: lc_weighting
!end type IOFname_t
!
!CONTAINS


function make_fname2(ichunker, ochunker, ents, lc_weighting, &
    laisource, cropsource, var, year, istep, ostep, ver, &
    doytype, idoy, varsuffix) &
result(fn)
    class(Chunker_t), intent(IN) :: ichunker, ochunker
    type(EntSet_t), intent(IN) :: ents
    logical, intent(IN) :: lc_weighting

    character*(*), intent(IN) :: laisource   ! M (MODIS, Nancy’s old version), BNU (Carl and Elizabeth)

    character*(*), intent(IN) :: cropsource   ! M (Monfreda et al. 2008)
    character*(*), intent(IN) :: var    ! lc, lai, laimax, height
    integer, intent(IN) :: year
    character*(*), intent(IN) :: istep  ! raw, pure, trimmed, trimmed_scaled, trimmed_scaled_nocrops, trimmed_scaled_nocrops_ext, trimmed_scaled_crops_ext (lai only)
    character*(*), intent(IN) :: ostep  ! raw, pure, trimmed, trimmed_scaled, trimmed_scaled_nocrops, trimmed_scaled_nocrops_ext, trimmed_scaled_crops_ext (lai only)
    character*(*) , intent(IN) :: ver  ! 1.1, 1.2, etc
    character*(*), intent(IN), OPTIONAL :: doytype ! ann,doy,month
    integer, intent(IN), OPTIONAL :: idoy
    character*(*), intent(IN), OPTIONAL :: varsuffix
    type(IOFname_t) :: fn       ! RESULT

    ! ------------- Locals
    type(FileInfo_t) :: iinfo, oinfo

    call ichunker%file_info(iinfo, ents, &
        laisource, cropsource, var, year, istep, ver, doytype, idoy, varsuffix)
    call ochunker%file_info(oinfo, ents, &
        laisource, cropsource, var, year, ostep, ver, doytype, idoy, varsuffix)

    fn = make_fname(LC_LAI_ENT_DIR, &
        iinfo%dir, iinfo%leaf, &
        oinfo%dir, oinfo%leaf, &
        lc_weighting, &
        iinfo%vname, iinfo%long_name, iinfo%units)
end function make_fname2


subroutine regrid_lais(esub, fname)
    type(EntSet_t), intent(IN) :: esub
    type(IOFname_t),intent(IN), target :: fname(:)
    ! ------------ Local Vars
    type(Chunker_t), target :: chunker, chunkerlr
    type(ChunkIO_t), target :: io_lc_pure(esub%ncover)
    type(ChunkIO_t) :: ioall_lai(size(fname,1)), io_lai(esub%ncover,size(fname,1))
    type(ChunkIO_t) :: ioall_laiout(size(fname,1)), io_laiout(esub%ncover,size(fname,1))

    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid


    integer :: jchunk, ichunk    ! Index of current chunk
    integer :: jc, ic    ! Index WITHIN current chunk
    integer :: jj, ii            ! Index in full space
    integer :: k
    integer :: ndoy,idoy
    type(IOFname_t), pointer :: fn
    type(Weighting_t) :: wgt
    logical :: need_lc

    ndoy = size(fname,1)
    do idoy=1,ndoy
        print *,'****************** BEGIN Regrid ',trim(fname(idoy)%ileaf),' --> ',trim(fname(idoy)%oleaf)
    end do

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 300, 320, 20)
    call chunkerlr%init(IMLR,JMLR, IMH*2,JMH*2, 'qxq', 300, 320, 20)

    ! Hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1), chunker%ngrid(2), 0d0, 180d0*60d0 / chunker%ngrid(2))
    spec_lr = hntr_spec(chunkerlr%chunk_size(1), chunkerlr%ngrid(2), 0d0, 180d0*60d0 / chunkerlr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, 0d0)   ! datmis=0

    ! Figure out if we need LC weighting
    need_lc = .false.
    do idoy=1,ndoy
        need_lc = need_lc.or.fname(idoy)%lc_weighting
    end do


    ! ------------- Input Files
    ! LC written by A04; in the esub indexing scheme
    if (need_lc) then
        call chunker%nc_open_set(esub, io_lc_pure, &
            'BNU', 'M', 'lc', 2004, 'pure', '1.1')
    end if

    ! LAI
    do idoy=1,ndoy
        fn => fname(idoy)

        call chunker%nc_open(ioall_lai(idoy), trim(fn%root), &
            trim(fn%idir), trim(fn%ileaf)//'.nc', trim(fn%vname), 0)
        do k = 1,esub%ncover
            call chunker%nc_reuse_var(ioall_lai(idoy), io_lai(k,idoy), (/1,1,k/))
        enddo

        ! ----------- Output Files
        ! No weighting needed because create_lr==.false.
        call chunkerlr%nc_create(ioall_laiout(idoy), &
            weighting(chunkerlr%wta1,1d0,0d0), &
            trim(fn%odir), trim(fn%oleaf), trim(fn%vname), &
            trim(fn%long_name), trim(fn%units), &
            esub%layer_names(), &
            create_lr=.false.)
        do k=1,esub%ncover
            call chunkerlr%nc_reuse_var(ioall_laiout(idoy), io_laiout(k,idoy), &
                (/1,1,k/), weighting(chunkerlr%wta1,1d0,0d0))
        end do   ! k

    end do    ! idoy

    call chunker%nc_check('A07_regrid_hr')
    call chunkerlr%nc_check('A07_regrid_lr')

!    call regrid_handles(chunker, chunkerlr, esub
!
!end subroutine regrid_lais
!
!subroutine regrid_handles()

#ifdef ENTGVSD_DEBUG
!    do jchunk = chunker%nchunk(2)*3/4,chunker%nchunk(2)*3/4+1
!    do ichunk = chunker%nchunk(1)*3/4,chunker%nchunk(1)*3/4+1
    ! Choose a smalle area to debug
    do jchunk = 11,12
    do ichunk = 5,7
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif

       call chunker%move_to(ichunk,jchunk)
       call chunkerlr%move_to(ichunk,jchunk)

        do idoy=1,ndoy
        do k=1,esub%ncover
            if (fn%lc_weighting) then
                wgt%buf => io_lc_pure(k)%buf
                wgt%MM = 1d0
                wgt%BB = 030
            else
                wgt%buf => chunker%wta1
                wgt%MM = 1d0
                wgt%BB = 030
            end if

            call hntr_lr%regrid4( &
                io_laiout(k,idoy)%buf, io_lai(k,idoy)%buf, &
                wgt%buf, wgt%MM, wgt%BB, &    ! weighting
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


subroutine do_regrid_all_lais

    use ent_labels_mod
    use gcm_labels_mod
    use chunkparams_mod

    type(EntSet_t) :: ent1
    type(GcmEntSet_t), target :: esub
    class(EntSet_t), pointer :: esub_p
    type(FileInfo_t) :: iinfo,oinfo

    integer, parameter :: MAXFILES = 20 !(1 + (1 + 2 + 12))   ! LC + LAI(annual + doy + monthly)
    type(IOFname_t) :: fname(MAXFILES),fnames1(1)
    integer :: nf,i0,i1,idoy,imonth
    type(Chunker_t) :: chunker, chunkerlr

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)
    esub_p => esub

    ! Chunkers just for make_fname2()
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 0,0,0)
    call chunkerlr%init(IMLR,JMLR, IMH*2,JMH*2, 'qxq', 0,0,0)

    nf = 0

    ! ----------- Annual
    nf = nf + 1
    fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
        .false., 'BNU', 'M', 'lc', 2004, 'pure', 'purelr', '1.1')


    nf = nf + 1
    fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
        .true., 'BNU', 'M', 'laimax', 2004, 'pure', 'purelr', '1.1')

    nf = nf + 1
    fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
        .true., 'BNU', 'M', 'hgt', 2004, 'pure', 'purelr', '1.1')

    ! ---------- DOY
    do idoy=1,NDOY
        nf = nf + 1
        fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
            .true., 'BNU', 'M', 'lai', 2004, 'pure', 'purelr', '1.1', &
            'doy', idoy)
    end do

    ! ---------- MONTH
    do imonth=1,NMONTH
        nf = nf + 1
        fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
            .true., 'BNU', 'M', 'lai', 2004, 'pure', 'purelr', '1.1', &
            'month', imonth)
    end do


    if (nf > MAXFILES) then
        write(ERROR_UNIT,*) 'MAXFILES too small',nf
        STOP
    end if




    ! ---------------------------------
    ! Create single-layer EntSet_t
    call ent1%allocate(1,1)
    call ent1%add_covertype('n', 'bs_brightratio', 'Bare Soil Brightness Ratio')

    ! Bare soil Bright Ratio
    call chunkerlr%file_info(oinfo, ent1, &
        'BNU', 'M', 'bs_brightratio', 2004, 'purelr', '1.1')
    fnames1(1) = &
        make_fname(LAI3G_INPUT, &
            'lc_lai_ent/', 'bs_brightratio', &   ! 
            oinfo%dir, 'bs_brightratio', & !oinfo%leaf, &
            .false., &
            oinfo%vname, oinfo%long_name, oinfo%units)

    call regrid_lais(ent1, fnames1)
    ! ---------------------------------

    esub_p => esub
    do i0=1,nf,10
        i1 = min(nf,i0 + 10 - 1)

        call regrid_lais(esub_p, fname(i0:i1))
    end do

end subroutine do_regrid_all_lais



end module a07_mod

! =========================================================
program regrid
    use a07_mod
implicit none
    call do_regrid_all_lais
end program regrid
