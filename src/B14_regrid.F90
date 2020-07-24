! Regrids the PURE files from 1km to 1/2 degree resolution
!    in preparation of further trimming.  
! Creates:
!    purelr with water_ice and 
!    purelr_noh2o without water_ice, as input for trimming
!
! Includes LC, LAIMAX, Height and monthly LAI files
!
! Author: Elizabeth Fischer, Nancy Y. Kiang

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module b14_mod
    use netcdf
    use chunker_mod
    use ent_labels_mod
    use gcm_labels_mod
    use ent_params_mod
    use hntr_mod

implicit none


type IOFname_t
    character*(512) :: root
    character*(100) :: idir,odir
    character*(100) :: ileaf,oleaf
    logical :: lc_weighting
    type(FileInfo_t) :: info   ! For metadata ONLY
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
    fn%info%vname = vname
    fn%info%long_name = long_name
    fn%info%units = units
    fn%info%modification = ''

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
    character*(*), intent(IN) :: istep  ! ent17, pure, trimmed, trimmed_scaled, trimmed_scaled_natveg, trimmed_scaled_natveg_ext, trimmed_scaled_crops_ext (lai only)
    character*(*), intent(IN) :: ostep  ! ent17, pure, trimmed, trimmed_scaled, trimmed_scaled_natveg, trimmed_scaled_natveg_ext, trimmed_scaled_crops_ext (lai only)
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


    fn%root = ochunker%outputs_dir
    fn%idir = iinfo%dir
    fn%ileaf = iinfo%leaf
    fn%odir = oinfo%dir
    fn%oleaf = oinfo%leaf
    fn%lc_weighting = lc_weighting
    fn%info = iinfo   ! Copy all metadata

end function make_fname2


! Regrids a file, using its own FillValue as a mask
subroutine regrid_selfmask(root, idir,ifname,ivname,oinfo, rw)
    character*(*), intent(IN) :: root,idir,ifname,ivname
    type(FileInfo_t), intent(IN) :: oinfo
    type(ReadWrites_t) :: rw
    ! ------------ Local Vars
    type(Chunker_t), target :: chunker, chunkerlr
    type(ChunkIO_t) :: io_valin, io_valout

    type(HntrSpec_t) :: spec_hr, spec_lr
    type(HntrCalc_t) :: hntr_lr    ! Preparation to regrid


    integer :: jchunk, ichunk    ! Index of current chunk
    integer :: jc, ic    ! Index WITHIN current chunk
    integer :: jj, ii            ! Index in full space
    real*4, dimension(:,:), allocatable :: mask
    integer :: k
    integer :: ndoy,idoy
    logical :: need_lc

    print *,'****************** BEGIN Regrid ',trim(idir),trim(ifname),' --> ',trim(oinfo%leaf)

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 300, 320, 20, (/6,5/), outputs_dir=THIS_OUTPUTS_DIR)
    call chunkerlr%init(IMLR,JMLR, IMH*2,JMH*2, 'forplot', 300, 320, 20, (/6,5/), outputs_dir=THIS_OUTPUTS_DIR)

    ! Hntr stuff
    spec_hr = hntr_spec(chunker%chunk_size(1), chunker%ngrid(2), 0d0, 180d0*60d0 / chunker%ngrid(2))
    spec_lr = hntr_spec(chunkerlr%chunk_size(1), chunkerlr%ngrid(2), 0d0, 180d0*60d0 / chunkerlr%ngrid(2))
    hntr_lr = hntr_calc(spec_lr, spec_hr, 0d0)   ! datmis=0


    ! ------------- Input Files
    call chunker%nc_open(io_valin, trim(root), &
        trim(idir), trim(ifname)//'.nc', trim(ivname), 1)

    ! ----------- Output Files
    ! No weighting needed because create_lr==.false.
    allocate(mask(chunker%chunk_size(1),chunker%chunk_size(2)))
    call chunkerlr%nc_create1(io_valout, &
            weighting(mask,1d0,0d0), &
            trim(oinfo%dir), trim(oinfo%leaf), oinfo, &
            create_lr=.false.)

    call chunker%nc_check(rw=rw)
    call chunkerlr%nc_check(rw=rw)

#ifdef JUST_DEPENDENCIES
    return
#endif


#ifdef ENTGVSD_DEBUG
    do jchunk = dbj0_lc,dbj1_lc
    do ichunk = dbi0_lc,dbi1_lc
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif

        call chunker%move_to(ichunk,jchunk)
        call chunkerlr%move_to(ichunk,jchunk)

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)
            mask(ic,jc) = 1
            if (io_valin%buf(ic,jc) == FillValue.or.io_valin%buf(ic,jc) /= io_valin%buf(ic,jc)) then
                mask(ic,jc) = 0.
                io_valin%buf(ic,jc) = FillValue
            end if
        end do
        end do


        call hntr_lr%regrid4( &
            io_valout%buf, io_valin%buf, &
            mask, 1d0, 0d0, &   ! weighting
            io_valout%startB(2), io_valout%chunker%chunk_size(2))

        ! call chunker%write_chunks
        call chunkerlr%write_chunks

    end do
    end do

    call chunker%close_chunks
    call chunkerlr%close_chunks

end subroutine regrid_selfmask




subroutine regrid_lais(esub, fname, rw)
    type(EntSet_t), intent(IN) :: esub
    type(IOFname_t),intent(IN), target :: fname(:)
    type(ReadWrites_t) :: rw
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
    logical :: need_lc

    ndoy = size(fname,1)
    do idoy=1,ndoy
        print *,'****************** BEGIN Regrid ',trim(fname(idoy)%ileaf),' --> ',trim(fname(idoy)%oleaf)
    end do

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 300, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)
    call chunkerlr%init(IMLR,JMLR, IMH*2,JMH*2, 'forplot', 300, 320, 20, outputs_dir=THIS_OUTPUTS_DIR)

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
            LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1')
    end if

    ! LAI
    do idoy=1,ndoy
        fn => fname(idoy)

        call chunker%nc_open(ioall_lai(idoy), trim(fn%root), &
            trim(fn%idir), trim(fn%ileaf)//'.nc', trim(fn%info%vname), 0)
        do k = 1,esub%ncover
            call chunker%nc_reuse_var(ioall_lai(idoy), io_lai(k,idoy), (/1,1,k/))
        enddo

        ! ----------- Output Files
        ! No weighting needed because create_lr==.false.
        call chunkerlr%nc_create1(ioall_laiout(idoy), &
            weighting(chunkerlr%wta1,1d0,0d0), &
            trim(fn%odir), trim(fn%oleaf), fn%info, &
            esub%layer_names(), esub%long_layer_names(), &
            create_lr=.false.)
        do k=1,esub%ncover
            call chunkerlr%nc_reuse_var(ioall_laiout(idoy), io_laiout(k,idoy), &
                (/1,1,k/), weighting(chunkerlr%wta1,1d0,0d0))
        end do   ! k

    end do    ! idoy

    call chunker%nc_check(rw=rw)
    call chunkerlr%nc_check(rw=rw)

#ifdef JUST_DEPENDENCIES
    return
#endif

#ifdef ENTGVSD_DEBUG
    do jchunk = dbj0,dbj1
    do ichunk = dbi0,dbi1
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif

       call chunker%move_to(ichunk,jchunk)
       call chunkerlr%move_to(ichunk,jchunk)

        do idoy=1,ndoy
        do k=1,esub%ncover
            if (fname(idoy)%lc_weighting) then
                call hntr_lr%regrid4( &
                    io_laiout(k,idoy)%buf, io_lai(k,idoy)%buf, &
                    io_lc_pure(k)%buf, 1d0, 0d0, &   ! weighting
                    io_laiout(k,idoy)%startB(2), io_laiout(k,idoy)%chunker%chunk_size(2))
            else
                call hntr_lr%regrid4( &
                    io_laiout(k,idoy)%buf, io_lai(k,idoy)%buf, &
                    chunker%wta1, 1d0, 0d0, &   ! weighting
                    io_laiout(k,idoy)%startB(2), io_laiout(k,idoy)%chunker%chunk_size(2))
            end if

        end do    ! k=1,esub%ncover
        end do    ! idoy

        call chunker%write_chunks
        call chunkerlr%write_chunks

    end do
    end do

    call chunker%close_chunks
    call chunkerlr%close_chunks
end subroutine regrid_lais


subroutine do_regrid_all_lais(rw)
    use ent_labels_mod
    use gcm_labels_mod
    use ent_params_mod
    type(ReadWrites_t) :: rw
    ! --------------- Locals
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
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 0,0,0, outputs_dir=THIS_OUTPUTS_DIR)
    call chunkerlr%init(IMLR,JMLR, IMH*2,JMH*2, 'forplot', 0,0,0, outputs_dir=THIS_OUTPUTS_DIR)

    nf = 0

    ! ----------- Annual
    nf = nf + 1
    fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
        .false., LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', 'purelr', '1.1')


    nf = nf + 1
    fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
        .true., LAI_SOURCE, 'M', 'laimax', LAI_YEAR, 'pure', 'purelr', '1.1')

    nf = nf + 1
    fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
        .true., LAI_SOURCE, 'M', 'hgt', LAI_YEAR, 'pure', 'purelr', '1.1')

#if 0
! We only need the doy files at the ent17 and pure steps.  We do not
! need them at later trimming steps. The doy files are only for 1km
! hi-res evaluation against MODIS albedo.

    ! ---------- DOY
    do idoy=1,NDOY
        nf = nf + 1
        fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
            .true., LAI_SOURCE, 'M', 'lai', LAI_YEAR, 'pure', 'purelr', '1.1', &
            'doy', idoy)
    end do
#endif

    ! ---------- MONTH
    do imonth=1,NMONTH
        nf = nf + 1
        fname(nf) = make_fname2(chunker, chunkerlr, esub_p, &
            .true., LAI_SOURCE, 'M', 'lai', LAI_YEAR, 'pure', 'purelr', '1.1', &
            'month', imonth)
    end do


    if (nf > MAXFILES) then
        write(ERROR_UNIT,*) 'MAXFILES too small',nf
        STOP
    end if




    ! ---------------------------------
    ! NOT NEEDED: Use file already created by B07_soilalbedo.F90.
    ! Create single-layer EntSet_t
    ! Bare soil Bright Ratio
    !call chunkerlr%file_info(oinfo, ent1, &
    !    LAI_SOURCE, 'M', 'bs_brightratio', LAI_YEAR, 'purelr', '1.1')
    !
    !oinfo%leaf = 'bs_brightratio'
    !call regrid_selfmask( &
    !    chunkerlr%outputs_dir, 'soilalbedo/', 'soilalbedo_1km_bs_brightratio_fill', 'bs_brightratio', &
    !    oinfo, rw)
    ! ---------------------------------

    esub_p => esub
    do i0=1,nf,10
        i1 = min(nf,i0 + 10 - 1)

        call regrid_lais(esub_p, fname(i0:i1), rw)
    end do

end subroutine do_regrid_all_lais



end module b14_mod

! =========================================================
program regrid
    use b14_mod
implicit none
    type(ReadWrites_t) :: rw

    call rw%init(THIS_OUTPUTS_DIR, "B14_regrid", 40,40)
    call do_regrid_all_lais(rw)
    call rw%write_mk

end program regrid
