! Reclass monthly LAI from 20-cover classes to 18-cover classes
!
! Takes 20-cover classes (Ent 17 PFTs + 3 non-veg) and converts to
! 18-cover classes (Ent 16 PFTs and bright and dark bare soil fractions).
! Merges C3 and C4 crops into one crop cover type for Ent 16 PFTs, excludes water
! and permanent ice, and converts barse/sparse cover into equivalent veg type and
! bare soil bright and dark fractions, preserving total LAI of grid cell.
!
! Author: Nancy Kiang, Carlo Monte, Elizabeth Fischer
!
! See slib/cropmerge_laisparse_splitbare.f90

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module b13_mod

    use netcdf
    use chunker_mod
    use ent_labels_mod
    use gcm_labels_mod
    use ent_params_mod
    use cropmerge_laisparse_splitbare_mod

implicit none

CONTAINS

subroutine do_reindex(esub,m0,m1,rw)
    type(GcmEntSet_t), intent(IN), target :: esub
    integer :: m0,m1    ! First and last month to work on
    type(ReadWrites_t) :: rw

    class(EntSet_t), pointer :: esub_p
    type(Chunker_t) :: chunker

    ! ------ Input Files
    type(ChunkIO_t) :: io_lc_ent17(NENT20)
    type(ChunkIO_t) :: io_lc_pure(esub%ncover)
    real*4, ALLOCATABLE :: sum_lc(:,:)
    type(ChunkIO_t) :: io_laiin(NENT20,m1-m0+1)
    type(ChunkIO_t) :: io_bs
    type(ChunkIO_t) :: io_TCinave
    ! ------ Output files
    type(ChunkIO_t) :: io_laiout(esub%ncover,m1-m0+1)
    type(ChunkIO_t) :: io_lclai_checksum(m1-m0+1)


    type(FileInfo_t) :: info
    integer :: k,ksub
    integer :: im,imonth


    esub_p => esub
    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 300, 300, 20, outputs_dir=THIS_OUTPUTS_DIR)
    allocate(sum_lc(chunker%chunk_size(1), chunker%chunk_size(2)))

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! --- ENTPFTLC: Open outputs written by A00
    call chunker%nc_open_set(ent20, io_lc_ent17, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')

    ! LC written by A04; in the esub indexing scheme
    call chunker%nc_open_set(esub_p, io_lc_pure, &
        LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'pure', '1.1')

    ! laiin
    do im = m0,m1
        imonth = im - m0 + 1

        call chunker%nc_open_set(ent20, io_laiin(:,imonth), &
            LAI_SOURCE, 'M', 'lai', LAI_YEAR, 'ent17', '1.1', &
            doytype='month', idoy=im)
    end do

    ! bs ratio
    call chunker%nc_open(io_bs, chunker%outputs_dir, 'soilalbedo/', &
        'soilalbedo_1km_bs_brightratio_fill.nc', 'bs_brightratio', 1)

    ! Climate statistics (we want TCinave = temperature [C])
    call chunker%nc_open_input(io_TCinave, &
        INPUTS_URL, INPUTS_DIR, &
        'climstats/CRU-TS3.22_GPCC-V6/', &
        'TCinave.nc', 'TCinave', 1)

    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    do im = m0,m1
        imonth = im - m0 + 1

        call chunker%nc_create_set( &
            esub_p, io_laiout(:,imonth), lc_weights(io_lc_pure, 1d0, 0d0), &
            LAI_SOURCE, 'M', 'lai', LAI_YEAR, 'pure', '1.1', &
            doytype='month', idoy=im)

        call chunker%file_info(info, esub_p, &
            LAI_SOURCE, 'M', 'lclai', LAI_YEAR, 'pure', '1.1', &
            varsuffix = '_checksum', &
            doytype='month', idoy=im)
        call chunker%nc_create(io_lclai_checksum(imonth), &
            weighting(sum_lc, 1d0, 0d0), &
            info%dir, info%leaf, info%vname, &
            'SUM(LC*LAI)', info%units)
    end do   ! imonth

    call chunker%nc_check(rw=rw)
    print *,'Done opening files: nreads',chunker%nreads,'nwrite',chunker%nwrites
#ifdef JUST_DEPENDENCIES
    return
#endif

    call cropmerge_laisparse_splitbare(esub, chunker, m1-m0+1, &

#ifdef ENTGVSD_DEBUG
        dbj0,dbj1, &
        dbi0,dbi1, &
#else
        1,chunker%nchunk(2), &
        1,chunker%nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lc_ent17, io_laiin, io_bs, io_TCinave, &
        io_laiout, &
        sum_lc = sum_lc, &
        io_lclai_checksum=io_lclai_checksum)

    call chunker%close_chunks

end subroutine do_reindex
end module b13_mod
! ====================================================================

program convert
    use b13_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t), target :: esub
    integer, parameter :: monthchunk = 6
    integer :: m0,m1
    type(ReadWrites_t) :: rw

    call rw%init(THIS_OUTPUTS_DIR, 'B13_reclass_monthly', 300,300)

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    do m0=1,nmonth,monthchunk
        m1 = min(m0+monthchunk-1, nmonth)

        call do_reindex(esub,m0,m1,rw)
    end do

    call rw%write_mk


end program convert
