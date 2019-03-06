module a06_mod

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use gcm_labels_mod
    use geom_mod
    use trim_tiny_mod

implicit none


    ! Combine C3 and C4 crops into one PFT, for ModelE
    logical, parameter :: combine_crops_c3_c4 = .true.
    ! Split the bare soil into dark and light, to get the right albedo
    logical, parameter :: split_bare_soil = .true.

CONTAINS

subroutine do_reindex(esub,m0,m1)
    type(GcmEntSet_t), intent(IN) :: esub
    integer :: m0,m1    ! First and last month to work on

    type(Chunker_t) :: chunker

    ! ------ Input Files
    type(ChunkIO_t) :: io_lcin(NENT20), io_bs
    type(ChunkIO_t) :: ioall_laiin(m1-m0+1), io_laiin(NENT20,m1-m0+1)
    ! ------ Output files
    type(ChunkIO_t) :: ioall_laiout(m1-m0+1), io_laiout(esub%ncover,m1-m0+1)


    integer :: k,ksub
    integer :: im,imonth


    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'qxq', 150, 150)

    !------------------------------------------------------------------------
    ! OPEN INPUT FILES

    ! lcin
    do k=1,NENT20
        ! TODO: LAI3g????
        ! PathFilepre= '../../LAI3g/lc_lai_ent/EntMM_lc_laimax_1kmx1km/'
        call chunker%nc_open(io_lcin(k), LC_LAI_ENT_DIR, &
            'EntMM_lc_laimax_1kmx1km/', &
            trim(itoa2(k))//'_'//trim(ent20%abbrev(k))//'_lc.nc', &
            trim(ent20%abbrev(k)), 1)
    enddo


    ! laiin
    do im = m0,m1
        imonth = im - m0 + 1
        call chunker%nc_open(ioall_laiin(imonth), LC_LAI_ENT_DIR, &
            'nc/', 'EntMM_lc_lai_'//trim(MONTH(imonth))//'_1kmx1km.nc', 'EntPFT', 0)

        do k=1,NENT20
            call chunker%nc_reuse_var(ioall_laiin(imonth), io_laiin(k,imonth), &
                (/1,1,k/), 'r', weighting(chunker%wta1,1d0,0d0))
        end do

    end do

    ! bs ratio
    call chunker%nc_open_gz(io_bs, LAI3G_DIR, LAI3G_INPUT, &
        'lc_lai_ent/', 'bs_brightratio.nc', 'bs_brightratio', 1)
    ! For now at least, use BS BrightRatio from LAI3G dataset
    ! call chunker%nc_open(io_bs, LC_LAI_ENT_DIR, &
    !     '', 'bs_brightratio.nc', 'bs_brightratio', 1)

    !------------------------------------------------------------------------
    ! CREATE OUTPUT NETCDF FILES

    do im = m0,m1
        imonth = im - m0 + 1
        call chunker%nc_create(ioall_laiout(imonth), &
            weighting(chunker%wta1,1d0,0d0), &
            '16/nc/', 'V1km_EntGVSDv1.1_BNU16_lai_'//MONTH(imonth)//'_pure')

        do k=1,esub%ncover
            call chunker%nc_reuse_file(ioall_laiout(imonth), io_laiout(k,imonth), &
                'lai_'//trim(esub%abbrev(k)), &
                trim(esub%title(k)), 'm2 m-2', &
                trim(esub%title(k)), &
                weighting(chunker%wta1,1d0,0d0))   ! TODO: What convert LC to 18

        end do  ! k
    end do   ! imonth

    call chunker%nc_check('A06_trim_lc_laimax_laimonth')
    print *,'Done opening files: nreads',chunker%nreads,'nwrite',chunker%nwrites
#ifdef JUST_DEPENDENCIES
    stop 0
#endif

    call trim_tiny(esub, chunker, m1-m0+1, &

#ifdef ENTGVSD_DEBUG
        nchunk(2)*3/4,nchunk(2)*3/4+1, &
        nchunk(1)*3/4,nchunk(1)*3/4+1, &
#else
        1,nchunk(2), &
        1,nchunk(1), &
#endif
        combine_crops_c3_c4, split_bare_soil, &
        io_lcin, io_laiin, io_bs, &
        io_laiout)

    call chunker%close_chunks

end subroutine do_reindex
end module a06_mod
! ====================================================================

program convert
    use a06_mod
    use ent_labels_mod
    use gcm_labels_mod
implicit none

    ! -------------------------------------------------------
    type(GcmEntSet_t) :: esub

    call init_ent_labels
    esub = make_ent_gcm_subset(combine_crops_c3_c4, split_bare_soil)

    call do_reindex(esub,1,6)
    call do_reindex(esub,7,12)


end program convert
