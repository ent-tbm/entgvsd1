! Computes min,max,mean,std of Carrer soil albedo
! Elizabeth Fischer <elizabeth.fischer@columbia.edu>
! August 15, 2019

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

module B05_mod

    use netcdf
    use chunker_mod
    use ent_labels_mod
    use ent_params_mod
    use assign_laimax_mod
    use carrer_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


CONTAINS

subroutine do_carrer_mean(rw, iband, ndates)
    type(ReadWrites_t) :: rw
    integer, intent(IN) :: iband
    integer, intent(IN) :: ndates
    ! ------------ Local Vars

    type(Chunker_t) :: chunker
    type(FileInfo_t) :: info
    ! Input files
    type(ChunkIO_t), target :: ioall_albin
    type(ChunkIO_t), target :: io_albin(ndates)
    ! Output files
    type(ChunkIO_t), target :: ioall_albout
    type(ChunkIO_t), target :: io_albout(NSTATS)

    integer :: ichunk,jchunk,ic,jc,ii,jj
    integer :: id,n
    real(real64) :: val,byN,mean,sum,sumsq,xmin,xmax

    real, dimension(:,:), allocatable, target :: wta    ! Surmised landmask

    call chunker%init(IMK, JMK, IMH*2,JMH*2, 'forplot', 100, 10, 10, (/6,5/), outputs_dir=THIS_OUTPUTS_DIR)

    allocate(wta(chunker%chunk_size(1), chunker%chunk_size(2)))

    ! ====================== Input files
    call chunker%nc_open_input(ioall_albin, &
        INPUTS_URL, INPUTS_DIR, &
        'soilalbedo/', 'Carrer2014_soilalbedo_VIS_NIR_'//sLAI_YEAR//'_8day_6km.nc', &
        'soilalb_'//trim(sbands_modis(iband)), 0)
    do id=1,ndates
        call chunker%nc_reuse_var(ioall_albin, io_albin(id), (/1,1,id/))
    end do

    ! ======================== Output Files
    call clear_file_info(info)
    info%vname = 'albmodis_'//trim(sbands_modis(iband))
    info%long_name = 'Carrer soil albedo ('//trim(sbands_modis(iband))//' band)'
    info%units = '1'
    info%file_metadata_type = 'soilalbedo' !carrer
    !Here get grid name from IM and JM in chunker%ngrid
    !Hack until function to do this is written: hard-coded 6km into file name.
    call chunker%nc_create1_n(ioall_albout, weighting(wta,1d0,0d0), &
        'tmp/carrer/', &
        'soilalbedo_6km_Carrer2014_'//sLAI_YEAR//'ann_modis_'//trim(sbands_modis(iband)), info, &
        sstats, sstats)
    do id=1,NSTATS
        call chunker%nc_reuse_var(ioall_albout, io_albout(id), (/1,1,id/), &
            weighting(wta, 1d0, 0d0))
    end do

    call chunker%nc_check(rw=rw)
#ifdef JUST_DEPENDENCIES
    return
#endif

    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)

        call chunker%move_to(ichunk,jchunk)
        wta = 0

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            ! Compute overall NetCDF index of current cell
            ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
            jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1


            ! --------------- Compute statistics across days
            sum = 0
            sumsq = 0
            xmin = 1e30
            xmax = -1e30
            n = 0

            do id = 1,ndates
                val = io_albin(id)%buf(ic,jc)
                if (val /= FillValue) then
                    sum = sum + val
                    sumsq = sumsq + val*val
                    xmin = min(xmin, val)
                    xmax = max(xmax, val)
                    n = n + 1
                    wta(ic,jc) = 1
                end if
            end do

            if (wta(ic,jc) == 1) then
                io_albout(SMIN)%buf(ic,jc) = xmin
                io_albout(SMAX)%buf(ic,jc) = xmax
                byN = 1d0 / n
                mean = sum * byN
                io_albout(SMEAN)%buf(ic,jc) = mean
                io_albout(SSTD)%buf(ic,jc) = sqrt(byN*sumsq - mean*mean)
                io_albout(SNUM)%buf(ic,jc) = n
            else
                io_albout(SMIN)%buf(ic,jc) = FillValue
                io_albout(SMAX)%buf(ic,jc) = FillValue
                io_albout(SMEAN)%buf(ic,jc) = FillValue
                io_albout(SSTD)%buf(ic,jc) = FillValue
                io_albout(SNUM)%buf(ic,jc) = 0
            end if

if (io_albout(SSTD)%buf(ic,jc) > 1e10) then
    print *,'xBIG'
end if
                
        end do
        end do
        call chunker%write_chunks
    end do
    end do


    call chunker%close_chunks
end subroutine do_carrer_mean

end module B05_mod
! ================================================================
program carrer_mean
    use carrer_mod
    use B05_mod

implicit none


    integer :: err,ncid,dates_dimid,ndates,iband
    character(len=NF90_MAX_NAME) :: xname
    type(ReadWrites_t) :: rw
    integer :: nerr

    call rw%init(THIS_OUTPUTS_DIR, 'B05_carrer_mean', 20,20)

    ! Dummy open file to make sure it's downloaded
    nerr = 0
    err = download_input_file(nerr, &
        INPUTS_URL, INPUTS_DIR, &
        'soilalbedo/', 'Carrer2014_soilalbedo_VIS_NIR_'//sLAI_YEAR//'_8day_6km.nc')

    err = nf90_open(INPUTS_DIR//'soilalbedo/Carrer2014_soilalbedo_VIS_NIR_'//sLAI_YEAR//'_8day_6km.nc', NF90_NOWRITE, ncid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening file at beginning'
        STOP
    end if

    err = nf90_inq_dimid(ncid, 'dates', dates_dimid)
    err = nf90_inquire_dimension(ncid, dates_dimid, xname, ndates)
    if (err /= NF90_NOERR) then
        STOP
    end if

    err = nf90_close(ncid)

    do iband=1,NBANDS_MODIS
        call do_carrer_mean(rw, iband, ndates)
    end do
    call rw%write_mk

end program carrer_mean

