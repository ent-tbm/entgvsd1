! Program computes min,max,mean,std of Carrer soil albedo
! Elizabeth Fischer <elizabeth.fischer@columbia.edu>
! August 15, 2019

module A01a_mod

    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use geom_mod
    use assign_laimax_mod
    use carrer_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


CONTAINS

subroutine do_carrer_mean(iband, ndates)
    integer, intent(IN) :: iband
    integer, intent(IN) :: ndates
    ! ------------ Local Vars

    type(Chunker_t) :: chunker
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

    call chunker%init(IMK, JMK, IMH*2,JMH*2, 'forplot', 100, 10, 10, (/6,5/))

    allocate(wta(chunker%chunk_size(1), chunker%chunk_size(2)))

    ! ====================== Input files
    call chunker%nc_open(ioall_albin, DATA_INPUT_MANUAL, 'carrer/', 'carrer.nc', &
        'soilalb_'//trim(sbands_modis(iband)), 0)
    do id=1,ndates
        call chunker%nc_reuse_var(ioall_albin, io_albin(id), (/1,1,id/))
    end do

    ! ======================== Output Files
    call chunker%nc_create(ioall_albout, weighting(wta,1d0,0d0), &
        'carrer/', &
        'albmodis_'//trim(sbands_modis(iband)), &
        'albmodis_'//trim(sbands_modis(iband)), &
        'Soil albedo ('//trim(sbands_modis(iband))//' band)', &
        '1', sstats, sstats)
    do id=1,NSTATS
        call chunker%nc_reuse_var(ioall_albout, io_albout(id), (/1,1,id/), &
            weighting(wta, 1d0, 0d0))
    end do

    call chunker%nc_check(trim(MAIN_PROGRAM_FILE)//'_'//trim(sbands_modis(iband)))

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

end module A01a_mod
! ================================================================
program carrer_mean
    use carrer_mod
    use A01a_mod

implicit none


    integer :: err,ncid,dates_dimid,ndates,iband
    character(len=NF90_MAX_NAME) :: xname

    MAIN_PROGRAM_FILE = 'A01a_carrer_mean'
    err = nf90_open(DATA_INPUT_MANUAL//'carrer/carrer.nc', NF90_NOWRITE, ncid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening file at beginning'
        return
    end if

    err = nf90_inq_dimid(ncid, 'dates', dates_dimid)
    err = nf90_inquire_dimension(ncid, dates_dimid, xname, ndates)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error reading ndates'
        return
    end if

    err = nf90_close(ncid)

    do iband=1,NBANDS_MODIS
        call do_carrer_mean(iband, ndates)
    end do

print *,'ndates',ndates

end program carrer_mean

