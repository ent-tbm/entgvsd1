
!  Program to assign 1kmx1km BNU LAI of selected DOY to EntPFTs

!------------------------------------------------------------------------

program lc_lai_doy

use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none


!***************************************************
!*      PREFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
character*3, parameter :: EntPFT_files1(19) = &
     (/ &
     '01_', &
     '02_', &
     '03_', &
     '04_', &
     '05_', &
     '06_', &
     '07_', &
     '08_', &
     '09_', &
     '10_', &
     '11_', &
     '12_', &
     '13_', &
     '14_', &
     '15_', &
     '16_', &
     '17_', &
     '18_', &
     '19_' &
     /)

!***************************************************
!*      SUFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
character*14, parameter :: EntPFT_files2(19) = &
     (/ &
     'ever_br_early ', &
     'ever_br_late  ', &
     'ever_nd_early ', &
     'ever_nd_late  ', &
     'cold_br_early ', &
     'cold_br_late  ', &
     'drought_br    ', &
     'decid_nd      ', &
     'cold_shrub    ', &
     'arid_shrub    ', &
     'c3_grass_per  ', &
     'c4_grass      ', &
     'c3_grass_ann  ', &
     'c3_grass_arct ', &
     'crops_c3_herb ', &
     'crops_c4_herb ', &
     'crops_woody   ', &
     'snow_ice      ', &
     'bare_sparse   ' &
     /)


integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees

integer, parameter :: X1km = 43200 !long at 1 km
integer, parameter :: Y1km = 21600 !lat at 1 km

integer, parameter :: IM1km = X1km !long at 1 km
integer, parameter :: JM1km = Y1km !lat at 1 km

integer, parameter :: longin = 1
integer, parameter :: latin = 1
integer, parameter :: longout = 1 !Should be same at longin
integer, parameter :: latout = 1

real*4, parameter :: undef = -1e30 !AViewer undef value      !## UPDATE ME

character*3, parameter :: DOY(2) = &
     (/ &
     "017","201" &
     /)

character*256 :: PathFilepre,PathFilepost

integer, parameter :: NUMLAYERS = 40  !## UPDATE ME

integer, parameter :: NUMLAYERSLC = 20  ! ONLY LC LAYERS
integer, parameter :: NUMTITLES = 40  ! ALL LC AND LAI TITLES

character*50 :: PATHin, PATHout
real*4 :: LAYERIN
real*4 :: LAYEROUT
real*4 :: OFFIA, DLATAH, DLATA, OFFIB, DLATB, DATMIS

character*256 :: filelai,fileout,filepft,filewater
character*20 :: inqvarin

real*4 :: LCIN_in, LAI, lai_lc

integer :: i, j, k, f, p, z

integer :: err,dimidx,dimidy,dimidz
integer :: ycoord, xcoord

integer, parameter :: ndoy = 2

type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t), target :: io_lai(ndoy),io_water
type(ChunkIO_t), target :: io_pft(19)
! Output files
type(ChunkIO_t) :: ioall_out(ndoy), io_out(NUMLAYERSLC,ndoy)


!integer :: startA(1),startB(2),countA(1),countB(2)
!integer :: startX(1),startY(1),countX(1),countY(1)
!integer :: start3d(4),count3d(4)
integer :: dd(4),varidx(12),varidy(12)
integer :: layer_indices(20)
character*17 :: layer_names(20)
integer :: ichunk,jchunk,ic,jc
real*4, dimension(:,:), pointer :: wbuf

call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 100, 120)

!* Input file.

!     Monthly LAI
!      PathFilepre= '../../data/LAI/LAI3gMonthly/'
!      PathFilepost = 'LAI3g_'//MONTH(7)//'_2004_1kmx1km.nc'
!      filelai = trim(PathFilepre)//trim(PathFilepost)
!      err = NF90_OPEN(filelai,NF90_WRITE,laifileid)
!      err = NF90_INQ_VARID(laifileid,'lai',laivarid)
!      err = NF90_INQ_VARID(laifileid,'lon',varidx)
!      err = NF90_INQ_VARID(laifileid,'lat',varidy)

!     DOY LAI

do k = 1,2
    call chunker%nc_open_gz(io_lai(k), DATA_DIR, DATA_INPUT, &
        'LAI/', 'global_30s_2004_'//DOY(k)//'.nc', 'lai', 1)
enddo

!     Water LC

!call chunker%nc_open_gz(io_water, LAI3G_DIR, LAI3G_INPUT, &
!    'EntMM_lc_laimax_1kmx1km/', 'water_lc.nc', 'water_lc')
call chunker%nc_open(io_water, LC_LAI_ENT_DIR, &
    'EntMM_lc_laimax_1kmx1km/', 'water_lc.nc', 'water',1)

!     ENTPFTLC
do k = 1,19
    call chunker%nc_open(io_pft(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(EntPFT_files1(k))//trim(EntPFT_files2(k))//'_lc.nc', &
        trim(EntPFT_files2(k)),1)
end do

! Cons up layer names and indices to write into our output
layer_indices(1) = 0
layer_names(1) = '0_water'
do k=1,19
    layer_indices(k+1) = k
    layer_names(k+1) = EntPFT_files1(k) // EntPFT_files2(k)
end do

!     Fileout
do k = 1,2
    call chunker%nc_create(ioall_out(k), &
        weighting(chunker%wta1, 1d0, 0d0), &    ! TODO: Scale by _lc; store an array of 2D array pointers
        'nc/', 'EntMM_lc_lai_'//DOY(k)//'_1kmx1km', 'EntPFT', &
        'LAI output of A02', 'm2 m-2', 'LAI', &
        layer_indices, layer_names)
    do p=1,20
        if (p==1) then
            wbuf => io_water%buf
        else
            wbuf => io_pft(p-1)%buf
        end if

        call chunker%nc_reuse_var(ioall_out(k), io_out(p,k), &
            (/1,1,p/), 'w', weighting(wbuf, 1d0,0d0))
    end do

enddo

call chunker%nc_check('A02_lc_lai_doy')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

!-----------------------------------------------------------------
!     Loop for every grid point

! Use these loop bounds for testing...
! it chooses a land area in Asia
#ifdef ENTGVSD_DEBUG
do jchunk = nchunk(2)*3/4,nchunk(2)*3/4+1
do ichunk = nchunk(1)*3/4,nchunk(1)*3/4+1
#else
do jchunk = 1,nchunk(2)
do ichunk = 1,nchunk(1)
#endif

    call chunker%move_to(ichunk,jchunk)

    do jc = 1,chunker%chunk_size(2)
    do ic = 1,chunker%chunk_size(1)

         
!**   LAI data ------------------------------------------------------------------------------

        do k = 1,2
        LAI = io_lai(k)%buf(ic,jc)
        
            do p = 1,NUMLAYERSLC


                if (p.eq.1) then
                    LCIN_in = io_water%buf(ic,jc)
                else
                    LCIN_in = io_pft(p-1)%buf(ic,jc)
                end if

                lai_lc = LCIN_in*LAI
                io_out(p,k)%buf(ic,jc) = lai_lc
            end do ! p=1,NUMLAYERSLC
        end do   ! k=1,2
    end do    ! ic
    end do    !jc

    call chunker%write_chunks

enddo    ! ichunk
enddo    ! jchunk

call chunker%close_chunks

end program lc_lai_doy
