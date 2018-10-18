!  Program to assign 1kmx1km BNU LAImax to EntPFTs lc

!------------------------------------------------------------------------

program lc_laimax

use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod

 ! Read in GISS layer 0.5x0.5 degree files, and use HNTRP* to 
 ! interpolate to coarser resolutions.
implicit none

!***************************************************
!*      PREFIX OF ENTPFTS FILES FOR LC AND LAI     *
!***************************************************
character*3, parameter :: EntPFT_files1(20) = &
     (/ &
     '   ', &
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
character*14, parameter :: EntPFT_files2(20) = &
     (/ &
     'water         ', &
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

    character*50, parameter :: EntPFT_title(20) = &
         (/ &
         '0 - water                                       ', &
         '1 - evergreen broadleaf early successional      ', &
         '2 - evergreen broadleaf late successional       ', &
         '3 - evergreen needleleaf early successional     ', &
         '4 - evergreen needleleaf late successional      ', &
         '5 - cold deciduous broadleaf early successional ', &
         '6 - cold deciduous broadleaf late successional  ', &
         '7 - drought deciduous broadleaf                 ', &
         '8 - deciduous needleleaf                        ', &
         '9 - cold adapted shrub                          ', &
         '10 - arid adapted shrub                         ', &
         '11 - C3 grass perennial                         ', &
         '12 - C4 grass                                   ', &
         '13 - C3 grass - annual                          ', &
         '14 - arctic C3 grass                            ', &
         '15 - crops C3 herb                              ', &
         '16 - crops C4 herb                              ', &
         '17 - crops woody                                ', &
         '18 - Permanent snow/ice                         ', &
         '19 - Bare or sparsely vegetated, urban          ' &
         /)

integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees

integer, parameter :: XM = 43200 !long at 1 km
integer, parameter :: YM = 21600 !lat at 1 km

integer, parameter :: IM = XM !long at 1 km
integer, parameter :: JM = YM !lat at 1 km

integer, parameter :: longi = XM
integer, parameter :: latin = YM

integer, parameter :: NUMLAYERS = 40  !## UPDATE ME

integer, parameter :: NUMLAYERSLC = 20  ! ONLY LC LAYERS
integer, parameter :: NUMTITLES = 40  ! ALL LC AND LAI TITLES

character*50 :: PATHin, PATHout
real*4 :: LAYERIN
real*4 :: LAYEROUT
real*4 :: OFFIA, DLATAH, DLATA, OFFIB, DLATB, DATMIS

!character*256 :: filelai,fileout,filepft,filewater
character*20 :: inqvarin

real*4 :: LC(IM,JM), LAI(IM,JM), lailc(IM,JM)
real*4 :: lon(IM),lat(JM)

integer :: ichunk,jchunk,ic,jc,ii,jj,k, f, p, z

integer :: err,dimidx,dimidy

type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t) :: io_lai
type(ChunkIO_t) :: io_lc(20)
! Output files
type(ChunkIO_t) :: io_out(20)

integer :: start2d(2),count2d(2)
integer :: startX(1),startY(1),countX(1),countY(1)

call chunker%init(IM, JM, IMH*2,JMH*2, 100, 120)

!      LAI max
call chunker%nc_open_gz(io_lai, DATA_DIR, DATA_INPUT, &
    'LAI/', 'global_30s_2004_max.nc', 'lai', 1)


! --- ENTPFTLC: Open outputs written by A00
do k = 1,20
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', trim(EntPFT_files1(k))//trim(EntPFT_files2(k))//'_lc.nc', &
        trim(EntPFT_files2(k)), 1)
enddo

!     Files out
do k = 1,20
    call chunker%nc_create(io_out(k),  weighting(io_lc(k)%buf,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(EntPFT_files1(k))//trim(EntPFT_files2(k))//'_lai', &
        trim(EntPFT_files2(k)), &
        EntPFT_title(k), 'm2 m-2', TITLE_LAI)
enddo

!-----------------------------------------------------------------
!     Loop for every pft and  grid point


! Quit if we had any problems opening files
call chunker%nc_check('A01_lc_laimax')
#ifdef JUST_DEPENDENCIES
stop 0
#endif

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

        ! Compute overall NetCDF index of current cell
        ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
        jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

        do p = 1,NUMLAYERSLC
            !write(*,*) 'processing ',EntPFT_files2(p), p

            !**   LAI data, lon lat from LAI file                

            io_out(p)%buf(ic,jc) = io_lc(p)%buf(ic,jc) * io_lai%buf(ic,jc)
        end do
    end do
    end do

    call chunker%write_chunks

end do
end do

end program lc_laimax
      
      
      

