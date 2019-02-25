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
integer, parameter :: NUMLAYERSLC = 20  ! ONLY LC LAYERS
character*3, parameter :: EntPFT_files1(NUMLAYERSLC) = &
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
character*14, parameter :: EntPFT_files2(NUMLAYERSLC) = &
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

    character*50, parameter :: EntPFT_title(NUMLAYERSLC) = &
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

integer, parameter :: NUMTITLES = 40  ! ALL LC AND LAI TITLES

character*50 :: PATHin, PATHout
real*4 :: LAYERIN
real*4 :: LAYEROUT
real*4 :: OFFIA, DLATAH, DLATA, OFFIB, DLATB, DATMIS

!character*256 :: filelai,fileout,filepft,filewater
character*20 :: inqvarin

real*4 :: LC(IM,JM), LAI(IM,JM), lailc(IM,JM)
real*4 :: lon(IM),lat(JM)


real*8 :: CHECKSUM

integer :: ichunk,jchunk,ic,jc,ii,jj,k, f, p, z

integer :: err,dimidx,dimidy

type(Chunker_t) :: chunker
! Input files
type(ChunkIO_t) :: io_lai
type(ChunkIO_t) :: io_lc(NUMLAYERSLC)
! Output files
type(ChunkIO_t) :: io_laiout(NUMLAYERSLC), io_err(NUMLAYERSLC)
type(ChunkIO_t) :: io_checksum_lclai

integer :: start2d(2),count2d(2)
integer :: startX(1),startY(1),countX(1),countY(1)

call chunker%init(IM, JM, IMH*2,JMH*2, 'qxq', 100, 120)

! ================= Input Files
!      LAI max
call chunker%nc_open_gz(io_lai, DATA_DIR, DATA_INPUT, &
    'LAI/', 'global_30s_2004_max.nc', 'lai', 1)


! --- ENTPFTLC: Open outputs written by A00
do k = 1,NUMLAYERSLC
    call chunker%nc_open(io_lc(k), LC_LAI_ENT_DIR, &
        'EntMM_lc_laimax_1kmx1km/', trim(EntPFT_files1(k))//trim(EntPFT_files2(k))//'_lc.nc', &
        trim(EntPFT_files2(k)), 1)
enddo

! ================= Output Files
do k = 1,NUMLAYERSLC
    call chunker%nc_create(io_laiout(k),  weighting(io_lc(k)%buf,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(EntPFT_files1(k))//trim(EntPFT_files2(k))//'_lai', &
        trim(EntPFT_files2(k)), &
        EntPFT_title(k), 'm2 m-2', TITLE_LAI)

    call chunker%nc_create(io_err(k),  weighting(io_lc(k)%buf,1d0,0d0), &
        'EntMM_lc_laimax_1kmx1km/', &
        trim(EntPFT_files1(k))//trim(EntPFT_files2(k))//'_err', &
        trim(EntPFT_files2(k)), &
        EntPFT_title(k), 'm2 m-2', TITLE_LAI)
enddo

call chunker%nc_create(io_checksum_lclai,  weighting(chunker%wta1,1d0,0d0), &
    'EntMM_lc_laimax_1kmx1km/checksum_lclai/', &
    'lclai', &
    'Sum(LC*LAI) - LAI_orig == 0', 'm2 m-2', 'Sum of LC*LAI')



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

        ! Compute <original lai> - sum(LC*LAI), should equal 0
        CHECKSUM = 0d0
        do p = 1,NUMLAYERSLC

            ! If this LC type does NOT participate in this 1km grid cell,
            ! then LAI needs to be zero.
            if ((io_lc(p)%buf(ic,jc) <= 0d0).or.(io_lc(p)%buf(ic,jc) == undef)) then
                ! If lc==0, fix LAImax=0 (enforce the invariant)
                io_laiout(p)%buf(ic,jc) = 0d0
            else    ! non-zero LC type
                ! Use the single LAI for all LC types in this gridcell.
                io_laiout(p)%buf(ic,jc) = io_lai%buf(ic,jc)

                ! Problem if lc>0 but LAImax==0
                if (io_laiout%buf(ic,jc)==0) then
                    io_err(p) = io_lc(p)%buf(ic,jc)
                end if
            end if

            CHECKSUM = CHECKSUM + &
                io_lc(p)%buf(ic,jc) * io_laiout(p)%buf(ic,jc)
        end do
        CHECKSUM = CHECKSUM - io_lai%buf(ic,jc)
        io_checksum_lclai%buf(ic,jc) = CHECKSUM

    end do
    end do

    call chunker%write_chunks

end do
end do

call chunker%close_chunks

end program lc_laimax
