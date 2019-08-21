module carrer_mod

    integer, parameter :: VIS_MODIS = 1
    integer, parameter :: NIR_MODIS = 2
    integer, parameter :: NBANDS_MODIS = 2
    character*3, parameter :: sbands_modis(NBANDS_MODIS) = &
        (/"VIS", "NIR"/)

    integer, parameter :: SMAX = 1
    integer, parameter :: SMEAN = 2
    integer, parameter :: SMIN = 3
    integer, parameter :: SSTD = 4
    integer, parameter :: SNUM = 5
    integer, parameter :: NSTATS = 5
    character*4, parameter :: sstats(NSTATS) = &
        (/"MAX ", "MEAN", "MIN ", "STD ", "N   "/)

end module carrer_mod
