! Dimensions used in the Carrer soilalbedo file
!     <inputs>/soilalbedo/Carrer2014_soilalbedo_VIS_NIR_<LAI_YEAR>_8day_6km.nc
! Used by B05_carrer_mean, to read that file.  These MUST match what
! is in the file.
module carrer_mod

    ! MODIS bands: VIS (visible) or NIR (near infrared)
    integer, parameter :: VIS_MODIS = 1
    integer, parameter :: NIR_MODIS = 2
    integer, parameter :: NBANDS_MODIS = 2
    character*3, parameter :: sbands_modis(NBANDS_MODIS) = &
        (/"VIS", "NIR"/)

    ! Type of staticical summary:
    !     max, mean, min, standard deviation, number of elements
    integer, parameter :: SMAX = 1
    integer, parameter :: SMEAN = 2
    integer, parameter :: SMIN = 3
    integer, parameter :: SSTD = 4
    integer, parameter :: SNUM = 5
    integer, parameter :: NSTATS = 5
    integer, parameter :: NSTATS_FILLABLE = 3   ! Stats that will receive Poisson filling
    character*4, parameter :: sstats(NSTATS) = &
        (/"MAX ", "MEAN", "MIN ", "STD ", "N   "/)

end module carrer_mod
