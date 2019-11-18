module ent_params_mod
    use entgvsd_config_mod
implicit none

! ============== Parameters that might be made more variable

! Combine C3 and C4 crops into one PFT, for ModelE?
logical, parameter :: combine_crops_c3_c4 = .true.
! Split the bare soil into dark and light, to get the right albedo
logical, parameter :: split_bare_soil = .true.


! Source of LAI data
! 'L' = LAI3g; 'B' = BNU (Beijing Normal University)
CHARACTER (LEN=*), PARAMETER :: LAI_SOURCE = 'B'

! Sa Simard RH100 | Sb Simard RH100-stdev | Sc Simard RH100-0.5*stdev | Ha GEDI+ICESat-2 data (some day) RH100 | Hb GEDI+ICESat-2 RH100-stdev | etc.
acter*(*), parameter :: XHEIGHTSOURCE = 'Sa'

! =================== Paths

! Donwload input files from here
character(*), parameter :: INPUTS_URL= &
    'https://portal.nccs.nasa.gov/datashare/GISS/Ent_TBM/inputs/'

! Put downloaded input files here
character(*), parameter :: INPUTS_DIR= &
    ENTGVSD_PROJECT_SOURCE_DIR//'inputs_new/'

! Any EntGVSD1-generated output files go here
character(*), parameter :: OUTPUTS_DIR= &
    ENTGVSD_PROJECT_SOURCE_DIR//'outputs/'


! ======== Some standard ModelE (and related) resolutions for lon/lat grids

! 1kmx1km grid
integer, parameter :: X1km = 43200 !long at 1 km
integer, parameter :: Y1km = 21600 !lat at 1 km
integer, parameter :: IM1km = X1km !long at 1 km
integer, parameter :: JM1km = Y1km !lat at 1 km

integer, parameter :: IMH = 720 !long at 0.5 degrees
integer, parameter :: JMH = 360 !lat at 0.5 degrees

integer, parameter :: IMK = 7200 !long at 0.05 degrees ~ 1 km
integer, parameter :: JMK = 3600 !lat at 0.05 degrees ~ 1 im

integer, parameter :: IM1 = 360 !long at 1 degrees
integer, parameter :: JM1 = 180 !lat at 1 degrees
integer, parameter :: IM2 = 144 !long at 2.5 degrees
integer, parameter :: JM2 = 90 !lat at 2 degrees
integer, parameter :: IM4X5 = 72 !long at 5 degrees
integer, parameter :: JM4X5 = 46 !lat at 4 degrees

! Lo-res version (what we trim on)
integer, parameter :: IMLR=IMH
integer, parameter :: JMLR=JMH

CONTAINS

end module ent_params_mod
