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
character*(*), parameter :: XHEIGHTSOURCE = 'Sa'


    ! Donwload input files from here
    character(*), parameter :: INPUTS_URL= &
        'https://portal.nccs.nasa.gov/datashare/GISS/Ent_TBM/inputs/'

    ! Put downloaded input files here
    character(*), parameter :: INPUTS_DIR= &
        ENTGVSD_PROJECT_SOURCE_DIR//'inputs_new/'

    ! Any EntGVSD1-generated output files go here
    character(*), parameter :: OUTPUTS_DIR= &
        ENTGVSD_PROJECT_SOURCE_DIR//'outputs/'

CONTAINS

end module ent_params_mod
