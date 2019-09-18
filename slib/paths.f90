module paths_mod
implicit none

! This file lists the input and output directories used by EntGVSD.
! It is intended that the user change these when setting up on a new
! system.


! Input files are read from one of four "input files roots", which are
! defined in `slib/paths.f90` (see below).  Each root has two locations:
! The `_DIR` version is where the original (compressed) files reside.
! The `_INPUT` version is where EntGVSD1 uncompresses those files to for
! its own use.

! ================ Input Files
! Original data input files (e.g. any file NOT created by some program
! in EntGVSD: Schaaf UMB files, scaled files, pre-processd.)  are
! never overwritten.  However, they come compressed and cannot be read
! directly.  When it needs to open a data file, EntGVSD decompresses
! it and stores a directory structure of decompressed input files.
! For every original input directory, there is a corresponding
! decompressed input directory.  The correspondence is described by
! pairs of path constants:
! 
!    XYZ_DIR is an original input directory
!    XYZ_INPUT is the corresponding decompressed input directory



CHARACTER (LEN=*), PARAMETER :: DATA_DIR= &
   '/home2/rpfische/entgvsd0/discover/Vegcover_1km/data/'
CHARACTER (LEN=*), PARAMETER :: DATA_INPUT= &
   '/home2/rpfische/git/entgvsd1/inputs/data/'
CHARACTER (LEN=*), PARAMETER :: DATA_INPUT_MANUAL= &
   '/home2/rpfische/git/entgvsd1/inputs_manual/'



CHARACTER (LEN=*), PARAMETER :: LAI3G_DIR= &
   '/home2/rpfische/entgvsd0/discover/Vegcover_1km/LAI3g/'
CHARACTER (LEN=*), PARAMETER :: LAI3G_INPUT= &
   '/home2/rpfische/git/entgvsd1/inputs/LAI3g/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR= &
   '/home2/rpfische/entgvsd0/discover/Vegcover_1km/BNU/lc_lai_giss/'
CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_INPUT= &
   '/home2/rpfische/git/entgvsd1/inputs/lc_lai_giss/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR= &
    '/home2/rpfische/entgvsd0/discover/Vegcover_1km/lc_lai_for_1kmx1km/'
CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_INPUT= &
   '/home2/rpfische/git/entgvsd1/inputs/lc_lai_giss/lc_lai_for_1kmx1km/'



! ================ Output Files
! Output files are written to multiple roots, listed here.

CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR= &
   '/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent/'

CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT16_DIR= &
   '/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent16/'

! ===========================================================

! 'L' or 'B'
CHARACTER (LEN=*), PARAMETER :: LAI_SOURCE = 'B'

! Name of main Fortran file
character (LEN=40) :: MAIN_PROGRAM_FILE

! ===========================================================

! Sa Simard RH100 | Sb Simard RH100-stdev | Sc Simard RH100-0.5*stdev | Ha GEDI+ICESat-2 data (some day) RH100 | Hb GEDI+ICESat-2 RH100-stdev | etc.
character*(*), parameter :: XHEIGHTSOURCE = 'Sa'


end module paths_mod
