module paths_mod


implicit none


! Hardcoded directories for input / output

! ---------- Original values (Carlos on discover)
!      --- Input (read-only)
!      CHARACTER (LEN=*), PARAMETER :: DATA_DIR='../../data/'
!      !CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR='../lc_lai_giss/'
!      CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR=
!     &     '../../lc_lai_for_1km1km/'
!      --- Output
!      CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR='../lc_lai_ent/'

! Values for Elizabeth's work on gibbs
!      --- Input (read-only)
CHARACTER (LEN=*), PARAMETER :: DATA_DIR= &
   '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/data/'

CHARACTER (LEN=*), PARAMETER :: DATA_INPUT= &
   '/home2/rpfische/git/entgvsd/inputs/data/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR= &
   '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/BNU/'// &
   'lc_lai_giss'

CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_INPUT= &
   '/home2/rpfische/git/entgvsd/inputs/lc_lai_giss/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR= &
    '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/'// &
    'lc_lai_for_1kmx1km/'
CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_INPUT= &
   '/home2/rpfische/git/entgvsd/inputs/lc_lai_giss/lc_lai_for_1kmx1km/'




!      --- Output
CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR= &
   '/home2/rpfische/git/entgvsd/bnu/lc_lai_ent/'

CHARACTER (LEN=*), PARAMETER :: TEMPLATE_DIR= &
   '/home2/rpfische/git/entgvsd/templates/'


! Place of output files originally written by carlos
! Used to create templates.

CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_ORIG= &
   '/home2/rpfische/entgvsd0_orig/discover/Vegcover_1km/BNU/lc_lai_ent'

! ===========================================================

end module paths_mod
