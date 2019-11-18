program test_ent
!     Read in MODIS GISS layer and Monfreda 0.5x0.5,1x1,or 4x5 degree files, 
!     and convert veg types to Ent PFTs as GISS layers.
!     To convert to Ent data structures with mixed canopies, use program
!     in Ent repository.
!     Output is files prefixed "EntMM" for Ent-MODIS-Monfreda.
use netcdf
use chunker_mod
use ent_params_mod
use paths_mod
use ent_labels_mod
use geom_mod

implicit none

character*(8) :: s
type(Chunker_t) :: chunker
integer :: err

    call chunker%init(IM1km, JM1km, IMH*2,JMH*2, 'forplot', 100, 120, 10)

    err = download_input_file(chunker, INPUTS_URL, INPUTS_DIR, &
        '', 'hello.txt')
    print *,err

    err = download_input_file(chunker, INPUTS_URL, INPUTS_DIR, &
        '', 'hellox.txt')
    print *,err

end program test_ent

