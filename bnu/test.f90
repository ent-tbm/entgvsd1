program test_ent
!     Read in MODIS GISS layer and Monfreda 0.5x0.5,1x1,or 4x5 degree files, 
!     and convert veg types to Ent PFTs as GISS layers.
!     To convert to Ent data structures with mixed canopies, use program
!     in Ent repository.
!     Output is files prefixed "EntMM" for Ent-MODIS-Monfreda.
use netcdf
use chunker_mod
use chunkparams_mod
use paths_mod
use ent_labels_mod

implicit none

character*(8) :: s

    s = '12345678'
    s = ''
    print *,'"',trim(s),'"'
    print *,s==''
    print *,'Test program succeeded!'

end program test_ent

