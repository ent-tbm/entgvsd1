program test_ent
! Program tests the chunker infrastructure, a little bit.

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

