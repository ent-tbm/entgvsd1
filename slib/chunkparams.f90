module chunkparams_mod
implicit none

    ! ============== Parameters that might be made more variable

    ! Combine C3 and C4 crops into one PFT, for ModelE
    logical, parameter :: combine_crops_c3_c4 = .true.
    ! Split the bare soil into dark and light, to get the right albedo
    logical, parameter :: split_bare_soil = .true.


CONTAINS

end module chunkparams_mod
