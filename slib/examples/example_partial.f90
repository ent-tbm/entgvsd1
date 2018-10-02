program example_hntr

use hntr_mod
implicit none

    type(HntrSpec_t) :: g4,g8
    type(HntrCalc_t) :: hc
    real*4 :: B8(16,8), A4(8,4)
    real*4 :: WTA4(8,4)

    integer :: j

    g4 = hntr_spec(3,3, 2d0, 180d0*60d0/3)
    g8 = hntr_spec(6, 6, 4d0, 180d0*60d0/6)

    hc = hntr_calc(g8, g4, 0d0)

    WTA4 = 1d0
    do j=1,g4%im
        a4(j,:)=j
    end do
    b8 = 0d0/0d0
    call hc%regrid4(B8, A4, WTA4)

    do j=1,g8%jm
        print *,B8(:,j)
    end do

    print *,'--------------'
end program example_hntr
