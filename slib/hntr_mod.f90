
module hntr4_mod
implicit none

    real(8), parameter :: M_PI  = 4 * atan (1d0)
contains

type HntrSpec_t
    integer :: im    ! Number of cells in east-west direction
    integer :: jm    ! Number of cells in north-south direction

    ! number (fraction) of cells in east-west direction from
    ! International Date Line (180) to western edge of cell IA=1
    real*8 :: offi

    ! minutes of latitude for non-polar cells on grid A
    real*8 :: dlat
end type HntrSpec_t

type HntrGrid_t
    type(HntrSpec_t) :: spec
    real*8, dimension(:), allocatable :: dxyp
end type HntrGrid_t

type HntrCalc_t
    type(HntrSpec_t) :: specA, specB
    real*8, dimension(:), allocatable :: dxypA, dxypB

    real*8, dimension(:), allocatable :: SINA,SINB
    real*8, dimension(:), allocatable :: FMIN,FMAX
    real*8, dimension(:), allocatable :: IMIN,IMAX
    real*8, dimension(:), allocatable :: GMIN,GMAX
    real*8, dimension(:), allocatable :: JMIN,JMAX

    ! DATMIS = missing data value inserted in output array B when
    ! cell (IB,JB) has integrated value 0 of WTA
    real*8 DATMIS;

    contains

    procedure :: regrid
    procedure :: partition_east_west
    procedure :: partition_north_south

end type HntrCalc_t

CONTAINS

function hntr_spec(im, jm, offi, dlat)
    integer, intent(IN) :: im,jm
    real*8, intent(IN) :: offi
    real*8, intent(IN) :: dlat
    type(HntrSpec_t), intent(OUT) :: hntr_spec
    ! ------ Locals

    ! ------ Construct the HntrSpec
    hntr_spec%im = im
    hntr_spec%jm = jm
    hntr_spec%offi = offi
    hntr_spec%dlat = dlat
end function hntr_spec
<
subroutine make_dxyp(dxyp, spec)
    real*8, double(:), allocatable, intent(INOUT) :: dxyp
    type(HntrSpec_t), intent(IN) :: spec
    ! -------- Locals
    real*8 :: dLON, dLAT
    real*8 :: SINS, SINN

    ! ------ Compute dxyp
    allocate(dxyp(jm)

    ! Calculate the sperical area of grid cells
    ! (on a radius=1 sphere)
    dLON = (2d0 * M_PI) / hntr_grid%spec%im;
    dLAT = M_PI / hntr_grid%spec%jm;
    do j=1,jm
        SINS = sin(dLAT*(j-jm/2-1));
        SINN = sin(dLAT*(j-jm/2));
        dxyp(j) = dLON * (SINN - SINS);
    end do
end subroutine make_dxyp

function hntr_calc(specB, specA, datmis)
    type(HntrSpec_t) :: specB, specA
    real*8 :: datmis
    ! -------- Locals

    allocate(hntr_calc%SINA(0:specA%jm))
    allocate(hntr_calc%SINB(0:specB%jm))
    allocate(hntr_calc%FMIN(specB%im))
    allocate(hntr_calc%FMAX(specB%im))
    allocate(hntr_calc%IMIN(specB%im))
    allocate(hntr_calc%IMAX(specB%im))
    allocate(hntr_calc%GMIN(specB%jm))
    allocate(hntr_calc%GMAX(specB%jm))
    allocate(hntr_calc%JMIN(specB%jm))
    allocate(hntr_calc%JMAX(specB%jm))

    hntr_calc%DATMIS = datmis

    this%specA = specA
    make_dxyp(this%dxypA, specA)
    this%specB = specB
    make_dxyp(this%dxypB, specB)


    call hntr_calc%partition_east_west
    call hntr_calc%partition_north_south
end function hntr_calc


! Partitions in east-west (I) direction
! Domain, around the globe, is scaled to fit from 0 to IMA*Bgrid.spec.im
subroutine partition_east_west(this)
    class(HntrCalc_t) :: this
    ! --------- Locals
    real*8 :: DIA
    integer :: IA, IB
    real*8 :: RIA, RIB
    integer :: IBp1

    DIA = specB%im  !  width of single A grid cell in scaled domain
    IA = 1
    double RIA = (IA+specA%offi - specA%im)*specB%im  !  scaled longitude of eastern edge
    IB  = specB%im
    do IBp1=1,specB%im
        RIB = (IBp1-1+specB%offi)*specA%im    !  scaled longitude of eastern edge
        do while (RIA < RIB)
            IA  = IA + 1
            RIA = RIA + DIA
        end do

        if (RIA == RIB) then
        ! Eastern edges of cells IA of grid A and IB of grid B coincide
            this%IMAX(IB) = IA
            this%FMAX(IB) = 0
            IA = IA + 1
            RIA = RIA + DIA
            this%IMIN(IBp1) = IA
            this%FMIN(IBp1) = 0
        else
        ! Cell IA of grid A contains western edge of cell IB of grid B
            this%IMAX(IB) = IA
            this%FMAX(IB) = (RIA-RIB)/DIA
            this%IMIN(IBp1) = IA
            this%FMIN(IBp1) = 1-FMAX(IB)
        end if
        IB = IBp1
    end do
    this%IMAX(specB%im) = this%IMAX(specB%im) + specA%im
end subroutine partition_east_west


subroutine partition_north_south(this)
    class(HntrCalc_t) :: this
    ! --------- Locals
    ! Convert minutes to radians
    real*8, parameter :: MIN_TO_RAD = (2. * M_PI) / (360d0*60d0)
    real*8 :: FJEQA, RJA
    real*8 :: FJEQB, RJB
    integer :: JA, JB

    ! ------------------------------------------------
    ! Partitions in the north-south (J) direction
    ! Domain is measured in minutes (1/60-th of a degree)
    FJEQA = .5*(1+specA%jm)
    do JA=1,specA%jm-1
        RJA = (JA + .5-FJEQA) * specA%dlat  !  latitude in minutes of northern edge
        this%SINA(JA) = sin(RJA * MIN_TO_RAD)
    end do
    this%SINA(0) = -1
    this%SINA(specA%jm)=  1

    ! -----------
    double FJEQB = .5*(1+specB%jm)
    do JB=1,specB%jm-1
        double RJB = (JB+.5-FJEQB)*specB%dlat  !  latitude in minutes of northern edge
        this%SINB(JB) = sin(RJB * MIN_TO_RAD)
    end do
    this%SINB(0) = -1
    this%SINB(specB%jm) = 1

    ! -----------
    this%JMIN(1) = 1
    this%GMIN(1) = 0
    JA = 1
    do JB=1,specB%jm-1
        do while (SINA(JA) < SINB(JB))
            JA = JA + 1
        end do

        if (SINA(JA) == SINB(JB)) then
            !  Northern edges of cells JA of grid A and JB of grid B coincide
            JMAX(JB) = JA
            GMAX(JB) = 0
            JA = JA + 1
            JMIN(JB+1) = JA
            GMIN(JB+1) = 0
        else
            ! Cell JA of grid A contains northern edge of cell JB of grid B
            JMAX(JB) = JA
            GMAX(JB) = SINA(JA) - SINB(JB)
            JMIN(JB+1) = JA
            GMIN(JB+1) = SINB(JB) - SINA(JA-1)
        end if
    end do
    JMAX(specB%jm) = specA%jm
    GMAX(specB%jm) = 0

end subroutine partition_north_south


subroutine regrid(this, WTA,B,A)
    class(HntrCalc_t) :: this
    real*4, dimension(:), intent(IN) :: WTA
    real*4, dimension(:), intent(INOUT) :: B
    real*4, dimension(:), intent(IN) :: A
    ! ------- Locals

!**** Local variables
      Integer :: IA,JA,IJA, IB,JB,IJB, IAMIN,IAMAX,JAMIN,JAMAX, IAREV
      Real*8  :: WEIGHT,VALUE,F,G
!****
!**** Interpolate the A grid onto the B grid
!****

    integer :: JB,IB,IJB
    real*8 :: JAMIN,JAMAX
    real*8 :: WEIGHT, VALUE
    integer :: IAMIN, IAMAX
    real*8 :: F,G

    Do 20 JB=1,specB%jm
        JAMIN = this%JMIN(JB)
        JAMAX = this%JMAX(JB)
        Do 20 IB=1,IMB
            IJB  = IB + IMB*(JB-1)
            WEIGHT= 0
            VALUE = 0
            IAMIN = this%IMIN(IB)
            IAMAX = this%IMAX(IB)
            Do JA=JAMIN,JAMAX
                G = SINA(JA)-SINA(JA-1)
                If (JA==JAMIN)  G = G - this%GMIN(JB)
                If (JA==JAMAX)  G = G - this%GMAX(JB)
                Do 10 IAREV=IAMIN,IAMAX
                    IA  = 1 + Mod(IAREV-1,IMA)
                    IJA = IA + IMA*(JA-1)
                    F = 1d0
                    If (IAREV==IAMIN) F = F - this%FMIN(IB)
                    If (IAREV==IAMAX) F = F - this%FMAX(IB)
                    WEIGHT = WEIGHT + F*G*WTA(IJA)
                    VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
                end do
            end do
            if (WEIGHT == 0) then
                B(IJB) = DATMIS
            else
                B(IJB) = VALUE/WEIGHT
            end if
        end do
   end do
end subroutine regrid

end module hntr4_mod
