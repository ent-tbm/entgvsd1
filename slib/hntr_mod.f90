
module hntr_mod
implicit none

    real(8), parameter :: M_PI  = 4 * atan (1d0)

type HntrSpec_t
    integer :: im    ! Number of cells in east-west direction
    integer :: jm    ! Number of cells in north-south direction

    ! number (fraction) of cells in east-west direction from
    ! International Date Line (180) to western edge of cell IA=1
    real*8 :: offi

    ! minutes of latitude for non-polar cells on grid A
    real*8 :: dlat
end type HntrSpec_t

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

    procedure :: regrid4
    procedure :: partition_east_west
    procedure :: partition_north_south

end type HntrCalc_t

    CONTAINS

function hntr_spec(im, jm, offi, dlat)
    integer, intent(IN) :: im,jm
    real*8, intent(IN) :: offi
    real*8, intent(IN) :: dlat
    type(HntrSpec_t) :: hntr_spec
    ! ------ Locals

    ! ------ Construct the HntrSpec
    hntr_spec%im = im
    hntr_spec%jm = jm
    hntr_spec%offi = offi
    hntr_spec%dlat = dlat
end function hntr_spec

subroutine make_dxyp(dxyp, spec)
    real*8, dimension(:), allocatable, intent(INOUT) :: dxyp
    type(HntrSpec_t), intent(IN) :: spec
    ! -------- Locals
    real*8 :: dLON, dLAT
    real*8 :: SINS, SINN
    integer :: j

    ! ------ Compute dxyp
    allocate(dxyp(spec%jm))

    ! Calculate the sperical area of grid cells
    ! (on a radius=1 sphere)
    dLON = (2d0 * M_PI) / spec%im;
    dLAT = M_PI / spec%jm;
    do j=1,spec%jm
        SINS = sin(dLAT*(j-spec%jm/2-1));
        SINN = sin(dLAT*(j-spec%jm/2));
        dxyp(j) = dLON * (SINN - SINS);
    end do
end subroutine make_dxyp

function hntr_calc(specB, specA, datmis)
    type(HntrSpec_t) :: specB, specA
    real*8 :: datmis
    type(HntrCalc_t) :: hntr_calc
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

    hntr_calc%specA = specA
    call make_dxyp(hntr_calc%dxypA, specA)
    hntr_calc%specB = specB
    call make_dxyp(hntr_calc%dxypB, specB)


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

    DIA = this%specB%im  !  width of single A grid cell in scaled domain
    IA = 1
    RIA = (IA+this%specA%offi - this%specA%im)*this%specB%im  !  scaled longitude of eastern edge
    IB  = this%specB%im
    do IBp1=1,this%specB%im
        RIB = (IBp1-1+this%specB%offi)*this%specA%im    !  scaled longitude of eastern edge
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
            this%FMIN(IBp1) = 1-this%FMAX(IB)
        end if
        IB = IBp1
    end do
    this%IMAX(this%specB%im) = this%IMAX(this%specB%im) + this%specA%im
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
    FJEQA = .5*(1+this%specA%jm)
    do JA=1,this%specA%jm-1
        RJA = (JA + .5-FJEQA) * this%specA%dlat  !  latitude in minutes of northern edge
        this%SINA(JA) = sin(RJA * MIN_TO_RAD)
    end do
    this%SINA(0) = -1
    this%SINA(this%specA%jm)=  1

    ! -----------
    FJEQB = .5*(1+this%specB%jm)
    do JB=1,this%specB%jm-1
        RJB = (JB+.5-FJEQB)*this%specB%dlat  !  latitude in minutes of northern edge
        this%SINB(JB) = sin(RJB * MIN_TO_RAD)
    end do
    this%SINB(0) = -1
    this%SINB(this%specB%jm) = 1

    ! -----------
    this%JMIN(1) = 1
    this%GMIN(1) = 0
    JA = 1
    do JB=1,this%specB%jm-1
        do while (this%SINA(JA) < this%SINB(JB))
            JA = JA + 1
        end do

        if (this%SINA(JA) == this%SINB(JB)) then
            !  Northern edges of cells JA of grid A and JB of grid B coincide
            this%JMAX(JB) = JA
            this%GMAX(JB) = 0
            JA = JA + 1
            this%JMIN(JB+1) = JA
            this%GMIN(JB+1) = 0
        else
            ! Cell JA of grid A contains northern edge of cell JB of grid B
            this%JMAX(JB) = JA
            this%GMAX(JB) = this%SINA(JA) - this%SINB(JB)
            this%JMIN(JB+1) = JA
            this%GMIN(JB+1) = this%SINB(JB) - this%SINA(JA-1)
        end if
    end do
    this%JMAX(this%specB%jm) = this%specA%jm
    this%GMAX(this%specB%jm) = 0

end subroutine partition_north_south


! Interpolate the A grid onto the B grid
subroutine regrid4(this, B,A,WTA,jb0,njb)
    class(HntrCalc_t) :: this
    real*4, dimension(:,:), intent(INOUT) :: B
    real*4, dimension(:,:), intent(IN) :: A
    real*4, dimension(:,:), intent(IN) :: WTA
    ! ------- Locals

    integer :: JB,IB
    integer :: JA,IA
    integer :: JAMIN,JAMAX
    real*8 :: WEIGHT, VALUE
    integer :: IAMIN, IAMAX
    real*8 :: F,G
    integer :: IAREV

!    Do JB=1,this%specB%jm
    ja0 = this%JMIN(jb0)
    Do JB=jb0,jbo+njb-1
        JAMIN = this%JMIN(JB)
        JAMAX = this%JMAX(JB)
        Do IB=1,this%specB%im
            WEIGHT= 0
            VALUE = 0
            IAMIN = this%IMIN(IB)
            IAMAX = this%IMAX(IB)
            Do JA=JAMIN,JAMAX
                G = this%SINA(JA)-this%SINA(JA-1)
                If (JA==JAMIN)  G = G - this%GMIN(JB)
                If (JA==JAMAX)  G = G - this%GMAX(JB)
                Do IAREV=IAMIN,IAMAX
                    IA  = 1 + Mod(IAREV-1,this%specA%im)
                    F = 1d0
                    If (IAREV==IAMIN) F = F - this%FMIN(IB)
                    If (IAREV==IAMAX) F = F - this%FMAX(IB)
                    WEIGHT = WEIGHT + F*G*WTA(IA,JA)
                    VALUE  = VALUE  + F*G*WTA(IA,JA)*A(IA,JA-ja0+1)
                end do
            end do
            print *,ib,jb,weight
            if (WEIGHT == 0) then
                B(IB,JB-jb0+1) = this%DATMIS
            else
                B(IB,JB-jb0+1) = VALUE/WEIGHT
            end if
        end do
   end do
end subroutine regrid4

end module hntr_mod
