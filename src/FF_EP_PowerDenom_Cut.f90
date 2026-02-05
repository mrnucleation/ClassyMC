module FF_EasyEP_PowerDenom_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_PowerDenom_Cut
    real(dp), allocatable :: ATable(:,:)
    real(dp), allocatable :: BTable(:,:)
    real(dp), allocatable :: DTable(:,:)
    real(dp), allocatable :: nTable(:,:)

    contains
      procedure, pass :: Constructor => Constructor_EP_PowerDenom_Cut
      procedure, pass :: PairFunction => PairFunction_EP_PowerDenom_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_PowerDenom_Cut
  end type

  contains
  subroutine Constructor_EP_PowerDenom_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_PowerDenom_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%ATable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%BTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%DTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%nTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%ATable = 1E0_dp
    self%BTable = 1E0_dp
    self%DTable = 1E0_dp
    self%nTable = 6E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the PowerDenom Cut Pair Style"

  end subroutine

  function PairFunction_EP_PowerDenom_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_PowerDenom_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: A, B, D, n, r, denom

    A = self % ATable(atmType2, atmType1)
    B = self % BTable(atmType2, atmType1)
    D = self % DTable(atmType2, atmType1)
    n = self % nTable(atmType2, atmType1)
    r = sqrt(rsq)
    denom = B + D * (r ** n)
    if(denom > 0E0_dp) then
      E_Pair = A / denom
    else
      E_Pair = 0E0_dp
    endif

  end function

  subroutine ProcessIO_EP_PowerDenom_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_PowerDenom_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat
    integer :: type1, type2
    real(dp) :: Aval, Bval, Dval, nval, rCut, rMin
  
    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut

      case("powerdenom")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) Aval
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) Bval
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) Dval
        call GetXCommand(line, command, 7, lineStat)
        read(command, *) nval
        call GetXCommand(line, command, 8, lineStat)
        read(command, *) rMin
        Aval = Aval * inEngUnit
        rMin = rMin * inLenUnit

        self%ATable(type1, type2) = Aval
        self%ATable(type2, type1) = Aval

        self%BTable(type1, type2) = Bval
        self%BTable(type2, type1) = Bval

        self%DTable(type1, type2) = Dval
        self%DTable(type2, type1) = Dval

        self%nTable(type1, type2) = nval
        self%nTable(type2, type1) = nval

        self%rMinTable(type1, type2) = rMin**2
        self%rMinTable(type2, type1) = rMin**2

      case default
    end select

  end subroutine
end module
