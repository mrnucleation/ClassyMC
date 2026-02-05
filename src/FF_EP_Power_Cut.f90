module FF_EasyEP_Power_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_Power_Cut
    real(dp), allocatable :: ATable(:,:)
    real(dp), allocatable :: nTable(:,:)

    contains
      procedure, pass :: Constructor => Constructor_EP_Power_Cut
      procedure, pass :: PairFunction => PairFunction_EP_Power_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_Power_Cut
  end type

  contains
  subroutine Constructor_EP_Power_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Power_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%ATable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%nTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%ATable = 0E0_dp
    self%nTable = 12E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the Power Cut Pair Style"

  end subroutine

  function PairFunction_EP_Power_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Power_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: A, n, r

    A = self % ATable(atmType2, atmType1)
    n = self % nTable(atmType2, atmType1)
    r = sqrt(rsq)
    if(r > 0E0_dp) then
      E_Pair = A / (r ** n)
    else
      E_Pair = 0E0_dp
    endif

  end function

  subroutine ProcessIO_EP_Power_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_Power_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat
    integer :: type1, type2
    real(dp) :: Aval, nval, rCut, rMin
  
    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut

      case("power")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) Aval
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) nval
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) rMin
        Aval = Aval * inEngUnit
        rMin = rMin * inLenUnit

        self%ATable(type1, type2) = Aval
        self%ATable(type2, type1) = Aval

        self%nTable(type1, type2) = nval
        self%nTable(type2, type1) = nval

        self%rMinTable(type1, type2) = rMin**2
        self%rMinTable(type2, type1) = rMin**2

      case default
    end select

  end subroutine
end module
