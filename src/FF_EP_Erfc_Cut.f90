module FF_EasyEP_Erfc_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_Erfc_Cut
    real(dp), allocatable :: scaleTable(:,:)
    real(dp), allocatable :: alphaTable(:,:)

    contains
      procedure, pass :: Constructor => Constructor_EP_Erfc_Cut
      procedure, pass :: PairFunction => PairFunction_EP_Erfc_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_Erfc_Cut
  end type

  contains
  subroutine Constructor_EP_Erfc_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Erfc_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%scaleTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%alphaTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%scaleTable = 1E0_dp
    self%alphaTable = 1E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the Erfc Cut Pair Style"

  end subroutine

  function PairFunction_EP_Erfc_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Erfc_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: scale, alpha, r

    scale = self % scaleTable(atmType2, atmType1)
    alpha = self % alphaTable(atmType2, atmType1)
    r = sqrt(rsq)
    E_Pair = scale * erfc(alpha * r) / r

  end function

  subroutine ProcessIO_EP_Erfc_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_Erfc_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat
    integer :: type1, type2
    real(dp) :: scaleval, alphaval, rCut, rMin
  
    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut

      case("erfc")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) scaleval
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) alphaval
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) rMin
        scaleval = scaleval * inEngUnit
        alphaval = alphaval / inLenUnit
        rMin = rMin * inLenUnit

        self%scaleTable(type1, type2) = scaleval
        self%scaleTable(type2, type1) = scaleval

        self%alphaTable(type1, type2) = alphaval
        self%alphaTable(type2, type1) = alphaval

        self%rMinTable(type1, type2) = rMin**2
        self%rMinTable(type2, type1) = rMin**2

      case default
    end select

  end subroutine
end module
