!================================================================================
!================================================================================
module FF_EasyEP_Table
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_Table
!    real(dp), allocatable :: rMin(:)
!    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    

    contains
      procedure, pass :: Constructor => Constructor_EP_Table
      procedure, pass :: PairFunction => PairFunction_EP_Table
      procedure, pass ::  => PairFunction_EP_Table
      procedure, pass :: ProcessIO => ProcessIO_EP_Table
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EP_Table(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Table), intent(inout) :: self
    integer :: AllocateStat

    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the LJ/Cut Pair Style"

  end subroutine
  !=============================================================================+
  function PairFunction_EP_Table(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Table), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    integer :: highBin, lowBin
    real(dp) :: dr, r, rLow, rHigh

    dr = self%dR(atmType1, atmType2)
    

    r = sqrt(rsq)
    lowbin = floor(dr*r)
    highbin = ceiling(dr*r)

    rLow = real(lowbin,dp)/dr
    rHigh = real(highbin,dp)/dr

    E_Pair = 

  end function
  !=====================================================================
  subroutine ProcessIO_EP_Table(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_Table), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    character(len=150) :: filename
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, rCut, rMin
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut
      case("pairfile")
        call GetXCommand(line, filename, 2, lineStat)
      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      call CountCommands(line, nPar)
      select case(nPar)
        case(2)
          read(line, *) type1, rMin
          rMin = rMin * inLenUnit
          self%rMin(type1) = rMin 


          do jType = 1, nAtomTypes
            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case default
          lineStat = -1
          write(0,*) line
          stop "Unknown Command Found"
      end select
    endif

  end subroutine
  !=====================================================================
  subroutine PairFunction_(self, pairfile, atmtype1, atmtype2)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Table), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    character(len=150), intent(in) :: filename



  end subroutine
  !=====================================================================
end module
!=====================================================================
