!================================================================================
! Easy Pair Template for general pair forcefields. The purpose of this class is
! to act as an inheritable object that can allow the user to quickly set up
! any standard pair distance based forcefields by simply creating a child object
! from this class. By overwriting the PairFunction class and parameter settings
! the procedures for all cut-off based forcefields are automatically defined.
!================================================================================
module FF_EasyEP_LJ_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_LJ_Cut
!    real(dp), allocatable :: rMin(:)
!    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    real(dp), allocatable :: eps(:)
    real(dp), allocatable :: sig(:)

    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)

    contains
      procedure, pass :: Constructor => Constructor_EP_LJ_Cut
      procedure, pass :: PairFunction => PairFunction_EP_LJ_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_LJ_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EP_LJ_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%eps(1:nAtomTypes), stat = AllocateStat)
    allocate(self%sig(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)

    allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%sigTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%eps = 4E0_dp
    self%sig = 1E0_dp
    self%rMin = 0.5E0_dp

    self%epsTable = 4E0_dp
    self%sigTable = 1E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the LJ/Cut Pair Style"

  end subroutine
  !=============================================================================+
  function PairFunction_EP_LJ_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: ep, sig_sq, LJ

    ep = self % epsTable(atmType2, atmType1)
    sig_sq = self % sigTable(atmType2, atmType1)
    LJ = (sig_sq/rsq)
    LJ = LJ * LJ * LJ
    E_Pair = ep * LJ * (LJ-1E0_dp)

  end function
  !=====================================================================
  subroutine ProcessIO_EP_LJ_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
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

      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      call CountCommands(line, nPar)
      select case(nPar)
        case(4)
          read(line, *) type1, ep, sig, rMin
          ep = ep * inEngUnit
          sig = sig * inLenUnit
          rMin = rMin * inLenUnit

          self%eps(type1) = ep 
          self%sig(type1) = sig 
          self%rMin(type1) = rMin 


          do jType = 1, nAtomTypes
            self%epsTable(type1, jType) = 4E0_dp * sqrt(ep * self%eps(jType))
            self%epsTable(jType, type1) = 4E0_dp * sqrt(ep * self%eps(jType))

            self%sigTable(type1, jType) = (0.5E0_dp * (sig + self%sig(jType)))**2
            self%sigTable(jType, type1) = (0.5E0_dp * (sig + self%sig(jType)))**2

            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case(5)
          read(line, *) type1, type2, ep, sig, rMin
          self%epsTable(type1, type2) = ep
          self%epsTable(type2, type1) = ep

          self%sigTable(type1, type2) = sig
          self%sigTable(type2, type1) = sig

          self%rMinTable(type1, type2) = rMin**2
          self%rMinTable(type2, type1) = rMin**2


        case default
          lineStat = -1
      end select
    endif

  end subroutine
  !=====================================================================
end module
!=====================================================================
