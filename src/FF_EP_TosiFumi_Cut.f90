!================================================================================
! Tosi-Fumi Potential (Born-Mayer-Huggins based) with dispersion terms.
! The potential form is:
!   E(r) = A * exp(-r/p) - C/r^6
! where:
!   A = pre-exponential repulsion coefficient
!   p = softness parameter (characteristic length)
!   C = dispersion coefficient (C6 term)
!================================================================================
module FF_EasyEP_TosiFumi_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_TosiFumi_Cut
    ! Parameters stored as 2D arrays indexed by atom type pairs
    real(dp), allocatable :: A(:,:)       ! Pre-exponential factor
    real(dp), allocatable :: rho(:,:)     ! Softness parameter (p in the formula)
    real(dp), allocatable :: C6(:,:)      ! Dispersion coefficient
    real(dp), allocatable :: C8(:,:)      ! Dispersion coefficient
    real(dp), allocatable :: EShift(:,:)  ! Energy shift at cutoff
    logical :: shifted = .false.
    contains
      procedure, pass :: Constructor => Constructor_EP_TosiFumi_Cut
      procedure, pass :: PairFunction => PairFunction_EP_TosiFumi_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_TosiFumi_Cut
      procedure, pass :: ComputeEShift => ComputeEShift_EP_TosiFumi_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EP_TosiFumi_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_TosiFumi_Cut), intent(inout) :: self
    integer :: AllocateStat


    ! Allocate parameter arrays
    allocate(self%A(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
    allocate(self%rho(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
    allocate(self%C6(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
    allocate(self%C8(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
    allocate(self%EShift(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
    allocate(self%rMin(1:nAtomTypes), stat=AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)

    IF (AllocateStat /= 0) error STOP "Allocation error in EP_TosiFumi_Cut"

    ! Initialize to zero
    self%A = 0E0_dp
    self%rho = 1E0_dp  ! Avoid division by zero
    self%C6 = 0E0_dp
    self%C8 = 0E0_dp
    self%EShift = 0E0_dp
    self%rMin = 0.5E0_dp
    self%rMinTable = 0.25E0_dp

    ! Default cutoff
    self%rCut = 10E0_dp
    self%rCutSq = self%rCut**2

  end subroutine
  !======================================================================+
  function PairFunction_EP_TosiFumi_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_TosiFumi_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: r, r6, r8, A_ij, rho_ij, C6_ij, C8_ij
    real(dp) :: E_Rep, E_Disp

    r = sqrt(rsq)
    r6 = rsq * rsq * rsq  ! r^6
    r8 = r6* rsq  ! r^8

    ! Get parameters for this atom type pair
    A_ij = self%A(atmtype1, atmtype2)
    rho_ij = self%rho(atmtype1, atmtype2)
    C6_ij = self%C6(atmtype1, atmtype2)
    C8_ij = self%C8(atmtype1, atmtype2)

    ! Tosi-Fumi potential: A * exp(-r/p) - C/r^6
    E_Rep = A_ij * exp(-r / rho_ij)
    E_Disp = -C6_ij / r6
    E_Disp = E_Disp - C8_ij / r8
    E_Pair = E_Rep + E_Disp

    ! Apply energy shift if enabled
    if (self%shifted) then
      E_Pair = E_Pair - self%EShift(atmtype1, atmtype2)
    endif

  end function
  !=====================================================================
  subroutine ProcessIO_EP_TosiFumi_Cut(self, line)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use Common_MolInfo, only: nAtomTypes
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_TosiFumi_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: logicVal
    integer :: jType, iType, lineStat
    real(dp) :: realVal, A_val, rho_val, C6_val, C8_val

    call GetXCommand(line, command, 1, lineStat)
    select case(trim(adjustl(command)))

      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) realVal
        self%rCut = realVal * inLenUnit
        self%rCutSq = self%rCut * self%rCut

      case("shifted")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) logicVal
        self%shifted = logicVal

      case("pair")
        ! Format: pair iType jType A rho C6
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) iType
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) jType
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) A_val
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) rho_val
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) C6_val
        call GetXCommand(line, command, 7, lineStat)
        read(command, *) C8_val
        ! Apply units and store symmetrically
        A_val = A_val * inEngUnit
        rho_val = rho_val * inLenUnit
        C6_val = C6_val * inEngUnit * (inLenUnit**6)
        C8_val = C8_val * inEngUnit * (inLenUnit**8)
        self%A(iType, jType) = A_val
        self%A(jType, iType) = A_val
        self%rho(iType, jType) = rho_val
        self%rho(jType, iType) = rho_val
        self%C6(iType, jType) = C6_val
        self%C6(jType, iType) = C6_val
        self%C8(iType, jType) = C8_val
        self%C8(jType, iType) = C8_val
        ! Compute energy shift at cutoff
        call self%ComputeEShift(iType, jType)

      case("rmin")
        ! Format: rmin iType jType value
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) iType
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) jType
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) realVal
        realVal = realVal * inLenUnit
        self%rMinTable(iType, jType) = realVal * realVal  ! Store as r^2
        self%rMinTable(jType, iType) = realVal * realVal

      case("tailcorrection")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) logicVal
        self%usetailcorrection = logicVal

      case default
        ! Unknown command, could add error handling here
    end select

  end subroutine
  !=====================================================================
  subroutine ComputeEShift_EP_TosiFumi_Cut(self, iType, jType)
    implicit none
    class(EP_TosiFumi_Cut), intent(inout) :: self
    integer, intent(in) :: iType, jType

    real(dp) :: r, r6, r8, E_Rep, E_Disp, EShift_val

    r = self%rCut
    r6 = r**6
    r8 = r6* r * r  ! r^8
    E_Rep = self%A(iType, jType) * exp(-r / self%rho(iType, jType))
    E_Disp = -self%C6(iType, jType) / r6
    E_Disp = E_Disp - self%C8(iType, jType) / r8

    EShift_val = E_Rep + E_Disp

    self%EShift(iType, jType) = EShift_val
    self%EShift(jType, iType) = EShift_val

  end subroutine
  !=====================================================================
end module FF_EasyEP_TosiFumi_Cut
!=====================================================================
