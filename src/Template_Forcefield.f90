module Template_ForceField
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(classyClass) :: forcefield
    real(dp) :: rCut, rCutSq
    logical :: isSubPair = .false.  ! True if this is a sub-forcefield in a hybrid
    contains
      procedure, pass :: Constructor 
      procedure, pass :: DetailedECalc 
      procedure, pass :: DiffECalc
      procedure, pass :: SinglePair
      procedure, pass :: SinglePair_Approx
      procedure, pass :: ManyBody
!      procedure, pass :: ShiftECalc_Single
!      procedure, pass :: ShiftECalc_Multi
!      procedure, pass :: NewECalc
!      procedure, pass :: OldECalc
!      procedure, pass :: VolECalc
      procedure, pass :: ProcessIO
      procedure, pass :: GetCutOff
      procedure, pass :: PairForceScale
      procedure, pass :: HasComputeForces
      procedure, pass :: ComputeForces
      procedure, pass :: ComputeForcesFiniteDiff
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    implicit none
    class(forcefield), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!============================================================================
  subroutine DiffECalc(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(inout), target :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp


  end subroutine
!=============================================================================+
  function SinglePair(self, rsq, atmtype1, atmtype2) result(E_Pair)
    implicit none
    class(forcefield), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    E_Pair = 30E0_dp

  end function
!=============================================================================+
  function SinglePair_Approx(self, rsq, atmtype1, atmtype2) result(E_Pair)
    implicit none
    class(forcefield), intent(in) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

  end function
!=============================================================================+
  subroutine ManyBody(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept)
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    integer, intent(in) :: atmtype1
    integer, intent(in) :: atmtypes(:)
    real(dp), intent(in) :: pos1(:)
    real(dp), intent(in) :: posN(:,:)
    logical, intent(out) :: accept
    real(dp), intent(out) :: E_Many

    E_Many = 0E0_dp
    accept = .true.

  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(forcefield), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
  function GetRMin(self) result(rMinOut)
    implicit none
    class(forcefield), intent(inout) :: self
    real(dp), pointer :: rMinOut(:)

    rMinOut => null()

  end function
!=============================================================================+
  function GetCutOff(self) result(rCut)
    implicit none
    class(forcefield), intent(inout) :: self
    real(dp) :: rCut

!    write(*,*) self%rCut
    rCut = self%rCut
  end function
!=============================================================================+
  function HasComputeForces(self) result(hasForce)
    ! True when this forcefield overrides ComputeForces with a real
    ! implementation (analytic pair derivative, etc.). The box uses
    ! finite-difference DiffECalc when this returns false.
    implicit none
    class(forcefield), intent(in) :: self
    logical :: hasForce

    hasForce = .false.
  end function
!=============================================================================+
  function PairForceScale(self, rsq, atmtype1, atmtype2) result(scale)
    ! Returns the scalar such that F_i = scale * (r_i - r_j) for a pair.
    implicit none
    class(forcefield), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: scale

    scale = 0E0_dp
  end function
!=============================================================================+
  subroutine ComputeForces(self, curbox, forces, coords)
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(out) :: forces(:,:)
    real(dp), intent(in) :: coords(:,:)

    forces = 0E0_dp
  end subroutine
!=============================================================================+
  subroutine ComputeForcesFiniteDiff(self, curbox, forces, coords)
    ! Finite-difference fallback: F_i ≈ (U(r_i + δ ê) - U(r)) / δ
    ! evaluated from DiffECalc. This uses the current box coordinates as
    ! the reference state; coords is accepted for interface consistency
    ! with ComputeForces.
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(out) :: forces(:,:)
    real(dp), intent(in) :: coords(:,:)
    logical :: accept
    integer :: iAtom, nMax
    type(Displacement) :: disp(1:1)
    real(dp) :: E_Diff
    integer :: tempList(1,1)
    integer :: tempNNei(1)

    forces = 0E0_dp
    tempList = 0
    tempNNei = 0
    nMax = min(curbox%nMaxAtoms, size(forces, 2), size(coords, 2), size(curbox%atoms, 2))
    do iAtom = 1, nMax
      if(.not. curbox%IsActive(iAtom)) cycle
      disp(1)%molType = curbox%MolType(iAtom)
      disp(1)%molIndx = curbox%MolIndx(iAtom)
      disp(1)%atmIndx = iAtom

      disp(1)%x_new = curbox%atoms(1, iAtom) + curbox%forcedelta
      disp(1)%y_new = curbox%atoms(2, iAtom)
      disp(1)%z_new = curbox%atoms(3, iAtom)
      call self%DiffECalc(curbox, disp(1:1), tempList, tempNNei, E_Diff, accept)
      forces(1, iAtom) = E_Diff / curbox%forcedelta

      disp(1)%x_new = curbox%atoms(1, iAtom)
      disp(1)%y_new = curbox%atoms(2, iAtom) + curbox%forcedelta
      call self%DiffECalc(curbox, disp(1:1), tempList, tempNNei, E_Diff, accept)
      forces(2, iAtom) = E_Diff / curbox%forcedelta

      disp(1)%y_new = curbox%atoms(2, iAtom)
      disp(1)%z_new = curbox%atoms(3, iAtom) + curbox%forcedelta
      call self%DiffECalc(curbox, disp(1:1), tempList, tempNNei, E_Diff, accept)
      forces(3, iAtom) = E_Diff / curbox%forcedelta
    enddo
  end subroutine
!=============================================================================+
end module
!=============================================================================+
