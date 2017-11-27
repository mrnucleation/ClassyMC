!=============================================================
module Constrain_HardWall
  use ConstraintTemplate, only: constraint
  use VarPrecision

  type, public, extends(constraint) :: hardwall
    integer, allocatable :: atmTypes(:)
    real(dp) :: xhi, xlo
    real(dp) :: yhi, ylo
    real(dp) :: zhi, zlo
    contains
      procedure, pass :: CheckInitialConstraint => HardWall_CheckInitialConstraint
      procedure, pass :: ShiftCheck => HardWall_ShiftCheck
      procedure, pass :: SwapInCheck => HardWall_SwapInCheck
      procedure, pass :: SwapOutCheck => HardWall_SwapOutCheck
  end type

!=============================================================
  contains
!=============================================================
  subroutine HardWall_CheckInitialConstraint(self)
    implicit none
    class(hardwall), intent(in) :: self
  end subroutine
!=============================================================
  subroutine HardWall_ShiftCheck(self, trialBox, disp, accept)
    use CoordinateTypes, only: Displacement
    use Template_SimBox, only: SimBox
    implicit none
    class(hardwall), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.

  end subroutine
!=============================================================
  subroutine HardWall_SwapInCheck(self)
    implicit none
    class(hardwall), intent(in) :: self
  end subroutine
!=============================================================
  subroutine HardWall_SwapOutCheck(self)
    implicit none
    class(hardwall), intent(in) :: self
  end subroutine
!=============================================================
end module
!=============================================================
