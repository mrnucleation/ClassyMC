!=============================================================
module ConstraintTemplate
  use VarPrecision

  type, public :: constraint
    contains
      procedure, pass :: Constructor
      procedure, pass :: CheckInitialConstraint
      procedure, pass :: ShiftCheck
      procedure, pass :: SwapInCheck
      procedure, pass :: SwapOutCheck
      procedure, pass :: ProcessIO
  end type

  type, public :: constrainArray
    class(constraint), allocatable :: method
  end type
!=============================================================
  contains
!=============================================================
  subroutine Constructor(self)
    implicit none
    class(constraint), intent(inout) :: self
  end subroutine
!=============================================================
  subroutine CheckInitialConstraint(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine
!=============================================================
  subroutine ShiftCheck(self, trialBox, disp, accept)
    use CoordinateTypes, only: Displacement
    use Template_SimBox, only: SimBox
    implicit none
    class(constraint), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.

  end subroutine
!=============================================================
  subroutine SwapInCheck(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine
!=============================================================
  subroutine SwapOutCheck(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine
!=============================================================
  subroutine ProcessIO(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine
!=============================================================
end module
!=============================================================
