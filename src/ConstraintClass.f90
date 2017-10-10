!=============================================================
module ConstraintTemplate
  use VarPrecision
  type, public :: constraint
    contains
      procedure, pass :: CheckInitialConstraint 
      procedure, pass :: ShiftCheck
      procedure, pass :: SwapInCheck
      procedure, pass :: SwapOutCheck
  end type
!=============================================================
  contains
!=============================================================
  subroutine CheckInitialConstraint(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine
!=============================================================
  subroutine ShiftCheck(self, disp, rejMove)
    use CoordinateTypes, only: Displacement
    implicit none
    class(constraint), intent(in) :: self
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: rejMove

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
end module
!=============================================================
