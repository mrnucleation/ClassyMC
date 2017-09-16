module DistanceCriteria
  use VarPrecision
  type, public, extends(constraint) :: distcriteria
    contains
      procedure, pass :: CheckInitialConstraint 
      procedure, pass :: ShiftCheck
      procedure, pass :: SwapInCheck
      procedure, pass :: SwapOutCheck
  end type

  contains

  subroutine CheckInitialConstraint(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine

  subroutine ShiftCheck(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine

  subroutine SwapInCheck(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine

  subroutine SwapOutCheck(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine



end module
