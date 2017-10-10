!====================================================================
module DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint

  type, public, extends(constraint) :: distcriteria
    integer :: neiList = -1
    contains
      procedure, pass :: CheckInitialConstraint => DistCrit_CheckInitialConstraint
      procedure, pass :: ShiftCheck => DistCrit_ShiftCheck
      procedure, pass :: SwapInCheck => DistCrit_SwapInCheck
      procedure, pass :: SwapOutCheck => DistCrit_SwapOutCheck
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine DistCrit_CheckInitialConstraint(self)
    implicit none
    class(distcriteria), intent(in) :: self
  end subroutine
!=====================================================================
  subroutine DistCrit_ShiftCheck(self, disp, rejMove)
    use CoordinateTypes, only: Displacement
    implicit none
    class(distcriteria), intent(in) :: self
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: rejMove

    write(*,*) "Here in my Dispy!"
    rejMove = .false.
  end subroutine
!=====================================================================
  subroutine DistCrit_SwapInCheck(self)
    implicit none
    class(distcriteria), intent(in) :: self
  end subroutine
!=====================================================================
  subroutine DistCrit_SwapOutCheck(self)
    implicit none
    class(distcriteria), intent(in) :: self
  end subroutine
!=====================================================================
end module
!=====================================================================
