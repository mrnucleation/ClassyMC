!====================================================================
module DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Displacement
  use Template_SimBox, only: SimBox

  type, public, extends(constraint) :: distcriteria
    integer :: neiList = -1
    integer :: molType = -1
    contains
      procedure, pass :: CheckInitialConstraint => DistCrit_CheckInitialConstraint
!      procedure, pass :: DiffCheck
      procedure, pass :: ShiftCheck => DistCrit_ShiftCheck
      procedure, pass :: NewCheck => DistCrit_NewCheck
      procedure, pass :: OldCheck => DistCrit_OldCheck
!      procedure, pass :: ProcessIO => DistCrit_ProcessIO
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine DistCrit_CheckInitialConstraint(self, accept)
    implicit none
    class(distcriteria), intent(in) :: self
    logical, intent(out) :: accept
  end subroutine
!=====================================================================
  subroutine DistCrit_ShiftCheck(self, trialBox, disp, accept)
    implicit none
    class(distcriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    
    integer :: indx2
    real(dp) :: rx, ry, rz, rsq
    

    accept = .true.

  end subroutine
!=====================================================================
  subroutine DistCrit_NewCheck(self, trialBox, disp, accept)
    implicit none
    class(distcriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=====================================================================
  subroutine DistCrit_OldCheck(self, trialBox, disp, accept)
    implicit none
    class(distcriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=====================================================================
end module
!=====================================================================
