!====================================================================
module DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint

  type, public, extends(constraint) :: distcriteria
    integer, allocatable :: neiListIndx(:)
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
  subroutine DistCrit_ShiftCheck(self, trialBox, disp, accept)
    use CoordinateTypes, only: Displacement
    use SimBoxDef, only: SimBox
    implicit none
    class(distcriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    
    integer :: indx2
    real(dp) :: rx, ry, rz, rsq
    

    accept = .true.
!    write(*,*) trialBox%atoms(1,1), trialBox%atoms(1,2)
    if(disp(1)%atmindx == 1) then
      indx2 = 2 
    else
      indx2 = 1
    endif
    rx = disp(1)%x_new - trialBox%atoms(1, indx2)
    ry = disp(1)%y_new - trialBox%atoms(2, indx2)
    rz = disp(1)%z_new - trialBox%atoms(3, indx2)
    rsq = rx*rx + ry*ry + rz*rz
    if(rsq > 2.0E0_dp**2) then
      accept = .false.
    endif


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
