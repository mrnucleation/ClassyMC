!=========================================================================
module MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public :: MCMove
    real(dp) :: atmps = 1E-30_dp
    real(dp) :: accpt = 0E0_dp

    contains
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance 
  end type

 contains
!=========================================================================
  subroutine GeneratePosition(self, disp)
    use CoordinateTypes, only: Displacement
    implicit none
    class(MCMove), intent(in) :: self
    type(Displacement), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine FullMove(self, trialBox)
    class(MCMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
  end subroutine
!=========================================================================
  function GetAcceptRate(self) result(rate)
    class(MCMove), intent(in) :: self
    real(dp) :: rate

    if(self%atmps > 0E0_dp) then
      rate = 1E2_dp*self%accpt/self%atmps
    else
      rate = 0E0_dp
    endif

    return
  end function
!=========================================================================
  subroutine Maintenance(self)
    class(MCMove), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
