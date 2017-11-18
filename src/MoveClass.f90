!=========================================================================
module MoveClassDef
use SimBoxDef, only: SimBox
use VarPrecision

  type, public :: MCMove
    real(dp) :: atmps, accpt
    contains
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
      procedure, pass :: GetAcceptRate
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
    class(SimBox), intent(inout) :: trialBox
  end subroutine
!=========================================================================
  function GetAcceptRate(self) result(rate)
    class(MCMove), intent(in) :: self
    real(dp) :: rate

    if(self%atmps > 0) then
      rate = self%accpt/self%atmps
    else
      rate = 0E0_dp
    endif

    return
  end function

!=========================================================================
end module
!=========================================================================
