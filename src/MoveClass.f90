module MoveClassDef
use VarPrecision

  type, extensible :: MCMove
    real(dp) :: atmps, accpt
    contains
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
      procedure, nopass, non_overridable  :: GetAcceptRate
  end type

 contains

  subroutine GeneratePosition(self)
    class(MCMove), intent(in) :: self
  end subroutine

  subroutine FullMove(self, trialBox)
    class(MCMove), intent(in) :: self
    class(SimBox), intent(inout) :: trialBox
  end subroutine

  function GetAcceptRate(self) result(rate)
    class(MCMove), intent(in) :: self
    real(dp) :: rate
    if(self%atmps > 0) then
      rate = self%accpt/self%atmpts
    else
      rate = 0E0_dp
    endif

    return
  end subroutine


end module
