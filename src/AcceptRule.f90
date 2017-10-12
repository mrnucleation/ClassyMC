!====================================================================
module AcceptRuleTemplate
  use VarPrecision
  use CoordinateTypes, only: Displacement

  type, public :: acceptrule
    contains
       procedure, pass :: makedecision
  end type
!====================================================================
  contains
!====================================================================
  function MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use SimBoxDef, only: SimBox
    implicit none
    class(acceptrule), intent(in) :: self
    class(simBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff, inProb
    logical :: accept

    accept = .true.
  end function
!====================================================================
end module
!====================================================================
