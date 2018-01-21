!====================================================================
module AcceptRuleTemplate
  use VarPrecision
  use CoordinateTypes, only: Displacement

  type, public :: acceptrule
    contains
       procedure, pass :: MakeDecision
       procedure, pass :: Maintenance
       procedure, pass :: ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(in) :: self
    class(simBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff, inProb
    logical :: accept

    accept = .true.
  end function
!====================================================================
  subroutine Maintenance(self)
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(in) :: self

  end subroutine
!====================================================================
  subroutine ProcessIO(self, line, lineStat)
    implicit none
    class(acceptrule), intent(in) :: self
    integer, intent(out) :: lineStat
    character(len=*), intent(in) :: line   

  end subroutine
!====================================================================
end module
!====================================================================
