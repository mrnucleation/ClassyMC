!====================================================================
module AcceptRuleTemplate
  use MasterTemplate, only: classyClass
  use VarPrecision
  use CoordinateTypes, only: Displacement, Perturbation

  type, public, extends(classyClass) :: acceptrule
    contains
       procedure, pass :: MakeDecision
!       procedure, pass :: Maintenance
       procedure, pass :: ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
!    type(Displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff, inProb
    logical :: accept

    accept = .true.
  end function
!====================================================================
!  subroutine Maintenance(self)
!    implicit none
!    class(acceptrule), intent(inout) :: self
!
!  end subroutine
!====================================================================
  subroutine ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(acceptrule), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    lineStat = 0
  end subroutine
!====================================================================
end module
!====================================================================
