module AcceptRuleTemplate
  use VarPrecision
  use CoordinateTypes
  type, public :: acceptrule
    contains
       procedure, pass :: makedecision
  end type

  contains

  function MakeDecision(self, E_Diff, disp) result(accept)
    implicit none
    class(constraint), intent(in) :: self
    type(Displacement), intent(in) :: disp
    real(dp), intent(in) :: E_Diff
    logical :: accept

    

  end subroutine


end module
