!====================================================================
module MinMetroRule
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use AcceptRuleTemplate, only: acceptrule
 
  type, public, extends(acceptrule) :: MinMetro
    contains
       procedure, pass :: MakeDecision => MinMetro_MakeDecision
!       procedure, pass :: Maintenance => minmetro_Maintenance
!       procedure, pass :: ProcessIO => minmetro_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function MinMetro_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(minmetro), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(in) :: inProb
    real(dp), intent(in) :: E_Diff
    logical :: accept


    if(E_Diff <= 0.0E0_dp) then
      accept = .true.
    else
      accept = .false.
    endif


  end function
!====================================================================
end module
!====================================================================
