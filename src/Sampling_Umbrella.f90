!====================================================================
module UmbrellaRule
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use AcceptRuleTemplate, only: acceptrule
 
  type, public, extends(acceptrule) :: umbrella
    contains
       procedure, pass :: MakeDecision => Umbrella_MakeDecision
       procedure, pass :: GetBias => Umbrella_GetBias
!       procedure, pass :: Maintenance => Umbrella_Maintenance
!       procedure, pass :: ProcessIO => Umbrella_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function Umbrella_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(metropolis), intent(in) :: self
    class(simBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(in) :: inProb
    real(dp), intent(in) :: E_Diff
    logical :: accept
    real(dp) :: biasE


    accept = .false.
    biasE = -trialBox%beta * E_Diff + log(inProb)
    if(biasE > 0.0E0_dp) then
      accept = .true.
    elseif(biasE > log(grnd())) then
      accept = .true.
    endif


  end function
!====================================================================
end module
!====================================================================
