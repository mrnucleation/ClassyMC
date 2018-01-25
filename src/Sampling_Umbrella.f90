!====================================================================
module UmbrellaRule
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use AcceptRuleTemplate, only: acceptrule
 
  type, public, extends(acceptrule) :: Umbrella
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UHistTotal(:)
    real(dp), allocatable :: UBinSize(:)

    integer :: nBiasVar = 0
    integer, allocatable :: AnalysisIndex(:)

    contains
       procedure, pass :: MakeDecision => Umbrella_MakeDecision
       procedure, pass :: GetBias => Umbrella_GetBias
!       procedure, pass :: Prologue => Umbrella_Prologue
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
    class(Umbrella), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
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
