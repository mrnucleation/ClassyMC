!====================================================================
module MinMetroRule
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: MinMetro
    contains
       procedure, pass :: MakeDecision => MinMetro_MakeDecision
!       procedure, pass :: Maintenance => minmetro_Maintenance
!       procedure, pass :: ProcessIO => minmetro_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function MinMetro_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb, extraIn) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(minmetro), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional:: inProb, logProb
    real(dp), intent(in), optional:: extraIn
    real(dp), intent(in) :: E_Diff
    logical :: accept
    real(dp) :: biasE, extraTerms, probTerm


    if(present(inProb)) then
      probTerm = log(inProb)
    elseif(present(logProb)) then
      probTerm = logProb
    else
      write(0,*) "Coding Error! Probability has not been passed into Sampling "
      error stop
    endif

    if(present(extraIn)) then
      extraTerms = extraIn
    else
      extraTerms = 0E0_dp
    endif



    biasE = -trialBox%beta * E_Diff + probTerm + extraTerms
    if(biasE <= 0.0E0_dp) then
      accept = .true.
    else
      accept = .false.
    endif


  end function
!====================================================================
end module
!====================================================================
