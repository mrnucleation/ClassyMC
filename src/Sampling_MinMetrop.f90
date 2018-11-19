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
  function MinMetro_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(minmetro), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional:: inProb, logProb
    real(dp), intent(in) :: E_Diff
    logical :: accept
    real(dp) :: biasE, extraTerms, probTerm

    extraTerms = 0E0_dp
    select type(disp)
      class is(Addition)
          extraTerms = extraTerms + trialBox%chempot(disp(1)%molType)
      class is(Deletion)
          extraTerms = extraTerms - trialBox%chempot(disp(1)%molType)
      class is(VolChange)
          extraTerms = extraTerms + (disp(1)%volNew -disp(1)%volOld)*trialBox%pressure*trialBox%beta
    end select

    if(present(inProb)) then
      probTerm = log(inProb)
    elseif(present(logProb)) then
      probTerm = logProb
    else
      write(*,*) "Coding Error! Probability has not been passed into Sampling "
      stop
    endif



    biasE = -trialBox%beta * E_Diff + probTerm + extraTerms
    if(biasE <= 0.0E0_dp) then
      accept = .true.
      write(*,*) biasE
    else
      accept = .false.
    endif


  end function
!====================================================================
end module
!====================================================================
