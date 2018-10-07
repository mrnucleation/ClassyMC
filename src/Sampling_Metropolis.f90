!====================================================================
module MetropolisRule
  use VarPrecision
  use CoordinateTypes, only: Displacement, Perturbation, Addition, Deletion, VolChange
  use AcceptRuleTemplate, only: acceptrule
 
  type, public, extends(acceptrule) :: metropolis
    contains
       procedure, pass :: MakeDecision => Metropolis_MakeDecision
!       procedure, pass :: Maintenance => Metropolis_Maintenance
!       procedure, pass :: ProcessIO => Metropolis_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function Metropolis_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(metropolis), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
!    type(Displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: inProb
    real(dp), intent(in) :: E_Diff
    logical :: accept
    integer :: iDisp
    real(dp) :: biasE, extraTerms, chemPot

    accept = .false.
    if(inProb <= 0E0_dp) then
      return
!      write(nout,*) "Probability:", inProb
!      stop "CRITICAL ERROR! Probability passed to the Metropolis Sampling Function is zero or negative!"
    endif


     ! The purpose of this section is to add any terms such as the isobaric or
     ! grand canonical ensemble terms (IE the PV or chemical potential) to the
     ! detailed balance condition. 
     ! Currently only coded for single molecule moves, need to generalize this
     ! in the future.
    extraTerms = 0E0_dp
    select type(disp)
      class is(Addition)
          extraTerms = extraTerms + trialBox%chempot(disp(1)%molType)
      class is(Deletion)
          extraTerms = extraTerms - trialBox%chempot(disp(1)%molType)
      class is(VolChange)
          extraTerms = extraTerms + (disp(1)%volNew -disp(1)%volOld)*trialBox%pressure*trialBox%beta
    end select

    biasE = -trialBox%beta * E_Diff + log(inProb) + extraTerms
    if(biasE > 0.0E0_dp) then
      accept = .true.
    elseif( biasE > log(grnd()) ) then
      accept = .true.
    endif

  end function
!====================================================================
end module
!====================================================================
