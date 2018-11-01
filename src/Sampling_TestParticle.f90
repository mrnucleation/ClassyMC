!====================================================================
! Sampling method for use in Widom Test Particle simulations.
!====================================================================
module TestParticleRule
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: TestParticle
    logical :: useDiff = .false.
    real(dp) :: E_Diff = 0E0_dp
    real(dp) :: EnsembAvg = 0E0_dp
    real(dp) :: cnt = 0E0_dp
    contains
       procedure, pass :: MakeDecision => TestParticle_MakeDecision
       procedure, pass :: UpdateStatistics => TestParticle_UpdateStatistics
!       procedure, pass :: Maintenance => TestParticle_Maintenance
!       procedure, pass :: ProcessIO => TestParticle_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function TestParticle_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(TestParticle), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
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
!      stop "CRITICAL ERROR! Probability passed to the TestParticle Sampling Function is zero or negative!"
    endif


     ! The purpose of this section is to add any terms such as the isobaric or
     ! grand canonical ensemble terms (IE the PV or chemical potential) to the
     ! detailed balance condition. 
     ! Currently only coded for single molecule moves, need to generalize this
     ! in the future.
    extraTerms = 0E0_dp
    select type(disp)
      class is(Addition)
         self%useDiff = .true.
         self%E_Diff = E_Diff
         accept = .false.
         return
      class is(Deletion)
         stop "TestParticle Sampling Can't handle deletion just yet"
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
  subroutine TestParticle_UpdateStatistics(self, accept)
    implicit none
    class(TestParticle), intent(inout) :: self
    logical, intent(in) :: accept

    if(self%useDiff) then
       self%useDiff = .false.
    else
    endif


  end subroutine

!====================================================================
end module
!====================================================================
