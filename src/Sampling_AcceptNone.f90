!====================================================================
! As it's name suggests, this sampling style rejects every move. Not sure why
! you will need it, but hey it's here if you do.
!====================================================================
module AcceptNoneRule
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: AcceptNone
    contains
       procedure, pass :: MakeDecision => AcceptNone_MakeDecision
       procedure, pass :: MakeDecision2Box => AcceptNone_MakeDecision2Box
!       procedure, pass :: Maintenance => AcceptNone_Maintenance
!       procedure, pass :: ProcessIO => AcceptNone_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function AcceptNone_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb, extraIn) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(AcceptNone), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional:: inProb, logProb, extraIn
    real(dp), intent(in) :: E_Diff
    logical :: accept
    
    accept = .false.
  end function
!====================================================================
  function AcceptNone_MakeDecision2Box(self, trialBox1,  trialBox2, E_Diff1, E_Diff2,&
                                      disp1, disp2, inProb, logProb, extraIn ) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(AcceptNone), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox1, trialBox2
    class(Perturbation), intent(in) :: disp1(:), disp2(:)
    real(dp), intent(in) :: E_Diff1, E_Diff2
    real(dp), intent(in), optional :: inProb, logProb, extraIn
    logical :: accept
    integer :: iDisp
    real(dp) :: biasE, chemPot, extraTerms, probTerm

    accept = .false.
  end function
!====================================================================
end module
!====================================================================
