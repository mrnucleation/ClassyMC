!====================================================================
module Template_AcceptRule
  use MasterTemplate, only: classyClass
  use VarPrecision
  use CoordinateTypes

  type, public, extends(classyClass) :: acceptrule
    contains
       procedure, pass :: MakeDecision
       procedure, pass :: MakeDecision2Box
       procedure, pass :: UpdateStatistics 
       procedure, pass :: GetExtraTerms
!       procedure, pass :: Maintenance
       procedure, pass :: ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb, extraIn ) result(accept)
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff
    real(dp), intent(in), optional :: inProb, logProb, extraIn
    logical :: accept

    accept = .true.
  end function
!====================================================================
  function MakeDecision2Box(self, trialBox1,  trialBox2, E_Diff1, E_Diff2, &
                           disp1, disp2, inProb, logProb, extraIn ) result(accept)
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(inout) :: self
    class(simBox), intent(in) :: trialBox1, trialBox2
    class(Perturbation), intent(in) :: disp1(:), disp2(:)
    real(dp), intent(in) :: E_Diff1, E_Diff2
    real(dp), intent(in), optional :: inProb, logProb, extraIn
    logical :: accept

    stop "This Sampling Procedure does not have a 2box acceptance rule defined"
    accept = .true.
  end function
!====================================================================
  subroutine UpdateStatistics(self, accept)
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(inout) :: self
    logical, intent(in) :: accept

  end subroutine
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
  function GetExtraTerms(self, disp, trialBox) result(extraTerms)
    use Input_Format, only: maxLineLen
!    use SimpleSimBox, only: SimpleBox
    use Template_SimBox, only: SimBox
    implicit none
    class(acceptrule), intent(in) :: self
    class(Perturbation), intent(in) :: disp(:)
    class(SimBox), intent(in) :: trialBox
    real(dp) :: extraTerms

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
            extraTerms = extraTerms - (disp(1)%volNew-disp(1)%volOld) * trialBox%pressure*trialBox%beta
       end select

  end function
!====================================================================
end module
!====================================================================
