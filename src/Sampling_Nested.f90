!====================================================================
module NestedSampling
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: Nested
    real(dp) :: E_Median
    real(dp) :: E_Min, E_Max, dE
    real(dp) :: E_Hist(0:1000) = 0E0_dp
    character(len=560) :: logfile = "Nested.dat"
    contains
       procedure, pass :: MakeDecision => Nested_MakeDecision
       procedure, pass :: MakeDecision2Box => Nested_MakeDecision2Box
       procedure, pass :: Maintenance => Nested_Maintenance
       procedure, pass :: ProcessIO => Nested_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function Nested_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb, extraIn) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(Nested), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional:: inProb, logProb
    real(dp), intent(in), optional:: extraIn
    real(dp), intent(in) :: E_Diff
    logical :: accept
    integer :: ebin
    real(dp) :: E_Total, chemPot, extraTerms, probTerm

    E_Total = trialBox%ETotal + E_Diff
    if(E_Total <= self%E_Median) then
      accept = .true.
      ebin = floor( (E_Total-self%EMin)*self%dE )
      E_Hist(ebin) = E_Hist(ebin) + 1E0_dp
    else
      accept = .false.
    endif

  end function
!====================================================================
  function Nested_MakeDecision2Box(self, trialBox1,  trialBox2, E_Diff1, E_Diff2, &
                                       disp1, disp2, inProb, &
                                       logProb, extraIn ) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(Nested), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox1, trialBox2
    class(Perturbation), intent(in) :: disp1(:), disp2(:)
    real(dp), intent(in) :: E_Diff1, E_Diff2
    real(dp), intent(in), optional :: inProb, logProb, extraIn
    logical :: accept
    integer :: iDisp
    real(dp) :: biasE, chemPot, extraTerms, probTerm

    stop



    accept = .false.
    if(present(inProb)) then
      if(inProb <= 0E0_dp) then
        return
!      write(nout,*) "Probability:", inProb
!      stop "CRITICAL ERROR! Probability passed to the Nested Sampling Function is zero or negative!"
      endif
    endif

    if(present(extraIn)) then
      extraTerms = extraIn
    else
      extraTerms = 0E0_dp
    endif

    if(present(inProb)) then
      probTerm = log(inProb)
    elseif(present(logProb)) then
      probTerm = logProb
    else
      write(0,*) "Coding Error! Probability has not been passed into Sampling "
      stop
    endif

    biasE = -trialBox1%beta*E_Diff1 - trialBox2%beta*E_Diff2 + probTerm + extraTerms
    if(biasE > 0.0E0_dp) then
      accept = .true.
    elseif( biasE > log(grnd()) ) then
      accept = .true.
    endif

  end function
!====================================================================
  subroutine Nested_ProcessIO(self, line, linestat) 
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(Nested), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line   
    integer, intent(out) :: lineStat

    integer :: i, intVal, intVal2
    real(dp) :: realVal
    character(len=30) :: command

    lineStat  = 0
    call GetXCommand(line, command, 3, lineStat)
    if(lineStat < 0) then
      return
    endif
    select case(trim(adjustl(command)))
      case("adjustfreq")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) realVal
          self%maintFreq = floor(realVal)

      case("emax")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) self%EMax
          self%dE = 1000.0E0_dp/(self%EMax-self%EMin)

      case("emin")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) self%EMin
          self%dE = 1000.0E0_dp/(self%EMax-self%EMin)

      case("filename")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) self%logfile

      case default
        lineStat = -1
    end select

  end subroutine
!====================================================================
  subroutine Nested_Maintenance(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(Nested), intent(inout) :: self
    integer :: nMedian
    real(dp) :: norm
    real(dp) :: sumint

    norm = sum(self%EHist)*0.5E0_dp

    sumint = 0.0E0_dp
    nMedian = -1
    do while(sumint < norm)
      nMedian = nMedian + 1
      sumint = sumint + self%EHist(nMedian)
    enddo

    self%E_Median = 0.5E0_dp*(self%EHist(nMedian) + self%EHist(nMedian-1))
    self%EHist = 0E0_dp

  end subroutine
!====================================================================
end module
!====================================================================
