!====================================================================
module NestedSampling
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: Nested
    integer :: binmiss = 0
    real(dp) :: EMedian
    real(dp) :: EMin, EMax, dE
    real(dp) :: EHist(0:1000)
    character(len=560) :: logfile = "Nested.dat"
    contains
       procedure, pass :: Prologue => Nested_Prologue
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
    real(dp) :: E_Total, E_PerAtom, chemPot, extraTerms, probTerm


    E_Total = trialBox%ETotal + E_Diff
    E_PerAtom = E_Total/trialBox%nAtoms
    if(trialBox%ETotal/trialBox%nAtoms > self%EMedian) then
      if(E_PerAtom < trialBox%ETotal/trialBox%nAtoms) then
        accept = .true.
        return
      endif
    endif


    if(E_PerAtom <= self%EMedian) then
      accept = .true.
      if(E_PerAtom > self%EMin .and. E_Total < self%EMax) then
        ebin = floor( (E_PerAtom-self%EMin)*self%dE )
        self%EHist(ebin) = self%EHist(ebin) + 1E0_dp
      else
        self%binmiss = self%binmiss + 1
      endif
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
    use Units, only: outEngUnit
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
          self%EMax = self%EMax*outEngUnit
          self%dE = 1000.0E0_dp/(self%EMax-self%EMin)

      case("emin")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) self%EMin
          self%EMin = self%EMin*outEngUnit
          self%dE = 1000.0E0_dp/(self%EMax-self%EMin)

      case("filename")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) self%logfile

      case default
        lineStat = -1
    end select

  end subroutine
!====================================================================
  subroutine Nested_Prologue(self)
    implicit none
    class(Nested), intent(inout) :: self

    self%EMedian = self%EMax
    self%EHist = 0E0_dp

  end subroutine
!====================================================================
  subroutine Nested_Maintenance(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    use Units, only: outEngUnit
#ifdef PARALLEL
    use MPI
#endif
    implicit none
    class(Nested), intent(inout) :: self
    integer :: nMedian, iError, arraySize
    real(dp) :: norm
    real(dp) :: sumint
#ifdef PARALLEL
    real(dp) :: temphist(0:1000)
#endif

    norm = sum(self%EHist)*0.5E0_dp
    write(nout,*) "Updating Nested Sampling..."

#ifdef PARALLEL
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 

    arraySize = size(self%EHist)     
    if(myid .eq. 0) then
      TempHist = 0E0_dp
    endif
    call MPI_REDUCE(self%EHist, TempHist, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)       
#endif
    if(myid == 0) then
      sumint = 0.0E0_dp
      nMedian = -1
      do   
        nMedian = nMedian + 1
        sumint = sumint + self%EHist(nMedian)
  !      write(nout,*) nMedian, norm, sumint, self%EHist(nMedian)
        if(sumint > norm) then
          exit
        endif

      enddo
  !    write(nout,*) nMedian, norm, sumint
    endif

    self%EMedian = 0.5E0_dp*( (nMedian/self%dE) + (nMedian-1)/self%dE + 2.0*self%EMin)
    self%EHist = 0E0_dp
    self%EMax = self%EMedian
    self%dE = 1000.0E0_dp/(self%EMax - self%EMin)
    write(nout,*) "New Median Value:", self%EMedian/outEngUnit

  end subroutine
!====================================================================
end module
!====================================================================
