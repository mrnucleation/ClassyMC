!====================================================================
module UmbrellaRule
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use AcceptRuleTemplate, only: AcceptRule
 
  type, public, extends(AcceptRule) :: Umbrella

    integer :: nBiasVar = 0
    integer, allocatable :: AnalysisIndex(:)

    integer :: umbrellaLimit = 0
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UBinSize(:)
    real(dp), allocatable :: indexCoeff(:)
    real(dp), allocatable :: binMin(:), binMax(:)
    real(dp), allocatable :: binIndx(:)

    character(len=50) :: fileName = "umbrella.dat"

    contains
       procedure, pass :: Constructor => Umbrella_Constructor
       procedure, pass :: MakeDecision => Umbrella_MakeDecision
       procedure, pass :: GetBiasIndex => Umbrella_GetBiasIndex
       procedure, pass :: GetNewBiasIndex => Umbrella_GetNewBiasIndex
       procedure, pass :: ReadInitialBias => Umbrella_ReadInitialBias
!       procedure, pass :: Epilogue => Umbrella_Epilogue
       procedure, pass :: Prologue => Umbrella_Prologue
       procedure, pass :: GetUIndexArray => Umbrella_GetUIndexArray
!       procedure, pass :: Maintenance => Umbrella_Maintenance
       procedure, pass :: ProcessIO => Umbrella_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  subroutine Umbrella_Constructor(self)
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: AllocationStat

    allocate( self%indexCoeff(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%binMax(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%binMin(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%binIndx(1:self%nBiasVar), stat=AllocationStat) 

  end subroutine
!====================================================================
  subroutine Umbrella_Prologue(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: i,j, indx, AllocateStatus

    if(.not. allocated(self%indexCoeff) ) then
      call self % Constructor
    endif

    do i = 1, self%nBiasVar
      indx = self%AnalysisIndex(i)
      AnalysisArray(indx)%func%usedInMove = .true.
    enddo

     ! Since the number of biasing variables is only known at run time, the bias matrix
     ! must be stored in a 1D array.  The index coefficient variables are used to emulate a N dimension matrix
     ! using the linear mapping equation of this form:
     ! U =  a1*x1 + a2*x2 + a3*x3 + .....
     ! Which maps a N dimension matrix onto a 1D array. 

    self%indexCoeff(1) = 1
    do i = 2, self%nBiasVar 
      self%indexCoeff(i) = 1
      do j = 1, i-1
        self%indexCoeff(i) = self%indexCoeff(i) + self%indexCoeff(j) * (self%binMax(j) - self%binMin(j))
      enddo
    enddo      

    self%umbrellaLimit = 1
    do i = 1, self%nBiasVar 
      self%umbrellaLimit = self%umbrellaLimit + self%indexCoeff(i) * (self%binMax(i) - self%binMin(i))
    enddo

    write(nout,*) "Sampling Style: Histogram based Umbrella Sampling"
    write(nout,*) "Number of Umbrella Bins:", self%umbrellaLimit
       
    allocate(self%UBias(1:self%umbrellaLimit+1), STAT = AllocateStatus)
    allocate(self%UHist(1:self%umbrellaLimit+1), STAT = AllocateStatus)

    self%UBias = 0E0_dp
    self%UHist = 0E0_dp


  end subroutine
!====================================================================
  function Umbrella_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use AnalysisData, only: AnalysisArray
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(Umbrella), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(in) :: inProb
    real(dp), intent(in) :: E_Diff

    logical :: accept
    integer :: iBias, oldIndx, newIndx, indx
    real(dp) :: biasE, biasOld, biasNew

    do iBias = 1, self%nBiasVar
      indx = self%AnalysisIndex(iBias)
      call AnalysisArray(indx)%func%CalcNewState(disp)
    enddo

    oldIndx = self%GetBiasIndex()
    call self%GetNewBiasIndex(newIndx, accept)
    if(.not. accept) then
      self%UHist(oldIndx) = self%UHist(oldIndx) + 1E0_dp
      return
    endif
    biasOld = self%UBias(oldIndx)
    biasNew = self%UBias(newIndx)
    

    accept = .false.
    biasE = -trialBox%beta * E_Diff + log(inProb) + (biasNew-biasOld)
    if(biasE > 0.0E0_dp) then
      accept = .true.
    elseif(biasE > log(grnd())) then
      accept = .true.
    endif

    if(accept) then
      self%UHist(newIndx) = self%UHist(newIndx) + 1E0_dp
    else
      self%UHist(oldIndx) = self%UHist(oldIndx) + 1E0_dp
    endif

  end function
!==========================================================================
  function Umbrella_GetBiasIndex(self)  result(biasIndx)
    use AnalysisData, only: AnalysisArray
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: biasIndx

    integer :: analyIndx
    integer :: iBias, bin
    real(dp) :: biasVal
    
   ! Get the variables that are used for the biasing and figure out which histogram bin they
   ! fall into. 
    do iBias = 1, self%nBiasVar
      analyIndx = self%AnalysisIndex(iBias)
      biasVal = AnalysisArray(analyIndx)%func%GetResult()
      self%binIndx(iBias) = nint(biasVal/self%UBinSize(iBias))
    enddo


   ! Using the bin values from each biasing variable, determine which
    biasIndx = 1
    do iBias = 1, self%nBiasVar
!       write(*,*) biasIndx,self%indexCoeff(iBias), binIndx(iBias), binMin(iBias) 
      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
    enddo

!     write(*,*) biasIndx

  end function
!==========================================================================
  subroutine Umbrella_GetNewBiasIndex(self, biasIndx, accept)
    use AnalysisData, only: AnalysisArray, analyCommon
    implicit none
    class(Umbrella), intent(inout) :: self
    logical, intent(out) :: accept
    integer :: biasIndx

    integer :: analyIndx
    integer :: iBias, bin
    real(dp) :: biasVal
    
   ! Get the variables that are used for the biasing and figure out which histogram bin they
   ! fall into. 
    accept = .true.
    do iBias = 1, self%nBiasVar
      analyIndx = self%AnalysisIndex(iBias)
      biasVal = analyCommon(analyIndx)
      self%binIndx(iBias) = nint(biasVal/self%UBinSize(iBias))

      if(self%binIndx(iBias) > self%binMax(iBias) ) then
        accept = .false.
        return
      endif

      if(self%binIndx(iBias) < self%binMin(iBias) ) then
        accept = .false.
        return
      endif
    enddo


   ! Using the bin values from each biasing variable, determine which
    biasIndx = 1
    do iBias = 1, self%nBiasVar
!       write(*,*) biasIndx,self%indexCoeff(iBias), binIndx(iBias), binMin(iBias) 
      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
    enddo

!     write(*,*) biasIndx

  end subroutine

!==========================================================================================
    subroutine Umbrella_ReadInitialBias(self)
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: AllocateStatus
    integer :: j, iInput, iBias, inStat, biasIndx
    real(dp), allocatable :: varValue(:)
    real(dp) :: curBias

    open(unit=36, file=trim(adjustl(self%filename)) )

    allocate(varValue(1:self%nBiasVar), STAT = AllocateStatus )
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    self%UBias = 0E0_dp
    do 
      read(36, *, IOSTAT=inStat) (varValue(j), j=1,self%nBiasVar), curBias

      if(inStat .lt. 0) then
        exit
      endif

      call self%GetUIndexArray(varValue, biasIndx, inStat) 
!      write(*, *) (varValue(j), j=1,self%nBiasVar), curBias
!      write(*,*) varValue, biasIndx
      if(inStat .eq. 1) then
        cycle
      endif

      self%UBias(biasIndx) = curBias
    enddo

    deallocate(varValue)

    close(36)

    end subroutine
!====================================================================
  subroutine Umbrella_GetUIndexArray(self, varArray, biasIndx, stat) 
    implicit none
    class(Umbrella), intent(inout) :: self
    real(dp), intent(in) :: varArray(:)
    integer, intent(out) :: biasIndx, stat
    integer :: iBias
      

    stat = 0
    do iBias = 1, self%nBiasVar
!       binIndx(iBias) = floor( varArray(iBias) / UBinSize(iBias) + 1E-8 )
      self%binIndx(iBias) = nint( varArray(iBias) / self%UBinSize(iBias) )
      if(self%binIndx(iBias) > self%binMax(iBias)) then
        stat = 1
        return
      endif
      if(self%binIndx(iBias) < self%binMin(iBias)) then
        stat = -1
        return
      endif
    enddo

    biasIndx = 1
    do iBias = 1, self%nBiasVar
      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
    enddo

   end subroutine
!====================================================================
  subroutine Umbrella_ProcessIO(self, line, linestat) 
    use Input_Format, only: GetXCommand
    implicit none
    class(Umbrella), intent(inout) :: self
    character(len=*), intent(in) :: line   
    integer, intent(out) :: lineStat

    integer :: intVal, intVal2
    character(len=30) :: command

    lineStat  = 0
    call GetXCommand(line, command, 3, lineStat)
    if(lineStat < 0) then
      return
    endif
    select case(trim(adjustl(command)))
      case("biasvariables")
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) intVal
        self%nBiasVar = intVal
        allocate(self%AnalysisIndex(1:intVal) ) 

      case("analysis")
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) intVal

        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal2
 
        self%AnalysisIndex(intVal) = intVal2

      case default
        lineStat = -1
    end select



   end subroutine
!====================================================================
end module
!====================================================================
