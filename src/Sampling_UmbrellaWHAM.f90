!====================================================================
module UmbrellaWHAMRule
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use AcceptRuleTemplate, only: AcceptRule
 
  type, public, extends(AcceptRule) :: UmbrellaWHAM

    integer :: nBiasVar = 0
    integer, allocatable :: AnalysisIndex(:)

    integer :: umbrellaLimit = 0
    integer, allocatable :: indexCoeff(:)
    integer, allocatable :: binMin(:), binMax(:)
    integer, allocatable :: binIndx(:)
    integer, allocatable :: UArray(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UBinSize(:)
    real(dp), allocatable :: varValues(:)
    character(len=50) :: fileName = "umbrella.dat"

    real(dp), allocatable :: refVals(:)
    integer :: refBin
    integer :: maxSelfConsist = 3000
    real(dp) :: tolLimit = 1E-5_dp
!    integer :: equilInterval
!    integer :: whamEstInterval

    integer :: nWhamItter, nCurWhamItter

    real(dp), allocatable :: WHAM_Numerator(:)
    real(dp), allocatable :: WHAM_Denominator(:,:)
    real(dp), allocatable :: HistStorage(:)
    real(dp), allocatable :: FreeEnergyEst(:)
    real(dp), allocatable :: BiasStorage(:,:)
    real(dp), allocatable :: NewBias(:)
    real(dp), allocatable :: ProbArray(:)      
    real(dp), allocatable :: TempHist(:)
 

    contains
       procedure, pass :: Constructor => UmbrellaWHAM_Constructor
       procedure, pass :: MakeDecision => UmbrellaWHAM_MakeDecision
       procedure, pass :: GetBiasIndex => UmbrellaWHAM_GetBiasIndex
       procedure, pass :: GetNewBiasIndex => UmbrellaWHAM_GetNewBiasIndex
       procedure, pass :: ReadInitialBias => UmbrellaWHAM_ReadInitialBias
       procedure, pass :: GetUIndexArray => UmbrellaWHAM_GetUIndexArray
       procedure, pass :: OutputUmbrellaHist => UmbrellaWHAM_OutputUmbrellaHist
       procedure, pass :: FindVarValues => UmbrellaWHAM_FindVarValues
       procedure, pass :: AdjustHist => UmbrellaWHAM_AdjustHist
       procedure, pass :: WHAMSimOutput => UmbrellaWHAM_WHAMSimOutput
       procedure, pass :: Maintenance => UmbrellaWHAM_Maintenance

       procedure, pass :: ProcessIO => UmbrellaWHAM_ProcessIO
       procedure, pass :: Epilogue => UmbrellaWHAM_Epilogue
       procedure, pass :: Prologue => UmbrellaWHAM_Prologue

  end type
!====================================================================
  contains
!====================================================================
  subroutine UmbrellaWHAM_Constructor(self)
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: AllocationStat

    allocate( self%refVals(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%AnalysisIndex(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%indexCoeff(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%binMax(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%binMin(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%binIndx(1:self%nBiasVar), stat=AllocationStat) 

    allocate( self%UBinSize(1:self%nBiasVar), STAT =AllocationStat )
    allocate( self%UArray(1:self%nBiasVar), STAT =AllocationStat )
    allocate( self%varValues(1:self%nBiasVar), STAT =AllocationStat )

  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_Prologue(self)
    use AnalysisData, only: AnalysisArray
    use SimControl, only: nCycles
    use ParallelVar, only: nout, myid
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: i,j, indx, AllocateStatus

    do i = 1, self%nBiasVar
      indx = self%AnalysisIndex(i)
      if((indx > size(AnalysisArray)) .or. (indx < 1) ) then
        write(nout, *) "ERROR! The Umbrella Sampling routine has been directed to an "
        write(nout, *) "invalid Analysis fucntion"
        write(nout, *) "Chosen Function:", indx
        write(nout, *) "Number of Analysis Functions:", size(AnalysisArray)
        stop "Error dectected in Umbrella Sampling"
      endif
    enddo

    do i = 1, self%nBiasVar
      if(self%binMin(i) > self%binMax(i) ) then
        write(nout,*) "ERROR! The given bounds for one of the umbrella variables does not make sense!"
        write(nout,*) "Smallest bin is larger than the largest bin" 
        write(nout,*) "Minimum Bin:", self%binMin(i)
        write(nout,*) "Maximum Bin:", self%binMax(i)
        stop
      endif
    enddo

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
!    write(*,*) self%indexCoeff
!    write(*,*) self%binMin
!    write(*,*) self%binMax
!    write(*,*) self%UBinSize

    self%umbrellaLimit = 1
    do i = 1, self%nBiasVar 
      self%umbrellaLimit = self%umbrellaLimit + self%indexCoeff(i) * (self%binMax(i) - self%binMin(i))
    enddo

    write(nout,*) "Sampling Style: Histogram based Umbrella Sampling"
    write(nout,*) "Number of Umbrella Bins:", self%umbrellaLimit
       
    allocate(self%UBias(1:self%umbrellaLimit+1), STAT = AllocateStatus)
    allocate(self%UHist(1:self%umbrellaLimit+1), STAT = AllocateStatus)

!    write(nout,*) self%binMax
    self%UBias = 0E0_dp
    self%UHist = 0E0_dp

    call self%ReadInitialBias
    self%nWhamItter = ceiling(dble(nCycles)/dble(self%maintFreq))
    ! Allocation of the WHAM variables
    if(myid .eq. 0) then
      allocate(self%WHAM_Numerator(1:self%umbrellaLimit), STAT = AllocateStatus)
      allocate(self%WHAM_Denominator(1:self%umbrellaLimit,1:self%nWhamItter+1), STAT = AllocateStatus)
      allocate(self%HistStorage(1:self%umbrellaLimit), STAT = AllocateStatus)
      allocate(self%BiasStorage(1:self%umbrellaLimit,1:self%nWhamItter+1), STAT = AllocateStatus)
      allocate(self%FreeEnergyEst(1:self%umbrellaLimit), STAT = AllocateStatus)
      allocate(self%ProbArray(1:self%umbrellaLimit), STAT = AllocateStatus)

      self%WHAM_Numerator = 0E0
      self%WHAM_Denominator = 0E0
      self%HistStorage = 0E0
      self%BiasStorage = 0E0
      self%nCurWhamItter = 1
!        tolLimit = 1d-2
      open(unit = 96, file="WHAM_.incomp")
      open(unit = 97, file="WHAM_Potential.incomp")
      open(unit = 98, file="WHAM_Mid_DG.incomp")
    endif

    allocate(self%TempHist(1:self%umbrellaLimit), STAT = AllocateStatus)      
    allocate(self%NewBias(1:self%umbrellaLimit), STAT = AllocateStatus)
    self%NewBias = 0E0
    self%TempHist = 0E0


  end subroutine
!====================================================================
  function UmbrellaWHAM_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use AnalysisData, only: AnalysisArray
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
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

!    write(*,*) oldIndx, self%UBias(oldIndx)
!    write(*,*) newIndx
    if(.not. accept) then
      self%UHist(oldIndx) = self%UHist(oldIndx) + 1E0_dp
      return
    endif

!    write(*,*) oldIndx, self%UBias(oldIndx)
!    write(*,*) newIndx, self%UBias(newIndx)
!    write(*,*) 
    biasOld = self%UBias(oldIndx)
    biasNew = self%UBias(newIndx)
 


    accept = .false.
    biasE = -trialBox%beta * E_Diff + log(inProb) + (biasNew-biasOld)
!    write(2,*) E_Diff, log(inProb),  (biasNew-biasOld)
!    write(2,*) biasE
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
  function UmbrellaWHAM_GetBiasIndex(self)  result(biasIndx)
    use AnalysisData, only: AnalysisArray
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
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
!      write(*,*) biasVal, self%binIndx(iBias), self%UBinSize(iBias)
    enddo


   ! Using the bin values from each biasing variable, determine which
    biasIndx = 1
    do iBias = 1, self%nBiasVar
      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
!      write(*,*) biasIndx, self%indexCoeff(iBias), self%binIndx(iBias), self%binMin(iBias) 
    enddo

!    write(*,*) 

  end function
!==========================================================================
  subroutine UmbrellaWHAM_GetNewBiasIndex(self, biasIndx, accept)
    use AnalysisData, only: AnalysisArray, analyCommon
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    logical, intent(out) :: accept
    integer, intent(out) :: biasIndx

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

!      write(*,*) biasVal, self%binIndx(iBias)
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
!       write(*,*) biasIndx, self%indexCoeff(iBias), self%binIndx(iBias), self%binMin(iBias) 
      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
    enddo

!    write(*,*) biasIndx

  end subroutine

!==========================================================================================
    subroutine UmbrellaWHAM_ReadInitialBias(self)
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
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
  subroutine UmbrellaWHAM_GetUIndexArray(self, varArray, biasIndx, stat) 
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
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
!==========================================================================================
   subroutine UmbrellaWHAM_OutputUmbrellaHist(self)
     implicit none
     class(UmbrellaWHAM), intent(inout) :: self
     integer :: iUmbrella, iBias, iBin
     character(len = 100) :: outputString

!     write(*,*) "Output"
     write(outputString, *) "(", ("F12.5, 2x", iBias =1,self%nBiasVar), "2x, F18.1)"
     open(unit=60, file="UmbrellaHist.txt")
      
     do iUmbrella = 1, self%umbrellaLimit
       call self%findVarValues(iUmbrella, self%UArray)  
       do iBias = 1, self%nBiasVar        
         self%varValues(iBias) = real(self% UArray(iBias), dp) * self%UBinSize(iBias)          
       enddo

       if(self%UHist(iUmbrella) .ne. 0E0_dp) then
         write(60,outputString) (self%varValues(iBias), iBias=1,self%nBiasVar), self%UHist(iUmbrella)
       endif
        
     enddo 
    
     flush(60)
     close(60)
      
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_FindVarValues(self, UIndx, UArray) 
    implicit none 
    class(UmbrellaWHAM), intent(inout) :: self
    integer, intent(in) :: UIndx 
    integer, intent(inout) :: UArray(:) 
    integer :: i, iBias 
    integer :: remainder, curVal 
                            
    remainder = UIndx - 1 
    do i = 1, self%nBiasVar 
      iBias = self%nBiasVar - i + 1 
      curVal = int( real(remainder, dp)/real(self%indexCoeff(iBias),dp) ) 
      self%UArray(iBias) = curVal + self%binMin(iBias) 
      remainder = remainder - curVal * self%indexCoeff(iBias) 
    enddo 
                                                                             
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_ProcessIO(self, line, linestat) 
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line   
    integer, intent(out) :: lineStat

    integer :: i, intVal, intVal2
    real(dp) :: realVal
    character(len=30) :: command

    lineStat  = 0
    call GetXCommand(line, command, 3, lineStat)
!    write(*,*) " Umbrella", command
    if(lineStat < 0) then
      return
    endif
    select case(trim(adjustl(command)))
      case("biasvariables")
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) intVal
        self%nBiasVar = intVal
        if(.not. allocated(self%AnalysisIndex) ) then
          call self%Constructor
        endif

      case("analysis")
        if(allocated(self%AnalysisIndex) ) then
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) intVal

          call GetXCommand(line, command, 5, lineStat)
          read(command, *) intVal2
 
          self%AnalysisIndex(intVal) = intVal2
        endif

      case("bounds")
        if(allocated(self%AnalysisIndex) ) then
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) intVal

          call GetXCommand(line, command, 5, lineStat)
          read(command, *) realVal
          self%UBinSize(intVal) = realVal

          call GetXCommand(line, command, 6, lineStat)
          read(command, *) realVal
          self%binMin(intVal) = nint(realVal/self%UBinSize(intVal))


          call GetXCommand(line, command, 7, lineStat)
          read(command, *) realVal
          self%binMax(intVal) = nint(realVal/self%UBinSize(intVal))
        endif

      case("reference")
        do i = 1, self%nBiasVar
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) realVal
          self%refVals(i) = realVal
        enddo
        call self%GetUIndexArray(self%refVals, intVal, lineStat) 
        self%refBin = intVal

      case("whamfreq")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) realVal
          self%maintFreq = nint(realVal)

      case default
        lineStat = -1
    end select



   end subroutine
!====================================================================
  subroutine UmbrellaWHAM_Maintenance(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self


    call self%Adjusthist
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_Epilogue(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self

    write(nout,*) "Writing Umbrella Sampling Histogram..."
  end subroutine
!==================================================================================
  subroutine UmbrellaWHAM_WHAMSimOutput(self)
    use ParallelVar
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self

    integer :: i,j
    integer :: UArray(1:self%nBiasVar)
    character(len = 100) :: outputString

    write(outputString, *) "(", ("2x, F10.8,", j =1,self%nBiasVar), "2x, F18.8)"

!        This block exports the calculated free energy to a file
    rewind(96)
    do i = 1, self%umbrellaLimit
      if(self%HistStorage(i) .ne. 0E0 ) then
        call self%findVarValues(i, self%UArray)
        write(96, outputString) (self%UArray(j)*self%UBinSize(j),j=1,self%nBiasVar), self%HistStorage(i)
      endif
    enddo
    flush(96)

    write(outputString, *) "(", ("2x, F10.8,", j =1,self%nBiasVar), "2x, F18.8)"

!        This block exports the current umbrella bias
    rewind(97)
    do i = 1, self%umbrellaLimit
      call self%findVarValues(i, self%UArray)
      write(97, outputString) (self%UArray(j)*self%UBinSize(j),j=1,self%nBiasVar), self%UBias(i)
    enddo
    flush(97)


!        This block exports the calculated free energy to a file
    rewind(98)
    do i = 1, self%umbrellaLimit
      if(self%ProbArray(i) .ne. 0E0 ) then
        call self%findVarValues(i, self%UArray)
        write(98, outputString) (self%UArray(j)*self%UBinSize(j),j=1,self%nBiasVar), self%FreeEnergyEst(i)
      endif
    enddo
    flush(98)



           
  end subroutine
!=========================================================================
!     This subroutine periodically adjusts the Umbrella Sampling Bias
!     by collecting histogram data from across 
  subroutine UmbrellaWHAM_AdjustHist(self)
    use ParallelVar, only: myid, nout, ierror
    use MPI
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: arraySize, i, j, cnt, maxbin, maxbin2
    real(dp) :: norm, maxBias, denomSum
    real(dp) :: F_Estimate(1:self%nWhamItter), F_Old(1:self%nWhamItter), fSum     
    real(dp) :: tol, refBias

    write(nout,*) "Halting for WHAM"
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
    
!      This block condences the histogram data from all the different processors
!      into one collective array on the root (myid = 0) processor.        
    arraySize = size(self%UHist)     
    if(myid .eq. 0) then
      self%TempHist = 0E0
    endif
    call MPI_REDUCE(self%UHist, self%TempHist, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)       

!      write(*,*) "UHIST!"
    if(myid .eq. 0) then
!        This block calculates the terms needed to 
      norm = sum(self%TempHist)
      do i = 1, self%umbrellaLimit
        self%BiasStorage(i, self%nCurWhamItter) = self%UBias(i) 
        if(self%TempHist(i) .ne. 0E0) then
          self%WHAM_Numerator(i) = self%WHAM_Numerator(i) + self%TempHist(i)*self%TempHist(i)/norm
          self%WHAM_Denominator(i, self%nCurWhamItter) = self%TempHist(i)*exp(self%UBias(i))
!            WHAM_Numerator(i) = WHAM_Numerator(i) + (i)
!            WHAM_Denominator(i, nCurWhamItter) = norm*exp(UBias(i))
          self%HistStorage(i) = self%HistStorage(i) + self%TempHist(i)
        endif
      enddo 

!        This block solves for the free energy terms required by WHAM.  This is done
!        self-consistently.
      self%ProbArray = 0E0
      F_Estimate = 0E0
      tol = self%tolLimit + 1E0
      cnt = 0
      do while(tol > self%tolLimit)
        cnt = cnt + 1
!          Infinite Loop protection
        if(cnt > self%maxSelfConsist) then 
          write(35,*) "Self Consistent Limit Hit"
          exit
        endif

        self%ProbArray = 0E0
        do j = 1, self%nCurWhamItter
          F_Old(j) = F_Estimate(j)
        enddo
!            If bin #i has been sampled at any point in the simulation, estimate the unbiased probability
!            based on the current guess value for F
        do i = 1, self%umbrellaLimit
          if(self%WHAM_Numerator(i) /= 0E0) then
            denomSum = 0E0
!              maxBias = minval(F_Estimate)
            do j = 1, self%nCurWhamItter
              if(self%WHAM_Denominator(i,j) > 0E0) then
                denomSum = denomSum + self%WHAM_Denominator(i,j)*exp(-F_Estimate(j))
              endif
            enddo
            if(denomSum /= 0E0) then
              self%ProbArray(i) = self%WHAM_Numerator(i)/denomSum
            endif
          else
            self%ProbArray(i) = 0E0
          endif
        enddo 

        norm = sum(self%ProbArray)
        do i = 1, self%umbrellaLimit
          self%ProbArray(i) = self%ProbArray(i)/norm
        enddo 
!          Once all the unbiased probabilities have been estimated, use these unbiased estimates
!          to calculate a new estimate for F
        do j = 1, self%nCurWhamItter
          fSum = 0E0
!            maxBias = maxval(BiasStorage(:,j))
          do i = 1, self%umbrellaLimit
            if(self%ProbArray(i) .ne. 0E0) then
              fSum = fSum + self%ProbArray(i)*exp(self%BiasStorage(i,j))
!                fSum = fSum + ProbArray(i)*exp(BiasStorage(i,j) - maxBias)
            endif
          enddo
!            F_Estimate(j) = log(fSum)
          F_Estimate(j) = log(fSum)
          F_Estimate(j) = (F_Estimate(j) + F_Old(j))*0.5E0
       enddo 
!         Calculate the average change in F from the previous estimate and determine 
!         if there has been a significant change to the F values.
       tol = 0E0
       do j = 1, self%nCurWhamItter
         tol = tol + abs(F_Estimate(j) - F_Old(j))
       enddo
     enddo

!        Using the new estimates for the unbiased probability, calculate the free energy of nucleation
!        and modify the umbrella sampling bias to
      self%NewBias = 0E0
      maxbin = maxloc(self%HistStorage,1)
      maxbin2 = maxloc(self%TempHist,1)
      do i = 1, self%umbrellaLimit
        if(self%ProbArray(i) > 0E0) then
          self%FreeEnergyEst(i) = -log(self%ProbArray(i)/self%ProbArray(maxbin))
        endif
        if(self%TempHist(i) >= 1E0) then
          self%NewBias(i) = self%UBias(i) - self%UBias(maxbin2) - log(self%TempHist(i)/self%TempHist(maxbin2))
        endif
      enddo
      do i = 1, self%umbrellaLimit
        if(self%TempHist(i) .lt. 1E0) then
          self%NewBias(i) = self%UBias(i) - self%UBias(maxbin2) + log(self%TempHist(maxbin2))
        endif
      enddo
!        Rescale the pontential such that the reference free energy is set to 0
      refBias = self%NewBias(self%refBin)
      do i = 1, self%umbrellaLimit
        self%NewBias(i) = self%NewBias(i) - refBias
      enddo
      refBias = self%FreeEnergyEst(self%refBin)
      do i = 1, self%umbrellaLimit
        self%FreeEnergyEst(i) = self%FreeEnergyEst(i) - refBias
      enddo
!      endif

    endif      !End of processor 0 only block


    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
!      Distribute the new free energy estimate to all threads so that they can continue the simulation
!      with the new free energy. 
    arraySize = size(self%NewBias)      
    call MPI_BCast(self%NewBias, arraySize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 

    do i = 1, self%umbrellaLimit
!        write(*,*) i, UBias(i), Newbias(i)
      self%UBias(i) = self%NewBias(i)
      self%UHist(i) = 0E0
    enddo 
    if(myid .eq. 0) then
      call self%WHAMSimOutput
    endif

    self%nCurWhamItter = self%nCurWhamItter + 1
     
  end subroutine
!====================================================================
end module
!====================================================================
