!====================================================================
module UmbrellaRule
  use VarPrecision
  use AcceptRuleTemplate, only: AcceptRule
  use CoordinateTypes, only: Displacement, Perturbation, Addition, Deletion, VolChange
 
  type, public, extends(AcceptRule) :: Umbrella

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

    contains
       procedure, pass :: Constructor => Umbrella_Constructor
       procedure, pass :: MakeDecision => Umbrella_MakeDecision
       procedure, pass :: GetBiasIndex => Umbrella_GetBiasIndex
       procedure, pass :: GetNewBiasIndex => Umbrella_GetNewBiasIndex
       procedure, pass :: ReadInitialBias => Umbrella_ReadInitialBias
       procedure, pass :: GetUIndexArray => Umbrella_GetUIndexArray
       procedure, pass :: OutputUmbrellaHist => Umbrella_OutputUmbrellaHist
       procedure, pass :: FindVarValues => Umbrella_FindVarValues
       procedure, pass :: OutBias => Umbrella_OutBias
!       procedure, pass :: Maintenance => Umbrella_Maintenance

       procedure, pass :: ProcessIO => Umbrella_ProcessIO
       procedure, pass :: Epilogue => Umbrella_Epilogue
       procedure, pass :: Prologue => Umbrella_Prologue

  end type
!====================================================================
  contains
!====================================================================
  subroutine Umbrella_Constructor(self)
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: AllocationStat

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
  subroutine Umbrella_Prologue(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout
    implicit none
    class(Umbrella), intent(inout) :: self
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


  end subroutine
!====================================================================
  function Umbrella_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use AnalysisData, only: AnalysisArray
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(Umbrella), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
!    type(Displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: inProb
    real(dp), intent(in) :: E_Diff

    logical :: accept
    integer :: iBias, oldIndx, newIndx, indx
    real(dp) :: biasE, biasOld, biasNew
    real(dp) :: extraTerms, chemPot

    do iBias = 1, self%nBiasVar
      indx = self%AnalysisIndex(iBias)
      call AnalysisArray(indx)%func%CalcNewState(disp)
    enddo

    oldIndx = self%GetBiasIndex()
    call self%GetNewBiasIndex(newIndx, accept)

    extraTerms = 0E0_dp
    select type(disp)
      class is(Addition)
          extraTerms = extraTerms + trialBox%chempot(disp(1)%molType)
      class is(Deletion)
          extraTerms = extraTerms - trialBox%chempot(disp(1)%molType)
      class is(VolChange)
          extraTerms = extraTerms + (disp(1)%volNew -disp(1)%volOld)*trialBox%pressure*trialBox%beta
    end select



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
      self%binIndx(iBias) = floor(biasVal/self%UBinSize(iBias))
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
  subroutine Umbrella_GetNewBiasIndex(self, biasIndx, accept)
    use AnalysisData, only: AnalysisArray, analyCommon
    implicit none
    class(Umbrella), intent(inout) :: self
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

      select type( biasVar => analyCommon(analyIndx)%val )
        type is(integer)
            biasVal = biasVar
        type is(real)
            biasVal = biasVar
      end select

      self%binIndx(iBias) = floor(biasVal/self%UBinSize(iBias))

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
      self%binIndx(iBias) = floor( varArray(iBias) / self%UBinSize(iBias) )
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
   subroutine Umbrella_OutputUmbrellaHist(self)
     implicit none
     class(Umbrella), intent(inout) :: self
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
  subroutine Umbrella_FindVarValues(self, UIndx, UArray) 
    implicit none 
    class(Umbrella), intent(inout) :: self
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
  subroutine Umbrella_ProcessIO(self, line, linestat) 
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(Umbrella), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line   
    integer, intent(out) :: lineStat

    integer :: intVal, intVal2
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
          self%binMin(intVal) = floor(realVal/self%UBinSize(intVal))


          call GetXCommand(line, command, 7, lineStat)
          read(command, *) realVal
          self%binMax(intVal) = floor(realVal/self%UBinSize(intVal))

        endif

      case default
        lineStat = -1
    end select



   end subroutine
!====================================================================
  subroutine Umbrella_ColapseHist(self)
    use ParallelVar, only: nout, myid, ierror
    use MPI
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: iUmbrella, arraySize
    real(dp), allocatable :: TempHist(:)



    call MPI_BARRIER(MPI_COMM_WORLD, ierror)  
    allocate( TempHist(1:self%umbrellaLimit+1) ) 
    TempHist = 0E0_dp

    arraySize = size(self%UHist)
    call MPI_REDUCE(self%UHist, TempHist, arraySize, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

    if(myid == 0) then
      do iUmbrella = 1, self%umbrellaLimit
        self%UHist(iUmbrella) = TempHist(iUmbrella)
      enddo
    endif

    deallocate(TempHist)

  end subroutine
!====================================================================
  subroutine Umbrella_Epilogue(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(Umbrella), intent(inout) :: self

    write(nout,*) "Writing Umbrella Sampling Histogram..."
    if(myid == 0) then
      call self%OutputUmbrellaHist
      call self%OutBias
    endif
  end subroutine
!====================================================================
  subroutine Umbrella_Outbias(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: refBin
    real(dp) :: newBias, refbias
     integer :: iUmbrella, iBias, iBin
     character(len = 100) :: outputString



!     write(nout,*) "Outbias"
     write(outputString, *) "(", ("F12.5, 2x", iBias =1,self%nBiasVar), "2x, F18.7)"
     open(unit=60, file="NewBias.txt")
      
     refbin = 1
     do iUmbrella = 1, self%umbrellaLimit
       if(self%UHist(iUmbrella) > self%UHist(refBin)) then
         refBin = iUmbrella
         refbias = self%UBias(refBin)
       endif
     enddo

     do iUmbrella = 1, self%umbrellaLimit
       call self%findVarValues(iUmbrella, self%UArray)  
       do iBias = 1, self%nBiasVar        
         self%varValues(iBias) = real(self% UArray(iBias), dp) * self%UBinSize(iBias)          
       enddo

       if(self%UHist(iUmbrella) /= 0.0E0_dp) then
         newBias = self%UBias(iUmbrella) - refBias - &
                   log(self%UHist(iUmbrella)/self%UHist(refBin))
       else
         newBias = self%UBias(iUmbrella) - refBias + log(self%UHist(refBin))

       endif

       write(60,outputString) (self%varValues(iBias), iBias=1,self%nBiasVar), newBias
        
     enddo 
    
     flush(60)
     close(60)
 
  end subroutine
!====================================================================
end module
!====================================================================
