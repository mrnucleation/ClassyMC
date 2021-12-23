!====================================================================
!
!
!
!====================================================================
module UmbrellaRule
  use VarPrecision
  use Template_AcceptRule, only: acceptrule
  use CoordinateTypes, only: Displacement, Perturbation, Addition, Deletion, VolChange
  use UmbrellaWHAMRule, only: UmbrellaWHAM
 
  !--------------------------------------------------------------
  !Since the two classes have a lot of opperations in common, it's
  !cleaner to extend the Umbrella class from the UmbrellaWHAM.
  !-----------------------------------------------------------
  type, public, extends(UmbrellaWHAM) :: Umbrella
    integer, private :: histunit = -1
    contains
!       procedure, pass :: Constructor => Umbrella_Constructor
!       procedure, pass :: MakeDecision => Umbrella_MakeDecision
!       procedure, pass :: UpdateStatistics => Umbrella_UpdateStatistics

!       procedure, pass :: GetBiasIndex => Umbrella_GetBiasIndex
!       procedure, pass, private :: GetNewBiasIndex => Umbrella_GetNewBiasIndex
!       procedure, pass, private :: ReadInitialBias => Umbrella_ReadInitialBias
!       procedure, pass, private :: GetUIndexArray => Umbrella_GetUIndexArray
!       procedure, pass, private :: FindVarValues => Umbrella_FindVarValues
!       procedure, pass :: OutputUmbrellaHist => Umbrella_OutputUmbrellaHist
!       procedure, pass :: AdjustHist => Umbrella_AdjustHist
       procedure, pass :: CollectHist => Umbrella_CollectHist
       procedure, pass :: Maintenance => Umbrella_Maintenance

!       procedure, pass :: ProcessIO => Umbrella_ProcessIO
       procedure, pass :: Epilogue => Umbrella_Epilogue
       procedure, pass :: Prologue => Umbrella_Prologue
!       procedure, pass :: Update => Umbrella_Update

  end type
!====================================================================
  contains
!====================================================================
  subroutine Umbrella_Prologue(self)
    use AnalysisData, only: AnalysisArray
    use SimControl, only: nCycles
    use ParallelVar, only: nout, myid
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: i,j, indx, stat, AllocateStatus

    do i = 1, self%nBiasVar
      indx = self%AnalysisIndex(i)
      if((indx > size(AnalysisArray)) .or. (indx < 1) ) then
        write(nout, *) "ERROR! The Umbrella Sampling routine has been directed to an "
        write(nout, *) "invalid Analysis fucntion"
        write(nout, *) "Chosen Function:", indx
        write(nout, *) "Number of Analysis Functions:", size(AnalysisArray)
        error stop "Error dectected in Umbrella Sampling"
      endif
    enddo

    do i = 1, self%nBiasVar
      if(self%valMin(i) > self%valMax(i) ) then
        write(nout,*) "ERROR! The given bounds for one of the umbrella variables does not make sense!"
        write(nout,*) "Smallest bin is larger than the largest bin" 
        write(nout,*) "Minimum Value:", self%valMin(i)
        write(nout,*) "Maximum Value:", self%valMax(i)
        error stop
      endif
    enddo

    do i = 1, self%nBiasVar
      indx = self%AnalysisIndex(i)
      AnalysisArray(indx)%func%usedInMove = .true.
      AnalysisArray(indx)%func%permove = .true.
    enddo

     ! Since the number of biasing variables is only known at run time, the bias matrix
     ! must be stored in a 1D array. The index coefficient variables are used to emulate a N dimensional matrix
     ! using the linear mapping equation of this form:
     ! U =  a1*x1 + a2*x2 + a3*x3 + .....
     ! Which maps a N dimension matrix onto a 1D array. 

    self%indexCoeff(1) = 1
    do i = 2, self%nBiasVar 
      self%indexCoeff(i) = 1
      do j = 1, i-1
        self%indexCoeff(i) = self%indexCoeff(i) + self%indexCoeff(j) * self%nBins(j) 
      enddo
    enddo      
    self%umbrellaLimit = 1
    do i = 1, self%nBiasVar 
      self%umbrellaLimit = self%umbrellaLimit + self%indexCoeff(i) * self%nBins(i)
    enddo

    write(nout,*) "Sampling Style: Histogram based Umbrella Sampling"
    write(nout,*) "Number of Umbrella Bins:", self%umbrellaLimit
       
    allocate(self%UBias(1:self%umbrellaLimit+1), STAT = AllocateStatus)
    allocate(self%UHist(1:self%umbrellaLimit+1), STAT = AllocateStatus)

!    write(nout,*) self%binMax
    self%UBias = 0E0_dp
    self%UHist = 0E0_dp

    call self%GetUIndexArray(self%refVals, i, stat)
    if(stat /= 0) then
      error stop
    endif
    self%refBin = i

    call self%ReadInitialBias
    self%nWhamItter = ceiling(dble(nCycles)/dble(self%maintFreq))
    ! Allocation of the WHAM variables
    open(newunit = self%histunit, file="Umbrella_Histogram.dat")


    allocate(self%TempHist(1:self%umbrellaLimit), STAT = AllocateStatus)      
    allocate(self%NewBias(1:self%umbrellaLimit), STAT = AllocateStatus)
    self%NewBias = 0E0
    self%TempHist = 0E0

    write(nout,*) self%refVals, i 
    write(nout,*) "Bin Size:", self%UBinSize
!    self%oldIndx = self % GetBiasIndex()
    self%oldIndx = 1
  end subroutine
!====================================================================
  subroutine Umbrella_Maintenance(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(Umbrella), intent(inout) :: self

    call self%Collecthist
  end subroutine
!=========================================================================
!     This subroutine periodically 
  subroutine Umbrella_CollectHist(self)
    use ParallelVar, only: myid, nout, ierror
#ifdef MPIPARALLEL
    use MPI
#endif
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: arraySize, i, j, cnt, maxbin, maxbin2
    real(dp) :: norm, maxBias, denomSum
    real(dp) :: tol, refBias
    character(len = 100) :: outputString

#ifdef MPIPARALLEL
    write(nout,*) "Halting for Histogram Collection"
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
#endif
    
!      This block condences the histogram data from all the different processors
!      into one collective array on the root (myid = 0) processor.        
#ifdef MPIPARALLEL
    arraySize = size(self%UHist)     
    if(myid .eq. 0) then
      self%TempHist = 0E0_dp
    endif
    call MPI_REDUCE(self%UHist, self%TempHist, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)       
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
#else
    self%TempHist = self%UHist
#endif

!    This block exports the histogram to a file
    write(outputString, *) "(", ("2x, F10.4,", j =1,self%nBiasVar), "2x, F22.1)"
    rewind(self%histunit)
    do i = 1, self%umbrellaLimit
      if(self%HistStorage(i) /= 0E0_dp ) then
        call self%FindVarValues(i, self%UArray)
        write(self%histunit, outputString) (self%UArray(j)*self%UBinSize(j)+self%valMin(j), j=1,self%nBiasVar), self%TempHist(i)
      endif
    enddo
    flush(self%histunit)
    
  end subroutine
!====================================================================
  subroutine Umbrella_Epilogue(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(Umbrella), intent(inout) :: self

    write(nout,*) "Writing Umbrella Sampling Histogram..."
    call self%CollectHist
  end subroutine
!====================================================================
end module
!====================================================================
