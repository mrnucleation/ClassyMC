!=========================================================================
module Anaylsis_BlockAvgAverage
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: BlockAvgAverage
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    integer :: boxNum = -1
    integer :: thermNum = -1
    character(len=30) :: varName

    real(dp) :: varSum = 0E0_dp
    real(dp) :: varSumSq = 0E0_dp
    real(dp) :: nSamp = 1E-40_dp
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Compute => BlockAvg_Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => BlockAvg_ProcessIO
      procedure, pass :: WriteInfo => BlockAvg_WriteInfo
      procedure, pass :: GetResult => BlockAvg_GetResult
      procedure, pass :: GetSTDev => BlockAvg_GetSTDev
      procedure, pass :: Finalize => BlockAvg_Finalize
      procedure, pass :: CastCommonType => BlockAvg_CastCommonType
  end type

 contains
!=========================================================================
  subroutine BlockAvg_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(BlockAvgAverage), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: thermVal

    thermVal = BoxArray(self%boxNum) % box % GetBlockAvg(self%thermNum)

    self%varSum = self%varSum + thermVal
    self%varSumSq = self%varSumSq + thermVal*thermVal
    self%nSamp = self%nSamp + 1E0_dp

  end subroutine
!=========================================================================
  subroutine BlockAvg_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(BlockAvgAverage), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal


    !Input format
    ! BlockAvgAverage (BoxNum) (BlockAvg Variable) (Avg Frequency) (FileName)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal

    call GetXCommand(line, command, 3, lineStat)
    self%thermNum = BoxArray(self%boxNum) % box % BlockAvgLookUp(command)
    self%varName = command
    
    call GetXCommand(line, command, 4, lineStat)
    read(command, *) intVal
    self%UpdateFreq = intVal

    do iCharacter = 1, len(self%filename)
      if(self%filename(iCharacter:iCharacter) == "&") then
        write(idString, *) myid
        self%filename = ReplaceText(self%filename, "&", trim(adjustl(idString)))
        exit
      endif
    enddo


  end subroutine
!=========================================================================
  subroutine BlockAvg_WriteInfo(self)
    use ParallelVar, only: nout
    implicit none
    class(BlockAvgAverage), intent(inout) :: self

    write(nout, *) trim(adjustl(self%varName)), " Average: ", self%GetResult(), "+\-", self%GetSTDev()
  end subroutine
!=========================================================================
  function BlockAvg_GetResult(self) result(var)
    implicit none
    class(BlockAvgAverage), intent(in) :: self
    real(dp) :: var

    var = self%varSum/self%nSamp
  end function
!=========================================================================
  subroutine BlockAvg_Finalize(self)
#ifdef PARALLEL
    use MPI
#endif
    use ParallelVar, only: myid
    implicit none
    integer :: ierror
    class(BlockAvgAverage), intent(inout) :: self
    real(dp) :: dummy

#ifdef PARALLEL
    dummy = self%varSum
    call MPI_REDUCE(self%varSum, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    dummy = self%varSumSq
    call MPI_REDUCE(self%varSumSq, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    dummy = self%nSamp
    call MPI_REDUCE(self%nSamp, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 
#endif

    if(myid .eq. 0) then
      call self % WriteInfo
    endif

  end subroutine
!=========================================================================
  subroutine BlockAvg_CastCommonType(self, anaVar)
    implicit none
    class(BlockAvgAverage), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
!      write(*,*) "Allocated as Real"
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
