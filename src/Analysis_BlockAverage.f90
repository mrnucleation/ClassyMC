!=========================================================================
module Anaylsis_BlockAverage
  use AnaylsisClassDef, only: Analysis
  use VarPrecision



  type, public, extends(Analysis):: BlockAverage
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

#ifdef MPIPARALLEL
    logical :: parallel = .true.
#else
    logical :: parallel = .false.
#endif

    integer :: fileunit = 300
    integer :: boxNum = -1
    integer :: thermNum = -1
    integer :: writeNum = 0
    character(len=30) :: varName = ""
    character(len=80) :: fileName = ""

    real(dp) :: varSum = 0E0_dp
    real(dp) :: varSumSq = 0E0_dp
    real(dp) :: nSamp = 1E-40_dp
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Compute => BlockAvg_Compute
      procedure, pass :: Maintenance => BlockAvg_Maintenance
      procedure, pass :: ProcessIO => BlockAvg_ProcessIO
!      procedure, pass :: WriteInfo => BlockAvg_WriteInfo
!      procedure, pass :: GetResult => BlockAvg_GetResult
!      procedure, pass :: Finalize => BlockAvg_Finalize
      procedure, pass :: CastCommonType => BlockAvg_CastCommonType
  end type

 contains
!=========================================================================
  subroutine BlockAvg_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(BlockAverage), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: thermVal

    thermVal = BoxArray(self%boxNum) % box % GetThermo(self%thermNum)

    self%varSum = self%varSum + thermVal
    self%varSumSq = self%varSumSq + thermVal**2
    self%nSamp = self%nSamp + 1E0_dp

  end subroutine
!=========================================================================
  subroutine BlockAvg_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText
    use ParallelVar, only: myid, nout
    implicit none
    class(BlockAverage), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal, iCharacter
    character(len=30) :: idString


    !Input format
    ! BlockAverage (BoxNum) (BlockAvg Variable) (Avg Frequency) (Block Size) (FileName)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal

    call GetXCommand(line, command, 3, lineStat)
    self%thermNum = BoxArray(self%boxNum) % box % ThermoLookUp(command)
    self%varName = command
!    write(nout,*) command
 
    call GetXCommand(line, command, 4, lineStat)
    read(command, *) intVal
    self%UpdateFreq = intVal

    call GetXCommand(line, command, 5, lineStat)
    read(command, *) intVal
    self%MaintFreq = intVal

    call GetXCommand(line, command, 6, lineStat)
    read(command, *) self%fileName
    do iCharacter = 1, len(self%filename)
      if(self%filename(iCharacter:iCharacter) == "&") then
        write(idString, *) myid
        self%filename = ReplaceText(self%filename, "&", trim(adjustl(idString)))
#ifdef MPIPARALLEL
        self%parallel = .false.
#endif
        exit
      endif
    enddo

    do iCharacter = 1, len(self%filename)
      if(self%filename(iCharacter:iCharacter) == '"') then
        self%fileName(iCharacter:iCharacter) = " "
      endif
    enddo


    open(newunit=self%fileunit, file=trim(adjustl(self%fileName)) )

  end subroutine
!=========================================================================
  subroutine BlockAvg_Maintenance(self)
#ifdef MPIPARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout
    implicit none
    integer :: ierror
    class(BlockAverage), intent(inout) :: self
    real(dp) :: dummy


    self%writeNum = self%writeNum + 1
#ifdef MPIPARALLEL
    
    if(self%parallel) then
        write(nout, *) "Stopping for block averaging"
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)       

        dummy = self%varSum
        call MPI_REDUCE(self%varSum, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

        dummy = self%varSumSq
        call MPI_REDUCE(self%varSumSq, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

        dummy = self%nSamp
        call MPI_REDUCE(self%nSamp, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    endif
#endif

    if(self%parallel) then
      if(myid /= 0 ) then
        self%varSum = 0E0_dp
        self%varSumSq = 0E0_dp
        self%nSamp = 1E-40_dp
        return
      endif
    endif
    write(self%fileunit, *) self%writeNum*self%maintFreq, self%varSum/self%nSamp, self%varSumSq/self%nSamp
    flush(self%fileunit)
    self%varSum = 0E0_dp
    self%varSumSq = 0E0_dp
    self%nSamp = 1E-40_dp
  end subroutine
!=========================================================================
!  subroutine BlockAvg_WriteInfo(self)
!    use ParallelVar, only: nout
!    implicit none
!    class(BlockAverage), intent(inout) :: self
!
!    write(nout, *) trim(adjustl(self%varName)), " Average: ", self%GetResult(), "+\-", self%GetSTDev()
!  end subroutine
!=========================================================================
  subroutine BlockAvg_CastCommonType(self, anaVar)
    implicit none
    class(BlockAverage), intent(inout) :: self
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
