!=========================================================================
module Anaylsis_ThermAverage
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: ThermAverage
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
      procedure, pass :: Compute => Therm_Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => Therm_ProcessIO
      procedure, pass :: WriteInfo => Therm_WriteInfo
      procedure, pass :: GetResult => Therm_GetResult
      procedure, pass :: GetSTDev => Therm_GetSTDev
      procedure, pass :: Finalize => Therm_Finalize
  end type

 contains
!=========================================================================
  subroutine Therm_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(ThermAverage), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: thermVal

    thermVal = BoxArray(self%boxNum) % box % GetThermo(self%thermNum)

    self%varSum = self%varSum + thermVal
    self%varSumSq = self%varSumSq + thermVal*thermVal
    self%nSamp = self%nSamp + 1E0_dp

  end subroutine
!=========================================================================
  subroutine Therm_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(ThermAverage), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal

    call GetXCommand(line, command, 3, lineStat)
    self%thermNum = BoxArray(self%boxNum) % box % ThermoLookUp(command)
    self%varName = command
    
    call GetXCommand(line, command, 4, lineStat)
    read(command, *) intVal
    self%UpdateFreq = intVal


  end subroutine
!=========================================================================
  subroutine Therm_WriteInfo(self)
    use ParallelVar, only: nout
    implicit none
    class(ThermAverage), intent(inout) :: self

    write(nout, *) trim(adjustl(self%varName)), " Average: ", self%GetResult(), "+\-", self%GetSTDev()
  end subroutine
!=========================================================================
  function Therm_GetResult(self) result(var)
    implicit none
    class(ThermAverage), intent(in) :: self
    real(dp) :: var

    var = self%varSum/self%nSamp
  end function
!=========================================================================
  function Therm_GetSTDev(self) result(var)
    implicit none
    class(ThermAverage), intent(in) :: self
    real(dp) :: var

    var = self%varSumSq/self%nSamp - (self%varSum/self%nSamp)**2
    var = sqrt(var)
  end function
!=========================================================================
  subroutine Therm_Finalize(self)
    use MPI
    use ParallelVar, only: myid
    implicit none
    integer :: ierror
    class(ThermAverage), intent(inout) :: self
    real(dp) :: dummy

    dummy = self%varSum
    call MPI_REDUCE(self%varSum, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    dummy = self%varSumSq
    call MPI_REDUCE(self%varSumSq, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    dummy = self%nSamp
    call MPI_REDUCE(self%nSamp, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    if(myid .eq. 0) then
      call self % WriteInfo
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
