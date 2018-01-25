!=========================================================================
module Anaylsis_ThermoIntegration
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: ThermoIntegration
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    integer :: boxNum = -1
    integer :: thermNum = -1

    contains
!      procedure, pass :: Initialize
      procedure, pass :: Compute => ThermoInt_Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => ThermoInt_ProcessIO
      procedure, pass :: WriteInfo => ThermoInt_WriteInfo
      procedure, pass :: GetResult => ThermoInt_GetResult
      procedure, pass :: Finalize => ThermoInt_Finalize
  end type

 contains
!=========================================================================
  subroutine ThermoInt_Compute(self, accept)
    use ForcefieldData, only: EnergyCalculator
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: thermVal


  end subroutine
!=========================================================================
  subroutine ThermoInt_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(ThermoIntegration), intent(inout) :: self
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
  subroutine ThermoInt_WriteInfo(self)
    use ParallelVar, only: nout
    implicit none
    class(ThermoIntegration), intent(inout) :: self

    write(nout, *) trim(adjustl(self%varName)), " Average: ", self%GetResult(), "+\-", self%GetSTDev()
  end subroutine
!=========================================================================
  function ThermoInt_GetResult(self) result(var)
    implicit none
    class(ThermoIntegration), intent(in) :: self
    real(dp) :: var

    var = self%varSum/self%nSamp
  end function
!=========================================================================
  subroutine ThermoInt_Finalize(self)
    use MPI
    use ParallelVar, only: myid
    implicit none
    integer :: ierror
    class(ThermoIntegration), intent(inout) :: self
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
