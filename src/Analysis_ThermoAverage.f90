!=========================================================================
module Anaylsis_ThermoAverage
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: ThermoAverage
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
      procedure, pass :: Compute => Thermo_Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => Thermo_ProcessIO
      procedure, pass :: WriteInfo => Thermo_WriteInfo
      procedure, pass :: GetResult => Thermo_GetResult
      procedure, pass :: GetSTDev => Thermo_GetSTDev
      procedure, pass :: Epilogue => Thermo_Epilogue
      procedure, pass :: Finalize => Thermo_Finalize
      procedure, pass :: CastCommonType => Thermo_CastCommonType
  end type

 contains
!=========================================================================
  subroutine Thermo_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(ThermoAverage), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: thermVal

    thermVal = BoxArray(self%boxNum) % box % GetThermo(self%thermNum)

    self%varSum = self%varSum + thermVal
    self%varSumSq = self%varSumSq + thermVal*thermVal
    self%nSamp = self%nSamp + 1E0_dp

  end subroutine
!=========================================================================
  subroutine Thermo_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(ThermoAverage), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal


    !Input format
    ! ThermoAverage (BoxNum) (Thermo Variable) (Avg Frequency)
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
  subroutine Thermo_WriteInfo(self)
    use ParallelVar, only: nout
    implicit none
    class(ThermoAverage), intent(inout) :: self

    write(nout, *) trim(adjustl(self%varName)), " Average: ", self%GetResult(), "+\-", self%GetSTDev()
  end subroutine
!=========================================================================
  function Thermo_GetResult(self) result(var)
    implicit none
    class(ThermoAverage), intent(in) :: self
    real(dp) :: var

    var = self%varSum/self%nSamp
  end function
!=========================================================================
  function Thermo_GetSTDev(self) result(var)
    implicit none
    class(ThermoAverage), intent(in) :: self
    real(dp) :: var

    if(self%nSamp < 1E0_dp) then
      var = 0.0E0_dp
      return
    endif

    var = self%varSumSq/self%nSamp - (self%varSum/self%nSamp)**2
    var = sqrt(var)
  end function
!=========================================================================
  subroutine Thermo_Epilogue(self)
    use ParallelVar, only: myid
    implicit none
    class(ThermoAverage), intent(inout) :: self

    call self % WriteInfo

  end subroutine
!=========================================================================
  subroutine Thermo_Finalize(self)
#ifdef PARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout
    implicit none
    integer :: ierror
    class(ThermoAverage), intent(inout) :: self
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

    if(myid .eq. 0) then
      write(nout, *) trim(adjustl(self%varName)), " Thread Average: ", self%GetResult(), "+\-", self%GetSTDev()
!      call self % WriteInfo
    endif
#endif

  end subroutine
!=========================================================================
  subroutine Thermo_CastCommonType(self, anaVar)
    implicit none
    class(ThermoAverage), intent(inout) :: self
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
