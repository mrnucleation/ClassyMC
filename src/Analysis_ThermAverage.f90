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
    real(dp) :: varSum = 0E0_dp
    real(dp) :: nSamp = 1E-40_dp
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Compute => Therm_Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => Therm_ProcessIO
      procedure, pass :: WriteInfo => Therm_WriteInfo
      procedure, pass :: GetResult => Therm_GetResult
      procedure, pass :: Finalize => Therm_Finalize
  end type

 contains
!=========================================================================
  subroutine Therm_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(ThermAverage), intent(in) :: self
    logical, intent(in) :: accept


  end subroutine
!=========================================================================
  subroutine Therm_ProcessIO(self)
    implicit none
    class(ThermAverage), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine Therm_WriteInfo(self)
    implicit none
    class(ThermAverage), intent(inout) :: self
  end subroutine
!=========================================================================
  function Therm_GetResult(self) result(var)
    implicit none
    class(ThermAverage), intent(in) :: self
    real(dp) :: var

    var = self%varSum/self%nSamp
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

    dummy = self%nSamp
    call MPI_REDUCE(self%nSamp, dummy, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    if(myid .eq. 0) then
      
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
