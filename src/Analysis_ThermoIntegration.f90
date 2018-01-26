!=========================================================================
module Anaylsis_ThermoIntegration
use AnaylsisClassDef, only: Analysis
use AnalysisData, only: analyCommon
use FF_ThermoIntegration, only: Pair_ThermoIntegration
use VarPrecision

  type, public, extends(Analysis):: ThermoIntegration
    
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1
!    integer :: analyID = -1

    integer :: ECalc = -1
    real(dp) :: lambda

    contains
      procedure, pass :: Initialize => ThermoInt_Initialize
      procedure, pass :: Compute => ThermoInt_Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: CalcNewState => ThermoInt_CalcNewState
!      procedure, pass :: ProcessIO => ThermoInt_ProcessIO
!      procedure, pass :: WriteInfo => ThermoInt_WriteInfo
      procedure, pass :: GetResult => ThermoInt_GetResult
!      procedure, pass :: Finalize => ThermoInt_Finalize
  end type

 contains
!=========================================================================
  subroutine ThermoInt_Initialize(self)
    use ForcefieldData, only: EnergyCalculator
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    integer :: i

    do i = 1, size(EnergyCalculator)
      select type(eng => EnergyCalculator(i)%Method)
        class is(pair_ThermoIntegration)
          self%ECalc = i
          exit
      end select
    enddo


  end subroutine
!=========================================================================
  subroutine ThermoInt_Compute(self, accept)
    use AnalysisData, only: analyCommon
    use ForcefieldData, only: EnergyCalculator
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: lambda


    if(accept) then
      select type(eng => EnergyCalculator(self%ECalc)%Method)
        class is(pair_ThermoIntegration)
           self%lambda = eng%GetLambda()
      end select
    endif

    
  end subroutine
!=========================================================================
  subroutine ThermoInt_CalcNewState(self, disp, newVal)
    use AnalysisData, only: analyCommon
    use CoordinateTypes, only: Displacement
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    type(Displacement), intent(in), optional :: disp(:)
    real(dp), intent(in), optional :: newVal

    if( present(newVal) ) then
      analyCommon(self%analyID) = newVal
    endif
  end subroutine
!=========================================================================
!  subroutine ThermoInt_ProcessIO(self, line)
!    use Input_Format, only: maxLineLen, GetXCommand
!    implicit none
!    class(ThermoIntegration), intent(inout) :: self
!    character(len=maxLineLen), intent(in) :: line
!    character(len=30) :: command
!    integer :: lineStat = 0
!    integer :: intVal
!
!
!
!  end subroutine
!=========================================================================
!  subroutine ThermoInt_WriteInfo(self)
!    use ParallelVar, only: nout
!    implicit none
!    class(ThermoIntegration), intent(inout) :: self
!
!  end subroutine
!=========================================================================
  function ThermoInt_GetResult(self) result(var)
    implicit none
    class(ThermoIntegration), intent(in) :: self
    real(dp) :: var

    var = self%lambda
  end function
!=========================================================================
!  subroutine ThermoInt_Finalize(self)
!    use MPI
!    use ParallelVar, only: myid
!    implicit none
!    integer :: ierror
!    class(ThermoIntegration), intent(inout) :: self
!    real(dp) :: dummy
!
!
!    if(myid .eq. 0) then
!      call self % WriteInfo
!    endif
!
!  end subroutine
!=========================================================================
end module
!=========================================================================
