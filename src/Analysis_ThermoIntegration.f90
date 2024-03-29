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

    logical :: pushedValue = .false.
    integer :: ECalc = -1
    real(dp) :: lambda, newLambda

    contains
!      procedure, pass :: Initialize => ThermoInt_Initialize
      procedure, pass :: Prologue => ThermoInt_Prologue
      procedure, pass :: PushLambda => ThermoInt_PushLambda
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
  subroutine ThermoInt_Prologue(self)
    use ForcefieldData, only: EnergyCalculator
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    integer :: i

    self%perMove = .true.
    do i = 1, size(EnergyCalculator)
      select type(eng => EnergyCalculator(i)%Method)
        class is(pair_ThermoIntegration)
          self%ECalc = i
          exit
      end select
    enddo

    if(self%ECalc < 1) then
      write(0,*) "ERROR! The ThermoIntegration Analysis function must be used"
      write(0,*) "with the corresponding energy function!"
      error stop
    endif

  end subroutine
!=========================================================================
  subroutine ThermoInt_Compute(self, accept)
    use AnalysisData, only: analyCommon
    use ForcefieldData, only: EnergyCalculator
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    logical, intent(in) :: accept
    real(dp) :: lambda


    select type(eng => EnergyCalculator(self%ECalc)%Method)
      class is(pair_ThermoIntegration)
         self%lambda = eng%GetLambda()
    end select

    
  end subroutine
!=========================================================================
  subroutine ThermoInt_PushLambda(self, newVal)
    use AnalysisData, only: analyCommon
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    real(dp), intent(in) :: newVal

    self%newLambda = newVal
    self%pushedValue = .true.

  end subroutine
!=========================================================================
  subroutine ThermoInt_CalcNewState(self, disp, accept, newVal)
    use AnalysisData, only: analyCommon
    use CoordinateTypes, only: Perturbation
    implicit none
    class(ThermoIntegration), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    real(dp), intent(in), optional :: newVal
    logical, intent(out) :: accept

    accept = .true.

    if( self%pushedValue ) then
      select type( anaVar => analyCommon(self%analyID)%val )
        type is(real)
          anaVar = self%newLambda
      end select

      self%pushedValue = .false.
    else
      select type( anaVar => analyCommon(self%analyID)%val )
        type is(real)
           anaVar = self%lambda
      end select
!      analyCommon(self%analyID)%val  = self%lambda
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
