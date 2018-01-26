!=========================================================================
module Move_ThermoLambda
  use SimpleSimBox, only: SimpleBox
  use CoordinateTypes, only: Displacement
  use VarPrecision
  use MoveClassDef

  use FF_ThermoIntegration, only: Pair_ThermoIntegration
  use ForcefieldData, only: EnergyCalculator
  use AnalysisData, only: AnalysisArray


  type, public, extends(MCMove) :: ThermoLambda
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    integer :: AnalyFunc = -1
    integer :: EFunc = -1
    type(Displacement) :: disp(1:1)
    contains
      procedure, pass :: Constructor => ThermoLambda_Constructor
!      procedure, pass :: GeneratePosition => ThermoLambda_GeneratePosition
      procedure, pass :: FullMove => ThermoLambda_FullMove
!      procedure, pass :: GetAcceptRate
!      procedure, pass :: Maintenance => ThermoLambda_Maintenance
      procedure, pass :: Epilogue => ThermoLambda_Epilogue
  end type

 contains
!========================================================
  subroutine ThermoLambda_Constructor(self)
    use AnalysisData, only: AnalysisArray
    use ForcefieldData, only: EnergyCalculator
    use Anaylsis_ThermoIntegration, only: ThermoIntegration
    implicit none
    class(ThermoLambda), intent(inout) :: self
    integer :: i

    do i = 1, size(EnergyCalculator)
      select type(eng => EnergyCalculator(i)%Method)
        class is(Pair_ThermoIntegration)
          self%EFunc = i
          exit
      end select
    enddo

    if(self%Efunc <= 0) then
      write(*,*) "ERROR! To perform thermodynamical integration the cooresponding"
      write(*,*) "forcefield function must be defined."
      stop
    endif

    do i = 1, size(AnalysisArray)
      select type(analy => AnalysisArray(i)%func)
        class is(ThermoIntegration)
          self%AnalyFunc = i
          exit
      end select
    enddo

    if(self%AnalyFunc <= 0) then
      write(*,*) "ERROR! To perform thermodynamical integration the cooresponding"
      write(*,*) "analysis function must be defined."
      stop
    endif


  end subroutine
!=========================================================================
  subroutine ThermoLambda_FullMove(self, trialBox, accept)

    use AnalysisData, only: AnalysisArray
    use CommonSampling, only: Sampling
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd

    implicit none
    class(ThermoLambda), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i
    integer :: nAtom, nAtomNew, reduIndx, newtype, oldtype
    real(dp) :: OldProb, NewProb, Prob
    real(dp) :: E_Diff, lambdaNew


    self % atmps = self % atmps + 1E0_dp
    lambdaNew = grnd()

    call AnalysisArray(self%AnalyFunc)%func%CalcNewState(newVal=lambdaNew)

    select type(eng => EnergyCalculator(self%EFunc)%Method)
      class is(Pair_ThermoIntegration)
        call eng%LambdaShift(lambdaNew, E_Diff)
    end select

    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      select type(eng => EnergyCalculator(self%EFunc)%Method)
        class is(Pair_ThermoIntegration)
          call eng%UpdateLambda(lambdaNew)
      end select


    endif

  end subroutine
!=========================================================================
  subroutine ThermoLambda_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(ThermoLambda), intent(inout) :: self
    real(dp) :: accptRate
      
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,I15)") "Thermo Lambda Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Thermo Lambda Moves Attempted: ", nint(self%atmps)
    write(nout,"(1x,A,F15.8)") "Thermo Lambda Acceptance Rate: ", accptRate

  end subroutine
!=========================================================================
end module
!=========================================================================
