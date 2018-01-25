!=========================================================================
module Move_ThermoLambda
  use SimpleSimBox, only: SimpleBox
  use CoordinateTypes, only: Displacement
  use VarPrecision
  use MoveClassDef

  use FF_ThermoIntegration, only: ThermoIntegration
  use ForcefieldData, only: EnergyCalculator
  use AnalysisData, only: AnalysisArray


  type, public, extends(MCMove) :: ThermoLambda
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    integer :: AnalyFunc = -1
    integer :: EFunc = -1
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
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(ThermoLambda), intent(inout) :: self


  end subroutine
!=========================================================================
  subroutine ThermoLambda_FullMove(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin

    implicit none
    class(ThermoLambda), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i
    integer :: nAtom, nAtomNew, reduIndx, newtype, oldtype
    real(dp) :: OldProb, NewProb, Prob
    real(dp) :: E_Diff


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
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
