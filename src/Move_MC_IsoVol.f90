!=========================================================================
module Move_IsoVol
use SimpleSimBox, only: SimpleBox
use VarPrecision
use MoveClassDef

  type, public, extends(MCMove) :: IsoVol
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp   
    real(dp) :: maxDv
    integer :: nDim = 3
    real(dp), allocatable :: scalar(:)
    type(Displacement) :: disp(1:1)
    contains
      procedure, pass :: Constructor => IsoVol_Constructor
      procedure, pass :: GeneratePosition => IsoVol_GeneratePosition
      procedure, pass :: FullMove => IsoVol_FullMove
!      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance => IsoVol_Maintenance
  end type

 contains
!========================================================
  subroutine IsoVol_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(IsoVol), intent(inout) :: self

    allocate( self%scalar(1:self%nDim) )

  end subroutine
!=========================================================================
  subroutine IsoVol_GeneratePosition(self, disp)
    use CoordinateTypes, only: Displacement
    implicit none
    class(IsoVol), intent(in) :: self
  end subroutine
!=========================================================================
  subroutine IsoVol_FullMove(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling

    implicit none
    class(IsoVol), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i
    real(dp) :: dV
    real(dp) :: OldProb, NewProb, Prob
    real(dp) :: E_Diff

    self % atmps = self % atmps + 1E0_dp
    dV = self%maxDv * (2E0_dp*grnd()-1E0_dp)

    do i = 1, size(self%scalar)
      self%scalar(i) = 1E0_dp + dV/trialBox%volume
    enddo

    !Check Constraint
!    accept = trialBox % CheckConstraint( self%disp(1:1) )
!    if(.not. accept) then
!      return
!    endif

    call trialbox % EFunc % Method % VolECalc(trialbox, self%scalars, E_Diff)

    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdateVol(self%scalars)
    endif

  end subroutine
!=========================================================================
  subroutine IsoVol_Maintenance(self)
    class(IsoVol), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
