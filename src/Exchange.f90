!========================================================
module MCM_Exchange
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimBoxDef, only: SimBox
use VarPrecision

  type, public, extends(MCMove) :: AtomExchange
    real(dp) :: max_dist = 0.1E0_dp
    type(Displacement) :: disp(1:1)
    contains
      procedure, pass :: Constructor => Exchange_Constructor
      procedure, pass :: GeneratePosition => Exchange_GeneratePosition
      procedure, pass :: FullMove => Exchange_FullMove
  end type
!========================================================
 contains
!========================================================
  subroutine Exchange_Constructor(self)
    implicit none
    class(AtomExchange), intent(inout) :: self

  end subroutine
!========================================================
  subroutine Exchange_GeneratePosition(self, disp)
    use RandomGen, only: grnd
    implicit none
    class(AtomExchange), intent(in) :: self
    type(Displacement), intent(inout) :: disp
    real(dp) :: dx, dy, dz

  end subroutine
!===============================================
  subroutine Exchange_FullMove(self, trialBox) 
    use ForcefieldData, only: EnergyCalculator
    use BoxData, only: Constrain
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    implicit none
    class(AtomExchange), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    logical :: accept
    integer :: nMove, iConstrain
    integer :: CalcIndex
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
       !Choose Molecule
    nMove = floor( trialBox%nAtoms * grnd() + 1E0_dp)

 
    self%disp(1)%atmIndx = nMove

    !Check Constraint
    if( size(Constrain) > 0 ) then
      do iConstrain = 1, size(Constrain)
        call Constrain(iConstrain) % method % ShiftCheck( trialBox, self%disp(1:1), accept )
      enddo
      if(.not. accept) then
        return
      endif
    endif


    !Energy Calculation
    CalcIndex = trialBox % ECalcer
    call EnergyCalculator(CalcIndex) % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)


    !Accept/Reject
    accept = .false.
    accept = sampling % MakeDecision(trialBox, E_Diff, 1E0_dp, self%disp(1:1))
    if(accept) then
      self % accpt = self % atmps + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:1))
    endif

  end subroutine

!========================================================
end module
!========================================================
