!========================================================
module AtomTranslation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimBoxDef, only: SimBox
use VarPrecision

  type, public, extends(MCMove) :: AtomMolTranslate
    real(dp) :: max_dist = 0.1E0_dp
    type(Displacement) :: disp(1:1)
    contains
      procedure, pass :: Constructor => AtomTrans_Constructor
      procedure, pass :: GeneratePosition => AtomTrans_GeneratePosition
      procedure, pass :: FullMove => AtomTrans_FullMove
  end type
!========================================================
 contains
!========================================================
  subroutine AtomTrans_Constructor(self)
    implicit none
    class(AtomMolTranslate), intent(inout) :: self
    
!    allocate(self%disp(1:
  end subroutine
!========================================================
  subroutine AtomTrans_GeneratePosition(self, disp)
    use RandomGen, only: grnd
    implicit none
    class(AtomMolTranslate), intent(in) :: self
    type(Displacement), intent(inout) :: disp
    real(dp) :: dx, dy, dz
      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
  end subroutine
!===============================================
  subroutine AtomTrans_FullMove(self, trialBox) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(AtomMolTranslate), intent(inout) :: self
    type(SimBox), intent(inout) :: trialBox
    logical :: accept
    integer :: nMove
    integer :: CalcIndex
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
       !Choose Molecule
    nMove = floor( trialBox%nAtoms * grnd() + 1d0)
    dx = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dy = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dz = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
 
    self%disp(1)%atmIndx = nMove
    self%disp(1)%x_new = trialBox%atoms(1, nMove) + dx
    self%disp(1)%y_new = trialBox%atoms(2, nMove) + dy
    self%disp(1)%z_new = trialBox%atoms(3, nMove) + dz

    write(*,*) self%disp(1)%x_new, self%disp(1)%y_new, self%disp(1)%z_new
    write(*,*) trialBox%atoms(1, nMove), trialBox%atoms(2, nMove), trialBox%atoms(3, nMove)


    !Check Constraint
!    if(.not. accept) then
!      return
!    endif


    !Energy Calculation
    CalcIndex = trialBox % ECalcer
    call EnergyCalculator(CalcIndex) % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)
    write(*,*) E_Diff    

    accept = .true.
    !Accept/Reject
    if(accept) then
      self % accpt = self % atmps + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:1))
    endif

  end subroutine

!========================================================
end module
!========================================================
