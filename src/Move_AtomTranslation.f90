!========================================================
module MCMove_AtomTranslation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: AtomTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: max_dist = 0.05E0_dp
    type(Displacement) :: disp(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => AtomTrans_Constructor
      procedure, pass :: GeneratePosition => AtomTrans_GeneratePosition
      procedure, pass :: FullMove => AtomTrans_FullMove
      procedure, pass :: Maintenance => AtomTrans_Maintenance

  end type
!========================================================
 contains
!========================================================
  subroutine AtomTrans_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(AtomTranslate), intent(inout) :: self
!    integer :: iType, maxAtoms



    allocate( self%tempNNei(1) )
    allocate( self%tempList(1000, 1) )
  end subroutine
!========================================================
  subroutine AtomTrans_GeneratePosition(self, disp)
    use RandomGen, only: grnd
    implicit none
    class(AtomTranslate), intent(in) :: self
    type(Displacement), intent(inout) :: disp
    real(dp) :: dx, dy, dz
      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
  end subroutine
!===============================================
  subroutine AtomTrans_FullMove(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Box_Utility, only: FindAtom
    implicit none
    class(AtomTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp)
    call FindAtom(trialbox, rawIndx, nMove)
!    write(*,*) nMove, rawIndx
    dx = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dy = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dz = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
 
    self%disp(1)%newatom = .true.
    self%disp(1)%oldatom = .true.
    self%disp(1)%molType = trialBox%MolType(nMove)
    self%disp(1)%molIndx = trialBox%MolIndx(nMove)
    self%disp(1)%atmIndx = nMove

    self%disp(1)%oldMolType = trialBox%MolType(nMove)
    self%disp(1)%oldMolIndx = trialBox%MolIndx(nMove)
    self%disp(1)%oldAtmIndx = nMove

    self%disp(1)%x_new = trialBox%atoms(1, nMove) + dx
    self%disp(1)%y_new = trialBox%atoms(2, nMove) + dy
    self%disp(1)%z_new = trialBox%atoms(3, nMove) + dz

    !If the particle moved a large distance get a temporary neighborlist
    if(any([dx,dy,dz] > neighSkin)) then
      call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%disp(1))
      self%disp(1)%newlist = .true.
    else
      self%disp(1)%newlist = .false.
      self%disp(1)%listIndex = nMove
    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialbox% EFunc % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff)
    write(*,*) E_Diff

    !Accept/Reject
    accept = .false.
    accept = sampling % MakeDecision(trialBox, E_Diff, 1E0_dp, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine AtomTrans_Maintenance(self)
    implicit none
    class(AtomTranslate), intent(inout) :: self
    real(dp), parameter :: limit = 3.0E0_dp
      
    if(self%atmps .lt. 0.5E0_dp) then
      return
    endif

    if(self%GetAcceptRate() .gt. 50E0_dp) then
      if(self%max_dist*1.01E0_dp .lt. limit) then
        self%max_dist = self%max_dist * 1.01E0_dp
      else 
        self%max_dist = limit       
      endif
    else
      self%max_dist = self%max_dist * 0.99E0_dp
    endif

 

  end subroutine
!========================================================
end module
!========================================================
