!=========================================================================
module Move_AtomExchange
use SimpleSimBox, only: SimpleBox
use VarPrecision
use MoveClassDef

  type, public, extends(MCMove) :: AtomExchange
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

    contains
      procedure, pass :: GeneratePosition => AtomExchange_GeneratePosition
      procedure, pass :: FullMove => AtomExchange_FullMove
!      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance => AtomExchange_Maintenance
  end type

 contains
!=========================================================================
  subroutine AtomExchange_GeneratePosition(self, disp)
    use CoordinateTypes, only: Displacement
    implicit none
    class(AtomExchange), intent(in) :: self
    type(Displacement), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine AtomExchange_FullMove(self, trialBox)
    use Common_MolDef, only: nAtomTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    implicit none
    class(AtomExchange), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical :: accept
    integer :: i
    integer :: nAtom, nAtomNew, reduIndx, newtype
    real(dp) :: OldProb, NewProb, Prob


    self % atmps = self % atmps + 1E0_dp
    accept = .true.
    ! Choose 
    reduIndx = floor( trialBox%nTotal * grnd() + 1E0_dp)
    call FindAtom(trialBox, reduIndx, nAtom)
    oldtype = trialBox%AtomType(nAtom)
    if(trialBox%NMolMin(oldtype) > trialBox%NMol(oldtype)-1) then
      return
    endif

    newtype = oldtype
    do while(newtype == oldtype)
      newtype = floor( trialBox% * grnd() + 1E0_dp)
    enddo
    if(trialBox%NMolMax(oldtype) < trialBox%NMol(oldtype)+1) then
      return
    endif
    call FindFirstEmptyMol(box, newtype, nAtomNew)

    self%disp(1)%newAtom = .true.
    self%disp(1)%MolType = newType
    self%disp(1)%MolIndx = nAtomNew
    self%disp(1)%atmIndx = nAtomNew
    self%disp(1)%x_new = trialBox%atoms(1, nAtom)
    self%disp(1)%y_new = trialBox%atoms(2, nAtom)
    self%disp(1)%z_new = trialBox%atoms(3, nAtom)

    self%disp(1)%oldAtom = .false.
    self%disp(1)%oldMolType = oldType
    self%disp(1)%oldMolIndx = nAtom
    self%disp(1)%oldAtmIndx = nAtom

    self%disp(1)%newlist = .false.
    self%disp(1)%listIndx = nAtom

    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif

    call trialbox % EFunc % Method % DiffECalc(trialBox, self%disp(1:1), E_Diff)

    NewProb = 1E0_dp / real(trialBox % NMol(newType) + 1, dp)
    OldProb = 1E0_dp / real(trialBox % NMol(oldType), dp)
    Prob = OldProb/NewProb

    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % AddMol( self%disp(1)%molType )
      call trialBox % DeleteMol( self%disp(1)%oldMolIndx )
      call trialBox % UpdatePosition( self%disp(1:1) )
      do i = 1, size(trialBox%NeighList)
        call trialBox % NeighList(i) % TransferList(nAtom, nAtomNew)
      enddo
    endif

  end subroutine
!=========================================================================
  subroutine AtomExchange_Maintenance(self)
    class(AtomExchange), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
