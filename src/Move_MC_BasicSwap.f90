!=========================================================================
module Move_BasicSwap
use SimpleSimBox, only: SimpleBox
use CoordinateTypes, only: Displacement
use VarPrecision
use MoveClassDef

  type, public, extends(MCMove) :: BasicSwap
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: inProb = 0.5E0_dp
    type(displacement), allocatable :: disp(:)
    contains
      procedure, pass :: Constructor => BasicSwap_Constructor
      procedure, pass :: GeneratePosition => BasicSwap_GeneratePosition
      procedure, pass :: FullMove => BasicSwap_FullMove
      procedure, pass :: Insertion => BasicSwap_Insertion
      procedure, pass :: Deletion => BasicSwap_Deletion
!      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance => BasicSwap_Maintenance
      procedure, pass :: Epilogue => BasicSwap_Epilogue
  end type

 contains
!========================================================
  subroutine BasicSwap_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(BasicSwap), intent(inout) :: self


    allocate( self%disp(1) )
    allocate( self%tempNNei(1) )
    allocate( self%tempList(100, 1) )
  end subroutine
!=========================================================================
  subroutine BasicSwap_GeneratePosition(self, disp)
    use CoordinateTypes, only: Displacement
    implicit none
    class(BasicSwap), intent(in) :: self
    type(Displacement), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine BasicSwap_FullMove(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin

    implicit none
    class(BasicSwap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    self % atmps = self % atmps + 1E0_dp
    if(grnd() < self%inProb) then
      call self%Insertion(trialBox, accept)
    else
      call self%Deletion(trialBox, accept)
    endif

  end subroutine
!=========================================================================
  subroutine BasicSwap_Insertion(self, trialBox, accept)
    implicit none
    class(BasicSwap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    integer :: i
    integer :: nAtom, nAtomNew, reduIndx, newtype, oldtype
    real(dp) :: OldProb, NewProb, Prob
    real(dp) :: E_Diff


    newType = floor(grnd()*nMolTypes +1E0_dp)
    if(trialBox%NMolMax(newtype) < trialBox%NMol(newtype)+1) then
      return
    endif

    nAtomNew = 1
    do i = 1, newtype - 1
      nAtomNew = nAtomNew + trialBox%NMolMax(i)
    enddo
    nAtomNew = nAtomNew + trialBox%NMol(newtype)

    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif

    call trialbox % EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif

    NewProb = 1E0_dp / trialBox % Volume
    OldProb = 1E0_dp / real(trialBox % NMol(newType)+1, dp)
    Prob = OldProb/NewProb

    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % AddMol( self%disp(1)%molType )
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList(:,:), self%tempNNei(:))
     endif
  end subroutine
!=========================================================================
  subroutine BasicSwap_Deletion(self, trialBox, accept)
    implicit none
    class(BasicSwap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i
    integer :: nAtom, nAtomNew, reduIndx, newtype, oldtype
    real(dp) :: OldProb, NewProb, Prob
    real(dp) :: E_Diff


    accept = .true.

    reduIndx = floor( trialBox%nTotal * grnd() + 1E0_dp)
    call FindAtom(trialBox, reduIndx, nAtom)
    oldtype = trialBox%AtomType(nAtom)
    if(trialBox%NMolMin(oldtype) > trialBox%NMol(oldtype)-1) then
      return
    endif

    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif

    call trialbox % EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif


    NewProb = 1E0_dp / real(trialBox % NMol(newType), dp)
    OldProb = 1E0_dp / trialBox % Volume
    Prob = OldProb/NewProb

    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol( self%disp(1)%oldMolIndx )
    endif

  end subroutine
!=========================================================================
  subroutine BasicSwap_Maintenance(self)
    class(BasicSwap), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine BasicSwap_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(BasicSwap), intent(inout) :: self
    real(dp) :: accptRate
      
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,I15)") "Basic Swap Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Basic Swap Moves Attempted: ", nint(self%atmps)
    write(nout,"(1x,A,F15.8)") "Basic Swap Acceptance Rate: ", accptRate

  end subroutine
!=========================================================================
end module
!=========================================================================
