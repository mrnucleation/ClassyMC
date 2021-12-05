!=========================================================================
module Move_AtomExchange
use CoordinateTypes, only: Perturbation, AtomExchange
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: MC_AtomExchange
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    type(AtomExchange) :: disp(1:1)
    contains
      procedure, pass :: Constructor => AtomExchange_Constructor
      procedure, pass :: GeneratePosition => AtomExchange_GeneratePosition
      procedure, pass :: FullMove => AtomExchange_FullMove
!      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance => AtomExchange_Maintenance
      procedure, pass :: Epilogue => AtomExchange_Epilogue
  end type

 contains
!========================================================
  subroutine AtomExchange_Constructor(self)
!    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(MC_AtomExchange), intent(inout) :: self


  end subroutine
!=========================================================================
  subroutine AtomExchange_GeneratePosition(self, disp)
    implicit none
    class(MC_AtomExchange), intent(in) :: self
    class(Perturbation), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine AtomExchange_FullMove(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindMolecule, FindFirstEmptyMol
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
!    use Common_NeighData, only: neighSkin

    implicit none
    class(MC_AtomExchange), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i
    integer :: nAtom, nAtomNew, reduIndx, newtype, oldtype
    integer :: nTarget
    real(dp) :: Prob
    real(dp) :: E_Diff, E_Inter, E_Intra
    real(dp) :: extraTerms


    accept = .true.
    call self%LoadBoxInfo(trialBox, self%disp)
    ! Choose an atom to remove
    reduIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, reduIndx, nTarget)
    call trialBox % GetMolData(nTarget, molStart=nAtom)

!    write(*,*) nAtom
    oldtype = trialBox%MolType(nTarget)
    if(trialBox%NMolMin(oldtype) > trialBox%NMol(oldtype)-1) then
      return
    endif

    newtype = oldtype
    do while(newtype == oldtype)
      newtype = floor( nMolTypes * grnd() + 1E0_dp)
    enddo

    if(newtype == oldtype) then
      return
    endif

    if(trialBox%NMolMax(newtype) < trialBox%NMol(newtype)+1) then
      return
    endif

    self % atmps = self % atmps + 1E0_dp
    nAtomNew = 1
    do i = 1, newtype - 1
      nAtomNew = nAtomNew + trialBox%NMolMax(i)
    enddo
    nAtomNew = nAtomNew + trialBox%NMol(newtype)

    self%disp(1)%newAtmIndx = trialbox%MolStartIndx(nAtomNew)
    self%disp(1)%newType = newType
    self%disp(1)%oldAtmIndx = nAtom
    self%disp(1)%oldType = oldType

    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif

!    call trialbox % EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    call trialBox%ComputeEnergyDelta(self%disp(1:1),&
                                     self%templist, &
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.false.)
    if(.not. accept) then
      return
    endif
    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:1), E_Diff )
    if(.not. accept) then
      return
    endif


!    NewProb = 1E0_dp / real(trialBox % NMol(newType) + 1, dp)
!    OldProb = 1E0_dp / real(trialBox % NMol(oldType), dp)
!    Prob = OldProb/NewProb
    Prob = 1.0E0_dp
    extraTerms = sampling % GetExtraTerms(self%disp(1:1), trialBox)
!    write(*,*) E_Diff, extraTerms

    accept = sampling % MakeDecision(trialBox, E_Diff,  self%disp(1:1), inProb=Prob, extraIn=extraTerms)
!    write(*,*) accept
    if(accept) then
!      write(*,*) "Atom Exchange"
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter)
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList(:,:), self%tempNNei(:))
     endif

  end subroutine
!=========================================================================
  subroutine AtomExchange_Maintenance(self)
    class(MC_AtomExchange), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine AtomExchange_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(MC_AtomExchange), intent(inout) :: self
    real(dp) :: accptRate
      
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,I15)") "Atom Exchange Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Atom Exchange Moves Attempted: ", nint(self%atmps)
    write(nout,"(1x,A,F15.8)") "Atom Exchange Acceptance Rate: ", accptRate

  end subroutine
!=========================================================================
end module
!=========================================================================
