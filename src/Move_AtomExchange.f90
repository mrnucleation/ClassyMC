!=========================================================================
module Move_AtomExchange
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public :: AtomExchange
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
    use Box_Utility, only: FindAtom
    implicit none
    class(AtomExchange), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical :: accept
    integer :: nAtom, reduIndx, newtype


    self % atmps = self % atmps + 1E0_dp
    accept = .true.
    !Choose 
    reduIndx = floor( trialBox%nTotal * grnd() + 1E0_dp)
    if(trialBox%NMolMin(oldtype) > trialBox%NMol(oldtype)-1) then
      return
    endif
    call FindAtom(trialBox, reduIndx, nAtom)
    oldtype = trialBox%AtomType(nAtom)

    newtype = oldtype
    do while(newtype == oldtype)
      newtype = floor( trialBox% * grnd() + 1E0_dp)
    enddo

    if( size(Constrain) > 0 ) then
      do iConstrain = 1, size(Constrain)
        call Constrain(iConstrain) % method % ShiftCheck( trialBox, self%disp(1:1), accept )
      enddo
      if(.not. accept) then
        return
      endif
    endif    

    call trialbox % EFunc % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)
  end subroutine
!=========================================================================
  subroutine AtomExchange_Maintenance(self)
    class(AtomExchange), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
