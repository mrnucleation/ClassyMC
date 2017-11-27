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
    implicit none
    class(AtomExchange), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical :: accept
    integer :: oldtype, newtype


    self % atmps = self % atmps + 1E0_dp
    accept = .true.
    !Choose 
    oldtype = floor( trialBox% * grnd() + 1E0_dp)
    if(trialBox%NMolMin(oldtype) > trialBox%NMol(oldtype)-1) then
      return
    endif

    newtype = oldtype
    do while(newtype == oldtype)
      newtype = floor( trialBox% * grnd() + 1E0_dp)
    enddo

    nMove = floor( trialBox%nAtoms * grnd() + 1E0_dp)
    

    call trialbox% EFunc % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)
  end subroutine
!=========================================================================
  subroutine AtomExchange_Maintenance(self)
    class(AtomExchange), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
