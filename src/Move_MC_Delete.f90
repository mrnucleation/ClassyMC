!========================================================
!Guarnteed delete move.  Only used for debuging. 
!
module MCMove_Delete
use CoordinateTypes, only: Deletion
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: MoveDelete
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    type(Deletion) :: disp(1:1)
    contains
      procedure, pass :: Constructor => MoveDelete_Constructor
!      procedure, pass :: GeneratePosition => MoveDelete_GeneratePosition
      procedure, pass :: FullMove => MoveDelete_FullMove
      procedure, pass :: Maintenance => MoveDelete_Maintenance

  end type
!========================================================
 contains
!========================================================
  subroutine MoveDelete_Constructor(self)
    implicit none
    class(MoveDelete), intent(inout) :: self

!    allocate( self%tempNNei(1) )
!    allocate( self%tempList(10, 1) )

  end subroutine
!===============================================
  subroutine MoveDelete_FullMove(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Box_Utility, only: FindAtom
    implicit none
    class(MoveDelete), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer, parameter :: delVal = 2
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex
    real(dp) :: E_Diff


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

!    self%disp(1)%newAtom = .false.
    self%disp(1)%MolType = 1
    self%disp(1)%MolIndx = delVal
!    self%disp(1)%atmIndx = delVal

!    self%disp(1)%oldAtom = .true.
!    self%disp(1)%oldMolType = 1
!    self%disp(1)%oldMolIndx = delVal
!    self%disp(1)%oldAtmIndx = delVal

!    self%disp(1)%newlist = .false.
!    self%disp(1)%listIndex = delVal

    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif


    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff, accept)

    write(*,*) "Delete"

    !Accept/Reject
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(delVal)
!      call trialBox % UpdatePosition(self%disp(1:1))
    endif

  end subroutine
!=========================================================================
  subroutine MoveDelete_Maintenance(self)
    implicit none
    class(MoveDelete), intent(inout) :: self
    real(dp), parameter :: limit = 3.0E0_dp
      
 

  end subroutine
!========================================================
end module
!========================================================
