!========================================================
!Guarnteed delete move.  Only used for debuging. 
!
module MCMove_Delete
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: MoveDelete
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
!    type(Displacement) :: disp(1:1)
    contains
      procedure, pass :: Constructor => MoveDelete_Constructor
      procedure, pass :: GeneratePosition => MoveDelete_GeneratePosition
      procedure, pass :: FullMove => MoveDelete_FullMove
      procedure, pass :: Maintenance => MoveDelete_Maintenance

  end type
!========================================================
 contains
!========================================================
  subroutine MoveDelete_Constructor(self)
    implicit none
    class(MoveDelete), intent(inout) :: self
    
  end subroutine
!========================================================
  subroutine MoveDelete_GeneratePosition(self, disp)
    use RandomGen, only: grnd
    implicit none
    class(MoveDelete), intent(in) :: self
    type(Displacement), intent(inout) :: disp
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
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex
    real(dp) :: E_Diff


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    E_Diff = 0E0_dp
    write(*,*) "Delete"

    !Accept/Reject
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(3)
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
