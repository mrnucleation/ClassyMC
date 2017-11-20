!=========================================================================
module MultiBoxMoveDef
use SimpleSimBox, only: SimpleBox
use VarPrecision
use MoveClassDef, only: MCMove

  type, public, extends(MCMove) :: MCMultiBoxMove
!    real(dp) :: atmps, accpt
    integer, allocatable :: allowedBox(:)
    contains
!      procedure, pass :: GeneratePosition 
!      procedure, pass :: FullMove
!      procedure, pass :: GetAcceptRate
      procedure, pass :: MultiBox
  end type

 contains
!=========================================================================
!  subroutine GeneratePosition(self, disp)
!    use CoordinateTypes, only: Displacement
!    implicit none
!    class(MCMove), intent(in) :: self
!    type(Displacement), intent(inout) :: disp
!  end subroutine
!=========================================================================
!  subroutine FullMove(self, trialBox)
!    class(MCMove), intent(inout) :: self
!    class(SimpleBox), intent(inout) :: trialBox
!  end subroutine
!=========================================================================
  subroutine MultiBox(self, trialBox)
    class(MCMultiBoxMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox(:)
  end subroutine
!=========================================================================
end module
!=========================================================================
