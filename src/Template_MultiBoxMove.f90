!=========================================================================
module MultiBoxMoveDef
use SimpleSimBox, only: SimpleBox
use VarPrecision
use MoveClassDef, only: MCMove

  type, public, extends(MCMove) :: MCMultiBoxMove
!    real(dp) :: atmps, accpt
!    integer, allocatable :: allowedBox(:)
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
  subroutine MultiBox(self, accept)
    class(MCMultiBoxMove), intent(inout) :: self
!    class(SimpleBox), intent(inout) :: trialBox(:)
    logical, intent(out) :: accept

    write(0,*) "WARNING! MultiBox has not been defined for this move type!"
    stop
    accept = .true.
  end subroutine
!=========================================================================
end module
!=========================================================================
