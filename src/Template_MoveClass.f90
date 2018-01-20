!=========================================================================
module MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public :: MCMove
    real(dp) :: atmps = 1E-30_dp
    real(dp) :: accpt = 0E0_dp

    !Temporary Neighborlist Variables
    integer, allocatable :: tempNnei(:)
    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
      procedure, pass :: GetAcceptRate
      procedure, pass :: ProcessIO
      procedure, pass :: Maintenance 
  end type

 contains
!=========================================================================
  subroutine Constructor(self)
    implicit none
    class(MCMove), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine GeneratePosition(self, disp)
    use CoordinateTypes, only: Displacement
    implicit none
    class(MCMove), intent(in) :: self
    type(Displacement), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine FullMove(self, trialBox, accept)
    class(MCMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
  end subroutine
!=========================================================================
  function GetAcceptRate(self) result(rate)
    class(MCMove), intent(in) :: self
    real(dp) :: rate

    if(self%atmps > 0E0_dp) then
      rate = 1E2_dp*self%accpt/self%atmps
    else
      rate = 0E0_dp
    endif

    return
  end function
!=========================================================================
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(MCMove), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=========================================================================
  subroutine Maintenance(self)
    implicit none
    class(MCMove), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
