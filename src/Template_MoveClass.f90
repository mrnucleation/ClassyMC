!=========================================================================
module MoveClassDef
  use MasterTemplate, only: classyClass
  use SimpleSimBox, only: SimpleBox
  use VarPrecision

  type, public, extends(classyClass) :: MCMove
    real(dp) :: atmps = 1E-30_dp
    real(dp) :: accpt = 0E0_dp
    real(dp), allocatable :: boxProb(:)

    !Temporary Neighborlist Variables
    integer, allocatable :: tempNnei(:)
    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
      procedure, pass :: GetAcceptRate
      procedure, pass :: GetBoxProb
      procedure, pass :: ProcessIO
!      procedure, pass :: Maintenance 
  end type

 contains
!=========================================================================
  subroutine Constructor(self)
    implicit none
    class(MCMove), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine GeneratePosition(self, disp)
    use CoordinateTypes, only: Perturbation
    implicit none
    class(MCMove), intent(in) :: self
    class(Perturbation), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine FullMove(self, trialBox, accept)
    class(MCMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    accept = .true.
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
  subroutine GetBoxProb(self, boxProb)
    implicit none
    class(MCMove), intent(inout) :: self
    real(dp), intent(inout) :: boxProb(:)

  end subroutine
!=========================================================================
  subroutine ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(MCMove), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat

    lineStat = 0
  end subroutine
!=========================================================================
!  subroutine Maintenance(self)
!    implicit none
!    class(MCMove), intent(inout) :: self
!  end subroutine
!=========================================================================
end module
!=========================================================================
