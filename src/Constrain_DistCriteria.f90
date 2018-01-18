!====================================================================
module Constrain_DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Displacement
  use Template_SimBox, only: SimBox

  type, public, extends(constraint) :: DistCriteria
    integer :: neighList = -1
    integer :: molType = -1
    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: CheckInitialConstraint => DistCrit_CheckInitialConstraint
!      procedure, pass :: DiffCheck
      procedure, pass :: ShiftCheck => DistCrit_ShiftCheck
      procedure, pass :: NewCheck => DistCrit_NewCheck
      procedure, pass :: OldCheck => DistCrit_OldCheck
!      procedure, pass :: ProcessIO => DistCrit_ProcessIO
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine DistCrit_CheckInitialConstraint(self, accept)
    implicit none
    class(DistCriteria), intent(in) :: self
    logical, intent(out) :: accept

  end subroutine
!=====================================================================
  subroutine DistCrit_ShiftCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    
    integer :: indx2
    real(dp) :: rx, ry, rz, rsq
    

    accept = .true.

  end subroutine
!=====================================================================
  subroutine DistCrit_NewCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=====================================================================
  subroutine DistCrit_OldCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    implicit none
    class(DistCriteria), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % neighList = intVal

      case("moleculetype")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % moltype = intVal

      case("rcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % rCut = realVal
        self % rCutSq = realVal**2

      case default
        lineStat = -1
    end select

  end subroutine
!=====================================================================
end module
!=====================================================================
