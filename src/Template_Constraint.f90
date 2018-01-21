!=============================================================
module ConstraintTemplate
  use VarPrecision
  use MasterTemplate, only: classyClass
  use CoordinateTypes, only: Displacement
  use Template_SimBox, only: SimBox

  type, public, extends(classyClass) :: constraint
    contains
      procedure, pass :: Constructor
      procedure, pass :: CheckInitialConstraint
      procedure, pass :: DiffCheck
      procedure, pass :: ShiftCheck
      procedure, pass :: NewCheck
      procedure, pass :: OldCheck
      procedure, pass :: VolCheck
      procedure, pass :: ProcessIO
  end type

  type, public :: constrainArray
    class(constraint), allocatable :: method
  end type
!=============================================================
  contains
!=============================================================
  subroutine Constructor(self, boxID)
    implicit none
    class(constraint), intent(inout) :: self
    integer, intent(in) :: boxId

  end subroutine
!=============================================================
  subroutine CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    accept = .true.

  end subroutine
!=============================================================
  subroutine DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    accept = .true.

    !Called when a particle is either moved or replaced with another particle
    if(disp(1)%newAtom .and. disp(1)%oldAtom) then
      call self % ShiftCheck(trialBox, disp, accept)
      return
    endif

    !Called when a particle is added to a system
    if(disp(1)%newAtom) then
      call self % NewCheck(trialBox, disp, accept)
      return
    endif

    !Called when a particle is removed from a system.
    if(disp(1)%oldAtom) then
      call self % OldCheck(trialBox, disp, accept)
      return
    endif


  end subroutine
!=============================================================
  subroutine ShiftCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.

  end subroutine
!=============================================================
  subroutine NewCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine OldCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine VolCheck(self, trialBox, scalar, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    real(dp), intent(in) :: scalar(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine ProcessIO(self, line, lineStat)
    implicit none
    class(constraint), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

  end subroutine
!=============================================================
end module
!=============================================================
