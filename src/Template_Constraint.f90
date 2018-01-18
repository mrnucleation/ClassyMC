!=============================================================
module ConstraintTemplate
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use Template_SimBox, only: SimBox

  type, public :: constraint
    contains
      procedure, pass :: Constructor
      procedure, pass :: CheckInitialConstraint
      procedure, pass :: DiffCheck
      procedure, pass :: ShiftCheck
      procedure, pass :: NewCheck
      procedure, pass :: OldCheck
      procedure, pass :: ProcessIO
  end type

  type, public :: constrainArray
    class(constraint), allocatable :: method
  end type
!=============================================================
  contains
!=============================================================
  subroutine Constructor(self)
    implicit none
    class(constraint), intent(inout) :: self
  end subroutine
!=============================================================
  subroutine CheckInitialConstraint(self, accept)
    implicit none
    class(constraint), intent(in) :: self
    logical, intent(out) :: accept

    accept = .true.

  end subroutine
!=============================================================
  subroutine DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(in) :: self
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
    class(constraint), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.

  end subroutine
!=============================================================
  subroutine NewCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine OldCheck(self, trialBox, disp, accept)
    implicit none
    class(constraint), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine ProcessIO(self)
    implicit none
    class(constraint), intent(in) :: self
  end subroutine
!=============================================================
end module
!=============================================================
