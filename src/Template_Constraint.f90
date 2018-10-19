!=============================================================
module ConstraintTemplate
  use VarPrecision
  use MasterTemplate, only: classyClass
  use CoordinateTypes, only: Perturbation
  use Template_SimBox, only: SimBox

  type, public, extends(classyClass) :: constraint
    contains
      procedure, pass :: Constructor
      procedure, pass :: CheckInitialConstraint
      procedure, pass :: DiffCheck
      procedure, pass :: ProcessIO
!      procedure, pass :: Update
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
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    accept = .true.

    select type(disp)

    end select


  end subroutine
!=============================================================
  subroutine ProcessIO(self, line, lineStat)
    implicit none
    class(constraint), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    lineStat = 0
  end subroutine
!=============================================================
!  subroutine Update(self, trialBox, disp)
!    implicit none
!    class(constraint), intent(inout) :: self
!    class(SimBox), intent(in) :: trialBox
!    class(Perturbation), intent(in) :: disp(:)
!
!  end subroutine
!=============================================================
end module
!=============================================================
