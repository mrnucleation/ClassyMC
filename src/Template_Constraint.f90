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
      procedure, pass :: PostEnergy
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
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    accept = .true.

  end subroutine
!=============================================================
  subroutine PostEnergy(self, trialBox, disp, E_Diff, accept)
    implicit none
    class(constraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff
    logical, intent(out) :: accept
    accept = .true.


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
end module
!=============================================================
