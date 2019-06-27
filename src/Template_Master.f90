!====================================================================
module MasterTemplate
  use VarPrecision

  type, public, abstract :: classyClass
    logical :: screenIO = .false.
    integer :: maintFreq = 10
    contains
       procedure, pass :: Destructor
       procedure, pass :: Epilogue
       procedure, pass :: Maintenance
       procedure, pass :: ModifyIO
       procedure, pass :: Prologue
       procedure, pass :: SafetyCheck
       procedure, pass :: ScreenOut
!       procedure, pass :: Report
       procedure, pass :: Update
 !      procedure, pass :: ProcessIO
  end type
!====================================================================
  contains
!====================================================================
! Called at the start of the simulation to ensure this object has
! everything it needs to properly function.  This can be used for example
! to force the user to define another object that this current object requires.
  subroutine SafetyCheck(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
! Called periodically during the simulation. Can be used for processes
! that need to be updated occationally, but not after every move/cycle.
  subroutine Maintenance(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
! Used to modify this object's parameters via the input script
  subroutine ModifyIO(self, line, lineStat)
    implicit none
    class(classyClass), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    lineStat = -1
    write(0,*) "This Object does not contain any modifiable parameters"
    write(0,*) line

  end subroutine
!====================================================================
! Used to periodically print information about this object to the screen mid simulation.
  subroutine ScreenOut(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
!  subroutine ProcessIO(self, line, lineStat)
!    implicit none
!    class(classyClass), intent(in) :: self
!    integer, intent(out) :: lineStat
!    character(len=*), intent(in) :: line   
!
!  end subroutine
!====================================================================
!  subroutine Report(self)
!    implicit none
!    class(classyClass), intent(inout) :: self
!
!  end subroutine
!====================================================================
  subroutine Destructor(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
! Used to update this object's information after a successful MC move.
  subroutine Update(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
! Function ran at the end of the simulation.  Can be used to print final
! information or perform end of simulation calculations
  subroutine Epilogue(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
! Function ran at the start of the simulation.  Can be used to print initial
! information or perform start of simulation calculations
  subroutine Prologue(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
end module
!====================================================================
