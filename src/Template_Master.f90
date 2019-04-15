!====================================================================
module MasterTemplate
  use VarPrecision

  type, public, abstract :: classyClass
    logical :: screenIO = .false.
    integer :: maintFreq = 10
    contains
       procedure, pass :: Epilogue
       procedure, pass :: SafetyCheck
       procedure, pass :: ScreenOut
       procedure, pass :: Maintenance
       procedure, pass :: ModifyIO
!       procedure, pass :: Report
       procedure, pass :: Update
 !      procedure, pass :: ProcessIO
       procedure, pass :: Prologue
  end type
!====================================================================
  contains
!====================================================================
  subroutine SafetyCheck(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
  subroutine Maintenance(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
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
  subroutine Update(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
  subroutine Epilogue(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
  subroutine Prologue(self)
    implicit none
    class(classyClass), intent(inout) :: self

  end subroutine
!====================================================================
end module
!====================================================================
