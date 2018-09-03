!====================================================================
module MasterTemplate
  use VarPrecision

  type, public :: classyClass
    integer :: maintFreq = 10
    contains
       procedure, pass :: Epilogue
       procedure, pass :: SafetyCheck
       procedure, pass :: GetScreenData
       procedure, pass :: Maintenance
       procedure, pass :: ModifyIO
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

    lineStat = 0

  end subroutine
!====================================================================
  subroutine GetScreenData(self, val)
    implicit none
    class(classyClass), intent(inout) :: self
    class(*), intent(out) :: val

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
