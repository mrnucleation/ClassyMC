!====================================================================
module MasterTemplate
  use VarPrecision

  type, public :: classyClass
    integer :: maintFreq = 100
    contains
       procedure, pass :: Epilogue
       procedure, pass :: SafetyCheck
       procedure, pass :: Maintenance
!       procedure, pass :: ProcessIO
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
!  subroutine ProcessIO(self, line, lineStat)
!    implicit none
!    class(classyClass), intent(in) :: self
!    integer, intent(out) :: lineStat
!    character(len=*), intent(in) :: line   
!
!  end subroutine
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
