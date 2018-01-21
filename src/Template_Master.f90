!====================================================================
module MasterTemplate
  use VarPrecision

  type, public :: classyClass
    integer :: maintFreq
    contains
       procedure, pass :: Maintenance
!       procedure, pass :: ProcessIO
  end type
!====================================================================
  contains
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
end module
!====================================================================
