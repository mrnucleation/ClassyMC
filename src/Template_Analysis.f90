!=========================================================================
module AnaylsisClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public :: Anaylsis
    integer :: IOUnit = -1
    contains
      procedure, pass :: Initialize
      procedure, pass :: Compute
      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO
      procedure, pass :: WriteInfo
  end type

 contains
!=========================================================================
  subroutine Initialize(self)
    use CoordinateTypes, only: Displacement
    implicit none
    class(Analysis), intent(in) :: self
  end subroutine
!=========================================================================
  subroutine Compute(self, accept)
    use CoordinateTypes, only: Displacement
    implicit none
    class(Analysis), intent(in) :: self
    logical, intent(in) :: accept
  end subroutine
!=========================================================================
  subroutine Maintenance(self)
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine ProcessIO(self)
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine WriteInfo(self)
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
