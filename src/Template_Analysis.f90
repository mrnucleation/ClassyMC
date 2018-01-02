!=========================================================================
module AnaylsisClassDef
use VarPrecision

  type, public :: Analysis
    logical :: perMove = .false.
    integer :: IOUnit = -1
    integer :: UpdateFreq = -1
    contains
      procedure, pass :: Initialize
      procedure, pass :: Compute
      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO
      procedure, pass :: WriteInfo
      procedure, pass :: GetResult
      procedure, pass :: Finalize
  end type

 contains
!=========================================================================
  subroutine Initialize(self)
    implicit none
    class(Analysis), intent(in) :: self
  end subroutine
!=========================================================================
  subroutine Compute(self, accept)
    implicit none
    class(Analysis), intent(inout) :: self
    logical, intent(in) :: accept
  end subroutine
!=========================================================================
  subroutine Maintenance(self)
    implicit none
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Analysis), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
  end subroutine
!=========================================================================
  subroutine WriteInfo(self)
    implicit none
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
  function GetResult(self) result(var)
    implicit none
    class(Analysis), intent(in) :: self
    real(dp) :: var
  end function
!=========================================================================
  subroutine Finalize(self)
    implicit none
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
