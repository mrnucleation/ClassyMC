!=========================================================================
module Anaylsis_RDF
use SimpleSimBox, only: SimpleBox
use Constants, only: pi
use VarPrecision

  type, public, extends(Analysis) :: rdf
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    integer :: boxNum
    integer :: type1, type2
    real(dp) :: dr = 0.01E0_dp
    real(dp), allocatable :: hist(:)
    real(dp), allocatable :: tempHist(:)
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
    implicit none
    class(Analysis), intent(in) :: self
  end subroutine
!=========================================================================
  subroutine Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(Analysis), intent(in) :: self
    logical, intent(in) :: accept
 
    integer :: iAtom, jAtom


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
