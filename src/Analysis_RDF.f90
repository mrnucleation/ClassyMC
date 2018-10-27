!=========================================================================
module Analysis_RDF

use AnaylsisClassDef, only: Analysis
use SimpleSimBox, only: SimpleBox
use ClassyConstants, only: pi
use VarPrecision

  type, public, extends(Analysis) :: rdf
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    integer :: boxNum
    integer :: type1, type2
    integer :: bins = 1000
    real(dp) :: dr = 0.01E0_dp
    real(dp), allocatable :: hist(:)
    real(dp), allocatable :: tempHist(:)
    contains
      procedure, pass :: Initialize => RDF_Initialize
      procedure, pass :: Compute => RDF_Compute
      procedure, pass :: ProcessIO => RDF_ProcessIO
      procedure, pass :: WriteInfo => RDF_WriteInfo
  end type

 contains
!=========================================================================
  subroutine RDF_Initialize(self)
    implicit none
    class(rdf), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine RDF_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(rdf), intent(inout) :: self
    logical, intent(in) :: accept
 
    integer :: iAtom, jAtom


  end subroutine
!=========================================================================
  subroutine RDF_ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(rdf), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
  end subroutine
!=========================================================================
  subroutine RDF_WriteInfo(self)
    class(rdf), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
