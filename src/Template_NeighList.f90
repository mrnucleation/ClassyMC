module NeighListDef
use VarPrecision
use CoordinateTypes, only: Displacement
use Template_SimBox, only: SimBox

  type, public :: NeighList
      logical :: Strict = .false.
      integer, allocatable :: list(:,:)
      integer, allocatable :: nNeigh(:)
      integer :: maxNei
      real(dp) :: rCut, rCutSq
!      class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor
      procedure, pass :: InitializeList
      procedure, pass :: GetNewList
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine Constructor(self, parentID, rCut)
    implicit none
    class(NeighList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut

  end subroutine
!===================================================================================
  subroutine InitializeList(self)
    implicit none
    class(NeighList), intent(inout) :: self

  end subroutine
!===================================================================================
  subroutine GetNewList(self, disp, newList)
    implicit none
    class(NeighList), intent(inout) :: self
    type(Displacement), intent(in) :: disp
    real(dp), intent(out) :: newList(:)
  end subroutine
!===================================================================================
end module
!===================================================================================
