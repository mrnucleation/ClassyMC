!===================================================================================
module Template_NeighList
use VarPrecision
use CoordinateTypes, only: Displacement
!use Template_SimBox, only: SimBox

  abstract interface
    subroutine buildfunc(trialBox)
     
    end subroutine
  end interface

  type, public :: NeighListDef
      logical :: Sorted = .false.
      logical :: Strict = .false.
      integer, allocatable :: list(:,:)
      integer, allocatable :: nNeigh(:)
      integer :: maxNei
      real(dp) :: rCut = -1E0_dp
      real(dp) :: rCutSq
!      class(SimBox), pointer :: parent => null()
      
    contains
      procedure, pass :: Constructor
      procedure, pass :: BuildList
      procedure, pass :: GetNewList
      procedure, pass :: DeleteMol
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine Constructor(self, parentID, rCut)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut

  end subroutine
!===================================================================================
  subroutine BuildList(self)
    implicit none
    class(NeighListDef), intent(inout) :: self

  end subroutine
!===================================================================================
  subroutine GetNewList(self, iDisp, tempList, tempNNei, disp)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: iDisp
    type(Displacement), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
 
  end subroutine
!===================================================================================
  subroutine DeleteMol(self, molIndx, topIndx)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx


  end subroutine
!===================================================================================
end module
!===================================================================================
