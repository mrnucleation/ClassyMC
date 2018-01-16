!===================================================================================
module Template_NeighList
use VarPrecision
use CoordinateTypes, only: Displacement
!use Template_SimBox, only: SimBox

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
      procedure, pass :: AddMol
      procedure, pass :: TransferList
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
  subroutine AddMol(self, disp, tempList, tempNNei)
    implicit none
    class(NeighListDef), intent(inout) :: self
    type(Displacement), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)


  end subroutine
!===================================================================================
  subroutine TransferList(self, indx1, indx2)
    use SearchSort, only: SimpleSearch, BinarySearch
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: indx1, indx2
    integer :: iNei, jAtom, bin, mNei, j

    self%nNeigh(indx2) = self%nNeigh(indx1)
    do iNei = 1, self%nNeigh(indx1)
      jAtom = self%list(iNei, indx1)
      self%list(iNei, indx2) = jAtom

      if(jAtom == indx2) then
        cycle
      endif
      mNei = self%nNeigh(jAtom)
      if(self%sorted) then
        bin = BinarySearch(indx1, self%list(1:mNei, jAtom))
        self%sorted = .false.
      else
        bin = SimpleSearch(indx1, self%list(1:mNei, jAtom))
      endif
      self%list(bin, jAtom) = indx2

    enddo

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
