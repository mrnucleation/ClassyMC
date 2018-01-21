!===================================================================================
module Template_NeighList
  use MasterTemplate, only: classyClass
  use VarPrecision
  use CoordinateTypes, only: Displacement
!use Template_SimBox, only: SimBox

  type, public, extends(classyClass) :: NeighListDef
      logical :: Sorted = .false.
      logical :: Strict = .false.
      integer, allocatable :: list(:,:)
      integer, allocatable :: nNeigh(:)
      integer :: maxNei
      real(dp) :: rCut = -1E0_dp
      real(dp) :: rCutSq
      logical :: restrictType = .false.
      logical, allocatable :: allowed(:)

      
    contains
      procedure, pass :: Constructor
      procedure, pass :: BuildList
      procedure, pass :: GetNewList
      procedure, pass :: AddMol
      procedure, pass :: TransferList
      procedure, pass :: DeleteMol
      procedure, pass :: ProcessIO
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
      if(jAtom == indx2) then
        cycle
      endif
      self%list(iNei, indx2) = jAtom

      bin = self%nNeigh(jAtom)
      self%list(bin+1, jAtom) = indx2
      self%nNeigh(jAtom) = self%nNeigh(jAtom) + 1
    enddo

  end subroutine
!===================================================================================
  subroutine DeleteMol(self, molIndx, topIndx)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx


  end subroutine
!====================================================================
  subroutine ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    lineStat = 0

  end subroutine
!===================================================================================

end module
!===================================================================================
