!===================================================================================
module Template_NeighList
  use MasterTemplate, only: classyClass
  use VarPrecision
  use CoordinateTypes, only: Perturbation, AtomExchange
!use Template_SimBox, only: SimBox

  type, public, extends(classyClass) :: NeighListDef
      logical :: Sorted = .false.
      logical :: Strict = .false.
      integer, allocatable :: list(:,:)
      integer, allocatable :: nNeigh(:)
      integer, allocatable :: templist(:,:)
      integer, allocatable :: tempNNeigh(:)
      integer :: maxNei
      real(dp) :: rCut = -1E0_dp
      real(dp) :: rCutSq
      logical :: restrictType = .false.
      logical, allocatable :: allowed(:)

      
    contains
      procedure, pass :: Constructor
      procedure, pass :: BuildList
      procedure, pass :: SortList
      procedure, pass :: GetListArray
      procedure, pass :: GetTempListBounds
      procedure, pass :: GetNewList
      procedure, pass :: GetMaxNei
      procedure, pass :: GetNeighCount
      procedure, pass :: GetRCut
      procedure, pass :: AddMol
      procedure, pass :: SwapAtomType
      procedure, pass :: DumpList
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
  subroutine BuildList(self, listindx)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: listindx

  end subroutine
!===================================================================================
  subroutine SortList(self, forcesort)
    use SearchSort, only: QSort
    implicit none
    class(NeighListDef), intent(inout) :: self
    logical, intent(in), optional :: forcesort


  end subroutine
!===================================================================================
  subroutine GetListArray(self,list,  nNeigh )
    implicit none
    class(NeighListDef), intent(inout), target :: self
    integer, pointer, intent(inout) :: nNeigh(:)
    integer, pointer, intent(inout) :: list(:,:)

    nNeigh => self%nNeigh
    list => self%list
 
  end subroutine
!===================================================================================
  subroutine GetTempListBounds(self, b1, b2)
    implicit none
    class(NeighListDef), intent(inout), target :: self
    integer, intent(out) :: b1, b2

    b1 = ubound(self%list, 1)
    b2 = ubound(self%list, 2)

  end subroutine
!===================================================================================
  function GetMaxNei(self) result(outval)
    implicit none
    class(NeighListDef), intent(inout), target :: self
    integer :: outval

    outval = self%maxnei
  end function
!===================================================================================
  function GetRCut(self) result(outval)
    implicit none
    class(NeighListDef), intent(inout), target :: self
    integer :: outval

    outval = self%rcut
  end function
!===================================================================================
  subroutine GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: iDisp
!    type(Displacement), intent(inout) :: disp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount


 
  end subroutine
!===================================================================================
  function GetNeighCount(self, nAtom, rCount) result(nCount)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: nAtom
    real(dp), intent(in) ,optional :: rCount
    integer :: nCount
 
  end function 
!===================================================================================
  subroutine AddMol(self, disp, tempList, tempNNei)
    implicit none
    class(NeighListDef), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)


  end subroutine
!===================================================================================
  subroutine SwapAtomType(self, disp, topIndx)
    implicit none
    class(NeighListDef), intent(inout) :: self
    class(AtomExchange), intent(in) :: disp(:)
    integer, intent(in) :: topIndx
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
!===================================================================================
  subroutine DumpList(self, filenum)
    implicit none
    class(NeighListDef), intent(inout) :: self
    integer, intent(in) :: filenum
    integer :: uB1, uB2, lB1, lB2
    integer :: iNei, iAtom

    lB1 = lbound(self%list, 1)
    lB2 = lbound(self%list, 2)

    uB1 = ubound(self%list, 1)
    uB2 = ubound(self%list, 2)

    write(filenum,*) lB1, lB2, uB1, uB2
    do iAtom = lB2, uB2
      write(filenum, *) iAtom, "|", (self%list(iNei, iAtom), iNei = lB1, uB1)
    enddo
    write(filenum, *) 


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
