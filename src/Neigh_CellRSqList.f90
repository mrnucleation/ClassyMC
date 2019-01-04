!===================================================================================
! This module contains a simple neighborlist
!===================================================================================
module CellRSqListDef
use VarPrecision
use CoordinateTypes
use Template_SimBox, only: SimBox
use SimpleSimBox, only: SimpleBox
use Template_NeighList, only: NeighListDef

  type, public, extends(NeighListDef) :: CellRSqList
!      logical :: Sorted = .false.
!      logical :: Strict = .false.
!      integer, allocatable :: list(:,:)
!      integer, allocatable :: nNeigh(:)
!      integer :: maxNei
!      real(dp) :: rCut, rCutSq
!      logical :: restrictType = .false.
!      integer, allocatable :: allowed(:)
!      integer :: safetyCheck = .false.
      integer, allocatable :: cellID(:)
      integer, allocatable :: cellList(:)

      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => CellRSqList_Constructor 
      procedure, pass :: BuildList => CellRSqList_BuildList 
      procedure, pass :: GetNewList => CellRSqList_GetNewList
      procedure, pass :: AddMol => CellRSqList_AddMol
      procedure, pass :: GetNeighCount => CellRSqList_GetNeighCount
      procedure, pass :: ProcessIO => CellRSqList_ProcessIO
!      procedure, pass :: TransferList
      procedure, pass :: DeleteMol => CellRSqList_DeleteMol
      procedure, pass :: Prologue => CellRSqList_Prologue
      procedure, pass :: Update => CellRSqList_Update
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine CellRSqList_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 0.65E0_dp  !Used to estimate an approximate volume of 
    integer :: AllocateStatus

    self%parent => BoxArray(parentID)%box
    if(.not. allocated(self%parent%atoms) ) then
      stop
    endif

!     If no rCut value is given by the subroutine call attempt to pull
!     the rSq value from the parent box's energy function. The assumption being
!     that the neighborlist is used for the energy calculation routines.
    if( present(rCut) ) then
      self % rCut = rCut
      self % rCutSq = rCut * rCut
      self % maxNei = ceiling(rCut**3/atomRadius**3)
    else
      if(self%rCut > 0E0_dp) then
        self % rCutSq = (self%rCut)**2
        self % maxNei = ceiling(self%rCut**3/atomRadius**3)

      else
        self % rCut = self % parent % EFunc % Method % GetCutOff() + neighSkin
        self % rCutSq = (self%rCut)**2
        self % maxNei = ceiling(self%rCut**3/atomRadius**3)
      endif
    endif
 
    write(nout,*) "Neighbor List CutOff:", self%rCut
    if(self%maxNei > self%parent%nMaxAtoms) then
      self%maxNei = self%parent%nMaxAtoms
    endif
    write(nout,*) "Neighbor List Maximum Neighbors:", self%maxNei

    allocate( self%cellID(1:self%parent%nMaxAtoms), stat=AllocateStatus )
    allocate( self%list(1:self%maxNei, 1:self%parent%nMaxAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:self%parent%nMaxAtoms), stat=AllocateStatus )

    if(.not. allocated(self%allowed) ) then
      allocate(self%allowed(1:nAtomTypes), stat=AllocateStatus )
      self%allowed = .true.
    endif

    self%list = 0
    self%nNeigh = 0 
    IF (AllocateStatus /= 0) STOP "*** NeighRSQList: Not enough memory ***"

    self%restrictType = .false.
  end subroutine
!===================================================================================
  subroutine CellRSqList_Prologue(self)
    implicit none
    class(CellRSqList), intent(inout) :: self


!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine CellRSqList_Update(self)
    implicit none
    class(CellRSqList), intent(inout) :: self


!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine CellRSqList_BuildList(self)
    implicit none
    class(CellRSqList), intent(inout) :: self

    real(dp) :: boxdim(1:9)


    self%parent%GetDimensions(boxdim)


!    call Builder_RSq(self%parent)
  end subroutine
!===================================================================================
  function CellRSqList_GetNeighCount(self, nAtom, rCount) result(nCount)
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: nAtom
    real(dp), intent(in) ,optional :: rCount
    integer :: iNei
    real(dp) :: rCut, rCutSq
    real(dp) :: rx, ry, rz, rsq
    integer :: nCount, jAtom

    if(present(rCount)) then
      rCut = rCount
      rCutSq = rCount*rCount
    else
      rCut = self%rCut
      rCutSq = self%rCutSq
    endif

    nCount = 0
    do iNei = 1, self % nNeigh(nAtom)
      jAtom = self%list(iNei, nAtom) 
      rx = self%parent%atoms(1, nAtom) - self%parent%atoms(1, jAtom)
      ry = self%parent%atoms(2, nAtom) - self%parent%atoms(2, jAtom)
      rz = self%parent%atoms(3, nAtom) - self%parent%atoms(3, jAtom)
      call self%parent%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < rCutSq) then
        nCount = nCount + 1
      endif
    enddo

  end function
!===================================================================================
  subroutine CellRSqList_AddMol(self, disp, tempList, tempNNei)
    implicit none
    class(CellRSqList), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)

    select type(disp)

      class is(Addition)
        call UpdateList_AddMol_RSq(self%parent, disp, tempList, tempNNei)
    end select

  end subroutine
!===================================================================================
  subroutine CellRSqList_DeleteMol(self, molIndx, topIndx)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: BinarySearch, SimpleSearch
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx
    integer :: iAtom, iNei, jNei, nType, j
    integer :: nStart, topStart
    integer :: nEnd, topEnd
    integer :: atmIndx, topAtom
    integer :: curNei, curIndx, nNei

    nStart = self % parent % MolStartIndx(molIndx)
    topStart = self % parent % MolStartIndx(topIndx)

    nEnd = self % parent % MolEndIndx(molIndx)
    topEnd = self % parent % MolEndIndx(topIndx)
    nType = self % parent % MolType(nStart)

!    write(2,*) "----------------------------"
!    write(2,*) "Delete"
!    write(2,*) "Removed Mol:", molIndx
!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "Sorted?:", self%sorted

    do iAtom = 1, MolData(nType)%nAtoms
!      atmIndx = nStart + iAtom - 1
!      topAtom = topStart + iAtom - 1


      atmIndx = nEnd - iAtom + 1
      topAtom = topEnd - iAtom + 1
!      write(2,*) "iAtom", atmIndx, "topAtom", topAtom
      !Remove the deleted from the list of it's neighbors
      do iNei = 1, self % nNeigh(atmIndx)
        curNei = self % list(iNei, atmIndx)
        nNei = self%nNeigh(curNei)
        if(nNei == 0) then
          cycle
        endif

        if(self%sorted) then
          curIndx = BinarySearch( atmIndx, self%list(1:nNei, curNei) )
        else
          curIndx = SimpleSearch( atmIndx, self%list(1:nNei, curNei) )
        endif

!        if(nNei <= 1) then
!          self%nNeigh(curNei) = 0
!          self%list(:, curNei) = 0
!          curIndx = 0
!          cycle
!        else
!
!        endif
!        write(2,*) "curNei", curNei, "Indexes", curIndx, atmIndx
        if(curIndx /= 0) then
          if(nNei > 2) then
            self%list(1:nNei, curNei ) = [self%list(1:curIndx-1, curNei), &
                                          self%list(curIndx+1:nNei, curNei) ]
          else
            if(curIndx == 1) then
              self%list(1, curNei) = self%list(2,curNei)
            endif
          endif
          self%nNeigh(curNei) = self%nNeigh(curNei) - 1 

        endif
      enddo
    enddo
!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo

    do iAtom = 1, MolData(nType)%nAtoms
      atmIndx = nEnd - iAtom + 1
      topAtom = topEnd - iAtom + 1     
      !Move the top atom into the spot formerly taken up by the deleted atom.
      do iNei = 1, self % nNeigh(topAtom)
        self%list(iNei, atmIndx) = self%list(iNei, topAtom)
      enddo
      self%nNeigh(atmIndx) = self%nNeigh(topAtom)

      !Re-index the neighbor's of the top atom to it's new array location
      do iNei = 1, self % nNeigh(topAtom)
        curNei = self % list(iNei, topAtom)
        nNei = self%nNeigh(curNei)
        if(nNei == 0 ) then
          cycle
        endif
        if(self%sorted) then
          curIndx = BinarySearch( topAtom, self%list(1:nNei, curNei) )
          self%sorted = .false.
        else
          curIndx = SimpleSearch( topAtom, self%list(1:nNei, curNei) )
        endif
!        write(2,*) "Second", "curNei", curNei, "Indexes", curIndx, topAtom
        if(curIndx /= 0) then
          self % list(curIndx, curNei) = atmIndx
        endif
      enddo
      self%nNeigh(topAtom) = 0
!      do iNei = 1, self%parent%nMaxatoms
!        write(2,"(I3,A,1000(I3))") iNei,"|", (self%list(j, iNei) ,j=1,self%nNeigh(iNei))
!      enddo
    enddo

!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "---------------------------------"


    self % sorted = .false.

  end subroutine
!===================================================================================
  subroutine CellRSqList_GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: iDisp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount
    integer :: jType, jAtom, j, iAtom
    integer :: jUp, jLow, molIndx, jMol
    real(dp) :: xn, yn, zn
    real(dp) :: rx, ry, rz, rsq

    if(present(nCount)) then
      nCount = 0
    endif
    select type(disp)
      class is (Addition)
!        disp % newlist = .true.
        disp % listIndex = iDisp
        iAtom = disp%atmIndx
        molIndx = self%parent%MolIndx(iAtom)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
      class is (Displacement)
        disp % newlist = .true.
        disp % listIndex = iDisp
        iAtom = disp%atmIndx
        molIndx = self%parent%MolIndx(iAtom)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
    end select

    templist(:, iDisp) = 0
    tempNNei(iDisp) = 0

!    molStart = 1
!    do jType = 1, nMolTypes
      do jAtom = 1, self%parent%nMaxAtoms
        if( self%parent%MolSubIndx(jAtom) == molIndx ) then
          cycle
        endif
        if( self%parent%MolSubIndx(jAtom) > self%parent%NMol(self%parent%MolType(jAtom)) ) then
          cycle
        endif
        rx = xn - self%parent%atoms(1, jAtom)
        ry = yn - self%parent%atoms(2, jAtom)
        rz = zn - self%parent%atoms(3, jAtom)
        call self%parent%Boundary(rx,ry,rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          tempNNei(iDisp) = tempNNei(iDisp) + 1
          templist(tempNNei(iDisp), iDisp) = jAtom
        endif
        if(present(rCount)) then
          if(rsq < rCount*rCount) then
            nCount = nCount + 1
          endif
        endif
      enddo
!      molStart = molStart + self%parent%NMolMax(jType)
!    enddo
!    write(2,"(A, 1000(I3))") "New",   templist(1:tempNNei(iDisp), iDisp)
  end subroutine
!====================================================================
  subroutine CellRSqList_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetAllCommands, GetXCommand,maxLineLen
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: i, intVal, nPar
    real(dp) :: realVal

    character(len=30) :: command 
    character(len=30), allocatable :: parlist(:)


    lineStat = 0
    call GetXCommand(line, command, 6, lineStat)
    select case( trim(adjustl(command)) )
      case("rcut")
        call GetXCommand(line, command, 7, lineStat)
        read(command,*) realVal
        self%rCut = realVal
        self%rCutSq = realVal * realVal

      case("restricttype")
        call GetAllCommands(line, parlist, nPar, lineStat)
        self%restrictType = .true.
        if(.not. allocated(self%allowed) ) then
          allocate(self%allowed(1:nAtomTypes) )
        endif
        self%allowed = .false.
        do i = 7, size(parList)
          read(parList(i), *) intVal
          self%allowed(intVal) = .true.
        enddo

      case default
        lineStat = -1
    end select

  end subroutine
!===================================================================================
! End Type Bound
!===================================================================================
  subroutine UpdateList_AddMol_RSq(trialBox, disp, tempList, tempNNei)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: trialBox
    class(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iList, iDisp, iAtom, iNei, nNei, neiIndx, j
    real(dp) :: rx, ry, rz, rsq




    do iList = 1, size(trialBox%NeighList)
      if(iList == 1) then
!        write(2,*) "----------------------------"
!        write(2,*) "Add"
!        do iAtom = 1, trialBox%nMaxatoms
!          write(2,"(I3,A,1000(I3))") iAtom,"|", (trialBox % NeighList(iList)%list(j, iAtom) ,j=1,trialBox % NeighList(iList)%nNeigh(iAtom))
!1        enddo
!        write(2,*)
        do iDisp = 1, size(disp)
!          write(2,"(A,A,1000(I3))") "NewList","|", (tempList(j, iDisp) ,j=1,tempNNei(iDisp))
          iAtom = disp(iDisp)%atmIndx
          trialBox % NeighList(iList) % nNeigh(iAtom) = tempNNei(iDisp)
          do iNei = 1, tempNNei(iDisp)
            neiIndx = tempList(iNei, iDisp)
            trialBox % NeighList(iList) % list(iNei, iAtom) =  neiIndx
            trialBox % NeighList(iList) % list( trialBox%NeighList(iList)%nNeigh(neiIndx)+1, neiIndx ) = iAtom
            trialBox%NeighList(iList)%nNeigh(neiIndx)= trialBox%NeighList(iList)%nNeigh(neiIndx) + 1
          enddo
        enddo
!        do iAtom = 1, trialBox%nMaxatoms
!          write(2,"(I3,A,1000(I3))") iAtom,"|", (trialBox % NeighList(iList)%list(j, iAtom) ,j=1,trialBox % NeighList(iList)%nNeigh(iAtom))
!        enddo
!        write(2,*)
      endif
!      write(2,*) "N", trialBox%NeighList(iList)%nNeigh(:)     
    enddo



  end subroutine
!===================================================================================
end module
!===================================================================================
