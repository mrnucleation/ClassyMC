!===================================================================================
! This module contains a Cell based Verlet List. 
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
      logical :: initialized = .false.
      integer :: maxAtoms
      integer :: nCells = 0
      integer :: integriyfails = 0
      integer, allocatable :: cellID(:)
      integer, allocatable :: nCellAtoms(:)
      integer, allocatable :: cellList(:, :)
      integer, allocatable :: atomList(:)


      integer :: nX, nY, nZ
      integer :: coeffX, coeffY, coeffZ
      real(dp) :: dX, dY, dZ
      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => CellRSqList_Constructor 
      procedure, pass :: BuildList => CellRSqList_BuildList 
      procedure, pass :: SortList => CellRSqList_SortList 
      procedure, pass :: GetNewList => CellRSqList_GetNewList
      procedure, pass :: AddMol => CellRSqList_AddMol
      procedure, pass :: PurgeAtom => CellRSqList_PurgeAtom
      procedure, pass :: SwapAtomLists => CellRSqList_SwapAtomLists
      procedure, pass :: SwapAtomType => CellRSqList_SwapAtomType
      procedure, pass :: GetNeighCount => CellRSqList_GetNeighCount
      procedure, pass :: GetCellIndex => CellRSqList_GetCellIndex
      procedure, pass :: GetBins => CellRSqList_GetBins
      procedure, pass :: GetCellAtoms => CellRSqList_GetCellAtoms
      procedure, pass :: IntegrityCheck => CellRSqList_IntegrityCheck
      procedure, pass :: ProcessIO => CellRSqList_ProcessIO
!      procedure, pass :: TransferList
      procedure, pass :: DeleteMol => CellRSqList_DeleteMol
      procedure, pass :: Prologue => CellRSqList_Prologue
      procedure, pass :: Epilogue => CellRSqList_Epilogue
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
    real(dp), parameter :: atomRadius = 0.5E0_dp  !Used to estimate an approximate volume of 
    integer :: AllocateStatus
    real(dp), pointer :: coords(:,:)



    self%parent => BoxArray(parentID)%box
    if(self%initialized) then
      return
    endif
    call self%parent%GetCoordinates(coords)
!     If no rCut value is given by the subroutine call attempt to pull
!     the rSq value from the parent box's energy function. The assumption being
!     that the neighborlist is used for the energy calculation routines.
    if( present(rCut) ) then
      self % rCut = rCut
      self % rCutSq = rCut * rCut
    else
      if(self%rCut > 0E0_dp) then
        self % rCutSq = (self%rCut)**2
      else
        self % rCut = self % parent % EFunc % Method % GetCutOff() + neighSkin
        self % rCutSq = (self%rCut)**2
      endif
    endif
    self % maxNei = ceiling(self%rCut**3/atomRadius**3) + 1
 
    write(nout,*) "Neighbor List CutOff:", self%rCut
    if(self%maxNei > self%parent%nMaxAtoms) then
      self%maxNei = self%parent%nMaxAtoms
    endif
    write(nout,*) "Neighbor List Maximum Neighbors:", self%maxNei

    if(.not. self%initialized) then
        allocate( self%cellID(1:self%parent%nMaxAtoms), stat=AllocateStatus )
        allocate( self%list(1:self%maxNei, 1:self%parent%nMaxAtoms), stat=AllocateStatus )
        allocate( self%nNeigh(1:self%parent%nMaxAtoms), stat=AllocateStatus )
        allocate( self%atomList(1:self%parent%nMaxAtoms), stat=AllocateStatus )

      if(.not. allocated(self%allowed) ) then
        allocate(self%allowed(1:nAtomTypes), stat=AllocateStatus )
        self%allowed = .true.
      endif
    endif

    self%list = 0
    self%nNeigh = 0 
    IF (AllocateStatus /= 0) then
        write(0,*) "Error Code:", AllocateStatus
        STOP "*** CellNeighRSQList: Memory Allocation Error! ***"
    endif


    self%initialized = .true.
    self%restrictType = .false.
  end subroutine
!===================================================================================
  subroutine CellRSqList_Prologue(self)
    implicit none
    class(CellRSqList), intent(inout) :: self


!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine CellRSqList_Epilogue(self)
    implicit none
    class(CellRSqList), intent(inout) :: self

!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine CellRSqList_Update(self)
    implicit none
    class(CellRSqList), intent(inout) :: self


!    call self%IntegrityCheck(1)
!    call self%DumpList(2)
  end subroutine
!===================================================================================
  function CellRSqList_GetCellIndex(self, binx, biny, binz) result(cellIndx)
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: binx, biny, binz
    integer :: cellIndx

    integer :: wrapX, wrapY, wrapZ

    !In the event the next cell over is across the boundary, wrap it back to the
    !other side of the box.
    wrapx = binx
    wrapy = biny
    wrapz = binz
    do while(wrapx < 0) 
      wrapx = wrapx + self%nx
    enddo
    do while(wrapy < 0) 
      wrapy = wrapy + self%ny
    enddo
    do while(wrapz < 0)
      wrapz = wrapz + self%nz
    enddo

    wrapx = mod(wrapx, self%nx) 
    wrapy = mod(wrapy, self%ny) 
    wrapz = mod(wrapz, self%nz) 



    cellIndx = 1
    cellIndx = cellIndx + self%coeffX * wrapx
    cellIndx = cellIndx + self%coeffY * wrapy
    cellIndx = cellIndx + self%coeffZ * wrapz

  end function
!===================================================================================
  subroutine CellRSqList_GetBins(self, cellIndx, binx, biny, binz) 
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: cellIndx
    integer, intent(out) :: binx, biny, binz

    integer :: remainder, curVal

    remainder = cellIndx - 1 

    curVal = int( real(remainder, dp)/real(self%coeffZ,dp) ) 
    binz = curVal 
    remainder = remainder - curVal * self%coeffZ

    curVal = int( real(remainder, dp)/real(self%coeffY,dp) ) 
    biny = curVal 
    remainder = remainder - curVal * self%coeffY

    curVal = int( real(remainder, dp)/real(self%coeffX,dp) ) 
    binx = curVal 
    remainder = remainder - curVal * self%coeffX

  end subroutine
!===================================================================================
  subroutine CellRSqList_BuildList(self, listindx)
    use SearchSort, only: QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: listindx

    integer :: iAtom, jNei, jAtom
!    integer :: maxAtoms
    integer :: nCells
    integer :: binx, biny, binz
    integer :: cellIndx
    integer :: nCellAtoms
    integer :: nx, ny, nz
    integer :: i,j,k, nNeigh

    real(dp) :: boxdim(1:2, 1:3)
    real(dp), pointer :: coords(:,:)
    real(dp) :: Lx, Ly, Lz
    real(dp) :: dx, dy, dz
    real(dp) :: rx, ry, rz, rsq


    !Compute the total number of cells and the size of each cell
    call self%parent%GetDimensions(boxdim)
    call self%parent%GetCoordinates(coords)
    self%maxAtoms = self%parent%nMaxAtoms
    Lx = boxdim(2,1) - boxdim(1,1)
    Ly = boxdim(2,2) - boxdim(1,2)
    Lz = boxdim(2,3) - boxdim(1,3)
    nx = floor(Lx/self%rCut)
    ny = floor(Ly/self%rCut)
    nz = floor(Lz/self%rCut)
    ! If nx < 1 then the box is small and there is only one cell
    if(nx < 1) nx = 1
    if(ny < 1) ny = 1
    if(nz < 1) nz = 1
    dx = Lx/real(nx, dp)
    dy = Ly/real(ny, dp)
    dz = Lz/real(nz, dp)
    nCells = nx * ny * nz
    self%nCells = nCells

    self%coeffX = 1
    self%coeffY = 1 + self%coeffX * (nx-1)
    self%coeffZ = 1 + self%coeffX * (nx-1)
    self%coeffZ = self%coeffZ + self%coeffY * (ny-1)

    self%dx = dx
    self%dy = dy
    self%dz = dz

    self%nx = nx
    self%ny = ny
    self%nz = nz


    !Check to see if the cell list is allocated and is large enough to accomodate
    !the number of cells in the system
    if(allocated(self%celllist)) then
      if(ubound(self%cellList, 2) < self%nCells) then
        deallocate(self%cellList)
        allocate(self%cellList(1:self%maxNei, 1:self%nCells))
        if(allocated(self%nCellAtoms)) then
          deallocate(self%nCellAtoms)
        endif
        allocate(self%nCellAtoms(1:nCells))

      endif

    else
!      allocate(self%cellID(1:maxAtoms))
      allocate(self%nCellAtoms(1:nCells))
      allocate(self%cellList(1:self%maxNei, 1:self%nCells))
    endif

    !Assign atoms to their respective cell
    self%cellID = 0
    self%nCellAtoms = 0
    self%cellList = 0
    do iAtom = 1, self%maxAtoms
      if( .not. self%parent%IsActive(iAtom) ) then
        cycle
      endif

      binx = floor((coords(1, iAtom) - boxdim(1,1))/dx)
      biny = floor((coords(2, iAtom) - boxdim(1,2))/dy)
      binz = floor((coords(3, iAtom) - boxdim(1,3))/dz)
      cellIndx = self % GetCellIndex(binx, biny, binz)
      self % cellID(iAtom) = cellIndx
      self % nCellAtoms(cellIndx) = self % nCellAtoms(cellIndx) + 1
      self % cellList(self%nCellAtoms(cellIndx), cellIndx) = iAtom
    enddo

    self%nNeigh = 0
    self%list = 0
    do iAtom = 1, self%maxAtoms
      if( self%parent%MolSubIndx(iAtom) > self%parent%NMol(self%parent%MolType(iAtom)) ) then
        cycle
      endif
      call self%GetCellAtoms(nCellAtoms, self%atomlist, atmIndx=iAtom)
      do jNei = 1, nCellAtoms
        jAtom = self%atomlist(jNei)
        if(jAtom <= iAtom) cycle
        if( self%parent%MolIndx(iAtom) == self%parent%MolIndx(jAtom)) cycle
        rx = coords(1, iAtom) - coords(1, jAtom)
        ry = coords(2, iAtom) - coords(2, jAtom)
        rz = coords(3, iAtom) - coords(3, jAtom)
        call self%parent%Boundary(rx,ry,rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          self%nNeigh(iAtom) = self%nNeigh(iAtom) + 1
          self%list( self%nNeigh(iAtom), iAtom ) = jAtom

          self%nNeigh(jAtom) = self%nNeigh(jAtom) + 1
          self%list( self%nNeigh(jAtom), jAtom ) = iAtom
        endif
      enddo
    enddo

    do iAtom = 1, self%parent%nMaxAtoms
      if( self%parent%IsActive(iAtom) ) then
        cycle
      endif
      nNeigh = self%nNeigh(iAtom)
      call QSort(self%list(1:nNeigh, iAtom))
    enddo

    self%sorted = .true.

  end subroutine
!===================================================================================
  subroutine CellRSqList_SortList(self)
    use SearchSort, only: QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer :: iAtom, nNeigh
    integer :: iCell, nCellAtoms

    if(.not. self%sorted) then
      do iAtom = 1, self%parent%nMaxAtoms
        if( self%parent%IsActive(iAtom) ) then
          cycle
        endif
        nNeigh = self%nNeigh(iAtom)
        call QSort(self%list(1:nNeigh, iAtom))
      enddo

      do iCell = 1, self%nCells
        nCellAtoms = self%nCellAtoms(iCell)
        call QSort(self%cellList(1:nCellAtoms, iCell))
      enddo

      self%sorted = .true.
    endif


  end subroutine
!===================================================================================
  subroutine CellRSqList_GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    use Common_MolInfo, only: nMolTypes
    use SearchSort, only: QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: iDisp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount
    integer :: jType, jAtom, j, iAtom, jNei
    integer :: jUp, jLow, molIndx, jMol
    integer :: binx, biny, binz, cellindx, nCellAtoms
    real(dp) :: xn, yn, zn
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: boxdim(1:2, 1:3)
    real(dp), pointer :: coords(:,:)

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

    !Get relevant box information from the parent.
    call self%parent%GetCoordinates(coords)
    call self%parent%GetDimensions(boxdim)

    binx = floor((xn - boxdim(1,1))/self%dx)
    biny = floor((yn - boxdim(1,2))/self%dy)
    binz = floor((zn - boxdim(1,3))/self%dz)
    cellIndx = self % GetCellIndex(binx, biny, binz)

    templist(:, iDisp) = 0
    tempNNei(iDisp) = 0

    call self%GetCellAtoms(nCellAtoms, self%atomlist, cellIndx=cellIndx)
    do jNei = 1, nCellAtoms
        jAtom = self%atomlist(jNei)
!        if(.not. self%parent%IsActive(jAtom) ) then
!          cycle
!        endif
        rx = xn - coords(1, jAtom)
        ry = yn - coords(2, jAtom)
        rz = zn - coords(3, jAtom)
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

      call QSort(templist(1:tempNNei(iDisp), iDisp))
!      molStart = molStart + self%parent%NMolMax(jType)
!    enddo
!    write(2,"(A, 1000(I3))") "New",   templist(1:tempNNei(iDisp), iDisp)
  end subroutine
!===================================================================================
  subroutine CellRSqList_GetCellAtoms(self, nAtoms, atomlist, atmindx, cellIndx) 
    use SearchSort, only: QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in), optional :: atmindx
    integer, intent(in) , optional:: cellindx
    integer, intent(out) :: nAtoms
    integer, intent(inout) :: atomlist(:)

    integer :: startCell, curCell
    integer :: jAtom
    integer :: i, j, k
    integer :: binx, biny, binz

    integer :: nVisits 
    integer :: visitedcells(1:30)
    

    if(present(atmindx)) then
      startCell = self%cellID(atmindx)
    elseif(present(cellindx)) then
      startCell = cellIndx
    endif
    call self%GetBins(startCell, binx, biny, binz)


    visitedCells = 0
    nVisits = 0

    !In 
    atomlist = 0
    nAtoms = 0
    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          curCell = self%GetCellIndex(binx+i, biny+j, binz+k)
          if( any( visitedCells(1:nVisits) == curCell) ) then
            cycle
          endif
          nVisits = nVisits + 1
          visitedCells(nVisits) = curCell
          do jAtom = 1, self%nCellAtoms(curCell)
            nAtoms = nAtoms + 1
            atomlist(nAtoms) = self%cellList(jAtom, curcell)
          enddo
        enddo
      enddo
    enddo

    call QSort(atomlist(1:nAtoms))




  end subroutine
!===================================================================================
  subroutine CellRSqList_SwapAtomLists(self, atmindx1, atmindx2)
    !--------------------------------
    ! This function switches
    !--------------------------------

    use SearchSort, only: BinarySearch, SimpleSearch, QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: atmindx1, atmindx2
    integer :: nNeigh1, nNeigh2
    integer :: templist(1:2*self%maxnei)
    integer :: jAtom, iNei, jNei
    integer :: cellID1, cellID2, nCellAtoms

    integer :: iSearch
    integer :: nSearch = 0
    integer :: neiIndx1, neiIndx2
    integer :: searchlist(1:2*self%maxnei)


    nNeigh1 = self%nNeigh(atmIndx1)
    nNeigh2 = self%nNeigh(atmIndx2)

    !If neither atom has a neighbor, nothing needs to be done. 
    if( (nNeigh1 == 0) .and. (nNeigh2 == 0) ) then
      return
    endif

    !The first task is to get a list of neighbors of both atoms in 
    !order to figure out which lists need to be modified. To prevent double
    !counting a unique list is 
    templist(1:nNeigh1+nNeigh2) = [self%list(1:nNeigh1, atmIndx1), &
                                     self%list(1:nNeigh2, atmIndx2)]
    call QSort(templist(1:nNeigh1+nNeigh2))
    nSearch = 1
    searchlist(1) = templist(1)
    do jNei = 2, nNeigh1+nNeigh2
      if(templist(jNei-1) /= templist(jNei)) then
        nSearch = nSearch + 1
        searchlist(nSearch) = templist(jNei)
      endif
    enddo


    !Now that we have the lists we need to update, we need to search the
    !lists to 

!    write(2,*) "---------------------------"
!    write(2,*) "Atom 1:", atmindx1
!    write(2,*) "Atom 2:", atmindx2
!    write(2,*) 
    call self%SortList
    do iSearch = 1, nSearch
      jAtom = searchlist(iSearch)
!      write(2,*) jAtom,"|",self%list(1:self%nNeigh(jAtom), jAtom)
!      call QSort(self%list(1:self%nNeigh(jAtom), jAtom))
!      write(2,*) jAtom,"|",self%list(1:self%nNeigh(jAtom), jAtom)
      neiIndx1 = BinarySearch(atmIndx1, self%list(1:self%nNeigh(jAtom), jAtom))
      neiIndx2 = BinarySearch(atmIndx2, self%list(1:self%nNeigh(jAtom), jAtom))
!      write(2,*) neiIndx1, neiIndx2
      if(neiIndx1 /= 0) then
        self%list(neiIndx1, jAtom) = atmIndx2
      endif
      if(neiIndx2 /= 0) then
        self%list(neiIndx2, jAtom) = atmIndx1
      endif     
      call QSort(self%list(1:self%nNeigh(jAtom), jAtom))
!      write(2,*) jAtom,"|",self%list(1:self%nNeigh(jAtom), jAtom)
!      write(2,*)
    enddo


    !Now that the other lists are reorganized, we need to swap the two lists
    if( nNeigh1 == 0 ) then
      self%list(1:nNeigh2, atmIndx1) = self%list(1:nNeigh2, atmIndx2)
    elseif( nNeigh2 == 0 ) then
      self%list(1:nNeigh1, atmIndx2) = self%list(1:nNeigh1, atmIndx1)
    else
      templist(1:nNeigh1) = self%list(1:nNeigh1, atmIndx1)
      self%list(1:nNeigh2, atmIndx1) = self%list(1:nNeigh2, atmIndx2)
      self%list(1:nNeigh1, atmIndx2) = templist(1:nNeigh1)
    endif


    !Swap the neighbor counters.
    self%nNeigh(atmIndx1) = nNeigh2
    self%nNeigh(atmIndx2) = nNeigh1


    !In addition swap cell IDs, to avoid problems with
    !improper swaping.  Both IDs must be updated at the same time.
    cellID1 = self%cellID(atmIndx1)
    cellID2 = self%cellID(atmIndx2)

    if( (cellID1 == 0) .and. (cellID2 == 0)) then
      write(0,*) "At least one atom must be initialized for this function to be called"
      error stop
    endif

    !Find both atom's locations in their respective cells.  If an atom
    !is not initialized then the cell ID will be 0
    neiIndx1 = 0
    neiIndx2 = 0
!    write(*,*) cellId1, cellId2
    if(cellID1 /= 0) then
      nCellAtoms = self%nCellAtoms(cellID1)
      neiIndx1 = BinarySearch(atmIndx1, self%cellList(1:nCellAtoms, cellID1))
    endif
    if(cellId2 /= 0) then
      nCellAtoms = self%nCellAtoms(cellID2)
      neiIndx2 = BinarySearch(atmIndx2, self%cellList(1:nCellAtoms, cellID2))
    endif
!    write(2,*) atmIndx1, atmIndx2
!    write(2,*) cellId1, cellId2
!    write(2,*) neiIndx1, neiIndx2
!    nCellAtoms = self%nCellAtoms(cellID1)
!    write(2,*) cellId1, "|", self%cellList(1:nCellAtoms, cellId1)
!    write(2,*) 

    !Swap the two atom IDs to their new location.
    if((cellID1 /= 0) .and. (neiIndx1 /= 0)) self%cellList(neiIndx1, cellId1) = atmIndx2
    if((cellID2 /= 0) .and. (neiIndx2 /= 0)) self%cellList(neiIndx2, cellId2) = atmIndx1
    self%cellID(atmIndx1) = cellID2
    self%cellID(atmIndx2) = cellID1

    if(cellID1 == cellID2) then
      nCellAtoms = self%nCellAtoms(cellID1)
      call QSort(self%cellList(1:nCellAtoms, cellId1))
    else
      if(cellId1 /= 0) then
        nCellAtoms = self%nCellAtoms(cellID1)
        call QSort(self%cellList(1:nCellAtoms, cellId1)) 
      endif
      if(cellId2 /= 0) then
        nCellAtoms = self%nCellAtoms(cellID2)
        call QSort(self%cellList(1:nCellAtoms, cellId2))
      endif
    endif
  end subroutine
!===================================================================================
  subroutine CellRSqList_PurgeAtom(self, atmindx1)
    !--------------------------------
    ! This function removes an atom's entry from
    !--------------------------------

    use SearchSort, only: BinarySearch, SimpleSearch, QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: atmindx1
    integer :: nNeigh
    integer :: templist(1:self%maxnei)
    integer :: jAtom, iNei, jNei
    integer :: cellID, nCellAtoms
    integer :: neiIndx



    !If the atom lacks a neighbor, nothing needs to be done. 
    if( (self%nNeigh(atmIndx1) == 0)  ) then
      return
    endif

    call self%SortList
!    write(2,*) "=================================================="
!    write(2,*) "Purge: ", atmIndx1
    do jNei = 1, self%nNeigh(atmIndx1)
      jAtom = self%list(jNei, atmIndx1)
      nNeigh = self%nNeigh(jAtom) 
!      write(2, *) "JAtom:", jAtom
!      write(2,"(999(1x, I3))") self%list(1:nNeigh, jAtom)
      neiIndx = BinarySearch(atmIndx1, self%list(1:nNeigh, jAtom))
      if(nNeigh > 1) then
        self%list(1:nNeigh-1, jAtom) = [self%list(1:neiIndx-1, jAtom), &
                                        self%list(neiIndx+1:nNeigh, jAtom)]
        self%nNeigh(jAtom) = self%nNeigh(jAtom) - 1
        nNeigh = nNeigh - 1
        call QSort(self%list(1:nNeigh, jAtom))
      else
        self%nNeigh(jAtom) = 0 
        nNeigh = 0
      endif

!      write(2,"(999(1x, I3))") self%list(1:nNeigh, jAtom)
!      write(2,*)
    enddo
    self%nNeigh(atmIndx1) = 0

    cellID = self%cellID(atmIndx1)
!    write(2,*) atmIndx1, cellID
!    write(2,*) "Cell List:", cellID
    nCellatoms = self%nCellAtoms(cellID)
!    write(2,"(999(1x, I3))") self%cellList(1:nCellatoms, cellID)
    neiIndx = BinarySearch(atmIndx1, self%cellList(1:nCellatoms, cellID))
    self%cellList(1:nCellAtoms-1, cellID) = [self%cellList(1:neiIndx-1, cellID), &
                                             self%cellList(neiIndx+1:nCellAtoms, cellID)]
    self%nCellAtoms(cellID) = nCellatoms - 1

    call QSort(self%cellList(1:nCellAtoms-1, cellID))

!    write(2,"(999(1x, I3))") self%cellList(1:nCellatoms, cellID)
    self%cellID(atmIndx1) = 0

!    write(2,*) "=================================================="


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
    real(dp), pointer :: coords(:,:)

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
      rx = coords(1, nAtom) - coords(1, jAtom)
      ry = coords(2, nAtom) - coords(2, jAtom)
      rz = coords(3, nAtom) - coords(3, jAtom)
      call self%parent%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < rCutSq) then
        nCount = nCount + 1
      endif
    enddo

  end function

!===================================================================================
  subroutine CellRSqList_AddMol(self, disp, tempList, tempNNei)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: QSort
    implicit none
    class(CellRSqList), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iDisp, iAtom, iNei, nNei, neiIndx, j
    integer :: binX, binY, binZ, cellIndx
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: boxdim(1:2, 1:3)


    write(*,*) "Add"
    select type(disp)
      class is(Addition)
        call self%parent%GetDimensions(boxdim)
        do iDisp = 1, size(disp)
!          write(2,"(A,A,1000(I3,1x))") "NewList","|", (tempList(j, iDisp) ,j=1,tempNNei(iDisp))
           iAtom = disp(iDisp)%atmIndx
           binx = floor((disp(iDisp)%x_new - boxdim(1,1))/self%dx)
           biny = floor((disp(iDisp)%y_new - boxdim(1,2))/self%dy)
           binz = floor((disp(iDisp)%z_new - boxdim(1,3))/self%dz)
           cellIndx = self % GetCellIndex(binx, biny, binz)
           self % cellID(iAtom) = cellIndx
           self % nCellAtoms(cellIndx) = self % nCellAtoms(cellIndx) + 1
           self % cellList(self%nCellAtoms(cellIndx), cellIndx) = iAtom

           self%nNeigh(iAtom) = tempNNei(iDisp)
           self%list(1:tempNNei(iDisp), iAtom ) = templist(1:tempNNei(iDisp), iDisp)
!           do iNei = 1, tempNNei(iDisp)
!               neiIndx = tempList(iNei, iDisp)
!               self % list(iNei, iAtom) =  neiIndx
!               self % list( self%nNeigh(neiIndx)+1, neiIndx ) = iAtom
!               self%nNeigh(neiIndx)= self%nNeigh(neiIndx) + 1
!           enddo
       enddo

   end select



  end subroutine
!===================================================================================
  subroutine CellRSqList_DeleteMol(self, molIndx, topIndx)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: BinarySearch, SimpleSearch
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx
    integer :: nStart, nEnd
    integer :: iAtom, iNei, jNei, nType, j
    integer :: topStart
    integer :: topEnd
    integer :: cellIndx
    integer :: atmIndx, topAtom
    integer :: curNei, curIndx, nNei

    nStart = self % parent % MolStartIndx(molIndx)
    topStart = self % parent % MolStartIndx(topIndx)

    nEnd = self % parent % MolEndIndx(molIndx)
    topEnd = self % parent % MolEndIndx(topIndx)
    nType = self % parent % MolType(nStart)

    call self%sortlist
    do iAtom = 1, MolData(nType)%nAtoms
      atmIndx = nEnd - iAtom + 1
      topAtom = topEnd - iAtom + 1
!      write(*,*) atmIndx, topAtom
      if(topAtom /= atmIndx) then
        call self%SwapAtomLists(atmIndx, topAtom)
      endif
!      call self%sortlist
      call self%PurgeAtom(topAtom)
    enddo
    call self%sortlist
!    write(*,*) "Deleted Mol"
!    call self%IntegrityCheck(1)
    return

!    write(2,*) "----------------------------"
!    write(2,*) "Delete"
!    write(2,*) "Removed Mol:", molIndx
!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3,1x))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "Sorted?:", self%sorted
    call self%sortlist

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
!      write(2,"(I3,A,1000(I3,1x))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
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
!        write(2,"(I3,A,1000(I3,1x))") iNei,"|", (self%list(j, iNei) ,j=1,self%nNeigh(iNei))
!      enddo


      !Delete atom's the entry from the cell list.
      cellIndx = self%cellID(atmIndx)
      nNei = self%nCellAtoms(cellIndx)
      curIndx = SimpleSearch( atmIndx, self%celllist(1:nNei, cellIndx) )
      if(curIndx /= 0) then
         self%celllist(1:nNei, cellIndx ) = [self%celllist(1:curIndx-1, cellIndx), &
                                             self%celllist(curIndx+1:nNei, cellIndx) ]
         self%nCellAtoms(cellIndx) = self%nCellAtoms(cellIndx) - 1
      endif

      !Re-index the top to the old from the cell list.
      if(atmIndx /= topAtom) then
          cellIndx = self%cellID(topAtom)
          nNei = self%nCellAtoms(cellIndx)
          if(nNei > 0) then
            curIndx = SimpleSearch( topAtom, self%celllist(1:nNei, cellIndx) )
            if(curIndx /= 0) then
              self%cellList(curIndx, cellIndx) = atmIndx
            endif
            self%cellID(atmIndx) = self%cellID(topAtom)
          endif

       endif
       self%cellID(topAtom) = 0



    enddo

!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3,1x))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "---------------------------------"


    self % sorted = .false.
    call self%SortList


  end subroutine
!===================================================================================
  subroutine CellRSqList_SwapAtomType(self, disp, topIndx)
    !----------------
    ! Routine used to change the type of an atom.  Primarily used for atom swap moves.
    ! Not currently designed to work with molecules. 
    ! Variables
    !    input
    !        disp => Displacement class variable which contains information related to
    !                what changed in the system.  For example this might contain the old/new
    !                positions of an atom that was shifted. This routine expects an
    !                AtomExchange class which contains the array location of 
    !                the atom being changed and it's new array location.
    !    
    !    function variables
    !       iAtomNew => Index of the atom's new position in the array
    !       iAtomOld => Index of the atom's old position in the array
    !       iNei => Loop integer for looping over the number of neighbors in the neighborlist list
    !       cellIndx => Cell ID of the current atom. 
    !       curIndx => Atom Index returned by the search algorithm
    !
    !---------
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: SimpleSearch
    implicit none
    class(CellRSqList), intent(inout) :: self
    class(AtomExchange), intent(in) :: disp(:)
    integer, intent(in) :: topIndx
    integer :: iAtomNew, iAtomOld, iNei, jAtom, iAtom, nNei
    integer :: topAtom
    integer :: cellIndx, curIndx
    integer :: molIndxNew, molIndxOld

    iAtomNew = disp(1)%newAtmIndx
    iAtomOld = disp(1)%oldAtmIndx
    topAtom = self%parent%MolStartIndx(topIndx)
!    write(2,*) "New: ", iAtomNew
!    write(2,*) "Old: ", iAtomOld
!    write(2,*) "Top: ", topAtom
    call self%SwapAtomLists(iAtomNew, iAtomOld)
    if(iAtomOld /= topAtom) then
      call self%SwapAtomLists(iAtomOld, topAtom)
    endif

    self%nNeigh(topAtom) = 0
    self%cellId(topAtom) = 0
    call self%SortList
!    call self%IntegrityCheck(1)

    return


!    write(2,*) "--------------------------------------"
!    write(2,*) "AtomExchange"
    iAtomNew = disp(1)%newAtmIndx
    iAtomOld = disp(1)%oldAtmIndx
    topAtom = self%parent%MolStartIndx(topIndx)
!    write(2,*) "OutIndex", iAtomOld
!    write(2,*) "InIndex", iAtomNew
!    write(2,*) "TopIndex", topAtom
!    write(2,*)
!    do iAtom = 1, self%parent%nMaxAtoms
!      nNei = self%nNeigh(iatom)
!      write(2,*) iAtom, ":", self % list(1:nNei, iAtom)
!    enddo
!    write(2,*)



    molIndxNew = self%parent%MolIndx(iAtomNew)
    molIndxOld = self%parent%MolIndx(iAtomOld)

    !Search through the neighborlists and replace the old atom index for
    !a new atom index
    nNei = self%nNeigh(iAtomOld)
    do iNei = 1, nNei
      jAtom = self % list(iNei, iAtomOld)
      curIndx = SimpleSearch( iAtomOld, self%list(1:nNei, jAtom) )
      if(curIndx /= 0) then
        self%list(curIndx, jAtom) = iAtomNew
      endif
    enddo

!    write(2,*)
!    do iAtom = 1, self%parent%nMaxAtoms
!      if(.not. self%parent%IsActive(iAtom)) cycle
!      nNei = self%nNeigh(iatom)
!      write(2,*) iAtom, ":", self % list(1:nNei, iAtom)
!    enddo
!    write(2,*) 

    nNei = self%nNeigh(topAtom)
    do iNei = 1, nNei
      jAtom = self % list(iNei, topAtom)
      if(jAtom == iAtomNew) then
        cycle
      endif
      curIndx = SimpleSearch( topAtom, self%list(1:nNei, jAtom) )
      if(curIndx /= 0) then
        self%list(curIndx, jAtom) = iAtomOld
      endif
    enddo
    if(nNei > 0) then
        curIndx = SimpleSearch( topAtom, self%list(1:nNei, iAtomOld) )
        if(curIndx /= 0) then
          self%list(curIndx, iAtomOld) = iAtomOld
        endif
    endif

!    write(2,*)
!    do iAtom = 1, self%parent%nMaxAtoms
!      nNei = self%nNeigh(iatom)
!      write(2,*) iAtom, ":", self % list(1:nNei, iAtom)
!    enddo
!    write(2,*) ".."
    !Copy the neighbor of the atom being swapped list from the old location
    !to the new location
    do iNei = 1, self%nNeigh(iAtomOld)
      self % list(iNei, iAtomNew) =  self % list(iNei, iAtomOld)
    enddo
    do iNei = 1, self%nNeigh(topAtom)
      self % list(iNei, iAtomOld) =  self % list(iNei, topAtom)
    enddo


    self%nNeigh(iAtomNew) = self%nNeigh(iAtomOld) 
    self%nNeigh(iAtomOld) = self%nNeigh(topAtom) 
    self%nNeigh(topAtom) = 0
!    write(2,*)
!    do iAtom = 1, self%parent%nMaxAtoms
!      nNei = self%nNeigh(iatom)
!      write(2,*) iAtom, ":", self % list(1:nNei, iAtom)
!    enddo
!    write(2,*)

!    write(2,*) "Cell:"
    cellIndx = self%cellID(iAtomOld)
    nNei = self%nCellAtoms(cellIndx)
!    write(2,*) self%cellList(1:nNei, cellIndx)
    curIndx = SimpleSearch( iAtomOld, self%cellList(1:nNei, cellIndx) )
    self%cellList(curIndx, cellIndx) = iAtomNew
    self%cellID(iAtomNew) = cellIndx
!    write(2,*) self%cellList(1:nNei, cellIndx)

    if(topAtom == iAtomOld) then
        self%cellID(topAtom) = 0
    else
        cellIndx = self%cellID(topAtom)
        nNei = self%nCellAtoms(cellIndx)
!        write(2,*) self%cellList(1:nNei, cellIndx)
        curIndx = SimpleSearch(topAtom, self%cellList(1:nNei, cellIndx) )
        self%cellList(curIndx, cellIndx) = iAtomOld
        self%cellID(iAtomOld) = cellIndx
        self%cellID(topAtom) = 0
!        write(2,*) self%cellList(1:nNei, cellIndx)
    endif

    self % sorted = .false.


!    call self%IntegrityCheck(1)

!    acall self%DeleteMol(molIndxOld, topIndx)
    




  end subroutine
!====================================================================
  subroutine CellRSqList_IntegrityCheck(self, listindx)
    !-------------------------------
    ! Performs some sanity checks on the neighborlist to ensure
    ! that it was properly updated.
    ! 
    ! Input:
    !    integer listindx => This list's array index within its parent
    !                        box class. 
    ! Local Variables:
    !    integer iAtom => Loop variable over the atoms in the system
    !    integer jNei => Loop variable for looping over the neighbors of atom iAtom
    !    integer jAtom => Atomic Index for the jNei-th neighbor of atom iAtom
    !    integer nNeigh => Number of Neighbors of atom iAtom
    !-------------------------------
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(CellRSqList), intent(inout) :: self
    integer, intent(in) :: listindx
    integer :: iAtom, jAtom, jAtom2, jNei
    integer :: nNeigh, nNeigh2

    integer :: templist(1:self%maxnei, 1:self%maxAtoms)
    integer :: tempNNeigh(1:self%maxAtoms)




    !This first loop scans the neighborlist to see if somehow
    !an atom's index ended up on it's own neighborlist.
    !This usually indicates a transfer error.
    atomloop: do iAtom = 1, self%maxAtoms
      if(.not. self%parent%IsActive(iAtom) ) then
        cycle
      endif
      nNeigh = self%nNeigh(iAtom)
      do jNei = 1, nNeigh
        jAtom = self%list(jNei, iAtom)
        if(iAtom == jAtom) then
          write(0,*) "ERROR: NeighList Integrity Failure!"
          write(0,*) "An atom's neighbor list contains itself!"
          write(0,*) "This is usually a sign of a problem during a swap move"
          error stop
          return
        endif
      enddo
    enddo atomloop


    templist = self%list
    tempNNeigh = self%nNeigh

    call self%BuildList(1)
    atomloop2: do iAtom = 1, self%maxAtoms
      if(.not. self%parent%IsActive(iAtom) ) then
        cycle
      endif
      nNeigh = self%nNeigh(iAtom)
      nNeigh2 = tempNNeigh(iAtom)
      if(nNeigh /= nNeigh2) then
        write(0,*) "ERROR: NeighList Integrity Failure!"
        write(0,*) "An atom's neighborlist changed significantly after"
        write(0,*) "a rebuild!"
        write(0,*) "N-Neigh:", nNeigh, nNeigh2
        error stop
      endif

      do jNei = 1, nNeigh
        jAtom = self%list(jNei, iAtom)
        jAtom2 = templist(jNei, iAtom)
        if(jAtom /= jAtom2) then
          write(0,*) "ERROR: NeighList Integrity Failure!"
          write(0,*) "An atom's neighborlist changed significantly after"
          write(0,*) "a rebuild!"
          write(0,*) jNei, "|",jAtom, jAtom2
          error stop
        endif
      enddo
    enddo atomloop2

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
end module
!===================================================================================
