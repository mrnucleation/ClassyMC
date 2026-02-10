!===================================================================================
! This module contains a pure Cell List for spatial decomposition.
! Unlike CellRSqList, this does not maintain a distance-based neighbor list.
! Instead, it only tracks which atoms are in which cells.
! This is useful for:
!   - Force fields that want to iterate over cells directly
!   - Spatial decomposition for parallelization
!   - Quick proximity queries without maintaining full neighbor lists
!===================================================================================
module CellListDef
use VarPrecision
use CoordinateTypes
use Template_SimBox, only: SimBox
use SimpleSimBox, only: SimpleBox
use Template_NeighList, only: NeighListDef

  type, public, extends(NeighListDef) :: CellList
      logical :: initialized = .false.
      integer :: maxAtoms = 0
      
      ! Cell structure
      integer :: nCells = 0
      integer, allocatable :: cellID(:)          ! Cell ID for each atom
      integer, allocatable :: nCellAtoms(:)      ! Number of atoms in each cell
      integer, allocatable :: cellList(:, :)     ! List of atoms in each cell
      integer, allocatable :: atomList(:)        ! Scratch array for atom queries

      ! Cell dimensions
      integer :: nX = 0, nY = 0, nZ = 0
      integer :: coeffX, coeffY, coeffZ
      real(dp) :: dX = -1E0_dp
      real(dp) :: dY = -1E0_dp
      real(dp) :: dZ = -1E0_dp
      real(dp) :: cellSize = 3.0E0_dp           ! Target cell size
      
      ! Box origin (lower corner)
      real(dp) :: originX = 0E0_dp
      real(dp) :: originY = 0E0_dp
      real(dp) :: originZ = 0E0_dp
      
      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => CellList_Constructor 
      procedure, pass :: BuildList => CellList_BuildList 
      procedure, pass :: BuildNeighborPairs => CellList_BuildNeighborPairs
      procedure, pass :: GetNewList => CellList_GetNewList
      procedure, pass :: GetCellIndex => CellList_GetCellIndex
      procedure, pass :: GetBins => CellList_GetBins
      procedure, pass :: GetCellAtoms => CellList_GetCellAtoms
      procedure, pass :: GetNeighborCells => CellList_GetNeighborCells
      procedure, pass :: UpdateCellDimensions => CellList_UpdateCellDimensions
      procedure, pass :: InsertAtom => CellList_InsertAtom
      procedure, pass :: RemoveAtom => CellList_RemoveAtom
      procedure, pass :: MoveAtom => CellList_MoveAtom
      procedure, pass :: AddMol => CellList_AddMol
      procedure, pass :: DeleteMol => CellList_DeleteMol
      procedure, pass :: ProcessIO => CellList_ProcessIO
      procedure, pass :: Prologue => CellList_Prologue
      procedure, pass :: Epilogue => CellList_Epilogue
      procedure, pass :: Update => CellList_Update
      procedure, pass :: PrintList => CellList_PrintList
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine CellList_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use SimpleSimBox, only: SimpleBox
    use Common_MolInfo, only: nAtomTypes
    use Common_NeighData, only: neighSkin
    use ParallelVar, only: nout
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 0.85E0_dp  ! Used to estimate neighbor volume
    integer :: AllocateStatus
    integer :: maxCellAtoms

    select type(parbox => BoxArray(parentID)%box)
      type is(SimpleBox)
        error stop "ERROR! To use a cell list you must have a box with a well defined boundary!"
    end select
    
    self%parent => BoxArray(parentID)%box

    if(self%initialized) then
      return
    endif

    ! Set cell size: use provided rCut, or existing value, or get from forcefield
    if(present(rCut)) then
      self%cellSize = rCut
      self%rCut = rCut
      self%rCutSq = rCut * rCut
    else if(self%rCut > 0E0_dp) then
      self%cellSize = self%rCut
      self%rCutSq = self%rCut * self%rCut
    else
      ! Get cutoff from forcefield + neighbor skin
      self%rCut = self%parent%EFunc%Method%GetCutOff() + neighSkin
      self%cellSize = self%rCut
      self%rCutSq = self%rCut * self%rCut
    endif

    self%maxAtoms = self%parent%nMaxAtoms
    
    ! Estimate max neighbors based on cutoff volume (same formula as RSqList/CellRSqList)
    self%maxNei = ceiling(self%rCut**3 / atomRadius**3) + 1
    if(self%maxNei > self%maxAtoms) then
      self%maxNei = self%maxAtoms + 1
    endif
    write(nout,*) "CellList Maximum Neighbors:", self%maxNei
    
    ! Estimate max atoms per cell
    maxCellAtoms = max(200, self%maxAtoms / 10)
    
    allocate( self%cellID(1:self%maxAtoms), stat=AllocateStatus )
    allocate( self%atomList(1:self%maxAtoms), stat=AllocateStatus )
    
    ! Allocate neighbor list arrays
    allocate( self%list(1:self%maxNei+1, 1:self%maxAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:self%maxAtoms), stat=AllocateStatus )
    allocate( self%templist(1:self%maxNei+1, 1:1), stat=AllocateStatus )
    allocate( self%tempNNeigh(1:1), stat=AllocateStatus )

    if(.not. allocated(self%allowed) ) then
      allocate(self%allowed(1:nAtomTypes), stat=AllocateStatus )
      self%allowed = .true.
    endif

    self%cellID = 0
    self%list = 0
    self%nNeigh = 0

    IF (AllocateStatus /= 0) then
        write(0,*) "Error Code:", AllocateStatus
        error STOP "*** CellList: Memory Allocation Error! ***"
    endif

    write(nout,*) "Cell List initialized with cell size:", self%cellSize
    write(nout,*) "WARNING: CellList does not maintain neighbor pairs."
    write(nout,*) "         Use CellRSqList for standard energy calculations."

    self%initialized = .true.
    self%restrictType = .false.
  end subroutine
!===================================================================================
  subroutine CellList_Prologue(self)
    implicit none
    class(CellList), intent(inout) :: self
  end subroutine
!===================================================================================
  subroutine CellList_Epilogue(self)
    implicit none
    class(CellList), intent(inout) :: self
  end subroutine
!===================================================================================
  subroutine CellList_Update(self)
    implicit none
    class(CellList), intent(inout) :: self
  end subroutine
!===================================================================================
  subroutine CellList_UpdateCellDimensions(self)
    ! Update cell dimensions based on current box size without rebuilding atom assignments
    ! This is called when we need current cell geometry for queries (like GetNewList)
    ! but don't want the overhead of a full rebuild
    implicit none
    class(CellList), intent(inout) :: self
    real(dp) :: boxdim(1:2, 1:3)
    real(dp) :: Lx, Ly, Lz
    
    ! Get current box dimensions
    call self%parent%GetDimensions(boxdim)
    
    ! Store origin (lower corner of box)
    self%originX = boxdim(1,1)
    self%originY = boxdim(1,2)
    self%originZ = boxdim(1,3)
    
    Lx = boxdim(2,1) - boxdim(1,1)
    Ly = boxdim(2,2) - boxdim(1,2)
    Lz = boxdim(2,3) - boxdim(1,3)
    
    ! Calculate number of cells in each direction
    self%nx = max(1, floor(Lx/self%cellSize))
    self%ny = max(1, floor(Ly/self%cellSize))
    self%nz = max(1, floor(Lz/self%cellSize))
    
    ! Calculate actual cell sizes
    self%dx = Lx/real(self%nx, dp)
    self%dy = Ly/real(self%ny, dp)
    self%dz = Lz/real(self%nz, dp)
    
    ! Calculate cell index coefficients
    self%coeffX = 1
    self%coeffY = 1 + self%coeffX * (self%nx - 1)
    self%coeffZ = 1 + self%coeffX * (self%nx - 1)
    self%coeffZ = self%coeffZ + self%coeffY * (self%ny - 1)
    
  end subroutine
!===================================================================================
  function CellList_GetCellIndex(self, binx, biny, binz) result(cellIndx)
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: binx, biny, binz
    integer :: cellIndx
    integer :: wrapX, wrapY, wrapZ

    ! Wrap bin indices to handle periodic boundaries
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
  subroutine CellList_GetBins(self, cellIndx, binx, biny, binz) 
    implicit none
    class(CellList), intent(inout) :: self
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

  end subroutine
!===================================================================================
  subroutine CellList_BuildList(self, listindx)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: listindx

    integer :: iAtom, iType
    integer :: binx, biny, binz
    integer :: typeStart, typeEnd
    integer :: cellIndx, nCells
    integer :: maxCellAtoms
    real(dp) :: boxdim(1:2, 1:3)
    real(dp) :: Lx, Ly, Lz
    real(dp), pointer :: coords(:,:)

    ! Get box dimensions
    call self%parent%GetDimensions(boxdim)
    call self%parent%GetCoordinates(coords)
    
    ! Store origin (lower corner of box)
    self%originX = boxdim(1,1)
    self%originY = boxdim(1,2)
    self%originZ = boxdim(1,3)
    
    Lx = boxdim(2,1) - boxdim(1,1)
    Ly = boxdim(2,2) - boxdim(1,2)
    Lz = boxdim(2,3) - boxdim(1,3)
    
    ! Calculate number of cells in each direction
    self%nx = max(1, floor(Lx/self%cellSize))
    self%ny = max(1, floor(Ly/self%cellSize))
    self%nz = max(1, floor(Lz/self%cellSize))
    
    ! Calculate actual cell sizes
    self%dx = Lx/real(self%nx, dp)
    self%dy = Ly/real(self%ny, dp)
    self%dz = Lz/real(self%nz, dp)
    
    nCells = self%nx * self%ny * self%nz
    self%nCells = nCells

    ! Calculate cell index coefficients
    self%coeffX = 1
    self%coeffY = 1 + self%coeffX * (self%nx - 1)
    self%coeffZ = 1 + self%coeffX * (self%nx - 1)
    self%coeffZ = self%coeffZ + self%coeffY * (self%ny - 1)

    ! Estimate max atoms per cell
    maxCellAtoms = max(50, (self%maxAtoms / nCells) * 4 + 10)

    ! Allocate or reallocate cell arrays
    if(allocated(self%cellList)) then
      if(ubound(self%cellList, 2) < nCells) then
        deallocate(self%cellList)
        deallocate(self%nCellAtoms)
        allocate(self%cellList(1:maxCellAtoms, 1:nCells))
        allocate(self%nCellAtoms(1:nCells))
      endif
    else
      allocate(self%nCellAtoms(1:nCells))
      allocate(self%cellList(1:maxCellAtoms, 1:nCells))
    endif

    ! Initialize
    self%cellID = 0
    self%nCellAtoms = 0
    self%cellList = 0

    ! Assign atoms to cells
    do iType = 1, nMolTypes
      call self%parent%GetTypeAtoms(iType, typeStart, typeEnd)
      if( (typeStart < 1) .or. (typeStart > typeEnd) ) cycle
      do iAtom = typeStart, typeEnd
        binx = floor((coords(1, iAtom) - boxdim(1,1))/self%dx)
        biny = floor((coords(2, iAtom) - boxdim(1,2))/self%dy)
        binz = floor((coords(3, iAtom) - boxdim(1,3))/self%dz)
        cellIndx = self%GetCellIndex(binx, biny, binz)
        self%cellID(iAtom) = cellIndx
        call self%InsertAtom(cellIndx, iAtom)
      enddo
    enddo

    ! Now build the neighbor list arrays (self%list and self%nNeigh)
    call self%BuildNeighborPairs(coords)

    self%sorted = .true.

  end subroutine
!===================================================================================
  subroutine CellList_BuildNeighborPairs(self, coords)
    ! Build the neighbor pair list (self%list, self%nNeigh) using cell structure
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(CellList), intent(inout) :: self
    real(dp), pointer, intent(in) :: coords(:,:)
    
    integer :: iAtom, jAtom, iType, jNei, cellNei
    integer :: typeStart, typeEnd
    integer :: cellIndx, nNeighCells, nCellAtoms
    integer :: molIndx, jMol
    integer :: neighborCells(1:30)
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: xi, yi, zi

    ! Initialize neighbor counts
    self%nNeigh = 0

    ! Loop over all atoms
    do iType = 1, nMolTypes
      call self%parent%GetTypeAtoms(iType, typeStart, typeEnd)
      if( (typeStart < 1) .or. (typeStart > typeEnd) ) cycle
      
      do iAtom = typeStart, typeEnd
        call self%parent%GetAtomData(iAtom, molIndx)
        xi = coords(1, iAtom)
        yi = coords(2, iAtom)
        zi = coords(3, iAtom)
        cellIndx = self%cellID(iAtom)
        
        ! Get neighboring cells
        call self%GetNeighborCells(cellIndx, nNeighCells, neighborCells)
        
        ! Loop over atoms in neighboring cells
        do cellNei = 1, nNeighCells
          call self%GetCellAtoms(neighborCells(cellNei), nCellAtoms, self%atomList)
          
          do jNei = 1, nCellAtoms
            jAtom = self%atomList(jNei)
            
            ! Only count each pair once (jAtom > iAtom)
            if(jAtom <= iAtom) cycle
            
            ! Skip atoms in the same molecule
            call self%parent%GetAtomData(jAtom, jMol)
            if(molIndx == jMol) cycle
            
            ! Calculate distance
            rx = xi - coords(1, jAtom)
            ry = yi - coords(2, jAtom)
            rz = zi - coords(3, jAtom)
            call self%parent%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            
            ! Add to neighbor list if within cutoff
            if(rsq < self%rCutSq) then
              ! Add jAtom to iAtom's list
              self%nNeigh(iAtom) = self%nNeigh(iAtom) + 1
              if(self%nNeigh(iAtom) <= self%maxNei) then
                self%list(self%nNeigh(iAtom), iAtom) = jAtom
              endif
              
              ! Add iAtom to jAtom's list (symmetric)
              self%nNeigh(jAtom) = self%nNeigh(jAtom) + 1
              if(self%nNeigh(jAtom) <= self%maxNei) then
                self%list(self%nNeigh(jAtom), jAtom) = iAtom
              endif
            endif
          enddo
        enddo
      enddo
    enddo

  end subroutine
!===================================================================================
  subroutine CellList_GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    use Common_MolInfo, only: nMolTypes
    use SearchSort, only: QSort, IsSorted
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: iDisp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount
    
    integer :: jAtom, iAtom, jNei, cellNei
    integer :: molIndx, jMol
    integer :: binx, biny, binz, cellIndx
    integer :: nNeighCells, nCellAtoms
    integer :: neighborCells(1:30)
    integer :: ub, lb
    real(dp) :: xn, yn, zn
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: boxdim(1:2, 1:3)
    real(dp) :: Lx, Ly, Lz
    real(dp) :: expected_Lx, expected_Ly, expected_Lz, dimension_change
    real(dp), pointer :: coords(:,:)

    if(present(nCount)) then
      nCount = 0
    endif
    
    select type(disp)
      class is (Addition)
        disp%listIndex = iDisp
        iAtom = disp%atmIndx
        call self%parent%GetAtomData(iAtom, molIndx)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
      class is (Displacement)
        disp%newlist = .true.
        disp%listIndex = iDisp
        iAtom = disp%atmIndx
        call self%parent%GetAtomData(iAtom, molIndx)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
      class default
        error stop "CellList_GetNewList: Unsupported displacement type"
    end select

    ! Check if cell dimensions need updating (e.g., after volume moves)
    ! Compare current box dimensions to stored cell structure dimensions
    call self%parent%GetDimensions(boxdim)
    Lx = boxdim(2,1) - boxdim(1,1)
    Ly = boxdim(2,2) - boxdim(1,2)
    Lz = boxdim(2,3) - boxdim(1,3)
    
    ! If box size has changed by more than a small threshold, update cell dimensions
    ! This is needed for correct cell index calculations after volume moves
    if (self%initialized) then
      expected_Lx = self%dx * self%nx
      expected_Ly = self%dy * self%ny
      expected_Lz = self%dz * self%nz
      dimension_change = max(abs(Lx - expected_Lx)/expected_Lx, &
                             abs(Ly - expected_Ly)/expected_Ly, &
                             abs(Lz - expected_Lz)/expected_Lz)
      if (dimension_change > 0.01E0_dp) then  ! 1% change threshold
        call self%UpdateCellDimensions()
      endif
    endif

    ! Get coordinates
    call self%parent%GetCoordinates(coords)

    ! Calculate which cell the new position falls into
    binx = floor((xn - self%originX)/self%dx)
    biny = floor((yn - self%originY)/self%dy)
    binz = floor((zn - self%originZ)/self%dz)
    cellIndx = self%GetCellIndex(binx, biny, binz)

    ! Validate bounds
    lb = lbound(tempList, 2)
    ub = ubound(tempList, 2)
    if((lb /= 1) .or. (iDisp > ub)) then
      write(0,*) "Largest Displacement Index:", iDisp
      write(0,*) "Bounds:", lb, ub
      error stop "Bounds Error: Array Pointer has invalid bounds"
    endif
    
    ! Initialize output
    tempList(:, iDisp) = 0
    tempNNei(iDisp) = 0

    ! Get all neighboring cells (including the cell itself)
    call self%GetNeighborCells(cellIndx, nNeighCells, neighborCells)

    ! Loop over all atoms in neighboring cells
    do cellNei = 1, nNeighCells
      call self%GetCellAtoms(neighborCells(cellNei), nCellAtoms, self%atomList)
      
      do jNei = 1, nCellAtoms
        jAtom = self%atomList(jNei)
        call self%parent%GetAtomData(jAtom, jMol)
        
        ! Skip atoms in the same molecule
        if(molIndx == jMol) cycle
        
        ! Calculate distance
        rx = xn - coords(1, jAtom)
        ry = yn - coords(2, jAtom)
        rz = zn - coords(3, jAtom)
        call self%parent%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        
        ! Add to neighbor list if within cutoff
        if(rsq < self%rCutSq) then
          tempNNei(iDisp) = tempNNei(iDisp) + 1
          if(tempNNei(iDisp) <= size(tempList, 1)) then
            tempList(tempNNei(iDisp), iDisp) = jAtom
          else
            error stop "CellList_GetNewList: tempList overflow"
          endif
        endif
        
        if(present(rCount)) then
          if(rsq < rCount*rCount) then
            nCount = nCount + 1
          endif
        endif
      enddo
    enddo

    ! Sort the neighbor list
    if(tempNNei(iDisp) > 1) then
      if(.not. IsSorted(tempList(1:tempNNei(iDisp), iDisp))) then
        call QSort(tempList(1:tempNNei(iDisp), iDisp))
      endif
    endif

  end subroutine
!===================================================================================
  subroutine CellList_GetCellAtoms(self, cellIndx, nAtoms, atomlist) 
    ! Returns the atoms in a single cell
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: cellIndx
    integer, intent(out) :: nAtoms
    integer, intent(inout) :: atomlist(:)
    integer :: i

    nAtoms = self%nCellAtoms(cellIndx)
    do i = 1, nAtoms
      atomlist(i) = self%cellList(i, cellIndx)
    enddo

  end subroutine
!===================================================================================
  subroutine CellList_GetNeighborCells(self, cellIndx, nCells, cellIds)
    ! Returns all neighboring cells (including the cell itself)
    ! Returns up to 27 cells in 3D (3x3x3 cube of cells)
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: cellIndx
    integer, intent(out) :: nCells
    integer, intent(out) :: cellIds(:)
    integer :: binx, biny, binz
    integer :: i, j, k, curCell
    integer :: nVisits
    integer :: visitedCells(1:30)

    call self%GetBins(cellIndx, binx, biny, binz)

    visitedCells = 0
    nVisits = 0
    nCells = 0

    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          curCell = self%GetCellIndex(binx+i, biny+j, binz+k)
          ! Avoid duplicates (can happen with small boxes)
          if( any( visitedCells(1:nVisits) == curCell) ) cycle
          nVisits = nVisits + 1
          visitedCells(nVisits) = curCell
          nCells = nCells + 1
          cellIds(nCells) = curCell
        enddo
      enddo
    enddo

  end subroutine
!===================================================================================
  subroutine CellList_InsertAtom(self, cellIndx, atmIndx)
    ! Insert an atom into a cell
    use SearchSort, only: BinaryInsert
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: cellIndx, atmIndx
    integer :: nAtoms

    if(atmIndx < 1) then
      error stop "Invalid atom index passed to CellList_InsertAtom"
    endif

    nAtoms = self%nCellAtoms(cellIndx)
    
    ! Simple append for unsorted cell lists
    self%nCellAtoms(cellIndx) = nAtoms + 1
    self%cellList(nAtoms + 1, cellIndx) = atmIndx
    self%cellID(atmIndx) = cellIndx

  end subroutine
!===================================================================================
  subroutine CellList_RemoveAtom(self, atmIndx)
    ! Remove an atom from its current cell
    use SearchSort, only: SimpleSearch
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: atmIndx
    integer :: cellIndx, nAtoms, idx, i

    cellIndx = self%cellID(atmIndx)
    if(cellIndx < 1) return
    
    nAtoms = self%nCellAtoms(cellIndx)
    
    ! Find and remove the atom
    do i = 1, nAtoms
      if(self%cellList(i, cellIndx) == atmIndx) then
        idx = i
        exit
      endif
    enddo
    
    ! Shift remaining atoms
    if(idx < nAtoms) then
      self%cellList(idx:nAtoms-1, cellIndx) = self%cellList(idx+1:nAtoms, cellIndx)
    endif
    self%nCellAtoms(cellIndx) = nAtoms - 1
    self%cellID(atmIndx) = 0

  end subroutine
!===================================================================================
  subroutine CellList_MoveAtom(self, atmIndx, x, y, z)
    ! Update an atom's cell assignment based on new position
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: atmIndx
    real(dp), intent(in) :: x, y, z
    integer :: oldCell, newCell
    integer :: binx, biny, binz
    real(dp) :: boxdim(1:2, 1:3)

    call self%parent%GetDimensions(boxdim)
    
    binx = floor((x - boxdim(1,1))/self%dx)
    biny = floor((y - boxdim(1,2))/self%dy)
    binz = floor((z - boxdim(1,3))/self%dz)
    newCell = self%GetCellIndex(binx, biny, binz)
    
    oldCell = self%cellID(atmIndx)
    
    if(oldCell /= newCell) then
      if(oldCell > 0) call self%RemoveAtom(atmIndx)
      call self%InsertAtom(newCell, atmIndx)
    endif

  end subroutine
!===================================================================================
  subroutine CellList_AddMol(self, disp, tempList, tempNNei)
    implicit none
    class(CellList), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iDisp, iAtom
    integer :: binX, binY, binZ, cellIndx
    real(dp) :: boxdim(1:2, 1:3)

    select type(disp)
      class is(Addition)
        call self%parent%GetDimensions(boxdim)
        do iDisp = 1, size(disp)
           iAtom = disp(iDisp)%atmIndx
           binx = floor((disp(iDisp)%x_new - boxdim(1,1))/self%dx)
           biny = floor((disp(iDisp)%y_new - boxdim(1,2))/self%dy)
           binz = floor((disp(iDisp)%z_new - boxdim(1,3))/self%dz)
           cellIndx = self%GetCellIndex(binx, biny, binz)
           call self%InsertAtom(cellIndx, iAtom)
       enddo
   end select

  end subroutine
!===================================================================================
  subroutine CellList_DeleteMol(self, molIndx, topIndx)
    use Common_MolInfo, only: MolData
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx
    integer :: nStart, nEnd
    integer :: iAtom, nType
    integer :: topStart, topEnd
    integer :: atmIndx, topAtom

    nStart = self%parent%MolStartIndx(molIndx)
    nEnd = self%parent%MolEndIndx(molIndx)
    nType = self%parent%MolType(nStart)
    
    topStart = self%parent%MolStartIndx(topIndx)
    topEnd = self%parent%MolEndIndx(topIndx)

    do iAtom = 1, MolData(nType)%nAtoms
      atmIndx = nEnd - iAtom + 1
      topAtom = topEnd - iAtom + 1
      
      ! Remove deleted atom from cell
      call self%RemoveAtom(atmIndx)
      
      ! If not the same as top atom, update top atom's cell reference
      if(topAtom /= atmIndx) then
        ! Top atom takes deleted atom's place in the array
        self%cellID(atmIndx) = self%cellID(topAtom)
        self%cellID(topAtom) = 0
      endif
    enddo

  end subroutine
!===================================================================================
  subroutine CellList_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(CellList), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   
    real(dp) :: realVal
    character(len=30) :: command 

    lineStat = 0
    call GetXCommand(line, command, 6, lineStat)
    select case( trim(adjustl(command)) )
      case("cellsize")
        call GetXCommand(line, command, 7, lineStat)
        read(command,*) realVal
        self%cellSize = realVal
        self%rCut = realVal
        self%rCutSq = realVal * realVal

      case default
        lineStat = -1
    end select

  end subroutine
!===================================================================================
  subroutine CellList_PrintList(self, writeunit, label)
    use ParallelVar, only: nout
    implicit none
    class(CellList), intent(in) :: self
    integer, intent(in), optional :: writeunit
    character(len=*), intent(in), optional :: label
    integer :: iCell, outunit, j

    if(present(writeunit)) then
      outunit = writeunit
    else
      outunit = nout
    endif

    write(outunit,*) "----------------------------"
    if(present(label)) then
      write(outunit,*) label
    endif
    write(outunit,*) "Cell dimensions:", self%nX, "x", self%nY, "x", self%nZ
    write(outunit,*) "Total cells:", self%nCells
    write(outunit,*) "Cell List:"
    do iCell = 1, self%nCells
      if(self%nCellAtoms(iCell) < 1) cycle
      write(outunit,"(I5,A,1000(I5,1x))") iCell,"|", &
           (self%cellList(j, iCell), j=1,self%nCellAtoms(iCell))
    enddo
    
  end subroutine
!===================================================================================
end module
!===================================================================================

