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
      integer :: nCells = 0
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
      procedure, pass :: GetNewList => CellRSqList_GetNewList
      procedure, pass :: AddMol => CellRSqList_AddMol
      procedure, pass :: GetNeighCount => CellRSqList_GetNeighCount
      procedure, pass :: GetCellIndex => CellRSqList_GetCellIndex
      procedure, pass :: GetBins => CellRSqList_GetBins
      procedure, pass :: GetCellAtoms => CellRSqList_GetCellAtoms
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
    self % maxNei = ceiling(self%rCut**3/atomRadius**3)
 
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
  subroutine CellRSqList_Update(self)
    implicit none
    class(CellRSqList), intent(inout) :: self


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
    integer :: maxAtoms
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
    maxAtoms = self%parent%nMaxAtoms
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
      if(ubound(self%cellList, 2) < self %nCells) then
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
    do iAtom = 1, maxAtoms
      if( self%parent%MolSubIndx(iAtom) > self%parent%NMol(self%parent%MolType(iAtom)) ) then
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
    do iAtom = 1, maxAtoms
      if( self%parent%MolSubIndx(iAtom) > self%parent%NMol(self%parent%MolType(iAtom)) ) then
        cycle
      endif
      call self%GetCellAtoms(nCellAtoms, self%atomlist, atmIndx=iAtom)
      do jNei = 1, nCellAtoms
        jAtom = self%atomlist(jNei)
        if(jAtom <= iAtom) cycle
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
      if( self%parent%MolSubIndx(iAtom) > self%parent%NMol(self%parent%MolType(iAtom)) ) then
        cycle
      endif
      nNeigh = self%nNeigh(iAtom)
      call QSort(self%list(1:nNeigh, iAtom))
    enddo

    self%sorted = .false.

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
        if( self%parent%MolSubIndx(jAtom) == molIndx ) then
          cycle
        endif
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
    implicit none
    class(CellRSqList), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iDisp, iAtom, iNei, nNei, neiIndx, j
    integer :: binX, binY, binZ, cellIndx
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: boxdim(1:2, 1:3)


    select type(disp)
      class is(Addition)

!        write(2,*) "----------------------------"
!        write(2,*) "Add"
!        do iAtom = 1, trialBox%nMaxatoms
!          write(2,"(I3,A,1000(I3))") iAtom,"|", (trialBox % NeighList(iList)%list(j, iAtom) ,j=1,trialBox % NeighList(iList)%nNeigh(iAtom))
!1        enddo
!        write(2,*)

    call self%parent%GetDimensions(boxdim)
    do iDisp = 1, size(disp)
!          write(2,"(A,A,1000(I3))") "NewList","|", (tempList(j, iDisp) ,j=1,tempNNei(iDisp))
      iAtom = disp(iDisp)%atmIndx
      binx = floor((disp(iDisp)%x_new - boxdim(1,1))/self%dx)
      biny = floor((disp(iDisp)%y_new - boxdim(1,2))/self%dy)
      binz = floor((disp(iDisp)%z_new - boxdim(1,3))/self%dz)
      cellIndx = self % GetCellIndex(binx, biny, binz)
      self % cellID(iAtom) = cellIndx
      self % nCellAtoms(cellIndx) = self % nCellAtoms(cellIndx) + 1
      self % cellList(self%nCellAtoms(cellIndx), cellIndx) = iAtom

      self%nNeigh(iAtom) = tempNNei(iDisp)
      do iNei = 1, tempNNei(iDisp)
        neiIndx = tempList(iNei, iDisp)
        self % list(iNei, iAtom) =  neiIndx
        self % list( self%nNeigh(neiIndx)+1, neiIndx ) = iAtom
        self%nNeigh(neiIndx)= self%nNeigh(neiIndx) + 1
      enddo
    enddo
   end select
!        do iAtom = 1, trialBox%nMaxatoms
!          write(2,"(I3,A,1000(I3))") iAtom,"|", (trialBox % NeighList(iList)%list(j, iAtom) ,j=1,trialBox % NeighList(iList)%nNeigh(iAtom))
!        enddo
!        write(2,*)
!      write(2,*) "N", trialBox%NeighList(iList)%nNeigh(:)     



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


      !Delete atom's the entry from the cell list.
      cellIndx = self%cellID(atmIndx)
      nNei = self%nCellAtoms(cellIndx)
      curIndx = SimpleSearch( atmIndx, self%celllist(1:nNei, cellIndx) )
      self%celllist(1:nNei, cellIndx ) = [self%celllist(1:curIndx-1, cellIndx), &
                                          self%celllist(curIndx+1:nNei, cellIndx) ]
      self%nCellAtoms(cellIndx) = self%nCellAtoms(cellIndx) - 1

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



    enddo

!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "---------------------------------"


    self % sorted = .false.


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
