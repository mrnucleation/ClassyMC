!===================================================================================
! This module contains a Neighborlist which will force
! all pair computations in the system. Obviously not efficient,
! you should only use this if you know what you're doing.
!
! I usually only use this for small cluster simulations where the neighborlist
! doesn't change much. -Troy
!===================================================================================
module NeighAllAtomsDef
use VarPrecision
use CoordinateTypes
use Template_SimBox, only: SimBox
use SimpleSimBox, only: SimpleBox
use Template_NeighList, only: NeighListDef

  type, public, extends(NeighListDef) :: NeighAllAtoms
!      logical :: Sorted = .false.
!      logical :: Strict = .false.
!      integer, allocatable :: list(:,:)
!      integer, allocatable :: nNeigh(:)
!      integer :: maxNei
!      real(dp) :: rCut, rCutSq
!      logical :: restrictType = .false.
!      integer, allocatable :: allowed(:)
!      integer :: safetyCheck = .false.
      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => NeighAllAtoms_Constructor 
      procedure, pass :: BuildList => NeighAllAtoms_BuildList 
      procedure, pass :: SortList => NeighAllAtoms_SortList 
      procedure, pass :: GetNewList => NeighAllAtoms_GetNewList
      procedure, pass :: AddMol => NeighAllAtoms_AddMol
      procedure, pass :: PurgeAtom => NeighAllAtoms_PurgeAtom
      procedure, pass :: SwapAtomLists => NeighAllAtoms_SwapAtomLists
      procedure, pass :: SwapAtomType => NeighAllAtoms_SwapAtomType
      procedure, pass :: GetNeighCount => NeighAllAtoms_GetNeighCount
      procedure, pass :: IntegrityCheck => NeighAllAtoms_IntegrityCheck
      procedure, pass :: ProcessIO => NeighAllAtoms_ProcessIO
!      procedure, pass :: TransferList
      procedure, pass :: PrintList => NeighAllAtoms_PrintList
      procedure, pass :: DeleteMol => NeighAllAtoms_DeleteMol
      procedure, pass :: Prologue => NeighAllAtoms_Prologue
      procedure, pass :: Epilogue => NeighAllAtoms_Epilogue
      procedure, pass :: Update => NeighAllAtoms_Update
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine NeighAllAtoms_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use SimpleSimBox, only: SimpleBox
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 0.5E0_dp  !Used to estimate an approximate volume of 
    integer :: AllocateStatus
    real(dp), pointer :: coords(:,:)




    if(.not. self%initialized) then
        allocate( self%cellID(1:self%parent%nMaxAtoms), stat=AllocateStatus )
        allocate( self%list(1:self%parent%nMaxAtoms, 1:self%parent%nMaxAtoms), stat=AllocateStatus )
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
        error STOP "*** CellNeighRSQList: Memory Allocation Error! ***"
    endif


    self%initialized = .true.
    self%restrictType = .false.
  end subroutine
!===================================================================================
  subroutine NeighAllAtoms_Prologue(self)
    implicit none
    class(NeighAllAtoms), intent(inout) :: self


!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine NeighAllAtoms_BuildList(self, listindx)
    use SearchSort, only: QSort
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
    integer, intent(in) :: listindx

    integer :: iAtom, jNei, jAtom

    self%nNeigh = 0
    self%list = 0
    do iAtom = 1, self%maxAtoms
      if( .not. self%parent%IsActive(iAtom) ) then
        cycle
      endif
      do jAtom = 1, self%maxAtoms
        if( .not. self%parent%IsActive(iAtom) ) then
        if( self%parent%MolIndx(iAtom) == self%parent%MolIndx(jAtom)) cycle
        self%nNeigh(iAtom) = self%nNeigh(iAtom) + 1
        self%list( self%nNeigh(iAtom), iAtom ) = jAtom
      enddo
    enddo

    self%sorted = .true.

  end subroutine
!===================================================================================
! Used to generate a neighborlist for a newly added molecule.  Primarily used
! To give a temporary neighborlist to addition type Monte Carlo Moves.
  subroutine NeighAllAtoms_GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    use Common_MolInfo, only: nMolTypes
    use SearchSort, only: QSort
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
    integer, intent(in) :: iDisp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount
    integer :: jType, jAtom, j, iAtom, jNei
    integer :: jUp, jLow, molIndx, jMol
    real(dp) :: xn, yn, zn

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
      class default
        error stop
    end select

    !Get relevant box information from the parent.

    call self%GetCellAtoms(nCellAtoms, self%atomlist, cellIndx=cellIndx)
    do jAtom = 1, self%maxAtoms
      if(.not. self%parent%IsActive(jAtom) ) cycle
      if(self%parent%MolType(iAtom) == self%parent%MolType(jAtom) ) cycle
      tempNNei(iDisp) = tempNNei(iDisp) + 1
      templist(tempNNei(iDisp), iDisp) = jAtom
    enddo

  end subroutine
!===================================================================================
  subroutine NeighAllAtoms_SwapAtomLists(self, atmindx1, atmindx2)
    !--------------------------------
    ! This function switches
    !--------------------------------

    use SearchSort, only: BinarySearch, SimpleSearch, QSort
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
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



    call self%SortList
    do iSearch = 1, nSearch
      jAtom = searchlist(iSearch)
      neiIndx1 = BinarySearch(atmIndx1, self%list(1:self%nNeigh(jAtom), jAtom))
      neiIndx2 = BinarySearch(atmIndx2, self%list(1:self%nNeigh(jAtom), jAtom))
      if(neiIndx1 /= 0) then
        self%list(neiIndx1, jAtom) = atmIndx2
      endif
      if(neiIndx2 /= 0) then
        self%list(neiIndx2, jAtom) = atmIndx1
      endif     
      call QSort(self%list(1:self%nNeigh(jAtom), jAtom))
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

  end subroutine
!===================================================================================
  subroutine NeighAllAtoms_PurgeAtom(self, atmindx1)
    !--------------------------------
    ! This function removes an atom's entry from
    !--------------------------------

    use SearchSort, only: BinarySearch, SimpleSearch, QSort
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
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
    do jNei = 1, self%nNeigh(atmIndx1)
      jAtom = self%list(jNei, atmIndx1)
      nNeigh = self%nNeigh(jAtom) 
      neiIndx = BinarySearch(atmIndx1, self%list(1:nNeigh, jAtom))
      if(neiIndx == 0) then
        error stop "Neighborlist Error! Index of neighbor atom not found!"
      endif
      if(nNeigh > 1) then
        if(neiIndx /= nNeigh) then
          self%list(neiIndx:nNeigh-1, jAtom) = self%list(neiIndx+1:nNeigh, jAtom)
        endif
        self%nNeigh(jAtom) = self%nNeigh(jAtom) - 1
        nNeigh = nNeigh - 1
        call QSort(self%list(1:nNeigh, jAtom))
      else
        self%nNeigh(jAtom) = 0 
        nNeigh = 0
      endif
    enddo
    self%nNeigh(atmIndx1) = 0

  end subroutine
!===================================================================================
  function NeighAllAtoms_GetNeighCount(self, nAtom, rCount) result(nCount)
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
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
  subroutine NeighAllAtoms_AddMol(self, disp, tempList, tempNNei)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: QSort
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iDisp, iAtom, iNei, nNei, neiIndx, j
    integer :: binX, binY, binZ, cellIndx
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: boxdim(1:2, 1:3)


    select type(disp)
      class is(Addition)
        call self%parent%GetDimensions(boxdim)
        do iDisp = 1, size(disp)
           iAtom = disp(iDisp)%atmIndx
           self%nNeigh(iAtom) = tempNNei(iDisp)
           self%list(1:tempNNei(iDisp), iAtom ) = templist(1:tempNNei(iDisp), iDisp)
           do iNei = 1, tempNNei(iDisp)
               neiIndx = tempList(iNei, iDisp)
               self%nNeigh(neiIndx)= self%nNeigh(neiIndx) + 1
               self % list( self%nNeigh(neiIndx), neiIndx ) = iAtom
           enddo
       enddo
   end select

   call self%sortlist

  end subroutine
!===================================================================================
  subroutine NeighAllAtoms_DeleteMol(self, molIndx, topIndx)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: BinarySearch, SimpleSearch
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
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
  subroutine NeighAllAtoms_SwapAtomType(self, disp, topIndx)
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
    class(NeighAllAtoms), intent(inout) :: self
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
    self % sorted = .false.


  end subroutine
!====================================================================
  subroutine NeighAllAtoms_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetAllCommands, GetXCommand,maxLineLen
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(NeighAllAtoms), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: i, intVal, nPar
    real(dp) :: realVal

    character(len=30) :: command 
    character(len=30), allocatable :: parlist(:)


    lineStat = 0

  end subroutine
!===================================================================================
end module
!===================================================================================
