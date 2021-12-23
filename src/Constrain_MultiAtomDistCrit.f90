!====================================================================
!This module contains the Stilinger distance criteria that is used to 
!enforce clustering. 
!====================================================================
module Constrain_MultiAtomDistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use CoordinateTypes, only: Displacement, Deletion, Addition
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout

  type, public, extends(constraint) :: MultiAtomDistCrit
    integer :: neighList = 1
    integer :: molType = 1
    real(dp) :: rCut, rCutSq
    integer :: boxID = 1

    integer :: nAtomMax
    integer :: nMolMax
    integer, allocatable :: TypeStart(:), TypeEnd(:)
    logical, allocatable :: clustMemb(:)

    integer, allocatable :: nTopoNei(:)
    integer, allocatable :: topoList(:, :)
    integer, allocatable :: nTopoNewNei(:)
    integer, allocatable :: newTopoList(:, :)

    integer, allocatable :: newList(:)
    integer, allocatable :: newList2(:)
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => MultiDistCrit_Constructor
      procedure, pass :: CheckInitialConstraint => MultiDistCrit_CheckInitialConstraint
      procedure, pass :: DiffCheck => MultiDistCrit_DiffCheck
      procedure, pass :: PurgeMolIndex => MultiDistCrit_PurgeMolIndex
      procedure, pass :: SortList => MultiDistCrit_SortList
      procedure, pass :: SwapMolLists => MultiDistCrit_SwapMolLists
      procedure, pass :: DeleteMol => MultiDistCrit_DeleteMol
      procedure, pass :: ProcessIO => MultiDistCrit_ProcessIO
      procedure, pass :: Maintenance => MultiDistCrit_Maintenance
      procedure, pass :: Update => MultiDistCrit_Update
      procedure, pass :: Epilogue => MultiDistCrit_Epilogue
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine MultiDistCrit_Constructor(self, boxID)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: AllocateStat
    integer :: nMolMax
    integer :: iType

    self%boxID = boxID
    self%parent => BoxArray(boxID) % box 

!    nMolMax = self % parent % NMolMax(self%molType)
    nMolMax = self % parent % maxMol
    self%nAtomMax = 1

    allocate(self%TypeStart(1:nMolTypes), stat = AllocateStat)
    allocate(self%TypeEnd(1:nMolTypes), stat = AllocateStat)
    self%TypeStart = 0
    self%TypeEnd = 0
    self%TypeStart(1) = 1
    self%TypeEnd(1) = self%parent%NMolMax(1)
    if(nMolTypes > 1) then
      do iType = 2, nMolTypes
        self%TypeStart(iType) = self%TypeEnd(iType-1) + 1
        self%TypeEnd(iType) = self%TypeEnd(iType-1) + self%parent%NMolMax(iType) 
      enddo
    endif

    do iType = 1, nMolTypes
      self%nAtomMax = max(self%nAtomMax, MolData(iType)%nAtoms)
    enddo
    
    allocate(self%clustMemb(1:nMolMax), stat = AllocateStat)

    allocate(self%nTopoNei(1:nMolMax), stat = AllocateStat)
    allocate(self%topoList(1:nMolMax, 1:nMolMax), stat = AllocateStat)
    allocate(self%nTopoNewNei(1:nMolMax), stat = AllocateStat)
    allocate(self%newTopoList(1:nMolMax, 1:nMolMax), stat = AllocateStat)

    allocate(self%newlist(1:nMolMax), stat = AllocateStat)
    allocate(self%newlist2(1:nMolMax), stat = AllocateStat)

    self%nMolMax = nMolMax

    IF (AllocateStat /= 0) STOP "Allocation Error in Distance Constraint"
  end subroutine
!=====================================================================
  subroutine MultiDistCrit_CheckInitialConstraint(self, trialBox, accept)
    use Common_MolInfo, only: MolData
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    integer :: totalMol, nNew, nClust, neiIndx
    integer :: iMol,jMol, iAtom, iAtoms, nIAtoms, nJAtoms, jAtom, iLimit
    integer :: molIndx, molType, nNext, nNext2, iNext
    integer :: newNNei, iStart, iEnd, jStart, jEnd, jNei
    real(dp) :: rx, ry, rz, rsq

    real(dp), pointer :: atoms(:,:) => null()

    accept = .true.
    self%clustMemb = .false.
    self%topoList = 0
    self%nTopoNei = 0
    
    call trialbox%GetCoordinates(atoms)

    totalMol = self%nMolMax
    if(trialBox%nMolTotal < 2) then
      accept = .true.
      return
    endif


    do iMol = 1, totalMol-1
      call trialbox%GetMolData(iMol, nAtoms=nIAtoms, molStart=iStart, molEnd=iEnd)
      if(.not. trialbox%IsActive(iStart) ) cycle

      do jMol = iMol+1, totalMol
        call trialbox%GetMolData(jMol, nAtoms=nJAtoms, molStart=jStart, molEnd=jEnd)
        if(.not. trialbox%IsActive(jStart) ) cycle
        atomloop: do iAtom = iStart, iEnd
          do jAtom = jStart, jEnd
            rx = atoms(1, iAtom) - atoms(1, jAtom)
            ry = atoms(2, iAtom) - atoms(2, jAtom)
            rz = atoms(3, iAtom) - atoms(3, jAtom)
            call trialBox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq ) then
              newNNei = self%nTopoNei(iMol) + 1
              self%nTopoNei(iMol) = newNNei
              self%topolist(newNnei, iMol) = jMol

              newNNei = self%nTopoNei(jMol) + 1
              self%nTopoNei(jMol) = newNNei
              self%topolist(newNnei, jMol) = iMol
              exit atomloop
            endif
          enddo
        enddo atomloop
      enddo
    enddo


!    do iMol = 1, totalMol
!      if(.not. trialBox%IsActive(iMol)) cycle
!      newNNei = self%nTopoNei(iMol)
!      write(*,*) iMol, "|", self%topolist(1:newNNei, iMol)
!    enddo
!
!    stop



     !Seed the initial cluter check by adding the first particle in the array
     !to the cluster
    self%clustMemb = .false.
    nClust = 0
    do iMol = 1, totalMol
      if(self%nTopoNewNei(iMol) > 0) then
        self%clustMemb(iMol) = .true.
        nNew = 1
        nClust = 1
        self%newlist(1) = iMol
        exit
      endif
    enddo

    
    do iLimit = 1, totalMol
      nNext = nNew
      nNew = 0
      do iNext = 1, nNext
        iMol = self%newlist(iNext)
        do jNei = 1, self%nTopoNei(iMol)
          jMol = self%topolist(jNei, iMol)
          if(.not. self%clustMemb(jMol) )then
            self%clustMemb(jMol) = .true.
            nNew = nNew + 1
            nClust = nClust + 1
            self%newlist2(nNew) = jMol
          endif
        enddo
      enddo

       !If every molecule has been added then no further calculations are
       !needed.
      if(nClust >= totalMol) then
        exit
      endif 
       ! If no new molecules were added, the algorithm has hit a dead end
       ! and the cluster is broken.  
      if(nNew <= 0) then
        exit
      endif

      self%newlist(1:nNew) = self%newlist2(1:nNew)
    enddo


     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
    if( nClust < trialBox%nMolTotal ) then
      accept = .false.
      write(nout,*) "Detailed Cluster Criteria Check Failed!"
      return
    endif

    write(nout,*) "Detailed Cluster Criteria Check Succeeded!"
    accept = .true.

  end subroutine
!=============================================================
  subroutine MultiDistCrit_DiffCheck(self, trialBox, disp, accept)
    use SearchSort, only: QSort
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept

    integer :: iDisp
    integer :: totalMol, nNew, nClust, neiIndx, startMol, topMol
    integer :: nAtoms, newNnei, atmIndx
    integer :: molStart, molEnd
    integer :: iStart, iEnd, jStart, jEnd, jNei
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    integer :: molIndx, molType, nNext, nNext2, iNext
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: temppos(1:3, self%nAtomMax)
    real(dp), pointer :: atoms(:,:) => null()

    accept = .true.
    totalMol = trialBox%nMolTotal

    call trialbox%GetCoordinates(atoms)
    !This section creates the topology list of the new state using information
    !based on the what kind of perturbation was performed.

    self%newTopoList = self%topoList 
    self%nTopoNewNei = self%nTopoNei


    select type(disp)
       !----------------------------------------------------------------------------
      class is(Displacement)
!        write(*,*) "Shift"
        if(totalMol == 1) then
          accept = .true.
          return
        endif
        molIndx = disp(1)%molIndx
!        write(*,*) "Shift", molIndx
        call trialbox%GetMolData(molIndx, molType=molType, nAtoms=nAtoms, molStart=molStart, molEnd=molEnd)
        topMol = trialBox % NMol(molType) + self%TypeStart(molType) - 1
        call self%PurgeMolIndex(molIndx, self%newTopoList, self%nTopoNewNei)
        temppos(1:3, 1:nAtoms) = atoms(1:3, molStart:molEnd)
        do iDisp = 1, size(disp)
          atmIndx = disp(iDisp)%atmIndx - molStart + 1
          temppos(1, atmIndx) = disp(iDisp)%x_new
          temppos(2, atmIndx) = disp(iDisp)%y_new
          temppos(3, atmIndx) = disp(iDisp)%z_new
        enddo

        do jMol = 1, self%nMolMax
          call trialbox%GetMolData(jMol, molStart=jStart, molEnd=jEnd)
          if(.not. trialbox%IsActive(jStart) ) cycle
          if(molIndx == jMol) cycle
          atomloop: do iAtom = 1, nAtoms
            do jAtom = jStart, jEnd
              rx = temppos(1, iAtom) - atoms(1, jAtom)
              ry = temppos(2, iAtom) - atoms(2, jAtom)
              rz = temppos(3, iAtom) - atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
              if(rsq < self%rCutSq ) then
                newNNei = self%nTopoNewNei(molIndx) + 1
                self%nTopoNewNei(molIndx) = newNNei
                self%newTopoList(newNnei, molIndx) = jMol


                newNNei = self%nTopoNewNei(jMol) + 1
                self%nTopoNewNei(jMol) = newNNei
                self%newTopoList(newNnei, jMol) = molIndx
!                call QSort(self%newTopoList(1:newNNei, jMol))
                exit atomloop
              endif
            enddo
          enddo atomloop
        enddo


        ! Early Rejection by checking if any molecule has zero neighbors.
        do jMol = 1, self%nMolMax
          call trialbox%GetMolData(jMol, molStart=jStart)
          if(.not. trialBox%IsActive(jStart) ) cycle
          if(self%nTopoNewNei(jMol) < 1) then
            accept = .false.
            return
          endif
        enddo
        newNNei = self%nTopoNewNei(molIndx) 
!        call QSort(self%newTopoList(1:newNnei, molIndx))


!        call self%sortlist(self%newTopoList, self%nTopoNewNei)

       !----------------------------------------------------------------------------
      class is(Addition)
!        write(*,*) "Add"
        molIndx = disp(1)%molIndx
        call trialbox%GetMolData(molIndx, nAtoms=nAtoms, molStart=molStart)
        temppos(1:3, 1:nAtoms) = 0E0_dp
        self%nTopoNewNei(molIndx) = 0 
!        write(*,*) "Add", molIndx, nAtoms
        do iDisp = 1, size(disp)
          atmIndx = disp(iDisp)%atmIndx - molStart + 1
          temppos(1, atmIndx) = disp(iDisp)%x_new
          temppos(2, atmIndx) = disp(iDisp)%y_new
          temppos(3, atmIndx) = disp(iDisp)%z_new
        enddo

!        write(*,*) self%nMolMax
        do jMol = 1, self%nMolMax
          if(molIndx == jMol) cycle
          call trialbox%GetMolData(jMol, molStart=jStart, molEnd=jEnd)
!          write(*,*) jMol, nAtoms, trialBox%IsActive(jMol)
!          write(*,*) jEnd-jStart+1
!          write(*,*)
          if(.not. trialbox%IsActive(jStart) ) cycle
          atomloopadd: do iAtom = 1, nAtoms
            do jAtom = jStart, jEnd
              rx = temppos(1, iAtom) - atoms(1, jAtom)
              ry = temppos(2, iAtom) - atoms(2, jAtom)
              rz = temppos(3, iAtom) - atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
!              write(*,*) iAtom, jAtom, sqrt(rsq), self%rCut
              if(rsq < self%rCutSq ) then
!                write(*,*) molIndx, jMol
                newNNei = self%nTopoNewNei(molIndx) + 1
                self%nTopoNewNei(molIndx) = newNNei
                self%newTopoList(newNnei, molIndx) = jMol
!                write(*,*) newNNei,",",self%newTopoList(newNnei, molIndx)

                newNNei = self%nTopoNewNei(jMol) + 1
                self%nTopoNewNei(jMol) = newNNei
                self%newTopoList(newNnei, jMol) = molIndx

!                call QSort(self%newTopoList(1:newNNei, jMol))
!                write(*,*) newNNei,",", self%newTopoList(newNnei, jMol)
                exit atomloopadd
              endif
            enddo
          enddo atomloopadd
        enddo

        if(self%nTopoNewNei(molIndx) < 1) then
          accept = .false.
!          write(*,*) accept
          return
        endif
        totalMol = totalMol + 1
        newNNei = self%nTopoNewNei(molIndx) 
!        call QSort(self%newTopoList(1:newNnei, molIndx))
!        call self%sortlist(self%newTopoList, self%nTopoNewNei)
       !----------------------------------------------------------------------------
      class is(Deletion)
!        write(*,*) "Delete"
        molType = disp(1)%molType
        molIndx = disp(1)%molIndx
        topMol = trialBox % NMol(molType) + self%TypeStart(molType) - 1

        if(totalMol == 2) then
          self%newTopoList = 0
          self%nTopoNewNei = 0
          accept = .true.
          return
        endif
        call self%DeleteMol(molIndx, topMol, self%newTopoList, self%nTopoNewNei)

        totalMol = totalMol - 1
       !----------------------------------------------------------------------------
      class default
        stop "Distance criteria is not compatiable with this perturbation type."
       !----------------------------------------------------------------------------
    end select

!    write(*,*) "NTotal:", totalMol
!    do iMol = 1, self%nMolMax
!      newNNei = self%nTopoNewNei(iMol)
!      write(*,*) iMol, "|", self%newTopoList(1:newNNei, iMol)
!    enddo

    !Find a starting point.
    self%clustMemb = .false.
    nClust = 0
    do iMol = 1, self%nMolMax
      if(self%nTopoNewNei(iMol) > 0) then
        self%clustMemb(iMol) = .true.
        nNew = 1
        nClust = 1
        self%newlist(1) = iMol
        exit
      endif
    enddo

    if(nClust < 1) then
      accept = .false.
!      write(*,*) accept
      return
    endif

    !Using the newly constructed topology list, check to see if the new cluster satisfies
    !the cluster criteria. 
    do iLimit = 1, totalMol
      nNext = nNew
      nNew = 0
      do iNext = 1, nNext
        iMol = self%newlist(iNext)
        if(self%nTopoNewNei(iMol) < 1) cycle
        do jNei = 1, self%nTopoNewNei(iMol)
          jMol = self%newTopoList(jNei, iMol)
          if(.not. self%clustMemb(jMol) )then
            self%clustMemb(jMol) = .true.
            nNew = nNew + 1
            nClust = nClust + 1
            self%newlist2(nNew) = jMol
          endif
        enddo
      enddo

       !If every molecule has been added then no further calculations are
       !needed.
      if(nClust >= totalMol) then
        exit
      endif 
       ! If no new molecules were added, the algorithm has hit a dead end
       ! and the cluster is broken.  
      if(nNew <= 0) then
        exit
      endif

      self%newlist(1:nNew) = self%newlist2(1:nNew)
    enddo

     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
!    write(*,*) nClust, totalMol
    if( nClust < totalMol ) then
      accept = .false.
    else
      accept = .true.
    endif

!    write(*,*) accept
!    write(*,*)



  end subroutine
!===================================================================================
  subroutine MultiDistCrit_SortList(self, topolist, ntopo)
    use SearchSort, only: QSort
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    integer, intent(inout) :: topolist(:,:)
    integer, intent(inout) :: ntopo(:)
    integer :: iMol
    integer :: nNeigh

    do iMol = 1, size(topolist,2)
      nNeigh = ntopo(iMol)
      if( ntopo(iMol) < 1 ) then
        cycle
      endif
!      call QSort( topolist(1:nNeigh, iMol) )
    enddo


  end subroutine
!===================================================================================
  subroutine MultiDistCrit_SwapMolLists(self, topolist, nTopo, molindx1, molindx2)
    use SearchSort, only: BinarySearch, SimpleSearch, QSort
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    integer, intent(in) :: molindx1, molindx2
    integer, intent(inout) :: topolist(:,:)
    integer, intent(inout) :: ntopo(:)
    integer :: nNeigh1, nNeigh2
    integer :: templist(1:2*self%nMolMax)
    integer :: jMol, iNei, jNei

    integer :: iSearch
    integer :: nSearch = 0
    integer :: neiIndx1, neiIndx2
    integer :: searchlist(1:2*self%nMolMax)


    nNeigh1 = nTopo(molIndx1)
    nNeigh2 = nTopo(molIndx2)

    !If neither atom has a neighbor, nothing needs to be done. 
    if( (nNeigh1 == 0) .and. (nNeigh2 == 0) ) then
      return
    endif

    !The first task is to get a list of neighbors of both atoms in 
    !order to figure out which lists need to be modified. To prevent double
    !counting a unique list is computed
    templist(1:nNeigh1+nNeigh2) = [topolist(1:nNeigh1, molIndx1), &
                                   topolist(1:nNeigh2, molIndx2)]
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
!    call self%SortList(topolist,nTopo)
    do iSearch = 1, nSearch
      jMol = searchlist(iSearch)
!      neiIndx1 = BinarySearch(molIndx1, topolist(1:nTopo(jMol), jMol))
!      neiIndx2 = BinarySearch(molIndx2, topolist(1:nTopo(jMol), jMol))
      neiIndx1 = SimpleSearch(molIndx1, topolist(1:nTopo(jMol), jMol))
      neiIndx2 = SimpleSearch(molIndx2, topolist(1:nTopo(jMol), jMol))
      if(neiIndx1 /= 0) then
        topolist(neiIndx1, jMol) = molIndx2
      endif
      if(neiIndx2 /= 0) then
        topolist(neiIndx2, jMol) = molIndx1
      endif     
!      call QSort(topolist(1:nTopo(jMol), jMol))
    enddo


    !Now that the other lists are reorganized, we need to swap the two lists
    if( nNeigh1 == 0 ) then
      topolist(1:nNeigh2, molIndx1) = topolist(1:nNeigh2, molIndx2)
    elseif( nNeigh2 == 0 ) then
      topolist(1:nNeigh1, molIndx2) = topolist(1:nNeigh1, molIndx1)
    else
      templist(1:nNeigh1) = topolist(1:nNeigh1, molIndx1)
      topolist(1:nNeigh2, molIndx1) = topolist(1:nNeigh2, molIndx2)
      topolist(1:nNeigh1, molIndx2) = templist(1:nNeigh1)
    endif


    !Swap the neighbor counters.
    nTopo(molIndx1) = nNeigh2
    nTopo(molIndx2) = nNeigh1

  end subroutine
!===================================================================================
  subroutine MultiDistCrit_PurgeMolIndex(self, molIndx, topolist, nTopo)
    ! Searches for the molIndx value given inside of the topolist
    ! and then removes any entry of this index within the rows of the entire list
    use SearchSort, only: BinarySearch, SimpleSearch
    use SearchSort, only: QSort
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    integer, intent(in) :: molIndx
    integer, intent(inout) :: topolist(:,:)
    integer, intent(inout) :: ntopo(:)

    integer :: iMol, nNei, neighIndx

!    call self%sortlist(topolist, nTopo)


    do iMol = 1, self%nMolMax
      if(iMol == molIndx) cycle
      if(nTopo(iMol) < 1) cycle


!      neighIndx = BinarySearch(molIndx, topolist(1:nTopo(iMol), iMol))
      neighIndx = SimpleSearch(molIndx, topolist(1:nTopo(iMol), iMol))
      if(neighIndx == 0) cycle

      nNei = nTopo(iMol)
      if(nNei > 1 ) then
        if(neighIndx /= nTopo(iMol)) then
          topolist(neighIndx:nNei-1, iMol) = topolist(neighIndx+1:nNei, iMol)
        endif
        nTopo(iMol) = nTopo(iMol) - 1
      else
        topolist(:, iMol) = 0
        nTopo(iMol) = 0
      endif
    enddo
    topolist(:,molIndx) = 0
    nTopo(molIndx) = 0

  end subroutine
!===================================================================================
  subroutine MultiDistCrit_DeleteMol(self, molIndx, topMol, topolist, nTopo)
    ! Removes all traces of molIndx and re-indexes the arrays
    ! Such that the top molecule
    use SearchSort, only: BinarySearch, SimpleSearch
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    integer, intent(in) :: molIndx, topMol
    integer, intent(inout) :: topolist(:,:)
    integer, intent(inout) :: ntopo(:)

    integer :: molType
    integer :: iMol, nNei, neighIndx

!    call self%sortlist(topolist, nTopo)
    if(topMol /= molIndx) then
      call self%SwapMolLists(topolist, nTopo,  molIndx, topMol)
    endif
    call self%PurgeMolIndex(topMol, topolist, nTopo)
!    call self%sortlist(topolist, nTopo)

  end subroutine
!=============================================================
  subroutine MultiDistCrit_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ParallelVar, only: nout
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) realVal
    self%rCut = realVal
    self%rCutSq = realVal*realVal
    write(nout, *) "MultiAtom Distance Criteria:",self%rCut
    write(nout, *) "MultiAtom Distance Criteria (SQ):",self%rCutSq

  end subroutine
!====================================================================
  subroutine MultiDistCrit_Maintenance(self)
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self

  end subroutine
!====================================================================
  subroutine MultiDistCrit_Epilogue(self)
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    logical :: accept

    call self % CheckInitialConstraint(self%parent, accept)

  end subroutine
!=============================================================
  subroutine MultiDistCrit_Update(self)
    implicit none
    class(MultiAtomDistCrit), intent(inout) :: self
    logical :: accept
    integer :: iMol, newNNei

    self%topoList = self%newTopoList
    self%nTopoNei = self%nTopoNewNei
!    write(*,*) "Jesus save me from this code!"
!    do iMol = 1, self%nMolMax
!      newNNei = self%nTopoNewNei(iMol)
!      write(*,*) iMol, "|", self%newTopoList(1:newNNei, iMol)
!    enddo
!    write!(*,*)
!    call self % CheckInitialConstraint(self%parent, accept)
!    if(.not. accept) then
!    do iMol = 1, self%nMolMax
!      newNNei = self%nTopoNei(iMol)
!      write(*,*) iMol, "|", self%topolist(1:newNNei, iMol)
!    enddo
!
!      stop
!    endif
  end subroutine
!=====================================================================
end module
!=====================================================================
