!====================================================================
! Stillinger distance criteria over all molecule types in the box.
! Input (per type i): reference atom index, then rCut(i).
! Pair cutoff uses arithmetic mixing: r_ij = (r_i + r_j) / 2
!====================================================================
module Constrain_DistanceCriteriaAllMol
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use CoordinateTypes, only: Displacement, Deletion, Addition, NewState_IsoMol
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout
  implicit none

  type, public, extends(constraint) :: DistCriteriaAllMol
    integer :: boxID = 1
    integer :: nMolMax = 0

    integer, allocatable :: atomNumType(:)
    real(dp), allocatable :: rCut(:)
    real(dp), allocatable :: rCutPairSq(:,:)

    logical, allocatable :: topoList(:,:)
    logical, allocatable :: newTopoList(:,:)
    logical, allocatable :: clustMemb(:)
    integer, allocatable :: newlist(:) ! Used to create the queue for BFS
    integer, allocatable :: newlist2(:)

    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => DistAllMol_Constructor
      procedure, pass :: CheckInitialConstraint => DistAllMol_CheckInitialConstraint
      procedure, pass :: DiffCheck => DistAllMol_DiffCheck
      procedure, pass :: ProcessIO => DistAllMol_ProcessIO
      procedure, pass :: Maintenance => DistAllMol_Maintenance
      procedure, pass :: Update => DistAllMol_Update
      procedure, pass :: Epilogue => DistAllMol_Epilogue
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine DistAllMol_Constructor(self, boxID)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: AllocateStat

    self%boxID = boxID
    self%parent => BoxArray(boxID) % box

    self%nMolMax = self % parent % maxMol

    if(.not. allocated(self%rCut)) then
      allocate(self%rCut(1:nMolTypes), stat = AllocateStat)
      allocate(self%rCutPairSq(1:nMolTypes, 1:nMolTypes), stat = AllocateStat)
      allocate(self%atomNumType(1:nMolTypes), stat = AllocateStat)
      self%rCut = 0E0_dp
      self%rCutPairSq = 0E0_dp
      self%atomNumType = 1
    endif

    allocate(self%topoList(1:self%nMolMax, 1:self%nMolMax), stat = AllocateStat)
    allocate(self%newTopoList(1:self%nMolMax, 1:self%nMolMax), stat = AllocateStat)
    allocate(self%clustMemb(1:self%nMolMax), stat = AllocateStat)
    allocate(self%newlist(1:self%nMolMax), stat = AllocateStat)
    allocate(self%newlist2(1:self%nMolMax), stat = AllocateStat)

    self%topoList = .false.
    self%newTopoList = .false.
    self%clustMemb = .false.

    IF (AllocateStat /= 0) ERROR STOP "Allocation Error in All-Mol Distance Constraint Constructor"
  end subroutine
!=====================================================================
  subroutine DistAllMol_RebuildPairCutSq(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    integer :: iType, jType
    real(dp) :: rMix

    do iType = 1, nMolTypes
      do jType = 1, nMolTypes
        rMix = 0.5E0_dp * (self%rCut(iType) + self%rCut(jType))
        self%rCutPairSq(iType, jType) = rMix * rMix
      enddo
    enddo
  end subroutine
!=====================================================================
  subroutine DistAllMol_CheckInitialConstraint(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    integer :: nActive, nNew, nClust
    integer :: iMol, jMol, iAtom, jAtom, iLimit
    integer :: iType, jType, iStart, jStart, nNext, iNext
    integer :: startMol
    real(dp) :: rx, ry, rz, rsq
    logical :: leave

    ! Initialize the variables
    accept = .true.
    self%topoList = .false.
    self%clustMemb = .false.


    ! First check if we even have more than one molecule.  If we don't then just skip this. 
    nActive = trialBox%nMolTotal
    if(nActive < 2) then
      return
    endif

    do iMol = 1, self%nMolMax - 1
      !Get the molecule type and first atom index for the first molecule. Then see if this molecule is active.
      call trialBox%GetMolData(iMol, molType=iType, molStart=iStart)
      if(.not. trialBox%IsActive(iStart)) cycle
      !Each molecule has a different designated atom that is used. This computes that index. 
      iAtom = iStart + self%atomNumType(iType) - 1

      ! Now loop over the pairs from i+1 to N
      do jMol = iMol + 1, self%nMolMax

        ! Get the molecule information for the j-molecule. 
        call trialBox%GetMolData(jMol, molType=jType, molStart=jStart)
        if(.not. trialBox%IsActive(jStart)) cycle
        jAtom = jStart + self%atomNumType(jType) - 1

        !Compute distances with the key atoms. 
        rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
        ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
        rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
        call trialBox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz

        ! If the two are within the cutoff, they are neighbored and have a graph connection.
        if(rsq < self%rCutPairSq(iType, jType)) then
          self%topoList(iMol, jMol) = .true.
          self%topoList(jMol, iMol) = .true.
        endif
      enddo
    enddo

    ! If we have zero connections despite having more than one molecule, there's no cluster.
    if(all(self%topoList .eqv. .false.)) then
      accept = .false.
      write(nout,*) "All-Mol cluster criteria check failed!"
      return
    endif


    !Find the first molecule that is active in the system.  Once we have one, exit the loop.
    startMol = 0
    leave = .false.
    do iMol = 1, self%nMolMax
      do jMol = 1, self%nMolMax
        if(self%topoList(jMol, iMol)) then
          startMol = iMol
          leave = .true.
          exit
        endif
      enddo
      if(leave) exit
    enddo

    ! Initialize the breadth first search by adding the startMol to the cluster. Then add it to the schedule.
    self%clustMemb(startMol) = .true.
    self%newlist(1) = startMol
    nNew = 1
    nClust = 1

    ! We now grow the cluster by expanding the graph using breadth first search. The loop only needs
    ! to be done N_mol times since that's the worst case scenario of a valid cluster.  
    do iLimit = 1, nActive
      nNext = nNew
      nNew = 0
      do iNext = 1, nNext
        ! Get the next molecule from the queue.
        iMol = self%newlist(iNext)
        ! Loop over all molecules and check if they are connected to the current molecule.
        do jMol = 1, self%nMolMax
          ! If the two are connected, check next if it's already in the cluster. If not, add it to the queue.
          if(self%topoList(jMol, iMol)) then
            if(.not. self%clustMemb(jMol)) then
              self%clustMemb(jMol) = .true.
              ! Add the molecule to the queue and increment the counters.
              nNew = nNew + 1
              nClust = nClust + 1
              self%newlist2(nNew) = jMol
            endif
          endif
        enddo
      enddo

      ! If we have found everything, no need to continue! We're already done.
      if(nClust >= nActive) then
        exit
      endif

      ! If we haven't added anything new, exit the cluster.  We either succeeded or failed to find everything
      if(nNew <= 0) then
        exit
      endif

      self%newlist(1:nNew) = self%newlist2(1:nNew)
    enddo

    ! If we didn't find everything, or we didn't add anything new, the cluster is invalid.
    if((nNew <= 0) .or. (nClust < nActive)) then
      accept = .false.
      write(nout,*) "All-Mol cluster criteria check failed!"
      return
    endif

    !write(nout,*) "All-Mol cluster criteria check succeeded!"
    accept = .true.

  end subroutine
!=============================================================
  subroutine DistAllMol_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    logical :: leave
    integer :: iDisp
    integer :: nActive, nNew, nClust, nNext, iNext, iLimit
    integer :: iMol, jMol, iAtom, jAtom
    integer :: molIndx, molType, molStart, jType, jStart
    integer :: startMol, topIndex
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    nActive = trialBox%nMolTotal

    select type(disp)
      !----------------------------------------------------------------
      class is(Displacement)
        ! Initialize the new topology list to the current topology list.
        self%newTopoList = self%topoList
        accept = .true.

        do iDisp = 1, size(disp)
          molIndx = disp(iDisp)%molIndx
          call trialBox%GetMolData(molIndx, molType=molType, molStart=molStart)
          iAtom = molStart + self%atomNumType(molType) - 1

          ! If the atom being displaced is the key atom, we need to check the new connections.
          if(disp(iDisp)%atmIndx == iAtom) then
            accept = .false.
            iMol = molIndx
            ! Remove all connections from the moved molecule to all other molecules.
            do jMol = 1, self%nMolMax
              if(iMol /= jMol) then
                self%newTopoList(jMol, iMol) = .false.
                self%newTopoList(iMol, jMol) = .false.
              endif
            enddo
            ! Loop over all molecules and check if they are within the cutoff. If they are, add a connection.
            do jMol = 1, self%nMolMax
              if(iMol == jMol) cycle
              call trialBox%GetMolData(jMol, molType=jType, molStart=jStart)
              if(.not. trialBox%IsActive(jStart)) cycle
              jAtom = jStart + self%atomNumType(jType) - 1

              ! Compute the distance between the two atoms.
              rx = disp(iDisp)%x_new - trialBox%atoms(1, jAtom)
              ry = disp(iDisp)%y_new - trialBox%atoms(2, jAtom)
              rz = disp(iDisp)%z_new - trialBox%atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
              ! If the distance is within the cutoff, add a connection.
              if(rsq < self%rCutPairSq(molType, jType)) then
                self%newTopoList(jMol, iMol) = .true.
                self%newTopoList(iMol, jMol) = .true.
              endif
            enddo
          endif
        enddo

        ! This triggers when the atom being displaced is not the key atom.
        if(accept) then
          return
        endif

        ! Initialize the breadth first search by adding the startMol to the cluster. Then add it to the schedule.
        self%clustMemb = .false.
        startMol = 0
        leave = .false.
        do iMol = 1, self%nMolMax
          do jMol = 1, self%nMolMax
            if(self%newTopoList(jMol, iMol)) then
              startMol = iMol
              leave = .true.
              exit
            endif
          enddo
          if(leave) exit
        enddo
        if(startMol == 0) then
          if(nActive >= 2) then
            accept = .false.
            return
          endif
          nNew = 0
          nClust = nActive
        else
          self%clustMemb(startMol) = .true.
          self%newlist(1) = startMol
          nNew = 1
          nClust = 1
        endif

      !----------------------------------------------------------------
      class is(Addition)
        self%newTopoList = self%topoList
        accept = .true.
        molIndx = disp(1)%molIndx
        call trialBox%GetMolData(molIndx, molType=molType, molStart=molStart)
        iAtom = molStart + self%atomNumType(molType) - 1
        do iDisp = 1, size(disp)
          if(disp(iDisp)%atmIndx == iAtom) then
            accept = .false.
            iMol = molIndx
            do jMol = 1, self%nMolMax
              if(iMol == jMol) cycle
              call trialBox%GetMolData(jMol, molType=jType, molStart=jStart)
              if(.not. trialBox%IsActive(jStart)) cycle
              jAtom = jStart + self%atomNumType(jType) - 1
              rx = disp(iDisp)%x_new - trialBox%atoms(1, jAtom)
              ry = disp(iDisp)%y_new - trialBox%atoms(2, jAtom)
              rz = disp(iDisp)%z_new - trialBox%atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
              if(rsq < self%rCutPairSq(molType, jType)) then
                self%newTopoList(jMol, iMol) = .true.
                self%newTopoList(iMol, jMol) = .true.
              endif
            enddo
            exit
          endif
        enddo

        if(accept) then
          return
        endif

        nActive = nActive + 1

        self%clustMemb = .false.
        startMol = 0
        leave = .false.
        do iMol = 1, self%nMolMax
          do jMol = 1, self%nMolMax
            if(self%newTopoList(jMol, iMol)) then
              startMol = iMol
              leave = .true.
              exit
            endif
          enddo
          if(leave) exit
        enddo
        if(startMol == 0) then
          if(nActive >= 2) then
            accept = .false.
            return
          endif
          nNew = 0
          nClust = nActive
        else
          self%clustMemb(startMol) = .true.
          self%newlist(1) = startMol
          nNew = 1
          nClust = 1
        endif

      !----------------------------------------------------------------
      class is(Deletion)
        self%newTopoList = self%topoList
        accept = .false.
        iMol = disp(1)%molIndx

        !First loop through and disconnect the removed molecule. 
        do jMol = 1, self%nMolMax
          if(self%newTopoList(jMol, iMol)) then
            self%newTopoList(jMol, iMol) = .false.
            self%newTopoList(iMol, jMol) = .false.
          endif
        enddo

        nActive = nActive - 1

        !Next we need to find one molecule that is still present in the system.
        !This will be the starting point for the breadth first search.
        !This is a spot for potential optimization since the algorithm is O(N^2) in the number of molecules.  
        startMol = 0
        leave = .false.
        do iMol = 1, self%nMolMax
          do jMol = iMol+1, self%nMolMax
            if(self%newTopoList(jMol, iMol)) then
              startMol = iMol
              leave = .true.
              exit
            endif
          enddo
          if(leave) then
            exit
          endif
        enddo

        if(startMol == 0) then
          if(nActive >= 2) then
            accept = .false.
            return
          endif
          self%clustMemb = .false.
          nNew = 0
          nClust = nActive
        else
          self%clustMemb = .false.
          self%clustMemb(startMol) = .true.
          self%newlist(1) = startMol
          nNew = 1
          nClust = 1
        endif

      !----------------------------------------------------------------
      class is(NewState_IsoMol)
        self%newTopoList = .false.
        if(nActive < 2) then
          accept = .true.
          return
        endif
        if(.not. allocated(disp(1)%newAtoms)) then
          accept = .false.
          return
        endif
        do iMol = 1, self%nMolMax - 1
          call trialBox%GetMolData(iMol, molType=molType, molStart=molStart)
          if(.not. trialBox%IsActive(molStart)) cycle
          iAtom = molStart + self%atomNumType(molType) - 1
          do jMol = iMol + 1, self%nMolMax
            call trialBox%GetMolData(jMol, molType=jType, molStart=jStart)
            if(.not. trialBox%IsActive(jStart)) cycle
            jAtom = jStart + self%atomNumType(jType) - 1
            rx = disp(1)%newAtoms(1, iAtom) - disp(1)%newAtoms(1, jAtom)
            ry = disp(1)%newAtoms(2, iAtom) - disp(1)%newAtoms(2, jAtom)
            rz = disp(1)%newAtoms(3, iAtom) - disp(1)%newAtoms(3, jAtom)
            call trialBox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutPairSq(molType, jType)) then
              self%newTopoList(iMol, jMol) = .true.
              self%newTopoList(jMol, iMol) = .true.
            endif
          enddo
        enddo
        startMol = 0
        leave = .false.
        do iMol = 1, self%nMolMax
          do jMol = 1, self%nMolMax
            if(self%newTopoList(jMol, iMol)) then
              startMol = iMol
              leave = .true.
              exit
            endif
          enddo
          if(leave) exit
        enddo
        if(startMol == 0) then
          if(nActive >= 2) then
            accept = .false.
            return
          endif
          self%clustMemb = .false.
          nNew = 0
          nClust = nActive
        else
          self%clustMemb = .false.
          self%clustMemb(startMol) = .true.
          self%newlist(1) = startMol
          nNew = 1
          nClust = 1
        endif

      !----------------------------------------------------------------
      class default
        error stop "All-molecule distance criteria is not compatible with this perturbation type."
      !----------------------------------------------------------------
    end select


    ! Now we grow the cluster by expanding the graph using breadth first search. The loop only needs
    ! to be done N_mol times since that's the worst case scenario of a valid cluster
    do iLimit = 1, nActive
      nNext = nNew
      nNew = 0
      do iNext = 1, nNext
        ! Get the next molecule from the queue.
        iMol = self%newlist(iNext)
        do jMol = 1, self%nMolMax
          ! If the two are connected, check next if it's already in the cluster. If not, add it to the queue.
          if(self%newTopoList(jMol, iMol)) then
            if(.not. self%clustMemb(jMol)) then
              self%clustMemb(jMol) = .true.
              nNew = nNew + 1
              nClust = nClust + 1
              self%newlist2(nNew) = jMol
            endif
          endif
        enddo
      enddo

      if(nClust >= nActive) then
        exit
      endif
      if(nNew <= 0) then
        exit
      endif

      self%newlist(1:nNew) = self%newlist2(1:nNew)
    enddo

    select type(disp)
      class is(Deletion)
        if(nClust < nActive) then
          accept = .false.
          return
        endif

        topIndex = trialBox%TypeMolFirst(disp(1)%molType) &
                 + trialBox%NMol(disp(1)%molType) - 1
        iMol = disp(1)%molIndx

        if(iMol == topIndex) then
          self%newTopoList(topIndex,:) = .false.
          self%newTopoList(:,topIndex) = .false.
        else
          do jMol = 1, self%nMolMax
            if(self%newTopoList(jMol, topIndex)) then
              if(jMol /= iMol) then
                self%newTopoList(jMol, iMol) = .true.
                self%newTopoList(iMol, jMol) = .true.
              endif
            else
              self%newTopoList(jMol, iMol) = .false.
              self%newTopoList(iMol, jMol) = .false.
            endif
          enddo
          self%newTopoList(iMol, iMol) = .false.
          self%newTopoList(topIndex,:) = .false.
          self%newTopoList(:,topIndex) = .false.
        endif

      class default
        if(nClust < nActive) then
          accept = .false.
          return
        endif
    end select
    accept = .true.

  end subroutine
!=============================================================
  subroutine DistAllMol_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand
    use Common_MolInfo, only: nMolTypes
    use ParallelVar, only: nout
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: iType, intVal, tokenPos, AllocateStat
    real(dp) :: realVal
    character(len=30) :: command

    lineStat = 0
    if(.not. allocated(self%rCut)) then
      allocate(self%rCut(1:nMolTypes), stat = AllocateStat)
      allocate(self%rCutPairSq(1:nMolTypes, 1:nMolTypes), stat = AllocateStat)
      allocate(self%atomNumType(1:nMolTypes), stat = AllocateStat)
      if(AllocateStat /= 0) then
        lineStat = -1
        return
      endif
      self%rCut = 0E0_dp
      self%rCutPairSq = 0E0_dp
      self%atomNumType = 1
    endif

    tokenPos = 2
    do iType = 1, nMolTypes
      call GetXCommand(line, command, tokenPos, lineStat)
      if(lineStat /= 0) return
      read(command, *) intVal
      self%atomNumType(iType) = intVal
      tokenPos = tokenPos + 1

      call GetXCommand(line, command, tokenPos, lineStat)
      if(lineStat /= 0) return
      read(command, *) realVal
      self%rCut(iType) = realVal
      tokenPos = tokenPos + 1
    enddo
    call DistAllMol_RebuildPairCutSq(self)

    write(nout, *) "All-molecule distance criteria:"
    do iType = 1, nMolTypes
      write(nout, *) "  type", iType, " atom =", self%atomNumType(iType), &
        " rCut =", self%rCut(iType)
    enddo
    write(nout, *) "  mixing: r_ij = (r_i + r_j) / 2"

  end subroutine
!====================================================================
  subroutine DistAllMol_Maintenance(self)
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
  end subroutine
!====================================================================
  subroutine DistAllMol_Epilogue(self)
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    logical :: accept

    call self%CheckInitialConstraint(self%parent, accept)

  end subroutine
!=============================================================
  subroutine DistAllMol_Update(self, accept)
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    logical, intent(in) :: accept

    if(.not. accept) then
      return
    endif

    self%topoList = self%newTopoList

  end subroutine
!=====================================================================
end module
!=====================================================================
