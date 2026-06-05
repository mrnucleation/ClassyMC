!====================================================================
! Stillinger distance criteria over all molecule types in the box.
! Input (per type i): reference atom index, then rCut(i).
! Pair cutoff uses arithmetic mixing: r_ij = (r_i + r_j) / 2
!====================================================================
module Constrain_DistanceCriteriaAllMol
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use CoordinateTypes, only: Displacement, Deletion, Addition
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout
  use Graph_Module, only: undirected_graph
  implicit none

  type, public, extends(constraint) :: DistCriteriaAllMol
    integer :: boxID = 1
    integer :: nMolMax = 0

    integer, allocatable :: atomNumType(:)
    real(dp), allocatable :: rCut(:)
    real(dp), allocatable :: rCutPairSq(:,:)

    type(undirected_graph) :: topoGraph, temp_topoGraph

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
    integer :: iType

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

    call self%topoGraph%init(self%nMolMax)
    call self%temp_topoGraph%init(self%nMolMax)


    IF (AllocateStat /= 0) STOP "Allocation Error in All-Mol Distance Constraint Constructor"
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
    integer :: iType, jType, molStart, jStart, nNext, nNext2, iNext
    integer :: iGraphID, jGraphID
    real(dp) :: rx, ry, rz, rsq

    accept = .true.

    call self%topoGraph%reset_graph()




    ! First we need to create the nodes for each of the active molecules in the system.  We can do this with a single loop through all the molecules and checking if they're active.
    do iMol = 1, self%nMolMax
      call trialBox%GetMolData(iMol, molType=iType, molStart=molStart)
      !write(*,*) "Checking molecule ", iMol, " of type ", iType, " starting at atom ", molStart
      if(.not. trialBox%IsActive(molStart)) cycle
      call self%topoGraph%Add_Node(nodeID=iMol)
    enddo

    ! First check if we even have more than one molecule in the first place
    nActive = trialBox%nMolTotal
    if(nActive < 2) then
      return
    endif


    ! Loop through all the molecules current in the system and build the initial topology graph based on the distance criteria.  
    ! We only need to check each pair once, so we can use a double loop with jMol starting from iMol + 1.
    do iMol = 1, self%nMolMax - 1
      call trialBox%GetMolData(iMol, molType=iType, molStart=molStart)
      if(.not. trialBox%IsActive(molStart)) cycle
      iAtom = molStart + self%atomNumType(iType) - 1

      do jMol = iMol + 1, self%nMolMax
        call trialBox%GetMolData(jMol, molType=jType, molStart=jStart)
        if(.not. trialBox%IsActive(jStart)) cycle
        jAtom = jStart + self%atomNumType(jType) - 1
        rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
        ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
        rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
        call trialBox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutPairSq(iType, jType)) then
          ! Add_Edge works on graph node positions, not molecule indices.  In a
          ! multi-type system the molecule indices are sparse, so map each one to
          ! its current graph position first.
          iGraphID = self%topoGraph%find_node_by_ID(iMol)
          jGraphID = self%topoGraph%find_node_by_ID(jMol)
          call self%topoGraph%Add_Edge(iGraphID, jGraphID)
        endif
      enddo
    enddo

    if(.not. self%topoGraph%is_connected()) then
      accept = .false.
      write(nout,*) "All-Mol cluster criteria check failed!"
      return
    endif

    write(nout,*) "All-Mol cluster criteria check succeeded!"
    accept = .true.

  end subroutine
!=============================================================
  subroutine DistAllMol_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteriaAllMol), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    integer :: iDisp
    integer :: molIndx, molType, molStart, iAtom
    integer :: iMol, jMol, jType, jStart, jAtom
    integer :: nodeIndex
    integer :: graphID, jGraphID, lastGraphID
    integer :: topIndex

    real(dp) :: rx, ry, rz, rsq

    accept = .true.

    !call self%temp_topoGraph%print_graph()
    !call self%topoGraph%print_graph()
    self%temp_topoGraph = self%topoGraph
    select type(disp)
        !.........
      class is(Displacement)
        !write(*,*) "Displacement"
         ! Check if we have more than one molecule in the system.  If not, then we can just accept the move since it won't change any of the relevant distances.
        if(trialBox%nMolTotal < 2) then
          accept = .true.
          return
        endif
         ! The variable is set initially to accept=.true. because if the molecule being displaced is not actually moving its reference atom, then we can just skip the distance checks 
         ! and accept the move since it won't change any of the relevant distances.
        accept = .true.
        !write(*,*) disp(1)%molIndx
        graphID = self%temp_topoGraph%find_node_by_ID(disp(1)%molIndx)
        do iDisp = 1, size(disp)
          molIndx = disp(iDisp)%molIndx
          call trialBox%GetMolData(molIndx, molType=molType, molStart=molStart)
          iAtom = molStart + self%atomNumType(molType) - 1

           ! If the reference atom for this molecule is being displaced, then we need to check all the distances for this molecule in its new position to see if it creates or breaks any edges in the temporary topology graph.
          if(disp(iDisp)%atmIndx == iAtom) then
            accept = .false.
            self%temp_topoGraph = self%topoGraph

              !Disconnect the old edges for this molecule in the temporary graph since we're going to check all the distances again for this molecule in its new position.  
            call self%temp_topoGraph%reset_single_nodes_edges(node1=graphID)

              !Now use the distance to create the new edges for this molecule in the temporary graph.  
            do jMol = 1, self%nMolMax 
              if(molIndx == jMol) cycle
              call trialBox%GetMolData(jMol, molType=jType, molStart=jStart)
              if(.not. trialBox%IsActive(jStart)) cycle
              jAtom = jStart + self%atomNumType(jType) - 1
              rx = disp(iDisp)%x_new - trialBox%atoms(1, jAtom)
              ry = disp(iDisp)%y_new - trialBox%atoms(2, jAtom)
              rz = disp(iDisp)%z_new - trialBox%atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
              !write(nout,*) "Distance squared between molecule ", molIndx, " and molecule ", jMol, " is ", rsq
              if(rsq < self%rCutPairSq(molType, jType)) then
                jGraphID = self%temp_topoGraph%find_node_by_ID(jMol)
                call self%temp_topoGraph%Add_Edge(graphID, jGraphID)
              endif
            enddo
          endif
        enddo

        if(accept) then
          return
        endif
        
        !.........
      class is(Addition) ! If there's a new molecule being added. 
        !write(*,*) "Add"

        call self%temp_topoGraph%Add_Node(nodeID=disp(1)%molIndx)
        graphID = self%temp_topoGraph%find_node_by_ID(disp(1)%molIndx)

         ! For an addition move, we only need to check the distances for the molecule being added since no existing edges can be broken by adding a new molecule.  
         ! If any of the distances from the new molecule to existing molecules violates the distance criteria, then we reject the move.  Otherwise we accept it.  
         ! We can exit as soon as we find one violation, so we set

        accept = .true.
        do iDisp = 1, size(disp)

          molIndx = disp(iDisp)%molIndx
          !write(nout,*) "Addition: molIndx =", molIndx
          call trialBox%GetMolData(molIndx, molType=molType, molStart=molStart)
          iAtom = molStart + self%atomNumType(molType) - 1

          if(disp(iDisp)%atmIndx == iAtom) then

            !write(nout,*) "Addition: iAtom =", iAtom
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
                jGraphID = self%temp_topoGraph%find_node_by_ID(jMol)
                call self%temp_topoGraph%Add_Edge(graphID, jGraphID)
              endif

            enddo
          endif
        enddo

        if(accept) then
          return
        endif

        !.........
      class is(Deletion)
        !write(*,*) "Deletion"
        ! ClassyMC's DeleteMol uses a swap-with-last scheme: the last molecule of this
        ! type (topIndex) is copied into the deleted slot (disp(1)%molIndx).  We must
        ! mirror this in the graph: remove the deleted node, then rename the lastMol
        ! node to the deleted slot's index and move it to the correct sorted position.
        !
        ! Use the same formula as DeleteMol itself: TypeMolFirst + NMol - 1.
        ! (TypeMolLast is a static capacity array, not the current active last index.)
        topIndex = trialBox%TypeMolFirst(disp(1)%molType) + trialBox%NMol(disp(1)%molType) - 1
        graphID  = self%temp_topoGraph%find_node_by_ID(disp(1)%molIndx)
        call self%temp_topoGraph%Remove_Node(node1=graphID)

        if (topIndex /= disp(1)%molIndx) then
          ! lastMol node still present; find its current position (shifted after Remove_Node)
          ! and move it to graphID with the deleted molecule's index as its new ID.
          lastGraphID = self%temp_topoGraph%find_node_by_ID(topIndex)
          call self%temp_topoGraph%move_node(fromPos=lastGraphID, toPos=graphID, &
                                             newID=disp(1)%molIndx)
        end if

        graphID = 1  ! Start connectivity check from node 1 (any remaining node is fine).
          
        !.........
      class default
        error stop "All-molecule distance criteria is not compatible with this perturbation type. Currently only displacement, addition, and deletion moves are allowed."
    end select


    accept = self%temp_topoGraph%is_connected(startNode=graphID)

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

    ! 

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

    if (.not. accept) then
      return
    endif

    self%topoGraph = self%temp_topoGraph

  end subroutine
!=====================================================================
end module
!=====================================================================
