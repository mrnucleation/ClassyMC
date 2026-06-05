!==================================================================
module Graph_Module
  implicit none
  private

  public :: undirected_graph

  !------------------------------------------------------------------------------
  type :: undirected_graph

    integer, private :: nNodes = 0  ! Number of nodes currently in the graph
    integer, private :: nEdges = 0 
    integer, private :: nMaxNodes = 0  !Maximum number of nodes allowed

    integer, allocatable :: neighbor(:,:)      ! Adjacency list: neighbor(i, :) holds neighbors of node i
    integer, allocatable :: currentNeighbor(:)  ! Tracks the current number of neighbors for each node
    integer, allocatable :: nodeIDs(:)
  contains
    procedure :: init 
    procedure :: reset_graph
    procedure :: reset_single_nodes_edges
    procedure :: add_node
    procedure :: remove_node
    procedure :: add_edge
    procedure :: remove_edge
    procedure :: has_edge
    procedure :: get_neighbors
    procedure :: sort_neighlists
    procedure :: reindex_nodes
    procedure :: move_node
    procedure :: is_connected
    procedure :: count_neighbors
    procedure :: get_nNodes
    procedure :: find_node_by_ID
    procedure :: copy_graph
    procedure :: print_graph
    procedure :: finalize => finalize_graph
    generic, public :: assignment(=) => copy_graph

  end type undirected_graph
  !------------------------------------------------------------------------------

contains

   ! 05/23/2026 - Where I left off
   ! 1. Need to  write the routine which disconnects a node completely from all its edges
   ! 2. Need to track nodeIDs separate from the indicies.
   ! 3. 
  !------------------------------------------------------------------------------
  ! Initialize graph with maximum number of nodes
  subroutine init(self, max_nodes)
    implicit none
    class(undirected_graph), intent(inout) :: self
    integer, intent(in) :: max_nodes
    integer :: STATUS

    self%nMaxNodes = max_nodes
    !write(*,*) "Initializing graph with max nodes = ", max_nodes

    self%nNodes = 0
    allocate(self%neighbor(max_nodes, max_nodes))   ! Over-allocate columns for neighbors
    self%neighbor = 0
    allocate(self%currentNeighbor(max_nodes))
    allocate(self%nodeIDs(max_nodes))
    self%currentNeighbor = 0
  end subroutine init

  !------------------------------------------------------------------------------
  subroutine reset_graph(self)
    implicit none
    class(undirected_graph), intent(inout) :: self
    self%nNodes = 0
    self%nEdges = 0
    self%neighbor = 0
    self%currentNeighbor = 0
    self%nodeIDs = 0
  end subroutine reset_graph
  !------------------------------------------------------------------------------
  subroutine reset_single_nodes_edges(self, node1)
     !Removes all edges from a single node in the graph without touching the rest of the graph structure.  This is useful for algorithms that need to temporarily disconnect a node from the graph without affecting the overall structure or other nodes.
    use SearchSort, only: BinarySearch
    implicit none
    class(undirected_graph), intent(inout) :: self
    integer, intent(in) :: node1
    integer :: iNode, jNode


    !This requires us to loop over the neighbors of node1 and remove the edge between node1 and each of its neighbors.  
    !Doing this from the end of the neighbor list to the beginning allows us to remove edges without needing to worry about shifting the neighbor list down as we remove edges since we won't be accessing the already removed neighbors again.

    do while (self%currentNeighbor(node1) > 0)
      iNode = self%neighbor(node1, self%currentNeighbor(node1))
      call self%remove_edge(node1, iNode)
    end do



  end subroutine reset_single_nodes_edges
  !------------------------------------------------------------------------------
  subroutine copy_graph(dest, src)
     !Clones the graph structure and data from src to dest.  Note that dest must already be initialized with the same max_nodes as src before calling this function.
    implicit none
    class(undirected_graph), intent(inout) :: dest
    class(undirected_graph), intent(in) :: src



    !write(*,*) "------- Copying graph with ", src%nNodes, " nodes and ", src%nEdges, " edges. --------"
    !write(*,*) "Source graph max nodes: ", src%nMaxNodes
    !write(*,*) "Destination graph max nodes: ", dest%nMaxNodes
    !Check that the destination graph has the same max_nodes as the source graph.  If not, error out to prevent memory issues.
    if (dest%nMaxNodes /= src%nMaxNodes) then
      error stop "Error: Destination graph must have the same max_nodes as the source graph for copy_graph to work."
    end if


    dest%nNodes = src%nNodes
    dest%nEdges = src%nEdges
    dest%neighbor = src%neighbor
    dest%currentNeighbor = src%currentNeighbor

    dest%nodeIDs = src%nodeIDs

  end subroutine copy_graph
  !------------------------------------------------------------------------------
  integer function find_node_by_ID(self, nodeID) result(nodeIndex)
    use SearchSort, only: BinarySearch
    implicit none
    class(undirected_graph), intent(in) :: self
    integer, intent(in) :: nodeID
    integer :: iNode

    !write(*,*) "Finding node with ID ", nodeID, " in graph with ", self%nNodes, " nodes."
    nodeIndex = BinarySearch(nodeID, self%nodeIDs(1:self%nNodes))

    !write(*,*) "Current ID list: ", self%nodeIDs(1:self%nNodes)
    if(nodeIndex == 0) then
      error stop "Error: Node ID not found in the graph."
    end if

  end function find_node_by_ID
  !------------------------------------------------------------------------------
  ! Add a new node and insert it in sorted nodeID order so BinarySearch remains valid.
  subroutine add_node(self, nodeID)
    use SearchSort, only: SimpleSearch
    implicit none
    class(undirected_graph), intent(inout) :: self
    integer, intent(in), optional :: nodeID
    integer :: newID, insertPos, iNode, jNode

    self%nNodes = self%nNodes + 1

    ! Determine the ID for the new node.
    if (present(nodeID)) then
      newID = nodeID
    else
      newID = self%nNodes
    end if

    ! Duplicate check via linear search (works for any ordering).
    if (SimpleSearch(newID, self%nodeIDs(1:self%nNodes-1)) > 0) then
      error stop "Error: Node ID already exists in the graph."
    end if

    ! Find the sorted insertion position (first existing ID > newID).
    insertPos = self%nNodes   ! default: append at the end
    do iNode = 1, self%nNodes - 1
      if (self%nodeIDs(iNode) > newID) then
        insertPos = iNode
        exit
      end if
    end do

    ! Shift existing nodes at [insertPos..nNodes-1] right by 1 to [insertPos+1..nNodes].
    do iNode = self%nNodes, insertPos + 1, -1
      self%nodeIDs(iNode)         = self%nodeIDs(iNode - 1)
      self%currentNeighbor(iNode) = self%currentNeighbor(iNode - 1)
      self%neighbor(iNode, :)     = self%neighbor(iNode - 1, :)
    end do

    ! Place the new (empty) node at insertPos.
    self%nodeIDs(insertPos)         = newID
    self%currentNeighbor(insertPos) = 0
    self%neighbor(insertPos, :)     = 0

    ! Update neighbor position references in all nodes: positions >= insertPos are now +1.
    if (insertPos < self%nNodes) then
      do iNode = 1, self%nNodes
        if (iNode == insertPos) cycle   ! new node has no edges yet
        do jNode = 1, self%currentNeighbor(iNode)
          if (self%neighbor(iNode, jNode) >= insertPos) then
            self%neighbor(iNode, jNode) = self%neighbor(iNode, jNode) + 1
          end if
        end do
      end do
    end if

  end subroutine add_node
  !------------------------------------------------------------------------------
  ! Remove a node.  This is more involved since it requires removing all edges connected to this node and then shifting all subsequent nodes down by one to fill the gap.  
  subroutine remove_node(self, node1, nodeId)
      use SearchSort, only: BinarySearch
      implicit none
      class(undirected_graph), intent(inout) :: self
      integer, intent(in), optional :: node1, nodeId
      integer :: iNode, jNode
      integer :: nodeToRemove
      !Disconnect the node from all its neighbors by removing the edges between them.  This will also shift all subsequent nodes down by one to fill the gap left by the removed node.  Note that this function assumes that the neighbor lists are already sorted so that it can use binary search to check for neighbors.

      if (present(nodeId)) then
        nodeToRemove = BinarySearch(nodeId, self%nodeIDs(1:self%nNodes))
      elseif (present(node1)) then
        nodeToRemove = node1
      else
        error stop "Error: Either node1 or nodeId must be provided to remove_node."
      end if

      do iNode = 1, self%nNodes
        if (iNode == nodeToRemove) cycle
        call self%remove_edge(nodeToRemove, iNode)
      end do

      ! Clean up the neighbor list for the removed node
      self%neighbor(nodeToRemove, :) = 0
      self%currentNeighbor(nodeToRemove) = 0

      !Shift all subsequent nodes down by one to fill the gap left by the removed node.
      do iNode = nodeToRemove + 1, self%nNodes
        self%neighbor(iNode - 1, :) = self%neighbor(iNode, :)
        self%currentNeighbor(iNode - 1) = self%currentNeighbor(iNode)
        self%nodeIDs(iNode - 1) = self%nodeIDs(iNode)
      end do
      self%nNodes = self%nNodes - 1

      !Reindex the neighbor lists to account for the removed node.
      do iNode = 1, self%nNodes
        do jNode = 1, self%currentNeighbor(iNode)
          if (self%neighbor(iNode, jNode) > nodeToRemove) then
            self%neighbor(iNode, jNode) = self%neighbor(iNode, jNode) - 1
          end if
        end do
      end do
  end subroutine remove_node
  !------------------------------------------------------------------------------
  subroutine sort_neighlists(self, bubble)
     ! Sort the neighborlist so binary search can be used in subsequent calls. 
     ! By default it uses quicksort, but if bubble is set to true it will use bubble sort instead. 
     ! Bubble sort is much faster for small lists or if they are nearly sorted which these lists usually are.
    use SearchSort, only: BubbleSort, QSort
    implicit none
    class(undirected_graph), intent(inout) :: self
    logical, intent(in), optional :: bubble
    integer :: iNode


    do iNode = 1, self%nNodes
      if (present(bubble) .and. bubble) then
        call BubbleSort(self%neighbor(iNode, 1:self%currentNeighbor(iNode)))
      else
        call QSort(self%neighbor(iNode, 1:self%currentNeighbor(iNode)))
      endif
    enddo


  end subroutine
  !------------------------------------------------------------------------------
  ! Add undirected edge between node1ID and node2ID (1-based)
  subroutine add_edge(self, node1ID, node2ID)

    class(undirected_graph), intent(inout) :: self
    integer, intent(in) :: node1ID, node2ID

    integer :: iNode

    if (node1ID < 1 .or. node1ID > self%nNodes &
        .or. node2ID < 1 .or. node2ID > self%nNodes) then
      print *, "Error: Invalid node indices"
      return
    end if

    if (node1ID == node2ID) return  ! No self-loops by default

    if (self%has_edge(node1ID, node2ID)) return  ! No duplicate edges

    ! Add node2ID to node1ID's neighbor list
    self%currentNeighbor(node1ID) = self%currentNeighbor(node1ID) + 1
    self%neighbor(node1ID, self%currentNeighbor(node1ID)) = node2ID

    ! Add node1ID to node2ID's neighbor list (undirected)
    self%currentNeighbor(node2ID) = self%currentNeighbor(node2ID) + 1
    self%neighbor(node2ID, self%currentNeighbor(node2ID)) = node1ID

    !write(*,*) "Added edge between node ", node1ID, " and node ", node2ID

    call self%sort_neighlists(.true.)
  end subroutine add_edge
  !------------------------------------------------------------------------------
  subroutine remove_edge(self, node1ID, node2ID)
    use SearchSort, only: BinarySearch
    class(undirected_graph), intent(inout) :: self
    integer, intent(in) :: node1ID, node2ID

    integer :: iNode, jNode

    if (node1ID < 1 .or. node1ID > self%nNodes &
        .or. node2ID < 1 .or. node2ID > self%nNodes) then
      print *, "Error: Invalid node indices"
      return
    end if

    ! Remove node2ID from node1ID's neighbor list
    ! Shift remaining neighbors left to fill the gap (IE 1 2 3 4 5 -> remove 3 -> 1 2 _ 4 5  -> 1 2 4 5 _)
    iNode = BinarySearch(node2ID,  self%neighbor(node1ID, 1:self%currentNeighbor(node1ID)))
    if (iNode > 0) then
      self%neighbor(node1ID, iNode:self%currentNeighbor(node1ID)-1) = self%neighbor(node1ID, iNode+1:self%currentNeighbor(node1ID))
      self%currentNeighbor(node1ID) = self%currentNeighbor(node1ID) - 1
    end if


    ! Remove node1ID from node2ID's neighbor list (undirected)
    jNode = BinarySearch(node1ID,  self%neighbor(node2ID, 1:self%currentNeighbor(node2ID)))
    if (jNode > 0) then
      self%neighbor(node2ID, jNode:self%currentNeighbor(node2ID)-1) = self%neighbor(node2ID, jNode+1:self%currentNeighbor(node2ID))
      self%currentNeighbor(node2ID) = self%currentNeighbor(node2ID) - 1
    end if

    call self%sort_neighlists(.true.)

  end subroutine remove_edge
  !------------------------------------------------------------------------------
  logical function has_edge(self, node1, node2) result(res)
    use SearchSort, only: BinarySearch
    class(undirected_graph), intent(in) :: self
    integer, intent(in) :: node1, node2
     ! Check if an edge exists between node1 and node2 using binary search
    integer :: i
    res = .false.
    if (node1 < 1 .or. node1 > self%nNodes &
        .or. node2 < 1 .or. node2 > self%nNodes) then
      print *, "Error: Invalid node indices"
      return
    end if

    i = BinarySearch( node2, self%neighbor(node1, 1:self%currentNeighbor(node1)) )
    if (i > 0) res = .true.

  end function has_edge
  !------------------------------------------------------------------------------
  logical function is_connected(self, startnode, target_list) result(reslt)
      ! Check if the graph is fully connected (i.e. there's a path between any two nodes)
      ! This is done by performing a breadth first search starting from the first node and seeing if we can reach all other nodes.  
      ! Note that this function assumes that the neighbor lists are already sorted so that it can use binary search to check for neighbors.
    class(undirected_graph), intent(in) :: self
    integer :: iNode, jNode
    integer, intent(in), optional :: startnode ! If not provided, the search will start from node 1 by default.
    integer, intent(in), optional :: target_list(:) ! If provided, end the search early if all nodes in this list have been visited.  

    integer :: nVisited ! How many nodes have been visited so far during the breadth first search
    integer :: nStack ! How many nodes are currently in the stack for breadth first search
    integer :: nextNode(1:self%nMaxNodes) ! Stack for breadth first search, determines the next nodes to visit
    logical :: visited(1:self%nMaxNodes)  ! Tracks which nodes have been visited already during the search


    reslt = .false.

    ! Using a stack, perform a breadth first search to see if all nodes can be reached from the first node.
    visited = .false.
    nVisited = 0
    nStack = 1

    if (present(startnode)) then
      !write(*,*) "Setting initial node: ", startnode, "with NNodes = ", self%nNodes
      nextNode(1) = startnode
      visited(startnode) = .true.
    else
      nextNode(1) = 1
      visited(1) = .true.
    end if

    !write(*,*) "Number of nodes in graph: ", self%nNodes
    !write(*,*) "Starting breadth first search for connectivity check from node ", nextNode(1)

    do while (nStack > 0)
      !write(*,*) nStack, " nodes in stack, ", nVisited, " nodes visited so far."
       ! Take the next node from the stack and mark it as visited.  Then add all of its unvisited neighbors to the stack.
      iNode = nextNode(nStack)
      !write(*,*) "Visiting node ", iNode
      nStack = nStack - 1
      nVisited = nVisited + 1

      !write(*,*) nVisited

        ! Add all unvisited neighbors of the current node to the stack.
      do jNode = 1, self%currentNeighbor(iNode)
        !write(*,*) "Checking neighbor ", jNode 
        if (.not. visited(self%neighbor(iNode, jNode))) then
          visited(self%neighbor(iNode, jNode)) = .true.
          nStack = nStack + 1
          nextNode(nStack) = self%neighbor(iNode, jNode)
          !write(*,*) "Adding node ", self%neighbor(iNode, jNode), " to stack"
        end if
      end do ! End of loop over neighbors of the current node

      if (present(target_list)) then
         ! If a target list is provided, check if all nodes in the target list have been visited.  If so, we can end the search early and return true.
        if (all( visited(target_list) )) then
          reslt = .true.
          !write(*,*) "All target nodes in list have been visited. Ending search early with result = true."
          return
        end if
      end if

        !Early break. If we've already visited all the nodes, we can stop and return true without needing to continue the search.
      if (nVisited == self%nNodes) then
        reslt = .true.
        !write(*,*) "All nodes have been visited. Ending search with result = true."
        return
      end if

    end do ! End of breadth first search loop

    if (nVisited == self%nNodes) then
      reslt = .true.
    else 
      reslt = .false.
    end if

  end function is_connected
  !------------------------------------------------------------------------------
  subroutine reindex_nodes(self, bottomrange, toprange)
    implicit none
    class(undirected_graph), intent(inout) :: self
    integer, intent(in) :: bottomrange, toprange
    integer :: iNode

    do iNode = 1, self%nNodes
      if (self%nodeIDs(iNode) >= bottomrange .and. self%nodeIDs(iNode) <= toprange) then
        self%nodeIDs(iNode) = self%nodeIDs(iNode) - 1
      end if
    end do    

  end subroutine reindex_nodes
  !------------------------------------------------------------------------------
  subroutine move_node(self, fromPos, toPos, newID)
    ! Move a node from graph position fromPos to toPos, renaming its molecule ID to newID.
    ! Used to mirror ClassyMC's swap-with-last deletion: when molecule lastMol is copied
    ! into the deleted molecule's slot, the graph must similarly relocate lastMol's node.
    ! Assumes toPos <= fromPos (moving toward a lower position index).
    implicit none
    class(undirected_graph), intent(inout) :: self
    integer, intent(in) :: fromPos, toPos, newID
    integer :: iNode, jNode
    integer :: tempNeighbor(1:self%nMaxNodes)
    integer :: tempCount

    if (fromPos == toPos) then
      self%nodeIDs(toPos) = newID
      return
    end if

    ! Save the data from fromPos (the node being moved).
    tempCount = self%currentNeighbor(fromPos)
    tempNeighbor(:) = 0
    tempNeighbor(1:tempCount) = self%neighbor(fromPos, 1:tempCount)

    ! Pre-adjust the saved neighbor refs to account for the right-shift of [toPos..fromPos-1].
    do jNode = 1, tempCount
      if (tempNeighbor(jNode) >= toPos .and. tempNeighbor(jNode) < fromPos) then
        tempNeighbor(jNode) = tempNeighbor(jNode) + 1
      end if
    end do

    ! Shift nodes at [toPos..fromPos-1] right by 1 to [toPos+1..fromPos].
    do iNode = fromPos, toPos + 1, -1
      self%nodeIDs(iNode)          = self%nodeIDs(iNode - 1)
      self%currentNeighbor(iNode)  = self%currentNeighbor(iNode - 1)
      self%neighbor(iNode, :)      = self%neighbor(iNode - 1, :)
    end do

    ! Place the saved node at toPos with its new ID.
    self%nodeIDs(toPos)            = newID
    self%currentNeighbor(toPos)    = tempCount
    self%neighbor(toPos, :)        = 0
    self%neighbor(toPos, 1:tempCount) = tempNeighbor(1:tempCount)

    ! Update neighbor position references in all OTHER nodes.
    ! fromPos moved → toPos; positions [toPos..fromPos-1] shifted right by 1.
    do iNode = 1, self%nNodes
      if (iNode == toPos) cycle   ! Already handled via tempNeighbor above.
      do jNode = 1, self%currentNeighbor(iNode)
        if (self%neighbor(iNode, jNode) == fromPos) then
          self%neighbor(iNode, jNode) = toPos
        else if (self%neighbor(iNode, jNode) >= toPos .and. &
                 self%neighbor(iNode, jNode) < fromPos) then
          self%neighbor(iNode, jNode) = self%neighbor(iNode, jNode) + 1
        end if
      end do
    end do

    ! Re-sort neighbor lists since position updates may have disturbed order.
    call self%sort_neighlists(bubble=.true.)

  end subroutine move_node
  !------------------------------------------------------------------------------
  function get_neighbors(self, node1) result(neigh)
     ! Gets the neighbors of node node1.  Note that the size of the returned array is equal to the number of neighbors, not the maximum number of neighbors.
    class(undirected_graph), intent(in) :: self
    integer, intent(in) :: node1
    integer, allocatable :: neigh(:)
    integer :: i, count

    count = count_neighbors(self, node1)
    allocate(neigh(count))
    count = 0
    do i = 1, self%currentNeighbor(node1)
      if (self%neighbor(node1, i) > 0) then
        count = count + 1
        neigh(count) = self%neighbor(node1, i)
      end if
    end do

  end function get_neighbors
  !------------------------------------------------------------------------------
  integer function count_neighbors(self, node1) result(cnt)
    class(undirected_graph), intent(in) :: self
    integer, intent(in) :: node1
    integer :: i

    cnt = self%currentNeighbor(node1)

  end function count_neighbors
  !------------------------------------------------------------------------------
  integer function get_nNodes(self) result(n)
    class(undirected_graph), intent(in) :: self
    n = self%nNodes
  end function get_nNodes
  !------------------------------------------------------------------------------
  subroutine print_graph(self)
    class(undirected_graph), intent(in) :: self
    integer :: i, j

    print *, "Undirected Graph with", self%nNodes, "nodes:"
    do i = 1, self%nNodes
      write(*, '(A,I0,A)', advance='no') "  ", i, " : "
      do j = 1, self%currentNeighbor(i)
        if (self%neighbor(i, j) > 0) write(*, '(I0," ")', advance='no') self%neighbor(i, j)
      end do
      print *
    end do
  end subroutine print_graph

  !------------------------------------------------------------------------------
  subroutine finalize_graph(self)
    implicit none
    class(undirected_graph), intent(inout) :: self
    if (allocated(self%neighbor)) deallocate(self%neighbor)
    if (allocated(self%currentNeighbor)) deallocate(self%currentNeighbor)
  end subroutine finalize_graph
  !------------------------------------------------------------------------------
!==================================================================
end module Graph_Module
!==================================================================