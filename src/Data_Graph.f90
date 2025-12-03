!==================================================================
module Data_Graph
   !-------------------------------------------------
  type, public :: node
    integer, allocatable :: edge(:)
    integer, allocatable :: bondid(:)

     !Mid angles are defined where this node is the center vertex
     !in the angle. IE Angle(3,1,2) is Node 1's mid angle 
     !In contrast end angles are where this vertex is on the end.
     !IE Angle(3,2,1) is Node 1's end angle 

    integer, allocatable :: midangle(:,:) 
    integer, allocatable :: midAngID(:) 
    integer, allocatable :: endangle(:,:)
    integer, allocatable :: endAngID(:) 

!    integer allocatable :: midtorsion(:,:) 
!    integer allocatable :: midtorsID(:) 
    integer, allocatable :: endtorsion(:,:)
    integer, allocatable :: endtorsID(:) 
    contains
!      procedure, pass :: Constructor => Node_Constructor
      procedure, pass :: ExpandArray => Node_ExpandArray
      procedure, pass :: Expand2DArray => Node_Expand2DArray
!      procedure, pass :: Destructor =>
      procedure, pass :: AddEdge => Node_AddEdge
      procedure, pass :: AddAngle => Node_AddAngle
      procedure, pass :: AddTorsion => Node_AddTorsion
  end type
   !-------------------------------------------------
  type, public :: graph
    integer :: nNodes = -1
    type(node), allocatable :: nodes(:)
    contains
       procedure, pass :: Constructor => Graph_Constructor
       procedure, pass :: IsConnected => Graph_IsConnected
       procedure, pass :: AddEdge => Graph_AddEdge
       procedure, pass :: AddAngle => Graph_AddAngle
       procedure, pass :: AddTorsion => Graph_AddTorsion
  end type
!==================================================================
  contains
!==================================================================
  subroutine Node_ExpandArray(self, arr, newSize)
    implicit none
    class(node), intent(inout) :: self
    integer, intent(in) :: newsize
    integer, allocatable, intent(inout) :: arr(:)
    integer, allocatable :: temparray(:)
    integer :: nElements

    if(allocated(arr)) then
      nElements = size(arr)
      allocate(temparray(1:nElements))
      temparray(1:nElements) = arr(1:nElements)
      deallocate(arr)
      allocate( arr(1:newSize) )
      arr(1:nElements) = temparray(1:nElements)
      deallocate(temparray)
    else
      allocate(arr(1:newsize))
    endif

  end subroutine
  !==================================================================
  subroutine Node_Expand2DArray(self, arr, ndim1, ndim2)
    implicit none
    class(node), intent(inout) :: self
    integer, intent(in) :: ndim1, ndim2
    integer, allocatable, intent(inout) :: arr(:, :)
    integer, allocatable :: temparray(:, :)
    integer :: nElem1, nElem2

!    write(*,*) ndim1, ndim2
    if(allocated(arr)) then
      nElem1 = size(arr, 1)
      nElem2 = size(arr, 2)
      allocate(temparray(1:nElem1, 1:nElem2))
      temparray(1:nElem1, 1:nElem2) = arr(1:nElem1, 1:nElem2)
      deallocate(arr)
      allocate( arr(1:ndim1, 1:ndim2) )
      arr(1:nElem1, 1:nElem2) = temparray(1:nElem1, 1:nElem2)
      deallocate(temparray)
    else
      allocate(arr(1:ndim1, 1:ndim2))
    endif

  end subroutine
  !==================================================================
  subroutine Node_AddEdge(self, nodeid, bondId)
    implicit none
    class(node), intent(inout) :: self
    integer, intent(in) :: nodeid
    integer, intent(in) :: bondId
    integer :: nEdge


    if(allocated(self%edge) .and. allocated(self%bondid)) then
      nEdge = size(self%edge) + 1
      call self%ExpandArray(self%edge, nEdge)
      call self%ExpandArray(self%bondid, nEdge)

    else if( allocated(self%edge) .or. allocated(self%bondid) ) then
      error stop "ERROR! The edge array and bondid arrays are not allocated together!"
    else
      nEdge = 1
      allocate(self%edge(1:1))
      allocate(self%bondid(1:1))
    endif

    self%edge(nEdge) = nodeID
    self%bondID(nEdge) = bondID

  end subroutine
  !==================================================================
  subroutine Node_AddAngle(self, angArray, idArray, nodeid1, nodeid2, angleId)
    implicit none
    class(node), intent(inout) :: self
    integer, intent(inout), allocatable :: angArray(:,:), idArray(:)
    integer, intent(in) :: nodeid1, nodeid2
    integer, intent(in) :: angleID
    integer :: nAngle
    integer, parameter :: nNeigh=2


    if( allocated(angArray) .and. allocated(idArray) ) then
      nAngle = size(angArray, 2) + 1
      call self%Expand2DArray(angArray, nNeigh, nAngle)
      call self%ExpandArray(idArray, nAngle)
    else if( allocated(angArray) .or. allocated(idArray) ) then
      error stop "ERROR! The edge array and bondid arrays are not allocated together!"
    else
      nAngle = 1
      call self%Expand2DArray(angArray, nNeigh, nAngle)
      call self%ExpandArray(idArray, nAngle)
    endif

!    write(*,*) "a", nNeigh, nAngle
!    write(*,*) "a1", size(angArray, 1), size(angArray, 2)

    angArray(1, nAngle) = nodeID1
    angArray(2, nAngle) = nodeID2
    idArray(nAngle) = angleID

  end subroutine
  !==================================================================
  subroutine Node_AddTorsion(self, torsArray, idArray, nodeid1, nodeid2, nodeid3, torsId)
    implicit none
    class(node), intent(inout) :: self
    integer, intent(inout), allocatable :: torsArray(:,:), idArray(:)
    integer, intent(in) :: nodeid1, nodeid2, nodeid3
    integer, intent(in) :: torsID
    integer :: nAngle
    integer, parameter :: nNeigh=3


    nAngle = 0
    if( allocated(torsArray) .and. allocated(idArray) ) then
      nAngle = size(torsArray, 2) + 1
      call self%Expand2DArray(torsArray, nNeigh, nAngle)
      call self%ExpandArray(idArray, nAngle)
    else if( allocated(torsArray) .or. allocated(idArray) ) then
      error stop "ERROR! The edge array and bondid arrays are not allocated together!"
    else
      nAngle = 1
      call self%Expand2DArray(torsArray, nNeigh, nAngle)
      call self%ExpandArray(idArray, nAngle)
    endif
!    write(*,*) "t", nNeigh, nAngle
!    write(*,*) "t1", size(torsArray, 1), size(torsArray, 2)
    torsArray(1, nAngle) = nodeID1
    torsArray(2, nAngle) = nodeID2
    torsArray(3, nAngle) = nodeID3
    idArray(nAngle) = torsID


  end subroutine
  !==================================================================
  ! Graph Class Functions
  !==================================================================
  subroutine Graph_Constructor(self, nNodes)
    implicit none
    class(graph), intent(inout) :: self
    integer, intent(in) :: nNodes

    if(allocated(self%nodes)) then
      error stop "Graph has already been created."
    endif

    self%nNodes = nNodes
    allocate(self%nodes(1:nNodes))



  end subroutine
  !==================================================================
  subroutine Graph_AddEdge(self, atm1, atm2, bondid)
    implicit none
    class(graph), intent(inout) :: self
    integer, intent(in) :: atm1, atm2, bondid

    call self%nodes(atm1)%AddEdge(atm2, bondid)
    call self%nodes(atm2)%AddEdge(atm1, bondid)
  end subroutine
  !==================================================================
  subroutine Graph_AddAngle(self, atm1, atm2, atm3, angID)
     !Adds
    implicit none
    class(graph), intent(inout), target :: self
    integer, intent(in) :: atm1, atm2, atm3, angID
    type(node), pointer :: curnode => null()


    curnode => self%nodes(atm2)
    call curnode%AddAngle(curnode%midAngle, curnode%midAngID, atm1, atm3, angID)

    curnode => self%nodes(atm1)
    call curnode%AddAngle(curnode%endAngle, curnode%endAngID, atm2, atm3, angID)

    curnode => self%nodes(atm3)
    call curnode%AddAngle(curnode%endAngle, curnode%endAngID, atm1, atm2, angID)

    curnode => null()
  end subroutine
  !==================================================================
  subroutine Graph_AddTorsion(self, atm1, atm2, atm3, atm4, torsID)
    implicit none
    class(graph), intent(inout), target :: self
    integer, intent(in) :: atm1, atm2, atm3, atm4, torsID
    type(node), pointer :: curnode => null()

    curnode => self%nodes(atm4)
    call curnode%AddTorsion(curnode%endTorsion, curnode%endTorsID, atm3, atm2, atm1, torsID)

    curnode => self%nodes(atm1)
    call curnode%AddTorsion(curnode%endTorsion, curnode%endTorsID, atm2, atm3, atm4, torsID)

    curnode => null()
  end subroutine
  !==================================================================
  function Graph_IsConnected(self) result(allconnected)
    implicit none
    class(graph), intent(inout), target :: self
    logical :: allconnected

    type(node), pointer :: curnode => null()
    logical :: isconnected(1:self%nNodes)
    integer :: iNode, iEdge, neiNode
    integer :: nNext, nConnected
    integer :: nextNode(1:self%nNodes)

!    write(*,*) "BLAHNL"
    if(.not. allocated(self%nodes)) then
      error stop "IsConnected has been called before the graph has been initialized."
    endif

     !If there's only one node. It's connected and thus don't waste your time here.
    if(size(self%nodes) == 1) then
      allconnected = .true.
      return
    endif

    isconnected = .false.
    allconnected = .false.

     !Start the graph walk at node 1 and expand out to see if all atoms can be found.
    nNext = 1
    nextNode = 0
    neiNode = 0
    nextNode(1) = 1
    isconnected(1) = .true.
    nConnected = 1
!    write(*,*) nNext, nextNode(1:nNext)
    do while( nNext > 0 )
!      write(*,*) nNext, nextNode(1:nNext)
      iNode = nextNode(nNext)
      
      nNext = nNext - 1
      curnode => self%nodes(iNode)
      do iEdge = 1, size(curnode%edge)
        neiNode = curnode%edge(iEdge)
        if(isconnected(neiNode)) cycle
        isconnected(neiNode) = .true.
        nConnected = nConnected + 1
        nNext = nNext + 1
        nextNode(nNext) = neiNode
      enddo
      curnode => null()
      if(nConnected == size(self%nodes)) exit
    enddo
    
    if(nConnected == size(self%nodes)) then
      allconnected = .true.
    endif

  end function
  !==================================================================
  subroutine Graph_CheckIntegrity(self)
    implicit none
    class(graph), intent(inout), target :: self
    type(node), pointer :: curnode => null()

!    do iNode = 1, size(self%nodes)
!      curnode => self%nodes(iNode)
!    enddo
  end subroutine
!==================================================================
end module
!==================================================================
