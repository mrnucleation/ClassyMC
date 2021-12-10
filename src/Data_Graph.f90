!==================================================================
module Graph
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
      procedure, pass :: Constructor => Node_Constructor
      procedure, pass :: ExpandArray => Node_ExpandArray
      procedure, pass :: Expand2DArray => Node_Expand2DArray
!      procedure, pass :: Destructor =>
      procedure, pass :: AddEdge => Node_AddEdge
  end type
   !-------------------------------------------------
  type, public :: graph
    type(node), allocatable :: nodes(:)
    contains
       procedure, pass :: AddEdge => Graph_AddEdge
       procedure, pass :: AddAngles => Graph_AddAngles
       procedure, pass :: AddTorsion => Graph_AddTorsion
  end type
!==================================================================
  contains
!==================================================================
subroutine Node_ExpandArray(self, arr, newSize)
  implicit none
  class(node), intent(inout) :: self
  integer, intent(in) :: newsize
  integer, allocatable :: temparray(:)
  integer :: nElements

  if(allocated(arr)) then
    nElements = size(arr)
    allocate(temparray(1:nElements))
    temparray(1:nElements) = arr(1:nElements)
    deallocate(arr)
    allocate( arr(1:newSize) )
    arr(1:nElements) = temparray(1:nElements)
    nEdges = nEdges + 1
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
  integer, allocatable :: temparray(:, :)
  integer :: nElem1, nElem2

  if(allocated(arr)) then
    nElem1 = size(arr, 1)
    nElem2 = size(arr, 2)
    allocate(temparray(1:nElem1, 1:nElem2))
    temparray(1:nElem1, 1:nElem2) = arr(1:nElem1, 1:nElem2)
    deallocate(arr)
    allocate( arr(1:ndim1, 1:ndim2) )
    arr(1:nElem1, 1:nElem2) = temparray(1:nElem1, 1:nElem2)
    nEdges = nEdges + 1
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
  integer, intent(in) :: inBondId
  integer :: nEdges


  if(allocated(self%edge) .and. allocated(self%bondid)) then
    nEdges = size(self%edge) + 1
    call self%ExpandAray(self%edge, nEdges)
    call self%ExpandAray(self%bondid, nEdges)

  else if( allocated(self%edge) .or. allocated(self%bondid) ) then
    error stop "ERROR! The edge array and bondid arrays are not allocated together!"
  else
    nEdges = 1
    allocate(self%edge(1:1))
    allocate(self%bondid(1:1))
  endif

  self%edge(nEdges) = nodeID
  self%bondID(nEdges) = bondID

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


  if( allocated(angArray) .neqv. allocated(idArray) ) then
    error stop "ERROR! The edge array and bondid arrays are not allocated together!"
  else
    if(allocated(angArray)) then
      nAngle = size(angArray, 2)
    else
      nAngle = 1
    endif
    call self%Expand2DArray(angArray, nNeigh, nAngle)
    call self%ExpandArray(idArray, nAngle)
  endif

  angArray(1, nAngles) = nodeID1
  angArray(2, nAngles) = nodeID2
  idArray(nAngles) = angleID

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


  if( allocated(torsArray) .neqv. allocated(idArray) ) then
    error stop "ERROR! The edge array and bondid arrays are not allocated together!"
  else
    if(allocated(torsArray)) then
      nAngle = size(torsArray, 2)
    else
      nAngle = 1
    endif
    call self%Expand2DArray(torsArray, nNeigh, nAngle)
    call self%ExpandArray(idArray, nAngle)
  endif

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

  allocate(self%nodes(1:nNodes))



end subroutine
!==================================================================
subroutine Graph_AddEdge(self, atm1, atm2)
  implicit none
  class(graph), intent(inout) :: self
  integer, intent(in) :: atm1, atm2

  self%node(atm1)%AddEdge(atm2)
  self%node(atm2)%AddEdge(atm1)
end subroutine
!==================================================================
subroutine Graph_AddAngle(self, atm1, atm2, atm3, angID)
  implicit none
  class(graph), intent(inout), target :: self
  integer, intent(in) :: atm1, atm2, atm3, angID
  type(node), pointer :: node => null()


  node => self%node(atm2)
  call node%AddAngle(node%midAngle, node%midAngID, atm1, atm3, angID)

  node => self%node(atm1)
  call node%AddAngle(node%endAngle, node%endAngID, atm2, atm3, angID)

  node => self%node(atm3)
  call node%AddAngle(node%endAngle, node%endAngID, atm1, atm2, angID)
end subroutine
!==================================================================
subroutine Graph_AddTorsion(self, atm1, atm2, atm3, atm4, torsID)
  implicit none
  class(graph), intent(inout), target :: self
  integer, intent(in) :: atm1, atm2, atm3, torsID
  type(node), pointer :: node => null()

  node => self%node(atm4)
  call node%AddTorsion(node%endTorsion, node%midTorsID, atm3, atm2, atm1, torsID)

  node => self%node(atm1)
  call node%AddTorsion(node%midAngle, node%endTorsID, atm2, atm3, atm4, torsID)

end subroutine
!==================================================================
end module
!==================================================================