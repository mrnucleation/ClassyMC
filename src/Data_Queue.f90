!==================================================================
! Queue data structure for general use in algorithms
! Provides FIFO (First-In-First-Out) queue operations
!==================================================================
module Data_Queue
  use VarPrecision
  implicit none
  private

  !-------------------------------------------------
  ! Integer Queue - useful for graph traversal, BFS, etc.
  !-------------------------------------------------
  type, public :: IntQueue
    private
    integer :: capacity = 0
    integer :: head = 1      ! Index of front element
    integer :: tail = 0      ! Index of last element (0 = empty)
    integer :: count = 0     ! Number of elements in queue
    integer, allocatable :: data(:)
  contains
    procedure, public, pass :: Initialize => IntQueue_Initialize
    procedure, public, pass :: Enqueue => IntQueue_Enqueue
    procedure, public, pass :: Dequeue => IntQueue_Dequeue
    procedure, public, pass :: Peek => IntQueue_Peek
    procedure, public, pass :: IsEmpty => IntQueue_IsEmpty
    procedure, public, pass :: IsFull => IntQueue_IsFull
    procedure, public, pass :: Size => IntQueue_Size
    procedure, public, pass :: Clear => IntQueue_Clear
    procedure, public, pass :: Contains => IntQueue_Contains
    procedure, public, pass :: Destroy => IntQueue_Destroy
  end type

  !-------------------------------------------------
  ! Integer Stack - useful for DFS, backtracking, etc.
  !-------------------------------------------------
  type, public :: IntStack
    private
    integer :: capacity = 0
    integer :: top = 0       ! Index of top element (0 = empty)
    integer, allocatable :: data(:)
  contains
    procedure, public, pass :: Initialize => IntStack_Initialize
    procedure, public, pass :: Push => IntStack_Push
    procedure, public, pass :: Pop => IntStack_Pop
    procedure, public, pass :: Peek => IntStack_Peek
    procedure, public, pass :: IsEmpty => IntStack_IsEmpty
    procedure, public, pass :: IsFull => IntStack_IsFull
    procedure, public, pass :: Size => IntStack_Size
    procedure, public, pass :: Clear => IntStack_Clear
    procedure, public, pass :: Contains => IntStack_Contains
    procedure, public, pass :: Destroy => IntStack_Destroy
  end type

  !-------------------------------------------------
  ! Regrowth Node - for storing regrowth schedule information
  ! Contains the atom to grow and its "grown" neighbors
  !-------------------------------------------------
  type, public :: RegrowthNode
    integer :: atomID = 0             ! The atom to be regrown
    integer :: nGrownNeighbors = 0    ! Number of already-grown neighbors
    integer :: grownNeighbors(4) = 0  ! IDs of grown neighbors (max 4 for tetrahedral)
  end type

  !-------------------------------------------------
  ! Regrowth Queue - queue of RegrowthNodes for CBMC scheduling
  !-------------------------------------------------
  type, public :: RegrowthQueue
    private
    integer :: capacity = 0
    integer :: head = 1
    integer :: tail = 0
    integer :: count = 0
    type(RegrowthNode), allocatable :: data(:)
  contains
    procedure, public, pass :: Initialize => RegrowthQueue_Initialize
    procedure, public, pass :: Enqueue => RegrowthQueue_Enqueue
    procedure, public, pass :: Dequeue => RegrowthQueue_Dequeue
    procedure, public, pass :: Peek => RegrowthQueue_Peek
    procedure, public, pass :: IsEmpty => RegrowthQueue_IsEmpty
    procedure, public, pass :: Size => RegrowthQueue_Size
    procedure, public, pass :: Clear => RegrowthQueue_Clear
    procedure, public, pass :: Destroy => RegrowthQueue_Destroy
  end type

!==================================================================
contains
!==================================================================
! Integer Queue Implementation
!==================================================================
  subroutine IntQueue_Initialize(self, capacity)
    implicit none
    class(IntQueue), intent(inout) :: self
    integer, intent(in) :: capacity

    if (allocated(self%data)) deallocate(self%data)
    
    self%capacity = capacity
    self%head = 1
    self%tail = 0
    self%count = 0
    allocate(self%data(1:capacity))
    self%data = 0

  end subroutine
!==================================================================
  subroutine IntQueue_Enqueue(self, value, success)
    implicit none
    class(IntQueue), intent(inout) :: self
    integer, intent(in) :: value
    logical, intent(out), optional :: success

    if (self%count >= self%capacity) then
      if (present(success)) success = .false.
      return
    endif

    ! Circular buffer: wrap around
    self%tail = mod(self%tail, self%capacity) + 1
    self%data(self%tail) = value
    self%count = self%count + 1

    if (present(success)) success = .true.

  end subroutine
!==================================================================
  function IntQueue_Dequeue(self, success) result(value)
    implicit none
    class(IntQueue), intent(inout) :: self
    logical, intent(out), optional :: success
    integer :: value

    value = 0
    if (self%count <= 0) then
      if (present(success)) success = .false.
      return
    endif

    value = self%data(self%head)
    ! Circular buffer: wrap around
    self%head = mod(self%head, self%capacity) + 1
    self%count = self%count - 1

    if (present(success)) success = .true.

  end function
!==================================================================
  function IntQueue_Peek(self, success) result(value)
    implicit none
    class(IntQueue), intent(in) :: self
    logical, intent(out), optional :: success
    integer :: value

    value = 0
    if (self%count <= 0) then
      if (present(success)) success = .false.
      return
    endif

    value = self%data(self%head)
    if (present(success)) success = .true.

  end function
!==================================================================
  function IntQueue_IsEmpty(self) result(empty)
    implicit none
    class(IntQueue), intent(in) :: self
    logical :: empty

    empty = (self%count <= 0)

  end function
!==================================================================
  function IntQueue_IsFull(self) result(full)
    implicit none
    class(IntQueue), intent(in) :: self
    logical :: full

    full = (self%count >= self%capacity)

  end function
!==================================================================
  function IntQueue_Size(self) result(n)
    implicit none
    class(IntQueue), intent(in) :: self
    integer :: n

    n = self%count

  end function
!==================================================================
  subroutine IntQueue_Clear(self)
    implicit none
    class(IntQueue), intent(inout) :: self

    self%head = 1
    self%tail = 0
    self%count = 0

  end subroutine
!==================================================================
  function IntQueue_Contains(self, value) result(found)
    implicit none
    class(IntQueue), intent(in) :: self
    integer, intent(in) :: value
    logical :: found
    integer :: i, idx

    found = .false.
    if (self%count <= 0) return

    do i = 1, self%count
      idx = mod(self%head + i - 2, self%capacity) + 1
      if (self%data(idx) == value) then
        found = .true.
        return
      endif
    enddo

  end function
!==================================================================
  subroutine IntQueue_Destroy(self)
    implicit none
    class(IntQueue), intent(inout) :: self

    if (allocated(self%data)) deallocate(self%data)
    self%capacity = 0
    self%head = 1
    self%tail = 0
    self%count = 0

  end subroutine

!==================================================================
! Integer Stack Implementation
!==================================================================
  subroutine IntStack_Initialize(self, capacity)
    implicit none
    class(IntStack), intent(inout) :: self
    integer, intent(in) :: capacity

    if (allocated(self%data)) deallocate(self%data)
    
    self%capacity = capacity
    self%top = 0
    allocate(self%data(1:capacity))
    self%data = 0

  end subroutine
!==================================================================
  subroutine IntStack_Push(self, value, success)
    implicit none
    class(IntStack), intent(inout) :: self
    integer, intent(in) :: value
    logical, intent(out), optional :: success

    if (self%top >= self%capacity) then
      if (present(success)) success = .false.
      return
    endif

    self%top = self%top + 1
    self%data(self%top) = value

    if (present(success)) success = .true.

  end subroutine
!==================================================================
  function IntStack_Pop(self, success) result(value)
    implicit none
    class(IntStack), intent(inout) :: self
    logical, intent(out), optional :: success
    integer :: value

    value = 0
    if (self%top <= 0) then
      if (present(success)) success = .false.
      return
    endif

    value = self%data(self%top)
    self%top = self%top - 1

    if (present(success)) success = .true.

  end function
!==================================================================
  function IntStack_Peek(self, success) result(value)
    implicit none
    class(IntStack), intent(in) :: self
    logical, intent(out), optional :: success
    integer :: value

    value = 0
    if (self%top <= 0) then
      if (present(success)) success = .false.
      return
    endif

    value = self%data(self%top)
    if (present(success)) success = .true.

  end function
!==================================================================
  function IntStack_IsEmpty(self) result(empty)
    implicit none
    class(IntStack), intent(in) :: self
    logical :: empty

    empty = (self%top <= 0)

  end function
!==================================================================
  function IntStack_IsFull(self) result(full)
    implicit none
    class(IntStack), intent(in) :: self
    logical :: full

    full = (self%top >= self%capacity)

  end function
!==================================================================
  function IntStack_Size(self) result(n)
    implicit none
    class(IntStack), intent(in) :: self
    integer :: n

    n = self%top

  end function
!==================================================================
  subroutine IntStack_Clear(self)
    implicit none
    class(IntStack), intent(inout) :: self

    self%top = 0

  end subroutine
!==================================================================
  function IntStack_Contains(self, value) result(found)
    implicit none
    class(IntStack), intent(in) :: self
    integer, intent(in) :: value
    logical :: found
    integer :: i

    found = .false.
    do i = 1, self%top
      if (self%data(i) == value) then
        found = .true.
        return
      endif
    enddo

  end function
!==================================================================
  subroutine IntStack_Destroy(self)
    implicit none
    class(IntStack), intent(inout) :: self

    if (allocated(self%data)) deallocate(self%data)
    self%capacity = 0
    self%top = 0

  end subroutine

!==================================================================
! Regrowth Queue Implementation
!==================================================================
  subroutine RegrowthQueue_Initialize(self, capacity)
    implicit none
    class(RegrowthQueue), intent(inout) :: self
    integer, intent(in) :: capacity

    if (allocated(self%data)) deallocate(self%data)
    
    self%capacity = capacity
    self%head = 1
    self%tail = 0
    self%count = 0
    allocate(self%data(1:capacity))

  end subroutine
!==================================================================
  subroutine RegrowthQueue_Enqueue(self, node, success)
    implicit none
    class(RegrowthQueue), intent(inout) :: self
    type(RegrowthNode), intent(in) :: node
    logical, intent(out), optional :: success

    if (self%count >= self%capacity) then
      if (present(success)) success = .false.
      return
    endif

    self%tail = mod(self%tail, self%capacity) + 1
    self%data(self%tail) = node
    self%count = self%count + 1

    if (present(success)) success = .true.

  end subroutine
!==================================================================
  function RegrowthQueue_Dequeue(self, success) result(node)
    implicit none
    class(RegrowthQueue), intent(inout) :: self
    logical, intent(out), optional :: success
    type(RegrowthNode) :: node

    node%atomID = 0
    node%nGrownNeighbors = 0
    node%grownNeighbors = 0

    if (self%count <= 0) then
      if (present(success)) success = .false.
      return
    endif

    node = self%data(self%head)
    self%head = mod(self%head, self%capacity) + 1
    self%count = self%count - 1

    if (present(success)) success = .true.

  end function
!==================================================================
  function RegrowthQueue_Peek(self, success) result(node)
    implicit none
    class(RegrowthQueue), intent(in) :: self
    logical, intent(out), optional :: success
    type(RegrowthNode) :: node

    node%atomID = 0
    node%nGrownNeighbors = 0
    node%grownNeighbors = 0

    if (self%count <= 0) then
      if (present(success)) success = .false.
      return
    endif

    node = self%data(self%head)
    if (present(success)) success = .true.

  end function
!==================================================================
  function RegrowthQueue_IsEmpty(self) result(empty)
    implicit none
    class(RegrowthQueue), intent(in) :: self
    logical :: empty

    empty = (self%count <= 0)

  end function
!==================================================================
  function RegrowthQueue_Size(self) result(n)
    implicit none
    class(RegrowthQueue), intent(in) :: self
    integer :: n

    n = self%count

  end function
!==================================================================
  subroutine RegrowthQueue_Clear(self)
    implicit none
    class(RegrowthQueue), intent(inout) :: self

    self%head = 1
    self%tail = 0
    self%count = 0

  end subroutine
!==================================================================
  subroutine RegrowthQueue_Destroy(self)
    implicit none
    class(RegrowthQueue), intent(inout) :: self

    if (allocated(self%data)) deallocate(self%data)
    self%capacity = 0
    self%head = 1
    self%tail = 0
    self%count = 0

  end subroutine
!==================================================================
end module
!==================================================================
