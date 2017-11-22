!===================================================================================
! This module contains a simplistic neighbor list implimentation primarily used for debugging.  
!===================================================================================
module RSqListDef
use VarPrecision
use CoordinateTypes
use Template_SimBox, only: SimBox
use SimpleSimBox, only: SimpleBox
use Template_NeighList, only: NeighListDef

  type, public, extends(NeighListDef) :: RSqList
!      logical :: Strict = .false.
!      integer, allocatable :: list(:,:)
!      integer, allocatable :: nNeigh(:)
!      integer :: maxNei
!      real(dp) :: rCut, rCutSq
      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => RSqList_Constructor 
      procedure, pass :: InitializeList => RSqList_InitializeList 
      procedure, pass :: GetNewList => RSqList_GetNewList
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine RSqList_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use NeighData, only: neighSkin
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 1.0E0_dp  !Used to estimate an approximate volume of 
    integer :: AllocateStatus

    self%parent => BoxArray(parentID)%box
    if(.not. allocated(self%parent%atoms) ) then
      stop
    endif

!     If no rCut value is given by the subroutine call attempt to pull
!     the rSq value from the parent box's energy function. The assumption being
!     that the neighborlist is used for the energy calculation routines.
    if( present(rCut) ) then
      self % rCut = rCut
      self % rCutSq = rCut * rCut
      self % maxNei = ceiling(rCut**3/atomRadius**3)
    else
      self % rCut = self % parent % EFunc % Method % GetCutOff() + neighSkin
      write(*,*) self % rCut
      self % rCutSq = (self%rCut)**2
      self % maxNei = ceiling(self%rCut**3/atomRadius**3)
    endif

    allocate( self%list(1:self%maxNei, 1:self%parent%nAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:self%parent%nAtoms), stat=AllocateStatus )

    self%list = 0

    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine
!===================================================================================
  subroutine RSqList_InitializeList(self)
    implicit none
    class(RSqList), intent(inout) :: self

  end subroutine
!===================================================================================
  subroutine RSqList_GetNewList(self, disp, newList)
    implicit none
    class(RSqList), intent(inout) :: self
    type(Displacement), intent(in) :: disp
    real(dp), intent(out) :: newList(:)

    newList = 0E0_dp
  end subroutine
!===================================================================================
end module
!===================================================================================
