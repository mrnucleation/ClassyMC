module NeighListDef
use VarPrecision
use CoordinateTypes

  type, public :: NeighList
      logical :: Strict = .false.
      integer, allocatable :: list(:,:)
      integer, allocatable :: nNeigh(:)
      integer :: maxNei
      real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor
      procedure, pass :: InitializeList
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine Constructor(self, nAtoms, rCut)
    implicit none
    class(NeighList), intent(inout) :: self
    integer, intent(in) :: nAtoms
    real(dp), intent(in) :: rCut
    real(dp), parameter :: atomRadius = 1.0  !Used to estimate an approximate volume of 
    integer :: AllocateStatus

    self % rCut = rCut
    self % rCutSq = rCut * rCut

    self% maxNei = ceiling(atomRadius**3/rCut**3)

    allocate( self%list(1:self%maxNei, 1:nAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:nAtoms), stat=AllocateStatus )

    self%list = 0

    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine
!===================================================================================
  subroutine InitializeList(self)
  implicit none
  class(NeighList), intent(inout) :: self

  end subroutine
!===================================================================================
end module
!===================================================================================
