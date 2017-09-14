module NeighListDef
use VarPrecision
use CoordinateTypes
use Constants, only: pi

  type, public :: NeighList
      
      integer, allocatable :: list(:,:)
      integer, allocatable :: nNeigh(:)
      integer :: maxNei
      real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor
      procedure, pass :: InitializeList
  end type


  contains

  !------------------------------------------------------------------------------
  subroutine Constructor(self, nAtoms, rCut)
  implicit none
  class(NeighList), intent(in) :: self
  integer, intent(in) :: nAtoms
  real(dp), intent(in) :: rCut
  real(dp), parameter :: atomRadius = 1.0  !Used to estimate an approximate volume of 
  integer :: AllocateStatus

  self % rCut = rCut
  self % rCutSq = rCut * rCut

  self% maxNei = ceiling(atomRadius**3/rCut**3)

  allocate(self%list(1:nAtoms, 1:maxNei), status=AllocateStatus)
  allocate(self%nNei(1:nAtoms), status=AllocateStatus)

  self%list = 0

  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !------------------------------------------------------------------------------
  subroutine InitializeList(self)
  implicit none
  class(NeighList), intent(in) :: self
  real(dp), intent(in) :: E_Diff

    self % E_Total = self % E_Total + E_Diff

  end subroutine
  !------------------------------------------------------------------------------

end module
