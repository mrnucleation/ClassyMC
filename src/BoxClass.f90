!==========================================================================================
module SimBoxDef
use VarPrecision
use NeighListDef

use CoordinateTypes

  type, public :: SimBox
    integer :: nAtoms
    real(dp), allocatable :: atoms(:,:)

    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp) :: beta, temperature
    real(dp) :: ETotal
    integer, allocatable :: NMolMin(:), NMolMax(:)
    integer, allocatable :: NMol(:)
    integer, allocatable :: AtomType(:)
    integer, allocatable :: MolIndx(:)
    integer :: nTotal, nAtoms
! Constraint Class
    contains
      procedure, pass :: UpdateEnergy
      procedure, pass :: UpdatePosition
      procedure, pass :: CreateNeighList

  end type


  contains


  !------------------------------------------------------------------------------
  subroutine UpdateEnergy(self, E_Diff)
  implicit none
  class(SimBox), intent(in) :: self
  real(dp), intent(in) :: E_Diff

    self % E_Total = self % E_Total + E_Diff
!    self % ETable = self % ETable + self % dETable

  end subroutine

  !------------------------------------------------------------------------------
  subroutine UpdatePosition(self, disp)
  implicit none
  class(SimBox), intent(in) :: self
  type(Displacement), intent(in) :: disp(:)
  integer :: iDisp, dispLen, dispIndx

  dispLen = len(disp)

  do iDisp = 1, dispLen
    dispIndx = disp(iDisp) % atmIndx
    self % atoms(1, dispIndx) = disp(iDisp)%x_new
    self % atoms(2, dispIndx) = disp(iDisp)%y_new
    self % atoms(3, dispIndx) = disp(iDisp)%z_new
  enddo

  end subroutine
  !------------------------------------------------------------------------------

end module
!==========================================================================================
