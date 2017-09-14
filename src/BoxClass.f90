module SimBoxDef
use VarPrecision
use CoordinateTypes

  type, public :: SimBox
    real(dp), allocatable :: atoms(:,:)
    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp) :: beta, temperature
    real(dp) :: ETotal
    integer, allocatable :: NMin(:), NMax(:)
    integer, allocatable :: NMol(:)
    integer, allocatable :: AtomType(:)
    integer, allocatable :: MolIndx(:)
    integer :: nTotal
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
