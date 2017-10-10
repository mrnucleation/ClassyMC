!==========================================================================================
module SimBoxDef
use VarPrecision
!use ForcefieldData, only: ECalcArray
!use NeighListDef

use CoordinateTypes

  type, public :: SimBox
    integer :: nAtoms
    real(dp), allocatable :: atoms(:,:)

    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp) :: beta, temperature
    real(dp) :: ETotal
    integer, allocatable :: AtomType(:)
!    type(ECalcArray), pointer :: ECalc

! Constraint Class
    contains
      procedure, pass :: Constructor
      procedure, pass :: LoadCoordinates
      procedure, pass :: UpdateEnergy
      procedure, pass :: UpdatePosition
!      procedure, pass :: CreateNeighList

  end type


  contains

  !------------------------------------------------------------------------------
  subroutine Constructor(self)
    implicit none
    class(SimBox), intent(inout) :: self
!    character(len=*), intent(in) :: fileName


  end subroutine

  !------------------------------------------------------------------------------
  subroutine LoadCoordinates(self, fileName)
  implicit none
  class(SimBox), intent(inout) :: self
  character(len=*), intent(in) :: fileName


  end subroutine
  !------------------------------------------------------------------------------
  subroutine UpdateEnergy(self, E_Diff)
  implicit none
  class(SimBox), intent(inout) :: self
  real(dp), intent(in) :: E_Diff

    self % ETotal = self % ETotal + E_Diff
!    self % ETable = self % ETable + self % dETable

  end subroutine

  !------------------------------------------------------------------------------
  subroutine UpdatePosition(self, disp)
  use CoordinateTypes
  implicit none
  class(SimBox), intent(inout) :: self
  type(Displacement), intent(in) :: disp(:)
  integer :: iDisp, dispLen, dispIndx

  dispLen = size(disp)

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
