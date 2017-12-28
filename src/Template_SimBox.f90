!==========================================================================================
! The purpose of this module is to provide the base class for the simulation box
! family of objects. This type is only intended to provide the basic
! structure for child classes and should not be used directly. 
!==========================================================================================
module Template_SimBox
  use VarPrecision
!  use ForcefieldData, only: ECalcArray
  use Template_NeighList, only: NeighListDef
  use CoordinateTypes, only: Displacement
!  use ConstraintTemplate, only: constrainArray


  !Sim Box Definition
  type :: SimBox
    character(len=15) :: boxStr = "Empty"
    integer :: nAtoms
    integer :: nDimension = 3
    real(dp) :: ETotal
    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp) :: beta, temperature
    real(dp), allocatable :: atoms(:,:)
    integer, allocatable :: AtomType(:)

    integer :: nLists
    class(NeighListDef), allocatable :: NeighList(:)

    contains
      procedure, public, pass :: Constructor
      procedure, public, pass :: LoadAtomCoord
      procedure, public, pass :: LoadDimension
      procedure, public, pass :: BuildNeighList
      procedure, public, pass :: Boundary
      procedure, public, pass :: UpdateEnergy
      procedure, public, pass :: UpdatePosition
      procedure, public, pass :: UpdateNeighLists
      procedure, public, pass :: ComputeEnergy
      procedure, public, pass :: IOProcess
  end type

  public :: Constructor, LoadCoordinates, BuildNeighList,Boundary
  public :: UpdateEnergy, UpdatePosition, UpdateNeighLists, DummyCoords, IOProcess
  public :: DumpXYZConfig
!==========================================================================================
  contains
!==========================================================================================
  subroutine Constructor(self)
    implicit none
    class(SimBox), intent(inout) :: self
  end subroutine

!==========================================================================================
  subroutine LoadAtomCoord(self, line, lineStat)
    implicit none
    class(SimBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
  end subroutine
!==========================================================================================
  subroutine LoadDimension(self, line, lineStat)
    implicit none
    class(SimBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
  end subroutine
!==========================================================================================
  subroutine BuildNeighList(self)
    implicit none
    class(SimBox), intent(inout) :: self

  end subroutine
!==========================================================================================
  subroutine Boundary(self, rx, ry, rz)
    implicit none
    class(SimBox), intent(in) :: self
    real(dp), intent(inout) :: rx, ry, rz 

  end subroutine
!==========================================================================================
  subroutine UpdateEnergy(self, E_Diff)
    implicit none
    class(SimBox), intent(inout) :: self
    real(dp), intent(in) :: E_Diff
  end subroutine
!==========================================================================================
  subroutine UpdatePosition(self, disp)
    implicit none
    class(SimBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
  end subroutine
!==========================================================================================
  subroutine UpdateNeighLists(self, disp)
    use CoordinateTypes
    implicit none
    class(SimBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
  end subroutine
!==========================================================================================
  subroutine ComputeEnergy(self)
    implicit none
    class(SimBox), intent(inout) :: self
  end subroutine
!==========================================================================================
  subroutine IOProcess(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(SimBox), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line   
    integer, intent(out) :: lineStat

    lineStat = 0
  end subroutine
!==========================================================================================
end module
!==========================================================================================
