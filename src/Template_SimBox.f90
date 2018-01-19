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
  type, public :: SimBox
    character(len=15) :: boxStr = "Empty"
    integer :: nAtoms, nMaxAtoms
    integer :: nMolTotal
    integer :: nDimension = 3

    !Thermodynamic Variables
    real(dp) :: ETotal, HTotal
    real(dp) :: beta, temperature, pressure, volume
    real(dp), allocatable :: chempot(:)

    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp), allocatable :: atoms(:,:)

    integer, allocatable :: NMolMin(:), NMolMax(:)
    integer, allocatable :: NMol(:), MolStartIndx(:), MolEndIndx(:)

    integer, allocatable :: AtomType(:), MolType(:)
    integer, allocatable :: MolIndx(:), MolSubIndx(:), SubIndx(:)

    integer, allocatable :: TypeFirst(:), TypeLast(:)

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
      procedure, public, pass :: DumpData
      procedure, public, pass :: GetThermo
      procedure, public, pass :: ThermoLookUp
  end type
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
  subroutine UpdatePosition(self, disp, tempList, tempNNei)
    implicit none
    class(SimBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
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
  subroutine DumpData(self, filename)
    use Input_Format, only: maxLineLen
    implicit none
    class(SimBox), intent(inout) :: self
    character(len=*), intent(in) :: filename

  end subroutine
!==========================================================================================
  function GetThermo(self, thermInt) result(thermVal)
    use Input_Format, only: maxLineLen
    implicit none
    class(SimBox), intent(in) :: self
    integer, intent(in) :: thermInt
    real(dp) :: thermVal

    select case(thermInt)
      case(1) !Energy
        thermVal = self % ETotal
      case(2) !Ethalpy
        thermVal = self % HTotal
      case(3) !Volume
        thermVal = self % volume
      case(4) !Temperature
        thermVal = self % temperature
      case(5) !Pressure
        thermVal = self % pressure
    end select

  end function

!==========================================================================================
  function ThermoLookup(self, thermoStr) result(thermInt)
    use Input_Format, only: maxLineLen
    implicit none
    class(SimBox), intent(in) :: self
    character(len=30), intent(in) :: thermoStr
    integer :: thermInt

    select case(trim(adjustl(thermoStr)))
      case("energy") !Energy
        thermInt = 1
      case("enthalpy") !Ethalpy
        thermInt = 2
      case("volume") !Volume
        thermInt = 3
      case("temperature") !Temperature
        thermInt = 4
      case("pressure") !Pressure
        thermInt = 5
    end select

  end function


!==========================================================================================
end module
!==========================================================================================
