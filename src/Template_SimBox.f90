!==========================================================================================
! The purpose of this module is to provide the base class for the simulation box
! family of objects. This type is only intended to provide the basic
! structure for child classes and should not be used directly. 
! Unlike other template classes, the root class used by most box types is the SimpleBox.
!==========================================================================================
module Template_SimBox
  use MasterTemplate, only: classyClass
  use VarPrecision
!  use ForcefieldData, only: ECalcArray
  use Template_NeighList, only: NeighListDef
  use CoordinateTypes, only: Perturbation
!  use ConstraintTemplate, only: constrainArray


  !Sim Box Definition
  type, public, extends(classyClass) :: SimBox
    character(len=15) :: boxStr = "Empty"
    integer :: boxID
    integer :: nAtoms, nMaxAtoms
    integer :: nMolTotal, maxMol
    integer :: nDimension = 3

    !Thermodynamic Variables
    real(dp) :: ETotal, HTotal
    real(dp) :: pressure = 0E0_dp
    real(dp) :: beta, temperature, volume

    real(dp), allocatable :: chempot(:)
    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp), allocatable :: atoms(:,:), centerMass(:,:)

    logical :: forceoutofdate = .true.
    real(dp) :: forcedelta = 1E-6_dp
    real(dp), allocatable :: forces(:,:)

    integer, allocatable :: NMolMin(:), NMolMax(:)
    integer, allocatable :: NMol(:), MolStartIndx(:), MolEndIndx(:)

    integer, allocatable :: AtomType(:), MolType(:)
    integer, allocatable :: MolIndx(:), MolSubIndx(:), AtomSubIndx(:)

    integer, allocatable :: MolGlobalIndx(:, :)
    integer, allocatable :: TypeFirst(:), TypeLast(:)

    integer :: nLists
    class(NeighListDef), allocatable :: NeighList(:)

    contains
      procedure, public, pass :: Constructor
      procedure, public, pass :: LoadAtomCoord
      procedure, public, pass :: LoadDimension
      procedure, public, pass :: GetCoordinates
      procedure, public, pass :: GetAtomTypes
      procedure, public, pass :: GetAtomData
      procedure, public, pass :: GetMolData
      procedure, public, pass :: BuildNeighList
      procedure, public, pass :: Boundary
      procedure, public, pass :: BoundaryNew
      procedure, public, pass :: ComputeEnergy
      procedure, public, pass :: ProcessIO
      procedure, public, pass :: DumpData
      procedure, public, pass :: GetThermo
      procedure, public, pass :: ThermoLookUp
      procedure, public, pass :: IsActive
      procedure, public, pass :: UpdateEnergy
      procedure, public, pass :: UpdatePosition
      procedure, public, pass :: UpdateVol
      procedure, public, pass :: UpdateNeighLists

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
    lineStat = 0
  end subroutine
!==========================================================================================
  subroutine GetCoordinates(self, atoms, slice)
    ! Used to return a pointer array to the atoms(:,:) array within the box.
    ! Can also be used to simply return a subset of the atoms(:,:) such as a single
    ! molecule's coordinates
    implicit none
    class(SimBox), intent(inout), target :: self
    real(dp), pointer :: atoms(:,:)
    integer, optional :: slice(1:2)
    integer :: lb, ub

    if(present(slice)) then
        lb = slice(1)
        ub = slice(2)
        atoms => self%atoms(1:3, lb:ub)
    else
        atoms => self%atoms
    endif

  end subroutine
!==========================================================================================
  subroutine GetAtomTypes(self, Atomtype)
    implicit none
    class(SimBox), intent(inout), target :: self
    integer, pointer, intent(out) :: AtomType(:)

    Atomtype => self%AtomType

  end subroutine
!==========================================================================================
  subroutine LoadDimension(self, line, lineStat)
    implicit none
    class(SimBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    lineStat = 0
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
    real(dp), intent(inout), optional :: rx, ry, rz 

  end subroutine
!==========================================================================================
  subroutine BoundaryNew(self, rx, ry, rz, disp)
    implicit none
    class(SimBox), intent(in) :: self
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inout), optional :: rx, ry, rz 

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
    class(Perturbation), intent(inout) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
  end subroutine
!==========================================================================================
  subroutine UpdateNeighLists(self, disp)
    use CoordinateTypes
    implicit none
    class(SimBox), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
  end subroutine
!==========================================================================================
  subroutine UpdateVol(self, scalar)
    use CoordinateTypes
    implicit none
    class(SimBox), intent(inout) :: self
    real(dp), intent(in) :: scalar(:)
  end subroutine
!==========================================================================================
  subroutine ComputeEnergy(self, tablecheck)
    implicit none
    class(SimBox), intent(inout) :: self
    logical, optional, intent(in) :: tablecheck
  end subroutine
!==========================================================================================
  subroutine ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(SimBox), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line   
    integer, intent(out) :: lineStat

    lineStat = 0
  end subroutine
!==========================================================================================
  subroutine GetMolData(self, globalIndx, nAtoms, molStart, molEnd, molType, atomSubIndx)
    implicit none
    class(SimBox), intent(inout) :: self
    integer, intent(in)  :: globalIndx
    integer, intent(inout), optional :: nAtoms, molStart, molEnd, molType, atomSubIndx


  end subroutine
!==========================================================================================
  subroutine GetAtomData(self, atomglobalIndx, molIndx, atomSubIndx, atomtype)
    implicit none
    class(SimBox), intent(inout) :: self
    integer, intent(in)  :: atomglobalIndx
    integer, intent(inout), optional :: molIndx, atomSubIndx, atomtype


  end subroutine

!==========================================================================================
!  Checks to see if the atom is present in the box.  This is required 
  function IsActive(self, atmIndx) result(active)
    implicit none
    class(SimBox), intent(inout) :: self
    logical :: active
    integer, intent(in) :: atmIndx

    active = .false.

  end function
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
        thermVal = self%ETotal
      case(2) !Ethalpy
        thermVal = self%ETotal + self%Pressure*self%Volume
      case(3) !Energy/Mol
        thermVal = self%ETotal / self%nMolTotal
      case(4) !Enthalpy/Mol
        thermVal = (self%ETotal + self%Pressure*self%Volume)/self%nMolTotal
      case(5) !Volume
        thermVal = self%volume
      case(6) !Temperature
        thermVal = self%temperature
      case(7) !Pressure
        thermVal = self%pressure
      case(8) !Density
        thermVal = self%nMolTotal / self%volume

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

      case("energy/mol") !Energy/Molecule
        thermInt = 3

      case("enthalpy/mol") !Enthalpy/Molecule
        thermInt = 4

      case("volume") !Volume
        thermInt = 5

      case("temperature") !Temperature
        thermInt = 6

      case("pressure") !Pressure
        thermInt = 7

      case("density") !Density
        thermInt = 8

    end select

  end function


!==========================================================================================
end module
!==========================================================================================
