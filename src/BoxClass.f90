!==========================================================================================
module SimBoxDef
use NeighListDef, only: NeighList
use VarPrecision
!use ForcefieldData, only: ECalcArray
!use NeighListDef
use CoordinateTypes

  type, public :: SimBox
    integer :: nTotal, nAtoms
    real(dp), allocatable :: atoms(:,:)

    real(dp), allocatable :: ETable(:), dETable(:)
    real(dp) :: beta, temperature
    real(dp) :: ETotal
    integer, allocatable :: NMolMin(:), NMolMax(:)
    integer, allocatable :: NMol(:)
    integer, allocatable :: AtomType(:)
    integer, allocatable :: MolIndx(:)

    integer :: ECalcer = -1

    type(NeighList), allocatable :: NeighList(:)

! Constraint Class
    contains
      procedure, pass :: Constructor
      procedure, pass :: LoadCoordinates
      procedure, pass :: UpdateEnergy
      procedure, pass :: UpdatePosition
      procedure, pass :: DummyCoords
      procedure, pass :: IOProcess
!      procedure, pass :: CreateNeighList

  end type


  contains

  !==========================================================================================
  subroutine Constructor(self)
    implicit none
    class(SimBox), intent(inout) :: self
    integer :: AllocateStatus
 
    allocate(self%atoms(1:3, 1:self%nAtoms), stat= AllocateStatus)
    allocate(self%ETable(1:self%nAtoms), stat= AllocateStatus)
    allocate(self%dETable(1:self%nAtoms), stat= AllocateStatus)
    allocate(self%AtomType(1:self%nAtoms), stat= AllocateStatus)


    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine

  !==========================================================================================
  subroutine LoadCoordinates(self, fileName)
  implicit none
  class(SimBox), intent(inout) :: self
  character(len=*), intent(in) :: fileName


  end subroutine
  !==========================================================================================
  subroutine UpdateEnergy(self, E_Diff)
  implicit none
  class(SimBox), intent(inout) :: self
  real(dp), intent(in) :: E_Diff

    self % ETotal = self % ETotal + E_Diff
!    self % ETable = self % ETable + self % dETable

  end subroutine

  !==========================================================================================
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
  !==========================================================================================
  subroutine DummyCoords(self)
    use CoordinateTypes
    implicit none
    class(SimBox), intent(inout) :: self

    self % nAtoms = 2

    self % atoms(1, 1) = 0.0
    self % atoms(2, 1) = 0.0
    self % atoms(3, 1) = 0.0
 
    self % atoms(1, 2) = 2.0**(1.0/6.0)
    self % atoms(2, 2) = 0.0
    self % atoms(3, 2) = 0.0

  end subroutine

  !==========================================================================================
  subroutine IOProcess(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    implicit none

    class(SimBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: intVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    call LowerCaseLine(command)
    write(*,*) command
    select case( trim(adjustl(command)) )
      case("energycalc")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % ECalcer = intVal
      case default
        lineStat = -1
    end select
  end subroutine
  !==========================================================================================

end module
!==========================================================================================
