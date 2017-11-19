!==========================================================================================
module CubicBoxDef
  use Template_SimBox
  use SimpleSimBox
  use NeighListDef, only: NeighList
  use VarPrecision
!  use ForcefieldData, only: ECalcArray
!  use NeighListDef
  use CoordinateTypes, only: Displacement
!  use ConstraintTemplate, only: constrainArray


  !Sim Box Definition
  type, public, extends(SimpleBox) :: CubeBox
    real(dp) :: boxL, boxL2

    contains
      procedure, pass :: Constructor => Cube_Constructor
      procedure, pass :: LoadCoordinates => Cube_LoadCoordinates
      procedure, pass :: UpdateEnergy => Cube_UpdateEnergy
      procedure, pass :: Boundary => Cube_Boundary
      procedure, pass :: UpdatePosition => Cube_UpdatePosition
      procedure, pass :: DummyCoords => Cube_DummyCoords
      procedure, pass :: IOProcess => Cube_IOProcess
      procedure, pass :: DumpXYZConfig => Cube_DumpXYZConfig
  end type

!==========================================================================================
  contains
!==========================================================================================
  subroutine Cube_Constructor(self)
    implicit none
    class(CubeBox), intent(inout) :: self
    integer :: AllocateStatus
 
    write(*,*) "Box Type: Cubic"

    allocate(self%atoms(1:3, 1:self%nAtoms), stat= AllocateStatus)
    allocate(self%ETable(1:self%nAtoms), stat= AllocateStatus)
    allocate(self%dETable(1:self%nAtoms), stat= AllocateStatus)
    allocate(self%AtomType(1:self%nAtoms), stat= AllocateStatus)


    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine

!==========================================================================================
  subroutine Cube_LoadCoordinates(self, fileName, fileType)
    implicit none
    class(CubeBox), intent(inout) :: self
    character(len=*), intent(in) :: fileName
    character(len=*), intent(in), optional :: fileType
    character(len=3) :: sym
    integer :: iAtom
    integer :: AllocateStatus, IOSt

    open(unit=50, file=fileName, status="OLD")
    read(50,*) self%nAtoms
    read(50,*) self%boxL

    self%boxL2 = self%boxL/2.0

    allocate( self%atoms(1:3, 1:self%nAtoms), stat= AllocateStatus)
    allocate( self%AtomType(1:self%nAtoms), stat= AllocateStatus )
    allocate( self%ETable(1:self%nAtoms), stat= AllocateStatus )
    allocate( self%dETable(1:self%nAtoms), stat= AllocateStatus )
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    do iAtom = 1, self%nAtoms
      read(50,*, iostat=IOSt) sym, self%atoms(1, iAtom), self%atoms(2, iAtom), self%atoms(3, iAtom)
      if(IOSt .ne. 0) then
        write(*,*) "ERROR! Could not properly read the configuration file."
        stop
      endif
    enddo

    close(50)

  end subroutine
!==========================================================================================
  subroutine Cube_Boundary(self, rx, ry, rz)
  implicit none
  class(CubeBox), intent(in) :: self
  real(dp), intent(inout) :: rx, ry, rz 
  real(dp) :: rx_new, ry_new, rz_new

  if(abs(rx) > self%boxL2) then
    rx = rx - sign(self%boxL, rx)
  endif

  if(abs(ry) > self%boxL2) then
    ry = ry - sign(self%boxL, ry)
  endif

  if(abs(rz) > self%boxL2) then
    rz = rz - sign(self%boxL, rz)
  endif

  end subroutine
!==========================================================================================
  subroutine Cube_UpdateEnergy(self, E_Diff)
  implicit none
  class(CubeBox), intent(inout) :: self
  real(dp), intent(in) :: E_Diff

    self % ETotal = self % ETotal + E_Diff
    self % ETable = self % ETable + self % dETable

  end subroutine

!==========================================================================================
  subroutine Cube_UpdatePosition(self, disp)
    use CoordinateTypes
    implicit none
    class(CubeBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
    integer :: iDisp, dispLen, dispIndx


!  write(*,*) "cube"

    dispLen = size(disp)

    do iDisp = 1, dispLen
      dispIndx = disp(iDisp) % atmIndx
      call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
      self % atoms(1, dispIndx) = disp(iDisp)%x_new
      self % atoms(2, dispIndx) = disp(iDisp)%y_new
      self % atoms(3, dispIndx) = disp(iDisp)%z_new
    enddo

  end subroutine
!==========================================================================================
  subroutine Cube_DummyCoords(self)
    use CoordinateTypes
    implicit none
    class(CubeBox), intent(inout) :: self

    self % nAtoms = 2

    self % atoms(1, 1) = 0.0
    self % atoms(2, 1) = 0.0
    self % atoms(3, 1) = 0.0
 
    self % atoms(1, 2) = 2.0**(1.0/6.0)
    self % atoms(2, 2) = 0.0
    self % atoms(3, 2) = 0.0

  end subroutine

!==========================================================================================
  subroutine Cube_IOProcess(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    implicit none

    class(CubeBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: intVal
    real(dp) :: realVal
    character(len=30) :: command, val

!    write(*,*) "Cube"
    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    call LowerCaseLine(command)
    write(*,*) command
    select case( trim(adjustl(command)) )
      case("energycalc")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % EFunc => EnergyCalculator(intVal)

!        self % ECalcer = intVal
      case("boxlength")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % boxL = realVal
        self % boxL2 = realVal/2.0E0_dp
        write(*,*) "Box Length:", self%boxL
      case default
        lineStat = -1
    end select

  end subroutine
!==========================================================================================
  subroutine Cube_DumpXYZConfig(self, fileName)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    implicit none
    class(CubeBox), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: fileName


  end subroutine
!==========================================================================================
end module
!==========================================================================================
