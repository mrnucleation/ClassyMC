!==========================================================================================
module SimpleSimBox
  use NeighListDef, only: NeighList
  use VarPrecision
!  use ForcefieldData, only: ECalcArray
!  use NeighListDef
  use CoordinateTypes, only: Displacement
  use ForcefieldData, only: ECalcArray
  use Template_SimBox
!  use ConstraintTemplate, only: constrainArray


  !Sim Box Definition
  type, public, extends(SimBox) :: SimpleBox
    integer :: boxID
    integer :: nTotal
    integer :: nMaxAtoms

!    real(dp), allocatable :: atoms(:,:)
!    real(dp), allocatable :: ETable(:), dETable(:)
!    real(dp) :: beta, temperature

!    real(dp) :: ETotal
    integer, allocatable :: NMolMin(:), NMolMax(:)
    integer, allocatable :: NMol(:)
!    integer, allocatable :: AtomType(:)
    integer, allocatable :: MolIndx(:)

    class(ECalcArray), pointer :: EFunc
    integer :: ECalcer = -1
    integer :: ENeiList = -1
    type(NeighList), allocatable :: NeighList(:)

    contains
      procedure, pass :: Constructor => SimpleBox_Constructor
      procedure, pass :: LoadCoordinates => SimpleBox_LoadCoordinates
      procedure, pass :: BuildNeighList => SimpleBox_BuildNeighList
      procedure, pass :: Boundary => SimpleBox_Boundary
      procedure, pass :: UpdateEnergy => SimpleBox_UpdateEnergy
      procedure, pass :: UpdatePosition => SimpleBox_UpdatePosition
      procedure, pass :: UpdateNeighLists => SimpleBox_UpdateNeighLists
      procedure, pass :: DummyCoords => SimpleBox_DummyCoords
      procedure, pass :: IOProcess => SimpleBox_IOProcess
      procedure, pass :: DumpXYZConfig => SimpleBox_DumpXYZConfig
  end type

!==========================================================================================
  contains
!==========================================================================================
  subroutine SimpleBox_Constructor(self)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: AllocateStatus
 
    allocate(self%atoms(1:3, 1:self%nMaxAtoms), stat= AllocateStatus)
    allocate(self%ETable(1:self%nMaxAtoms), stat= AllocateStatus)
    allocate(self%dETable(1:self%nMaxAtoms), stat= AllocateStatus)
    allocate(self%AtomType(1:self%nMaxAtoms), stat= AllocateStatus)


    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine

!==========================================================================================
  subroutine SimpleBox_LoadCoordinates(self, fileName, fileType)
    implicit none
    class(SimpleBox), intent(inout) :: self
    character(len=*), intent(in) :: fileName
    character(len=*), intent(in), optional :: fileType
    character(len=3) :: sym
    integer :: iAtom
    integer :: AllocateStatus, IOSt

    open(unit=50, file=fileName, status="OLD")
    read(50,*) self%nAtoms
    read(50,*) 

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
  subroutine SimpleBox_BuildNeighList(self)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: iList
    integer :: iAtom, jAtom
    real(dp) :: rx, ry, rz, rsq

    do iAtom = 1, self%nAtoms-1
      do jAtom = iAtom+1, self%nAtoms
        rx = self%atoms(1, iAtom) - self%atoms(1, jAtom)
        ry = self%atoms(2, iAtom) - self%atoms(2, jAtom)
        rz = self%atoms(3, iAtom) - self%atoms(3, jAtom)
        call self%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        do iList = 1, size(self%NeighList)
          if( rsq <= self%NeighList(iList)%rCutSq ) then 
            self%NeighList(iList)%nNeigh(iAtom) = self%NeighList(iList)%nNeigh(iAtom) + 1
            self%NeighList(iList)%list( self%NeighList(iList)%nNeigh(iAtom), iAtom ) = jAtom

            self%NeighList(iList)%nNeigh(jAtom) = self%NeighList(iList)%nNeigh(jAtom) + 1
            self%NeighList(iList)%list( self%NeighList(iList)%nNeigh(jAtom), jAtom ) = iAtom
          endif
        enddo        
      enddo
    enddo

  end subroutine
!==========================================================================================
  subroutine SimpleBox_Boundary(self, rx, ry, rz)
    implicit none
    class(SimpleBox), intent(in) :: self
    real(dp), intent(inout) :: rx, ry, rz 

  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdateEnergy(self, E_Diff)
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: E_Diff

    self % ETotal = self % ETotal + E_Diff
    self % ETable = self % ETable + self % dETable

  end subroutine

!==========================================================================================
  subroutine SimpleBox_UpdatePosition(self, disp)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
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
  subroutine SimpleBox_UpdateNeighLists(self, disp)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
    integer :: iDisp, iList
    integer :: atmIndx, iAtom, jAtom
    real(dp) :: rx, ry, rz, rsq

    do iDisp = 1, size(disp)
      atmIndx = disp(iDisp)%atmIndx
      rx = disp(iDisp)%x_new - self%atoms(1, atmIndx)
      ry = disp(iDisp)%y_new - self%atoms(2, atmIndx)
      rz = disp(iDisp)%z_new - self%atoms(3, atmIndx)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < 1.0E0_dp) then
        cycle
      endif

      do jAtom = 1, self%nAtoms
        rx = self%atoms(1, atmIndx) - self%atoms(1, jAtom)
        ry = self%atoms(2, atmIndx) - self%atoms(2, jAtom)
        rz = self%atoms(3, atmIndx) - self%atoms(3, jAtom)
        call self%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        do iList = 1, size(self%NeighList)
          if( rsq <= self%NeighList(iList)%rCutSq ) then 
            self%NeighList(iList)%nNeigh(iAtom) = self%NeighList(iList)%nNeigh(iAtom) + 1
            self%NeighList(iList)%list( self%NeighList(iList)%nNeigh(iAtom), iAtom ) = jAtom

            self%NeighList(iList)%nNeigh(jAtom) = self%NeighList(iList)%nNeigh(jAtom) + 1
            self%NeighList(iList)%list( self%NeighList(iList)%nNeigh(jAtom), jAtom ) = iAtom
          endif
        enddo        
      enddo
    enddo

  end subroutine

!==========================================================================================
  subroutine SimpleBox_DummyCoords(self)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self

    self % nAtoms = 2

    self % atoms(1, 1) = 0.0
    self % atoms(2, 1) = 0.0
    self % atoms(3, 1) = 0.0
 
    self % atoms(1, 2) = 2.0**(1.0/6.0)
    self % atoms(2, 2) = 0.0
    self % atoms(3, 2) = 0.0

  end subroutine

!==========================================================================================
  subroutine SimpleBox_IOProcess(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    implicit none

    class(SimpleBox), intent(inout) :: self
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
        self % EFunc => EnergyCalculator(intVal)
      case default
        lineStat = -1
    end select
  end subroutine
!==========================================================================================
  subroutine SimpleBox_DumpXYZConfig(self, fileName)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    implicit none
    class(SimpleBox), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: fileName

  end subroutine
!==========================================================================================
end module
!==========================================================================================
