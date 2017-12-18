!==========================================================================================
module SimpleSimBox
!  use NeighListDef, only: NeighList
  use VarPrecision
!  use ForcefieldData, only: ECalcArray
!  use NeighListDef
  use CoordinateTypes, only: Displacement
  use ForcefieldData, only: ECalcArray
  use ConstraintTemplate, only: constrainArray
  use Template_SimBox, only: SimBox


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
    integer, allocatable :: NMol(:), MolStartIndx(:)
!    integer, allocatable :: AtomType(:)
    integer, allocatable :: MolIndx(:), SubIndx(:)

    type(constrainArray), allocatable :: Constrain(:)
    class(ECalcArray), pointer :: EFunc
!    class(NeighList), allocatable :: NeighList(:)

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
    use Common_MolInfo
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: AllocateStatus
    integer :: iType, iMol, iAtom, atmIndx, molIndx, maxMol
 
    if( .not. allocated(self%NMolMin) ) then
      write(*,*) "ERROR! The maximum and minimum molecules allowed in the box must be defined"
      write(*,*) "prior to box initialization!"
      stop 
    endif

    !First begin by computing the maximium number of atoms that the box can potentially contain
    self%nMaxAtoms = 0
    maxMol = 0
    do iType = 1, nMolTypes
      self%nMaxAtoms = self%nMaxAtoms + self%NMolMax(iType)*MolData(iType)%nAtoms
      maxMol = maxMol + self%NMolMax(iType)
    enddo

    !Allocate the position and energy related arrays. 
    allocate(self%atoms(1:3, 1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%ETable(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%dETable(1:self%nMaxAtoms), stat=AllocateStatus)

    !Allocate the arrays which contain the atom type and quick look up information.
    allocate(self%AtomType(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%MolIndx(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%SubIndx(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%MolStartIndx(1:self%nMaxAtoms), stat=AllocateStatus)

    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    self%AtomType = 0
    self%MolIndx = 0
    self%SubIndx = 0
    self%MolStartIndx = 0

    atmIndx = 0
    molIndx = 0
    do iType = 1, nMolTypes
      do iMol = 1, self%NMolMax(iType)
        molIndx = molIndx + 1
        self%MolStartIndx(molIndx) = atmIndx + 1
        do iAtom = 1, MolData(iType)%nAtoms
          atmIndx = atmIndx + 1
          self%AtomType(atmIndx) = MolData(iType)%atomType(iAtom)
          self%MolIndx(atmIndx)  = molIndx
          self%SubIndx(atmIndx)  = iAtom
        enddo
      enddo 
    enddo


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
    read(50,*) self%nMaxAtoms
    read(50,*) 

    if( .not. allocated(self%nAtoms) ) then
      call self%Constructor
    endif

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

!    write(*,*) "Here"
    do iList = 1, size(self%NeighList)
      self%NeighList(iList)%nNeigh = 0
      self%NeighList(iList)%list = 0
!      write(*,*) iList, size(self%NeighList(iList)%nNeigh)
    enddo

    do iAtom = 1, self%nAtoms-1
      do jAtom = iAtom+1, self%nAtoms
        rx = self%atoms(1, iAtom) - self%atoms(1, jAtom)
        ry = self%atoms(2, iAtom) - self%atoms(2, jAtom)
        rz = self%atoms(3, iAtom) - self%atoms(3, jAtom)
        call self%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) iAtom, jAtom, rsq
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
!    write(*,*) "Here"


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
      case("molmin")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal

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
