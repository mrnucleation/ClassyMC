!==========================================================================================
module SimpleSimBox
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use ForcefieldData, only: ECalcArray
  use ConstraintTemplate, only: constrainArray
  use Template_SimBox, only: SimBox


  !Sim Box Definition
  type, public, extends(SimBox) :: SimpleBox
    integer :: boxID
    integer :: nTotal
!    integer :: nMaxAtoms

!    real(dp), allocatable :: atoms(:,:)
!    real(dp), allocatable :: ETable(:), dETable(:)

!    real(dp) :: pressure = 0E0_dp
!    real(dp) :: beta, temperature, volume
!    real(dp), allocatable :: chempot(:)

!    real(dp) :: ETotal
!    integer, allocatable :: NMolMin(:), NMolMax(:)
!    integer, allocatable :: NMol(:), MolStartIndx(:), MolEndIndx(:)

!    integer, allocatable :: AtomType(:)
!    integer, allocatable :: MolIndx(:), SubIndx(:)


    type(constrainArray), allocatable :: Constrain(:)
    class(ECalcArray), pointer :: EFunc
!    class(NeighList), allocatable :: NeighList(:)


    contains
      procedure, pass :: Constructor => SimpleBox_Constructor
      procedure, pass :: AllocateMolBound => SimpleBox_AllocateMolBound
      procedure, pass :: LoadAtomCoord => Simplebox_LoadAtomCoord
      procedure, pass :: LoadDimension => Simplebox_LoadDimension
      procedure, pass :: BuildNeighList => SimpleBox_BuildNeighList
      procedure, pass :: Boundary => SimpleBox_Boundary
      procedure, pass :: ComputeEnergy => SimpleBox_ComputeEnergy
      procedure, pass :: IOProcess => SimpleBox_IOProcess
      procedure, pass :: CheckConstraint => SimpleBox_CheckConstraint
      procedure, pass :: DumpData => SimpleBox_DumpData

      procedure, pass :: AddMol => SimpleBox_AddMol
      procedure, pass :: DeleteMol => SimpleBox_DeleteMol
      procedure, pass :: UpdateEnergy => SimpleBox_UpdateEnergy
      procedure, pass :: UpdatePosition => SimpleBox_UpdatePosition
      procedure, pass :: UpdateNeighLists => SimpleBox_UpdateNeighLists

!      procedure, public, pass :: GetThermo
!      procedure, public, pass :: ThermoLookUp

      procedure, pass :: Maintenence => SimpleBox_Maintenence
      procedure, pass :: Prologue => SimpleBox_Prologue
      procedure, pass :: Epilogue => SimpleBox_Epilogue

!      GENERIC, PUBLIC :: ASSIGNMENT(=) => SimpleBox_CopyBox

  end type

  interface assignment(=)
    module procedure SimpleBox_CopyBox
  end interface

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

    self%boxStr = "NoBox"
    !First begin by computing the maximium number of atoms that the box can potentially contain
    self%nMaxAtoms = 0
    maxMol = 0
    do iType = 1, nMolTypes
      self%nMaxAtoms = self%nMaxAtoms + self%NMolMax(iType)*MolData(iType)%nAtoms
      maxMol = maxMol + self%NMolMax(iType)
    enddo

    self%nAtoms = 0
    do iType = 1, nMolTypes
      self%nAtoms = self%nAtoms + self%NMol(iType)*MolData(iType)%nAtoms
    enddo

    !Allocate the position and energy related arrays. 
    allocate(self%atoms(1:3, 1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%ETable(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%dETable(1:self%nMaxAtoms), stat=AllocateStatus)

    !Allocate the arrays which contain the atom type and quick look up information.
    allocate(self%AtomType(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%MolType(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%MolIndx(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%MolSubIndx(1:self%nMaxAtoms), stat=AllocateStatus)
    allocate(self%SubIndx(1:self%nMaxAtoms), stat=AllocateStatus)

    allocate(self%MolStartIndx(1:maxMol), stat=AllocateStatus)
    allocate(self%MolEndIndx(1:maxMol), stat=AllocateStatus)

    allocate(self%TypeFirst(1:nMolTypes), stat=AllocateStatus)
    allocate(self%TypeLast(1:nMolTypes), stat=AllocateStatus)

    maxMol = maxval(self%NMolMax(:))
    allocate(self%MolGlobalIndx(1:nMolTypes, 1:maxMol), stat=AllocateStatus)

    allocate(self%chempot(1:nMolTypes), stat=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    self%AtomType = 0
    self%MolType = 0
    self%MolIndx = 0
    self%MolSubIndx = 0
    self%SubIndx = 0
    self%MolStartIndx = 0
    self%MolEndIndx = 0

    self%TypeFirst = 0
    self%TypeLast = 0

    self%MolGlobalIndx = 0

    self%chempot = 0E0_dp

    atmIndx = 0
    molIndx = 0
    do iType = 1, nMolTypes
      self%TypeFirst(iType) = atmIndx + 1
      do iMol = 1, self%NMolMax(iType)
        molIndx = molIndx + 1
        self%MolGlobalIndx(iType, iMol) = molIndx
        self%MolStartIndx(molIndx) = atmIndx + 1
        self%MolEndIndx(molIndx) = atmIndx + MolData(iType)%nAtoms 
        do iAtom = 1, MolData(iType)%nAtoms
          atmIndx = atmIndx + 1
          self%MolType(atmIndx) = iType
          self%AtomType(atmIndx) = MolData(iType)%atomType(iAtom)
          self%MolIndx(atmIndx)  = molIndx
          self%MolSubIndx(atmIndx)  = iMol
          self%SubIndx(atmIndx)  = iAtom
        enddo
      enddo 
      self%TypeLast(iType) = atmIndx 
    enddo


  end subroutine
!==========================================================================================
  subroutine SimpleBox_AllocateMolBound(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: AllocateStatus
 
    allocate(self%NMol(1:nMolTypes), stat=AllocateStatus)
    allocate(self%NMolMax(1:nMolTypes), stat=AllocateStatus)
    allocate(self%NMolMin(1:nMolTypes), stat=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  end subroutine
!==========================================================================================
  subroutine Simplebox_LoadDimension(self, line, lineStat)
    use Input_Format, only: GetXCommand
    implicit none
    class(SimpleBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    lineStat = 0

  end subroutine
!==========================================================================================
  subroutine Simplebox_LoadAtomCoord(self, line, lineStat)
!    use Box_Utility, only: FindMolecule
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
    integer :: molType, molIndx, atmIndx
    integer :: arrayIndx, subIndx
    integer :: iType, iMol, iAtom
    real(dp) :: x,y,z

    if( .not. allocated(self%atoms) ) then
      call self%Constructor
    endif

    read(line, *) molType, molIndx, atmIndx, x, y ,z

    if((molType > nMolTypes) .or. (molType < 1)) then
      write(*,*) "ERROR! Type Index out of bounds!"
      write(*,*) molType, molIndx, atmIndx
      lineStat = -1
      return
    endif

    if( (molIndx > self%NMolMax(molType)) .or. (molIndx < 1) ) then
      write(*,*) "ERROR! Molecule Index out of bounds!"
      write(*,*) molType, molIndx, atmIndx
      lineStat = -1
      return
    endif

    subIndx = 0
    do iType = 1, molType-1
      subIndx = self%NMolMax(iType)
    enddo
    subIndx = subIndx + molIndx
    arrayIndx = self%MolStartIndx(subIndx)
    arrayIndx = arrayIndx + atmIndx - 1

    self%atoms(1, arrayIndx) = x
    self%atoms(2, arrayIndx) = y
    self%atoms(3, arrayIndx) = z

  end subroutine
!==========================================================================================
  subroutine SimpleBox_BuildNeighList(self)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: iList
    integer :: iAtom, jAtom
    real(dp) :: rx, ry, rz, rsq

  end subroutine
!==========================================================================================
  subroutine SimpleBox_Boundary(self, rx, ry, rz)
    implicit none
    class(SimpleBox), intent(in) :: self
    real(dp), intent(inout) :: rx, ry, rz 

  end subroutine
!==========================================================================================
subroutine SimpleBox_ComputeEnergy(self)
  implicit none
  class(SimpleBox), intent(inout) :: self
  logical :: accept

  call self % EFunc % Method % DetailedECalc( self, self%ETotal, accept )
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
  subroutine SimpleBox_UpdateNeighLists(self, disp)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
    integer :: iDisp, iList
    integer :: atmIndx, iAtom, jAtom
    real(dp) :: rx, ry, rz, rsq

  end subroutine
!==========================================================================================
  function SimpleBox_CheckConstraint(self, disp) result(accept)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    type(Displacement), intent(in) :: disp(:)
    logical :: accept
    integer :: nDisp, iConstrain

    accept = .true.
    if( .not. allocated(self%Constrain) ) then
      return
    endif

    nDisp = size(disp)
    if( size(self%Constrain) > 0 ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % DiffCheck( self, disp(1:nDisp), accept )
      enddo
      if(.not. accept) then
        return
      endif
    endif     

  end function
!==========================================================================================
  subroutine SimpleBox_IOProcess(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    implicit none

    class(SimpleBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    logical :: logicVal
    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("buildfreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % buildfreq = intVal
        self % maintFreq = intVal

      case("chempot")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self % chempot(intVal) = realVal

      case("energycalc")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % EFunc => EnergyCalculator(intVal)

      case("neighcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self%NeighList(intVal)%rCut = realVal

      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call self%NeighList(intVal)%ProcessIO(line, lineStat)

      case("pressure")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % pressure = realVal

      case("temperature")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % temperature = realVal
        self % beta = 1E0_dp/realVal


!      case("recenter")
!        call GetXCommand(line, command, 5, lineStat)
!        read(command, *) logicVal
!        self % zeroCoords = logicVal

      case default
        lineStat = -1
    end select
  end subroutine
!==========================================================================================
  subroutine SimpleBox_DumpData(self, filename)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: iType, iMol, iAtom, jType, subIndx, arrayIndx

    open(unit=50, file=trim(adjustl(filename)))

    write(50,*) "boxtype nobox"
    write(50,*) 
    write(50,*) "molmin", (self%NMolMin(iType), iType=1,nMolTypes)
    write(50,*) "molmax", (self%NMolMax(iType), iType=1,nMolTypes)
    write(50,*) "mol", (self%NMol(iType), iType=1,nMolTypes)
    write(50,*) "# MolType,  MolNumber, AtomNumber, x1, x2......."

    do iType = 1, nMolTypes
      do iMol = 1, self%NMol(iType)
        do iAtom = 1, MolData(iType)%nAtoms
          subIndx = 0
          do jType = 1, iType-1
            subIndx = self%NMolMax(jType)
          enddo
          subIndx = subIndx + iMol
          arrayIndx = self%MolStartIndx(subIndx)
          arrayIndx = arrayIndx + iAtom - 1

          write(50,*) iType, iMol, iAtom, self%atoms(1,arrayIndx), &
                       self%atoms(2,arrayIndx), self%atoms(3,arrayIndx)
        enddo
      enddo
    enddo


    close(50)

  end subroutine
!==========================================================================================
  subroutine SimpleBox_AddMol(self, molType)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in) :: molType

    self % NMol(molType) = self % NMol(molType) + 1
    self % nAtoms = self % nAtoms + MolData(molType)%nAtoms
  end subroutine
!==========================================================================================
  subroutine SimpleBox_DeleteMol(self, molIndx)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in) :: molIndx
    integer :: iList, iDimn, iAtom
    integer :: lastMol, iType, nType
    integer :: nStart, jStart

    nStart = self % MolStartIndx(molIndx)
    nType = self % MolType(nStart)

    lastMol = 0
    do iType = 1, nType-1
      lastMol = lastMol + self%NMolMax(iType) 
    enddo
    lastMol = lastMol + self%NMol(nType)
!    if(molIndx == lastMol) then
!      self % NMol(nType) = self % NMol(nType) - 1 
!      self % nAtoms = self % nAtoms - MolData(nType)%nAtoms
!      return
!    endif
    jStart = self%MolStartIndx(lastMol)

!     Take the top molecule from the atom array and move it's position in the deleted
!     molecule's slot.
    do iAtom = 1, MolData(nType) % nAtoms
      do iDimn = 1, self%nDimension
        self % atoms(iDimn, nStart+iAtom-1 ) = self % atoms(iDimn, jStart+iAtom-1 )
      enddo
    enddo

    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % DeleteMol(molIndx, lastMol)
    enddo
    self % NMol(nType) = self % NMol(nType) - 1 
    self % nAtoms = self % nAtoms - MolData(nType)%nAtoms

  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdatePosition(self, disp, tempList, tempNNei)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    type(Displacement), intent(inout) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iDisp, dispLen, dispIndx

    dispLen = size(disp)

    do iDisp = 1, dispLen
      if( disp(iDisp)%newAtom ) then 
        dispIndx = disp(iDisp) % atmIndx
        call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
        self % atoms(1, dispIndx) = disp(iDisp)%x_new
        self % atoms(2, dispIndx) = disp(iDisp)%y_new
        self % atoms(3, dispIndx) = disp(iDisp)%z_new
      endif
    enddo

    if(disp(iDisp)%newlist) then
      call self % NeighList(1) % AddMol(disp, tempList, tempNNei)
    endif

  end subroutine

!==========================================================================================
  subroutine SimpleBox_Maintenence(self)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self

    call self % NeighList(1) % BuildList


  end subroutine
!==========================================================================================
  subroutine SimpleBox_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: self

    call self % ComputeEnergy
    write(nout,*) "NeighList"
    call self % NeighList(1) % BuildList
    write(nout,*) "NeighList"

    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Initial Energy: ", self % ETotal
    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol


  end subroutine
!==========================================================================================
  subroutine SimpleBox_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp) :: E_Culm

    E_Culm = self%ETotal

    write(nout,*) "--------Box", self%boxID , "Energy---------"
    call self % ComputeEnergy
    call self % NeighList(1) % BuildList

    write(nout, *) "Final Energy:", self % ETotal
    if(self%ETotal /= 0) then
      if( abs((E_Culm-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm
        write(nout, *) "Final Energy: ", self%ETotal
        write(nout, *) "Difference: ", self%ETotal-E_Culm
      endif
    else
      if( abs(E_Culm-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm
        write(nout, *) "Final Energy: ", self%ETotal
        write(nout, *) "Difference: ", self%ETotal-E_Culm
      endif
    endif


  end subroutine
!==========================================================================================
  pure subroutine SimpleBox_CopyBox(box1, box2) 
    implicit none
    type(SimpleBox), intent(out) :: box1
    type(SimpleBox), intent(in) :: box2
    integer :: iAtom, iDim

    do iAtom = 1,box1%nMaxAtoms
      do iDim = 1, box1%nDimension
        box1%atoms(iDim, iatom) = box2%atoms(iDim, iatom)
      enddo
      box1%ETable(iAtom) = box2%ETable(iAtom)
    enddo
    box1%ETotal = box2%ETotal
    box1%Volume = box2%Volume

  end subroutine
!==========================================================================================
end module
!==========================================================================================
