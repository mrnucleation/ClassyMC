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
      procedure, pass :: AllocateMolBound => SimpleBox_AllocateMolBound
      procedure, pass :: LoadAtomCoord => Simplebox_LoadAtomCoord
      procedure, pass :: LoadDimension => Simplebox_LoadDimension
      procedure, pass :: BuildNeighList => SimpleBox_BuildNeighList
      procedure, pass :: Boundary => SimpleBox_Boundary
      procedure, pass :: UpdateEnergy => SimpleBox_UpdateEnergy
      procedure, pass :: UpdatePosition => SimpleBox_UpdatePosition
      procedure, pass :: UpdateNeighLists => SimpleBox_UpdateNeighLists
      procedure, pass :: ComputeEnergy => SimpleBox_ComputeEnergy
      procedure, pass :: IOProcess => SimpleBox_IOProcess
      procedure, pass :: CheckConstraint => SimpleBox_CheckConstraint
      procedure, pass :: DumpData => SimpleBox_DumpData
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

    self%boxStr = "NoBox"
    !First begin by computing the maximium number of atoms that the box can potentially contain
    self%nMaxAtoms = 0
    maxMol = 0
    do iType = 1, nMolTypes
      self%nMaxAtoms = self%nMaxAtoms + self%NMolMax(iType)*MolData(iType)%nAtoms
      maxMol = maxMol + self%NMolMax(iType)
    enddo
    write(*,*) "maxatoms", self%nMaxAtoms

    self%nAtoms = 0
    do iType = 1, nMolTypes
      self%nAtoms = self%nAtoms + self%NMol(iType)*MolData(iType)%nAtoms
    enddo
    write(*,*) "atoms", self%nAtoms


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


  end subroutine
!==========================================================================================
  subroutine Simplebox_LoadAtomCoord(self, line, lineStat)
!    use Box_Utility, only: FindMolecule
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
      write(*,*) self%nAtoms
    endif

    read(line, *) molType, molIndx, atmIndx, x, y ,z

    if( molIndx > self%NMolMax(molType) ) then
      write(*,*) "ERROR! Index out of bounds!"
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

!    write(*,*) molType, molIndx, atmIndx, arrayIndx
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

  call self % EFunc % Method % DetailedECalc( self, self%ETotal )
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
  function SimpleBox_CheckConstraint(self, disp) result(accept)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(in) :: self
    type(Displacement), intent(in) :: disp(:)
    logical :: accept
    integer :: nDisp, iConstrain

    accept = .true.
    if( .not. allocated(self%Constrain) ) then
      return
    endif

    if( size(self%Constrain) > 0 ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % ShiftCheck( self, disp(1:nDisp), accept )
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

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
!    write(*,*) command
    select case( trim(adjustl(command)) )
      case("energycalc")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % EFunc => EnergyCalculator(intVal)

      case("temperature")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % temperature = realVal
        self % beta = 1E0_dp/realVal

      case("neighcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self%NeighList(intVal)%rCut = realVal

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
end module
!==========================================================================================
