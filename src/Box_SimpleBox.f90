!========================================================================================
#define __StdErr__ 0
!========================================================================================
module SimpleSimBox
  use VarPrecision
  use ForcefieldData, only: ECalcArray
  use ConstraintTemplate, only: constrainArray
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes


  !Sim Box Definition
  type, public, extends(SimBox) :: SimpleBox
!    ----------------------------
!   Inherited Variables from SimBox
!    character(len=15) :: boxStr = "Empty"
!    integer :: boxID
!    integer :: nAtoms, nMaxAtoms
!    integer :: nMolTotal
!    integer :: nDimension = 3
!
!    !Thermodynamic Variables
!    real(dp) :: ETotal, HTotal
!    real(dp) :: pressure = 0E0_dp
!    real(dp) :: beta, temperature, volume
!    real(dp), allocatable :: chempot(:)
!
!    real(dp), allocatable :: ETable(:), dETable(:)
!    real(dp), allocatable :: atoms(:,:)
!
!    Molecule Based Indexing and Census Arrays 
!    integer, allocatable :: NMolMin(:), NMolMax(:)
!    integer, allocatable :: NMol(:), MolStartIndx(:), MolEndIndx(:)
!
!    Atom Based Indexing Arrays 
!    integer, allocatable :: AtomType(:), MolType(:)
!    integer, allocatable :: MolIndx(:), MolSubIndx(:), AtomSubIndx(:)
!
!      
!    integer, allocatable :: MolGlobalIndx(:, :)
!    integer, allocatable :: TypeFirst(:), TypeLast(:)
!
!    integer :: nLists
!    ----------------------------
!    
    integer :: nTotal
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
      procedure, pass :: EnergySafetyCheck => SimpleBox_EnergySafetyCheck
      procedure, pass :: ComputeCM => SimpleBox_ComputeCM


      !IO Functions
      procedure, pass :: ProcessIO => SimpleBox_ProcessIO
      procedure, pass :: ScreenOut => SimpleBox_ScreenOut
!      procedure, pass :: ProcessIOCommon => SimpleBox_ProcessIOCommon

      procedure, pass :: CheckConstraint => SimpleBox_CheckConstraint
      procedure, pass :: CheckPostEnergy => SimpleBox_CheckPostEnergy
      procedure, pass :: DumpData => SimpleBox_DumpData

      !Coordinate Processing Functions
      procedure, pass :: GetNeighborList => Simplebox_GetNeighborList
      procedure, pass :: GetNewNeighborList => Simplebox_GetNewNeighborList
      procedure, pass :: GetDimensions => Simplebox_GetDimensions
      procedure, pass :: GetIndexData => Simplebox_GetIndexData
      procedure, pass :: GetMolData => SimpleBox_GetMolData
      procedure, pass :: GetMolEnergy => SimpleBox_GetMolEnergy
      procedure, pass :: GetMaxAtoms => SimpleBox_GetMaxAtoms
      procedure, pass :: CountAtoms => SimpleBox_CountAtoms
!      procedure, pass :: GetCoordinates => SimpleBox_GetCoordinates

      !New Property Gathering Functions
      procedure, pass :: GetNewEnergy => Simplebox_GetNewEnergy
      procedure, pass :: GetNewMolCount => Simplebox_GetNewMolCount
!      procedure, pass :: GetNewAtomCount => Simplebox_GetNewAtomCount



      !Update Functions
      procedure, pass :: IsActive => SimpleBox_IsActive
      procedure, pass :: AddMol => SimpleBox_AddMol
      procedure, pass :: DeleteMol => SimpleBox_DeleteMol
      procedure, pass :: Update => SimpleBox_Update
      procedure, pass :: UpdateVolume => SimpleBox_UpdateVolume
      procedure, pass :: UpdateEnergy => SimpleBox_UpdateEnergy
      procedure, pass :: UpdatePosition => SimpleBox_UpdatePosition
      procedure, pass :: UpdateNeighLists => SimpleBox_UpdateNeighLists

      procedure, pass :: GetReducedCoords => SimpleBox_GetReducedCoords
      procedure, pass :: GetRealCoords => SimpleBox_GetRealCoords

!      procedure, public, pass :: GetThermo
!      procedure, public, pass :: ThermoLookUp

      procedure, pass :: Maintenance => SimpleBox_Maintenance
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
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: AllocateStatus
    integer :: iType, iMol, iAtom, atmIndx, molIndx, maxMol, maxSingleMol
 
    if( .not. allocated(self%NMolMin) ) then
      write(*,*) "ERROR! The maximum and minimum molecules allowed in the box must be defined"
      write(*,*) "prior to box initialization!"
      stop 
    endif
    AllocateStatus = 0
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
    allocate(self%AtomSubIndx(1:self%nMaxAtoms), stat=AllocateStatus)

    allocate(self%MolStartIndx(1:maxMol), stat=AllocateStatus)
    allocate(self%MolEndIndx(1:maxMol), stat=AllocateStatus)
    allocate(self%centerMass(1:3, 1:maxMol), stat=AllocateStatus)

    allocate(self%TypeFirst(1:nMolTypes), stat=AllocateStatus)
    allocate(self%TypeLast(1:nMolTypes), stat=AllocateStatus)

    maxSingleMol = maxval(self%NMolMax(:))
    allocate(self%MolGlobalIndx(1:nMolTypes, 1:maxSingleMol), stat=AllocateStatus)

    allocate(self%chempot(1:nMolTypes), stat=AllocateStatus)
    self%chempot = 0E0_dp
    IF (AllocateStatus /= 0) STOP "Allocation Error in Simulation Box Def"

    self%maxMol = maxMol
    self%AtomType = 0
    self%MolType = 0
    self%MolIndx = 0
    self%MolSubIndx = 0
    self%AtomSubIndx = 0
    self%MolStartIndx = 0
    self%MolEndIndx = 0

    self%TypeFirst = 0
    self%TypeLast = 0

    self%MolGlobalIndx = 0

    self%chempot = 0E0_dp


    !This block creates the indexing arrays that can be used to find features
    !such as what molecule an atom belongs to, indexs relative to an atoms position in the
    !molecule, the molecules type, atom type, the indicies of other atoms in the molecule.
    !
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
          self%AtomSubIndx(atmIndx)  = iAtom
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
    IF (AllocateStatus /= 0) STOP "*** SimpleBox: Unable to allocate Mol Bounds ***"
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
  subroutine Simplebox_ComputeCM(self, molIndx)
    use Common_MolInfo, only: AtomData 
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in) :: molIndx
    integer :: iAtom, atmType
    integer :: molStart, molEnd
    real(dp) :: xcm, ycm, zcm, totalMass
    real(dp) :: x1, y1, z1
    real(dp) :: xn, yn, zn

    call self%GetMolData(molIndx, molStart=molStart, molEnd=molEnd)
    totalMass = 0E0_dp
    xcm = 0E0_dp
    ycm = 0E0_dp
    zcm = 0E0_dp

    x1 = self%atoms(1, molStart)
    y1 = self%atoms(2, molStart)
    z1 = self%atoms(3, molStart)

    do iAtom = molStart, molEnd
      atmType = self % AtomType(iAtom) 
      xn = self%atoms(1, iAtom) - x1
      yn = self%atoms(2, iAtom) - y1
      zn = self%atoms(3, iAtom) - z1
      call self%Boundary(xn, yn, zn)
      xcm = xcm + AtomData(atmType)%mass * xn
      ycm = ycm + AtomData(atmType)%mass * yn
      zcm = zcm + AtomData(atmType)%mass * zn
      totalMass = totalMass + AtomData(atmType)%mass
    enddo

    do iAtom = molStart, molEnd
      xcm = xcm/totalMass + x1
      ycm = ycm/totalMass + y1
      zcm = zcm/totalMass + z1
    enddo
    call self%Boundary(xcm, ycm, zcm)

    self % centerMass(1, molIndx) = xcm
    self % centerMass(2, molIndx) = ycm
    self % centerMass(3, molIndx) = zcm

  end subroutine
!==========================================================================================
  subroutine Simplebox_GetDimensions(self, list)
    use Input_Format, only: GetXCommand
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(inout) :: list(:, :)

    list = 0E0_dp

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
      call self % Constructor
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
! This subroutine recomputes the entire system energy from scratch. If
! tablecheck is passed the code will compare the current energy table
! to a recalculated version to see if it has bene properly kept.
  subroutine SimpleBox_ComputeEnergy(self, tablecheck)
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical, optional, intent(in) :: tablecheck
    logical :: accept
    integer :: iAtom
    real(dp), allocatable :: tempETable(:)
    real(dp) :: ECul, ECalc, EDiff

    if(present(tablecheck))then
      if(tablecheck) then
        allocate(tempETable(1:size(self%ETable)))
      endif
      tempETable = self%ETable
    endif

    call self % EFunc % Method % DetailedECalc( self, self%ETotal, accept )

    if(present(tablecheck))then
      do iAtom = 1, self%nMaxAtoms
        if(.not. self%IsActive(iAtom)) cycle
        ECul = tempETable(iAtom)
        ECalc = self%ETable(iAtom)
        EDiff = abs(ECul-ECalc)
        if(ECalc /= 0E0_dp) then
          if( EDiff/abs(ECalc)  > 1E-7_dp ) then
            write(nout, *) "Table Discrepancy: ", iAtom, ECul/outEngUnit, ECalc/outEngUnit,  EDiff/outEngUnit, engStr
          endif
        else
          if( EDiff > 1E-7_dp ) then
            write(nout, *) "Table Discrepancy: ", iAtom, ECul/outEngUnit, ECalc/outEngUnit,  EDiff/outEngUnit, engStr
          endif
        endif
      enddo
      if(tablecheck) then
        deallocate(tempETable)
      endif
      write(nout,*) "Energy Table Check complete!"
    endif
  end subroutine
!==========================================================================================
  subroutine SimpleBox_EnergySafetyCheck(self)
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    real(dp) :: E_Current


    E_Current = self%ETotal
    call self % EFunc % Method % DetailedECalc( self, self%ETotal, accept )

    if(self%ETotal /= 0E0_dp) then
      if( abs((E_Current-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "!!!!!!!!!!!!!!!!!!ERROR! Large energy drift detected!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Current/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Current)/outEngUnit, engStr
      endif
    else
      if( abs(E_Current-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "!!!!!!!!!!!!!!!!!!ERROR! Large energy drift detected!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Current/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Current)/outEngUnit, engStr
      endif
    endif



  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdateEnergy(self, E_Diff)
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: E_Diff

    self % ETotal = self % ETotal + E_Diff
    self % ETable = self % ETable + self % dETable
    self%dETable = 0E0_dp

  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdateVolume(self, disp)
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdateNeighLists(self, disp)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    integer :: iDisp, iList
    integer :: atmIndx, iAtom, jAtom
    real(dp) :: rx, ry, rz, rsq

  end subroutine
!==========================================================================================
  function SimpleBox_CheckConstraint(self, disp) result(accept)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    logical :: accept
    integer :: nDisp, iConstrain

    accept = .true.
    if( .not. allocated(self%Constrain) ) then
      return
    endif

    nDisp = size(disp)
    if( allocated(self%Constrain) ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % DiffCheck( self, disp(1:nDisp), accept )
      enddo
      if(.not. accept) then
        return
      endif
    endif

  end function
!==========================================================================================
  function SimpleBox_CheckPostEnergy(self, disp, E_Diff) result(accept)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff
    logical :: accept
    integer :: nDisp, iConstrain

    accept = .true.
    if( .not. allocated(self%Constrain) ) then
      return
    endif


    nDisp = size(disp)
    if( allocated(self%Constrain) ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % PostEnergy( self, disp(1:nDisp), E_Diff, accept )
      enddo
      if(.not. accept) then
        return
      endif
    endif

  end function

!==========================================================================================
  subroutine SimpleBox_ProcessIO(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    use Units, only: inPressUnit
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
        self % maintFreq = intVal

      case("chempotential")
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
        self % pressure = realVal*inPressUnit

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

          write(50,*) iType, iMol, iAtom, self%atoms(1,arrayIndx), self%atoms(2,arrayIndx), self%atoms(3,arrayIndx)
        enddo
      enddo
    enddo


    close(50)

  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetIndexData(self, MolIndx, MolSubIndx, AtomSubIndx)
    implicit none
    class(SimpleBox), intent(inout), target :: self
    integer, intent(inout), pointer, optional :: MolIndx(:), MolSubIndx(:), AtomSubIndx(:)

    if(present(MolIndx)) then
      MolIndx => self%MolIndx
    endif

    if(present(MolSubIndx)) then
      MolSubIndx => self%MolSubIndx
    endif

    if(present(AtomSubIndx)) then
      AtomSubIndx => self%AtomSubIndx
    endif



  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetMolData(self, globalIndx, molStart, molEnd, molType, atomSubIndx)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in)  :: globalIndx
    integer, intent(inout), optional :: molStart, molEnd, molType, atomSubIndx


    if((size(self%molStartIndx) < globalIndx) .or. (globalIndx < 1)) then
      write(__StdErr__, *) "Error in Sim Box class: GetMolData has been given an invalid index"
      write(__StdErr__, *) "Given Index:", globalIndx 
      write(__StdErr__, *) "Array Size:", size(self%molStartIndx)
    endif

    if( present(molStart) ) then
      molStart = self % MolStartIndx(globalIndx)
    endif

    if( present(molEnd) ) then
      molEnd = self % MolEndIndx(globalIndx)
    endif

    if( present(molType) ) then
      molType = self % MolType(globalIndx)
    endif

    if( present(atomSubIndx) ) then
      atomSubIndx = self%AtomSubIndx(globalIndx)
    endif

  end subroutine

!==========================================================================================
  function SimpleBox_GetMaxAtoms(self) result(nMax)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: nMax
    integer :: maxMol, iType
    

    self%nMaxAtoms = 0
    maxMol = 0
    do iType = 1, nMolTypes
      self%nMaxAtoms = self%nMaxAtoms + self%NMolMax(iType)*MolData(iType)%nAtoms
      maxMol = maxMol + self%NMolMax(iType)
    enddo
    nMax = self%nMaxAtoms


  end function
!==========================================================================================
! Returns the total number of atoms in the system.  Alternatively can be used to
! get the atoms of a single type by passing a value to atmtype
  function SimpleBox_CountAtoms(self, atmtype) result(nCount)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, optional, intent(in) :: atmtype
    integer :: iAtom, atmType1
    integer :: nCount

    nCount = 0

    if(present(atmtype)) then
       do iAtom = 1, self%nMaxAtoms
         if(.not. self%IsActive(iAtom)) then
           cycle
         endif
         atmType1 = self%AtomType(iAtom)
         if(atmType1 == atmType) then
           nCount = nCount + 1
         endif
      enddo
    else
       do iAtom = 1, self%nMaxAtoms
         if(self%IsActive(iAtom)) then
           nCount = nCount + 1
         endif
      enddo
    endif


  end function
!==========================================================================================
  subroutine SimpleBox_GetNeighborList(self, listindx, neighlist, nNei)
    implicit none
    class(SimpleBox), intent(inout), target :: self
    integer, intent(in) :: listIndx
    integer, pointer, intent(inout) :: neighlist(:,:), nNei(:)

    call self%NeighList(listindx)%GetListArray(neighlist, nNei)

  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetNewNeighborList(self, listindx, iAtom, tempList, tempNNei, disp)
    implicit none
    class(SimpleBox), intent(inout), target :: self
    integer, intent(in) :: listIndx
    integer, intent(in) :: iAtom
    integer, intent(inout) :: tempNnei(:)
    integer, intent(inout) :: tempList(:, :)
    class(Perturbation), intent(inout) :: disp

    call self%NeighList(listindx)%GetNewList(iAtom, tempList, tempNNei, disp)
  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetEnergyTable(self, etable, detable)
    implicit none
    class(SimpleBox), intent(inout), target :: self
    real(dp), pointer, optional :: ETable(:)
    real(dp), pointer, optional :: dETable(:)

    if(present(etable)) then
      ETable => self%ETable
    endif

    if(present(dETable)) then
      dETable => self%dETable
    endif



  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetMolEnergy(self, molIndx, EMol, newstate)
    implicit none
    class(SimpleBox), intent(inout), target :: self
    integer, intent(in) :: molIndx
    logical, intent(in), optional :: newstate
    real(dp), intent(out) :: EMol
    logical :: calcDE
    integer :: molStart, molEnd
    integer :: iAtom

    if(present(newstate)) then
      calcDE = newstate
    else
      calcDE =.false.
    endif

    call self%GetMolData(molIndx, molStart=molStart, molEnd=molEnd)

    EMol = 0E0_dp
    do iAtom = molStart, molEnd
      if(.not. calcDE) then
        EMol = EMol + self%ETable(iAtom)
      else
        EMol = EMol + self%ETable(iAtom) + self%dETable(iAtom)
      endif
    enddo


  end subroutine
!==========================================================================================
  function SimpleBox_GetNewEnergy(self, E_Diff) result(E_New)
    implicit none
    class(SimpleBox), intent(in) :: self
    real(dp), intent(in)  :: E_Diff
    real(dp)  :: E_New

    E_New = self % ETotal + E_Diff

  end function
!==========================================================================================
  function SimpleBox_GetNewMolCount(self, disp) result(nNew)
    implicit none
    class(SimpleBox), intent(in) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer :: nNew, molType

    nNew = self%nMolTotal

    select type(disp)
      class is(Addition)
        nNew = nNew + 1
      class is(Deletion)
        nNew = nNew - 1
    end select

  end function
!==========================================================================================
!  Checks to see if the atom is present in the box.  This is required 
  function SimpleBox_IsActive(self, atmIndx) result(active)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: active
    integer, intent(in) :: atmIndx

    if( self%MolSubIndx(atmIndx) > self%NMol(self%MolType(atmIndx)) ) then
      active = .false.
    else
      active = .true.
    endif

  end function
!==========================================================================================
  subroutine SimpleBox_AddMol(self, molType)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in) :: molType

    self % NMol(molType) = self % NMol(molType) + 1
    self % nAtoms = self % nAtoms + MolData(molType)%nAtoms
    self % nMolTotal = self % nMolTotal + 1
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
    self%centermass(1:3,molIndx) = self%centermass(1:3, lastMol)
!    if(molIndx == lastMol) then
!      self % NMol(nType) = self % NMol(nType) - 1 
!      self % nAtoms = self % nAtoms - MolData(nType)%nAtoms
!      return
!    endif
!    write(*,*) lastMol, nStart, nType, molIndx
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
    self % nMolTotal = self % nMolTotal - 1
    self % nAtoms = self % nAtoms - MolData(nType)%nAtoms

  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdatePosition(self, disp, tempList, tempNNei)
    use CoordinateTypes
    use Common_MolInfo, only: MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iDisp, dispLen, dispIndx
    integer :: iAtom, iMol, molStart, molEnd
    real(dp) :: dx, dy, dz

    select type(disp)
       !-------------------------------------------------
      class is(Displacement)
        dispLen = size(disp)
        do iDisp = 1, dispLen
          dispIndx = disp(iDisp) % atmIndx
          call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
          self % atoms(1, dispIndx) = disp(iDisp)%x_new
          self % atoms(2, dispIndx) = disp(iDisp)%y_new
          self % atoms(3, dispIndx) = disp(iDisp)%z_new
        enddo
        call self%ComputeCM(disp(1)%molIndx)


       !-------------------------------------------------
      class is(Addition)
        dispLen = size(disp)
        do iDisp = 1, dispLen
          dispIndx = disp(iDisp) % atmIndx
          call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
          self % atoms(1, dispIndx) = disp(iDisp)%x_new
          self % atoms(2, dispIndx) = disp(iDisp)%y_new
          self % atoms(3, dispIndx) = disp(iDisp)%z_new
        enddo
        call self % NeighList(1) % AddMol(disp, tempList, tempNNei)
        call self % AddMol(disp(1)%molType)
        call self % ComputeCM(disp(1)%molIndx)

       !-------------------------------------------------
      class is(OrthoVolChange)
        do iMol = 1, self%maxMol
          call self%GetMolData(iMol, molStart=molStart, molEnd=molEnd)
          if( self%MolSubIndx(molStart) <= self%NMol(self%MolType(molStart)) ) then
            dx = (disp(1)%xScale-1E0_dp) * self%centerMass(1, iMol)
            dy = (disp(1)%yScale-1E0_dp) * self%centerMass(2, iMol)
            dz = (disp(1)%zScale-1E0_dp) * self%centerMass(3, iMol)
            do iAtom = molStart, molEnd
!              write(*,*) self%atoms(1:3, iAtom)
              self%atoms(1, iAtom) = self%atoms(1, iAtom) + dx
              self%atoms(2, iAtom) = self%atoms(2, iAtom) + dy
              self%atoms(3, iAtom) = self%atoms(3, iAtom) + dz
!              write(*,*) self%atoms(1:3, iAtom)
            enddo
            self%centerMass(1, iMol) = self%centerMass(1, iMol) * disp(1)%xScale
            self%centerMass(2, iMol) = self%centerMass(2, iMol) * disp(1)%yScale
            self%centerMass(3, iMol) = self%centerMass(3, iMol) * disp(1)%zScale
          endif
        enddo
        call self%UpdateVolume(disp)

       !-------------------------------------------------
      class default
        stop "The code does not know how to update coordinates for this perturbation type."
    end select


  end subroutine
!==========================================================================================
  subroutine SimpleBox_Maintenance(self)
    use CoordinateTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    integer :: iList
    real(dp) :: tempE

    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % BuildList(iList)
    enddo
!    call self % ComputeEnergy

!    call self % EFunc % Method % DetailedECalc( self, tempE, accept )

!    if(abs(self%ETotal-tempE)/abs(tempE) > 1e-4) then
!   write(2,*) self%ETotal, tempE, self%ETotal-tempE
!    endif

  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetReducedCoords(self,realCoords,reducedCoords )
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: realCoords(:)
    real(dp), intent(out) :: reducedCoords(1:3)

    reducedCoords = 0E0_dp
    stop "GetReducedCoords routine has been called on a system where no box is defined."

  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetRealCoords(self, reducedCoords, realCoords)
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: reducedCoords(:)
    real(dp), intent(out) :: realCoords(1:3)

    realCoords = 0E0_dp
    stop "GetReducedCoords routine has been called on a system where no box is defined."
  end subroutine
!==========================================================================================
  subroutine SimpleBox_ScreenOut(self)
    use Common_MolInfo, only: nMolTypes
    use ParallelVar, only: nout
    use Input_Format, only: ReplaceText
    use Units, only: outEngUnit
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    integer :: iConstrain, iMol
    integer :: iType, molStart
    character(len=200) :: tempStr
    character(len=80) :: tempStr2

    write(tempStr, "(A)") "      ************************ Box %s ************************ "
    write(tempStr2, "(I50)") self%boxID
    tempStr = ReplaceText(tempStr, "%s", trim(adjustl(tempStr2)))
    write(nout, "(A)") trim(tempStr)

    write(tempStr, "(A)") "        Total Energy: %s1   Number of Molecules: %s2"
    write(tempStr2, "(F40.8)") self%ETotal/outEngUnit
    tempStr = ReplaceText(tempStr, "%s1", trim(adjustl(tempStr2)))
    write(tempStr2, "(I40)") self%nMolTotal
    tempStr = ReplaceText(tempStr, "%s2", trim(adjustl(tempStr2)))
    write(nout, "(A)") trim(tempStr)

    write(tempStr, "(A)") "        Per Mol Energy: %s1"
    write(tempStr2, "(F40.8)") self%ETotal/(outEngUnit*self%nMolTotal)
    tempStr = ReplaceText(tempStr, "%s1", trim(adjustl(tempStr2)))
    write(nout, "(A)") trim(tempStr)





  end subroutine
!==========================================================================================
  subroutine SimpleBox_Prologue(self)
    use Common_MolInfo, only: nMolTypes
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    integer :: iConstrain, iMol, iList
    integer :: iType, molStart

    self%nMolTotal = 0
    do iType = 1, nMolTypes    
      self%nMolTotal = self%nMolTotal + self % NMol(iType)
    enddo
    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout,*) "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal
    write(nout,*) "Box ", self%boxID, " Temperature: ", self % temperature

    call self % ComputeEnergy
    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % BuildList(iList)
    enddo

    if( size(self%Constrain) > 0 ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % Prologue
        call self%Constrain(iConstrain) % method % CheckInitialConstraint(self, accept)
      enddo
      if(.not. accept) then
        write(nout,*) "Initial Constraints are not statisfied!"
        stop
      endif
    endif


    write(nout, "(1x,A,I2,A,E15.8,1x,A)") "Box ", self%boxID, " Initial Energy: ", self % ETotal/outEngUnit, engStr
    write(nout, "(1x,A,I2,A,E15.8,1x,A)") "Box ", self%boxID, " Initial Energy (Per Mol): ", &
                                           self % ETotal/(outEngUnit*self%nMolTotal), engStr

    do iMol = 1, self%maxMol
      call self%GetMolData(iMol, molStart=molStart)
      if( self%MolSubIndx(molStart) <= self%NMol(self%MolType(molStart)) ) then
        call self%ComputeCM(iMol)
      endif
    enddo

  end subroutine
!==========================================================================================
  subroutine SimpleBox_Epilogue(self)
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: iConstrain, iList
    real(dp) :: E_Culm

    E_Culm = self%ETotal

    write(nout,*) "--------Box", self%boxID , "Energy---------"
    call self % ComputeEnergy(tablecheck=.true.)
    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % BuildList(iList)
    enddo

    write(nout, *) "Final Energy:", self % ETotal/outEngUnit, engStr
    if(self%ETotal /= 0) then
      if( abs((E_Culm-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit, engStr
      endif
    else
      if( abs(E_Culm-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit, engStr
      endif
    endif
    write(nout, "(1x,A,I2,A,E15.8,1x,A)") "Box ", self%boxID, " Final Energy (Per Mol): ", &
                                           self % ETotal/(outEngUnit*self%nMolTotal), engStr
    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout,*) "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal
    if( allocated(self%Constrain) ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % Epilogue
      enddo
    endif
  end subroutine
!==========================================================================================
  subroutine SimpleBox_Update(self)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer :: iConstrain, iList

    if( allocated(self%Constrain) ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % Update
      enddo
    endif

    if( allocated(self%NeighList) ) then
      do iList = 1, size(self%NeighList)
        call self%NeighList(iList) % Update
      enddo
    else
      stop "No Neighbor List has been defined!"
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
