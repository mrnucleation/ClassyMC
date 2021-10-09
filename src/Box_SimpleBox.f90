!========================================================================================
! Wall-less simulation box which acts as the parent class to other simulaion box types.
! Unlike the base template this box can be used for either a boundless problem
! or in conjunction with a hard wall constrain to make a non-periodic condition.
!========================================================================================
#define __StdErr__ 0
module SimpleSimBox
  use VarPrecision
  use ForcefieldData, only: ECalcArray
  use ConstraintTemplate, only: constrainArray
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes


  !Sim Box Definition
  !Commented variables are inherieted from the parent class. Listed here for
  !convience.
  type, public, extends(SimBox) :: SimpleBox
!    ----------------------------
!   Inherited Variables from SimBox are shown as commented variables
!      integer :: boxID => This box's array index within the BoxArray
!      integer :: nAtoms => Number of current atoms within this box
!      integer :: nMaxAtoms => Largest number of atoms allowed within this box. 
!                              The box is full when nAtoms == nMaxAtoms
!      integer :: nMolTotal => Number of total molecules currently in the box
!      integer :: nDimension => Number of coordinates per atom. Default 
!                               is 3 for a 3-dimensional box with (x,y,z) 
!                               coordinates for each atom. 
!
!    character(len=15) :: boxStr = "Empty" => Box Type String used for output purposes.
!
!    !Thermodynamic Variables
!    real(dp) :: ETotal => Total Energy of the box, the sum of all intra and inter contributions
!    real(dp) :: HTotal => Total Enthalpy of the Box

!    real(dp) :: pressure => Pressure of the current simulation box. Constant for NPT 
!                            and computed for NVT
!    real(dp) :: beta => 1/kT value precomputed for convience
!    real(dp) :: temperature => Simulation temperature
!    real(dp) :: volume =>  Total volume of the box. 
!    real(dp), allocatable :: chempot(:) => Chemical potential of each molecule type 
!
!    real(dp), allocatable :: ETable(:) => Table of inter-molecular energies per molecule
!    real(dp), allocatable :: dETable(:) => Temporary Placeholder Array. Contains Change in 
!                                           ETable as a result of a given move.
!    real(dp), allocatable :: atoms(:,:) => Positional array of all the atoms including
!                                           placeholders for molecules that can be added.
!                                           First index is for the atom-subcoorindate (x,y,z..)
!                                           Second index is for the atom index
!                                           Atoms of a given molecule are sequential
!                                           Information for how to access a given
!                                           molecules positions are kept in the other arrays.
!                                                  
!                                           
!
!    logical :: forceoutofdate (default:.true.) => Used for when the system has changed
!                                                  to signal that the forces must be recomputed
!    real(dp) :: forcedelta => (Force delta used for finite difference method)
!    real(dp), allocatable :: forces(:,:) => Similar to the atoms(:,:) array, but keeps
!                                            the force gradient instead of the positional
!                                            components. 

    !forceERecompute => Forces the box to recompute the energy of the system.
    !                   Setting this to true is expensive and should only be used for debug
    !                   purposes when there's an error with the delta calculations of a system.
    logical :: forceERecompute = .false.

    !rebuilds => Book keeping array which tracks the number of times a neighbor list has
    !            been rebuilt
    !dangerbuilds => Keeps track of the number of times an atom moved a much larger distance
    !                than expected and as such an incorrect neighborlist might have been
    !                tabulated. It's bit buggy at the moment so this number may not be
    !                accurate.
    !dr => Tracks how much an atom has moved since the last neighborlist update.
    !      Used to determine if a neighborlist should be recomputed.
    !largestdr => Used for holding the largest displacement found in the dr array. 
    integer, private :: rebuilds = 0
    integer, private :: dangerbuilds = 0
    real(dp), allocatable, private :: dr(:,:)
    real(dp), private :: largestdr = 0E0_dp

    ! Temporary Molecule Position Storage Arrays. 
    real(dp), allocatable, private :: newpos(:,:)
    real(dp), allocatable, private :: temppos(:,:)

!
!    Molecule Based Indexing and Census Arrays 
!    These arrays are used to quickly look up location and informatics
!    in arrays such as the atoms(:,:), forces(:,:), dr(:,:), etc.
!    in order to get the spacial information required to perform
!    all manners of calculations.
! 
!     - Count and Bounds arrays for the number of
!     - (Indexing is based on number of molecule types)
!    NMolMin => Lower Bound for a molecule type. System is not allowed to go below this
!    NMolMax => Upper Bound for a molecule type. System is not allowed to go above this
!    NMol => Current molecule count current in the system of a given molecule type
!
!     - Quick look up for finding the location of a molecule in the atoms(:,:) arrays.
!     - (Indexing is based on global molecule index for this box)
!    MolStartIndx => atom(:,:) array index of the first atom of this molecule  
!    MolEndIndx => atom(:,:) array index of the last atom of this molecule 
!
!     -Reverse Look up array for molecules
!     -(First index is the molecule's type number and second is the molecule relative
!       position index)
!    MolGlobalIndx => Finds a molecule's global index from the Molecule Type (first index)
!                     and relative molecule index (second index)
!
!    integer, allocatable :: NMol(:), NMolMin(:), NMolMax(:)
!    integer, allocatable :: MolStartIndx(:), MolEndIndx(:)
!    integer, allocatable :: MolGlobalIndx(:, :)
!
!
!    - Look up for information about this Atom's type and what molecule it is part of.
!    - (Indexing is based on global atom index for this box)
!    MolType => Contains the type of molecule this atom belongs to. 
!    AtomType => Contains the type of an atom in the system.
!    MolIndx => Global Mol Index for the molecule this atom belongs to.
!    MolSubIndx => Relative Mol Index for the molecule this atom belongs to.
!    AtomSubIndx => Relative Atom Index of this atom.  This gives information about it's
!                   position with respect to it's parent molecule.  For example
!                   the first defined atom in the molecule will have Relative Index of 1,
!                   the second atom defined will have a relative index of 2, etc.
!
!    integer, allocatable :: AtomType(:), MolType(:)
!    integer, allocatable :: MolIndx(:), MolSubIndx(:), AtomSubIndx(:)
!
!
!     -Molecule Type look up
!     -(Indexing is based on molecule type ID)
!    TypeFirst => The atom global index of the first atom of the first molecule of this type
!    TypeLast => The atom global index of the kast atom of the last molecule of this type
!
!     nLists => Number of neighborlists, depreciated.
!
!     integer, allocatable :: TypeFirst(:), TypeLast(:)
!     integer :: nLists
!    ----------------------------
!    
    integer :: volnum
!    integer :: nTotal

    ! Constraint Class Array, used to reject configurations that
    ! deviate from the user specified criterias. 
    ! See Constraint object types for more info.
    ! This only contains constraints applied to this specific box.
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
      procedure, pass :: CheckLists => SimpleBox_CheckLists

      procedure, pass :: ComputeCM => SimpleBox_ComputeCM
      procedure, pass :: ComputeEnergy => SimpleBox_ComputeEnergy
      procedure, pass :: ComputeMolIntra => SimpleBox_ComputeMolIntra
      procedure, pass :: ComputeIntraEnergy => SimpleBox_ComputeIntraEnergy
      procedure, pass :: ComputeIntraEnergyDelta => SimpleBox_ComputeIntraEnergyDelta
      procedure, pass :: ComputeEnergyDelta => SimpleBox_ComputeEnergyDelta
      procedure, pass :: ComputeForces => SimpleBox_ComputeForces
      procedure, pass :: EnergySafetyCheck => SimpleBox_EnergySafetyCheck

      !IO Functions
      procedure, pass :: ProcessIO => SimpleBox_ProcessIO
      procedure, pass :: ScreenOut => SimpleBox_ScreenOut
!      procedure, pass :: ProcessIOCommon => SimpleBox_ProcessIOCommon

      procedure, pass :: CheckConstraint => SimpleBox_CheckConstraint
      procedure, pass :: CheckPostEnergy => SimpleBox_CheckPostEnergy
      procedure, pass :: DumpData => SimpleBox_DumpData

      !Coordinate Processing Functions
!      procedure, pass :: GetAtomTypes => Simplebox_GetAtomTypes
      procedure, pass :: GetNeighborList => Simplebox_GetNeighborList
      procedure, pass :: GetNewNeighborList => Simplebox_GetNewNeighborList
      procedure, pass :: GetLargestNNei => SimpleBox_GetLargestNNei
      procedure, pass :: GetDimensions => Simplebox_GetDimensions
      procedure, pass :: GetForceArray => SimpleBox_GetForceArray
      procedure, pass :: GetIndexData => Simplebox_GetIndexData
      procedure, pass :: GetNeighSkin => SimpleBox_GetNeighSkin
      procedure, pass :: GetAtomData => SimpleBox_GetAtomData
      procedure, pass :: GetMolData => SimpleBox_GetMolData
      procedure, pass :: GetMolEnergy => SimpleBox_GetMolEnergy
      procedure, pass :: GetMaxAtoms => SimpleBox_GetMaxAtoms
      procedure, pass :: CountAtoms => SimpleBox_CountAtoms
!      procedure, pass :: GetCoordinates => SimpleBox_GetCoordinates
      procedure, pass :: GetEFunc => SimpleBox_GetEFunc

      !New Property Gathering Functions
      procedure, pass :: GetNewEnergy => Simplebox_GetNewEnergy
      procedure, pass :: GetNewMolCount => Simplebox_GetNewMolCount
!      procedure, pass :: GetNewAtomCount => Simplebox_GetNewAtomCount

      !Search Functions
      procedure, pass :: FindMolByTypeIndex => SimpleBox_FindMolByTypeIndex

      !Update Functions
      procedure, pass :: IsActive => SimpleBox_IsActive
      procedure, pass :: AddMol => SimpleBox_AddMol
      procedure, pass :: DeleteMol => SimpleBox_DeleteMol
      procedure, pass :: SwapAtomType => SimpleBox_SwapAtomType


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

      procedure, pass :: Destructor => SimpleBox_Destructor

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
      write(0,*) "ERROR! The maximum and minimum molecules allowed in the box must be defined"
      write(0,*) "prior to box initialization!"
      error stop 
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
    allocate(self%dr(1:3, 1:self%nMaxAtoms), stat=AllocateStatus)
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
    IF (AllocateStatus /= 0) error STOP "Allocation Error in Simulation Box Def"

    self%volume = 0E0_dp
    self%pressure = 0E0_dp
    self%HTotal = 0E0_dp
    self%atoms = 0.0E0_dp
    self%dr = 0E0_dp
    self%maxMol = maxMol
    self%AtomType = 0
    self%MolType = 0
    self%MolIndx = 0
    self%MolSubIndx = 0
    self%AtomSubIndx = 0
    self%MolStartIndx = 0
    self%MolEndIndx = 0
    self%nMolTotal = 0

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
!          write(*,*) iType, molIndx, atmIndx
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
  subroutine SimpleBox_Destructor(self)
    implicit none
    class(SimpleBox), intent(inout) :: self

    deallocate(self%atoms)
    deallocate(self%dr)
    deallocate(self%ETable)
    deallocate(self%dETable)

    deallocate(self%AtomType)
    deallocate(self%MolType)
    deallocate(self%MolIndx)
    deallocate(self%MolSubIndx)
    deallocate(self%AtomSubIndx)

    deallocate(self%MolStartIndx)
    deallocate(self%MolEndIndx)
    deallocate(self%centerMass)

    deallocate(self%TypeFirst)
    deallocate(self%TypeLast)

    deallocate(self%MolGlobalIndx)
    deallocate(self%chempot)

    deallocate(self%NMol)
    deallocate(self%NMolMax)
    deallocate(self%NMolMin)


    if( allocated(self%Constrain) ) then
      deallocate(self%Constrain)
    endif

    if( allocated(self%NeighList) ) then
      deallocate(self%NeighList)
    endif

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
    IF (AllocateStatus /= 0) error STOP "*** SimpleBox: Unable to allocate Mol Bounds ***"
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
  subroutine SimpleBox_GetEFunc(self, epointer)
    ! Used to return a pointer array to the atoms(:,:) array within the box.
    ! Can also be used to simply return a subset of the atoms(:,:) such as a single
    ! molecule's coordinates
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(ECalcArray), pointer :: epointer

    epointer => self%EFunc
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


    !Safety block to ensure the user didn't reference something that doesn't exist
    !or is outside of the acceptable bounds.
    if((molType > nMolTypes) .or. (molType < 1)) then
      write(0,*) "ERROR! Type Index out of bounds!"
      write(0,*) molType, molIndx, atmIndx
      lineStat = -1
      return
    endif

    if( (molIndx > self%NMolMax(molType)) .or. (molIndx < 1) ) then
      write(0,*) "ERROR! Molecule Index out of bounds!"
      write(0,*) molType, molIndx, atmIndx
      lineStat = -1
      return
    endif

    if( (molIndx > self%NMol(molType)) .or. (molIndx < 1) ) then
      write(0,*) "ERROR! Molecule Index references a molecule that doesn't exist!"
      write(0,*) molType, molIndx, atmIndx
      lineStat = -1
      return
    endif

    if( (atmIndx > MolData(molType)%nAtoms) .or. (atmIndx < 1) ) then
      write(0,*) "ERROR! Atom Index out of bounds!"
      write(0,*) molType, molIndx, atmIndx
      lineStat = -1
      return
    endif

    call self%Boundary(x,y,z)


    subIndx = 0
    do iType = 1, molType-1
      subIndx = subIndx + self%NMolMax(iType)
    enddo
    subIndx = subIndx + molIndx
    arrayIndx = self%MolStartIndx(subIndx)
    arrayIndx = arrayIndx + atmIndx - 1
!    write(*,*) arrayIndx, trim(adjustl(line))
!    write(*,*) molType, subindx, self%MolStartIndx(subIndx)

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
    real(dp), intent(inout), optional :: rx, ry, rz 

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
    real(dp) :: E_Intra, E_Inter

    if(present(tablecheck))then
      if(tablecheck) then
        allocate(tempETable(1:size(self%ETable)))
      endif
      tempETable = self%ETable
    endif

    call self % ComputeIntraEnergy(E_Intra, accept)
    if(.not. accept) then
      return
    endif
    self%E_Intra = E_Intra

    write(nout,*) "Intra Energy:", E_Intra
    call self % EFunc % Method % DetailedECalc( self, E_Inter, accept )
    self%E_Inter = E_Inter
    self%ETotal = E_Inter + E_Intra

    if(present(tablecheck))then
      do iAtom = 1, self%nMaxAtoms
        if(.not. self%IsActive(iAtom)) cycle
        ECul = tempETable(iAtom)
        ECalc = self%ETable(iAtom)
        EDiff = abs(ECul-ECalc)
        if(ECalc > 1.0E0_dp) then
          if(ECalc /= 0E0_dp) then
            if( EDiff/abs(ECalc)  > 1E-4_dp ) then
              write(nout, *) "Table Discrepancy: ", iAtom, ECul/outEngUnit, ECalc/outEngUnit,  EDiff/outEngUnit, engStr
            endif
          else
            if( EDiff > 1E-5_dp ) then
              write(nout, *) "2 Table Discrepancy: ", iAtom, ECul/outEngUnit, ECalc/outEngUnit,  EDiff/outEngUnit, engStr
            endif
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
! This subroutine recomputes the intra contribution of entire system 
! energy from scratch. 
  subroutine SimpleBox_ComputeIntraEnergyDelta(self, disp, E_Intra)
    use Common_MolInfo, only:MolData
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(out) :: E_Intra
    logical :: accept
    integer :: iType, iMol, iAtom, iBond, iAngle, iTors
    integer :: iDisp
    integer :: molStart, molEnd, molIndx, molType
    integer :: intraType, nAtoms
    integer :: atomsubindx
    integer :: mem1, mem2, mem3, mem4
    logical, allocatable :: changed
    real(dp) :: E_Temp

    E_Temp = 0E0_dp
    E_Intra = 0E0_dp
    select type(disp)
      class is(Displacement)
         molType = disp(1)%molType 
         molIndx = disp(1)%molIndx 
         call self % GetMolData(molIndx, nAtoms=nAtoms, molStart=molStart, molEnd=molEnd)
         self%newpos(1:3, 1:nAtoms) = self%atoms(1:3, molStart:molEnd)
         do iDisp = 1, size(disp)
             atomsubindx = self%AtomSubIndx(disp(iDisp)%atmIndx)
             self%newpos(1,atomsubindx) = disp(iDisp)%x_new
             self%newpos(2,atomsubindx) = disp(iDisp)%y_new
             self%newpos(3,atomsubindx) = disp(iDisp)%z_new
         enddo
         call self%ComputeMolIntra(molType, molIndx, E_Temp, accept, self%newpos)
         if(.not. accept) then
           return
         endif
         E_Intra = E_Temp
         call self%ComputeMolIntra(molType, molIndx, E_Temp, accept)
         E_Intra = E_Intra - E_Temp

      class is(Addition)
         molType = disp(1)%molType 
         molIndx = disp(1)%molIndx 
         call self % GetMolData(molIndx, nAtoms=nAtoms, molStart=molStart, molEnd=molEnd)
         self%newpos(1:3, 1:nAtoms) = self%atoms(1:3, molStart:molEnd)
         do iDisp = 1, size(disp)
             atomsubindx = self%AtomSubIndx(disp(iDisp)%atmIndx)
             self%newpos(1,atomsubindx) = disp(iDisp)%x_new
             self%newpos(2,atomsubindx) = disp(iDisp)%y_new
             self%newpos(3,atomsubindx) = disp(iDisp)%z_new
         enddo
         call self%ComputeMolIntra(molType, molIndx, E_Temp, accept, self%newpos)
         E_Intra = E_Temp


      class is(Deletion)
         molType = disp(1)%molType 
         molIndx = disp(1)%molIndx 
         call self%ComputeMolIntra(molType, molIndx, E_Intra, accept)
         E_Intra = -E_Intra

      class is(OrthoVolChange)
         !Vol Change does not modify internal degrees of freedom
         !As the molecules are all translated by their center of mass.
         return

      class is(AtomExchange)
         !Atom Exchange assumes mono-atomic system
         return
    end select

  end subroutine
!==========================================================================================
! This subroutine computes the intra contribution of a single molecule
  subroutine SimpleBox_ComputeMolIntra(self, molType, molIndx, E_Intra, accept, newpos)
    use Common_MolInfo, only: nMolTypes, MolData, BondData, AngleData, &
                              TorsionData, mostAtoms
    use ErrorChecking, only: IsNan, IsInf
    implicit none
    class(SimpleBox), intent(inout), target :: self
    integer, intent(in) :: molType, molIndx
    real(dp), intent(out) :: E_Intra
    logical, intent(out) :: accept
    real(dp), intent(in), optional, target :: newpos(:,:)
    real(dp), pointer :: molpos(:,:)
    integer :: iAtom, iBond, iAngle, iTors, iMisc
    integer :: molStart, molEnd
    integer :: intraType, nAtoms
    integer :: mem1, mem2, mem3, mem4
    real(dp) :: E_Bond, E_Angle, E_Torsion, E_Misc


    E_Intra = 0E0_dp
    E_Bond = 0E0_dp
    E_Angle = 0E0_dp
    E_Torsion = 0E0_dp
    E_Misc = 0E0_dp
    if(MolData(molType)%ridgid) then
      return
    endif
    if( .not. allocated(self%temppos) ) then
      allocate( self%newpos(1:3, 1:mostAtoms)  )
      allocate( self%temppos(1:3, 1:mostAtoms) )
    endif
    nAtoms = MolData(molType)%nAtoms
    call self % GetMolData(molIndx, molStart=molStart, molEnd=molEnd)
    if(present(newpos)) then
        molpos => newpos(1:3, 1:nAtoms)
    else
        molpos => self%atoms(1:3, molStart:molEnd)
    endif

    do iBond = 1, MolData(molType)%nBonds
       intraType = MolData(molType)%bond(iBond)%bondType

       mem1 = MolData(molType)%bond(iBond)%mem1
       mem2 = MolData(molType)%bond(iBond)%mem2
       self % temppos(1:3, 1) = molpos(1:3, mem1)
       self % temppos(1:3, 2) = molpos(1:3, mem2)

       call BondData(intraType) % bondFF % DetailedECalc(self, self%temppos(1:3, 1:2), &
                                                         E_Bond, accept)
       if(.not. accept) then
         return
       endif
       if(IsNan(E_Bond) ) then
         write(0,*) "NaN Error in Bond Calculation!"
         error stop
       endif
       E_Intra = E_Intra + E_Bond
     enddo

     do iAngle = 1, MolData(molType)%nAngles
       intraType = MolData(molType)%angle(iAngle)%angleType

       mem1 = MolData(molType)%angle(iAngle)%mem1
       mem2 = MolData(molType)%angle(iAngle)%mem2
       mem3 = MolData(molType)%angle(iAngle)%mem3
       self%temppos(1:3, 1) = molpos(1:3, mem1)
       self%temppos(1:3, 2) = molpos(1:3, mem2)
       self%temppos(1:3, 3) = molpos(1:3, mem3)

       call AngleData(intraType) % angleFF % DetailedECalc(self, self%temppos(1:3,1:3), &
                                                         E_Angle, accept)

       if(.not. accept) then
         return
       endif
       if(IsNan(E_Angle) ) then
         write(0,*) "NaN Error in Angle Calculation!"
         error stop
       endif
       E_Intra = E_Intra + E_Angle
     enddo

     do iTors = 1, MolData(molType)%nTors
       intraType = MolData(molType)%torsion(iTors)%torsType
       mem1 = MolData(molType)%torsion(iTors)%mem1
       mem2 = MolData(molType)%torsion(iTors)%mem2
       mem3 = MolData(molType)%torsion(iTors)%mem3
       mem4 = MolData(molType)%torsion(iTors)%mem4
       self%temppos(1:3, 1) = molpos(1:3, mem1)
       self%temppos(1:3, 2) = molpos(1:3, mem2)
       self%temppos(1:3, 3) = molpos(1:3, mem3)
       self%temppos(1:3, 4) = molpos(1:3, mem4)
       call TorsionData(intraType) % torsionFF % DetailedECalc(self, &
                                                         self%temppos(1:3,1:4), &
                                                         E_Torsion, &
                                                         accept)
       if(.not. accept) then
         return
       endif
       if(IsNan(E_Torsion) ) then
         write(0,*) "NaN Error in Torsion Calculation!"
         error stop
       endif
       E_Intra = E_Intra + E_Torsion
     enddo


     self%temppos(1:3, 1:nAtoms) = molpos(1:3, 1:nAtoms)
     do iMisc = 1, MolData(molType)%nMisc
       call MolData(molType) % miscdata(iMisc) % miscFF % DetailedECalc(self, &
                                                         self%temppos(1:3,1:nAtoms), &
                                                         E_Misc, &
                                                         accept)
       if(.not. accept) then
         return
       endif
       if(IsNan(E_Misc) ) then
         write(0,*) "NaN Error in Misc Calculation!"
         error stop
       endif
       E_Intra = E_Intra + E_Misc
      enddo

  end subroutine

!==========================================================================================
! This subroutine recomputes the intra contribution of entire system 
! energy from scratch. 
  subroutine SimpleBox_ComputeIntraEnergy(self, E_Intra, accept)
    use Common_MolInfo, only: nMolTypes, MolData, BondData, AngleData, &
                              TorsionData
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(out) :: E_Intra
    logical, intent(out) :: accept
    integer :: iType, iMol, iAtom, iBond, iAngle, iTors
    integer :: molStart, molEnd, molIndx, molType
    real(dp) :: E_Mol

    E_Intra = 0E0_dp
    do iType = 1, nMolTypes
      if(MolData(iType)%ridgid) then
        cycle
      endif
      do iMol = 1, self%NMol(iType)
         molIndx = self % MolGlobalIndx(iType, iMol)
         call self%ComputeMolIntra(iType, molIndx, E_Mol, accept)
         if(.not. accept) then
           return
         endif
         E_Intra = E_Intra + E_Mol
      enddo
    enddo

  end subroutine
!==========================================================================================
! This subroutine computes the energy delta of a system primarily used to compute
! The change in the system that is a result of the various Monte Carlo moves.
!
!  Input
!    disp => Displacement vector which describes the result of the MC Move
!    templist => Temporary Neighborlist
!    tempNNei => Number of
!  Output
!    E_Inter => The energy of interaction between molecules
!    E_Itra => The energy of interaction within a molecule
!    accept => Rejection flag used for early rejection.
!    computeintra (Optional) => If true the intramolecular components
!                               will be recomputed.  For moves
!                               that do not change internal degrees of freedom
!                               set this to false or omit it. 
!     
  subroutine SimpleBox_ComputeEnergyDelta(self, disp, templist, tempNNei, E_Inter, E_Intra, E_Total, accept, computeintra)

    use Common_MolInfo, only: mostAtoms
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    logical, intent(in), optional :: computeintra
    real(dp), intent(inout) :: E_Inter, E_Intra, E_Total
    logical, intent(out) :: accept
    integer :: iAtom

    accept = .true.
    if( .not. allocated(self%temppos) ) then
      allocate(self%newpos(1:3, 1:mostAtoms))
      allocate(self%temppos(1:3, 1:mostAtoms))
    endif
    E_Total = 0E0_dp
    E_Intra = 0E0_dp
    E_Inter = 0E0_dp

    if(present(computeintra)) then
      if(computeintra) then
        call self%ComputeIntraEnergyDelta(disp, E_Intra)
      endif
      if(.not. accept) then
        return
      endif
    endif

    call self % EFunc % Method % DiffECalc(self, &
                                           disp, &
                                           tempList, &
                                           tempNNei, &
                                           E_Inter, &
                                           accept)



    E_Total = E_Inter + E_Intra

  end subroutine
!==========================================================================================
! This subroutine recomputes the forces for the system.
  subroutine SimpleBox_ComputeForces(self)
    use ParallelVar, only: nout
    use CoordinateTypes, only: Displacement
    use Units, only: outEngUnit, engStr
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    integer :: iAtom, atmType
    type(Displacement) :: disp(1:1)
    real(dp) :: E_Diff
    integer :: tempNnei(1)
    integer :: tempList(1, 1)


    if(.not. allocated(self%forces)) then
      allocate(self%forces(1:3, 1:self%nMaxAtoms))
    endif

    !Since force computation is very expensive, only do it if absolutely
    !nessisary.
    if(.not. self%forceoutofdate) then
      return
    endif

    self%forces = 0E0_dp


    do iAtom = 1, self%nMaxAtoms
      if(.not. self%isActive(iAtom)) cycle
      disp(1)%molType = self%MolType(iAtom)
      disp(1)%molIndx = self%MolIndx(iAtom)
      disp(1)%atmIndx = iAtom
      !Approximate the gradient by making a very small displacement on the atom position
      !and computing the slope (E2-E2)/(x2-x1).  Do this in each direction for the fx, fy, and fz components
      disp(1)%x_new = self%atoms(1, iAtom) + self%forcedelta
      disp(1)%y_new = self%atoms(2, iAtom) 
      disp(1)%z_new = self%atoms(3, iAtom)
      call self%EFunc%Method%DiffECalc(self, disp(1:1), tempList, tempNNei, E_Diff, accept)
      self%forces(1, iAtom) = E_diff/self%forcedelta

      disp(1)%x_new = self%atoms(1, iAtom) 
      disp(1)%y_new = self%atoms(2, iAtom) + self%forcedelta
      call self%EFunc%Method%DiffECalc(self, disp(1:1), tempList, tempNNei, E_Diff, accept)
      self%forces(2, iAtom) = E_diff/self%forcedelta

      disp(1)%y_new = self%atoms(2, iAtom) 
      disp(1)%z_new = self%atoms(3, iAtom) + self%forcedelta
      call self%EFunc%Method%DiffECalc(self, disp(1:1), tempList, tempNNei, E_Diff, accept)
      self%forces(3, iAtom) = E_diff/self%forcedelta
    enddo
    self%forceoutofdate = .false.
  end subroutine
!==========================================================================================
  subroutine SimpleBox_EnergySafetyCheck(self)
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    real(dp) :: E_Current, E_InterCur, E_IntraCur


    E_Current = self%ETotal
    E_InterCur = self%E_Inter
    E_IntraCur = self%E_Intra
    call self % ComputeEnergy(tablecheck=.false.)
!    call self % EFunc % Method % DetailedECalc( self, self%ETotal, accept )

    if(self%ETotal > 1E0_dp) then

      if( abs((E_Current-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "!!!!!!!!!!!!!!!!!!ERROR! Large energy drift detected!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Current/outEngUnit, engStr
        write(nout, *) "Culmative (Inter) Energy: ", E_InterCur/outEngUnit, engStr
        write(nout, *) "Culmative (Intra) Energy: ", E_IntraCur/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Final Energy(Inter): ", self%E_Inter/outEngUnit, engStr
        write(nout, *) "Final Energy(Intra): ", self%E_Intra/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Current)/outEngUnit, engStr
      endif
    else
      if( abs(E_Current-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "!!!!!!!!!!!!!!!!!!ERROR! Large energy drift detected!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Current/outEngUnit, engStr
        write(nout, *) "Culmative (Inter) Energy: ", E_InterCur/outEngUnit, engStr
        write(nout, *) "Culmative (Intra) Energy: ", E_IntraCur/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Final Energy(Inter): ", self%E_Inter/outEngUnit, engStr
        write(nout, *) "Final Energy(Intra): ", self%E_Intra/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Current)/outEngUnit, engStr
      endif
    endif



  end subroutine
!==========================================================================================
  subroutine SimpleBox_UpdateEnergy(self, E_Diff, E_Inter, E_Intra)
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: E_Diff, E_Inter
    real(dp), intent(in), optional :: E_Intra

    self % ETotal = self % ETotal + E_Diff
    self % ETable = self % ETable + self % dETable
    self % dETable = 0E0_dp
    self % E_Inter = self % E_Inter + E_Inter
    if(present(E_Intra)) then
      self % E_Intra = self % E_Intra + E_Intra
    endif

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
        if(.not. accept) then
          return
        endif
      enddo
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
        if(.not. accept) then
          return
        endif
      enddo
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

      case("energyrecompute")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self % forceERecompute = logicVal

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
  subroutine SimpleBox_GetMolData(self, globalIndx, nAtoms, molStart, molEnd, molType, &
                                  slice)
    !Routine takes a
    !
    !Inputs
    !    globalIndx => Global Molecule Index irrespective of molecule type.
    !                  This may not be populated depending on the state of the system
    !Outputs
    !    nAtoms => Number of Atoms in this molecule
    !    molStart => Array Index of the first Atom of this molecule in the atoms(:,:)
    !    molEnd => Array Index of the last Atom of this molecule in the atoms(:,:)
    !    slice => Returns an array with (molStart, molEnd) as data indexes. Used
    !             for getting array slices of individual molecules and such.
    !    molType => Molecule Type ID
    implicit none
    class(SimpleBox), intent(in) :: self
    integer, intent(in)  :: globalIndx
    integer, intent(inout), optional :: molStart, molEnd, molType
    integer, intent(inout), optional :: nAtoms
    integer, intent(inout), optional :: slice(1:2)
    integer :: tempInt



    if((size(self%molStartIndx) < globalIndx) .or. (globalIndx < 1)) then
      write(__StdErr__, *) "Error in Sim Box class: GetMolData has been given an invalid index"
      write(__StdErr__, *) "Given Index:", globalIndx 
      write(__StdErr__, *) "Array Size:", size(self%molStartIndx)
      error stop
    endif

    if( present(molStart) ) then
      molStart = self % MolStartIndx(globalIndx)
    endif

    if( present(molEnd) ) then
      molEnd = self % MolEndIndx(globalIndx)
    endif

    if(present(slice)) then
      slice(1) = self % MolStartIndx(globalIndx) 
      slice(2) = self % MolEndIndx(globalIndx) 
    endif

    if( present(nAtoms) ) then
      nAtoms = self % MolEndIndx(globalIndx) - self % MolStartIndx(globalIndx) + 1
    endif

    if( present(molType) ) then
      tempInt = self % MolStartIndx(globalIndx)
      molType = self % MolType(tempInt)
    endif

!    if( present(atomSubIndx) ) then
!      atomSubIndx = self%AtomSubIndx(globalIndx)
!    endif

  end subroutine
!====================================================================================
  function SimpleBox_GetNeighSkin(self, listindex) result(outval)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in)  :: listindex
    real(dp) :: outval

  

  end function
!====================================================================================
! Returns various indexing information 
  subroutine SimpleBox_GetAtomData(self, atomglobalIndx, molIndx, atomSubIndx, atomtype)
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in)  :: atomglobalIndx
    integer, intent(inout), optional :: molIndx, atomSubIndx, atomtype

!    integer, allocatable :: AtomType(:), MolType(:)
!    integer, allocatable :: MolIndx(:), MolSubIndx(:), AtomSubIndx(:)

    if((self%nMaxAtoms < atomglobalIndx) .or. (atomglobalIndx < 1)) then
      write(__StdErr__, *) "Error in Sim Box class: GetAtomData has been given an invalid index"
      write(__StdErr__, *) "Given Index:", atomglobalIndx 
      write(__StdErr__, *) "Array Size:", self%nMaxAtoms
      error stop
    endif

    if( present(molIndx) ) then
      molIndx = self % MolIndx(atomglobalindx)
    endif

    if( present(atomSubIndx) ) then
      atomSubIndx = self%AtomSubIndx(atomglobalindx)
    endif

    if( present(atomtype) ) then
      atomtype = self%AtomType(atomglobalindx)
    endif

  end subroutine

!==========================================================================================
function SimpleBox_FindMolByTypeIndex(self, MolType, nthofMolType, nthAtom) result(molstartindx)
    !------------------------------------------------------------------------------
    ! Returns the atom index of the first atom of the specified molecule by default.
    ! Can also return the n-th atom of the molecule 
    ! Search input requires the molecule type id and molecule sub-index to be given
    !
    ! Input:
    !   Integer: MolType => The molecule type ID to be searched for
    !   Integer: nthofMolType => The nth molecule of this type. Also known as a sub-index
    !   Integer, Optional: nthAtom => Can be specified if the user wants a different atom b
    !
    ! Output:
    !   Integer: molStartIndx => First atom of the molecule in this index
    ! Function Variables:
    !   Integer: iType => Loop Variable to iterate over the number of molecule types
    !   Integer: globalMolIndx => The value
    !   Integer: atomIndx => Local version of nthAtom. Set to 1 if nthAtom is not given.
    !------------------------------------------------------------------------------
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(SimpleBox), intent(inout) :: self
    integer, intent(in) :: MolType, nthofMolType
    integer, intent(in), optional :: nthAtom
    integer :: molStartIndx    

    integer :: iType
    integer :: globalMolIndx, atomIndx


    if(present(nthAtom)) then
      atomIndx = nthAtom
    else
      atomIndx = 1
    endif

    !Safety block to ensure the user didn't reference something that doesn't exist
    !or is outside of the acceptable bounds.
    if((molType > nMolTypes) .or. (molType < 1)) then
      write(0,*) "ERROR! Type Index out of bounds!"
      write(0,*) molType, nthofMolType, atomIndx
      error stop
    endif

    if( (nthofMolType > self%NMolMax(molType)) .or. (nthofMolType < 1) ) then
      write(0,*) "ERROR! Molecule Index out of bounds!"
      write(0,*) molType, nthofMolType, atomIndx
      error stop
    endif

    if( (nthofMolType > self%NMol(molType)) .or. (nthofMolType < 1) ) then
      write(0,*) "ERROR! Molecule Index references a molecule that doesn't exist!"
      write(0,*) molType, nthofMolType, atomIndx
      error stop
    endif

    globalMolIndx = 0
    do iType = 1, molType-1
      globalMolIndx = globalMolIndx + self%NMolMax(iType)
    enddo
    globalMolIndx = globalMolIndx + nthofMolType
    MolStartIndx = self%MolStartIndx(globalMolIndx)
    MolStartIndx = MolStartIndx + atomIndx - 1


  end function
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
  subroutine SimpleBox_GetForceArray(self, forces)
    implicit none
    class(SimpleBox), intent(inout), target :: self
    real(dp), pointer, intent(inout) :: forces(:,:)


    forces => self%forces


  end subroutine
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
  function SimpleBox_GetLargestNNei(self) result(maxnei)
    !-----------------
    ! Returns the Largest "Max Neighbor Count" for a single atom across all lists
    !-----------------
    implicit none
    class(SimpleBox), intent(inout), target :: self
    integer :: maxnei 
    integer :: curmaxnei = 0
    integer :: iList

    maxnei = 0

    if(allocated(self%NeighList)) then
      do iList = 1, size(self%NeighList)
         curmaxnei = self%NeighList(iList)%GetMaxNei()
         if(maxnei < curmaxnei) then
           maxnei = curmaxnei
         endif
      enddo
    endif

  end function 
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
    class(SimpleBox), intent(in) :: self
    logical :: active
    integer, intent(in) :: atmIndx

    if( self%MolSubIndx(atmIndx) > self%NMol(self%MolType(atmIndx)) ) then
      active = .false.
    else
      active = .true.
    endif

  end function
!==========================================================================================
  subroutine SimpleBox_SwapAtomType(self, disp)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: self
    class(AtomExchange), intent(inout) :: disp(:)
    integer :: iAtom, iAtomNew, iAtomOld
    integer :: iList, iType, jStart, iDimn
    integer :: molIndxNew, molIndxOld
    integer :: nStartNew, nStartOld, lastMol, nTypeNew, nTypeOld

    iAtomNew = disp(1)%newAtmIndx
    iAtomOld = disp(1)%oldAtmIndx

    molIndxOld = self%MolStartIndx(iAtomOld)
    molIndxNew = self%MolStartIndx(iAtomNew)

    nStartNew = self % MolStartIndx(molIndxNew)
    nStartOld = self% MolStartIndx(molIndxOld)
    nTypeNew = self % MolType(nStartNew)
    nTypeOld = self % MolType(nStartOld)
    if(MolData(nTypeNew)%nAtoms > 1) then
      write(0,*) "INVALID. An AtomExchange Perturbation was passed into the update function containing"
      write(0,*) "the index of a molecule type with more than one atom present.  This perturbation was"
      write(0,*) "only designed for single atom molecules"
      error stop
    endif
    if(MolData(nTypeOld)%nAtoms > 1) then
      write(0,*) "INVALID. An AtomExchange Perturbation was passed into the update function containing"
      write(0,*) "the index of a molecule type with more than one atom present.  This perturbation was"
      write(0,*) "only designed for single atom molecules"
      error stop
    endif


    lastMol = 0
    do iType = 1, nTypeOld-1
      lastMol = lastMol + self%NMolMax(iType) 
    enddo
    lastMol = lastMol + self%NMol(nTypeOld)
    jStart = self % MolStartIndx(lastMol)


    self%centermass(1:3,molIndxNew) = self%centermass(1:3, molIndxOld)
    self%centermass(1:3,molIndxOld) = self%centermass(1:3, lastMol)

!     Take the top molecule from the atom array and move it's position in the deleted
!     molecule's slot.

    do iDimn = 1, self%nDimension
      self % atoms(iDimn, nStartNew ) = self % atoms(iDimn, nStartOld)
      self % dr(iDimn, nStartNew ) = self % dr(iDimn, nStartOld )
    enddo
    self % ETable(nStartNew) = self % ETable(nStartOld)

    do iDimn = 1, self%nDimension
      self % atoms(iDimn, nStartOld) = self % atoms(iDimn, jStart)
      self % dr(iDimn, nStartOld) = self % dr(iDimn, jStart)
    enddo
    self % ETable(nStartOld) = self % ETable(jStart)
    self % ETable(jStart) = 0E0_dp


    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % SwapAtomType(disp, lastMol)
    enddo

    self % NMol(nTypeNew) = self % NMol(nTypeNew) + 1 
    self % NMol(nTypeOld) = self % NMol(nTypeOld) - 1 
!    self % nMolTotal = self % nMolTotal 
!    self % nAtoms = self % nAtoms - MolData(nType)%nAtoms


  end subroutine
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

    if(self%NMol(nType) == self%NMolMin(nType)) then
      write(0,*) "ERROR! Box attempted to delete a molecule, but this"
      write(0,*) "       would take it below the minimum bounds for this box" 
      write(0,*) "Box Number:", self%boxID
      write(0,*) "MolType:", nType
      error stop 
    endif


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
!    write(*,*) molindx
!    write(*,*) "e", self%ETable
    do iAtom = 1, MolData(nType) % nAtoms
      do iDimn = 1, self%nDimension
        self % atoms(iDimn, nStart+iAtom-1 ) = self % atoms(iDimn, jStart+iAtom-1 )
        self % dr(iDimn, nStart+iAtom-1 ) = self % dr(iDimn, jStart+iAtom-1 )
      enddo
      self % ETable(nStart+iAtom-1 ) = self % ETable(jStart+iAtom-1 )
      self % ETable(jStart+iAtom-1 ) = 0E0_dp
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
    integer :: molType
    integer :: iDisp, dispLen, dispIndx
    integer :: iAtom, iMol, molStart, molEnd
    real(dp) :: dx, dy, dz
    
    self%forceoutofdate = .true.
    select type(disp)
       !-------------------------------------------------
      class is(Displacement)
        dispLen = size(disp)
        do iDisp = 1, dispLen
          dispIndx = disp(iDisp) % atmIndx
          dx = disp(iDisp)%x_new - self%atoms(1, dispIndx)
          dy = disp(iDisp)%y_new - self%atoms(2, dispIndx)
          dz = disp(iDisp)%z_new - self%atoms(3, dispIndx)
          call self%Boundary( dx, dy, dz )
          self % dr(1, dispIndx) = self % dr(1, dispIndx) + dx
          self % dr(2, dispIndx) = self % dr(2, dispIndx) + dy 
          self % dr(3, dispIndx) = self % dr(3, dispIndx) + dz
          call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
          self % atoms(1, dispIndx) = disp(iDisp)%x_new
          self % atoms(2, dispIndx) = disp(iDisp)%y_new
          self % atoms(3, dispIndx) = disp(iDisp)%z_new
        enddo
        call self%ComputeCM(disp(1)%molIndx)


       !-------------------------------------------------
      class is(Addition)
        dispLen = size(disp)
        molType = disp(1)%molType
        if(self%NMol(molType) == self%NMolMax(molType)) then
          write(0,*) "ERROR! Box attempted to add a molecule, but this"
          write(0,*) "       would take it above the maximum bounds for this box" 
          write(0,*) "Box Number:", self%boxID
          write(0,*) "MolType:", molType
          error stop 
        endif
        do iDisp = 1, dispLen
          dispIndx = disp(iDisp) % atmIndx
          call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
          self % atoms(1, dispIndx) = disp(iDisp)%x_new
          self % atoms(2, dispIndx) = disp(iDisp)%y_new
          self % atoms(3, dispIndx) = disp(iDisp)%z_new
          self % dr(1:3, dispIndx) = 0E0_dp
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
              self % dr(1, iAtom) = self % dr(1, iAtom) + dx
              self % dr(2, iAtom) = self % dr(2, iAtom) + dy
              self % dr(3, iAtom) = self % dr(3, iAtom) + dz
!              write(*,*) self%atoms(1:3, iAtom)
            enddo
            self%centerMass(1, iMol) = self%centerMass(1, iMol) * disp(1)%xScale
            self%centerMass(2, iMol) = self%centerMass(2, iMol) * disp(1)%yScale
            self%centerMass(3, iMol) = self%centerMass(3, iMol) * disp(1)%zScale
          endif
        enddo
        call self%UpdateVolume(disp)
      class is(AtomExchange)
        call self%SwapAtomType(disp)

       !-------------------------------------------------
      class default
        error stop "The code does not know how to update coordinates for this perturbation type."
    end select


  end subroutine
!==========================================================================================
  subroutine SimpleBox_Maintenance(self)
    use CoordinateTypes
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: self
    logical :: accept
    integer :: iList
    integer :: iAtom

!    return
!    do iList = 1, size(self%NeighList)
!      call self % NeighList(iList) % BuildList(iList)
!    enddo


    if(self%forceERecompute) then
      call self % ComputeEnergy(tablecheck=.false.)
    endif

!    call self % EFunc % Method % DetailedECalc( self, tempE, accept )

!    if(abs(self%ETotal-tempE)/abs(tempE) > 1e-4) then
!   write(2,*) self%ETotal, tempE, self%ETotal-tempE
!    endif

  end subroutine
!==========================================================================================
!  Checks the Neighbor List criteria periodically to see if the lists need to be rebuilt.
  subroutine SimpleBox_CheckLists(self)
    implicit none
    class(SimpleBox), intent(inout) :: self

    integer :: iList
    integer :: iAtom, iMol
    integer :: molStart, molEnd
    real(dp) :: tempE, drsq
    real(dp) :: maxdr, maxdr2, molMax
    real(dp) :: neighSkin, E_RCut, Nei_RCut


    !Check to see if the particles have shifted enough from their original position to justify a
    !neighorlist rebuild.
    maxdr = 0E0_dp
    maxdr2 = 0E0_dp
    do iMol = 1, self%maxMol
      molmax = 0E0_dp
      call self%GetMolData(iMol, molStart=molStart, molEnd=molEnd)
      do iAtom = molStart, molEnd
        if(.not. self%IsActive(iAtom)) cycle
        drsq = self%dr(1, iAtom)*self%dr(1, iAtom) + self%dr(2, iAtom)*self%dr(2, iAtom) + &
               self%dr(3, iAtom)*self%dr(3, iAtom)
!        write(*,*) iAtom, drsq
        molmax = max(drsq, molmax)
      enddo
      if(molmax > maxdr) then
        maxdr2 = maxdr
        maxdr = molmax
      else
        if(molmax > maxdr2) then
          maxdr2 = molmax
        endif
      endif
    enddo

    maxdr = sqrt(maxdr)
    maxdr2 = sqrt(maxdr2)

!    write(*,*) maxdr, maxdr2
    !Keep track of the largest displacement between rebuilds seen through the course of the simulation
    E_rcut = self%EFunc%Method%GetCutOff()
    do iList = 1, size(self%NeighList)
      Nei_RCut = self%NeighList(iList)%GetRCut()
      if( (iList == 1) .and. (E_rcut > Nei_RCut)) then
        error stop "ERROR! User set the first neighbor list's RCut shorter than required by the Forcefield!!"
      endif
      neighSkin = Nei_RCut - E_rcut
      if( (maxdr + maxdr2) > neighSkin ) then
!        write(*,*) maxdr, maxdr2, neighSkin
        self%rebuilds = self%rebuilds + 1
        if(maxdr > neighskin*0.7E0_dp) then
          self%dangerbuilds = self%dangerbuilds + 1
  !        write(__StdErr__, *) "Warning, Dangerous Neighborlist Build Detected!"
        endif
        if(maxdr > self%largestdr) then
          self%largestdr = maxdr
        endif
        call self % NeighList(iList) % BuildList(iList)
        self%dr = 0E0_dp
      endif
    enddo
  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetReducedCoords(self,realCoords,reducedCoords )
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: realCoords(:)
    real(dp), intent(out) :: reducedCoords(1:3)

    reducedCoords = 0E0_dp
    error stop "GetReducedCoords routine has been called on a system where no box is defined."

  end subroutine
!==========================================================================================
  subroutine SimpleBox_GetRealCoords(self, reducedCoords, realCoords)
    implicit none
    class(SimpleBox), intent(inout) :: self
    real(dp), intent(in) :: reducedCoords(:)
    real(dp), intent(out) :: realCoords(1:3)

    realCoords = 0E0_dp
    error stop "GetReducedCoords routine has been called on a system where no box is defined."
  end subroutine
!==========================================================================================
  subroutine SimpleBox_ScreenOut(self)
    use Common_MolInfo, only: nMolTypes
    use ParallelVar, only: nout
    use Input_Format, only: ReplaceText
    use Units, only: outEngUnit, outLenUnit
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

    write(tempStr, "(A)") "        Per Mol Energy: %s1    Volume: %s2"
    write(tempStr2, "(F40.8)") self%ETotal/(outEngUnit*self%nMolTotal)
    tempStr = ReplaceText(tempStr, "%s1", trim(adjustl(tempStr2)))
    write(tempStr2, "(F40.8)") self%Volume/(outLenUnit**3)
    tempStr = ReplaceText(tempStr, "%s2", trim(adjustl(tempStr2)))
    write(nout, "(A)") trim(tempStr)

    write(nout, "(A, 999(1x,I10))") "        Number of Molecules By Type: ", self%nMol(1:nMolTypes)




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
    write(nout,*) "Box ", self%boxID, " Molecule Count: ", (self % NMol(iMol), iMol=1,nMolTypes)
    write(nout,*) "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal
    write(nout,*) "Box ", self%boxID, " Temperature: ", self % temperature

    call self % ComputeEnergy
    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % BuildList(iList)
    enddo
    if( allocated(self%Constrain) ) then
      if( size(self%Constrain) > 0 ) then
        do iConstrain = 1, size(self%Constrain)
          call self%Constrain(iConstrain) % method % Prologue
          call self%Constrain(iConstrain) % method % CheckInitialConstraint(self, accept)
        enddo
        if(.not. accept) then
          write(nout,*) "Initial Constraints are not statisfied!"
          error stop
        endif
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
    real(dp) :: E_Culm, E_Culm_Inter, E_Culm_Intra

    E_Culm = self%ETotal
    E_Culm_Inter = self%E_Inter
    E_Culm_Intra = self%E_Intra



    write(nout,*) "--------Box", self%boxID , "Energy---------"
    if(isnan(E_Culm)) then
      write(nout, *) "ERROR! Final Culmative Energy is not a number!"
    endif
    call self % ComputeEnergy(tablecheck=.false.)
    do iList = 1, size(self%NeighList)
      call self % NeighList(iList) % BuildList(iList)
    enddo

    write(nout, *) "Final Energy:", self % ETotal/outEngUnit, engStr
    if(self%ETotal /= 0) then
      if( abs((E_Culm-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit, engStr
        write(nout, *) "Culmative (Inter) Energy: ", E_Culm_Inter/outEngUnit, engStr
        write(nout, *) "Culmative (Intra) Energy: ", E_Culm_Intra/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Final (Inter) Energy: ", self%E_Inter/outEngUnit, engStr
        write(nout, *) "Final (Intra) Energy: ", self%E_Intra/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit, engStr
      endif
    else
      if( abs(E_Culm-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit, engStr
        write(nout, *) "Culmative (Inter) Energy: ", E_Culm_Inter/outEngUnit, engStr
        write(nout, *) "Culmative (Intra) Energy: ", E_Culm_Intra/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Final (Inter) Energy: ", self%E_Inter/outEngUnit, engStr
        write(nout, *) "Final (Intra) Energy: ", self%E_Intra/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit, engStr
      endif
    endif
    write(nout, "(1x,A4,I2,A,E15.8,1x,A)") "Box ", self%boxID, " Final Energy (Per Mol): ", &
                                           self % ETotal/(outEngUnit*self%nMolTotal), engStr
    write(nout, "(1x,A4,I2,A,I8)") "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout, "(1x,A4,I2,A,I8)") "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal
    write(nout, "(1x,A4,I2,A,I8)") "Box ", self%boxID, " Neigh Rebuilds: ", self % rebuilds
    write(nout, "(1x,A4,I2,A,I8)") "Box ", self%boxID, " Dangerous Builds: ", self % dangerbuilds
    write(nout, "(1x,A4,I2,A,E15.8)") "Box ", self%boxID, " Largest Atom Displacement Between Neigh Builds: ", self % largestdr

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
      error stop "No Neighbor List has been defined!"
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
