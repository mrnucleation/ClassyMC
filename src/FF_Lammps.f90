!=============================================================================+
! LAMMPS forcefield interface module. This module provides an interface to
! LAMMPS' energy calculation routines through the LAMMPS C library API.
! 
! BUILDING:
! ---------
! To add LAMMPS to the compilation:
!   1. Build LAMMPS as a shared library (liblammps.so):
!      cd lammps/build && cmake -D BUILD_SHARED_LIBS=on ../cmake && make
!   2. Compile ClassyMC with the LAMMPS target:
!      make lammps LAMMPS_LIB_PATH=/path/to/lammps/build
!   3. Set LD_LIBRARY_PATH to include the LAMMPS library path
!
! The LAMMPS library documentation can be found at:
! https://docs.lammps.org/Library.html
!
! USAGE IN INPUT FILE:
! --------------------
! In your ClassyMC input file, specify the LAMMPS forcefield:
!
!   forcefield 1
!     lammps
!     rcut 10.0
!     lammps_cmd create_box 2 box
!     lammps_cmd mass 1 12.0
!     lammps_cmd mass 2 16.0
!     lammps_cmd pair_style lj/cut 10.0
!     lammps_cmd pair_coeff 1 1 0.1 3.4
!     lammps_cmd pair_coeff 1 2 0.15 3.2
!     lammps_cmd pair_coeff 2 2 0.2 3.0
!     energy_conversion 0.0433641  # kcal/mol to eV
!     1 0.5  # atom type 1 rMin = 0.5
!     2 0.5  # atom type 2 rMin = 0.5
!   end
!
! AVAILABLE COMMANDS:
! -------------------
!   rcut <value>              - Set the cutoff distance
!   lammps_cmd <command>      - Execute a LAMMPS command during initialization
!   inputscript <filename>    - Load commands from a LAMMPS input script file
!   energy_conversion <value> - Factor to convert LAMMPS energy to ClassyMC units
!   length_conversion <value> - Factor to convert ClassyMC length to LAMMPS units
!   typemap <classy_type> <lammps_type> - Map atom types between systems
!   <type> <rMin>             - Set minimum distance for atom type
!   <type1> <type2> <rMin>    - Set minimum distance between two atom types
!
! C LIBRARY INTERFACE:
! --------------------
! This module uses ISO_C_BINDING to interface with the LAMMPS C library.
! The key LAMMPS library functions used are:
!   - lammps_open_no_mpi: Create LAMMPS instance without MPI
!   - lammps_close: Close LAMMPS instance
!   - lammps_command: Execute a single LAMMPS command
!   - lammps_scatter_atoms: Set atomic positions
!   - lammps_gather_atoms: Get atomic positions
!   - lammps_extract_compute: Get compute results (e.g., potential energy)
!   - lammps_reset_box: Update simulation box dimensions
!
!=============================================================================+
module FF_Ext_LAMMPS
  use Template_ForceField, only: ForceField
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use Input_Format, only: maxLineLen
  use VarPrecision
  use CoordinateTypes
  use iso_c_binding

  implicit none
  private

#ifdef LAMMPS
  !---------------------------------------------------------------------------
  ! LAMMPS C Library Interface Declarations
  !---------------------------------------------------------------------------
  interface
    ! Create a LAMMPS instance without MPI (returns handle directly)
    function lammps_open_no_mpi(argc, argv, ptr) bind(C, name='lammps_open_no_mpi')
      import :: c_ptr, c_int
      integer(c_int), value :: argc
      type(c_ptr), value :: argv
      type(c_ptr), value :: ptr  ! Pass c_null_ptr to use return value
      type(c_ptr) :: lammps_open_no_mpi
    end function lammps_open_no_mpi

    ! Close a LAMMPS instance
    subroutine lammps_close(ptr) bind(C, name='lammps_close')
      import :: c_ptr
      type(c_ptr), value :: ptr
    end subroutine lammps_close

    ! Execute a single LAMMPS command
    subroutine lammps_command(ptr, cmd) bind(C, name='lammps_command')
      import :: c_ptr, c_char
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: cmd
    end subroutine lammps_command

    ! Execute multiple LAMMPS commands from a string
    subroutine lammps_commands_string(ptr, str) bind(C, name='lammps_commands_string')
      import :: c_ptr, c_char
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: str
    end subroutine lammps_commands_string

    ! Get number of atoms
    function lammps_get_natoms(ptr) bind(C, name='lammps_get_natoms')
      import :: c_ptr, c_double
      type(c_ptr), value :: ptr
      real(c_double) :: lammps_get_natoms
    end function lammps_get_natoms

    ! Scatter atoms - set positions from Fortran array to LAMMPS
    subroutine lammps_scatter_atoms(ptr, name, type, count, data) bind(C, name='lammps_scatter_atoms')
      import :: c_ptr, c_char, c_int, c_double
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: name
      integer(c_int), value :: type
      integer(c_int), value :: count
      real(c_double), dimension(*), intent(in) :: data
    end subroutine lammps_scatter_atoms

    ! Gather atoms - get positions from LAMMPS to Fortran array
    subroutine lammps_gather_atoms(ptr, name, type, count, data) bind(C, name='lammps_gather_atoms')
      import :: c_ptr, c_char, c_int, c_double
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: name
      integer(c_int), value :: type
      integer(c_int), value :: count
      real(c_double), dimension(*), intent(out) :: data
    end subroutine lammps_gather_atoms

    ! Extract a compute value
    function lammps_extract_compute(ptr, id, style, type) bind(C, name='lammps_extract_compute')
      import :: c_ptr, c_char, c_int
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: id
      integer(c_int), value :: style
      integer(c_int), value :: type
      type(c_ptr) :: lammps_extract_compute
    end function lammps_extract_compute

    ! Extract a global value (like thermo quantities)
    function lammps_extract_global(ptr, name) bind(C, name='lammps_extract_global')
      import :: c_ptr, c_char
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: name
      type(c_ptr) :: lammps_extract_global
    end function lammps_extract_global

    ! Extract variable
    function lammps_extract_variable(ptr, name, group) bind(C, name='lammps_extract_variable')
      import :: c_ptr, c_char, c_double
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: name
      character(kind=c_char), dimension(*), intent(in) :: group
      real(c_double) :: lammps_extract_variable
    end function lammps_extract_variable

    ! Get thermo value by keyword (modern API for getting PE, KE, etc.)
    function lammps_get_thermo(ptr, name) bind(C, name='lammps_get_thermo')
      import :: c_ptr, c_char, c_double
      type(c_ptr), value :: ptr
      character(kind=c_char), dimension(*), intent(in) :: name
      real(c_double) :: lammps_get_thermo
    end function lammps_get_thermo

    ! Reset the simulation box
    subroutine lammps_reset_box(ptr, boxlo, boxhi, xy, yz, xz) bind(C, name='lammps_reset_box')
      import :: c_ptr, c_double
      type(c_ptr), value :: ptr
      real(c_double), dimension(3), intent(in) :: boxlo
      real(c_double), dimension(3), intent(in) :: boxhi
      real(c_double), value :: xy, yz, xz
    end subroutine lammps_reset_box

    ! Create atoms in LAMMPS
    subroutine lammps_create_atoms(ptr, n, id, type, x, v, image, shrinkexceed) &
        bind(C, name='lammps_create_atoms')
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: ptr
      integer(c_int), value :: n
      type(c_ptr), value :: id
      type(c_ptr), value :: type
      real(c_double), dimension(*), intent(in) :: x
      type(c_ptr), value :: v
      type(c_ptr), value :: image
      integer(c_int), value :: shrinkexceed
    end subroutine lammps_create_atoms

  end interface
#endif

  !---------------------------------------------------------------------------
  ! LAMMPS Forcefield Type Definition
  !---------------------------------------------------------------------------
  type, public, extends(EasyPair_Cut) :: Pair_LAMMPS
    logical :: initialized = .false.
    logical :: lammps_created = .false.
    type(c_ptr) :: lmp = c_null_ptr
    
    ! Storage for LAMMPS input script commands
    character(len=maxLineLen), allocatable :: setupCommands(:)
    integer :: nSetupCommands = 0
    character(len=maxLineLen) :: inputScript = ""
    
    ! Atom type mapping (Classy type -> LAMMPS type)
    integer, allocatable :: typeMap(:)
    
    ! Temporary coordinate storage
    real(dp), allocatable :: tempcoords(:)
    integer, allocatable :: atomTypes(:)
    
    ! Box dimensions for LAMMPS
    real(dp) :: boxlo(3) = 0.0_dp
    real(dp) :: boxhi(3) = 10.0_dp
    
    ! Energy conversion factor (LAMMPS units -> Classy units)
    real(dp) :: energyConversion = 1.0_dp
    real(dp) :: lengthConversion = 1.0_dp
    
    contains
#ifdef LAMMPS
      procedure, pass :: Constructor => Constructor_LAMMPS
      procedure, pass :: DetailedECalc => DetailedECalc_LAMMPS
      procedure, pass :: DiffECalc => DiffECalc_LAMMPS
      procedure, pass :: ProcessIO => ProcessIO_LAMMPS
      procedure, pass :: Prologue => Prologue_LAMMPS
      procedure, pass :: Epilogue => Epilogue_LAMMPS
      procedure, pass :: InitLAMMPS => InitLAMMPS_LAMMPS
      procedure, pass :: UpdatePositions => UpdatePositions_LAMMPS
      procedure, pass :: GetEnergy => GetEnergy_LAMMPS
#endif
  end type

contains

#ifdef LAMMPS
!=============================================================================+
  subroutine Constructor_LAMMPS(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    integer :: AllocateStat

    write(nout, *) "Initializing LAMMPS Forcefield Interface"
    
    ! Allocate rMin tables (inherited from EasyPair_Cut)
    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%typeMap(1:nAtomTypes), stat = AllocateStat)
    
    self%rMin = 0.0_dp
    self%rMinTable = 0.0_dp
    
    ! Default type mapping: Classy type = LAMMPS type
    do AllocateStat = 1, nAtomTypes
      self%typeMap(AllocateStat) = AllocateStat
    enddo
    
    ! Allocate setup commands array
    allocate(self%setupCommands(100), stat = AllocateStat)
    self%nSetupCommands = 0
    
  end subroutine
!=============================================================================+
  subroutine InitLAMMPS_LAMMPS(self, curbox)
    use ParallelVar, only: nout
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    
    integer :: i, iAtom, nAtoms, idx
    integer(c_int), allocatable, target :: atomIDs(:), atomTypesArr(:)
    character(len=512) :: cmd
    real(dp) :: tempdim(3,3)
    real(dp), pointer :: atoms(:,:) => null()
    logical :: hasUnits, hasAtomStyle
    
    ! For suppressing LAMMPS output
    integer, parameter :: nargs = 5
    character(len=16), target :: args(nargs)
    type(c_ptr), target :: argv(nargs)
    
    if (self%lammps_created) return
    
    write(nout, *) "Creating LAMMPS instance..."
    
    ! Setup command-line arguments to suppress LAMMPS screen/log output
    args(1) = "lammps" // c_null_char
    args(2) = "-screen" // c_null_char
    args(3) = "none" // c_null_char
    args(4) = "-log" // c_null_char
    args(5) = "lmp.log" // c_null_char
    
    do i = 1, nargs
      argv(i) = c_loc(args(i))
    enddo
    
    ! Create LAMMPS instance with suppressed output (use return value, pass c_null_ptr for 3rd arg)
    self%lmp = lammps_open_no_mpi(int(nargs, c_int), c_loc(argv), c_null_ptr)
    
    if (.not. c_associated(self%lmp)) then
      error stop "Failed to create LAMMPS instance"
    endif
    
    self%lammps_created = .true.
    
    ! Get box dimensions from simulation box
    select type(curbox)
      class is(CubeBox)
        call curbox%GetDimensions(tempdim)
        self%boxlo(1) = tempdim(1,1) * self%lengthConversion
        self%boxlo(2) = tempdim(1,2) * self%lengthConversion
        self%boxlo(3) = tempdim(1,3) * self%lengthConversion
        self%boxhi(1) = tempdim(2,1) * self%lengthConversion
        self%boxhi(2) = tempdim(2,2) * self%lengthConversion
        self%boxhi(3) = tempdim(2,3) * self%lengthConversion
      class is(OrthoBox)
        call curbox%GetDimensions(tempdim)
        self%boxlo(1) = tempdim(1,1) * self%lengthConversion
        self%boxlo(2) = tempdim(1,2) * self%lengthConversion
        self%boxlo(3) = tempdim(1,3) * self%lengthConversion
        self%boxhi(1) = tempdim(2,1) * self%lengthConversion
        self%boxhi(2) = tempdim(2,2) * self%lengthConversion
        self%boxhi(3) = tempdim(2,3) * self%lengthConversion
    end select
    
    !write(nout, *) "DEBUG: ClassyMC box dimensions (raw tempdim):"
    !write(nout, *) "  X: ", tempdim(1,1), " to ", tempdim(2,1)
    !write(nout, *) "  Y: ", tempdim(1,2), " to ", tempdim(2,2)
    !write(nout, *) "  Z: ", tempdim(1,3), " to ", tempdim(2,3)
    !write(nout, *) "  Length conversion: ", self%lengthConversion
    !write(nout, *) "DEBUG: LAMMPS box will be:"
    !write(nout, *) "  X: ", self%boxlo(1), " to ", self%boxhi(1)
    !write(nout, *) "  Y: ", self%boxlo(2), " to ", self%boxhi(2)
    !write(nout, *) "  Z: ", self%boxlo(3), " to ", self%boxhi(3)
    
    ! Execute early setup commands (units, atom_style) from user first, or use defaults
    ! These MUST come before region/box definition in LAMMPS
    hasUnits = .false.
    hasAtomStyle = .false.
    
    ! Check if user specified units or atom_style, and execute them first
    do i = 1, self%nSetupCommands
      if (index(self%setupCommands(i), "units ") == 1) then
        call lammps_command(self%lmp, trim(self%setupCommands(i)) // c_null_char)
        hasUnits = .true.
      else if (index(self%setupCommands(i), "atom_style ") == 1) then
        call lammps_command(self%lmp, trim(self%setupCommands(i)) // c_null_char)
        hasAtomStyle = .true.
      endif
    enddo
    
    ! Apply defaults only if user didn't specify
    if (.not. hasUnits) then
    call lammps_command(self%lmp, "units real" // c_null_char)
    endif
    if (.not. hasAtomStyle) then
    call lammps_command(self%lmp, "atom_style atomic" // c_null_char)
    endif
    call lammps_command(self%lmp, "atom_modify map array sort 0 0.0" // c_null_char)
    
    ! Create simulation box region
    write(cmd, '(A,6(F12.4,1X))') "region box block ", &
        self%boxlo(1), self%boxhi(1), &
        self%boxlo(2), self%boxhi(2), &
        self%boxlo(3), self%boxhi(3)
    call lammps_command(self%lmp, trim(cmd) // c_null_char)
    
    ! Execute remaining user-defined setup commands (create_box, mass, pair_style, pair_coeff)
    !write(nout, *) "DEBUG: Executing user LAMMPS commands:"
    do i = 1, self%nSetupCommands
      ! Skip units and atom_style commands (already executed above)
      if (index(self%setupCommands(i), "units ") == 1) cycle
      if (index(self%setupCommands(i), "atom_style ") == 1) cycle
      write(nout, *) "  CMD: ", trim(self%setupCommands(i))
      call lammps_command(self%lmp, trim(self%setupCommands(i)) // c_null_char)
    enddo
    
    ! Now create atoms in LAMMPS based on ClassyMC simulation box
    ! This must be done AFTER create_box but BEFORE run/scatter_atoms
    call curbox%GetCoordinates(atoms)
    
    ! Count active atoms and prepare position/type arrays
    nAtoms = 0
    do iAtom = 1, curbox%nMaxAtoms
      if (curbox%IsActive(iAtom)) then
        nAtoms = nAtoms + 1
      endif
    enddo
    
    write(nout, *) "Creating", nAtoms, "atoms in LAMMPS..."
    
    if (nAtoms > 0) then
      ! Ensure tempcoords is allocated (normally done in Prologue)
      if (.not. allocated(self%tempcoords)) then
        allocate(self%tempcoords(3 * curbox%nMaxAtoms))
      endif
      
      ! Allocate temporary arrays for atom creation
      allocate(atomIDs(nAtoms))
      allocate(atomTypesArr(nAtoms))
      
      ! Fill arrays with atom data
      idx = 0
      do iAtom = 1, curbox%nMaxAtoms
        if (curbox%IsActive(iAtom)) then
          idx = idx + 1
          atomIDs(idx) = idx  ! Consecutive IDs starting from 1
          atomTypesArr(idx) = self%typeMap(curbox%AtomType(iAtom))
          ! Store positions in tempcoords (allocated in Prologue)
          self%tempcoords((idx-1)*3 + 1) = atoms(1, iAtom) * self%lengthConversion
          self%tempcoords((idx-1)*3 + 2) = atoms(2, iAtom) * self%lengthConversion
          self%tempcoords((idx-1)*3 + 3) = atoms(3, iAtom) * self%lengthConversion
          ! Debug: print first few atom positions
          !if (idx <= 5) then
            !write(nout, *) "DEBUG: Atom", idx, "type", atomTypesArr(idx), "pos:", &
            !    self%tempcoords((idx-1)*3 + 1), &
            !    self%tempcoords((idx-1)*3 + 2), &
            !    self%tempcoords((idx-1)*3 + 3)
          !endif
        endif
      enddo
      
      ! Create atoms using LAMMPS API
      call lammps_create_atoms(self%lmp, int(nAtoms, c_int), &
                               c_loc(atomIDs), c_loc(atomTypesArr), &
                               self%tempcoords, c_null_ptr, c_null_ptr, 0_c_int)
      
      ! Verify LAMMPS received the atoms
      !write(nout, *) "DEBUG: LAMMPS reports", int(lammps_get_natoms(self%lmp)), "atoms"
      
      deallocate(atomIDs)
      deallocate(atomTypesArr)
    endif
    
    ! Setup compute for potential energy
    ! call lammps_command(self%lmp, "compute pe_compute all etotal pe" // c_null_char)  ! Invalid: no such compute style 'etotal'; causing segfault in thermo init
    call lammps_command(self%lmp, "thermo_style custom step pe" // c_null_char)
    call lammps_command(self%lmp, "thermo 1" // c_null_char)
    
    ! DEBUG: Dump atom positions to verify they were created correctly
    !call lammps_command(self%lmp, "dump debug_dump all custom 1 lammps_debug.dump id type x y z" // c_null_char)
    !call lammps_command(self%lmp, "run 0 no pre" // c_null_char)
    !call lammps_command(self%lmp, "undump debug_dump" // c_null_char)
    !write(nout, *) "DEBUG: Atom positions written to lammps_debug.dump"
    
    self%initialized = .true.
    write(nout, *) "LAMMPS instance created successfully with", nAtoms, "atoms"
    
  end subroutine
!=============================================================================+
  subroutine UpdatePositions_LAMMPS(self, curbox)
    use CubicBoxDef, only: CubeBox
    use OrthoBoxDef, only: OrthoBox
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    
    integer :: iAtom, nAtoms, idx
    real(dp), pointer :: atoms(:,:) => null()
    real(dp) :: tempdim(1:2, 1:3)
    
    ! Update box dimensions in LAMMPS to match curbox
    select type(curbox)
      class is(CubeBox)
        call curbox%GetDimensions(tempdim)
        self%boxlo(1) = tempdim(1,1) * self%lengthConversion
        self%boxlo(2) = tempdim(1,2) * self%lengthConversion
        self%boxlo(3) = tempdim(1,3) * self%lengthConversion
        self%boxhi(1) = tempdim(2,1) * self%lengthConversion
        self%boxhi(2) = tempdim(2,2) * self%lengthConversion
        self%boxhi(3) = tempdim(2,3) * self%lengthConversion
        call lammps_reset_box(self%lmp, self%boxlo, self%boxhi, 0.0_c_double, 0.0_c_double, 0.0_c_double)
      class is(OrthoBox)
        call curbox%GetDimensions(tempdim)
        self%boxlo(1) = tempdim(1,1) * self%lengthConversion
        self%boxlo(2) = tempdim(1,2) * self%lengthConversion
        self%boxlo(3) = tempdim(1,3) * self%lengthConversion
        self%boxhi(1) = tempdim(2,1) * self%lengthConversion
        self%boxhi(2) = tempdim(2,2) * self%lengthConversion
        self%boxhi(3) = tempdim(2,3) * self%lengthConversion
        call lammps_reset_box(self%lmp, self%boxlo, self%boxhi, 0.0_c_double, 0.0_c_double, 0.0_c_double)
    end select
    
    ! Update atom positions
    call curbox%GetCoordinates(atoms)
    
    ! Count active atoms and prepare position array
    nAtoms = 0
    do iAtom = 1, curbox%nMaxAtoms
      if (curbox%IsActive(iAtom)) then
        nAtoms = nAtoms + 1
        idx = (nAtoms - 1) * 3
        self%tempcoords(idx + 1) = atoms(1, iAtom) * self%lengthConversion
        self%tempcoords(idx + 2) = atoms(2, iAtom) * self%lengthConversion
        self%tempcoords(idx + 3) = atoms(3, iAtom) * self%lengthConversion
      endif
    enddo
    
    ! Scatter positions to LAMMPS
    call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
    
  end subroutine
!=============================================================================+
  function GetEnergy_LAMMPS(self) result(E_Total)
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    real(dp) :: E_Total
    integer(c_int) :: n_atoms
    real(c_double) :: pe_value
    
    ! Run a single step to evaluate energy (run 0 just computes)
    call lammps_command(self%lmp, "run 0 pre no post no" // c_null_char)
    
    ! Get potential energy from custom compute via thermo keyword

    pe_value = lammps_get_thermo(self%lmp, "pe" // c_null_char)
    n_atoms = lammps_get_natoms(self%lmp)
    !write(6, *) "DEBUG GetEnergy: raw pe=", pe_value, " converted=", pe_value * self%energyConversion
    E_Total = pe_value * self%energyConversion * n_atoms
    
  end function
!=============================================================================+
  subroutine DetailedECalc_LAMMPS(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    
    integer :: iAtom, jAtom
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, rmin_ij
    real(dp), pointer :: atoms(:,:) => null()
    
    accept = .true.
    E_T = 0.0_dp
    
    ! Initialize LAMMPS if not done
    if (.not. self%initialized) then
      call self%InitLAMMPS(curbox)
    endif
    
    call curbox%GetCoordinates(atoms)
    
    ! First check for overlapping atoms using rMin criteria
    do iAtom = 1, curbox%nMaxAtoms - 1
      if (.not. curbox%IsActive(iAtom)) cycle
      atmType1 = curbox%AtomType(iAtom)
      
      do jAtom = iAtom + 1, curbox%nMaxAtoms
        if (.not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
        
        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        
        atmType2 = curbox%AtomType(jAtom)
        rmin_ij = self%rMinTable(atmType1, atmType2)
        
        if (rsq < rmin_ij) then
          write(*,*) "ERROR! Overlapping atoms found:", iAtom, jAtom, sqrt(rsq)
          accept = .false.
          return
        endif
      enddo
    enddo
    
    ! Update positions in LAMMPS
    call self%UpdatePositions(curbox)
    
    ! Get energy from LAMMPS
    E_T = self%GetEnergy()
    
    write(nout, *) "Total LAMMPS Energy:", E_T
    
  end subroutine
!=============================================================================+
  subroutine DiffECalc_LAMMPS(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(inout), target :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inout) :: E_Diff
    logical, intent(out) :: accept
    
    integer :: iDisp, iAtom, jAtom, jNei, idx
    integer :: atmType1, atmType2, molStart, molEnd
    real(dp) :: rx, ry, rz, rsq, rmin_ij
    real(dp) :: E_Old, E_New, E_Circle
    real(dp), pointer :: atoms(:,:) => null()
    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    character(len=256) :: cmd
    
    type(Displacement), pointer :: move(:) => null()
    type(Addition), pointer :: add(:) => null()
    type(Deletion), pointer :: del(:) => null()
    type(OrthoVolChange), pointer :: ortho(:) => null()
    
    ! Locals for trial state management in perturbation calculations
    real(c_double), allocatable :: vol_old_coords(:), vol_trial_coords(:)
    real(c_double), allocatable :: disp_old_coords(:)
    real(dp) :: vol_old_boxlo(3), vol_old_boxhi(3)
    real(dp) :: vol_xScale, vol_yScale, vol_zScale
    integer :: vol_iAtom, vol_lmp_idx, vol_molIndx, vol_n_active_atoms
    real(dp) :: vol_com_x, vol_com_y, vol_com_z, vol_shift_x, vol_shift_y, vol_shift_z
    real(dp) :: vol_pos_sim_x, vol_pos_sim_y, vol_pos_sim_z, vol_rel_x, vol_rel_y, vol_rel_z
    
    accept = .true.
    E_Diff = 0.0_dp
    curbox%dETable = 0.0_dp
    
    ! Ensure LAMMPS is synchronized with curbox accepted state before computing trial
    ! (Previous rejected moves might have left LAMMPS in a trial state)
    call self%UpdatePositions(curbox)
    call lammps_command(self%lmp, "run 0" // c_null_char)
    
    call curbox%GetCoordinates(atoms)
    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    
    select type(disp)
      class is(Displacement)
        move => disp
        ! Check rMin for the new positions
        do iDisp = 1, size(move)
          iAtom = move(iDisp)%atmIndx
          atmType1 = curbox%AtomType(iAtom)
          !write(6, *) move(iDisp)%x_new, move(iDisp)%y_new, move(iDisp)%z_new
          !write(6, *) atoms(1, iAtom), atoms(2, iAtom), atoms(3, iAtom)
          
          do jNei = 1, nNeigh(iAtom)
            jAtom = neighlist(jNei, iAtom)
            atmType2 = curbox%AtomType(jAtom)
            rmin_ij = self%rMinTable(atmType1, atmType2)
            
            rx = move(iDisp)%x_new - atoms(1, jAtom)
            ry = move(iDisp)%y_new - atoms(2, jAtom)
            rz = move(iDisp)%z_new - atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            
            if (rsq < rmin_ij) then
              accept = .false.
              return
            endif
          enddo
        enddo
        
        ! Calculate old energy (already stored from previous calculation)
        E_Old = curbox%E_Inter
        
        ! Gather current positions from LAMMPS, modify, then scatter back
        ! This ensures we match LAMMPS's expected data format exactly
        call lammps_gather_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        
        ! Save old coordinates for revert (using predeclared disp_old_coords)
        allocate(disp_old_coords(size(self%tempcoords)))
        disp_old_coords = self%tempcoords
        
        ! Update the gathered positions with new trial positions
        do iDisp = 1, size(move)
          iAtom = move(iDisp)%atmIndx
          ! Find the LAMMPS index for this atom (account for inactive atoms)
          idx = 0
          do jAtom = 1, iAtom
            if (curbox%IsActive(jAtom)) idx = idx + 1
          enddo
          ! Update tempcoords with new position (idx is 1-based, so atom 1 is at indices 1,2,3)
          self%tempcoords((idx-1)*3 + 1) = move(iDisp)%x_new * self%lengthConversion
          self%tempcoords((idx-1)*3 + 2) = move(iDisp)%y_new * self%lengthConversion
          self%tempcoords((idx-1)*3 + 3) = move(iDisp)%z_new * self%lengthConversion
        enddo
        
        ! Scatter modified positions back to LAMMPS
        call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        
        ! Force neighbor list rebuild and energy computation
        call lammps_command(self%lmp, "run 0" // c_null_char)
        
        E_New = self%GetEnergy()
        E_Diff = E_New - E_Old
        !write(6, *) "DEBUG DiffECalc: E_Old=", E_Old, " E_New=", E_New, " E_Diff=", E_Diff
        
        ! Revert LAMMPS positions to old state
        self%tempcoords = disp_old_coords
        call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        call lammps_command(self%lmp, "run 0" // c_null_char)
        deallocate(disp_old_coords)
        
      class is(Addition)
        add => disp
        ! Addition moves - compute trial energy then revert
        E_Old = curbox%E_Inter
        
        ! Gather old state
        call lammps_gather_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        allocate(disp_old_coords(size(self%tempcoords)))
        disp_old_coords = self%tempcoords
        
        ! Apply trial (UpdatePositions includes the added particle from curbox)
        call self%UpdatePositions(curbox)
        call lammps_command(self%lmp, "run 0" // c_null_char)
        E_New = self%GetEnergy()
        E_Diff = E_New - E_Old
        
        ! Revert LAMMPS to old state (without the added particle)
        self%tempcoords = disp_old_coords
        call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        call lammps_command(self%lmp, "run 0" // c_null_char)
        deallocate(disp_old_coords)
        
      class is(Deletion)
        del => disp
        ! Deletion moves - compute trial energy then revert
        E_Old = curbox%E_Inter
        
        ! Gather old state
        call lammps_gather_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        allocate(disp_old_coords(size(self%tempcoords)))
        disp_old_coords = self%tempcoords
        
        ! Apply trial (UpdatePositions excludes the deleted particle from curbox)
        call self%UpdatePositions(curbox)
        call lammps_command(self%lmp, "run 0" // c_null_char)
        E_New = self%GetEnergy()
        E_Diff = E_New - E_Old
        
        ! Revert LAMMPS to old state (with the deleted particle)
        self%tempcoords = disp_old_coords
        call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        call lammps_command(self%lmp, "run 0" // c_null_char)
        deallocate(disp_old_coords)
        
      class is(OrthoVolChange)
        ortho => disp
        ! Handle volume change - update box dimensions and shift positions by COM adjustment in LAMMPS
        E_Old = curbox%E_Inter
        
        ! Save old box bounds (using predeclared vol_old_box*)
        vol_old_boxlo = self%boxlo
        vol_old_boxhi = self%boxhi
        
        ! Compute new scaled box bounds (scaling around coordinate origin 0)
        vol_xScale = ortho(1)%xScale
        vol_yScale = ortho(1)%yScale
        vol_zScale = ortho(1)%zScale
        
        self%boxlo(1) = vol_old_boxlo(1) * vol_xScale
        self%boxlo(2) = vol_old_boxlo(2) * vol_yScale
        self%boxlo(3) = vol_old_boxlo(3) * vol_zScale
        self%boxhi(1) = vol_old_boxhi(1) * vol_xScale
        self%boxhi(2) = vol_old_boxhi(2) * vol_yScale
        self%boxhi(3) = vol_old_boxhi(3) * vol_zScale
        
        ! Gather current (old) positions from LAMMPS
        call lammps_gather_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        
        ! Allocate and save old/trial coords (using predeclared vol_*_coords)
        allocate(vol_old_coords(size(self%tempcoords)), vol_trial_coords(size(self%tempcoords)))
        vol_old_coords = self%tempcoords

        !Check Boundary Conditions for the old positions 


        vol_trial_coords = self%tempcoords
        
        ! Apply molecule COM shifts to trial positions with re-centering (matching FF_EasyPair_Cut and box logic; using predeclared vol_* vars)
        vol_n_active_atoms = 0
        do vol_iAtom = 1, curbox%nMaxAtoms
          if (curbox%IsActive(vol_iAtom)) then
            vol_n_active_atoms = vol_n_active_atoms + 1
            vol_lmp_idx = (vol_n_active_atoms - 1) * 3
            vol_molIndx = curbox%MolIndx(vol_iAtom)
            
            ! Convert current position to sim units for re-centering
            vol_pos_sim_x = vol_trial_coords(vol_lmp_idx + 1) / self%lengthConversion
            vol_pos_sim_y = vol_trial_coords(vol_lmp_idx + 2) / self%lengthConversion
            vol_pos_sim_z = vol_trial_coords(vol_lmp_idx + 3) / self%lengthConversion
            
            ! Re-center relative to mol COM (subtract COM, boundary rel pos, add back COM)
            vol_rel_x = vol_pos_sim_x - curbox%centerMass(1, vol_molIndx)
            vol_rel_y = vol_pos_sim_y - curbox%centerMass(2, vol_molIndx)
            vol_rel_z = vol_pos_sim_z - curbox%centerMass(3, vol_molIndx)
            call curbox%Boundary(vol_rel_x, vol_rel_y, vol_rel_z)
            vol_pos_sim_x = vol_rel_x + curbox%centerMass(1, vol_molIndx)
            vol_pos_sim_y = vol_rel_y + curbox%centerMass(2, vol_molIndx)
            vol_pos_sim_z = vol_rel_z + curbox%centerMass(3, vol_molIndx)
            
            ! Compute and add COM shift for volume scaling
            vol_com_x = curbox%centerMass(1, vol_molIndx)
            vol_com_y = curbox%centerMass(2, vol_molIndx)
            vol_com_z = curbox%centerMass(3, vol_molIndx)
            vol_shift_x = vol_com_x * (vol_xScale - 1.0_dp)
            vol_shift_y = vol_com_y * (vol_yScale - 1.0_dp)
            vol_shift_z = vol_com_z * (vol_zScale - 1.0_dp)
            vol_pos_sim_x = vol_pos_sim_x + vol_shift_x
            vol_pos_sim_y = vol_pos_sim_y + vol_shift_y
            vol_pos_sim_z = vol_pos_sim_z + vol_shift_z
            
            ! Convert back to LAMMPS units and update trial coords
            vol_trial_coords(vol_lmp_idx + 1) = vol_pos_sim_x * self%lengthConversion
            vol_trial_coords(vol_lmp_idx + 2) = vol_pos_sim_y * self%lengthConversion
            vol_trial_coords(vol_lmp_idx + 3) = vol_pos_sim_z * self%lengthConversion
          end if
        end do
        
        ! Reset box to trial dimensions
        call lammps_reset_box(self%lmp, self%boxlo, self%boxhi, 0.0_c_double, 0.0_c_double, 0.0_c_double)
        ! Apply trial positions to LAMMPS
        self%tempcoords = vol_trial_coords
        call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        
        
        ! Update LAMMPS internal state
        call lammps_command(self%lmp, "run 0" // c_null_char)
        
        ! Get trial energy
        E_New = self%GetEnergy()
        E_Diff = E_New - E_Old
        
        ! Revert LAMMPS to old state
        self%tempcoords = vol_old_coords
        call lammps_scatter_atoms(self%lmp, "x" // c_null_char, 1_c_int, 3_c_int, self%tempcoords)
        self%boxlo = vol_old_boxlo
        self%boxhi = vol_old_boxhi
        call lammps_reset_box(self%lmp, self%boxlo, self%boxhi, 0.0_c_double, 0.0_c_double, 0.0_c_double)
        call lammps_command(self%lmp, "run 0" // c_null_char)
        
        E_Circle = self%GetEnergy()
        write(6, *) "DEBUG DiffECalc: E_Old=", E_Old, " E_New=", E_New, " E_Diff=", E_Diff, " E_Circle=", E_Circle
        deallocate(vol_old_coords, vol_trial_coords)
        
      class default
        write(0,*) "Unknown Perturbation Type Encountered by LAMMPS Pair Style."
        error stop
        
    end select
    
  end subroutine
!=============================================================================+
  subroutine ProcessIO_LAMMPS(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Units, only: inEngUnit, inLenUnit
    use ParallelVar, only: nout
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=100) :: command
    character(len=maxLineLen) :: restOfLine
    integer :: lineStat, nPar, type1, type2, jType
    real(dp) :: rMin, rCut
    logical :: param = .false.
    
    call GetXCommand(line, command, 1, lineStat)
    
    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self%rCut = rCut * inLenUnit
        self%rCutSq = self%rCut * self%rCut
        write(nout, *) "LAMMPS cutoff set to:", self%rCut
        
      case("lammps_cmd")
        ! Add a LAMMPS command to be executed during initialization
        self%nSetupCommands = self%nSetupCommands + 1
        ! Extract everything after "lammps_cmd" by finding its position in the line
        type1 = index(line, "lammps_cmd")
        if (type1 > 0) then
          restOfLine = adjustl(line(type1 + 10:))  ! 10 = len("lammps_cmd")
        else
          restOfLine = ""
        endif
        self%setupCommands(self%nSetupCommands) = trim(restOfLine)
        write(nout, *) "Added LAMMPS command: ", trim(restOfLine)
        
      case("inputscript")
        ! Specify a LAMMPS input script file
        call GetXCommand(line, command, 2, lineStat)
        self%inputScript = trim(adjustl(command))
        write(nout, *) "LAMMPS input script: ", trim(self%inputScript)
        
      case("energy_conversion")
        ! Conversion factor from LAMMPS energy units to Classy units
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) self%energyConversion
        write(nout, *) "Energy conversion factor:", self%energyConversion
        
      case("length_conversion")
        ! Conversion factor from Classy length units to LAMMPS units
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) self%lengthConversion
        write(nout, *) "Length conversion factor:", self%lengthConversion
        
      case("typemap")
        ! Map Classy atom type to LAMMPS atom type
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        self%typeMap(type1) = type2
        write(nout, *) "Type mapping: Classy", type1, "-> LAMMPS", type2
        
      case default
        param = .true.
    end select
    
    ! Handle rMin parameters (similar to other forcefields)
    if (param) then
      call CountCommands(line, nPar)
      select case(nPar)
        case(2)
          read(line, *) type1, rMin
          rMin = rMin * inLenUnit
          self%rMin(type1) = rMin
          do jType = 1, nAtomTypes
            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
          
        case(3)
          read(line, *) type1, type2, rMin
          rMin = rMin * inLenUnit
          self%rMinTable(type1, type2) = rMin**2
          self%rMinTable(type2, type1) = rMin**2
          
        case default
          lineStat = -1
      end select
    endif
    
  end subroutine
!=============================================================================+
  subroutine Prologue_LAMMPS(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    integer :: atomLimit, iBox
    integer :: AllocateStat
    
    ! Determine maximum number of atoms
    atomLimit = 0
    do iBox = 1, size(BoxArray)
      if (atomLimit < boxArray(iBox)%box%GetMaxAtoms()) then
        atomLimit = boxArray(iBox)%box%GetMaxAtoms()
      endif
    enddo
    
    ! Allocate temporary coordinate storage (x, y, z for each atom)
    allocate(self%tempcoords(3 * atomLimit), stat = AllocateStat)
    allocate(self%atomTypes(atomLimit), stat = AllocateStat)
    
    if (AllocateStat /= 0) then
      error stop "Memory allocation error in LAMMPS Prologue"
    endif
    
    write(nout, *) "LAMMPS Forcefield Prologue Complete"
    write(nout, *) "  Max atoms:", atomLimit
    write(nout, *) "  Cutoff distance:", self%rCut
    
  end subroutine
!=============================================================================+
  subroutine Epilogue_LAMMPS(self)
    use ParallelVar, only: nout
    implicit none
    class(Pair_LAMMPS), intent(inout) :: self
    
    ! Close LAMMPS instance
    if (self%lammps_created .and. c_associated(self%lmp)) then
      write(nout, *) "Closing LAMMPS instance..."
      call lammps_close(self%lmp)
      self%lmp = c_null_ptr
      self%lammps_created = .false.
      self%initialized = .false.
    endif
    
    ! Deallocate arrays
    if (allocated(self%tempcoords)) deallocate(self%tempcoords)
    if (allocated(self%atomTypes)) deallocate(self%atomTypes)
    if (allocated(self%setupCommands)) deallocate(self%setupCommands)
    if (allocated(self%typeMap)) deallocate(self%typeMap)
    
    write(nout, *) "LAMMPS Forcefield Epilogue Complete"
    
  end subroutine
!=============================================================================+
#endif
! End LAMMPS Safety Block
!=============================================================================+
end module FF_Ext_LAMMPS
!=============================================================================+

