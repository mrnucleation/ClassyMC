!=========================================================================
!  Energy Module that is designed to call a Python Script 
!  in order to compute energies for Monte Carlo simulations.
!  This module is designed to pass Classy style
!  object information into Python.  
!
!  For a valid Energy script the following functions must be defined:
!
!  def prologue(boxlist):
!      """Called once at the start of the simulation to initialize.
!      
!      Args:
!          boxlist: List of box dictionaries containing simulation state
!      """
!      pass
!
!  def compute_total(boxlist):
!      """Compute the total energy of the system.
!      
!      Args:
!          boxlist: List of box dictionaries containing:
!              - 'boxtype': str - type of box ('cube', 'ortho', etc.)
!              - 'temperature': float - temperature
!              - 'pressure': float - pressure
!              - 'volume': float - volume
!              - 'boxdimensions': numpy array (2 x ndim) - box bounds
!              - 'atomtype': numpy array - atom types (1-indexed)
!              - 'raw_atoms': numpy array (3 x nMaxAtoms) - atomic coordinates
!              - 'moleculecount': numpy array - molecules per type
!              - 'moltype': numpy array - molecule type per atom
!              - 'molsubindx': numpy array - molecule sub-index per atom
!              - 'energytable': numpy array - per-atom energies
!              - 'neighlists': list of dicts (one per neighbor list), each containing:
!                  - 'list': numpy array (maxNei x nMaxAtoms) - neighbor indices per atom
!                  - 'nneigh': numpy array (nMaxAtoms) - number of neighbors per atom
!                  Note: For atom i, valid neighbors are list[0:nneigh[i], i]
!                        (Fortran column-major order, so atom index is the 2nd axis)
!      
!      Returns:
!          dict: Dictionary containing:
!              - 'energy': float - total system energy
!              - 'accept': bool - True if configuration is valid (optional, default True)
!              - 'energytable': numpy array - per-atom energies (optional)
!      """
!      return {'energy': 0.0, 'accept': True}
!
!  def compute_diff(boxlist, displist):
!      """Compute the energy difference for a proposed move.
!      
!      Args:
!          boxlist: List of box dictionaries (current state)
!          displist: List of displacement dictionaries, each containing:
!              - 'type': str - 'displacement', 'addition', 'deletion', etc.
!              - 'moltype': int - molecule type
!              - 'molindex': int - molecule index
!              - 'atomindex': int - atom index
!              - 'x_new', 'y_new', 'z_new': float - new coordinates (if applicable)
!              - For 'newstate_isomol':
!                    'new_atoms': numpy array (3 x nMaxAtoms)
!      
!      Returns:
!          dict: Dictionary containing:
!              - 'energy': float - energy difference (E_new - E_old)
!              - 'accept': bool - True if move is valid (no overlaps, etc.)
!              - 'denergytable': numpy array - per-atom energy differences (optional)
!      """
!      return {'energy': 0.0, 'accept': True}
!
!  def compute_forces(boxlist, coords):
!      """Optional. Forces at the given coordinates (current or trial state).
!      Needed for native force-bias / smart MC. If omitted, Fortran falls
!      back to finite-difference of compute_diff on the current box.
!
!      Args:
!          boxlist: List of box dictionaries (types, PBC, occupancy).
!                   Do not use box['raw_atoms'] as the configuration;
!                   that is the current box, not necessarily coords.
!          coords: numpy array (3 x nMaxAtoms), Fortran order — the
!                  configuration to evaluate, which may be a trial state.
!
!      Returns:
!          dict with 'forces': numpy array (3 x nMaxAtoms), Fortran order.
!          F = -dU/dr. Inactive atoms should be zero.
!      """
!      return {'forces': numpy.zeros_like(coords)}
!
!  def update(boxlist):
!      """Called after an accepted move to update internal state.
!      
!      Args:
!          boxlist: List of box dictionaries (now updated)
!      """
!      pass
!
!=========================================================================
#define errcheck_macro if(ierror/=0) then;call err_print;stop;endif
!=========================================================================
module FF_PythonEnergy
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use Input_Format, only: maxLineLen
  use VarPrecision
  use CoordinateTypes

#ifdef EMBPYTHON
  use forpy_mod, only: dict, dict_create, list, list_create, call_py, &
                       module_py, import_py, object, call_py_noret, tuple, &
                       tuple_create, ndarray, ndarray_create, ndarray_create_nocopy, &
                       cast, err_print, err_clear, get_sys_path, is_none
!------------------------------------------------------------------------------
  type, public, extends(forcefield) :: Pair_PythonEnergy
    logical, private :: initialized = .false.
    logical, private :: hasPyForces = .false.
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    type(module_py) :: pyenergy
    type(list) :: boxlist
    type(tuple) :: args
    type(tuple) :: diff_args
    type(dict), allocatable :: boxdicts(:)
    character(len=100) :: pymodule
    contains
      procedure, pass :: Constructor => Constructor_PythonEnergy
      procedure, pass :: DetailedECalc => DetailedECalc_PythonEnergy
      procedure, pass :: DiffECalc => DiffECalc_PythonEnergy
      procedure, pass :: HasComputeForces => HasComputeForces_PythonEnergy
      procedure, pass :: ComputeForces => ComputeForces_PythonEnergy
      procedure, pass :: ProcessIO => ProcessIO_PythonEnergy
      procedure, pass :: Prologue => Prologue_PythonEnergy
      procedure, pass :: Epilogue => Epilogue_PythonEnergy
      procedure, pass :: GetCutOff => GetCutOff_PythonEnergy
      procedure, pass :: Update => Update_PythonEnergy
  end type

 contains
!=========================================================================
  subroutine Constructor_PythonEnergy(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 0E0_dp
    self%rMinTable = 0E0_dp
    self%rCut = 1000E0_dp  ! Large default cutoff
    self%rCutSq = self%rCut * self%rCut

  end subroutine
!=========================================================================
  subroutine Prologue_PythonEnergy(self)
    use BoxData, only: BoxArray
    use ClassyPyObj, only: createboxdict_nocopy
    use ParallelVar, only: nout
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    integer :: ierror, iBox, nBoxes
    type(list) :: paths
    type(object) :: forceattr

    if(self%initialized) then
      return
    endif

    nBoxes = size(BoxArray)
    allocate(self%boxdicts(1:nBoxes))

    ierror = get_sys_path(paths)
    errcheck_macro
    ierror = paths%append(".")
    errcheck_macro
    call paths%destroy

    write(nout,"(1x,A,A)") "Loading Python Energy Module: ", trim(adjustl(self%pymodule))
    ierror = import_py(self%pyenergy, trim(adjustl(self%pymodule)))
    errcheck_macro

    self%hasPyForces = .false.
    ierror = self%pyenergy%getattribute(forceattr, "compute_forces")
    if(ierror == 0) then
      self%hasPyForces = .true.
      call forceattr%destroy
    else
      call err_clear()
    endif

    ierror = list_create(self%boxlist)
    errcheck_macro
    do iBox = 1, nBoxes
      self%boxdicts(iBox) = createboxdict_nocopy(iBox)
      ierror = self%boxlist%append(self%boxdicts(iBox))
      errcheck_macro
    enddo

    ierror = tuple_create(self%args, 1)
    errcheck_macro
    ierror = self%args%setitem(0, self%boxlist)
    errcheck_macro

    ierror = tuple_create(self%diff_args, 4)
    errcheck_macro
    ierror = self%diff_args%setitem(0, self%boxlist)
    errcheck_macro

    ierror = call_py_noret(self%pyenergy, "prologue", args=self%args)
    errcheck_macro

    self%initialized = .true.

    write(nout,"(1x,A,A)") "(Python Energy) Loaded module: ", trim(adjustl(self%pymodule))
    if(self%hasPyForces) then
      write(nout,"(1x,A)") "(Python Energy) compute_forces found; native force-bias will use Python forces"
    else
      write(nout,"(1x,A)") "(Python Energy) no compute_forces; forces fall back to finite difference"
    endif

  end subroutine
!=========================================================================
  subroutine Epilogue_PythonEnergy(self)
    use ClassyPyObj, only: refreshboxdicts
    use ParallelVar, only: nout
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    integer :: ierror, iBox

    if(.not. self%initialized) then
      return
    endif

    call refreshboxdicts(self%boxdicts)
    ierror = call_py_noret(self%pyenergy, "epilogue", args=self%args)
    errcheck_macro

    call self%args%destroy
    call self%diff_args%destroy
    call self%boxlist%destroy
    do iBox = 1, size(self%boxdicts)
      call self%boxdicts(iBox)%destroy
    enddo
    deallocate(self%boxdicts)
    call self%pyenergy%destroy
    self%initialized = .false.

    write(nout,"(1x,A)") "(Python Energy) Epilogue completed"
  end subroutine

!=========================================================================
  subroutine DetailedECalc_PythonEnergy(self, curbox, E_T, accept)
    use ClassyPyObj, only: refreshboxdicts
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    
    type(object) :: returnobj, item
    type(dict) :: resultdict
    type(ndarray) :: np_etable
    integer :: ierror
    real(dp) :: energy
    real(dp), pointer :: etable_ptr(:) => null()

    accept = .true.
    E_T = 0E0_dp

    call refreshboxdicts(self%boxdicts)

    ierror = call_py(returnobj, self%pyenergy, "compute_total", args=self%args)
    errcheck_macro

    ierror = cast(resultdict, returnobj)
    errcheck_macro

    ierror = resultdict%getitem(item, "energy")
    if(ierror /= 0) then
      write(0,*) "ERROR! Python compute_total must return dict with 'energy' key"
      error stop
    endif
    ierror = cast(energy, item)
    errcheck_macro
    call item%destroy
    E_T = energy

    ierror = resultdict%getitem(item, "accept")
    if(ierror == 0) then
      ierror = cast(accept, item)
      call item%destroy
    endif

    ierror = resultdict%getitem(item, "energytable")
    if(ierror == 0) then
      ierror = cast(np_etable, item)
      if(ierror == 0) then
        ierror = np_etable%get_data(etable_ptr)
        if(ierror == 0 .and. associated(etable_ptr)) then
          if(size(etable_ptr) >= curbox%nMaxAtoms) then
            curbox%ETable(1:curbox%nMaxAtoms) = etable_ptr(1:curbox%nMaxAtoms)
          endif
        endif
        call np_etable%destroy
      endif
      call item%destroy
    endif

    call resultdict%destroy
    call returnobj%destroy

    if(isnan(E_T)) then
      write(0,*) "ERROR! NaN energy returned from Python compute_total"
      error stop
    endif
  end subroutine
!=========================================================================
  subroutine DiffECalc_PythonEnergy(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    use ClassyPyObj, only: createdisplist, refreshboxdicts
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(inout), target :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    type(object) :: returnobj, item
    type(dict) :: resultdict
    type(list) :: displist
    type(ndarray) :: np_detable
    type(ndarray) :: np_nlist, np_nneigh
    integer :: ierror
    real(dp) :: energy
    real(dp), pointer :: detable_ptr(:) => null()

    accept = .true.
    E_Diff = 0E0_dp
    curbox%dETable = 0E0_dp

    call refreshboxdicts(self%boxdicts)

    displist = createdisplist(disp)
    ierror = self%diff_args%setitem(1, displist)
    errcheck_macro
    call displist%destroy

    ierror = ndarray_create_nocopy(np_nlist, tempList)
    errcheck_macro
    ierror = self%diff_args%setitem(2, np_nlist)
    errcheck_macro
    call np_nlist%destroy

    ierror = ndarray_create_nocopy(np_nneigh, tempNNei)
    errcheck_macro
    ierror = self%diff_args%setitem(3, np_nneigh)
    errcheck_macro
    call np_nneigh%destroy

    ierror = call_py(returnobj, self%pyenergy, "compute_diff", args=self%diff_args)
    errcheck_macro

    ierror = cast(resultdict, returnobj)
    errcheck_macro

    ierror = resultdict%getitem(item, "energy")
    if(ierror /= 0) then
      write(0,*) "ERROR! Python compute_diff must return dict with 'energy' key"
      error stop
    endif
    ierror = cast(energy, item)
    errcheck_macro
    call item%destroy
    E_Diff = energy

    ierror = resultdict%getitem(item, "accept")
    if(ierror == 0) then
      ierror = cast(accept, item)
      call item%destroy
    endif

    ierror = resultdict%getitem(item, "denergytable")
    if(ierror == 0) then
      ierror = cast(np_detable, item)
      if(ierror == 0) then
        ierror = np_detable%get_data(detable_ptr)
        if(ierror == 0 .and. associated(detable_ptr)) then
          if(size(detable_ptr) >= curbox%nMaxAtoms) then
            curbox%dETable(1:curbox%nMaxAtoms) = detable_ptr(1:curbox%nMaxAtoms)
          endif
        endif
        call np_detable%destroy
      endif
      call item%destroy
    endif

    call resultdict%destroy
    call returnobj%destroy

    if(isnan(E_Diff)) then
      write(0,*) "ERROR! NaN energy returned from Python compute_diff"
      error stop
    endif
  end subroutine
!=========================================================================
  subroutine ProcessIO_PythonEnergy(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use ParallelVar, only: nout
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=100) :: command
    logical :: param = .false.
    integer :: jType, lineStat, nPar
    integer :: type1, type2
    real(dp) :: rCut, rMin

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("module")
        call GetXCommand(line, command, 2, lineStat)
        self%pymodule = trim(adjustl(command))
        write(nout,"(1x,A,A)") "(Python Energy) Module set to: ", trim(self%pymodule)

      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self%rCut = rCut
        self%rCutSq = rCut * rCut
        write(nout,"(1x,A,F12.4)") "(Python Energy) Cutoff set to: ", rCut

      case("rmin")
        call CountCommands(line, nPar)
        select case(nPar)
          case(2)
            ! rmin type1 value
            call GetXCommand(line, command, 2, lineStat)
            read(command, *) type1
            call GetXCommand(line, command, 3, lineStat)
            read(command, *) rMin
            self%rMin(type1) = rMin
            do jType = 1, nAtomTypes
              self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
              self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
            enddo
          case(3)
            ! rmin type1 type2 value
            call GetXCommand(line, command, 2, lineStat)
            read(command, *) type1
            call GetXCommand(line, command, 3, lineStat)
            read(command, *) type2
            call GetXCommand(line, command, 4, lineStat)
            read(command, *) rMin
            self%rMinTable(type1, type2) = rMin**2
            self%rMinTable(type2, type1) = rMin**2
          case default
            write(0,*) "ERROR! Invalid rmin format. Use: rmin type value OR rmin type1 type2 value"
        end select

      case default
        write(0,*) "ERROR! Unknown Python energy command: ", trim(adjustl(command))
        write(0,*) "Valid commands: module, rcut, rmin"
    end select

  end subroutine
!=========================================================================
  function GetCutOff_PythonEnergy(self) result(rCut)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
!=========================================================================
  subroutine Update_PythonEnergy(self, accept)
    use ClassyPyObj, only: refreshboxdicts
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    logical, intent(in) :: accept

    type(object) :: returnobj
    integer :: ierror

    if(.not. self%initialized) then
      return
    endif

    call refreshboxdicts(self%boxdicts)
    ierror = call_py(returnobj, self%pyenergy, "update", args=self%args)
    errcheck_macro

    call returnobj%destroy
  end subroutine
!=========================================================================
  function HasComputeForces_PythonEnergy(self) result(hasForce)
    implicit none
    class(Pair_PythonEnergy), intent(in) :: self
    logical :: hasForce

    hasForce = self%hasPyForces
  end function
!=========================================================================
  subroutine ComputeForces_PythonEnergy(self, curbox, forces, coords)
    use ClassyPyObj, only: refreshboxdicts
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(out) :: forces(:,:)
    real(dp), intent(in) :: coords(:,:)

    type(object) :: returnobj, item
    type(dict) :: resultdict
    type(tuple) :: force_args
    type(ndarray) :: np_coords, np_forces
    real(dp), pointer :: fptr(:,:) => null()
    integer :: ierror, nCopy

    forces = 0E0_dp
    if(.not. self%initialized) then
      call self%Prologue()
    endif
    if(.not. self%hasPyForces) then
      call self%ComputeForcesFiniteDiff(curbox, forces, coords)
      return
    endif

    call refreshboxdicts(self%boxdicts)

    ierror = tuple_create(force_args, 2)
    errcheck_macro
    ierror = force_args%setitem(0, self%boxlist)
    errcheck_macro
    ierror = ndarray_create(np_coords, coords)
    errcheck_macro
    ierror = force_args%setitem(1, np_coords)
    errcheck_macro
    call np_coords%destroy

    ierror = call_py(returnobj, self%pyenergy, "compute_forces", args=force_args)
    errcheck_macro
    call force_args%destroy

    ierror = cast(resultdict, returnobj)
    if(ierror /= 0) then
      write(0,*) "ERROR! Python compute_forces must return a dict with a 'forces' key"
      error stop
    endif

    ierror = resultdict%getitem(item, "forces")
    if(ierror /= 0) then
      write(0,*) "ERROR! Python compute_forces must return dict with 'forces' key"
      write(0,*) "  'forces' must be a numpy array of shape (3, nAtoms), Fortran order"
      error stop
    endif
    ierror = cast(np_forces, item)
    errcheck_macro
    call item%destroy

    ierror = np_forces%get_data(fptr)
    if(ierror /= 0 .or. .not. associated(fptr)) then
      write(0,*) "ERROR! Python compute_forces 'forces' array could not be read."
      write(0,*) "  Use a contiguous (3, nAtoms) float64 array with order='F'."
      error stop
    endif
    if(size(fptr, 1) /= 3) then
      write(0,*) "ERROR! Python compute_forces 'forces' must have shape (3, nAtoms), got", &
        size(fptr, 1), size(fptr, 2)
      error stop
    endif
    nCopy = min(size(forces, 2), size(fptr, 2))
    forces(1:3, 1:nCopy) = fptr(1:3, 1:nCopy)

    call np_forces%destroy
    call resultdict%destroy
    call returnobj%destroy
  end subroutine
!=========================================================================
#else
! Stub type when EMBPYTHON is not defined
  type, public, extends(forcefield) :: Pair_PythonEnergy
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    contains
      procedure, pass :: Constructor => Constructor_PythonEnergy_Stub
      procedure, pass :: DetailedECalc => DetailedECalc_PythonEnergy_Stub
      procedure, pass :: DiffECalc => DiffECalc_PythonEnergy_Stub
      procedure, pass :: ProcessIO => ProcessIO_PythonEnergy_Stub
      procedure, pass :: Prologue => Prologue_PythonEnergy_Stub
      procedure, pass :: Epilogue => Epilogue_PythonEnergy_Stub
      procedure, pass :: GetCutOff => GetCutOff_PythonEnergy_Stub
  end type

 contains
!=========================================================================
  subroutine Constructor_PythonEnergy_Stub(self)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self

    write(0,*) "ERROR! Python energy requires compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  subroutine DetailedECalc_PythonEnergy_Stub(self, curbox, E_T, accept)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    write(0,*) "ERROR! Python energy requires compilation with -DEMBPYTHON flag"
    error stop
    accept = .false.
  end subroutine
!=========================================================================
  subroutine DiffECalc_PythonEnergy_Stub(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(inout), target :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    write(0,*) "ERROR! Python energy requires compilation with -DEMBPYTHON flag"
    error stop
    accept = .false.
  end subroutine
!=========================================================================
  subroutine ProcessIO_PythonEnergy_Stub(self, line)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

    write(0,*) "ERROR! Python energy requires compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  subroutine Prologue_PythonEnergy_Stub(self)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self

    write(0,*) "ERROR! Python energy requires compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  function GetCutOff_PythonEnergy_Stub(self) result(rCut)
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    real(dp) :: rCut

    write(0,*) "ERROR! Python energy requires compilation with -DEMBPYTHON flag"
    error stop
    rCut = 0E0_dp
  end function
!=========================================================================
#endif
end module
!=========================================================================
