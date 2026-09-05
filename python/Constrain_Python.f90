!=========================================================================
!  Constraint Module that is designed to call a Python Script 
!  in order to enforce constraints on Monte Carlo moves.
!  This module is designed to pass Classy style
!  object information into Python.  
!
!  For a valid Constraint script the following functions must be defined:
!
!  def prologue(boxinfo):
!      """Called once at the start of the simulation to initialize.
!      
!      Args:
!          boxinfo: Dictionary containing box information (see createboxdict)
!      """
!      pass
!
!  def check_initial(boxinfo):
!      """Check if the initial configuration satisfies the constraint.
!      
!      Args:
!          boxinfo: Standardized box dictionary from ClassyPyObj containing:
!              - 'boxtype': str - type of box ('cube', 'ortho', 'nobox', etc.)
!              - 'thread_id': int - parallel thread id
!              - 'box_id': int - box identifier
!              - 'energy': float - total energy
!              - 'temperature': float - temperature
!              - 'pressure': float - pressure
!              - 'volume': float - volume
!              - 'boxdimensions': numpy array (2 x ndim) - box bounds
!              - 'chemicalpotential': numpy array - chemical potentials
!              - 'atomtype': numpy array - atom types (1-indexed)
!              - 'raw_atoms': numpy array (ndim x nMaxAtoms) - atomic coordinates
!              - 'moleculecount': numpy array - molecules per type
!              - 'energytable': numpy array - per-atom energies
!              - 'moltype': numpy array - molecule type per atom
!              - 'molsubindx': numpy array - molecule sub-index per atom
!              - 'molindx': numpy array - molecule index per atom
!              - 'neighlists': list of neighbor-list dicts (if present)
!      
!      Returns:
!          bool: True if constraint satisfied, False otherwise
!      """
!      return True
!
!  def diff_check(boxinfo, displist):
!      """Check if a proposed move satisfies the constraint.
!      
!      Args:
!          boxinfo: Dictionary with current box state
!          displist: List of displacement dictionaries, each containing:
!              - 'type': str - 'displacement', 'addition', 'deletion', etc.
!              - 'moltype': int - molecule type
!              - 'molindex': int - molecule index
!              - 'atomindex': int - atom index
!              - 'x_new', 'y_new', 'z_new': float - new coordinates (if applicable)
!      
!      Returns:
!          bool: True if constraint satisfied, False otherwise
!      """
!      return True
!
!  def post_energy(boxinfo, displist, energy_diff):
!      """Check constraint after energy calculation (optional).
!      
!      Args:
!          boxinfo: Dictionary with current box state
!          displist: List of displacement dictionaries
!          energy_diff: float - energy difference from the move
!      
!      Returns:
!          bool: True if constraint satisfied, False otherwise
!      """
!      return True
!
!  def update(boxinfo):
!      """Called after a move is accepted to update internal state.
!      
!      Args:
!          boxinfo: Dictionary with updated box state
!      """
!      pass
!
!  def epilogue(boxinfo):
!      """Called at the end of the simulation for cleanup.
!      
!      Args:
!          boxinfo: Dictionary with final box state
!      """
!      pass
!
!=========================================================================
#define errcheck_macro if(ierror/=0) then;call err_print;stop;endif
!=========================================================================
module Constrain_Python
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout

#ifdef EMBPYTHON
  use forpy_mod, only: dict, list, call_py, module_py, import_py, object, &
                       call_py_noret, tuple, tuple_create, ndarray, &
                       ndarray_create, cast, err_print, get_sys_path
!------------------------------------------------------------------------------
  type, public, extends(constraint) :: PythonConstraint
    logical, private :: initialized = .false.
    integer :: boxID = 1
    type(module_py) :: pyconstraint
    type(dict) :: boxinfo
    type(tuple) :: args_init
    type(tuple) :: args_diff
    type(tuple) :: args_post
    type(tuple) :: args_update
    character(len=100) :: pymodule
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => PythonConstraint_Constructor
      procedure, pass :: CheckInitialConstraint => PythonConstraint_CheckInitialConstraint
      procedure, pass :: DiffCheck => PythonConstraint_DiffCheck
      procedure, pass :: PostEnergy => PythonConstraint_PostEnergy
      procedure, pass :: ProcessIO => PythonConstraint_ProcessIO
      procedure, pass :: Maintenance => PythonConstraint_Maintenance
      procedure, pass :: Update => PythonConstraint_Update
      procedure, pass :: Epilogue => PythonConstraint_Epilogue
  end type

 contains
!=========================================================================
  subroutine PythonConstraint_Constructor(self, boxID)
    use BoxData, only: BoxArray
    use ClassyPyObj, only: createboxdict_nocopy
    implicit none
    class(PythonConstraint), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: ierror
    type(list) :: paths

    if(self%initialized) then
      return
    endif

    self%boxID = boxID
    self%parent => BoxArray(boxID)%box

    ! Load Python module
    write(nout,"(1x,A,A)") "Loading Python Constraint Module: ", trim(adjustl(self%pymodule))
    ierror = get_sys_path(paths)
    errcheck_macro
    ierror = paths%append(".")
    errcheck_macro
    call paths%destroy

    ierror = import_py(self%pyconstraint, trim(adjustl(self%pymodule)))
    errcheck_macro

    ! Use the standardized ClassyPyObj box dictionary
    self%boxinfo = createboxdict_nocopy(boxID)

    ! Create argument tuples
    ierror = tuple_create(self%args_init, 1)
    errcheck_macro
    ierror = self%args_init%setitem(0, self%boxinfo)
    errcheck_macro

    ierror = tuple_create(self%args_diff, 2)
    errcheck_macro
    ierror = self%args_diff%setitem(0, self%boxinfo)
    errcheck_macro

    ierror = tuple_create(self%args_post, 3)
    errcheck_macro
    ierror = self%args_post%setitem(0, self%boxinfo)
    errcheck_macro

    ierror = tuple_create(self%args_update, 1)
    errcheck_macro
    ierror = self%args_update%setitem(0, self%boxinfo)
    errcheck_macro

    ! Call Python prologue
    ierror = call_py_noret(self%pyconstraint, "prologue", args=self%args_init)
    errcheck_macro

    self%initialized = .true.

    write(nout,"(1x,A,A)") "(Python Constraint) Loaded module: ", trim(adjustl(self%pymodule))

  end subroutine
!=========================================================================
  subroutine PythonConstraint_CheckInitialConstraint(self, trialBox, accept)
    use ClassyPyObj, only: refreshboxdict
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    type(object) :: returnobj
    integer :: ierror

    accept = .true.

    call refreshboxdict(self%boxinfo, self%boxID)

    ierror = call_py(returnobj, self%pyconstraint, "check_initial", args=self%args_init)
    errcheck_macro
    ierror = cast(accept, returnobj)
    errcheck_macro
    call returnobj%destroy

    if(accept) then
      write(nout,"(1x,A)") "(Python Constraint) Initial check: PASSED"
    else
      write(nout,"(1x,A)") "(Python Constraint) Initial check: FAILED"
    endif

  end subroutine
!=========================================================================
  subroutine PythonConstraint_DiffCheck(self, trialBox, disp, accept)
    use ClassyPyObj, only: createdisplist
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    type(object) :: returnobj
    type(list) :: displist
    integer :: ierror

    accept = .true.

    ! Create displacement list for Python
    displist = createdisplist(disp)

    ! Set displacement list in args
    ierror = self%args_diff%setitem(1, displist)
    errcheck_macro

    ! Call Python diff_check
    ierror = call_py(returnobj, self%pyconstraint, "diff_check", args=self%args_diff)
    errcheck_macro
    ierror = cast(accept, returnobj)
    errcheck_macro
    call returnobj%destroy
    call displist%destroy

  end subroutine
!=========================================================================
  subroutine PythonConstraint_PostEnergy(self, trialBox, disp, E_Diff, accept)
    use ClassyPyObj, only: createdisplist
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff
    logical, intent(out) :: accept
    type(object) :: returnobj
    type(list) :: displist
    integer :: ierror

    accept = .true.

    ! Create displacement list for Python
    displist = createdisplist(disp)

    ! Set displacement list and energy in args
    ierror = self%args_post%setitem(1, displist)
    errcheck_macro
    ierror = self%args_post%setitem(2, E_Diff)
    errcheck_macro

    ! Call Python post_energy
    ierror = call_py(returnobj, self%pyconstraint, "post_energy", args=self%args_post)
    errcheck_macro
    ierror = cast(accept, returnobj)
    errcheck_macro
    call returnobj%destroy
    call displist%destroy

  end subroutine
!=========================================================================
  subroutine PythonConstraint_Update(self, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    logical, intent(in) :: accept
    integer :: ierror

    !Need to fix the argument list to include accept, and then pass that to Python.  
    !This will allow users to have an update function that only runs on accepted moves if desired.  
    !If accept is false, we can skip the Python call entirely to save time for rejected moves.

    ! Call Python update
    ierror = call_py_noret(self%pyconstraint, "update", args=self%args_update)
    ! Don't error check - update is optional

  end subroutine
!=========================================================================
  subroutine PythonConstraint_Maintenance(self)
    implicit none
    class(PythonConstraint), intent(inout) :: self

    ! No-op for now, could add Python callback if needed

  end subroutine
!=========================================================================
  subroutine PythonConstraint_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(PythonConstraint), intent(inout) :: self
    logical :: accept
    integer :: ierror

    ! Check final constraint
    call self%CheckInitialConstraint(self%parent, accept)

    ! Call Python epilogue
    ierror = call_py_noret(self%pyconstraint, "epilogue", args=self%args_init)
    ! Don't error check - epilogue is optional

    write(nout,"(1x,A)") "(Python Constraint) Epilogue completed"

  end subroutine
!=========================================================================
  subroutine PythonConstraint_ProcessIO(self, line, lineStat)
    !-------
    ! IO Format: python module_name
    ! Example: python my_constraint
    !------
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(PythonConstraint), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=100) :: command

    lineStat = 0

    ! Get Python module name (2nd argument - 1st is "python")
    call GetXCommand(line, command, 2, lineStat)
    if(lineStat /= 0) then
      write(0,*) "ERROR! Python constraint requires module name"
      return
    endif
    self%pymodule = trim(adjustl(command))

  end subroutine
!=========================================================================
#else
! Stub type when EMBPYTHON is not defined
  type, public, extends(constraint) :: PythonConstraint
    integer :: boxID = 1
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => PythonConstraint_Stub_Constructor
      procedure, pass :: CheckInitialConstraint => PythonConstraint_Stub_CheckInitialConstraint
      procedure, pass :: DiffCheck => PythonConstraint_Stub_DiffCheck
      procedure, pass :: PostEnergy => PythonConstraint_Stub_PostEnergy
      procedure, pass :: ProcessIO => PythonConstraint_Stub_ProcessIO
      procedure, pass :: Maintenance => PythonConstraint_Stub_Maintenance
      procedure, pass :: Update => PythonConstraint_Stub_Update
      procedure, pass :: Epilogue => PythonConstraint_Stub_Epilogue
  end type

 contains
!=========================================================================
  subroutine PythonConstraint_Stub_Constructor(self, boxID)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    integer, intent(in) :: boxID

    write(0,*) "ERROR! Python constraints require compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    write(0,*) "ERROR! Python constraints require compilation with -DEMBPYTHON flag"
    error stop
    accept = .false.
  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept

    write(0,*) "ERROR! Python constraints require compilation with -DEMBPYTHON flag"
    error stop
    accept = .false.
  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_PostEnergy(self, trialBox, disp, E_Diff, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff
    logical, intent(out) :: accept

    write(0,*) "ERROR! Python constraints require compilation with -DEMBPYTHON flag"
    error stop
    accept = .false.
  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(PythonConstraint), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    write(0,*) "ERROR! Python constraints require compilation with -DEMBPYTHON flag"
    error stop
    lineStat = -1
  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_Maintenance(self)
    implicit none
    class(PythonConstraint), intent(inout) :: self

  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_Update(self, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    logical, intent(in) :: accept

  end subroutine
!=========================================================================
  subroutine PythonConstraint_Stub_Epilogue(self)
    implicit none
    class(PythonConstraint), intent(inout) :: self

  end subroutine
!=========================================================================
#endif
end module
!=========================================================================
