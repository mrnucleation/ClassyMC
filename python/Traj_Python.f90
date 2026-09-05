!=========================================================================
!  Trajectory Module that is designed to call a Python Script 
!  in order to write trajectory output.
!  This module is designed to pass Classy style
!  object information into Python for trajectory output.
!
!  This is particularly useful for interfacing with ASE (Atomic Simulation
!  Environment) or other Python-based analysis/visualization tools.
!
!  For a valid Trajectory script the following functions must be defined:
!
!  def prologue(boxinfo):
!      """Called once at the start of the simulation to initialize output.
!      
!      Args:
!          boxinfo: Dictionary containing:
!              - 'filename': str - output filename
!              - 'boxnum': int - box number
!              - 'outfreq': int - output frequency
!              - 'boxtype': str - type of box ('cube', 'ortho', etc.)
!              - 'temperature': float - temperature
!              - 'pressure': float - pressure
!              - 'volume': float - volume
!              - 'boxdimensions': numpy array (2 x ndim) - box bounds
!              - 'atomtype': numpy array - atom types (1-indexed)
!              - 'raw_atoms': numpy array (3 x nMaxAtoms) - atomic coordinates
!              - 'moleculecount': numpy array - molecules per type
!              - 'natoms': int - current number of atoms
!              - 'nmaxatoms': int - maximum atoms (for GC ensembles)
!              - 'moltype': numpy array - molecule type per atom
!              - 'molsubindx': numpy array - molecule sub-index per atom
!              - 'atomsymbols': list of str - atomic symbols by type
!      """
!      pass
!
!  def write_frame(boxinfo, cycle):
!      """Called each time a frame should be written.
!      
!      Args:
!          boxinfo: Dictionary with current box state (same as prologue)
!          cycle: int - current simulation cycle number
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
module Traj_Python
  use VarPrecision
  use TrajectoryTemplate, only: trajectory

#ifdef EMBPYTHON
  use forpy_mod, only: dict, dict_create, list, list_create, call_py, &
                       module_py, import_py, object, call_py_noret, tuple, &
                       tuple_create, ndarray, ndarray_create, ndarray_create_nocopy, &
                       cast, err_print, get_sys_path
!------------------------------------------------------------------------------
  type, public, extends(trajectory) :: PythonTraj
    logical, private :: initialized = .false.
    type(module_py) :: pytraj
    type(dict) :: boxinfo
    type(tuple) :: args_prologue
    type(tuple) :: args_frame
    type(tuple) :: args_epilogue
    character(len=100) :: pymodule
    contains
      procedure, pass :: WriteFrame => PythonTraj_WriteFrame
      procedure, pass :: Prologue => PythonTraj_Prologue
      procedure, pass :: Epilogue => PythonTraj_Epilogue
      procedure, pass :: ProcessIO => PythonTraj_ProcessIO
  end type

 contains
!=========================================================================
  subroutine PythonTraj_Prologue(self)
    use BoxData, only: BoxArray
    use CubicBoxDef, only: CubeBox
    use OrthoBoxDef, only: OrthoBox
    use SimpleSimBox, only: SimpleBox
    use Common_MolInfo, only: AtomData, nAtomTypes
    use Input_Format, only: ReplaceText
    use ParallelVar, only: nout, myid
    implicit none
    class(PythonTraj), intent(inout) :: self
    integer :: ierror
    integer :: boxNum, iType, iCharacter
    type(list) :: paths, atomsymbols
    type(ndarray) :: np_atoms, np_boxdim, np_atomtype, np_nmol
    type(ndarray) :: np_moltype, np_molsubindx
    real(dp), allocatable :: boxDim(:,:)
    character(len=30) :: idString

    if(self%initialized) then
      return
    endif

    boxNum = self%boxNum

    ! Handle filename substitution (& -> process id)
    do iCharacter = 1, len(self%filename)
      if(self%filename(iCharacter:iCharacter) == "&") then
        write(idString, *) myid
        self%filename = ReplaceText(self%filename, "&", trim(adjustl(idString)))
        exit
      endif
    enddo

    ! Load Python module
    write(nout,"(1x,A,A)") "Loading Python Trajectory Module: ", trim(adjustl(self%pymodule))
    ierror = get_sys_path(paths)
    errcheck_macro
    ierror = paths%append(".")
    errcheck_macro
    call paths%destroy

    ierror = import_py(self%pytraj, trim(adjustl(self%pymodule)))
    errcheck_macro

    ! Create the boxinfo dictionary
    ierror = dict_create(self%boxinfo)
    errcheck_macro

    ! Add filename and settings
    ierror = self%boxinfo%setitem("filename", trim(adjustl(self%fileName)))
    errcheck_macro
    ierror = self%boxinfo%setitem("boxnum", boxNum)
    errcheck_macro
    ierror = self%boxinfo%setitem("outfreq", self%outFreq)
    errcheck_macro

    ! Add box type
    select type(box => BoxArray(boxNum)%box)
      type is(SimpleBox)
        ierror = self%boxinfo%setitem("boxtype", "nobox")
      class is(CubeBox)
        ierror = self%boxinfo%setitem("boxtype", "cube")
      class is(OrthoBox)
        ierror = self%boxinfo%setitem("boxtype", "ortho")
      class default
        ierror = self%boxinfo%setitem("boxtype", "unknown")
    end select
    errcheck_macro

    ! Add thermodynamic properties
    ierror = self%boxinfo%setitem("temperature", BoxArray(boxNum)%box%temperature)
    errcheck_macro
    ierror = self%boxinfo%setitem("pressure", BoxArray(boxNum)%box%pressure)
    errcheck_macro
    ierror = self%boxinfo%setitem("volume", BoxArray(boxNum)%box%volume)
    errcheck_macro
    ierror = self%boxinfo%setitem("natoms", BoxArray(boxNum)%box%nAtoms)
    errcheck_macro
    ierror = self%boxinfo%setitem("nmaxatoms", BoxArray(boxNum)%box%nMaxAtoms)
    errcheck_macro
    ierror = self%boxinfo%setitem("ndimension", BoxArray(boxNum)%box%nDimension)
    errcheck_macro

    ! Add box dimensions
    allocate(boxDim(1:2, 1:BoxArray(boxNum)%box%nDimension))
    call BoxArray(boxNum)%box%GetDimensions(boxDim)
    ierror = ndarray_create(np_boxdim, boxDim)
    errcheck_macro
    ierror = self%boxinfo%setitem("boxdimensions", np_boxdim)
    errcheck_macro
    call np_boxdim%destroy

    ! Add atomic data with no-copy for efficiency
    ierror = ndarray_create_nocopy(np_atoms, BoxArray(boxNum)%box%atoms)
    errcheck_macro
    ierror = self%boxinfo%setitem("raw_atoms", np_atoms)
    errcheck_macro
    call np_atoms%destroy

    ierror = ndarray_create_nocopy(np_atomtype, BoxArray(boxNum)%box%AtomType)
    errcheck_macro
    ierror = self%boxinfo%setitem("atomtype", np_atomtype)
    errcheck_macro
    call np_atomtype%destroy

    ierror = ndarray_create_nocopy(np_nmol, BoxArray(boxNum)%box%NMol)
    errcheck_macro
    ierror = self%boxinfo%setitem("moleculecount", np_nmol)
    errcheck_macro
    call np_nmol%destroy

    ierror = ndarray_create_nocopy(np_moltype, BoxArray(boxNum)%box%MolType)
    errcheck_macro
    ierror = self%boxinfo%setitem("moltype", np_moltype)
    errcheck_macro
    call np_moltype%destroy

    ierror = ndarray_create_nocopy(np_molsubindx, BoxArray(boxNum)%box%MolSubIndx)
    errcheck_macro
    ierror = self%boxinfo%setitem("molsubindx", np_molsubindx)
    errcheck_macro
    call np_molsubindx%destroy

    ! Add atom symbols list
    ierror = list_create(atomsymbols)
    errcheck_macro
    do iType = 1, nAtomTypes
      ierror = atomsymbols%append(trim(adjustl(AtomData(iType)%Symb)))
      errcheck_macro
    enddo
    ierror = self%boxinfo%setitem("atomsymbols", atomsymbols)
    errcheck_macro
    call atomsymbols%destroy

    ! Create argument tuples
    ierror = tuple_create(self%args_prologue, 1)
    errcheck_macro
    ierror = self%args_prologue%setitem(0, self%boxinfo)
    errcheck_macro

    ierror = tuple_create(self%args_frame, 2)
    errcheck_macro
    ierror = self%args_frame%setitem(0, self%boxinfo)
    errcheck_macro

    ierror = tuple_create(self%args_epilogue, 1)
    errcheck_macro
    ierror = self%args_epilogue%setitem(0, self%boxinfo)
    errcheck_macro

    ! Call Python prologue
    ierror = call_py_noret(self%pytraj, "prologue", args=self%args_prologue)
    errcheck_macro

    self%initialized = .true.

    write(nout,"(1x,A,A)") "(Python Trajectory) Loaded module: ", trim(adjustl(self%pymodule))
    write(nout,"(1x,A,A)") "(Python Trajectory) Output file: ", trim(adjustl(self%fileName))

    ! Write initial frame
    call self%WriteFrame(0_8)

  end subroutine
!=========================================================================
  subroutine PythonTraj_WriteFrame(self, iCycle)
    use BoxData, only: BoxArray
    use ClassyPyObj, only: refreshboxdict
    implicit none
    class(PythonTraj), intent(inout) :: self
    integer(kind=8), intent(in) :: iCycle
    integer :: ierror, boxNum

    boxNum = self%boxNum

    call refreshboxdict(self%boxinfo, boxNum)
    ierror = self%boxinfo%setitem("natoms", BoxArray(boxNum)%box%nAtoms)
    ierror = self%boxinfo%setitem("cycle", iCycle)

    ierror = self%args_frame%setitem(1, iCycle)
    errcheck_macro

    ierror = call_py_noret(self%pytraj, "write_frame", args=self%args_frame)
    errcheck_macro

  end subroutine
!=========================================================================
  subroutine PythonTraj_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(PythonTraj), intent(inout) :: self
    integer :: ierror

    ! Call Python epilogue for cleanup
    ierror = call_py_noret(self%pytraj, "epilogue", args=self%args_epilogue)
    errcheck_macro

    write(nout,"(1x,A)") "(Python Trajectory) Epilogue completed"

  end subroutine
!=========================================================================
  subroutine PythonTraj_ProcessIO(self, line, lineStat)
    !-------
    ! IO Format: python boxnum freq filename module_name
    ! Example: python 1 1000 traj.xyz ase_trajectory
    !------
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(PythonTraj), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=100) :: command
    integer :: i, intVal

    lineStat = 0

    ! Get box number (2nd argument - 1st is "python")
    call GetXCommand(line, command, 2, lineStat)
    if(lineStat /= 0) then
      write(0,*) "ERROR! Python trajectory requires box number"
      return
    endif
    read(command, *) intVal
    call self%SetBox(intVal)

    ! Get output frequency (3rd argument)
    call GetXCommand(line, command, 3, lineStat)
    if(lineStat /= 0) then
      write(0,*) "ERROR! Python trajectory requires output frequency"
      return
    endif
    read(command, *) intVal
    call self%SetFreq(intVal)

    ! Get output filename (4th argument)
    call GetXCommand(line, command, 4, lineStat)
    if(lineStat /= 0) then
      write(0,*) "ERROR! Python trajectory requires output filename"
      return
    endif
    do i = 1, len(command)
      if(command(i:i) == '"') then
        command(i:i) = " "
      endif
    enddo
    call self%SetFileName(command)

    ! Get Python module name (5th argument)
    call GetXCommand(line, command, 5, lineStat)
    if(lineStat /= 0) then
      write(0,*) "ERROR! Python trajectory requires module name"
      return
    endif
    self%pymodule = trim(adjustl(command))

  end subroutine
!=========================================================================
#else
! Stub type when EMBPYTHON is not defined
  type, public, extends(trajectory) :: PythonTraj
    contains
      procedure, pass :: WriteFrame => PythonTraj_Stub_WriteFrame
      procedure, pass :: Prologue => PythonTraj_Stub_Prologue
      procedure, pass :: Epilogue => PythonTraj_Stub_Epilogue
      procedure, pass :: ProcessIO => PythonTraj_Stub_ProcessIO
  end type

 contains
!=========================================================================
  subroutine PythonTraj_Stub_WriteFrame(self, iCycle)
    implicit none
    class(PythonTraj), intent(inout) :: self
    integer(kind=8), intent(in) :: iCycle

    write(0,*) "ERROR! Python trajectories require compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  subroutine PythonTraj_Stub_Prologue(self)
    implicit none
    class(PythonTraj), intent(inout) :: self

    write(0,*) "ERROR! Python trajectories require compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  subroutine PythonTraj_Stub_Epilogue(self)
    implicit none
    class(PythonTraj), intent(inout) :: self

    write(0,*) "ERROR! Python trajectories require compilation with -DEMBPYTHON flag"
    error stop
  end subroutine
!=========================================================================
  subroutine PythonTraj_Stub_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(PythonTraj), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    write(0,*) "ERROR! Python trajectories require compilation with -DEMBPYTHON flag"
    error stop
    lineStat = -1
  end subroutine
!=========================================================================
#endif
end module
!=========================================================================
