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
!      
!      Returns:
!          dict: Dictionary containing:
!              - 'energy': float - energy difference (E_new - E_old)
!              - 'accept': bool - True if move is valid (no overlaps, etc.)
!              - 'denergytable': numpy array - per-atom energy differences (optional)
!      """
!      return {'energy': 0.0, 'accept': True}
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
                       cast, err_print, get_sys_path, is_none
!------------------------------------------------------------------------------
  type, public, extends(forcefield) :: Pair_PythonEnergy
    logical, private :: initialized = .false.
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
      procedure, pass :: ProcessIO => ProcessIO_PythonEnergy
      procedure, pass :: Prologue => Prologue_PythonEnergy
      procedure, pass :: GetCutOff => GetCutOff_PythonEnergy
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

    if(self%initialized) then
      return
    endif

    nBoxes = size(BoxArray)
    allocate(self%boxdicts(1:nBoxes))

    ! Add current directory to Python path
    ierror = get_sys_path(paths)
    errcheck_macro
    ierror = paths%append(".")
    errcheck_macro

    ! Load Python module
    write(nout,"(1x,A,A)") "Loading Python Energy Module: ", trim(adjustl(self%pymodule))
    ierror = import_py(self%pyenergy, trim(adjustl(self%pymodule)))
    errcheck_macro

    ! Create boxlist with no-copy references to Fortran arrays
    ierror = list_create(self%boxlist)
    errcheck_macro
    do iBox = 1, nBoxes
      self%boxdicts(iBox) = createboxdict_nocopy(iBox)
      ierror = self%boxlist%append(self%boxdicts(iBox))
      errcheck_macro
    enddo

    ! Create argument tuples
    ierror = tuple_create(self%args, 1)
    errcheck_macro
    ierror = self%args%setitem(0, self%boxlist)
    errcheck_macro

    ierror = tuple_create(self%diff_args, 2)
    errcheck_macro
    ierror = self%diff_args%setitem(0, self%boxlist)
    errcheck_macro

    ! Call Python prologue
    ierror = call_py_noret(self%pyenergy, "prologue", args=self%args)
    errcheck_macro

    self%initialized = .true.

    write(nout,"(1x,A,A)") "(Python Energy) Loaded module: ", trim(adjustl(self%pymodule))

  end subroutine
!=========================================================================
  subroutine DetailedECalc_PythonEnergy(self, curbox, E_T, accept)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(Pair_PythonEnergy), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    
    type(object) :: returnobj, item
    type(dict) :: resultdict
    type(ndarray) :: np_etable
    integer :: ierror, boxID
    real(dp) :: energy
    real(dp), pointer :: etable_ptr(:) => null()

    accept = .true.
    E_T = 0E0_dp
    boxID = curbox%boxID

    ! Call Python compute_total
    ierror = call_py(returnobj, self%pyenergy, "compute_total", args=self%args)
    errcheck_macro

    ! Cast return value to dict
    ierror = cast(resultdict, returnobj)
    errcheck_macro

    ! Get energy from result
    ierror = resultdict%getitem(item, "energy")
    if(ierror /= 0) then
      write(0,*) "ERROR! Python compute_total must return dict with 'energy' key"
      error stop
    endif
    ierror = cast(energy, item)
    errcheck_macro
    call item%destroy
    E_T = energy

    ! Check if accept flag is present
    ierror = resultdict%getitem(item, "accept")
    if(ierror == 0) then
      ierror = cast(accept, item)
      call item%destroy
    endif

    ! Check if energytable is provided (optional)
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
    use BoxData, only: BoxArray
    use ClassyPyObj, only: createdisplist
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
    integer :: ierror, iDisp
    real(dp) :: energy
    real(dp), pointer :: detable_ptr(:) => null()

    accept = .true.
    E_Diff = 0E0_dp
    curbox%dETable = 0E0_dp

    ! Create displacement list for Python
    displist = createdisplist(disp)

    ! Set displacement list in args
    ierror = self%diff_args%setitem(1, displist)
    errcheck_macro

    ! Call Python compute_diff
    ierror = call_py(returnobj, self%pyenergy, "compute_diff", args=self%diff_args)
    errcheck_macro

    ! Cast return value to dict
    ierror = cast(resultdict, returnobj)
    errcheck_macro

    ! Get energy difference from result
    ierror = resultdict%getitem(item, "energy")
    if(ierror /= 0) then
      write(0,*) "ERROR! Python compute_diff must return dict with 'energy' key"
      error stop
    endif
    ierror = cast(energy, item)
    errcheck_macro
    call item%destroy
    E_Diff = energy

    ! Check if accept flag is present
    ierror = resultdict%getitem(item, "accept")
    if(ierror == 0) then
      ierror = cast(accept, item)
      call item%destroy
    endif

    ! Check if denergytable is provided (optional)
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
      endif
      call item%destroy
    endif

    ! Clean up displacement list
    do iDisp = 0, size(disp)-1
      ierror = displist%getitem(item, iDisp)
      call item%destroy
    enddo
    call displist%destroy

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
