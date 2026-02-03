!=========================================================================
!  Move Module that is designed to call a Python Script 
!  in order to carry out Monte Carlo moves.
!  This module is designed to pass Classy style
!  object information into Python.  
!
!  For a valid Move script the following functions must be defined:
!
!  def get_move_type():
!      """Return the displacement type this move will generate.
!      
!      Returns:
!          str: One of:
!              - 'displacement' - single atom/molecule translation
!              - 'deletion' - remove a molecule
!              - 'addition' - add a molecule
!              - 'atomexchange' - exchange atom types
!              - 'volchange' - isotropic volume change
!              - 'orthovolchange' - anisotropic volume change
!      """
!      return 'displacement'
!
!  def prologue(boxlist):
!      """Called once at the start of the simulation to initialize the move."""
!      pass
!
!  def generate_move(boxlist):
!      """Generate a trial move and return displacement information.
!      
!      Returns depend on move type (from get_move_type):
!      
!      For 'displacement':
!          {'moltype': int, 'molindex': int, 'atomindex': int,
!           'x_new': float, 'y_new': float, 'z_new': float}
!      
!      For 'deletion':
!          {'moltype': int, 'molindex': int, 'atomindex': int}
!      
!      For 'addition':
!          {'moltype': int, 'molindex': int, 'atomindex': int,
!           'x_new': float, 'y_new': float, 'z_new': float}
!      
!      For 'atomexchange':
!          {'oldatomindex': int, 'oldtype': int, 
!           'newatomindex': int, 'newtype': int}
!      
!      For 'volchange':
!          {'volnew': float, 'volold': float}
!      
!      For 'orthovolchange':
!          {'volnew': float, 'volold': float,
!           'xscale': float, 'yscale': float, 'zscale': float}
!      
!      Any type can also return:
!          {'reject': True} - to immediately reject the move
!      """
!      return {'reject': True}
!
!  def accept_move(boxlist, accepted, forward_prob, reverse_prob):
!      """Called after the move decision to allow state updates.
!      
!      Args:
!          boxlist: List of box dictionaries
!          accepted: bool - True if move was accepted
!          forward_prob: float - forward generation probability used
!          reverse_prob: float - reverse generation probability used
!      """
!      pass
!
!  def compute_reverse_prob(boxlist, move_dict):
!      """Compute the reverse generation probability for detailed balance.
!      
!      This function is called after generate_move() returns, before the
!      acceptance decision. It receives the same boxlist (still in the OLD
!      configuration) and the move_dict returned by generate_move().
!      
!      The reverse probability is the probability of generating the reverse
!      move (going from the trial state back to the current state).
!      
!      Args:
!          boxlist: List of box dictionaries (current/old configuration)
!          move_dict: Dictionary returned by generate_move()
!      
!      Returns:
!          float: Reverse generation probability. Return 1.0 if not needed
!                 (symmetric moves) or if forward_prob already accounts for it.
!      """
!      return 1.0
!
!=========================================================================
#define errcheck_macro if(ierror/=0) then;call err_print;stop;endif
!=========================================================================
module MCMove_Python
  use MoveClassDef
  use CoordinateTypes, only: Perturbation, Displacement, Deletion, Addition, &
                             AtomExchange, VolChange, OrthoVolChange
  use SimpleSimBox, only: SimpleBox
  use VarPrecision
  use iso_c_binding, only: C_CHAR
#ifdef EMBPYTHON
  use forpy_mod, only: dict, dict_create, get_sys_path, list, call_py, module_py, import_py, &
                       object, call_py_noret, tuple, &
                       tuple_create, list_create, cast, err_print, is_none
#endif
  
  ! Move type constants
  integer, parameter :: MOVE_DISPLACEMENT = 1
  integer, parameter :: MOVE_DELETION = 2
  integer, parameter :: MOVE_ADDITION = 3
  integer, parameter :: MOVE_ATOMEXCHANGE = 4
  integer, parameter :: MOVE_VOLCHANGE = 5
  integer, parameter :: MOVE_ORTHOVOLCHANGE = 6

#ifdef EMBPYTHON
!------------------------------------------------------------------------------
  type, public, extends(MCMove) :: PythonMove
    logical, private :: initialized = .false.
    integer :: moveType = MOVE_DISPLACEMENT
    type(module_py) :: pymove
    type(list) :: boxlist
    type(tuple) :: args
    type(tuple) :: accept_args
    type(tuple) :: reverse_prob_args
    type(dict), allocatable :: boxdicts(:)
    character(len=100) :: filename
    
    ! Polymorphic perturbation array - can hold any perturbation type
    class(Perturbation), allocatable :: perturb(:)
    
    !Rejection Counters
    integer :: ovlaprej = 0 
    integer :: constrainrej = 0 
    integer :: detailedrej = 0 

    contains
      procedure, pass :: Constructor => PythonMove_Constructor
      procedure, pass :: FullMove => PythonMove_FullMove
      procedure, pass :: Prologue => PythonMove_Prologue
      procedure, pass :: Epilogue => PythonMove_Epilogue
      procedure, pass :: ProcessIO => PythonMove_ProcessIO
  end type

 contains
!=========================================================================
  subroutine PythonMove_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    use BoxData, only: BoxArray
    implicit none
    class(PythonMove), intent(inout) :: self
    integer :: iType, maxAtoms, nBoxes

    nBoxes = size(boxArray)
    if(.not. allocated(self%boxProb)) then
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo

    call self%CreateTempArray(maxAtoms)

  end subroutine
!=========================================================================
  subroutine PythonMove_Prologue(self)
    use BoxData, only: BoxArray
    use ClassyPyObj, only: createboxdict_nocopy
    use Common_MolInfo, only: MolData, nMolTypes
    use Input_Format, only: ReplaceText
    use ParallelVar, only: nout
    implicit none
    class(PythonMove), intent(inout) :: self
    integer :: ierror
    integer :: iBox, nBoxes, iType, maxAtoms
    type(object) :: returnobj
    character(kind=C_CHAR, len=:), allocatable :: moveTypeStr

    if(self%initialized) then
      return
    endif

    nBoxes = size(BoxArray)
    allocate( self%boxdicts(1:nBoxes) )

    !Open up the Python Script 
    write(nout,*) "Loading Python Move Module: ", trim(adjustl(self%filename))
    ierror = import_py(self%pymove, trim(adjustl(self%filename)))
    errcheck_macro

    ierror = list_create(self%boxlist)
    errcheck_macro
    do iBox = 1, size(BoxArray)
      self%boxdicts(iBox) = createboxdict_nocopy(iBox)
      ierror = self%boxlist%append( self%boxdicts(iBox) )
      errcheck_macro
    enddo

    ierror = tuple_create(self%args, 1)
    ierror = self%args%setitem(0, self%boxlist)
    errcheck_macro

    ! accept_args: (boxlist, accepted, forward_prob, reverse_prob)
    ierror = tuple_create(self%accept_args, 4)
    errcheck_macro

    ! reverse_prob_args: (boxlist, move_dict)
    ierror = tuple_create(self%reverse_prob_args, 2)
    errcheck_macro

    ! Query Python for the move type
    ierror = call_py(returnobj, self%pymove, "get_move_type")
    errcheck_macro
    ierror = cast(moveTypeStr, returnobj)
    errcheck_macro
    call returnobj%destroy

    ! Parse the move type string and allocate appropriate perturbation type
    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo

    select case(trim(adjustl(moveTypeStr)))
      case("displacement")
        self%moveType = MOVE_DISPLACEMENT
        allocate(Displacement :: self%perturb(1:maxAtoms))
        write(nout,"(1x,A)") "(Python Move) Move type: displacement"
        
      case("deletion")
        self%moveType = MOVE_DELETION
        allocate(Deletion :: self%perturb(1:maxAtoms))
        write(nout,"(1x,A)") "(Python Move) Move type: deletion"
        
      case("addition")
        self%moveType = MOVE_ADDITION
        allocate(Addition :: self%perturb(1:maxAtoms))
        write(nout,"(1x,A)") "(Python Move) Move type: addition"
        
      case("atomexchange")
        self%moveType = MOVE_ATOMEXCHANGE
        allocate(AtomExchange :: self%perturb(1:1))
        write(nout,"(1x,A)") "(Python Move) Move type: atomexchange"
        
      case("volchange")
        self%moveType = MOVE_VOLCHANGE
        allocate(VolChange :: self%perturb(1:1))
        write(nout,"(1x,A)") "(Python Move) Move type: volchange"
        
      case("orthovolchange")
        self%moveType = MOVE_ORTHOVOLCHANGE
        allocate(OrthoVolChange :: self%perturb(1:1))
        write(nout,"(1x,A)") "(Python Move) Move type: orthovolchange"
        
      case default
        write(0,*) "ERROR! Unknown move type from Python: ", trim(adjustl(moveTypeStr))
        write(0,*) "Valid types: displacement, deletion, addition, atomexchange, volchange, orthovolchange"
        error stop
    end select

    ! Call Python prologue function
    ierror = call_py_noret(self%pymove, "prologue", args=self%args)
    errcheck_macro

    self%initialized = .true.
    
    ! Also call Constructor for templist allocation
    call self%Constructor

    write(nout,"(1x,A,A)") "(Python Move) Loaded module: ", trim(adjustl(self%filename))

  end subroutine
!=========================================================================
  subroutine PythonMove_FullMove(self, trialBox, accept) 
    use CommonSampling, only: sampling
    use Common_MolInfo, only: MolData, nMolTypes
    use ParallelVar, only: nout
    implicit none
    class(PythonMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    
    type(object) :: returnobj, item
    type(dict) :: dispdict
    integer :: ierror
    integer :: boxID
    integer :: molType, molIndx, atmIndx
    integer :: oldAtmIndx, oldType, newAtmIndx, newType
    real(dp) :: x_new, y_new, z_new
    real(dp) :: volNew, volOld, xScale, yScale, zScale
    real(dp) :: E_Diff, E_Inter, E_Intra
    real(dp) :: Prob, forward_prob, reverse_prob
    logical :: reject_flag

    boxID = trialBox % boxID
    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    ! Call Python generate_move function
    ierror = call_py(returnobj, self%pymove, "generate_move", args=self%args)
    errcheck_macro

    ! Cast the return value to a dict
    ierror = cast(dispdict, returnobj)
    errcheck_macro

    ! Check if the move should be immediately rejected
    ierror = dispdict%getitem(item, "reject")
    if(ierror == 0) then
      ierror = cast(reject_flag, item)
      call item%destroy
      if(reject_flag) then
        accept = .false.
        call dispdict%destroy
        call returnobj%destroy
        return
      endif
    endif

    ! Parse forward_prob from the return dict (default 1.0)
    forward_prob = 1E0_dp
    ierror = dispdict%getitem(item, "forward_prob")
    if(ierror == 0) then
      ierror = cast(forward_prob, item)
      call item%destroy
    endif

    ! Parse displacement based on move type
    select case(self%moveType)
      case(MOVE_DISPLACEMENT)
        call ParseDisplacement(dispdict, accept)
        
      case(MOVE_DELETION)
        call ParseDeletion(dispdict, accept)
        
      case(MOVE_ADDITION)
        call ParseAddition(dispdict, accept)
        
      case(MOVE_ATOMEXCHANGE)
        call ParseAtomExchange(dispdict, accept)
        
      case(MOVE_VOLCHANGE)
        call ParseVolChange(dispdict, accept)
        
      case(MOVE_ORTHOVOLCHANGE)
        call ParseOrthoVolChange(dispdict, accept)
    end select

    if(.not. accept) then
      call dispdict%destroy
      call returnobj%destroy
      return
    endif

    ! Call Python compute_reverse_prob function to get reverse probability
    reverse_prob = 1E0_dp
    ierror = self%reverse_prob_args%setitem(0, self%boxlist)
    ierror = self%reverse_prob_args%setitem(1, dispdict)
    call returnobj%destroy  ! Reuse returnobj
    ierror = call_py(returnobj, self%pymove, "compute_reverse_prob", args=self%reverse_prob_args)
    if(ierror == 0) then
      ierror = cast(reverse_prob, returnobj)
      if(ierror /= 0) reverse_prob = 1E0_dp
    endif
    call returnobj%destroy
    call dispdict%destroy

    ! Compute generation probability ratio for detailed balance
    ! acc = min(1, (reverse_prob / forward_prob) * exp(-beta * dE))
    if(forward_prob > 0E0_dp) then
      Prob = reverse_prob / forward_prob
    else
      ! Invalid forward probability - reject move
      accept = .false.
      return
    endif

    ! Set box ID on perturbation
    self%perturb(1)%boxID = boxID

    !Check Constraint
    accept = trialBox % CheckConstraint( self%perturb(1:1) )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      call NotifyPython(.false.)
      return
    endif

    !Energy Calculation
    call trialBox%ComputeEnergyDelta(self%perturb(1:1),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      call NotifyPython(.false.)
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%perturb(1:1), E_Diff )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      call NotifyPython(.false.)
      return
    endif

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%perturb(1:1), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%perturb(1:1), self%tempList, self%tempNNei)
    else
      self%detailedrej = self%detailedrej + 1
    endif

    ! Notify Python of the acceptance/rejection
    call NotifyPython(accept)

  contains
    !-----------------------------------------------------------------------
    subroutine NotifyPython(accepted)
      logical, intent(in) :: accepted
      integer :: ierr

      ierr = self%accept_args%setitem(0, self%boxlist)
      ierr = self%accept_args%setitem(1, accepted)
      ierr = self%accept_args%setitem(2, forward_prob)
      ierr = self%accept_args%setitem(3, reverse_prob)
      ierr = call_py_noret(self%pymove, "accept_move", args=self%accept_args)
    end subroutine
    !-----------------------------------------------------------------------
    subroutine ParseDisplacement(ddict, ok)
      type(dict), intent(inout) :: ddict
      logical, intent(out) :: ok
      integer :: ierr
      
      ok = .true.
      
      ierr = ddict%getitem(item, "moltype")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(molType, item)
      call item%destroy

      ierr = ddict%getitem(item, "molindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(molIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "atomindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(atmIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "x_new")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(x_new, item)
      call item%destroy

      ierr = ddict%getitem(item, "y_new")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(y_new, item)
      call item%destroy

      ierr = ddict%getitem(item, "z_new")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(z_new, item)
      call item%destroy

      select type(p => self%perturb(1))
        type is(Displacement)
          p%molType = molType
          p%molIndx = molIndx
          p%atmIndx = atmIndx
          p%x_new = x_new
          p%y_new = y_new
          p%z_new = z_new
          p%newlist = .false.
          p%listIndex = 1
      end select
    end subroutine
    !-----------------------------------------------------------------------
    subroutine ParseDeletion(ddict, ok)
      type(dict), intent(inout) :: ddict
      logical, intent(out) :: ok
      integer :: ierr
      
      ok = .true.
      
      ierr = ddict%getitem(item, "moltype")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(molType, item)
      call item%destroy

      ierr = ddict%getitem(item, "molindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(molIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "atomindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(atmIndx, item)
      call item%destroy

      select type(p => self%perturb(1))
        type is(Deletion)
          p%molType = molType
          p%molIndx = molIndx
          p%atmIndx = atmIndx
      end select
    end subroutine
    !-----------------------------------------------------------------------
    subroutine ParseAddition(ddict, ok)
      type(dict), intent(inout) :: ddict
      logical, intent(out) :: ok
      integer :: ierr
      
      ok = .true.
      
      ierr = ddict%getitem(item, "moltype")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(molType, item)
      call item%destroy

      ierr = ddict%getitem(item, "molindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(molIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "atomindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(atmIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "x_new")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(x_new, item)
      call item%destroy

      ierr = ddict%getitem(item, "y_new")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(y_new, item)
      call item%destroy

      ierr = ddict%getitem(item, "z_new")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(z_new, item)
      call item%destroy

      select type(p => self%perturb(1))
        type is(Addition)
          p%molType = molType
          p%molIndx = molIndx
          p%atmIndx = atmIndx
          p%x_new = x_new
          p%y_new = y_new
          p%z_new = z_new
          p%listIndex = 1
      end select
    end subroutine
    !-----------------------------------------------------------------------
    subroutine ParseAtomExchange(ddict, ok)
      type(dict), intent(inout) :: ddict
      logical, intent(out) :: ok
      integer :: ierr
      
      ok = .true.
      
      ierr = ddict%getitem(item, "oldatomindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(oldAtmIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "oldtype")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(oldType, item)
      call item%destroy

      ierr = ddict%getitem(item, "newatomindex")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(newAtmIndx, item)
      call item%destroy

      ierr = ddict%getitem(item, "newtype")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(newType, item)
      call item%destroy

      select type(p => self%perturb(1))
        type is(AtomExchange)
          p%oldAtmIndx = oldAtmIndx
          p%oldType = oldType
          p%newAtmIndx = newAtmIndx
          p%newType = newType
      end select
    end subroutine
    !-----------------------------------------------------------------------
    subroutine ParseVolChange(ddict, ok)
      type(dict), intent(inout) :: ddict
      logical, intent(out) :: ok
      integer :: ierr
      
      ok = .true.
      
      ierr = ddict%getitem(item, "volnew")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(volNew, item)
      call item%destroy

      ierr = ddict%getitem(item, "volold")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(volOld, item)
      call item%destroy

      select type(p => self%perturb(1))
        type is(VolChange)
          p%volNew = volNew
          p%volOld = volOld
      end select
    end subroutine
    !-----------------------------------------------------------------------
    subroutine ParseOrthoVolChange(ddict, ok)
      type(dict), intent(inout) :: ddict
      logical, intent(out) :: ok
      integer :: ierr
      
      ok = .true.
      
      ierr = ddict%getitem(item, "volnew")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(volNew, item)
      call item%destroy

      ierr = ddict%getitem(item, "volold")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(volOld, item)
      call item%destroy

      ierr = ddict%getitem(item, "xscale")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(xScale, item)
      call item%destroy

      ierr = ddict%getitem(item, "yscale")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(yScale, item)
      call item%destroy

      ierr = ddict%getitem(item, "zscale")
      if(ierr /= 0) then; ok = .false.; return; endif
      ierr = cast(zScale, item)
      call item%destroy

      select type(p => self%perturb(1))
        type is(OrthoVolChange)
          p%volNew = volNew
          p%volOld = volOld
          p%xScale = xScale
          p%yScale = yScale
          p%zScale = zScale
      end select
    end subroutine
    !-----------------------------------------------------------------------

  end subroutine
!=========================================================================
  subroutine PythonMove_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(PythonMove), intent(inout) :: self
    real(dp) :: accptRate
 
    write(nout,*) 
    write(nout,"(1x,A,I15)") "Python Move Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Python Move Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Python Move Acceptance Rate: ", accptRate
    write(nout, "(1x,A,I15)") "Python Move, Rejections due to overlap:", self%ovlaprej
    write(nout, "(1x,A,I15)") "Python Move, Rejections due to constraint:", self%constrainrej
    write(nout, "(1x,A,I15)") "Python Move, Rejections due to detailed balance:", self%detailedrej

  end subroutine
!=========================================================================
  subroutine PythonMove_ProcessIO(self, line, lineStat)
    !-------
    ! IO Format: python (probability) (python module name) [options]
    !------
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(PythonMove), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=100) :: command

    lineStat = 0

    ! Get the Python module name (3rd argument after "python" and probability)
    call GetXCommand(line, command, 3, lineStat)
    if(lineStat /= 0) then
      write(0,*) "ERROR! Python move requires a module name"
      lineStat = -1
      return
    endif
    self%filename = trim(adjustl(command))

  end subroutine
!=========================================================================
#else
! Stub type when EMBPYTHON is not defined
  type, public, extends(MCMove) :: PythonMove
    contains
      procedure, pass :: FullMove => PythonMove_Stub_FullMove
      procedure, pass :: ProcessIO => PythonMove_Stub_ProcessIO
  end type

 contains
!=========================================================================
  subroutine PythonMove_Stub_FullMove(self, trialBox, accept)
    implicit none
    class(PythonMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    write(0,*) "ERROR! Python moves require compilation with -DEMBPYTHON flag"
    error stop
    accept = .false.
  end subroutine
!=========================================================================
  subroutine PythonMove_Stub_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(PythonMove), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat

    write(0,*) "ERROR! Python moves require compilation with -DEMBPYTHON flag"
    error stop
    lineStat = -1
  end subroutine
!=========================================================================
#endif
end module
!=========================================================================
