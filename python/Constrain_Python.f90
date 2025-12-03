!====================================================================
!  Constraint Module that is designed to call a Python Script 
!  in order to carry out analytical computations.
!  This module is designed to pass Classy style
!  object information into Python.  
!
!  def initialcheck(boxlist):
!
!  def diffcheck(boxlist, disp):
!====================================================================
module Constrain_Python
  use CoordinateTypes, only: Displacement, Deletion, Addition, OrthoVolChange
  use ParallelVar, only: nout

#ifdef EMBPYTHON
  use forpy_mod, only: dict,get_sys_path, list, call_py, module_py, import_py, &
                       object, call_py_noret, tuple, &
                       tuple_create, list_create, cast, err_print

  type, public, extends(constraint) :: PythonConstraint
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => PythonConstraint_Constructor
      procedure, pass :: CheckInitialConstraint => PythonConstraint_CheckInitialConstraint
      procedure, pass :: DiffCheck => PythonConstraint_DiffCheck
      procedure, pass :: ProcessIO => PythonConstraint_ProcessIO
      procedure, pass :: Maintenance => PythonConstraint_Maintenance
      procedure, pass :: Update => PythonConstraint_Update
      procedure, pass :: Epilogue => PythonConstraint_Epilogue

  end type
!=====================================================================
  contains
!=====================================================================
  subroutine PythonConstraint_Constructor(self, boxID)
    use BoxData, only: BoxArray
    implicit none
    class(PythonConstraint), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: AllocateStat
    integer :: nMolMax

    self%boxID = boxID
    self%parent => BoxArray(boxID) % box 


    IF (AllocateStat /= 0) STOP "Allocation Error in Distance Constraint"
  end subroutine
!=====================================================================
  subroutine PythonConstraint_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept


    accept = .true.
    write(nout,*) "Detailed Cluster Criteria Check Succeeded!"
    accept = .true.

  end subroutine
!=============================================================
  subroutine PythonConstraint_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept


    accept = .true.
  end subroutine
!=============================================================
  subroutine PythonConstraint_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ParallelVar, only: nout
    implicit none
    class(PythonConstraint), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%molType = intVal

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) realVal


  end subroutine
!====================================================================
  subroutine PythonConstraint_Maintenance(self)
    implicit none
    class(PythonConstraint), intent(inout) :: self

  end subroutine
!====================================================================
  subroutine PythonConstraint_Epilogue(self)
    implicit none
    class(PythonConstraint), intent(inout) :: self
    logical :: accept

    call self % CheckInitialConstraint(self%parent, accept)

  end subroutine
!=============================================================
  subroutine PythonConstraint_Update(self)
    implicit none
    class(PythonConstraint), intent(inout) :: self

  end subroutine
!=====================================================================
#endif
end module
!=====================================================================
