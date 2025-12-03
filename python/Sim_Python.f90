!===========================================================================
#define __StdErr__ 0
!===========================================================================
module SimPython
  use ParallelVar, only: myid, ierror, nout
  use VarPrecision
  use forpy_mod, only: dict,get_sys_path, list, call_py, module_py, import_py, &
                       object, call_py_noret, tuple, &
                       tuple_create

!===========================================================================
contains
!===========================================================================
  subroutine RunPythonTest(filename)
#ifdef PARALLEL
    use MPI
#endif

    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use ClassyPyObj, only: createboxdict
    use MCMoveData, only: Moves, MoveProb
    use MoveClassDef, only: MCMove
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use SimControl, only: nMoves, nCycles, screenfreq, configfreq, energyCheck
    use SimMonteCarlo, only: ScreenOut, Prologue, Epilogue, SafetyCheck, Trajectory
    use Units, only: outEngUnit

    implicit none
 
    logical :: accept
    integer :: i, j, nBoxes
    integer :: iAtom, moveNum, boxNum
    integer(kind=8) :: iCycle, iMove
    character(len=50) :: fileName
!    class(MCMove), pointer :: curMove
    type(list) :: paths
    type(module_py) :: testscript
    type(object) :: classyfunc
    type(tuple) :: args
    type(dict) :: boxdict
    real(dp), allocatable :: boxProb(:)

    iCycle = 0
    iMove = 0

    call Prologue
    call SafetyCheck
    ierror = get_sys_path(paths)
    ierror = paths%append(".")
    write(*,*) "Starting Python Unit Test"
    write(*,*) "Filename", filename
    ierror = import_py(testscript, trim(adjustl(filename)))

    boxdict = createboxdict(1)
    ierror = tuple_create(args, 1)
    ierror = args%setitem(0, boxdict)
    ierror = call_py_noret(testscript, "classyunit", args)


    call Epilogue

      
  end subroutine
!===========================================================================
end module
!===========================================================================
