!=========================================================================
!  Analysis Module that is designed to call a Python Script 
!  in order to carry out analytical computations.
!  This module is designed to pass Classy style
!  object information into Python.  
!
!  def myfunction(boxlist):
!=========================================================================
#define errcheck_macro if(ierror/=0) then;call err_print;stop;endif
!=========================================================================
module Anaylsis_PythonFunc
  use AnaylsisClassDef, only: Analysis
  use AnalysisData, only: analyCommon
  use VarPrecision
  use SimpleSimBox, only: SimpleBox
#ifdef EMBPYTHON
  use forpy_mod, only: dict,get_sys_path, list, call_py, module_py, import_py, &
                       object, call_py_noret, tuple, &
                       tuple_create, list_create, cast, err_print
!------------------------------------------------------------------------------
  type, public, extends(Analysis):: PythonFunc
    logical, private :: new
    real(dp) :: lastvalue
    real(dp) :: newvalue
    type(module_py) :: pyanalysis
    type(list) :: boxlist
    type(tuple) :: args, newargs
    type(dict), allocatable :: boxdicts(:)
    character(len=30) :: filename
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Prologue => PythonFunc_Prologue
      procedure, pass :: Compute => PythonFunc_Compute
      procedure, pass :: CalcNewState => PythonFunc_CalcNewState
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => PythonFunc_ProcessIO
!      procedure, pass :: WriteInfo => PythonFunc_WriteInfo
      procedure, pass :: GetResult => PythonFunc_GetResult
      procedure, pass :: CastCommonType => PythonFunc_CastCommonType
!      procedure, pass :: Finalize => PythonFunc_Finalize
  end type

 contains
!=========================================================================
  subroutine PythonFunc_Prologue(self)
    use BoxData, only: BoxArray
    use ClassyPyObj, only: createboxdict,createboxdict_nocopy
    use Input_Format, only: ReplaceText
    use ParallelVar, only: nout
    implicit none
    class(PythonFunc), intent(inout) :: self
    logical :: accept
    integer :: ierror
    integer :: iBox, nBoxes

    if(.not. allocated(self%boxdicts)) then
      nBoxes = size(BoxArray)
      allocate( self%boxdicts(1:nBoxes) )
    else
      return
    endif
    !Open up the Python Script and 
    write(nout,*) "Loading Python Analysis Module: ", self%filename
    ierror = import_py(self%pyanalysis, trim(adjustl(self%filename)))
    errcheck_macro



    ierror = list_create(self%boxlist)
!    errcheck_macro
    do iBox = 1, size(BoxArray)
!      self%boxdicts(iBox) = createboxdict(iBox)
      self%boxdicts(iBox) = createboxdict_nocopy(iBox)
      ierror = self%boxlist%append( self%boxdicts(iBox) )
      errcheck_macro
    enddo

    ierror = tuple_create(self%args, 1)
    ierror = self%args%setitem(0, self%boxlist)

    ierror = tuple_create(self%newargs, 2)
    errcheck_macro
    ierror = self%newargs%setitem(0, self%boxlist)

    self%perMove = .true.
    accept = .true.
    call self%Compute(accept)

  end subroutine
!=========================================================================
!
!    Epilogue code
!    do iBox = 1, size(BoxArray)
!      call self%boxdicts(iBox)%destroy
!      errcheck_macro
!    enddo
!    call self%boxlist%destroy
!=========================================================================
  subroutine PythonFunc_Compute(self, accept)
    use BoxData, only: BoxArray
    use ClassyPyObj, only: createboxdict
    implicit none
    class(PythonFunc), intent(inout) :: self
    logical, intent(in) :: accept
    type(object) :: returnval
    integer :: iBox
    integer :: ierror



    if(.not. accept) then
      return
    endif

    if(self%new) then
      self%lastvalue = self%newvalue
      self%new = .false.
!      ierror = call_py_noret(self%pyanalysis, "update", args=self%args)
      return
    endif

    ierror = call_py(returnval, self%pyanalysis, "compute", args=self%args)
    errcheck_macro
    ierror = cast(self%lastvalue, returnval)    
    errcheck_macro

    if(isnan(self%lastvalue)) then
      error stop "NaN value returned from Python compute function!"
    endif



  end subroutine
!=========================================================================
  subroutine PythonFunc_CalcNewState(self, disp, newVal)
    use CoordinateTypes, only: Displacement, Perturbation
    use ClassyPyObj, only: createdisplist
    implicit none
    class(PythonFunc), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    integer :: iDisp
    real(dp), intent(in), optional :: newVal
    type(object) :: returnobj
    type(list) :: displist
    integer :: ierror
    real(dp) :: returnval

    displist = createdisplist(disp)



    ierror = self%newargs%setitem(1, displist)
    errcheck_macro

    ierror = call_py(returnobj, self%pyanalysis, "compute_new", args=self%newargs)
    errcheck_macro
    ierror = cast(returnval, returnobj)    
    errcheck_macro

    if(isnan(returnval)) then
      error stop "NaN value returned from Python compute_new function!"
    endif

!    call self%newargs%destroy
    call displist%destroy

    self%new = .true.
    self%newvalue = returnval
    select type(anaVar =>  analyCommon(self%analyID)%val)
      type is(real(dp))
        anaVar = returnval
    end select

  end subroutine
!=========================================================================
  function PythonFunc_GetResult(self) result(var)
    implicit none
    class(PythonFunc), intent(in) :: self
    real(dp) :: var

    var = self%lastvalue
  end function
!=========================================================================
  subroutine PythonFunc_ProcessIO(self, line)
    !-------
    ! IO Format: python (python module)
    !------
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(PythonFunc), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0

    call GetXCommand(line, command, 2, lineStat)
    self%filename = command
   
!    call GetXCommand(line, command, 3, lineStat)
  end subroutine
!=========================================================================
  subroutine PythonFunc_CastCommonType(self, anaVar)
    implicit none
    class(PythonFunc), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
      write(*,*) "Allocated as Real"
    endif

  end subroutine
!=========================================================================
#endif
end module
!=========================================================================
