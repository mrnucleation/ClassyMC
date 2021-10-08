!=========================================================================
! Analyzes the number of molecules of a specific type in the cluster
! Primarily used to bias the addition algorithm in conjunction with
! umbrella sampling to compute the free energy as a function of
! molecule count.
!=========================================================================
module Anaylsis_TotalSize
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: TotalSize
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    integer, private :: boxNum = -1
    integer, private :: nMol 
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Prologue => TotalSize_Prologue
      procedure, pass :: Compute => TotalSize_Compute
      procedure, pass :: CalcNewState => TotalSize_CalcNewState
      procedure, pass :: CastCommonType => TotalSize_CastCommonType

!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => TotalSize_ProcessIO
      procedure, pass :: GetResult => TotalSize_GetResult
  end type

 contains
!=========================================================================
  subroutine TotalSize_Prologue(self)
    use BoxData, only: BoxArray
    implicit none
    class(TotalSize), intent(inout) :: self

    self%UpdateFreq = 1


  end subroutine
!=========================================================================
  subroutine TotalSize_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(TotalSize), intent(inout) :: self
    logical, intent(in) :: accept

    self%nMol = BoxArray(self%boxNum) % box % nMolTotal

  end subroutine
!=========================================================================
  function TotalSize_GetResult(self) result(var)
    implicit none
    class(TotalSize), intent(in) :: self
    real(dp) :: var

    var = self%nMol
  end function
!=========================================================================
  subroutine TotalSize_CalcNewState(self, disp, newVal)
    use AnalysisData, only: analyCommon
    use CoordinateTypes, only: Perturbation, Deletion, Addition, AtomExchange
    implicit none
    class(TotalSize), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    real(dp), intent(in), optional :: newVal
    integer :: Diff
    integer :: molNew, molOld
    integer :: typeNew, typeOld

    Diff = 0
    select type(disp)
      class is(Deletion)
        Diff = Diff - 1

      class is(Addition)
        Diff = Diff + 1

    end select

    select type(anaVar => analyCommon(self%analyID)%val)
      type is(integer)
        anaVar = self%nMol + Diff
    end select
  end subroutine
!=========================================================================
  subroutine TotalSize_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(TotalSize), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal

    !Format = (TotalSize) (Box Number) (Mol Type)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal

  end subroutine
!=========================================================================
  subroutine TotalSize_CastCommonType(self, anaVar)
    implicit none
    class(TotalSize), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    integer :: def

    def = 0
    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
