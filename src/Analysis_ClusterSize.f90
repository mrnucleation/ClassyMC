!=========================================================================
module Anaylsis_ClusterSize
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: ClusterSize
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    integer :: boxNum = -1
    integer :: molType = -1

    integer :: nMol 
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Prologue => ClusterSize_Prologue
      procedure, pass :: Compute => ClusterSize_Compute
      procedure, pass :: CalcNewState => ClusterSize_CalcNewState
      procedure, pass :: CastCommonType => ClusterSize_CastCommonType

!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => ClusterSize_ProcessIO
      procedure, pass :: GetResult => ClusterSize_GetResult
  end type

 contains
!=========================================================================
  subroutine ClusterSize_Prologue(self)
    use BoxData, only: BoxArray
    implicit none
    class(ClusterSize), intent(inout) :: self

    self%UpdateFreq = 1


  end subroutine
!=========================================================================
  subroutine ClusterSize_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(ClusterSize), intent(inout) :: self
    logical, intent(in) :: accept

    self%nMol = BoxArray(self%boxNum) % box % NMol(self%molType)

  end subroutine
!=========================================================================
  function ClusterSize_GetResult(self) result(var)
    implicit none
    class(ClusterSize), intent(in) :: self
    real(dp) :: var

    var = self%nMol
  end function
!=========================================================================
  subroutine ClusterSize_CalcNewState(self, disp, newVal)
    use AnalysisData, only: analyCommon
    use CoordinateTypes, only: Perturbation, Deletion, Addition
    implicit none
    class(ClusterSize), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    real(dp), intent(in), optional :: newVal
    integer :: Diff

    Diff = 0
    select type(disp)
      class is(Deletion)
          if(disp(1)%MolType == self%molType) then
            Diff = Diff - 1
          endif

      class is(Addition)
          if(disp(1)%MolType == self%molType) then
            Diff = Diff + 1
          endif
        
    end select
!    write(*,*) "Newcalc", Diff, self%nMol

    select type(anaVar => analyCommon(self%analyID)%val)
      type is(integer)
!        analyCommon(self%analyID)%val = real(self%nMol + Diff, dp) + 1e-10
        anaVar = self%nMol + Diff
    end select
  end subroutine
!=========================================================================
  subroutine ClusterSize_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(ClusterSize), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal

    !Format = (ClusterSize) (Box Number) (Mol Type)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) intVal
    self%molType = intVal

  end subroutine
!=========================================================================
  subroutine ClusterSize_CastCommonType(self, anaVar)
    implicit none
    class(ClusterSize), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    integer :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
      write(*,*) "Allocated as Integer"
    endif

  end subroutine

!=========================================================================
end module
!=========================================================================
