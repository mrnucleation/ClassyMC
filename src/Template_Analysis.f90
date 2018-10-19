!=========================================================================
module AnaylsisClassDef
  use MasterTemplate, only: classyClass
  use CoordinateTypes, only: Perturbation
  use VarPrecision

  type, public, extends(classyClass) :: Analysis
    logical :: perMove = .false.
    logical :: usedInMove = .false.
    integer :: IOUnit = -1
    integer :: UpdateFreq = -1
    integer :: analyID = -1
    contains
      procedure, pass :: Initialize
      procedure, pass :: CalcNewState
      procedure, pass :: Compute
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO
      procedure, pass :: CastCommonType
      procedure, pass :: WriteInfo
      procedure, pass :: GetResult
      procedure, pass :: Finalize
  end type

 contains
!=========================================================================
  subroutine Initialize(self)
    implicit none
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine CalcNewState(self, disp, newVal)
    use CoordinateTypes, only: Perturbation
    implicit none
    class(Analysis), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    real(dp), intent(in), optional :: newVal
  end subroutine
!=========================================================================
  subroutine Compute(self, accept)
    implicit none
    class(Analysis), intent(inout) :: self
    logical, intent(in) :: accept
  end subroutine
!=========================================================================
!  subroutine Maintenance(self)
!    implicit none
!    class(Analysis), intent(inout) :: self
!  end subroutine
!=========================================================================
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Analysis), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
  end subroutine
!=========================================================================
  subroutine CastCommonType(self, anaVar)
    implicit none
    class(Analysis), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
      write(*,*) "Allocated as Real"
    endif

  end subroutine
!=========================================================================
  subroutine WriteInfo(self)
    implicit none
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
  function GetResult(self) result(var)
    implicit none
    class(Analysis), intent(in) :: self
    real(dp) :: var
  end function
!=========================================================================
  subroutine Finalize(self)
    implicit none
    class(Analysis), intent(inout) :: self
  end subroutine
!=========================================================================
end module
!=========================================================================
