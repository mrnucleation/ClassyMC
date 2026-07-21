!=========================================================================
module Anaylsis_PsuedoDistPair
use AnaylsisClassDef, only: Analysis
use AnalysisData, only: analyCommon
use VarPrecision
use SimpleSimBox, only: SimpleBox

  type, public, extends(Analysis):: PsuedoDistPair
!    logical :: perMove = .false.
!    logical :: usedInMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1
!    integer :: analyID = -1

    integer :: boxNum = 1
    integer :: atom1, atom2
    real(dp) :: dist
    class(SimpleBox), pointer :: box => null()

    contains
!      procedure, pass :: Initialize
      procedure, pass :: Prologue => PsuedoDistPair_Prologue
      procedure, pass :: Compute => PsuedoDistPair_Compute
      procedure, pass :: CalcNewState => PsuedoDistPair_CalcNewState
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => PsuedoDistPair_ProcessIO
      procedure, pass :: WriteInfo => PsuedoDistPair_WriteInfo
      procedure, pass :: GetResult => PsuedoDistPair_GetResult
      procedure, pass :: CastCommonType => PsuedoDistPair_CastCommonType
!      procedure, pass :: Finalize => PsuedoDistPair_Finalize
  end type

 contains
!=========================================================================
  subroutine PsuedoDistPair_Prologue(self)
    use BoxData, only: BoxArray
    implicit none
    class(PsuedoDistPair), intent(inout) :: self
    logical :: accept

    self%box => BoxArray(self%boxNum)%box

    self%perMove = .true.
    accept = .true.
    call self%Compute(accept)

  end subroutine
!=========================================================================
  subroutine PsuedoDistPair_Compute(self, accept)
    implicit none
    class(PsuedoDistPair), intent(inout) :: self
    logical, intent(in) :: accept

    real(dp) :: rx, rz, rsq, r

!    accept = .true.
    if(.not. accept) then
      return
    endif

    rx = self%box % atoms(1, self%atom1) - self%box % atoms(1, self%atom2)
    rz = self%box % atoms(3, self%atom1) - self%box % atoms(3, self%atom2)
    call self%box% Boundary(rx=rx, rz=rz)
    rsq = rx*rx + rz*rz
    r = sqrt(rsq)
 
!    write(*,*) r
    self%dist = r

  end subroutine
!=========================================================================
  subroutine PsuedoDistPair_CalcNewState(self, disp, accept, newVal)
    use CoordinateTypes, only: Displacement, Perturbation
    implicit none
    class(PsuedoDistPair), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    integer :: iDisp
    real(dp), intent(in), optional :: newVal
    logical, intent(out) :: accept
    real(dp) :: rx, rz, rsq, r

    accept = .true.
    r = self%dist
    select type(disp)
      class is(Displacement)
        do iDisp = 1, size(disp)
          if(disp(iDisp)%atmIndx == self%atom1 ) then
            rx = disp(iDisp)%x_new - self%box % atoms(1, self%atom2)
            rz = disp(iDisp)%z_new - self%box % atoms(3, self%atom2)
            call self%box% Boundary(rx=rx, rz=rz)
            rsq = rx*rx + rz*rz
            r = sqrt(rsq)
           
          elseif(disp(iDisp)%atmIndx == self%atom2 ) then
            rx = disp(iDisp)%x_new - self%box % atoms(1, self%atom1)
            rz = disp(iDisp)%z_new - self%box % atoms(3, self%atom1)
            call self%box% Boundary(rx=rx, rz=rz)
            rsq = rx*rx + rz*rz
            r = sqrt(rsq)
          endif
       enddo
    end select

!    write(*,*) r
    select type(anaVar =>  analyCommon(self%analyID)%val)
      type is(real(dp))
        anaVar = r 
    end select

  end subroutine
!=========================================================================
  subroutine PsuedoDistPair_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(PsuedoDistPair), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%atom1 = intVal

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) intVal
    self%atom2 = intVal

   

  end subroutine
!=========================================================================
  subroutine PsuedoDistPair_WriteInfo(self)
    use ParallelVar, only: nout
    implicit none
    class(PsuedoDistPair), intent(inout) :: self

  end subroutine
!=========================================================================
  function PsuedoDistPair_GetResult(self) result(var)
    implicit none
    class(PsuedoDistPair), intent(in) :: self
    logical :: accept
    real(dp) :: var

    var = self%dist
  end function
!=========================================================================
  subroutine PsuedoDistPair_CastCommonType(self, anaVar)
    implicit none
    class(PsuedoDistPair), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
      write(*,*) "Allocated as Real"
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
