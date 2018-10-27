!=========================================================================
module Anaylsis_DistPair
use AnaylsisClassDef, only: Analysis
use AnalysisData, only: analyCommon
use VarPrecision
use SimpleSimBox, only: SimpleBox

  type, public, extends(Analysis):: DistPair
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
      procedure, pass :: Prologue => DistPair_Prologue
      procedure, pass :: Compute => DistPair_Compute
      procedure, pass :: CalcNewState => DistPair_CalcNewState
!      procedure, pass :: Maintenance 
      procedure, pass :: ProcessIO => DistPair_ProcessIO
      procedure, pass :: WriteInfo => DistPair_WriteInfo
      procedure, pass :: GetResult => DistPair_GetResult
      procedure, pass :: CastCommonType => DistPair_CastCommonType
!      procedure, pass :: Finalize => DistPair_Finalize
  end type

 contains
!=========================================================================
  subroutine DistPair_Prologue(self)
    use BoxData, only: BoxArray
    implicit none
    class(DistPair), intent(inout) :: self
    logical :: accept

    self%box => BoxArray(self%boxNum)%box

    self%perMove = .true.
    accept = .true.
    call self%Compute(accept)

  end subroutine
!=========================================================================
  subroutine DistPair_Compute(self, accept)
    implicit none
    class(DistPair), intent(inout) :: self
    logical, intent(in) :: accept

    real(dp) :: rx, ry, rz, rsq, r

!    accept = .true.
    if(.not. accept) then
      return
    endif

    rx = self%box % atoms(1, self%atom1) - self%box % atoms(1, self%atom2)
    ry = self%box % atoms(2, self%atom1) - self%box % atoms(2, self%atom2)
    rz = self%box % atoms(3, self%atom1) - self%box % atoms(3, self%atom2)
    call self%box% Boundary(rx, ry, rz)
    rsq = rx*rx + ry*ry + rz*rz
    r = sqrt(rsq)
 
!    write(*,*) r
    self%dist = r

  end subroutine
!=========================================================================
  subroutine DistPair_CalcNewState(self, disp, newVal)
    use CoordinateTypes, only: Displacement, Perturbation
    implicit none
    class(DistPair), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    integer :: iDisp
    real(dp), intent(in), optional :: newVal
    real(dp) :: rx, ry, rz, rsq, r

    r = self%dist
    select type(disp)
      class is(Displacement)
        do iDisp = 1, size(disp)
          if(disp(iDisp)%atmIndx == self%atom1 ) then
            rx = disp(iDisp)%x_new - self%box % atoms(1, self%atom2)
            ry = disp(iDisp)%y_new - self%box % atoms(2, self%atom2)
            rz = disp(iDisp)%z_new - self%box % atoms(3, self%atom2)
            call self%box% Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            r = sqrt(rsq)
           
          elseif(disp(iDisp)%atmIndx == self%atom2 ) then
            rx = disp(iDisp)%x_new - self%box % atoms(1, self%atom1)
            ry = disp(iDisp)%y_new - self%box % atoms(2, self%atom1)
            rz = disp(iDisp)%z_new - self%box % atoms(3, self%atom1)
            call self%box% Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
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
  subroutine DistPair_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(DistPair), intent(inout) :: self
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
  subroutine DistPair_WriteInfo(self)
    use ParallelVar, only: nout
    implicit none
    class(DistPair), intent(inout) :: self

  end subroutine
!=========================================================================
  function DistPair_GetResult(self) result(var)
    implicit none
    class(DistPair), intent(in) :: self
    logical :: accept
    real(dp) :: var

    var = self%dist
  end function
!=========================================================================
  subroutine DistPair_CastCommonType(self, anaVar)
    implicit none
    class(DistPair), intent(inout) :: self
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
