!=============================================================================+
module IntraAngle_Ridgid
  use Template_IntraAngle, only: Angle_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  type, public, extends(Angle_FF) :: RidgidAngle
!    real(dp) :: r0
!    real(dp) :: k0
    contains
      procedure, pass :: Constructor => RidgidAngle_Constructor
      procedure, pass :: DetailedECalc => RidgidAngle_DetailedECalc
      procedure, pass :: DiffECalc => RidgidAngle_DiffECalc
      procedure, pass :: GenerateDist => RidgidAngle_GenerateDist
      procedure, pass :: ProcessIO => RidgidAngle_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine RidgidAngle_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RidgidAngle), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine RidgidAngle_DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(RidgidAngle), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    E_T = 0E0_dp
    accept = .true.
  end subroutine
!============================================================================
  subroutine RidgidAngle_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(RidgidAngle), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine RidgidAngle_GenerateDist(self, val, probgen)
    implicit none
    class(RidgidAngle), intent(inout) :: self
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = self%theta0
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  subroutine RidgidAngle_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(RidgidAngle), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%theta0
    self%theta0 = self%theta0*inAngUnit

    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the ridgid Angle style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
