!=============================================================================+
module IntraBond_Ridgid
  use Template_IntraBond, only: Bond_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(Bond_FF) :: RidgidBond
!    real(dp) :: r0
    real(dp) :: k0
    contains
      procedure, pass :: Constructor => RidgidBond_Constructor
      procedure, pass :: DetailedECalc => RidgidBond_DetailedECalc
      procedure, pass :: DiffECalc => RidgidBond_DiffECalc
      procedure, pass :: GenerateDist => RidgidBond_GenerateDist
      procedure, pass :: ProcessIO => RidgidBond_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine RidgidBond_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RidgidBond), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine RidgidBond_DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(RidgidBond), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    E_T = 0E0_dp
    accept = .true.
  end subroutine
!============================================================================
  subroutine RidgidBond_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(RidgidBond), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine RidgidBond_GenerateDist(self, val, probgen)
    implicit none
    class(RidgidBond), intent(inout) :: self
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = self%r0
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  subroutine RidgidBond_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(RidgidBond), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%r0

    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the ridgid bond style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
