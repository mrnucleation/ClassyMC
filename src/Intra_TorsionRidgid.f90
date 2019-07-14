!=============================================================================+
! Ridgid Torsion Angle Forcefield.  Not sure why you might need it, but it's
! here if you do.  Fixes a Torsion angle to an exact value. May not work
! well with Improper Torsion Angles.
!=============================================================================+
module IntraTorsion_Ridgid
  use Template_IntraTorsion, only: Torsion_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  type, public, extends(Torsion_FF) :: RidgidTorsion
    real(dp) :: ang0
    contains
      procedure, pass :: Constructor => RidgidTorsion_Constructor
      procedure, pass :: DetailedECalc => RidgidTorsion_DetailedECalc
      procedure, pass :: DiffECalc => RidgidTorsion_DiffECalc
      procedure, pass :: GenerateDist => RidgidTorsion_GenerateDist
      procedure, pass :: ProcessIO => RidgidTorsion_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine RidgidTorsion_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RidgidTorsion), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine RidgidTorsion_DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(RidgidTorsion), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    E_T = 0E0_dp
    accept = .true.
  end subroutine
!============================================================================
  subroutine RidgidTorsion_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(RidgidTorsion), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine RidgidTorsion_GenerateDist(self, beta, val, probgen)
    implicit none
    class(RidgidTorsion), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = self%ang0
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  subroutine RidgidTorsion_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(RidgidTorsion), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%ang0
    self%ang0 = self%ang0*inAngUnit

    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the ridgid Torsion style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
