!=============================================================================+
module IntraBond_Harmonic
  use Template_IntraBond, only: Bond_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(Bond_FF) :: HarmonicBond
!    real(dp) :: r0
    real(dp) :: k0
    contains
      procedure, pass :: Constructor => HarmonicBond_Constructor
      procedure, pass :: DetailedECalc => HarmonicBond_DetailedECalc
      procedure, pass :: DiffECalc => HarmonicBond_DiffECalc
      procedure, pass :: GenerateDist => HarmonicBond_GenerateDist
      procedure, pass :: ProcessIO => HarmonicBond_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine HarmonicBond_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(HarmonicBond), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine HarmonicBond_DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(HarmonicBond), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!============================================================================
  subroutine HarmonicBond_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(HarmonicBond), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine HarmonicBond_GenerateDist(self, val, probgen)
    implicit none
    class(HarmonicBond), intent(inout) :: self
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = 0E0_dp
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  subroutine HarmonicBond_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(HarmonicBond), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%r0

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) self%k0

    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the harmonic bond style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
