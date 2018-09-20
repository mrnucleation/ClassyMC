!=============================================================================+
module IntraBond_Harmonic
  use Template_IntraBond, only: Bond_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(Bond_FF) :: IntraBond_Harmonic
    real(dp) :: r0, k0
    contains
      procedure, pass :: Constructor 
      procedure, pass :: DetailedECalc 
      procedure, pass :: DiffECalc
      procedure, pass :: GenerateDist
      procedure, pass :: ProcessIO
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Bond_FF), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(Bond_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!============================================================================
  subroutine DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Bond_FF), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine GenerateDist(self, val, probgen)
    implicit none
    class(Bond_FF), intent(inout) :: self
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = 0E0_dp
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Bond_FF), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
end module
!=============================================================================+
