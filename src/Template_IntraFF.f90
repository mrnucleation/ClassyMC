!=============================================================================+
module Template_Intra_FF
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(classyClass) :: Intra_FF
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
    implicit none
    class(Intra_FF), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(Intra_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!============================================================================
  subroutine DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Intra_FF), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine GenerateDist(self, beta, val, probgen)
    implicit none
    class(Intra_FF), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = 0E0_dp
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Intra_FF), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
end module
!=============================================================================+
