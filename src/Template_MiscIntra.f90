!=============================================================================+
module Template_IntraMiscIntra
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes

  type, public, extends(Intra_FF) :: MiscIntra_FF
    contains
      procedure, pass :: GenerateDist
  end type

  contains
!==========================================================================
  subroutine GenerateDist(self, beta, val, probgen, E_T)
    implicit none
    class(MiscIntra_FF), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T

    val = 0E0_dp
    E_T = 0E0_dp
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
end module
!=============================================================================+
