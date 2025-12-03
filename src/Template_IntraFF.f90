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
      procedure, pass :: GenerateTrial
      procedure, pass :: GenerateDist
      procedure, pass :: GenerateReverseDist
!      procedure, pass :: ComputeProb
      procedure, pass :: ProcessIO
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    implicit none
    class(Intra_FF), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, atompos, E_T, accept)
    implicit none
    class(Intra_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    real(dp), intent(in) :: atompos(:, :)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!==========================================================================
  subroutine GenerateTrial(self, beta, val, bounds)
    implicit none
    class(Intra_FF), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(in), optional :: bounds(1:2)
   
    val = 0E0_dp
  end subroutine
!==========================================================================
  subroutine GenerateDist(self, beta, val, probgen, E_T)
    implicit none
    class(Intra_FF), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T

    val = 0E0_dp
    probgen = 1E0_dp
    E_T = 0E0_dp

  end subroutine
!==========================================================================
  subroutine GenerateReverseDist(self, curbox, atompos, probgen)
    implicit none
    class(Intra_FF), intent(inout) :: self
    class(SimBox), intent(inout) :: curBox
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(out) :: probgen


    probgen = 1E0_dp

  end subroutine
!==========================================================================
  function ComputeProb(self, beta, val) result(probgen)
    implicit none
    class(Intra_FF), intent(in) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: val
    real(dp) :: probgen
!    real(dp) :: E_Val

!    E_Val = self%EFunc(val)
!    probgen = exp(-beta*E_Val)
  end function
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
