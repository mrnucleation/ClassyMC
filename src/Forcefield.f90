module ForceFieldTemplate
  use VarPrecision
  type, public :: forcefield
    contains
      procedure, pass :: DetailedECalc 
      procedure, pass :: ShiftECalc
      procedure, pass :: SwapInECalc
      procedure, pass :: SwapOutECalc
      procedure, pass :: SetParameter
  end type

  contains

  subroutine DetailedECalc(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine

  subroutine ShiftCheck(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine

  subroutine SwapInCheck(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine

  subroutine SwapOutCheck(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine


  subroutine SetParameter(self, parIndex,  parVal)
    implicit none
    class(forcefield), intent(in) :: self
    integer, intent(in) :: parIndex
    real(dp), intent(in) :: parVal
  end subroutine


end module
