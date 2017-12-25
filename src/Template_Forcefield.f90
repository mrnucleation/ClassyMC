module Template_ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public :: forcefield
    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor 
      procedure, pass :: DetailedECalc 
      procedure, pass :: NewECalc
      procedure, pass :: ShiftECalc_Single
      procedure, pass :: ShiftECalc_Multi
      procedure, pass :: SwapInECalc
      procedure, pass :: SwapOutECalc
      procedure, pass :: ExchangeECalc
      procedure, pass :: ProcessIO
      procedure, pass :: GetCutOff
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(forcefield), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, E_T)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T

  end subroutine
!=============================================================================+
  subroutine NewECalc(self, curbox, disp, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff
  end subroutine
!=============================================================================+
  subroutine ShiftECalc_Single(self, curbox, disp, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff
  end subroutine
!=============================================================================+
  subroutine ShiftECalc_Multi(self, curbox, disp, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
  end subroutine
!=============================================================================+
  subroutine SwapInECalc(self, curbox, disp, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff

  end subroutine
!=============================================================================+
  subroutine SwapOutECalc(self, curbox, atmIndx, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      real(dp), intent(inOut) :: E_Diff
      integer, intent(in) :: atmIndx(:)
  end subroutine
!=============================================================================+
  subroutine ExchangeECalc(self, curbox, atmIndx, newType, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      real(dp), intent(inOut) :: E_Diff
      integer, intent(in) :: atmIndx, newType
  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    implicit none
    class(forcefield), intent(inout) :: self
    character(len=*), intent(in) :: line

  end subroutine
!=============================================================================+
  function GetCutOff(self) result(rCut)
    implicit none
    class(forcefield), intent(inout) :: self
    real(dp) :: rCut

    write(*,*) self%rCut
    rCut = self%rCut
  end function
!=============================================================================+
end module
!=============================================================================+
