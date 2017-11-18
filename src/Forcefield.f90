module ForceFieldTemplate
  use VarPrecision
  use SimBoxDef, only: SimBox
  use CoordinateTypes

  type, public :: forcefield
    contains
      procedure, pass :: Constructor 
      procedure, pass :: DetailedECalc 
      procedure, pass :: ShiftECalc_Single
      procedure, pass :: ShiftECalc_Multi
      procedure, pass :: SwapInECalc
      procedure, pass :: SwapOutECalc
      procedure, pass :: Exchange
      procedure, pass :: SetParameter
      procedure, pass :: ReadParFile
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    use Common_MolDef, only: nMolTypes
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
  subroutine SwapOutECalc(self, curbox, atmIndx, E_Diff)
    implicit none
      class(forcefield), intent(in) :: self
      class(simBox), intent(inout) :: curbox
      real(dp), intent(inOut) :: E_Diff
      integer, intent(in) :: atmIndx(:)
  end subroutine
!=============================================================================+
  subroutine SetParameter(self, parIndex,  parVal)
    implicit none
    class(forcefield), intent(inout) :: self
    integer, intent(in) :: parIndex(:)
    real(dp), intent(in) :: parVal
  end subroutine
!=============================================================================+
  subroutine ReadParFile(self, fileName)
    implicit none
    class(forcefield), intent(inout) :: self
    character(len=*), intent(in) :: fileName
  end subroutine
!=============================================================================+

end module
