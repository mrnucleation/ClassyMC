!==========================================================================================
! The purpose of this module is to provide the base class for the molecule
! constructor class.
!==========================================================================================
module MolCon_SimpleRegrowth
  use CoordinateTypes, only: Displacement, Perturbation
  use MasterTemplate, only: classyClass
  use VarPrecision


  type, public, extends(MolConstructor) :: SimpleRegrowth
    contains
      procedure, public, pass :: Constructor => SimpleRegrwoth_Constructor
      procedure, public, pass :: GenerateConfig => SimpleRegrwoth_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine SimpleRegrowth_Constructor(self)
    implicit none
    class(MolConstructor), intent(inout) :: self
  end subroutine
!==========================================================================================
  subroutine SimpleRegrowth_GenerateConfig(self, trialBox, disp, probconstruct)
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimpleBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct = 1E0_dp


  end subroutine
!==========================================================================================
end module
!==========================================================================================
