!==========================================================================================
! The purpose of this module is to provide the base class for the molecule
! constructor class.
!==========================================================================================
module Template_MolConstructor
  use CoordinateTypes, only: Displacement, Perturbation
  use MasterTemplate, only: classyClass
!  use SimpleSimBox, only: SimpleBox
  use Template_SimBox, only: SimBox
  use VarPrecision


  type, public, extends(classyClass) :: MolConstructor
    contains
      procedure, public, pass :: Constructor
      procedure, public, pass :: GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine Constructor(self, molType)
    implicit none
    class(MolConstructor), intent(inout) :: self
    integer, intent(in) :: molType
  end subroutine
!==========================================================================================
  subroutine GenerateConfig(self, trialBox, disp, probconstruct)
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct 

    probconstruct = 1E0_dp
  end subroutine
!==========================================================================================
end module
!==========================================================================================
