!==========================================================================================
! The purpose of this module is to provide the base class for the simulation box
! family of objects. This type is only intended to provide the basic
! structure for child classes and should not be used directly. 
!==========================================================================================
module Template_MolConstructor
  use CoordinateTypes, only: Displacement, Perturbation
  use MasterTemplate, only: classyClass
  use VarPrecision


  !Sim Box Definition
  type, public, extends(classyClass) :: MolConstructor
    contains
      procedure, public, pass :: Constructor
      procedure, public, pass :: GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine Constructor(self)
    implicit none
    class(MolConstructor), intent(inout) :: self
  end subroutine
!==========================================================================================
  subroutine GenerateConfig(self, disp)
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
  end subroutine
!==========================================================================================
end module
!==========================================================================================
