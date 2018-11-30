!==========================================================================================
! The purpose of this module is to provide the base class for the molecule
! constructor class.
!==========================================================================================
module Template_MolConstructor
  use CoordinateTypes, only: Perturbation
  use MasterTemplate, only: classyClass
!  use SimpleSimBox, only: SimpleBox
  use Template_SimBox, only: SimBox
  use VarPrecision


  type, public, extends(classyClass) :: MolConstructor
    integer :: insPoints = 1
    integer :: molType = -1
    contains
      procedure, public, pass :: Constructor
      procedure, public, pass :: SetMolType
      procedure, public, pass :: GenerateConfig
      procedure, public, pass :: ReverseConfig
      procedure, public, pass :: GasConfig
      procedure, public, pass :: GetNInsertPoints
      procedure, public, pass :: ProcessIO
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine Constructor(self)
    implicit none
    class(MolConstructor), intent(inout) :: self
!    integer, intent(in) :: molType
  end subroutine

!==========================================================================================
  subroutine SetMolType(self, molType)
    implicit none
    class(MolConstructor), intent(inout) :: self
    integer, intent(in) :: molType

    self % molType = molType
  end subroutine
!==========================================================================================
  subroutine GenerateConfig(self, trialBox, disp, probconstruct, insPoint)
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)
    real(dp), intent(out) :: probconstruct 

    probconstruct = 1E0_dp
  end subroutine
!==========================================================================================
  subroutine ReverseConfig(self, trialBox, probconstruct, accept)
    implicit none
    class(MolConstructor), intent(inout) :: self
!    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct 
    logical, intent(out) :: accept

    accept = .true.
    probconstruct = 1E0_dp
  end subroutine
!==========================================================================================
  subroutine GasConfig(self, molType, probGas)
    implicit none
    class(MolConstructor), intent(inout) :: self
    integer, intent(in) :: molType
    real(dp), intent(out) :: probGas

    probGas = 1E0_dp
  end subroutine
!==========================================================================================
  function GetNInsertPoints(self) result(nPoints)
    implicit none
    class(MolConstructor), intent(inout) :: self
    integer :: nPoints

    nPoints = self%insPoints

  end function
!==========================================================================================
  subroutine ProcessIO(self, line, linestat)
    implicit none
    class(MolConstructor), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: linestat

    linestat = 0
  end subroutine
!==========================================================================================
end module
!==========================================================================================
