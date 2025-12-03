!======================================================
module ForcefieldData
use Template_ForceField, only: ForceField

  integer :: nForceFields = 0
  type ECalcArray 
    class(forcefield), allocatable :: Method
  end type

  type(ECalcArray), allocatable, target :: EnergyCalculator(:)

end module
!======================================================
