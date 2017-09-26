!======================================================
module ForcefieldData
use ForceFieldTemplate, only: ForceField

  integer :: nForceFields = 0
  type ECalcArray 
    class(forcefield), allocatable :: Method
  end type

  type(ECalcArray), allocatable, target :: EnergyCalculator(:)

end module
!======================================================
