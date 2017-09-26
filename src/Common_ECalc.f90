!======================================================
module ForcefieldData
use ForceFieldTemplate, only: ForceField

  type ECalcArray 
    type(forcefield), allocatable :: Method
  end type

  type(ECalcArray), allocatable, target :: EnergyCalculator(:)

end module
!======================================================
