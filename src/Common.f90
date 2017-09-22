!======================================================
      module CoordinateTypes
      use VarPrecision


      type Displacement
        integer(kind=atomIntType) :: molType, atmIndx
        real(dp) :: x_new, y_new, z_new
      end type

      end module

!======================================================
      module EnergyCalculators
      use ForceFieldTemplate


      type(forcefield), allocatable :: ECalcs(:)

      end module
!======================================================

