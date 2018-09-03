!==============================================================
module ScreenData
  use AnaylsisClassDef, only: Analysis
  use VarPrecision

  type ScreenVals
    class(*), allocatable:: val
  end type


!  real(dp), allocatable :: analyCommon(:)
  type(ScreenVals), allocatable :: screenCommon(:)

end module
!==============================================================
