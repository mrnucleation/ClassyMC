!==============================================================
module AnalysisData
  use AnaylsisClassDef, only: Analysis
  use VarPrecision

  type AnalysisArr
    class(Analysis), allocatable:: func
  end type

  type(AnalysisArr), allocatable, target  :: AnalysisArray(:)

  real(dp), allocatable :: analyCommon(:)

end module
!==============================================================
