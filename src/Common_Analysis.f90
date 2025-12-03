!==============================================================
module AnalysisData
  use AnaylsisClassDef, only: Analysis
  use VarPrecision

  type AnalysisArr
    class(Analysis), allocatable:: func
  end type

  type(AnalysisArr), allocatable, target  :: AnalysisArray(:)

  type AnalysisVals
    class(*), allocatable:: val
  end type


!  real(dp), allocatable :: analyCommon(:)
  type(AnalysisVals), allocatable :: analyCommon(:)

end module
!==============================================================
