!==============================================================
module AnalysisData
  use AnaylsisClassDef, only: Analysis

  type AnalysisArr
    class(Analysis), allocatable:: func
  end type

  type(AnalysisArr), allocatable, target  :: AnalysisArray(:)

end module
!==============================================================
