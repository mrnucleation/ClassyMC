!==============================================================
module AnalysisData
  use AnaylsisClassDef, only: Anaylsis

  type AnalysisArr
    class(Analysis), allocatable:: func
  end type

  type(AnalysisArr), allocatable, target  :: AnalysisArray(:)

end module
!==============================================================
