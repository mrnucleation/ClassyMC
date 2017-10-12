!====================================================================
module CommonSampling
  use VarPrecision
  use MetropolisRule, only: metropolis
  use AcceptRuleTemplate, only: acceptrule

  class(acceptrule), allocatable :: sampling

end module
!====================================================================
