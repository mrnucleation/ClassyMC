!====================================================================
module CommonSampling
  use VarPrecision
  use MetropolisRule, only: metropolis
  use Template_AcceptRule, only: acceptrule

  class(acceptrule), allocatable :: sampling

end module
!====================================================================
