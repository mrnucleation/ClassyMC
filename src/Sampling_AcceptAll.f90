!====================================================================
! As it's name suggests, this sampling style accepts every move. It is useful
! for generating random configurations for training things such as NeuroNets,
! forcefield training, etc.  All other constraints, energy calculations, etc.
! are still used with this sampling style on.
!
! It's like the Cleveland Brown's defense, everything gets past it. But unlike
! the Brown's defense it's actually useful. 
!====================================================================
module AcceptAllRule
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: AcceptAll
    contains
       procedure, pass :: MakeDecision => AcceptAll_MakeDecision
!       procedure, pass :: Maintenance => AcceptAll_Maintenance
!       procedure, pass :: ProcessIO => AcceptAll_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function AcceptAll_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(AcceptAll), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional:: inProb, logProb
    real(dp), intent(in) :: E_Diff
    logical :: accept
    
    accept = .true.
  end function
!====================================================================
end module
!====================================================================
