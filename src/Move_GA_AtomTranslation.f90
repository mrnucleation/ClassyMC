!========================================================
module GAMove_AtomTranslation
use CoordinateTypes, only: Displacement
use MultiBoxMoveDef, only: MCMultiBoxMove
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMultiBoxMove) :: GA_AtomTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: max_dist = 0.05E0_dp

    integer :: movedAtoms = 1

    contains
      procedure, pass :: Constructor => AtomTrans_Constructor
!      procedure, pass :: FullMove => AtomTrans_FullMove
      procedure, pass :: MultiBox => AtomTrans_MultiBox
      procedure, pass :: Maintenance => AtomTrans_Maintenance
!      procedure, pass :: Prologue => AtomTrans_Prologue
!      procedure, pass :: Epilogue => AtomTrans_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine AtomTrans_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(GA_AtomTranslate), intent(inout) :: self


  end subroutine
!========================================================
  subroutine MultiBox(self, trialBox)
    class(GA_AtomTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox(:)

  end subroutine
!===============================================
  subroutine AtomTrans_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use BoxData, only: BoxArray
    use Box_Utility, only: FindAtom
    implicit none
    class(GA_AtomTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex
    integer :: parent, child
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp)
    call FindAtom(trialbox, rawIndx, nMove)
    dx = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dy = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dz = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
 

    !Check Constraint
!    accept = trialBox % CheckConstraint( self%disp(1:1) )
!    if(.not. accept) then
!      return
!    endif

    BoxArray(1) % box =  BoxArray(2) % box
    !Energy Calculation
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff, accept)
!    if(.not. accept) then
!      return
!    endif


    !Accept/Reject
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
    endif

  end subroutine
!=========================================================================
  subroutine AtomTrans_Maintenance(self)
    implicit none
    class(GA_AtomTranslate), intent(inout) :: self
    real(dp), parameter :: limit = 3.0E0_dp
      
 

  end subroutine
!=========================================================================
  subroutine AtomTrans_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(GA_AtomTranslate), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Atom Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Atom Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Atom Translation Acceptance Rate: ", accptRate
 

  end subroutine
!========================================================
end module
!========================================================
