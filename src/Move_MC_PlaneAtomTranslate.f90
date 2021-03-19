!========================================================
! Standard Molecular Translation move with the ability to restrict the
! degrees of motion by turning off movement in the x,y, or z directions. 
! This can be used to constraint to the system to 2D or 1D gemoetries. 
!========================================================
module MCMove_PlaneAtomTranslation

use MCMove_PlaneTranslation, only: PlaneTranslate
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(PlaneTranslate) :: PlaneAtomTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
!    real(dp), allocatable :: boxatmps(:)
!    real(dp) , allocatable:: boxaccpt(:)
!    logical :: proportional = .true.
!    logical :: tuneMax = .true.
!    logical :: xdir = .true.
!    logical :: ydir = .true.
!    logical :: zdir = .true.
!    real(dp) :: limit = 3.00E0_dp
!    real(dp) :: targAccpt = 50E0_dp
!    real(dp) :: max_dist = 0.05E0_dp
!    real(dp), allocatable :: boxlimit(:)
!    real(dp), allocatable :: boxtargAccpt(:)
!    real(dp), allocatable :: boxmax_dist(:)
!    type(Displacement), allocatable :: disp(:)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: FullMove => PlaneAtomTranslate_FullMove
      procedure, pass :: Prologue => PlaneAtomTranslate_Prologue
      procedure, pass :: Epilogue => PlaneAtomTranslate_Epilogue
  end type
!========================================================
 contains
!========================================================
!  subroutine PlaneAtomTranslate_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(PlaneAtomTranslate), intent(in) :: self
!    type(Displacement), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine PlaneAtomTranslate_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(PlaneAtomTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: boxID, iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: molIndx, molStart, molEnd, molType
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, E_Inter, E_Intra, biasE
    real(dp), parameter :: Prob = 1E0_dp

    boxID = trialBox % boxID
    self % atmps = self % atmps + 1E0_dp
    self % boxatmps(boxID) = self % boxatmps(boxID) + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp )
    call FindAtom(trialbox, rawIndx, nMove)
    call trialBox % GetAtomData(nMove, molIndx=molIndx)
    call trialBox % GetMolData(molIndx, molStart=molStart, molEnd=molEnd, &
                               molType=molType)

    dx = 0E0_dp
    dy = 0E0_dp
    dz = 0E0_dp

    if(self%xdir) then
      dx = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    endif
    if(self%ydir) then
      dy = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    endif
    if(self%zdir) then
      dz = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    endif
 
    self%disp(1)%molType = molType
    self%disp(1)%molIndx = molIndx
    self%disp(1)%atmIndx = nMove

    self%disp(1)%x_new = trialBox%atoms(1, nMove) + dx
    self%disp(1)%y_new = trialBox%atoms(2, nMove) + dy
    self%disp(1)%z_new = trialBox%atoms(3, nMove) + dz

    self%disp(1)%newlist = .false.
    self%disp(1)%listIndex = 1

    !If the particle moved a large distance get a temporary neighborlist
!    if(any([dx,dy,dz] > neighSkin)) then
!      call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%disp(1))
!      self%disp(1)%newlist = .true.
!    else

!    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialBox%ComputeEnergyDelta(self%disp(1:1),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)

    if(.not. accept) then
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:1), E_Diff )
    if(.not. accept) then
      return
    endif


    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:1), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % boxaccpt(boxID) = self % boxaccpt(boxID) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine PlaneAtomTranslate_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(PlaneAtomTranslate), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      

    write(nout,"(1x,A,F15.8)") "(Plane Atom Translate) Maximum Displacement: ", self%max_dist

  end subroutine
!=========================================================================
  subroutine PlaneAtomTranslate_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(PlaneAtomTranslate), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Plane Atom Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Plane Atom Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Plane Atom Translation Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,100F15.8)") "Final Maximum Displacement: ", self%boxmax_dist(1:)
    endif
 

  end subroutine
!========================================================
end module
!========================================================
