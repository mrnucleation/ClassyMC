!========================================================
module MCMove_AtomTranslation

use MCMove_MolTranslation, only: MolTranslate
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  !In order to prevent having to recode this again and again
  !The atomtranslation will be treated a subset of the moltranslation
  !subroutine. 
  type, public, extends(MolTranslate) :: AtomTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

!    real(dp), allocatable :: boxatmps(:)
!    real(dp) , allocatable:: boxaccpt(:)
!    logical :: verbose = .true.
!    logical :: proportional = .true.
!    logical :: tuneMax = .true.
!    real(dp) :: limit = 3.00E0_dp
!    real(dp) :: targAccpt = 50E0_dp
!    real(dp) :: max_dist = 0.05E0_dp
!    real(dp), allocatable :: boxlimit(:)
!    real(dp), allocatable :: boxtargAccpt(:)
!    real(dp), allocatable :: boxmax_dist(:)
!    type(Displacement), allocatable :: disp(:)

    !Rejection Counters
!    integer :: ovlaprej = 0 
!    integer :: constrainrej = 0 
!    integer :: detailedrej = 0 

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)


    contains
!      procedure, pass :: Constructor => AtomTrans_Constructor
!      procedure, pass :: GeneratePosition => AtomTrans_GeneratePosition
      procedure, pass :: FullMove => AtomTrans_FullMove
!      procedure, pass :: Maintenance => AtomTrans_Maintenance
!      procedure, pass :: Prologue => AtomTrans_Prologue
      procedure, pass :: Epilogue => AtomTrans_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine AtomTrans_GeneratePosition(self, disp)
    use RandomGen, only: grnd
    implicit none
    class(AtomTranslate), intent(in) :: self
    type(Displacement), intent(inout) :: disp
    real(dp) :: dx, dy, dz
      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
  end subroutine
!===============================================
  subroutine AtomTrans_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Box_Utility, only: FindAtom
    implicit none
    class(AtomTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: boxID, molIndx
    integer :: nMove, rawIndx
    integer :: molStart, molEnd, molType
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, E_Inter, E_Intra
    real(dp), parameter :: Prob = 1E0_dp

    boxID = trialBox % boxID
    call self%LoadBoxInfo(trialBox, self%disp)
    self % atmps = self % atmps + 1E0_dp
    self % boxatmps(boxID) = self % boxatmps(boxID) + 1E0_dp
    accept = .true.

    !Propose move by randomly selecting one atom in the system.
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp )
    
    call FindAtom(trialbox, rawIndx, nMove)
    call trialBox % GetAtomData(nMove, molIndx=molIndx)
    call trialBox % GetMolData(molIndx, molStart=molStart, molEnd=molEnd, &
                               molType=molType)

    dx = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    dy = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    dz = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
 

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
      self%constrainrej = self%constrainrej + 1
      return
    endif

    !Energy Calculation
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:nAtoms), self%tempList, self%tempNNei, E_Diff, accept)
    call trialBox%ComputeEnergyDelta(self%disp(1:1),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:1), E_Diff )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif


    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:1), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % boxaccpt(boxID) = self % boxaccpt(boxID) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList, self%tempNNei)
    else
      self%detailedrej = self%detailedrej + 1
!      write(*,*) E_Diff, trialBox%beta, Prob
    endif


  end subroutine
!=========================================================================
  subroutine AtomTrans_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(AtomTranslate), intent(inout) :: self
    real(dp) :: accptRate
 

    write(nout,*) 
    write(nout,"(1x,A,I15)") "Atom Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Atom Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Atom Translation Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,100F15.8)") "Final Maximum Displacement: ", self%boxmax_dist(1:)
    endif

    if(self%verbose) then
      write(nout, "(1x,A,I15)") "Atom Translation, Rejections due to overlap:", self%ovlaprej
      write(nout, "(1x,A,I15)") "Atom Translation, Rejections due to constraint:", self%constrainrej
      write(nout, "(1x,A,I15)") "Atom Translation, Rejections due to detailed balance:", self%detailedrej
    endif



  end subroutine
!========================================================
end module
!========================================================
