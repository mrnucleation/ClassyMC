!========================================================
module MCMove_ParticleExchange
use CoordinateTypes, only: Addition, Deletion
use MultiBoxMoveDef, only: MCMultiBoxMove
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMultiBoxMove) :: ParticleExchange
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    type(Addition), allocatable, private :: newPart(:)
    type(Deletion), private :: oldPart(1:1)

    real(dp), private, allocatable :: insPoint(:, :)
    real(dp), private, allocatable :: insProb(:)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: SafetyCheck => ParticleExchange_SafetyCheck
      procedure, pass :: Constructor => ParticleExchange_Constructor
      procedure, pass :: MultiBox => ParticleExchange_MultiBox
      procedure, pass :: AllocateProb => ParticleExchange_AllocateProb
!      procedure, pass :: GeneratePosition => ParticleExchange_GeneratePosition
!      procedure, pass :: FullMove => ParticleExchange_FullMove
!      procedure, pass :: Maintenance => ParticleExchange_Maintenance
      procedure, pass :: Prologue => ParticleExchange_Prologue
      procedure, pass :: Epilogue => ParticleExchange_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine ParticleExchange_SafetyCheck(self)
    use BoxData, only: boxArray
    use SimpleSimBox, only: SimpleBox
    implicit none
    class(ParticleExchange), intent(inout) :: self
    integer :: iBox

    !Ensure a boundary condition is set.
    do iBox = 1, size(BoxArray)
      select type(box => BoxArray(iBox)%box)
        type is(SimpleBox)
            write(*,*) "WARNING! Basic Exchange move is not designed to work without"
            write(*,*) "a bounding box!"
            write(*,*) "Box Number:", iBox
      end select
    enddo



  end subroutine
!========================================================
  subroutine ParticleExchange_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(ParticleExchange), intent(inout) :: self
!    integer :: iType, maxAtoms



!    allocate( self%tempNNei(1) )
!    allocate( self%tempList(200, 1) )
  end subroutine
!========================================================
!  subroutine ParticleExchange_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(ParticleExchange), intent(in) :: self
!    type(DisplacementNew), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!========================================================
  subroutine  ParticleExchange_AllocateProb(self, nInsPoints)
    implicit none
    class(ParticleExchange), intent(inout) :: self
    integer, intent(in) :: nInsPoints

    if(.not. allocated(self%insPoint) ) then
      allocate(self%insPoint(1:3, 1:nInsPoints))
      allocate(self%insProb(1:nInsPoints))
    else if(size(self%insPoint, 2) < nInsPoints) then
      deallocate(self%insPoint)
      deallocate(self%insProb)
      allocate(self%insPoint(1:3, 1:nInsPoints))
      allocate(self%insProb(1:nInsPoints))
    endif

  end subroutine
!=========================================================================
  subroutine ParticleExchange_MultiBox(self, accept)
    use Box_Utility, only: FindMolecule
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, MolData
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use RandomGen, only: grnd, ListRNG
    use CommonSampling, only: sampling

    implicit none
    class(ParticleExchange), intent(inout) :: self
    logical, intent(out) :: accept
    integer :: i, iAtom, iPoint, boxNum, nAtoms
    integer :: rawIndx, atomIndx, nMoveOut, nMoveIn
    integer :: molType, molStart, molEnd
    integer :: nInsPoints
    class(SimpleBox), pointer :: box1, box2
    real(dp) :: Prob, ProbSubIn, ProbSubOut, Norm, extraTerms, half
    real(dp) :: vol
    real(dp) :: reduced(1:3)
    real(dp) :: E_Diff1, E_Diff2
    real(dp) :: E_Inter1, E_Intra1
    real(dp) :: E_Inter2, E_Intra2
    real(dp) :: rescale(1:size(self%boxprob))
    character(len=30), parameter :: volume="volume"



    !Randomly choose which boxes will exchange particles
    boxNum = ListRNG(self%boxProb)
    box1 => BoxArray(boxNum)%box

    !To avoid picking the same box twice, rescale the probability such that
    !box1's probability is equal to 0.
    rescale = self%boxprob
    rescale(boxNum) = 0E0_dp
    norm = 0E0_dp
    do i = 1, size(rescale)
      norm = norm + rescale(i)
    enddo
    if(norm == 0E0_dp) then
      write(0,*) "WARNING! To use ParticleExchange a nonzero probability must"
      write(0,*) "be set for more than one box!"
      stop 
    endif
    do i = 1, size(rescale)
      rescale(i) = rescale(i)/norm
    enddo
    boxNum = ListRNG(rescale)
    box2 => BoxArray(boxNum)%box


    !Increment attempt counter.
    self % atmps = self % atmps + 1E0_dp

    !Choose a particle to remove from Box1
    rawIndx = floor( box1%nMolTotal * grnd() + 1E0_dp)
    if(rawIndx > box1%nMolTotal) then
      rawIndx = box1%nMolTotal
    endif
    call FindMolecule(box1, rawIndx, nMoveOut)
    call box1 % GetMolData(nMoveOut, molType=molType, molStart=molStart)

    !Ensure Box1 is above it's minimum molecule limit
    if(box1%NMol(molType) - 1 < box1%NMolMin(molType)) then
      accept = .false.
      return
    endif
    self%oldPart(1)%molType = molType
    self%oldPart(1)%molIndx = nMoveOut

    !Ensure Box2 is below it's maximum molecule limit
    if(box2%NMol(molType) + 1 > box2%NMolMax(molType)) then
      accept = .false.
      return
    endif

    !Check Box1's Constraints
    accept = box1 % CheckConstraint( self%oldPart(1:1) )
    if(.not. accept) then
      return
    endif

    !Get the first empty molecule position for box2
    nMoveIn = box2%NMol(molType) + 1
    nMoveIn = box2%MolGlobalIndx(molType, nMoveIn)
    call box2 % GetMolData(nMoveIn, molStart=molStart, molEnd=molEnd)

    !Generate an Insertion Point for the new particle
    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    do iPoint = 1, nInsPoints
      call self%AllocateProb(nInsPoints)
      do i = 1, 3
        reduced(i) = grnd()
      enddo
      call box2%GetRealCoords(reduced, self%insPoint(1:3, iPoint))
      self%insProb(iPoint) = 1E0_dp
    enddo

    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      self%newPart(iAtom)%molType = molType
      self%newPart(iAtom)%molIndx = nMoveIn
      self%newPart(iAtom)%atmIndx = atomIndx
    enddo


    call MolData(molType) % molConstruct % GenerateConfig(box2, self%newPart(1:nAtoms), ProbSubIn, accept ,self%insPoint, self%insProb)
    if(.not. accept) then
      return
    endif

    do iAtom = 1, nAtoms
      call box2 % NeighList(1) % GetNewList(iAtom, self%tempList, self%tempNNei, &
                                                self%newPart(iAtom))
      self%newPart(iAtom)%listIndex = iAtom
    enddo 

    !Check Box 2's Constraints
    accept = box2 % CheckConstraint( self%newPart(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Box1's Energy Calculation
!    call box1 % EFunc % Method % DiffECalc(box1, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff1, accept)
    call box1%ComputeEnergyDelta(self%oldpart(1:1),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter1, &
                                     E_Intra1, &
                                     E_Diff1, &
                                     accept, &
                                     computeintra=.true.)
    if(.not. accept) then
!      write(*,*) "Energy Rejection"
      return
    endif

    !Box2's Energy Calculation
    call box2%ComputeEnergyDelta(self%newpart(1:nAtoms),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter2, &
                                     E_Intra2, &
                                     E_Diff2, &
                                     accept, &
                                     computeintra=.true.)
    if(.not. accept) then
      return
    endif


    !Compute the Reverse Generation Probability For Swapping the Same Particle Back into
    !It's Original Position in Box 1
    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    do iPoint = 1, nInsPoints
      call self%AllocateProb(nInsPoints)
      do i = 1, 3
        reduced(i) = grnd()
      enddo
      call box1%GetRealCoords(reduced, self%insPoint(1:3, iPoint))
      self%insProb(iPoint) = 1E0_dp
    enddo
    call MolData(molType) % molConstruct % ReverseConfig(self%oldpart(1:1), &
                                                         box1, &
                                                         ProbSubOut, &
                                                         accept, &
                                                         self%insPoint(1:3, 1:nInsPoints), &
                                                         self%insProb(1:nInsPoints)) 


    !Forward Probability = 1/N_box1 * 1/vol_Box2 * ForwardGenProbility
    !Reverse Probability = 1/(N_box2+1) * 1/vol_Box1 * ReverseGenProbility
    !Rev/For = N_box1 * vol_box2 / ( (N_box2+1) * vol_Box1 ) * GenProb
    vol = box2 % GetThermo_String(volume)
    Prob = real(box1%nMolTotal, dp) * vol * ProbSubIn

    vol = box1 % GetThermo_String(volume)
    Prob = Prob/(real(box2%nMolTotal+1, dp) * vol * ProbSubOut)

    extraTerms = sampling % GetExtraTerms(self%oldpart(1:1), box1)
    half = sampling % GetExtraTerms(self%newpart(1:nAtoms), box2)
    extraTerms = extraTerms + half


    !Accept/Reject
    accept = sampling % MakeDecision2Box(box1,  box2, E_Diff1, E_Diff2, &
                            self%oldPart(1:1), self%newPart(1:nAtoms), inprob=prob, &
                            extraIn=extraTerms )
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call Box1 % UpdateEnergy(E_Diff1, E_Inter1, E_Intra1)
      call Box1 % DeleteMol(self%oldPart(1)%molIndx)

      call Box2 % UpdateEnergy(E_Diff2, E_Inter2, E_Intra2)
      call Box2 % UpdatePosition(self%newpart(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine ParticleExchange_Maintenance(self)
    implicit none
    class(ParticleExchange), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine ParticleExchange_Prologue(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    use ClassyConstants, only: pi
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(ParticleExchange), intent(inout) :: self
    integer :: iType, maxAtoms, maxPoints
    integer :: nBoxes

    maxAtoms = 0
    maxPoints = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
      if(MolData(iType)%molConstruct%GetNInsertPoints() > maxPoints) then
        maxPoints = MolData(iType)%molConstruct%GetNInsertPoints()
      endif
    enddo

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif
!    write(*,*) self%ubVol

    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(2000,maxAtoms ) )
    allocate( self%newPart(1:maxAtoms) )
  end subroutine
!=========================================================================
  subroutine ParticleExchange_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(ParticleExchange), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Basic Exchange Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Basic Exchange Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Basic Exchange Acceptance Rate: ", accptRate


  end subroutine
!========================================================
end module
!========================================================
