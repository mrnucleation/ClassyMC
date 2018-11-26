!========================================================
module MCMove_ParticleExchange
use CoordinateTypes, only: Addition, Deletion
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: ParticleExchange
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    type(Addition), allocatable :: newPart(:)
    type(Deletion) :: oldPart(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
    procedure, pass :: SafetyCheck => ParticleExchange_SafetyCheck
      procedure, pass :: Constructor => ParticleExchange_Constructor
!      procedure, pass :: GeneratePosition => ParticleExchange_GeneratePosition
!      procedure, pass :: FullMove => ParticleExchange_FullMove
      procedure, pass :: SwapIn => ParticleExchange_SwapIn
      procedure, pass :: SwapOut => ParticleExchange_SwapOut
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
!=========================================================================
  subroutine ParticleExchange_MultiBox(self, accept)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use RandomGen, only: grnd, ListRNG
    use CommonSampling, only: sampling
    implicit none
    class(ParticleExchange), intent(inout) :: self
    logical, intent(out) :: accept

  end subroutine
!=========================================================================
  subroutine ParticleExchange_MultiBox(self, accept)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use RandomGen, only: grnd, ListRNG
    use CommonSampling, only: sampling

    implicit none
    class(ParticleExchange), intent(inout) :: self
    logical, intent(out) :: accept
    integer :: i, boxNum
    integer :: nMoveOut, nMoveIn
    class(SimpleBox), pointer :: box1, box2
    real(dp) :: dV
    real(dp) :: Prob, Norm, extraTerms, half
    real(dp) :: E_Diff1, E_Diff2, scaleFactor
    real(dp) :: rescale(1:size(self%boxprob))



    !Randomly choose which boxes will exchange particles
    boxNum = ListRNG(self%boxProb)
    box1 => BoxArray(boxNum)%box

    !To avoid picking the same box twice, rescale the probability such that
    !box1's probability is equal to 0.
    rescale = self%boxprob
    rescale(boxNum) = 0E0_dp
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
    call FindMolecule(trialbox, rawIndx, nMoveOut)
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
    do i = 1, 3
      reduced(i) = grnd()
    enddo
    call box2%GetRealCoords(reduced, insPoint)


    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      self%newPart(iAtom)%molType = molType
      self%newPart(iAtom)%molIndx = nMoveIn
      self%newPart(iAtom)%atmIndx = atomIndx
    enddo


    call MolData(molType) % molConstruct % GenerateConfig(box2, self%newPart(1:nAtoms), ProbSub , insPoint)

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
    call box1 % EFunc % Method % DiffECalc(box1, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff1, accept)
    if(.not. accept) then
!      write(*,*) "Energy Rejection"
      return
    endif

    !Box2's Energy Calculation
    call box2% EFunc % Method % DiffECalc(box2, self%newPart(1:nAtoms), self%tempList, &
                                              self%tempNNei, E_Diff2, accept)
    if(.not. accept) then
      return
    endif


    !Compute the Generation Probability
    vol = box1 % GetThermo(3)
    Prob = real(box1%nMolTotal-1, dp)/vol

    vol = box2 % GetThermo(3)
    Prob = Prob*vol/real(box2%nMolTotal+1, dp) 

    extraTerms = sampling % GetExtraTerms(self%oldpart(1:1), box1)
    half = sampling % GetExtraTerms(self%newpart(1:nAtoms), box2)
    extraTerms = extraTerms + half


    !Accept/Reject
    accept = sampling % MakeDecision2Box(box1,  box2, E_Diff1, E_Diff2, &
                            self%oldPart(1:1), self%newPart(1:nAtoms), prob=prob, &
                            extraIn=extraTerms )
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call box1 % UpdateEnergy(E_Diff1)
      call box1 % UpdatePosition(self%disp1(1:1), self%tempList, self%tempNNei)

      call box2 % UpdateEnergy(E_Diff2)
      call box2 % DeleteMol(self%oldPart(1)%molIndx)
    endif

  end subroutine

!===============================================
  subroutine ParticleExchange_SwapIn(self, trialBox, accept) 
    use Box_Utility, only: FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(ParticleExchange), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    logical :: relative
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove, nCount
    integer :: i, iAtom, iDisp
    integer :: molType, molStart, molEnd, atomIndx, nAtoms
    integer :: targStart
    real(dp) :: reduced(1:3)
    real(dp) :: insPoint(1:3)
    real(dp) :: dx, dy, dz, vol
    real(dp) :: E_Diff, biasE, radius
    real(dp) :: Prob = 1E0_dp
    real(dp) :: ProbSub

    self % atmps = self % atmps + 1E0_dp
    self % inatmps = self % inatmps + 1E0_dp
    accept = .true.

!    integer(kind=atomIntType) :: molType, atmIndx, molIndx
!    real(dp) :: x_new, y_new, z_new
!    integer :: listIndex = -1
    nType = floor(nMolTypes * grnd() + 1E0_dp)
    if(trialBox%NMol(nType) + 1 > trialBox%NMolMax(nType)) then
      accept = .false.
      return
    endif

    nMove = trialBox%NMol(nType) + 1
    nMove = trialbox%MolGlobalIndx(nType, nMove)
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart, &
                               molEnd=molEnd)

    do i=1,3
      reduced(i) = grnd()
    enddo
    call trialBox%GetRealCoords(reduced, insPoint)


    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      self%newPart(iAtom)%molType = molType
      self%newPart(iAtom)%molIndx = nMove
      self%newPart(iAtom)%atmIndx = atomIndx
    enddo


    call MolData(molType) % molConstruct % GenerateConfig(trialBox, self%newPart(1:nAtoms), ProbSub , insPoint)

    do iAtom = 1, nAtoms
      call trialBox % NeighList(1) % GetNewList(iAtom, self%tempList, self%tempNNei, &
                                                self%newPart(iAtom))
      self%newPart(iAtom)%listIndex = iAtom
    enddo 

    !Check Constraint
    accept = trialBox % CheckConstraint( self%newPart(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%newPart(1:nAtoms), self%tempList, &
                                              self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif

    !Compute the generation probability
    vol = trialBox % GetThermo(3)
    Prob = vol/real(trialBox%nMolTotal+1, dp) 
!    write(*,*) "Prob In", Prob, E_Diff, trialBox%nMolTotal, self%ubVol, nCount, trialBox%nMolTotal+1

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%newPart(1:nAtoms), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%newPart(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine ParticleExchange_SwapOut(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData
    use Box_Utility, only: FindMolecule
    use ParallelVar, only: nout
    implicit none
    class(ParticleExchange), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: molType, molStart, molEnd
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, nNei, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, vol
    real(dp) :: Prob = 1E0_dp
    real(dp) :: Probconstruct = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.

    !Energy Calculation
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
!      write(*,*) "Energy Rejection"
      return
    endif

    call MolData(molType) % molConstruct % ReverseConfig( trialBox, probconstruct, accept)

    vol = trialBox % GetThermo(3)
    Prob = real(trialBox%nMolTotal, dp)/vol
!    write(*,*) "Prob Out:", Prob, trialBox%nMolTotal, self%ubVol, nNei, trialBox%nMolTotal-1

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,self%oldPart(1:1),inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(self%oldPart(1)%molIndx)
    endif


  end subroutine
!=========================================================================
  subroutine ParticleExchange_Maintenance(self)
    implicit none
    class(ParticleExchange), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine ParticleExchange_Prologue(self)
    use ParallelVar, only: nout
    use ClassyConstants, only: pi
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(ParticleExchange), intent(inout) :: self
    integer :: iType, maxAtoms, maxPoints

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


!    write(*,*) self%ubVol

    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(200,maxAtoms ) )
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

    write(nout,"(1x,A,I15)") "Basic Exchange Out Moves Accepted: ", nint(self%outaccpt)
    write(nout,"(1x,A,I15)") "Basic Exchange Out Moves Attempted: ", nint(self%outatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "Basic Exchange Out Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "Basic Exchange In Moves Accepted: ", nint(self%inaccpt)
    write(nout,"(1x,A,I15)") "Basic Exchange In Moves Attempted: ", nint(self%inatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "Basic Exchange In Acceptance Rate: ", accptRate
 

  end subroutine
!========================================================
end module
!========================================================
