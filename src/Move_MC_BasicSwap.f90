!========================================================
module MCMove_Basic_Swap
use CoordinateTypes, only: Addition, Deletion
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: Basic_Swap
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: inatmps = 1E-30_dp
    real(dp) :: inaccpt = 0E0_dp
    real(dp) :: outatmps = 1E-30_dp
    real(dp) :: outaccpt = 0E0_dp
    type(Addition), allocatable :: newPart(:)
    type(Deletion) :: oldPart(1:1)

    real(dp), allocatable :: insPoint(:,:)
    real(dp), allocatable :: insProb(:)
!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
    procedure, pass :: SafetyCheck => Basic_Swap_SafetyCheck
      procedure, pass :: Constructor => Basic_Swap_Constructor
!      procedure, pass :: GeneratePosition => Basic_Swap_GeneratePosition
      procedure, pass :: FullMove => Basic_Swap_FullMove
      procedure, pass :: AllocateProb => Basic_Swap_AllocateProb
      procedure, pass :: SwapIn => Basic_Swap_SwapIn
      procedure, pass :: SwapOut => Basic_Swap_SwapOut
!      procedure, pass :: Maintenance => Basic_Swap_Maintenance
      procedure, pass :: Prologue => Basic_Swap_Prologue
      procedure, pass :: Epilogue => Basic_Swap_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine Basic_Swap_SafetyCheck(self)
    use BoxData, only: boxArray
    use SimpleSimBox, only: SimpleBox
    implicit none
    class(Basic_Swap), intent(inout) :: self
    integer :: iBox

    !Ensure a boundary condition is set.
    do iBox = 1, size(BoxArray)
      select type(box => BoxArray(iBox)%box)
        type is(SimpleBox)
            write(*,*) "WARNING! Basic Swap move is not designed to work without"
            write(*,*) "a bounding box!"
            write(*,*) "Box Number:", iBox
      end select
    enddo



  end subroutine
!========================================================
  subroutine Basic_Swap_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(Basic_Swap), intent(inout) :: self
!    integer :: iType, maxAtoms



  end subroutine
!========================================================
!  subroutine Basic_Swap_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(Basic_Swap), intent(in) :: self
!    type(DisplacementNew), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!========================================================
  subroutine Basic_Swap_AllocateProb(self, nInsPoints)
    implicit none
    class(Basic_Swap), intent(inout) :: self
    integer, intent(in) :: nInsPoints
    integer :: errstat=0

    if((.not. allocated(self%insPoint)) .and. (.not. allocated(self%insProb)) ) then
      allocate(self%insPoint(1:3, 1:nInsPoints), stat=errstat)
      if(errstat /= 0) then
        error stop "Allocation Error trying to allocate newPart!"
      endif
      allocate(self%insProb(1:nInsPoints), stat=errstat)
      if(errstat /= 0) then
        error stop "Allocation Error trying to allocate insProb!"
      endif
    else if(size(self%insPoint, 2) < nInsPoints) then
      deallocate(self%insPoint)
      deallocate(self%insProb)
      allocate(self%insPoint(1:3, 1:nInsPoints), stat=errstat)
      if(errstat /= 0) then
        error stop "Allocation Error trying to allocate newPart!"
      endif
      allocate(self%insProb(1:nInsPoints), stat=errstat)
      if(errstat /= 0) then
        error stop "Allocation Error trying to allocate insProb!"
      endif
    else
      return
    endif
    self%insPoint = 0E0_dp
    self%insProb = 1E0_dp

  end subroutine

!===============================================
  subroutine Basic_Swap_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    implicit none
    class(Basic_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    if(grnd() > 0.5E0_dp) then
      call self % SwapIn(trialBox, accept)
    else
      call self % SwapOut(trialBox, accept)
    endif

  end subroutine
!===============================================
  subroutine Basic_Swap_SwapIn(self, trialBox, accept) 
    use Box_Utility, only: FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(Basic_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    logical :: relative
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove, nCount
    integer :: i, iAtom, iPoint, iDisp, nInsPoints
    integer :: molType, molStart, molEnd, atomIndx, nAtoms
    real(dp) :: reduced(1:3)
    real(dp) :: dx, dy, dz, vol
    real(dp) :: E_Diff, biasE, radius
    real(dp) :: E_Inter, E_Intra
    real(dp) :: Prob = 1E0_dp
    real(dp) :: GenProb, GasProb 
    character(len=30), parameter :: volume = "volume"

    self % atmps = self % atmps + 1E0_dp
    self % inatmps = self % inatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialBox, self%newPart)

    nType = floor(nMolTypes * grnd() + 1E0_dp)
    if(trialBox%NMol(nType) + 1 > trialBox%NMolMax(nType)) then
      accept = .false.
      return
    endif

    nMove = trialBox%NMol(nType) + 1
    nMove = trialbox%MolGlobalIndx(nType, nMove)
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart, &
                               molEnd=molEnd)

    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    call self%AllocateProb(nInsPoints)
    do iPoint = 1, nInsPoints
      do i=1,3
        reduced(i) = grnd()
      enddo
      call trialBox%GetRealCoords(reduced, self%insPoint(1:3, iPoint))
      self%insprob(iPoint) = 1E0_dp 
    enddo


    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      self%newPart(iAtom)%molType = molType
      self%newPart(iAtom)%molIndx = nMove
      self%newPart(iAtom)%atmIndx = atomIndx
    enddo

    !Generate the full atomic positions for 
    call MolData(molType) % molConstruct % GenerateConfig(trialBox, &
                                                      self%newPart(1:nAtoms),&
                                                      GenProb, &
                                                      accept, &
                                                      self%insPoint, &
                                                      self%insProb)
    if(.not. accept) then
        return
    endif


    call trialBox % NeighList(1) % GetTempListArray(self%tempList, self%tempNNei)
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
    E_Inter = 0E0_dp
    E_Intra = 0E0_dp
    E_Diff = 0E0_dp
    call trialBox%ComputeEnergyDelta(self%newpart(1:nAtoms),&
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
    accept = trialBox % CheckPostEnergy( self%newPart(1:nAtoms), E_Diff )
    if(.not. accept) then
      return
    endif
    call MolData(molType) % molConstruct % GasConfig(GasProb)

    !Compute the generation probability
    vol = trialBox % GetThermo(5)
!    vol = trialBox % GetThermo_String(volume)
    Prob = vol/real(trialBox%nMolTotal+1, dp) 
    Prob = GasProb*Prob/GenProb

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%newPart(1:nAtoms), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%newPart(1:nAtoms), self%tempList(:,:), self%tempNNei(:))
    endif


  end subroutine
!===============================================
  subroutine Basic_Swap_SwapOut(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData
    use Box_Utility, only: FindMolecule
    use ParallelVar, only: nout
    implicit none
    class(Basic_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i, iPoint, nInsPoints
    integer :: molType, molStart, molEnd
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, nNei, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, vol
    real(dp) :: E_Inter, E_Intra
    real(dp) :: Prob = 1E0_dp
    real(dp) :: GenProb, GasProb
    real(dp) :: reduced(1:3)
    character(len=30), parameter :: volume = "volume"

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialBox, self%oldPart)

    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    if(rawIndx > trialBox%nMolTotal) then
      rawIndx = trialBox%nMolTotal
    endif
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart)

    if(trialBox%NMol(molType) - 1 < trialBox%NMolMin(molType)) then
      accept = .false.
      return
    endif


    self%oldPart(1)%molType = molType
    self%oldPart(1)%molIndx = nMove

    !Check Constraint
    accept = trialBox % CheckConstraint( self%oldPart(1:1) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    call trialBox%ComputeEnergyDelta(self%oldpart(1:1),&
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


    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    call self%AllocateProb(nInsPoints)
    do iPoint = 1, nInsPoints
      do i=1,3
        reduced(i) = grnd()
      enddo
      call trialBox%GetRealCoords(reduced, self%insPoint(1:3, iPoint))
      self%insprob(iPoint) = 1E0_dp 
    enddo

    call MolData(molType) % molConstruct % ReverseConfig(self%oldpart(1:1), &
                                                         trialBox, &
                                                         GenProb, &
                                                         accept, &
                                                         self%insPoint(1:3, 1:nInsPoints), &
                                                         self%insProb(1:nInsPoints)) 

    call MolData(molType) % molConstruct % GasConfig(GasProb)
    vol = trialBox % GetThermo_String(volume)
    Prob = real(trialBox%nMolTotal, dp)/vol
    Prob = Prob*GenProb/GasProb

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,self%oldPart(1:1),inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % DeleteMol(self%oldPart(1)%molIndx)
    endif


  end subroutine
!=========================================================================
  subroutine Basic_Swap_Maintenance(self)
    implicit none
    class(Basic_Swap), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine Basic_Swap_Prologue(self)
    use ParallelVar, only: nout
    use ClassyConstants, only: pi
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(Basic_Swap), intent(inout) :: self
    integer :: iType, maxAtoms, maxPoints, errstat
    integer :: iBox, maxnei

    maxAtoms = 0
    maxPoints = 0
    maxnei = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
      if(MolData(iType)%molConstruct%GetNInsertPoints() > maxPoints) then
        maxPoints = MolData(iType)%molConstruct%GetNInsertPoints()
      endif
    enddo





    allocate( self%newPart(1:maxAtoms) )
    if(errstat /= 0) then
      error stop "Allocation Error trying to allocate newPart!"
    endif
  end subroutine
!=========================================================================
  subroutine Basic_Swap_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(Basic_Swap), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Basic Swap Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Basic Swap Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Basic Swap Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "Basic Swap Out Moves Accepted: ", nint(self%outaccpt)
    write(nout,"(1x,A,I15)") "Basic Swap Out Moves Attempted: ", nint(self%outatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "Basic Swap Out Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "Basic Swap In Moves Accepted: ", nint(self%inaccpt)
    write(nout,"(1x,A,I15)") "Basic Swap In Moves Attempted: ", nint(self%inatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "Basic Swap In Acceptance Rate: ", accptRate
 

  end subroutine
!========================================================
end module
!========================================================
