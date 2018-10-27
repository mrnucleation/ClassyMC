!========================================================
module MCMove_UB_Swap
use CoordinateTypes, only: Addition, Deletion
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: UB_Swap
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: inatmps = 1E-30_dp
    real(dp) :: inaccpt = 0E0_dp
    real(dp) :: outatmps = 1E-30_dp
    real(dp) :: outaccpt = 0E0_dp
    real(dp) :: ubRad = 4.0E0_dp
    real(dp) :: ubRadSq = 4.0E0_dp**2
    real(dp) :: ubVol = 0E0_dp
    type(Addition), allocatable :: newPart(:)
    type(Deletion) :: oldPart(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => UB_Swap_Constructor
!      procedure, pass :: GeneratePosition => UB_Swap_GeneratePosition
      procedure, pass :: FullMove => UB_Swap_FullMove
      procedure, pass :: SwapIn => UB_Swap_SwapIn
      procedure, pass :: SwapOut => UB_Swap_SwapOut
      procedure, pass :: CountSites => UB_Swap_CountSites
!      procedure, pass :: Maintenance => UB_Swap_Maintenance
      procedure, pass :: Prologue => UB_Swap_Prologue
      procedure, pass :: Epilogue => UB_Swap_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine UB_Swap_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(UB_Swap), intent(inout) :: self
!    integer :: iType, maxAtoms



!    allocate( self%tempNNei(1) )
!    allocate( self%tempList(200, 1) )
  end subroutine
!========================================================
!  subroutine UB_Swap_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(UB_Swap), intent(in) :: self
!    type(DisplacementNew), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine UB_Swap_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    implicit none
    class(UB_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    if(grnd() > 0.5E0_dp) then
      call self % SwapIn(trialBox, accept)
    else
      call self % SwapOut(trialBox, accept)
    endif

  end subroutine
!===============================================
  subroutine UB_Swap_SwapIn(self, trialBox, accept) 
    use Box_Utility, only: FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(UB_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    logical :: relative
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove, nCount
    integer :: iAtom, iDisp
    integer :: molType, molStart, molEnd, atomIndx, nAtoms
    integer :: targStart
    real(dp) :: insPoint(1:3)
    real(dp) :: dx, dy, dz
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
    !Choose an atom to serve as the target for the new molecule.
!    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp)
!    call FindAtom(trialbox, rawIndx, nTarget)
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nTarget)
    call trialBox % GetMolData(nTarget, molStart=targStart)


!    call MolData(molType) % molConstruct % ReverseConfig( trialBox, ProbSub, accept)

    !Choose the position relative to the target atom 
    call Generate_UnitSphere(dx, dy, dz)
    radius = self % ubRad * grnd()**(1.0E0_dp/3.0E0_dp)
    dx = radius * dx
    dy = radius * dy
    dz = radius * dz
    insPoint(1) = trialBox%atoms(1, targStart) + dx
    insPoint(2) = trialBox%atoms(2, targStart) + dy
    insPoint(3) = trialBox%atoms(3, targStart) + dz


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
      if(iAtom == 1) then
        call self % CountSites(trialBox, self%newPart(1)%x_new, self%newPart(1)%y_new, &
                               self%newPart(1)%z_new, self%tempNNei(1), self%tempList(:,1), &
                               nCount )
      endif
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
    Prob = real(trialBox%nMolTotal, dp) * self%ubVol
    Prob = Prob/(real(nCount, dp) * real(trialBox%nMolTotal+1, dp))
!    write(*,*) "Prob In", Prob, E_Diff, trialBox%nMolTotal, self%ubVol, nCount, trialBox%nMolTotal+1

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%newPart(1:nAtoms))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%newPart(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine UB_Swap_SwapOut(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData
    use Box_Utility, only: FindMolecule
    use ParallelVar, only: nout
    implicit none
    class(UB_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: molType, molStart, molEnd
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, nNei, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE
    real(dp) :: Prob = 1E0_dp
    real(dp) :: Probconstruct = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart)


    if(trialBox%NMol(molType) - 1 < trialBox%NMolMin(molType)) then
!      write(*,*) "Bounds Rejection"
      accept = .false.
      return
    endif


    self%oldPart(1)%molType = molType
    self%oldPart(1)%molIndx = nMove

    !Check Constraint
    accept = trialBox % CheckConstraint( self%oldPart(1:1) )
    if(.not. accept) then
!      write(*,*) "Constraint Rejection"
!      write(*,*) "============================================"
      return
    endif

    !Energy Calculation
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
!      write(*,*) "Energy Rejection"
      return
    endif

    call MolData(molType) % molConstruct % ReverseConfig( trialBox, probconstruct, accept)
    call self % CountSites(trialBox, &
                           trialBox%atoms(1,molStart), &
                           trialBox%atoms(2,molStart), &
                           trialBox%atoms(3,molStart), &
                           trialBox%NeighList(1)%nNeigh(molStart), &
                           trialBox%NeighList(1)%list(:,molStart), &
                           nNei  )

    Prob = real(nNei, dp) * real(trialBox%nMolTotal, dp)
    Prob = Prob/(real(trialBox%nMolTotal-1, dp) * self%ubVol)
!    write(*,*) "Prob Out:", Prob, trialBox%nMolTotal, self%ubVol, nNei, trialBox%nMolTotal-1

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%oldPart(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(self%oldPart(1)%molIndx)
    endif


  end subroutine
!========================================================
  subroutine UB_Swap_CountSites(self, trialBox, x,y,z, nNei, tempList, nCount)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(UB_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    integer, intent(in) :: NNei
    integer, intent(in) :: tempList(:)
    real(dp), intent(in) :: x,y,z
    integer, intent(out) :: NCount
    integer :: iNei, iAtom, molIndx
    real(dp) :: rx, ry,rz, rsq


    nCount = 0
    do iNei = 1, NNei
      iAtom = tempList(iNei)
      molIndx = trialBox%MolIndx(iAtom)
      if(iAtom == trialBox%MolStartIndx(molIndx)) then
        rx = x - trialBox % atoms(1, iAtom)
        ry = y - trialBox % atoms(2, iAtom)
        rz = z - trialBox % atoms(3, iAtom)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%ubRadSq) then
          nCount = nCount + 1
        endif
      endif
    enddo

  end subroutine
!=========================================================================
  subroutine UB_Swap_Maintenance(self)
    implicit none
    class(UB_Swap), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine UB_Swap_Prologue(self)
    use ParallelVar, only: nout
    use ClassyConstants, only: pi
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(UB_Swap), intent(inout) :: self
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


    self%ubVol = (4E0_dp/3E0_dp)*pi*self%ubRad**3
    self%ubRadSq = self%ubRad * self%ubRad
!    write(*,*) self%ubVol

    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(200,maxAtoms ) )
    allocate( self%newPart(1:maxAtoms) )
  end subroutine
!=========================================================================
  subroutine UB_Swap_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(UB_Swap), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "UB Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "UB Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "UB Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "UB Out Moves Accepted: ", nint(self%outaccpt)
    write(nout,"(1x,A,I15)") "UB Out Moves Attempted: ", nint(self%outatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "UB Out Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "UB In Moves Accepted: ", nint(self%inaccpt)
    write(nout,"(1x,A,I15)") "UB In Moves Attempted: ", nint(self%inatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "UB In Acceptance Rate: ", accptRate
 

  end subroutine
!========================================================
end module
!========================================================
