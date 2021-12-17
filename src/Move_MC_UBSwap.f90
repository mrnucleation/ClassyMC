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

    real(dp), allocatable :: insPoint(:,:)
    real(dp), allocatable :: insProb(:)
    type(Addition), allocatable :: newPart(:)
    type(Deletion) :: oldPart(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => UB_Swap_Constructor
!      procedure, pass :: GeneratePosition => UB_Swap_GeneratePosition
      procedure, pass :: AllocateProb => UB_Swap_AllocateProb
      procedure, pass :: FullMove => UB_Swap_FullMove
      procedure, pass :: SwapIn => UB_Swap_SwapIn
      procedure, pass :: SwapOut => UB_Swap_SwapOut
      procedure, pass :: CheckReversibility => UB_Swap_CheckReversibility
      procedure, pass :: CountSites => UB_Swap_CountSites
      procedure, pass :: CreateForward => UB_Swap_CreateForward
!      procedure, pass :: Maintenance => UB_Swap_Maintenance
      procedure, pass :: ProcessIO => UB_Swap_ProcessIO
      procedure, pass :: Prologue => UB_Swap_Prologue
      procedure, pass :: Epilogue => UB_Swap_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine UB_Swap_Constructor(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes, mostAtoms
    implicit none
    class(UB_Swap), intent(inout) :: self
    integer :: nBoxes



    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif


  end subroutine
!========================================================
  subroutine UB_Swap_AllocateProb(self, nInsPoints)
    implicit none
    class(UB_Swap), intent(inout) :: self
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
!========================================================
  subroutine UB_Swap_CreateForward(self, targetatoms, inspoint)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(UB_Swap), intent(inout) :: self
    real(dp), intent(in) :: targetatoms(:,:) 
    real(dp), intent(out) :: inspoint(1:3)

    real(dp) :: dx, dy, dz, radius


    call Generate_UnitSphere(dx, dy, dz)
    radius = self % ubRad * grnd()**(1.0E0_dp/3.0E0_dp)
    dx = radius * dx
    dy = radius * dy
    dz = radius * dz
    insPoint(1) = targetatoms(1,1) + dx
    insPoint(2) = targetatoms(2,1) + dy
    insPoint(3) = targetatoms(3,1) + dz



  end subroutine
!===============================================
  subroutine UB_Swap_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    implicit none
    class(UB_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

#ifdef DETAILED
    call self%CheckReversibility(trialBox, accept)
    return
#endif

    if(grnd() > 0.5E0_dp) then
      call self % SwapIn(trialBox, accept)
    else
      call self % SwapOut(trialBox, accept)
    endif

  end subroutine
!===============================================
  subroutine UB_Swap_CheckReversibility(self, trialBox, accept) 
!   Function which performs a swap move and then undoes it
!   to see if the forward and reverse probabilties are consistent.
    use ParallelVar, only: nout
    implicit none
    class(UB_Swap), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove
    real(dp) :: A12, A21
    real(dp) :: E1, E2
    real(dp) :: P12, P21
    real(dp) :: detailbalance

    E1 = trialBox%ETotal
    call self % SwapIn(trialBox, accept, moveid=nMove, ProbIn=A12)
    if(.not. accept) then
      return
    endif

    E2 = trialBox%ETotal
    call self % SwapOut(trialBox, accept, forceid=nMove, ProbOut=A21)
    if(.not. accept) then
      return
    endif

    !Doing this in log space to avoid underflow and overflow issues.
    P12 = log(A12) - trialBox%beta*(E2-E1)
    if(P12 > 0E0_dp) P12 = 0E0_dp
    P21 = log(A21) - trialBox%beta*(E1-E2)
    if(P21 > 0E0_dp) P21 = 0E0_dp


    !Our A12 and A21 in this case are actually the ratio of the two

!    detailbalance = P12/P21
!    detailbalance = A21*detailbalance*exp(-trialBox%beta*(E1-E2))

    detailbalance = log(A21) + P12 - P21 - trialBox%beta*(E1-E2)
    if(abs(detailbalance) > 1E-6_dp) then
      write(nout,*) "E1, E2:", E1, E2
      write(nout,*) "A12, A21, A12*A21:", A12, A21, A12*A21
      write(nout,*) "ln(P12), ln(P21):", P12, P21
      write(nout,*) "Detailed:", detailbalance
      write(nout,*) "Violation of Detailed Balance Detected!"
      error stop
    endif



  end subroutine
!===============================================
  subroutine UB_Swap_SwapIn(self, trialBox, accept, moveid, ProbIn) 
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
    integer, intent(out), optional :: moveid
    real(dp), intent(out), optional :: ProbIn

    logical :: relative
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove, nCount
    integer :: iIns, nInsPoints
    integer :: iAtom, iDisp
    integer :: molType, molStart, molEnd, atomIndx, nAtoms
    integer :: targStart,targEnd, targType
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, radius, extraTerms
    real(dp) :: E_Inter, E_Intra
    real(dp) :: Prob = 1E0_dp
    real(dp) :: ProbSub, GenProb, GasProb
    integer :: slice(1:2)
    real(dp), pointer :: targetatoms(:,:) => null()


    self % atmps = self % atmps + 1E0_dp
    self % inatmps = self % inatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialbox, self%newPart)


    nType = floor(nMolTypes * grnd() + 1E0_dp)
    if(trialBox%NMol(nType) + 1 > trialBox%NMolMax(nType)) then
      accept = .false.
      return
    endif
    nMove = trialBox%NMol(nType) + 1
    nMove = trialbox%MolGlobalIndx(nType, nMove)
    call trialBox % GetMolData(nMove, molStart=molStart, &
                               molEnd=molEnd)

    !Choose an atom to serve as the target for the new molecule.
!    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
!    call FindMolecule(trialbox, rawIndx, nTarget)
    nTarget = self%UniformMoleculeSelect(trialBox)
    call trialBox % GetMolData(nTarget, molType=molType, molStart=targStart, molEnd=targEnd)

    slice(1) = targStart
    slice(2) = targEnd
    call trialbox%GetCoordinates(targetatoms, slice=slice)
    !Choose the position relative to the target atom 
    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    call self%AllocateProb(nInsPoints)
    do iIns = 1, nInsPoints
      call self%CreateForward(targetatoms, self%inspoint(1:3, iIns))
      self%insprob(iIns) = 1E0_dp 
    enddo

    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      self%newPart(iAtom)%molType = nType
      self%newPart(iAtom)%molIndx = nMove
      self%newPart(iAtom)%atmIndx = atomIndx
    enddo

    call MolData(molType) % molConstruct % GenerateConfig(trialBox, &
                                                          self%newPart(1:nAtoms),&
                                                          ProbSub, &
                                                          accept, &
                                                          self%insPoint(:, 1:nInsPoints), &
                                                          self%insProb(1:nInsPoints))
    if(.not. accept) then
      return
    endif
    GenProb = ProbSub

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
    if(nCount == 0) then
      return
    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%newPart(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
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





    !Compute the generation probability
    Prob = real(trialBox%nMolTotal, dp) * self%ubVol
    Prob = Prob/(real(nCount, dp) * real(trialBox%nMolTotal+1, dp))

!    write(*,*) "Prob In", trialBox%nMolTotal, self%ubVol, nCount
!    write(*,*) "Prob In", Prob, GenProb, E_Inter, E_Intra, E_Diff
!    write(*,*) "Ratio:", exp(-trialBox%beta*E_Intra)/GenProb
!    write(*,*) 

    call MolData(molType) % molConstruct % GasConfig(GasProb)
    Prob = GasProb*Prob/GenProb

    !Get the chemical potential term for GCMC
    extraTerms = sampling % GetExtraTerms(self%newpart(1:nAtoms), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%newPart(1:nAtoms), inProb=Prob, extraIn=extraTerms)
    if( present(moveid) ) then
      accept = .true.
    endif

    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      if( present(moveid) )  moveid = self%newPart(1)%molIndx
      if( present(ProbIn) )  ProbIn = Prob 
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
!      write(*,*) E_Diff, E_Inter, E_Intra
      call trialBox % UpdatePosition(self%newPart(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine UB_Swap_SwapOut(self, trialBox, accept, forceid, ProbOut) 
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
    integer, intent(in), optional :: forceid
    real(dp), intent(out), optional :: ProbOut
    integer :: molType, molStart, molEnd
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, nNei, nCount
    integer :: nInsPoints, iIns
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, extraTerms
    real(dp) :: E_Inter, E_Intra
    real(dp) :: GenProb, ProbSub, GasProb
    real(dp) :: Prob = 1E0_dp
    real(dp) :: Probconstruct = 1E0_dp

    integer :: slice(1:2)
    real(dp), pointer :: atoms(:,:) => null()




    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialbox, self%oldPart)

    !Propose move
    if(.not. present(forceid) ) then
      rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
      call FindMolecule(trialbox, rawIndx, nMove)
    else
      nMove = forceid
    endif
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart, molEnd=molEnd)
    slice(1) = molStart
    slice(2) = molEnd
    call trialbox%GetCoordinates(atoms, slice=slice)

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
      return
    endif

    !Energy Calculation
    call trialBox%ComputeEnergyDelta(self%oldpart(1:1),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)

    if(.not. accept) then
!      write(*,*) "Energy Rejection"
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%oldPart(1:1), E_Diff )
    if(.not. accept) then
      return
    endif


    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    call self%AllocateProb(nInsPoints)
    do iIns = 1, nInsPoints-1
      call self%CreateForward(atoms, self%inspoint(1:3, iIns))
      self%insprob(iIns) = 1E0_dp 
    enddo
    self%insProb(nInsPoints) = 1E0_dp

    call MolData(molType) % molConstruct % ReverseConfig(self%oldpart(1:1), &
                                                         trialBox, &
                                                         ProbSub, &
                                                         accept, &
                                                         self%insPoint(1:3, 1:nInsPoints), &
                                                         self%insProb(1:nInsPoints)) 
    GenProb = ProbSub

    call self % CountSites(trialBox, &
                           atoms(1,1), &
                           atoms(2,1), &
                           atoms(3,1), &
                           trialBox%NeighList(1)%nNeigh(molStart), &
                           trialBox%NeighList(1)%list(:,molStart), &
                           nNei  )

    call MolData(molType) % molConstruct % GasConfig(GasProb)

    Prob = real(nNei, dp) * real(trialBox%nMolTotal, dp)
    Prob = Prob/(real(trialBox%nMolTotal-1, dp) * self%ubVol)

!    write(*,*) "Prob Out", trialBox%nMolTotal, self%ubVol, nNei
!    write(*,*) "Prob Out", Prob, GenProb, E_Inter, E_Intra, E_Diff
!    write(*,*) "Ratio:", exp(-trialBox%beta*E_Intra)*GenProb
!    write(*,*) 

!    write(*,*) GasProb
    Prob = Prob*GenProb/GasProb

    !Get chemical potential term
    extraTerms = sampling % GetExtraTerms(self%oldpart(1:1), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%oldPart(1:1), inProb=Prob, extraIn=extraTerms)
    if( present(forceid) ) then
      accept = .true.
    endif

    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
      if( present(ProbOut) )  ProbOut = Prob
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
!      write(*,*) E_Diff, E_Inter, E_Intra
!      write(*,*) 
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
        call trialbox%Boundary(rx,ry,rz)
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

    allocate( self%newPart(1:maxAtoms) )
    call self%CreateTempArray(maxAtoms)
  end subroutine
!=========================================================================
  subroutine UB_Swap_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(UB_Swap), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,*) 
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
!=========================================================================
  subroutine UB_Swap_ProcessIO(self, line, lineStat)
    use ClassyConstants, only: pi
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(UB_Swap), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("radius")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%ubRad = realVal
        self%ubRadSq = realVal*realVal
        self%ubVol = (4E0_dp/3E0_dp)*pi*self%ubRad**3

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine

!========================================================
end module
!========================================================
