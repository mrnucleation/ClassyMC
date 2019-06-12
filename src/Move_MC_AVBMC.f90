!========================================================
module MCMove_AVBMC
use CoordinateTypes, only: Addition, Deletion
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: AVBMC
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    logical :: ebias = .false.
    integer :: neilistindx = 1
    real(dp) :: inatmps = 1E-30_dp
    real(dp) :: inaccpt = 0E0_dp
    real(dp) :: outatmps = 1E-30_dp
    real(dp) :: outaccpt = 0E0_dp
    real(dp) :: avbmcRad = 4.0E0_dp
    real(dp) :: avbmcRadSq = 4.0E0_dp**2
    real(dp) :: avbmcVol = 0E0_dp
    type(Addition), allocatable :: newPart(:)
    type(Deletion) :: oldPart(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => AVBMC_Constructor
!      procedure, pass :: GeneratePosition => AVBMC_GeneratePosition
      procedure, pass :: FullMove => AVBMC_FullMove
      procedure, pass :: SwapIn => AVBMC_SwapIn
      procedure, pass :: SwapOut => AVBMC_SwapOut
      procedure, pass :: CountSites => AVBMC_CountSites
      procedure, pass :: SelectNeigh => AVBMC_SelectNeigh
      procedure, pass :: SelectNeighReverse => AVBMC_SelectNeighReverse
!      procedure, pass :: Maintenance => AVBMC_Maintenance
      procedure, pass :: ProcessIO => AVBMC_ProcessIO
      procedure, pass :: Prologue => AVBMC_Prologue
      procedure, pass :: Epilogue => AVBMC_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine AVBMC_Constructor(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(AVBMC), intent(inout) :: self
    integer :: nBoxes

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif




!    allocate( self%tempNNei(1) )
!    allocate( self%tempList(200, 1) )
  end subroutine
!===============================================
  subroutine AVBMC_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    if(grnd() > 0.5E0_dp) then
      call self % SwapIn(trialBox, accept)
    else
      call self % SwapOut(trialBox, accept)
    endif

  end subroutine
!===============================================
  subroutine AVBMC_SwapIn(self, trialBox, accept) 
    use Box_Utility, only: FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(AVBMC), intent(inout) :: self
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
    real(dp) :: E_Diff, biasE, radius, extraTerms
    real(dp) :: Prob = 1E0_dp
    real(dp) :: ProbSub, ProbSel

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
    radius = self % avbmcRad * grnd()**(1.0E0_dp/3.0E0_dp)
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

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%newPart(1:nAtoms), E_Diff )
    if(.not. accept) then
      return
    endif


    !Compute the probability of reversing the move
    call self%SelectNeighReverse(trialBox, nTarget, ProbSel)


    !Compute the generation probability
    Prob = real(trialBox%nMolTotal, dp) * self%avbmcVol
    Prob = Prob*ProbSel/real(trialBox%nMolTotal+1, dp)
!    write(*,*) "Prob In", Prob, E_Diff, trialBox%nMolTotal, self%avbmcVol, nCount, trialBox%nMolTotal+1, ProbSel

    !Get the chemical potential term for GCMC
    extraTerms = sampling % GetExtraTerms(self%newpart(1:nAtoms), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%newPart(1:nAtoms), inProb=Prob, extraIn=extraTerms)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%newPart(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine AVBMC_SwapOut(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData
    use Box_Utility, only: FindMolecule
    use ParallelVar, only: nout
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: molType, molStart, molEnd
    integer :: nTarget, nMove, rawIndx, iConstrain
    integer :: CalcIndex, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, extraTerms
    real(dp) :: ProbSel = 1E0_dp
    real(dp) :: Prob = 1E0_dp
    real(dp) :: Probconstruct = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nTarget)
    call trialBox % GetMolData(nTarget, molType=molType, molStart=molStart)
    call self%SelectNeigh(trialBox, nTarget, nMove, ProbSel)

    if(nMove < 1) then
      accept = .false.
      return
    endif
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

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%oldPart(1:1), E_Diff )
    if(.not. accept) then
      return
    endif



    call MolData(molType) % molConstruct % ReverseConfig( trialBox, probconstruct, accept)
    Prob = real(trialBox%nMolTotal, dp)/ProbSel
    Prob = Prob/(real(trialBox%nMolTotal-1, dp) * self%avbmcVol)
!    write(*,*) "Prob Out:", Prob, trialBox%nMolTotal, self%avbmcVol, trialBox%nMolTotal-1, ProbSel

    !Get chemical potential term
    extraTerms = sampling % GetExtraTerms(self%oldpart(1:1), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%oldPart(1:1), inProb=Prob, extraIn=extraTerms)

    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(self%oldPart(1)%molIndx)
    endif


  end subroutine
!========================================================
  subroutine AVBMC_CountSites(self, trialBox, targetIndx, nNei)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    integer, intent(in) :: targetIndx
    integer, intent(out) :: NNei
    integer :: iNei, iAtom, molIndx
    integer :: targetStart
    real(dp) :: rx, ry,rz, rsq

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call trialbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call trialbox%GetCoordinates(atoms)

    nNei = 0
    do iNei = 1, nNeigh(targetIndx)
      iAtom = neighlist(iNei, targetIndx)
      molIndx = trialBox%MolIndx(iAtom)
      if(iAtom == trialBox%MolStartIndx(molIndx)) then
        rx = atoms(1, iAtom) - atoms(1, targetIndx)
        ry = atoms(2, iAtom) - atoms(2, targetIndx)
        rz = atoms(3, iAtom) - atoms(3, targetIndx)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%avbmcRadSq) then
          nNei = nNei + 1
        endif
      endif
    enddo

  end subroutine
!========================================================
  subroutine AVBMC_SelectNeigh(self, trialBox, targetIndx, removeMolIndx, ProbOut)
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: ListRNG, grnd
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    integer, intent(in) :: targetIndx
    integer, intent(out) :: removeMolIndx
    real(dp), intent(out) :: ProbOut
    integer :: removeIndx
    integer :: iNei, iAtom, molIndx
    integer :: nNei
    integer :: targetStart
    integer :: removelist(1:60)
    real(dp) :: rx, ry,rz, rsq
    real(dp) :: weights(1:60), norm, maxweight, EMol

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call trialbox%GetNeighborList(self%neilistindx, neighlist, nNeigh)
    call trialbox%GetCoordinates(atoms)

    nNei = 0
    removelist = 0
    do iNei = 1, nNeigh(targetIndx)
      iAtom = neighlist(iNei, targetIndx)
      molIndx = trialBox%MolIndx(iAtom)
      if(iAtom == trialBox%MolStartIndx(molIndx)) then
        rx = atoms(1, iAtom) - atoms(1, targetIndx)
        ry = atoms(2, iAtom) - atoms(2, targetIndx)
        rz = atoms(3, iAtom) - atoms(3, targetIndx)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%avbmcRadSq) then
          nNei = nNei + 1
          removelist(nNei) = molIndx
        endif
      endif
    enddo


    if(nNei < 1) then
      removeMolIndx = -1
      ProbOut = 0E0_dp
      return
    endif

    if(self%ebias) then
      weights = 0E0_dp
      do iNei = 1, nNei
        call trialbox%GetMolEnergy(removelist(iNei), EMol, newstate=.false.)
        weights(iNei) = EMol
      enddo
      maxweight = maxval(weights(1:nNei))
      do iNei = 1, nNei
        EMol = weights(iNei)
        weights(iNei) = exp((EMol-maxweight)*trialBox%beta)
      enddo     
      norm = sum(weights(1:nNei))
      removeIndx = ListRNG(weights, norm)
      removeMolIndx = removelist(removeIndx)
      ProbOut = weights(removeIndx)/norm
    else
      removeIndx = floor(grnd()*nNei + 1E0_dp)
      removeMolIndx = removelist(removeIndx)
      ProbOut = 1E0_dp/real(nNei,dp)
    endif

  end subroutine
!========================================================
  subroutine AVBMC_SelectNeighReverse(self, trialBox, targetIndx, ProbOut)
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: ListRNG, grnd
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    integer, intent(in) :: targetIndx
    real(dp), intent(out) :: ProbOut
    integer :: iNei, iAtom, molIndx
    integer :: nNei
    integer :: targetStart
    integer :: removelist(1:60)
    real(dp) :: rx, ry,rz, rsq
    real(dp) :: weights(1:60), norm, maxweight, EMol

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call trialbox%GetNeighborList(self%neilistindx, neighlist, nNeigh)
    call trialbox%GetCoordinates(atoms)

    nNei = 0
    do iNei = 1, nNeigh(targetIndx)
      iAtom = neighlist(iNei, targetIndx)
      molIndx = trialBox%MolIndx(iAtom)
      if(iAtom == trialBox%MolStartIndx(molIndx)) then
        rx = atoms(1, iAtom) - atoms(1, targetIndx)
        ry = atoms(2, iAtom) - atoms(2, targetIndx)
        rz = atoms(3, iAtom) - atoms(3, targetIndx)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%avbmcRadSq) then
          nNei = nNei + 1
          removelist(nNei) = molIndx
        endif
      endif
    enddo
    nNei = nNei + 1
    removelist(nNei) = self%newPart(1)%molIndx

    if(nNei < 1) then
      ProbOut = 0E0_dp
      return
    endif

    if(self%ebias) then
      weights = 0E0_dp
      do iNei = 1, nNei
        call trialbox%GetMolEnergy(removelist(iNei), EMol, newstate=.true.)
        weights(iNei) = EMol
      enddo
      maxweight = maxval(weights(1:nNei))
      do iNei = 1, nNei
        EMol = weights(iNei)
        weights(iNei) = exp((EMol-maxweight)*trialBox%beta)
      enddo     
      norm = sum(weights(1:nNei))
      ProbOut = weights(nNei)/norm
    else
      ProbOut = 1E0_dp/real(nNei,dp)
    endif

  end subroutine
!=========================================================================
  subroutine AVBMC_Maintenance(self)
    implicit none
    class(AVBMC), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine AVBMC_Prologue(self)
    use ParallelVar, only: nout
    use ClassyConstants, only: pi
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(AVBMC), intent(inout) :: self
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


    self%avbmcVol = (4E0_dp/3E0_dp)*pi*self%avbmcRad**3
    self%avbmcRadSq = self%avbmcRad * self%avbmcRad
!    write(*,*) self%avbmcVol

    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(2000,maxAtoms ) )
    allocate( self%newPart(1:maxAtoms) )
  end subroutine
!=========================================================================
  subroutine AVBMC_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(AVBMC), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "AVBMC Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "AVBMC Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "AVBMC Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "AVBMC Out Moves Accepted: ", nint(self%outaccpt)
    write(nout,"(1x,A,I15)") "AVBMC Out Moves Attempted: ", nint(self%outatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "AVBMC Out Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "AVBMC In Moves Accepted: ", nint(self%inaccpt)
    write(nout,"(1x,A,I15)") "AVBMC In Moves Attempted: ", nint(self%inatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "AVBMC In Acceptance Rate: ", accptRate
 

  end subroutine
!=========================================================================
  subroutine AVBMC_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(AVBMC), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("energybias")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%ebias = logicVal

      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%neilistindx = intVal

      case("radius")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%avbmcRad = realVal
        self%avbmcRadSq = realVal*realVal

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine

!========================================================
end module
!========================================================
