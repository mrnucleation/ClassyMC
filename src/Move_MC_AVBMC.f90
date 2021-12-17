!========================================================
!The AVBMC move is the intra-box variant of the move created by Bin Chen et al. 
!This is a biased swap move that increase the chances of creating and breaking 
!bonded configurations by using a small volume around the molecules present in the 
!system as a "target" location for a newly inserted molecule.
!This move type is suitable for both bulk and cluster simulations as it can be conformed 
!to Distance Criteria. For sparse bulk simulations it is suggested that the user still
!couple this move with addition swap move types that are suitable for low density simulations.
!
! Modifiable Parameters
!   radius (float) => Sets the maxmimum distance in angstroms that the first atom
!                     can be placed at
!                             
!     
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

    real(dp), allocatable :: insPoint(:,:)
    real(dp), allocatable :: insProb(:)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => AVBMC_Constructor
!      procedure, pass :: GeneratePosition => AVBMC_GeneratePosition
      procedure, pass :: FullMove => AVBMC_FullMove
      procedure, pass :: CheckReversibility => AVBMC_CheckReversibility
      procedure, pass :: AllocateProb => AVBMC_AllocateProb
      procedure, pass :: SwapIn => AVBMC_SwapIn
      procedure, pass :: SwapOut => AVBMC_SwapOut
      procedure, pass :: CountSites => AVBMC_CountSites
      procedure, pass :: CreateForward => AVBMC_CreateForward
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
    use Common_MolInfo, only: MolData, nMolTypes, mostAtoms
    implicit none
    class(AVBMC), intent(inout) :: self
    integer :: iBox, nBoxes
    integer :: iMol, nMaxPerMol, nMaxPerBox, maxnei

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    !To allocate the temporary neighborlist, we need to find out how big of a list we need.
    !In order to do this we need the number of atoms in the largest molecule and
    !the largest number of neighbors a single atom can encounter.

    nMaxPerMol = 0
    do iMol = 1, nMolTypes
      if(MolData(iMol)%natoms > nMaxPerMol) then
        nMaxPerMol = MolData(iMol)%natoms
      endif
    enddo
    
    nMaxPerBox = 0
    do iBox = 1, nBoxes
      maxnei = BoxArray(iBox)%box%GetLargestNNei()
      if(nMaxPerBox < maxnei) then
        nMaxPerBox = maxnei
      endif
    enddo


  end subroutine
!===============================================
  subroutine AVBMC_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    implicit none
    class(AVBMC), intent(inout) :: self
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
  subroutine AVBMC_CheckReversibility(self, trialBox, accept) 
!   Function which performs a swap move and then undoes it
!   to see if the forward and reverse probabilties are consistent.
    use ParallelVar, only: nout
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove, targid
    real(dp) :: A12, A21
    real(dp) :: E1, E2
    real(dp) :: P12, P21
    real(dp) :: detailbalance

    E1 = trialBox%ETotal
    call self % SwapIn(trialBox, accept, moveid=nMove, targid=targid, ProbIn=A12)
    if(.not. accept) then
      return
    endif

    E2 = trialBox%ETotal
    call self % SwapOut(trialBox, accept, forcetargid=targid, forceid=nMove, ProbOut=A21)
    if(.not. accept) then
      return
    endif

    if( ((A12 == 0E0_dp) .or. (A21 == 0E0_dp)) .and. (A21 /= A12) ) then
      write(0,*) "Zero Probability observed in move that should be reverisbile"
      write(0,*) "This implies a calculation error in the detailed balance."
      write(0,*) "A12, A21:", A12, A21
      error stop
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
!    write(*,*) "Log Balance", detailbalance
    if(abs(detailbalance) > 1E-6_dp) then
      write(nout,*) "E1, E2:", E1, E2
      write(nout,*) "A12, A21, A12*A21:", A12, A21, A12*A21
      write(nout,*) "ln(P12), ln(P21):", P12, P21
      write(nout,*) "Detailed:", detailbalance
      write(nout,*) "Violation of Detailed Balance Detected!"
      error stop
    endif



  end subroutine

!========================================================
  subroutine AVBMC_AllocateProb(self, nInsPoints)
    implicit none
    class(AVBMC), intent(inout) :: self
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
  subroutine AVBMC_CreateForward(self, targetatoms, inspoint)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(AVBMC), intent(inout) :: self
    real(dp), intent(in) :: targetatoms(:,:) 
    real(dp), intent(out) :: inspoint(1:3)

    real(dp) :: dx, dy, dz, radius


    call Generate_UnitSphere(dx, dy, dz)
    radius = self % avbmcRad * grnd()**(1.0E0_dp/3.0E0_dp)
    dx = radius * dx
    dy = radius * dy
    dz = radius * dz
    insPoint(1) = targetatoms(1,1) + dx
    insPoint(2) = targetatoms(2,1) + dy
    insPoint(3) = targetatoms(3,1) + dz

  end subroutine
!===============================================
  subroutine AVBMC_SwapIn(self, trialBox, accept, moveid, targid, ProbIn) 
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
    integer, intent(out), optional :: moveid, targid
    real(dp), intent(out), optional :: ProbIn
    logical :: relative
    integer :: nInsPoints, iIns
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove, nCount
    integer :: iAtom, iDisp
    integer :: molType, molStart, molEnd, atomIndx, nAtoms
    integer :: targStart, targEnd
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, radius, extraTerms
    real(dp) :: E_Inter, E_Intra
    real(dp) :: Prob = 1E0_dp
    real(dp) :: ProbSub, ProbSel, GenProb
    real(dp) :: GasProb = 1E0_dp

    integer :: slice(1:2)
    real(dp), pointer :: targetatoms(:,:) => null()

    self % atmps = self % atmps + 1E0_dp
    self % inatmps = self % inatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialbox, self%newPart)

!    integer(kind=atomIntType) :: molType, atmIndx, molIndx
!    real(dp) :: x_new, y_new, z_new
!    integer :: listIndex = -1
    nType = floor(nMolTypes * grnd() + 1E0_dp)
    if(trialBox%NMol(nType) + 1 > trialBox%NMolMax(nType)) then
      accept = .false.
      if(present(ProbIn)) ProbIn = 0E0_dp
      return
    endif

    nMove = trialBox%NMol(nType) + 1
    nMove = trialbox%MolGlobalIndx(nType, nMove)
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart, &
                               molEnd=molEnd)
    !Choose an atom to serve as the target for the new molecule.
!    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
!    call FindMolecule(trialbox, rawIndx, nTarget)
    nTarget = self%UniformMoleculeSelect(trialBox)
    call trialBox % GetMolData(nTarget, molStart=targStart, molEnd=targEnd)

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
                                                      self%insPoint, &
                                                      self%insProb)
    if(.not. accept) then
      if(present(ProbIn)) ProbIn = 0E0_dp
      return
    endif
    GenProb = ProbSub

    do iAtom = 1, nAtoms
      call trialBox % NeighList(1) % GetNewList(iAtom, self%tempList, self%tempNNei, &
                                                self%newPart(iAtom))
      self%newPart(iAtom)%listIndex = iAtom
    enddo 

    !Check Constraint
    accept = trialBox % CheckConstraint( self%newPart(1:nAtoms) )
    if(.not. accept) then
      if(present(ProbIn)) ProbIn = 0E0_dp
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
      if(present(ProbIn)) ProbIn = 0E0_dp
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%newPart(1:nAtoms), E_Diff )
    if(.not. accept) then
      if(present(ProbIn)) ProbIn = 0E0_dp
      return
    endif

    !Compute the probability of reversing the move
    ProbSel = 0E0_dp
    call self%SelectNeighReverse(trialBox, nTarget, ProbSel)

    !Compute the generation probability
    call MolData(molType) % molConstruct % GasConfig(GasProb)

    Prob = real(trialBox%nMolTotal, dp) * self%avbmcVol
    Prob = Prob*ProbSel/real(trialBox%nMolTotal+1, dp)
    Prob = GasProb*Prob/GenProb
!    write(*,*) "In", Prob, real(trialBox%nMolTotal, dp), 1E0_dp/ProbSel, GenProb
!    write(*,*) "  ", real(trialBox%nMolTotal+1, dp), self%avbmcVol, GasProb, molType

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
      if( present(targid) )  targid = nTarget
      if( present(ProbIn) )  ProbIn = Prob 
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%newPart(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine AVBMC_SwapOut(self, trialBox, accept, forcetargid, forceid, ProbOut) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData
    use Box_Utility, only: FindMolecule
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer, intent(in), optional :: forcetargid, forceid
    real(dp), intent(out), optional :: ProbOut
    integer :: iIns, nInsPoints
    integer :: molType, molStart, molEnd, targStart, targEnd
    integer :: nTarget, nMove, rawIndx, iConstrain
    integer :: CalcIndex, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, extraTerms
    real(dp) :: E_Inter, E_Intra
    real(dp) :: ProbSel = 1E0_dp
    real(dp) :: Prob = 1E0_dp
    real(dp) :: GenProb = 1E0_dp
    real(dp) :: GasProb = 1E0_dp
    real(dp) :: ProbSub = 1E0_dp

    integer :: slice(1:2)
    real(dp), pointer :: targetatoms(:,:) => null()

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialbox, self%oldPart)

    !Propose move
    call FindMolecule(trialbox, rawIndx, nTarget)
    if(.not. present(forcetargid) ) then
!      rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
!      call FindMolecule(trialbox, rawIndx, nMove)
      nTarget = self%UniformMoleculeSelect(trialBox)
    else
      nTarget = forcetargid
    endif

    call trialBox % GetMolData(nTarget, molType=molType, molStart=molStart)
    if( present(forceid) ) then
      call self%SelectNeigh(trialBox, nTarget, nMove, ProbSel, forceid)
    else
      call self%SelectNeigh(trialBox, nTarget, nMove, ProbSel)
    endif

    if(nMove < 1) then
      accept = .false.
      if(present(ProbOut)) ProbOut = 0E0_dp
      return
    endif
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart)


    if(trialBox%NMol(molType) - 1 < trialBox%NMolMin(molType)) then
      accept = .false.
      if(present(ProbOut)) ProbOut = 0E0_dp
      return
    endif


    self%oldPart(1)%molType = molType
    self%oldPart(1)%molIndx = nMove

    !Check Constraint
    accept = trialBox % CheckConstraint( self%oldPart(1:1) )
    if(.not. accept) then
      if(present(ProbOut)) ProbOut = 0E0_dp
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
      if(present(ProbOut)) ProbOut = 0E0_dp
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%oldPart(1:1), E_Diff )
    if(.not. accept) then
      if(present(ProbOut)) ProbOut = 0E0_dp
      return
    endif

    call trialBox % GetMolData(nTarget, molStart=targStart, molEnd=targEnd)
    slice(1) = targStart
    slice(2) = targEnd
    call trialbox%GetCoordinates(targetatoms, slice=slice)


    nInsPoints = MolData(molType) % molConstruct % GetNInsertPoints()
    call self%AllocateProb(nInsPoints)
    do iIns = 1, nInsPoints
      call self%CreateForward(targetatoms, self%inspoint(1:3, iIns))
      self%insprob(iIns) = 1E0_dp 
    enddo


    call MolData(molType) % molConstruct % ReverseConfig(self%oldpart(1:1), &
                                                         trialBox, &
                                                         ProbSub, &
                                                         accept, &
                                                         self%insPoint(1:3, 1:nInsPoints), &
                                                         self%insProb(1:nInsPoints)) 
    GenProb = ProbSub

    call MolData(molType) % molConstruct % GasConfig(GasProb)

    !Compute Acceptance Probability  
    Prob = real(trialBox%nMolTotal, dp)
    Prob = Prob/(real(trialBox%nMolTotal-1, dp) * self%avbmcVol * ProbSel)
    Prob = Prob*GenProb/GasProb
!    write(*,*) "Out", Prob, real(trialBox%nMolTotal, dp), 1E0_dp/ProbSel, GenProb
!    write(*,*) "   ", real(trialBox%nMolTotal-1, dp), self%avbmcVol
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
  subroutine AVBMC_SelectNeigh(self, trialBox, targetIndx, removeMolIndx, ProbOut, forceid)
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: ListRNG, grnd
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    integer, intent(in) :: targetIndx
    integer, intent(out) :: removeMolIndx
    real(dp), intent(out) :: ProbOut
    integer, intent(in), optional :: forceid
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
        call trialbox%Boundary(rx,ry,rz)
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
      if(present(forceid)) then
        removeMolIndx = forceid
      else
        removeIndx = floor(grnd()*nNei + 1E0_dp)
        removeMolIndx = removelist(removeIndx)
      endif
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
        call trialbox%Boundary(rx,ry,rz)
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

      !Note: Because we add the newly inserted molecule into the removelist() above
      !      we've already accounted for the (N_nei+1) factor in the traditional AVBMC
      !      formula.  Thus for this function we only need to return the 1/N_nei as
      !      computed by this function. 
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
    use Common_MolInfo, only: MolData, nMolTypes, mostAtoms
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


    call self%CreateTempArray(mostAtoms)
    self%avbmcVol = (4E0_dp/3E0_dp)*pi*self%avbmcRad**3
    self%avbmcRadSq = self%avbmcRad * self%avbmcRad
    write(nout,*) "AVBMC Radius:", self%avbmcRad
    write(nout,*) "AVBMC Volume:", self%avbmcVol

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
    use ClassyConstants, only: pi
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
        self%avbmcVol = (4E0_dp/3E0_dp)*pi*self%avbmcRad**3

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine

!========================================================
end module
!========================================================
