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
    use Common_MolInfo, only: MolData, nMolTypes
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


    allocate( self%tempNNei(nMaxPerMol) )
    allocate( self%tempList(nMaxPerBox, nMaxPerMol) )
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
    real(dp) :: ProbSub, ProbSel, GenProb, GasProb

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
      return
    endif

    nMove = trialBox%NMol(nType) + 1
    nMove = trialbox%MolGlobalIndx(nType, nMove)
    call trialBox % GetMolData(nMove, molType=molType, molStart=molStart, &
                               molEnd=molEnd)
    !Choose an atom to serve as the target for the new molecule.
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nTarget)

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
      self%newPart(iAtom)%molType = molType
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

    !Compute the probability of reversing the move
    call self%SelectNeighReverse(trialBox, nTarget, ProbSel)

    !Compute the generation probability
    Prob = real(trialBox%nMolTotal, dp) * self%avbmcVol
    Prob = Prob*ProbSel/real(trialBox%nMolTotal+1, dp)
    Prob = Prob/ProbSub

    call MolData(molType) % molConstruct % GasConfig(GasProb)
    Prob = GasProb*Prob/GenProb

    !Get the chemical potential term for GCMC
    extraTerms = sampling % GetExtraTerms(self%newpart(1:nAtoms), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%newPart(1:nAtoms), inProb=Prob, extraIn=extraTerms)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
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
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(AVBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: iIns, nInsPoints
    integer :: molType, molStart, molEnd, targStart, targEnd
    integer :: nTarget, nMove, rawIndx, iConstrain
    integer :: CalcIndex, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, extraTerms
    real(dp) :: E_Inter, E_Intra
    real(dp) :: ProbSel = 1E0_dp
    real(dp) :: Prob = 1E0_dp
    real(dp) :: ProbSub, GenProb, GasProb

    integer :: slice(1:2)
    real(dp), pointer :: targetatoms(:,:) => null()

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialbox, self%oldPart)

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
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%oldPart(1:1), E_Diff )
    if(.not. accept) then
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



    call MolData(molType) % molConstruct % GasConfig(GasProb)

    !Compute Acceptance Probability  
    Prob = real(trialBox%nMolTotal, dp)/ProbSel
    Prob = Prob/(real(trialBox%nMolTotal-1, dp) * self%avbmcVol)
    Prob = Prob*GenProb/GasProb
    !Get chemical potential term
    extraTerms = sampling % GetExtraTerms(self%oldpart(1:1), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff,  self%oldPart(1:1), inProb=Prob, extraIn=extraTerms)

    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
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
    write(nout,*) "AVBMC Radius:", self%avbmcRad
    write(nout,*) "AVBMC Volume:", self%avbmcVol

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
