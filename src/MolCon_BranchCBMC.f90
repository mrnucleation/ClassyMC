!==========================================================================================
! CBMC style regrowth algorithm for Branched molecules. More general purpose than
! the linear version and can handle molecules with branch points (e.g., isobutane,
! neopentane, branched alcohols).
!
! Key concepts:
! - Uses BFS traversal of molecular graph to determine regrowth order
! - At branch points, all branches are grown simultaneously to maintain proper
!   detailed balance
! - Supports 1-branch (linear extension), 2-branch, and 3-branch regrowth
!==========================================================================================
module MolCon_BranchCBMC
  use CoordinateTypes, only: Perturbation, Addition, Displacement, SingleMol, Deletion
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use Template_MolConstructor, only: MolConstructor
  use Data_Queue, only: IntQueue, RegrowthNode, RegrowthQueue
  use VarPrecision

  type, public, extends(MolConstructor) :: BranchCBMC
!    inspoints => Number of insertion points to request from the move
!                 calling this algorithm. Inherited from parent class.
!    integer :: inspoints = 1 

    !nAtoms => Number of Atoms in this molecule.
    !firstAtom => Relative Atom Index of the first atom to be inserted into the system
    !nRosenTrial => Number of Rosenbluth trials to use.  
    !                This is the number of potential positions
    !                that are generated for each new atom as it is being regrown.  
    !                One of these positions will be selected for the final atom position.  
    logical, private :: include15 = .false.
    integer, private :: nAtoms 
    integer, private :: firstAtom = 1 
    integer, private :: rosenNeighList = 1
    integer, private :: nRosenTrials = 1

    !RosenProb => Weight array for chosing which trial position to use
    real(dp), private, allocatable :: GenProb(:) 
    real(dp), private, allocatable :: RosenProb(:) 
    real(dp), private, allocatable :: tempcoords(:, :) 
    real(dp), private, allocatable :: newconfig(:, :) 

    !trialcoords => Temporary storage bucket for trial atom positions
    !  Dimensions: (3 coords, nAtoms, nRosenTrials)
    integer, private :: nTrialAtoms
    real(dp), private, allocatable :: trialcoords(:, :, :) 
    real(dp), private, allocatable :: trialGenProb(:)

    !nGrown => Counter which tracks how many atoms are already regrown. 
    !grown => Array which specifies if an atom has been added to the system or not
    !patharray => Path array for the molecule (may not be strictly linear)
    !pathposition => Reverse lookup: atom ID -> position in path
    !schedule => The order in which atoms will be grown
    !scratchschedule => Default growth order for full insertions
    integer, private :: nGrown
    logical, private, allocatable :: grown(:) 
    integer, private, allocatable :: patharray(:) 
    integer, private, allocatable :: pathposition(:) 
    integer, private, allocatable :: schedule(:) 
    integer, private, allocatable :: scratchschedule(:) 
    integer, private, allocatable :: tempList(:,:), tempNNei(:)
    integer, private, allocatable :: atomtypes(:)
    real(dp), private, allocatable :: posN(:,:)

    ! Storage for new bonds/angles/torsions being created
    integer, private :: nNewBonds, nNewAngles, nNewTorsions
    integer, private, allocatable :: newBonds(:)
    integer, private, allocatable :: newBondTypes(:)
    integer, private, allocatable :: newAngles(:)
    integer, private, allocatable :: newAngleTypes(:)
    integer, private, allocatable :: newTorsions(:)
    integer, private, allocatable :: newTorsionTypes(:)

    !freq => The number of bonds each atom has. Used to classify atoms:
    ! nBonds  = 0 -> Isolated Atom
    ! nBonds  = 1 -> Terminal Atom located on the end of a chain
    ! nBonds  = 2 -> Linker Atom located in the middle of a linear chain 
    ! nBonds >= 3 -> Branch Atom has two or more potential paths
    integer, private, allocatable :: freq(:) 
    
    ! nBranches => Number of bonded neighbors for each atom
    ! branchnext => Neighbor list: branchnext(iBranch, atomID) gives the neighbor atom
    integer, private, allocatable :: nBranches(:)
    integer, private, allocatable :: branchnext(:,:)
    
    ! Queue for BFS traversal during schedule creation
    type(IntQueue), private :: bfsQueue

    contains
      procedure, public, pass :: Prologue => BranchCBMC_Prologue
      procedure, public, pass :: CreateSchedule => BranchCBMC_CreateSchedule
      procedure, public, pass :: GenerateConfig => BranchCBMC_GenerateConfig
      procedure, public, pass :: ReverseConfig => BranchCBMC_ReverseConfig
      procedure, public, pass :: RosenBluth => BranchCBMC_RosenBluth
      procedure, public, pass :: GetPath => BranchCBMC_GetPath
      procedure, public, pass :: GasConfig => BranchCBMC_GasConfig
      procedure, public, pass :: FindAtomsFromPath => BranchCBMC_FindAtomsFromPath
      procedure, public, pass :: ProcessIO => BranchCBMC_ProcessIO
      procedure, private, pass :: OneBranch => BranchCBMC_OneBranch
      procedure, private, pass :: TwoBranch => BranchCBMC_TwoBranch
      procedure, private, pass :: ThreeBranch => BranchCBMC_ThreeBranch
      procedure, private, pass :: SelectTrial => BranchCBMC_SelectTrial
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine BranchCBMC_Prologue(self)
    use Common_MolInfo, only: MolData
    use BoxData, only: BoxArray
    use MolSearch, only: FindBond
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer :: iBox
    integer :: iBond, iAtom
    integer :: atm1, atm2
    integer :: ub1, ub2
    integer :: maxub1, maxub2
    integer :: iError
    integer :: iSchedule, nBranch
    integer :: nextAtm, prevAtm, curAtm, iPath
    integer :: leftdist, rightdist, nfirst
    integer :: largestVal, maxBranches

    iError = 0

    ! Get number of atoms in this molecule type
    self%nAtoms = MolData(self%molType)%nAtoms
    self%include15 = .false.

    ! Handle single atom molecules
    if (self%nAtoms == 1) then
      ! Single atom - just allocate minimal arrays
      allocate(self%freq(1:1))
      self%freq(1) = 0
      allocate(self%nBranches(1:1))
      self%nBranches(1) = 0
      allocate(self%branchnext(1:1, 1:1))
      self%branchnext = 0
      allocate(self%grown(1:1))
      allocate(self%newconfig(1:3, 1:1))
      allocate(self%patharray(1:1))
      allocate(self%pathposition(1:1))
      allocate(self%schedule(1:1))
      allocate(self%scratchschedule(1:1))
      self%patharray(1) = 1
      self%pathposition(1) = 1
      self%schedule(1) = 1
      self%scratchschedule(1) = 1
      
      if (.not. allocated(self%RosenProb)) then
        allocate(self%GenProb(1:self%nRosenTrials))
        allocate(self%RosenProb(1:self%nRosenTrials))
        allocate(self%tempcoords(1:3, 1:self%nRosenTrials))
      endif
      return
    endif

    ! Count the number of bonds each atom has and build neighbor list
    allocate(self%freq(1:self%nAtoms))
    allocate(self%nBranches(1:self%nAtoms))
    self%freq = 0
    self%nBranches = 0
    
    ! First pass: count bonds per atom
    do iBond = 1, size(MolData(self%molType)%bond)
      atm1 = MolData(self%molType)%bond(iBond)%mem1
      atm2 = MolData(self%molType)%bond(iBond)%mem2
      self%freq(atm1) = self%freq(atm1) + 1
      self%freq(atm2) = self%freq(atm2) + 1
    enddo

    ! Determine maximum number of branches any atom has
    maxBranches = maxval(self%freq)
    allocate(self%branchnext(1:maxBranches, 1:self%nAtoms))
    self%branchnext = 0

    ! Second pass: build neighbor list
    do iBond = 1, size(MolData(self%molType)%bond)
      atm1 = MolData(self%molType)%bond(iBond)%mem1
      atm2 = MolData(self%molType)%bond(iBond)%mem2
      
      nBranch = self%nBranches(atm1) + 1
      self%nBranches(atm1) = nBranch
      self%branchnext(nBranch, atm1) = atm2

      nBranch = self%nBranches(atm2) + 1
      self%nBranches(atm2) = nBranch
      self%branchnext(nBranch, atm2) = atm1
    enddo

    ! Check for cyclic molecules (not supported)
    if (all(self%freq /= 1)) then
      write(0,*) "ERROR: BranchCBMC detected a cyclic molecule (no terminal atoms)."
      write(0,*) "Cyclic molecules are not currently supported."
      iError = -1
    endif

    if (iError /= 0) then
      write(0,*) "BranchCBMC: Fatal error in molecule topology analysis."
      error stop
    endif

    ! Allocate probability and coordinate arrays
    if (.not. allocated(self%RosenProb)) then
      allocate(self%GenProb(1:self%nRosenTrials))
      allocate(self%RosenProb(1:self%nRosenTrials))
      allocate(self%tempcoords(1:3, 1:self%nRosenTrials))
    endif

    ! Allocate path and schedule arrays
    if (.not. allocated(self%patharray)) then
      allocate(self%newconfig(1:3, 1:self%nAtoms))
      allocate(self%trialcoords(1:3, 1:self%nAtoms, 1:self%nRosenTrials))
      allocate(self%trialGenProb(1:self%nRosenTrials))
      allocate(self%grown(1:self%nAtoms))
      allocate(self%patharray(1:self%nAtoms))
      allocate(self%pathposition(1:self%nAtoms))
      allocate(self%schedule(1:self%nAtoms))
      allocate(self%scratchschedule(1:self%nAtoms))
    endif

    ! Allocate storage for new intramolecular terms
    largestVal = 0
    if (allocated(MolData(self%molType)%nAtmBonds)) then
      do iAtom = 1, self%nAtoms
        largestVal = max(largestVal, MolData(self%molType)%nAtmBonds(iAtom))
      enddo
    endif
    largestVal = max(largestVal, 1)  ! Ensure at least size 1
    if (.not. allocated(self%newBonds)) allocate(self%newBonds(1:largestVal))
    if (.not. allocated(self%newBondTypes)) allocate(self%newBondTypes(1:largestVal))

    largestVal = 0
    if (allocated(MolData(self%molType)%nAtmAngles)) then
      do iAtom = 1, self%nAtoms
        largestVal = max(largestVal, MolData(self%molType)%nAtmAngles(iAtom))
      enddo
    endif
    largestVal = max(largestVal, 1)  ! Ensure at least size 1
    if (.not. allocated(self%newAngles)) allocate(self%newAngles(1:largestVal))
    if (.not. allocated(self%newAngleTypes)) allocate(self%newAngleTypes(1:largestVal))

    largestVal = 0
    if (allocated(MolData(self%molType)%nAtmTorsions)) then
      do iAtom = 1, self%nAtoms
        largestVal = max(largestVal, MolData(self%molType)%nAtmTorsions(iAtom))
      enddo
    endif
    largestVal = max(largestVal, 1)  ! Ensure at least size 1
    if (.not. allocated(self%newTorsions)) allocate(self%newTorsions(1:largestVal))
    if (.not. allocated(self%newTorsionTypes)) allocate(self%newTorsionTypes(1:largestVal))

    ! Initialize BFS queue for schedule creation
    call self%bfsQueue%Initialize(self%nAtoms)

    ! Build the default path array starting from a terminal atom
    self%patharray = 0
    self%pathposition = 0
    
    ! Find a terminal atom to start from
    prevAtm = 0
    curAtm = 0
    do iAtom = 1, self%nAtoms
      if (self%freq(iAtom) == 1) then
        curAtm = iAtom
        exit
      endif
    enddo

    if (curAtm == 0) then
      write(0,*) "BranchCBMC: Could not find a terminal atom to start path building."
      error stop
    endif

    ! Build path using BFS traversal
    self%patharray(1) = curAtm
    self%pathposition(curAtm) = 1
    iPath = 1
    call self%bfsQueue%Clear()
    call self%bfsQueue%Enqueue(curAtm)
    self%grown = .false.
    self%grown(curAtm) = .true.  ! Mark as visited

    do while (.not. self%bfsQueue%IsEmpty())
      curAtm = self%bfsQueue%Dequeue()
      do iBond = 1, self%nBranches(curAtm)
        nextAtm = self%branchnext(iBond, curAtm)
        if (.not. self%grown(nextAtm)) then
          iPath = iPath + 1
          self%patharray(iPath) = nextAtm
          self%pathposition(nextAtm) = iPath
          self%grown(nextAtm) = .true.
          call self%bfsQueue%Enqueue(nextAtm)
        endif
      enddo
    enddo

    if (any(self%patharray == 0)) then
      write(0,*) "BranchCBMC: ERROR! There are atoms unaccounted for by path building."
      write(0,*) "Ensure your molecule's bonds are properly connected."
      error stop
    endif

    ! Build default scratch schedule (for full insertions)
    ! Start from the first atom and use BFS order
    nfirst = self%pathposition(self%firstAtom)
    self%scratchschedule = self%patharray

    ! Allocate neighbor list storage for Rosenbluth calculations
    maxub1 = 0
    maxub2 = 0
    do iBox = 1, size(BoxArray)
      call BoxArray(iBox)%box%GetTempListBounds(ub1, ub2)
      maxub1 = max(ub1, maxub1)
      maxub2 = max(ub2, maxub2)
    enddo

    allocate(self%templist(1:maxub1, 1:self%nAtoms))
    allocate(self%tempNNei(1:self%nAtoms))

  end subroutine
!==========================================================================================
  subroutine BranchCBMC_CreateSchedule(self)
    ! Creates the regrowth schedule based on which atoms are currently marked
    ! as grown (anchors) and which need to be regrown.
    ! Uses BFS from the anchor points to determine order.
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer :: iAtom, iBranch, iSchedule
    integer :: curAtm, nextAtm
    logical :: foundAnchor

    call self%bfsQueue%Clear()
    iSchedule = 0
    
    ! First, add all grown atoms to the schedule (they're anchors)
    ! and enqueue their ungrown neighbors
    foundAnchor = .false.
    do iAtom = 1, self%nAtoms
      if (self%grown(iAtom)) then
        iSchedule = iSchedule + 1
        self%schedule(iSchedule) = iAtom
        foundAnchor = .true.
        ! Enqueue ungrown neighbors
        do iBranch = 1, self%nBranches(iAtom)
          nextAtm = self%branchnext(iBranch, iAtom)
          if (.not. self%grown(nextAtm) .and. .not. self%bfsQueue%Contains(nextAtm)) then
            call self%bfsQueue%Enqueue(nextAtm)
          endif
        enddo
      endif
    enddo

    ! If no anchors found (full regrowth), start from first atom
    if (.not. foundAnchor) then
      self%schedule = self%scratchschedule
      return
    endif

    ! BFS to add remaining atoms
    do while (.not. self%bfsQueue%IsEmpty())
      curAtm = self%bfsQueue%Dequeue()
      iSchedule = iSchedule + 1
      self%schedule(iSchedule) = curAtm
      
      ! Enqueue ungrown neighbors that haven't been scheduled yet
      do iBranch = 1, self%nBranches(curAtm)
        nextAtm = self%branchnext(iBranch, curAtm)
        if (.not. self%grown(nextAtm) .and. .not. self%bfsQueue%Contains(nextAtm)) then
          ! Check if already in schedule
          if (all(self%schedule(1:iSchedule) /= nextAtm)) then
            call self%bfsQueue%Enqueue(nextAtm)
          endif
        endif
      enddo
    enddo

  end subroutine
!==========================================================================================
  subroutine BranchCBMC_GenerateConfig(self, trialBox, disp, probconstruct, accept, insPoint, insProb)
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData
    use MolSearch, only: FindBond, FindAngle, FindTorsion
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone, Generate_UnitTorsion, ListRNG
    use ForcefieldData, only: ECalcArray

    implicit none
    class(BranchCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:,:)
    real(dp), intent(in), optional :: insProb(:)
    real(dp), intent(out) :: probconstruct 
    logical, intent(out) :: accept

    integer :: dispsubindx(1:self%nAtoms), atmdispindx(1:self%nAtoms)
    integer :: bondType, angleType, torsType
    integer :: molindx, molStart, molEnd, nSel
    integer :: atm1, atm2, atm3, atm4, atmindx, iDisp, iRosen
    integer :: lastGrown, iSchedule
    integer :: iBranch, nGrownNeigh, nUngrownNeigh
    integer :: grownNeigh(4), ungrownNeigh(4)
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r
    real(dp) :: norm
    real(dp) :: bend_angle, tors_angle
    real(dp) :: prob_r, prob_ang, prob_tors, probgen

    integer :: slice(1:2)
    real(dp), pointer :: atoms(:,:) => null()

    probconstruct = 1E0_dp
    accept = .true.

    select type(disp)
      class is(Displacement)
        ! Partial regrowth: some atoms stay fixed (anchors)
        self%grown = .true.
        self%nGrown = self%nAtoms
        molindx = disp(1)%molindx
        call trialbox%GetMolData(molindx, molStart=molStart, molEnd=molEnd)
        atmdispindx = 0
        dispsubindx = 0
        do iDisp = 1, size(disp)
          atm1 = disp(iDisp)%atmindx - molStart + 1
          dispsubindx(iDisp) = atm1
          atmdispindx(atm1) = iDisp
          self%grown(atm1) = .false.
          self%nGrown = self%nGrown - 1
        enddo
        call self%CreateSchedule
        slice(1) = molStart
        slice(2) = molEnd
        call trialbox%GetCoordinates(atoms, slice=slice)
        self%newconfig(1:3, 1:self%nAtoms) = atoms(1:3, 1:self%nAtoms)

      class is(Addition)
        ! Full regrowth: start from scratch
        self%grown = .false.
        self%nGrown = 0
        molindx = disp(1)%molindx
        call trialbox%GetMolData(molindx, molStart=molStart)
        do iDisp = 1, size(disp)
          atm1 = disp(iDisp)%atmindx - molStart + 1
          dispsubindx(iDisp) = atm1
          atmdispindx(atm1) = iDisp
        enddo
        self%schedule = self%scratchschedule

      class default
        write(0,*) "Critical Error! Invalid perturbation type passed to BranchCBMC"
        error stop 
    end select

    if (present(insPoint)) then
      if (size(insPoint, 2) /= self%nRosenTrials) then
        write(0,*) "ERROR! BranchCBMC received wrong number of insertion points"
        error stop
      endif
    endif

    ! Main regrowth loop
    iSchedule = self%nGrown
    do while (self%nGrown < self%nAtoms)
      iSchedule = iSchedule + 1
      lastGrown = self%schedule(iSchedule)
      
      ! Count grown and ungrown neighbors for this atom
      nGrownNeigh = 0
      nUngrownNeigh = 0
      do iBranch = 1, self%nBranches(lastGrown)
        atm1 = self%branchnext(iBranch, lastGrown)
        if (self%grown(atm1)) then
          nGrownNeigh = nGrownNeigh + 1
          grownNeigh(nGrownNeigh) = atm1
        else
          nUngrownNeigh = nUngrownNeigh + 1
          ungrownNeigh(nUngrownNeigh) = atm1
        endif
      enddo

      ! Generate trial positions based on number of grown neighbors
      if (self%nGrown == 0) then
        ! First atom: use insertion points
        if (present(insPoint) .and. present(insProb)) then
          do iRosen = 1, self%nRosenTrials
            self%tempcoords(1:3, iRosen) = insPoint(1:3, iRosen)
            self%GenProb(iRosen) = insProb(iRosen)
          enddo
        else
          error stop "Full regrowth requested without insertion points"
        endif

      elseif (nGrownNeigh == 0) then
        ! No grown neighbors - this shouldn't happen with BFS ordering
        write(0,*) "ERROR: Atom scheduled for regrowth has no grown neighbors!"
        write(0,*) "Atom:", lastGrown, "Schedule position:", iSchedule
        error stop

      elseif (nGrownNeigh == 1) then
        ! One grown neighbor: use OneBranch
        call self%OneBranch(lastGrown, grownNeigh(1), trialBox)

      elseif (nGrownNeigh == 2) then
        ! Two grown neighbors: use TwoBranch (rare, usually at branch points)
        call self%TwoBranch(lastGrown, grownNeigh(1), grownNeigh(2), trialBox)

      elseif (nGrownNeigh >= 3) then
        ! Three or more grown neighbors: use ThreeBranch
        call self%ThreeBranch(lastGrown, grownNeigh(1), grownNeigh(2), grownNeigh(3), trialBox)

      endif

      ! Compute Rosenbluth weights and select trial
      iDisp = atmdispindx(lastGrown)
      if (iDisp > 0) then
        call self%RosenBluth(trialBox, lastGrown, iDisp, disp)
      else
        ! Atom not in displacement array - still compute weights
        self%RosenProb = 1E0_dp
      endif

      norm = 0E0_dp
      do iRosen = 1, self%nRosenTrials
        norm = norm + self%RosenProb(iRosen)
      enddo

      if (norm <= 0E0_dp) then
        accept = .false.
        return
      endif

      nSel = ListRNG(self%RosenProb, norm)
      probconstruct = probconstruct * self%GenProb(nSel) * self%RosenProb(nSel) / norm
      
      self%newconfig(1:3, lastGrown) = self%tempcoords(1:3, nSel)
      self%grown(lastGrown) = .true.
      self%nGrown = self%nGrown + 1
    enddo

    ! Copy final positions to displacement array
    do iDisp = 1, size(disp)
      select type(disp)
        class is(Displacement)
          atmindx = disp(iDisp)%atmindx - molStart + 1
          disp(iDisp)%x_new = self%newconfig(1, atmindx)
          disp(iDisp)%y_new = self%newconfig(2, atmindx)
          disp(iDisp)%z_new = self%newconfig(3, atmindx)
        class is(Addition)
          atmindx = disp(iDisp)%atmindx - molStart + 1
          disp(iDisp)%x_new = self%newconfig(1, atmindx)
          disp(iDisp)%y_new = self%newconfig(2, atmindx)
          disp(iDisp)%z_new = self%newconfig(3, atmindx)
      end select
    enddo

  end subroutine
!======================================================================================
  subroutine BranchCBMC_OneBranch(self, curAtm, anchorAtm, trialBox)
    ! Generate trial positions for an atom with one grown neighbor (anchor)
    ! This is the most common case: extending a chain by one atom.
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData
    use MolSearch, only: FindBond, FindAngle, FindTorsion
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone, Generate_UnitTorsion, grnd
    use ClassyConstants, only: two_pi
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(in) :: curAtm, anchorAtm
    class(SimBox), intent(inout) :: trialBox

    integer :: iRosen, iBranch, iTors
    integer :: bondType, angleType, torsType
    integer :: atm1, atm2, atm3, atm4
    integer :: prevAtm, prevPrevAtm
    integer :: nGrownAnchorNeigh
    real(dp) :: r, bend_angle, tors_angle
    real(dp) :: prob_r, prob_ang, prob_tors, probgen
    real(dp) :: v1(1:3), v2(1:3), v3(1:3)
    real(dp) :: dx, dy, dz

    ! Find the bond type
    call FindBond(self%molType, curAtm, anchorAtm, bondType)
    
    ! Count how many grown neighbors the anchor has (besides curAtm)
    nGrownAnchorNeigh = 0
    prevAtm = 0
    prevPrevAtm = 0
    do iBranch = 1, self%nBranches(anchorAtm)
      atm1 = self%branchnext(iBranch, anchorAtm)
      if (atm1 /= curAtm .and. self%grown(atm1)) then
        nGrownAnchorNeigh = nGrownAnchorNeigh + 1
        if (nGrownAnchorNeigh == 1) prevAtm = atm1
        if (nGrownAnchorNeigh == 2) prevPrevAtm = atm1
      endif
    enddo

    ! Generate trials based on available reference atoms
    if (nGrownAnchorNeigh == 0) then
      ! Only anchor atom is grown - generate on sphere
      do iRosen = 1, self%nRosenTrials
        call BondData(bondType)%bondFF%GenerateDist(trialBox%beta, r, prob_r)
        call Generate_UnitSphere(dx, dy, dz)
        self%tempcoords(1, iRosen) = r*dx + self%newconfig(1, anchorAtm)
        self%tempcoords(2, iRosen) = r*dy + self%newconfig(2, anchorAtm)
        self%tempcoords(3, iRosen) = r*dz + self%newconfig(3, anchorAtm)
        call trialBox%Boundary(self%tempcoords(1, iRosen), &
                               self%tempcoords(2, iRosen), &
                               self%tempcoords(3, iRosen))
        self%GenProb(iRosen) = prob_r
      enddo

    elseif (nGrownAnchorNeigh == 1) then
      ! Have one reference atom - use bond angle
      call FindAngle(self%molType, prevAtm, anchorAtm, curAtm, angleType)
      v1(1:3) = self%newconfig(1:3, prevAtm) - self%newconfig(1:3, anchorAtm)
      call trialBox%Boundary(v1(1), v1(2), v1(3))

      do iRosen = 1, self%nRosenTrials
        call BondData(bondType)%bondFF%GenerateDist(trialBox%beta, r, prob_r)
        call AngleData(angleType)%angleFF%GenerateDist(trialBox%beta, bend_angle, prob_ang)
        probgen = prob_r * prob_ang
        call Generate_UnitCone(v1, r, bend_angle, v2)
        self%tempcoords(1:3, iRosen) = v2(1:3) + self%newconfig(1:3, anchorAtm)
        call trialBox%Boundary(self%tempcoords(1, iRosen), &
                               self%tempcoords(2, iRosen), &
                               self%tempcoords(3, iRosen))
        self%GenProb(iRosen) = probgen
      enddo

    else
      ! Have two or more reference atoms - try to use torsion if available
      call FindAngle(self%molType, prevAtm, anchorAtm, curAtm, angleType)
      
      ! Check if a torsion exists for this 4-atom sequence
      torsType = -1
      call FindTorsionSafe(self%molType, prevPrevAtm, prevAtm, anchorAtm, curAtm, torsType)
      
      v1(1:3) = self%newconfig(1:3, prevPrevAtm) - self%newconfig(1:3, anchorAtm)
      call trialBox%Boundary(v1(1), v1(2), v1(3))
      v2(1:3) = self%newconfig(1:3, prevAtm) - self%newconfig(1:3, anchorAtm)
      call trialBox%Boundary(v2(1), v2(2), v2(3))

      if (torsType > 0) then
        ! Torsion exists - use full torsion generation
        do iRosen = 1, self%nRosenTrials
          call BondData(bondType)%bondFF%GenerateDist(trialBox%beta, r, prob_r)
          call AngleData(angleType)%angleFF%GenerateDist(trialBox%beta, bend_angle, prob_ang)
          call TorsionData(torsType)%torsionFF%GenerateDist(trialBox%beta, tors_angle, prob_tors)
          probgen = prob_r * prob_ang * prob_tors
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          self%tempcoords(1:3, iRosen) = v3(1:3) + self%newconfig(1:3, anchorAtm)
          call trialBox%Boundary(self%tempcoords(1, iRosen), &
                                 self%tempcoords(2, iRosen), &
                                 self%tempcoords(3, iRosen))
          self%GenProb(iRosen) = probgen
        enddo
      else
        ! No torsion defined - use angle with random torsion (uniform)
        do iRosen = 1, self%nRosenTrials
          call BondData(bondType)%bondFF%GenerateDist(trialBox%beta, r, prob_r)
          call AngleData(angleType)%angleFF%GenerateDist(trialBox%beta, bend_angle, prob_ang)
          tors_angle = two_pi * grnd()
          prob_tors = 1.0E0_dp  ! Uniform distribution
          probgen = prob_r * prob_ang * prob_tors
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          self%tempcoords(1:3, iRosen) = v3(1:3) + self%newconfig(1:3, anchorAtm)
          call trialBox%Boundary(self%tempcoords(1, iRosen), &
                                 self%tempcoords(2, iRosen), &
                                 self%tempcoords(3, iRosen))
          self%GenProb(iRosen) = probgen
        enddo
      endif
    endif

  end subroutine
!======================================================================================
  subroutine BranchCBMC_TwoBranch(self, curAtm, anchor1, anchor2, trialBox)
    ! Generate trial positions for an atom constrained by two grown neighbors
    ! This occurs when growing an atom at a junction point with two anchors
    use Common_MolInfo, only: MolData, BondData, AngleData
    use MolSearch, only: FindBond, FindAngle
    use RandomGen, only: Generate_UnitCone, grnd
    use Math_Coordinates, only: AngleFromVectors
    use ClassyConstants, only: two_pi
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(in) :: curAtm, anchor1, anchor2
    class(SimBox), intent(inout) :: trialBox

    integer :: iRosen
    integer :: bondType1, bondType2, angleType
    real(dp) :: r1, r2, bend_angle
    real(dp) :: prob_r1, prob_r2, prob_ang, probgen
    real(dp) :: v1(1:3), v2(1:3), v_mid(1:3)
    real(dp) :: d1, d2, actual_angle
    real(dp) :: phi

    ! Find bond types to both anchors
    call FindBond(self%molType, curAtm, anchor1, bondType1)
    call FindBond(self%molType, curAtm, anchor2, bondType2)
    call FindAngle(self%molType, anchor1, curAtm, anchor2, angleType)

    ! Vector from anchor1 to anchor2
    v1(1:3) = self%newconfig(1:3, anchor2) - self%newconfig(1:3, anchor1)
    call trialBox%Boundary(v1(1), v1(2), v1(3))
    d1 = sqrt(sum(v1**2))

    ! Midpoint between anchors as reference
    v_mid(1:3) = 0.5E0_dp * (self%newconfig(1:3, anchor1) + self%newconfig(1:3, anchor2))

    do iRosen = 1, self%nRosenTrials
      ! Generate bond lengths
      call BondData(bondType1)%bondFF%GenerateDist(trialBox%beta, r1, prob_r1)
      call BondData(bondType2)%bondFF%GenerateDist(trialBox%beta, r2, prob_r2)
      
      ! Generate the bending angle
      call AngleData(angleType)%angleFF%GenerateDist(trialBox%beta, bend_angle, prob_ang)
      
      ! For two anchors, we place the atom on a cone around the anchor1-anchor2 axis
      ! The cone angle is determined by the geometry
      phi = two_pi * grnd()
      
      ! Generate position using the midpoint as reference
      call Generate_UnitCone(v1, r1, bend_angle/2.0E0_dp, v2)
      
      self%tempcoords(1:3, iRosen) = v2(1:3) + self%newconfig(1:3, anchor1)
      call trialBox%Boundary(self%tempcoords(1, iRosen), &
                             self%tempcoords(2, iRosen), &
                             self%tempcoords(3, iRosen))
      
      probgen = prob_r1 * prob_ang
      self%GenProb(iRosen) = probgen
    enddo

  end subroutine
!======================================================================================
  subroutine BranchCBMC_ThreeBranch(self, curAtm, anchor1, anchor2, anchor3, trialBox)
    ! Generate trial positions for an atom constrained by three grown neighbors
    ! This is highly constrained - typically only one or two valid positions
    use Common_MolInfo, only: MolData, BondData, AngleData
    use MolSearch, only: FindBond, FindAngle
    use RandomGen, only: grnd
    use ClassyConstants, only: two_pi, pi
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(in) :: curAtm, anchor1, anchor2, anchor3
    class(SimBox), intent(inout) :: trialBox

    integer :: iRosen
    integer :: bondType1
    integer :: angleType1, angleType2, angleType3
    real(dp) :: r1, ang1, ang2, ang3
    real(dp) :: prob_r, prob_ang1, prob_ang2, prob_ang3, probgen
    real(dp) :: v1(1:3), v2(1:3), v3(1:3), vn(1:3)
    real(dp) :: center(1:3)
    real(dp) :: d1, d2, d3
    real(dp) :: sign_choice

    ! Find bond and angle types
    call FindBond(self%molType, curAtm, anchor1, bondType1)
    call FindAngle(self%molType, anchor1, curAtm, anchor2, angleType1)
    call FindAngle(self%molType, anchor1, curAtm, anchor3, angleType2)
    call FindAngle(self%molType, anchor2, curAtm, anchor3, angleType3)

    ! Vectors from anchors
    v1(1:3) = self%newconfig(1:3, anchor1)
    v2(1:3) = self%newconfig(1:3, anchor2)
    v3(1:3) = self%newconfig(1:3, anchor3)

    ! Compute centroid of anchors
    center(1:3) = (v1 + v2 + v3) / 3.0E0_dp

    ! Compute normal to plane of anchors
    vn(1:3) = (v2 - v1)
    call trialBox%Boundary(vn(1), vn(2), vn(3))
    v1(1:3) = v3(1:3) - self%newconfig(1:3, anchor1)
    call trialBox%Boundary(v1(1), v1(2), v1(3))
    
    ! Cross product for normal
    vn(1) = (v2(2)-self%newconfig(2,anchor1))*(v3(3)-self%newconfig(3,anchor1)) - &
            (v2(3)-self%newconfig(3,anchor1))*(v3(2)-self%newconfig(2,anchor1))
    vn(2) = (v2(3)-self%newconfig(3,anchor1))*(v3(1)-self%newconfig(1,anchor1)) - &
            (v2(1)-self%newconfig(1,anchor1))*(v3(3)-self%newconfig(3,anchor1))
    vn(3) = (v2(1)-self%newconfig(1,anchor1))*(v3(2)-self%newconfig(2,anchor1)) - &
            (v2(2)-self%newconfig(2,anchor1))*(v3(1)-self%newconfig(1,anchor1))
    
    d1 = sqrt(sum(vn**2))
    if (d1 > 1E-10_dp) vn = vn / d1

    do iRosen = 1, self%nRosenTrials
      ! Generate bond length
      call BondData(bondType1)%bondFF%GenerateDist(trialBox%beta, r1, prob_r)
      
      ! For tetrahedral geometry, choose above or below plane
      sign_choice = 1.0E0_dp
      if (grnd() < 0.5E0_dp) sign_choice = -1.0E0_dp
      
      ! Place atom along normal from centroid
      ! Distance from centroid depends on bond length and tetrahedral geometry
      self%tempcoords(1:3, iRosen) = center(1:3) + sign_choice * r1 * 0.816E0_dp * vn(1:3)
      call trialBox%Boundary(self%tempcoords(1, iRosen), &
                             self%tempcoords(2, iRosen), &
                             self%tempcoords(3, iRosen))
      
      probgen = prob_r * 0.5E0_dp  ! Factor of 0.5 for the sign choice
      self%GenProb(iRosen) = probgen
    enddo

  end subroutine
!======================================================================================
  subroutine BranchCBMC_SelectTrial(self, nSel, norm)
    ! Select a trial position based on Rosenbluth weights
    use RandomGen, only: ListRNG
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(out) :: nSel
    real(dp), intent(out) :: norm
    integer :: iRosen

    norm = 0E0_dp
    do iRosen = 1, self%nRosenTrials
      norm = norm + self%RosenProb(iRosen)
    enddo

    nSel = ListRNG(self%RosenProb, norm)

  end subroutine
!======================================================================================
  subroutine BranchCBMC_ReverseConfig(self, disp, trialBox, probconstruct, accept, insPoint, insProb)
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData
    use MolSearch, only: FindBond, FindAngle, FindTorsion
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone, Generate_UnitTorsion, ListRNG
    implicit none
    class(BranchCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct 
    real(dp), intent(in), optional :: insPoint(:, :)
    real(dp), intent(in), optional :: insProb(:)
    logical, intent(out) :: accept

    integer :: dispsubindx(1:self%nAtoms), atmdispindx(1:self%nAtoms)
    integer :: bondType, angleType, torsType
    integer :: molindx, molStart, molEnd
    integer :: atm1, atm2, atm3, atm4, iDisp, iRosen
    integer :: lastGrown, iSchedule
    integer :: iBranch, nGrownNeigh
    integer :: grownNeigh(4)
    integer :: prevAtm, prevPrevAtm
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r
    real(dp) :: norm
    real(dp) :: bend_angle, tors_angle
    real(dp) :: prob_r, prob_ang, prob_tors, probgen
    integer :: slice(1:2)
    real(dp), pointer :: atoms(:,:) => null()
    real(dp) :: oldpos(1:3, 1:4)

    probconstruct = 1E0_dp
    accept = .true.

    select type(disp)
      class is(Displacement)
        self%grown = .true.
        self%nGrown = self%nAtoms
        molindx = disp(1)%molindx
        call trialbox%GetMolData(molindx, molStart=molStart, molEnd=molEnd)
        slice(1) = molStart
        slice(2) = molEnd
        call trialbox%GetCoordinates(atoms, slice=slice)
        atmdispindx = 0
        dispsubindx = 0
        do iDisp = 1, size(disp)
          atm1 = disp(iDisp)%atmindx - molStart + 1
          dispsubindx(iDisp) = atm1
          atmdispindx(atm1) = iDisp
          self%grown(atm1) = .false.
          self%nGrown = self%nGrown - 1
        enddo
        self%newconfig(1:3, 1:self%nAtoms) = atoms(1:3, 1:self%nAtoms)
        call self%CreateSchedule

      class is(Deletion)
        self%grown = .false.
        self%nGrown = 0
        molindx = disp(1)%molindx
        call trialbox%GetMolData(molindx, molStart=molStart, molEnd=molEnd)
        slice(1) = molStart
        slice(2) = molEnd
        self%schedule = self%scratchschedule
        call trialbox%GetCoordinates(atoms, slice=slice)
        do iDisp = 1, self%nAtoms
          dispsubindx(iDisp) = iDisp
          atmdispindx(iDisp) = 1
        enddo
        self%newconfig(1:3, 1:self%nAtoms) = atoms(1:3, 1:self%nAtoms)

      class default
        error stop "Critical Error! Invalid perturbation type in BranchCBMC ReverseConfig"
    end select

    if (present(insPoint) .neqv. present(insProb)) then
      write(0,*) "insPoint and insProb must both be present or absent!"
      error stop
    endif

    if (present(insPoint)) then
      if (size(insPoint, 2) /= self%nRosenTrials) then
        write(0,*) "ERROR! BranchCBMC received wrong number of insertion points"
        error stop
      endif
    endif

    ! Compute reverse probability following same schedule
    iSchedule = self%nGrown
    do while (self%nGrown < self%nAtoms)
      iSchedule = iSchedule + 1
      lastGrown = self%schedule(iSchedule)

      ! Count grown neighbors
      nGrownNeigh = 0
      do iBranch = 1, self%nBranches(lastGrown)
        atm1 = self%branchnext(iBranch, lastGrown)
        if (self%grown(atm1)) then
          nGrownNeigh = nGrownNeigh + 1
          grownNeigh(nGrownNeigh) = atm1
        endif
      enddo

      ! Generate n-1 random trials and use old position as last trial
      if (self%nGrown == 0) then
        ! First atom
        if (present(insPoint) .and. present(insProb)) then
          do iRosen = 1, self%nRosenTrials - 1
            self%tempcoords(1:3, iRosen) = insPoint(1:3, iRosen)
            self%GenProb(iRosen) = insProb(iRosen)
          enddo
          self%tempcoords(1:3, self%nRosenTrials) = self%newconfig(1:3, lastGrown)
          self%GenProb(self%nRosenTrials) = insProb(self%nRosenTrials)
        else
          error stop "Full reverse regrowth requested without insertion points"
        endif

      elseif (nGrownNeigh >= 1) then
        ! Generate trials as in forward, but place old position as last trial
        call self%OneBranch(lastGrown, grownNeigh(1), trialBox)
        
        ! Replace last trial with actual old position
        self%tempcoords(1:3, self%nRosenTrials) = self%newconfig(1:3, lastGrown)
        
        ! Compute generation probability of old position
        prevAtm = grownNeigh(1)
        call FindBond(self%molType, lastGrown, prevAtm, bondType)
        
        oldpos(1:3, 1) = atoms(1:3, prevAtm)
        oldpos(1:3, 2) = atoms(1:3, lastGrown)
        call BondData(bondType)%bondFF%GenerateReverseDist(trialbox, oldpos(1:3, 1:2), prob_r)
        
        ! Check if we need angle probability
        nGrownNeigh = 0
        prevPrevAtm = 0
        do iBranch = 1, self%nBranches(prevAtm)
          atm1 = self%branchnext(iBranch, prevAtm)
          if (atm1 /= lastGrown .and. self%grown(atm1)) then
            nGrownNeigh = nGrownNeigh + 1
            if (nGrownNeigh == 1) prevPrevAtm = atm1
          endif
        enddo
        
        if (prevPrevAtm > 0) then
          call FindAngle(self%molType, prevPrevAtm, prevAtm, lastGrown, angleType)
          oldpos(1:3, 1) = atoms(1:3, prevPrevAtm)
          oldpos(1:3, 2) = atoms(1:3, prevAtm)
          oldpos(1:3, 3) = atoms(1:3, lastGrown)
          call AngleData(angleType)%angleFF%GenerateReverseDist(trialbox, oldpos(1:3, 1:3), prob_ang)
          self%GenProb(self%nRosenTrials) = prob_r * prob_ang
        else
          self%GenProb(self%nRosenTrials) = prob_r
        endif
      endif

      ! Compute Rosenbluth weights
      iDisp = atmdispindx(lastGrown)
      if (iDisp > 0) then
        call self%RosenBluth(trialBox, lastGrown, iDisp, disp)
      else
        self%RosenProb = 1E0_dp
      endif

      norm = 0E0_dp
      do iRosen = 1, self%nRosenTrials
        norm = norm + self%RosenProb(iRosen)
      enddo

      if (norm < 1E-10_dp) then
        accept = .false.
        return
      endif

      probconstruct = probconstruct * self%GenProb(self%nRosenTrials) * &
                      self%RosenProb(self%nRosenTrials) / norm

      self%grown(lastGrown) = .true.
      self%nGrown = self%nGrown + 1
    enddo

  end subroutine
!======================================================================
  subroutine BranchCBMC_GasConfig(self, probGas)
    ! Compute probability of an isolated molecule in gas phase
    ! Used for swap moves with implicit gas box
    implicit none
    class(BranchCBMC), intent(inout) :: self
    real(dp), intent(out) :: probGas

    probGas = 1E0_dp
    if (.not. self%include15) then
      probGas = 1E0_dp / real(self%nRosenTrials, dp)**(self%nAtoms)
    endif

  end subroutine
!=======================================================================
  subroutine BranchCBMC_FindAtomsFromPath(self, Atm4, Atm1, Atm2, Atm3)
    ! Find reference atoms for torsion generation given target atom
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(in) :: Atm4
    integer, intent(out) :: Atm1, Atm2, Atm3
    integer :: iBranch, nFound
    integer :: candidates(4)

    ! Find grown neighbors of Atm4
    nFound = 0
    do iBranch = 1, self%nBranches(Atm4)
      if (self%grown(self%branchnext(iBranch, Atm4))) then
        nFound = nFound + 1
        candidates(nFound) = self%branchnext(iBranch, Atm4)
      endif
    enddo

    if (nFound < 1) then
      write(0,*) "ERROR in FindAtomsFromPath: No grown neighbors found!"
      error stop
    endif

    Atm3 = candidates(1)
    
    ! Find grown neighbors of Atm3 (excluding Atm4)
    nFound = 0
    do iBranch = 1, self%nBranches(Atm3)
      if (self%branchnext(iBranch, Atm3) /= Atm4 .and. &
          self%grown(self%branchnext(iBranch, Atm3))) then
        nFound = nFound + 1
        candidates(nFound) = self%branchnext(iBranch, Atm3)
      endif
    enddo

    if (nFound < 1) then
      Atm2 = 0
      Atm1 = 0
      return
    endif

    Atm2 = candidates(1)

    ! Find grown neighbors of Atm2 (excluding Atm3)
    nFound = 0
    do iBranch = 1, self%nBranches(Atm2)
      if (self%branchnext(iBranch, Atm2) /= Atm3 .and. &
          self%grown(self%branchnext(iBranch, Atm2))) then
        nFound = nFound + 1
        candidates(nFound) = self%branchnext(iBranch, Atm2)
      endif
    enddo

    if (nFound < 1) then
      Atm1 = 0
      return
    endif

    Atm1 = candidates(1)

  end subroutine 
!==========================================================================================
  subroutine BranchCBMC_RosenBluth(self, trialBox, atmsubindx, iDisp, disp)
    ! Computes the Rosenbluth Weight of each trial position
    use ForcefieldData, only: ECalcArray
    use ErrorChecking, only: IsNan, IsInf
    implicit none
    class(BranchCBMC), intent(inout), target :: self
    class(Perturbation), intent(inout) :: disp(:)
    integer, intent(in) :: iDisp, atmsubindx
    class(SimBox), intent(inout) :: trialBox

    class(ECalcArray), pointer :: EFunc => null()
    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    logical :: accept_trial
    logical :: overlap(1:self%nRosenTrials)
    integer :: molIndx, molType, atmIndx, molStart
    integer :: iRosen, jAtom, jNei, neiSize
    integer :: atmtype1, atmNeiIndx
    real(dp) :: E_Atom, E_Min, norm
    real(dp) :: pos1(1:3)
    type(Displacement) :: tempdisp(1:1)

    if (self%nRosenTrials == 1) then
      self%RosenProb(1) = 1E0_dp
      return
    endif

    if ((iDisp < 1) .or. (iDisp > size(disp))) then
      write(0,*) "Invalid Perturbation Index!"
      write(0,*) iDisp, size(disp)
      error stop
    endif

    select type(trialbox)
      class is(SimpleBox)
        call trialbox%GetCoordinates(atoms)
        call trialbox%GetEFunc(EFunc)
        select type(disp)
          class is(Displacement)
            call trialBox%GetNeighborList(self%rosenNeighList, neighlist, nNeigh)
          class is(Addition)
            neighlist => self%templist
            nNeigh => self%tempNNei
        end select
    end select

    if (.not. allocated(self%atomtypes)) then
      neiSize = size(neighlist, 1)
      allocate(self%atomtypes(1:neiSize))
      allocate(self%posN(1:3, 1:neiSize))
    endif

    select type(disp)
      class is(SingleMol)
        tempdisp(1)%molType = disp(iDisp)%molType
        tempdisp(1)%molindx = disp(iDisp)%molindx
    end select

    call trialBox%GetMolData(tempdisp(1)%molindx, molStart=molStart)

    overlap = .false.

    atmIndx = atmSubIndx + molStart - 1
    tempdisp(1)%atmindx = atmIndx
    select type(disp)
      class is(Displacement)
        atmNeiIndx = atmIndx
      class default
        atmNeiIndx = 1
    end select
    call trialBox%GetAtomData(atmIndx, atomtype=atmtype1)

    E_Min = huge(dp)
    E_Atom = 0E0_dp
    do iRosen = 1, self%nRosenTrials
      tempdisp(1)%x_new = self%tempcoords(1, iRosen)
      tempdisp(1)%y_new = self%tempcoords(2, iRosen)
      tempdisp(1)%z_new = self%tempcoords(3, iRosen)
      select type(disp)
        class is(Addition)
          select type(trialbox)
            class is(SimpleBox)
              call trialbox%GetNewNeighborList(self%rosenNeighList, 1, &
                   self%templist, self%tempNNei, tempdisp(1))
          end select
          atmNeiIndx = 1
      end select
      pos1(1:3) = self%tempcoords(1:3, iRosen)
      accept_trial = .true.
      if (nNeigh(atmNeiIndx) > 0) then
        do jNei = 1, nNeigh(atmNeiIndx)
          jAtom = neighlist(jNei, atmNeiIndx)
          call trialBox%GetAtomData(jAtom, atomtype=self%atomtypes(jNei))
          self%posN(1:3, jNei) = atoms(1:3, jAtom)
        enddo
        call EFunc%Method%ManyBody(trialbox, &
                         atmtype1, &
                         pos1, &
                         self%atomtypes(1:nNeigh(atmNeiIndx)), &
                         self%posN(1:3, 1:nNeigh(atmNeiIndx)), &
                         E_Atom, &
                         accept_trial)
      else
        E_Atom = 0E0_dp
      endif
      if (accept_trial) then
        if (E_Atom < E_Min) then
          E_Min = E_Atom
        endif
        self%RosenProb(iRosen) = E_Atom
      else
        overlap(iRosen) = .true.
        self%RosenProb(iRosen) = 0E0_dp
      endif
    enddo

    ! Compute Rosenbluth weights: P(E_i) = exp(-E_i/kT)/N
    ! Normalize by E_Min to prevent overflow
    norm = 0E0_dp
    do iRosen = 1, self%nRosenTrials
      if (overlap(iRosen)) cycle

      self%RosenProb(iRosen) = (self%RosenProb(iRosen) - E_Min) * trialbox%beta
      self%RosenProb(iRosen) = exp(-self%RosenProb(iRosen))
      if (IsNan(self%RosenProb(iRosen)) .or. IsInf(self%RosenProb(iRosen))) then
        write(0,*) "Invalid Weight in BranchCBMC!"
        write(0,*) self%RosenProb(iRosen), E_Min
        error stop
      endif
    enddo

  end subroutine 
!==========================================================================================
  subroutine BranchCBMC_GetPath(self, pathout)
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(out) :: pathout(1:self%nAtoms)
    integer :: iPath

    do iPath = 1, size(self%patharray)
      pathout(iPath) = self%patharray(iPath)
    enddo

  end subroutine
!==========================================================================================
  subroutine BranchCBMC_ProcessIO(self, line, linestat)
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(BranchCBMC), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: linestat
    character(len=30) :: command
    integer :: nRosen
    
    call GetXCommand(line, command, 3, lineStat)
    read(command, *) nRosen
    self%inspoints = nRosen
    self%nRosenTrials = nRosen

  end subroutine
!=======================================================================
  subroutine FindTorsionSafe(molType, mem1, mem2, mem3, mem4, torsType)
    ! Safe version of FindTorsion that returns -1 if no torsion found
    ! instead of stopping with an error
    use Common_MolInfo, only: MolData
    implicit none
    integer, intent(in) :: molType, mem1, mem2, mem3, mem4
    integer, intent(out) :: torsType
    integer :: iTors

    torsType = -1  ! Default: not found
    
    ! Check if torsion array is allocated and has elements
    if (.not. allocated(MolData(molType)%torsion)) return
    if (size(MolData(molType)%torsion) < 1) return

    do iTors = 1, size(MolData(molType)%torsion)
      if ((MolData(molType)%torsion(iTors)%mem1 == mem1) .or. &
          (MolData(molType)%torsion(iTors)%mem2 == mem1) .or. &
          (MolData(molType)%torsion(iTors)%mem3 == mem1) .or. &
          (MolData(molType)%torsion(iTors)%mem4 == mem1)) then
        if ((MolData(molType)%torsion(iTors)%mem1 == mem2) .or. &
            (MolData(molType)%torsion(iTors)%mem2 == mem2) .or. &
            (MolData(molType)%torsion(iTors)%mem3 == mem2) .or. &
            (MolData(molType)%torsion(iTors)%mem4 == mem2)) then
          if ((MolData(molType)%torsion(iTors)%mem1 == mem3) .or. &
              (MolData(molType)%torsion(iTors)%mem2 == mem3) .or. &
              (MolData(molType)%torsion(iTors)%mem3 == mem3) .or. &
              (MolData(molType)%torsion(iTors)%mem4 == mem3)) then
            if ((MolData(molType)%torsion(iTors)%mem1 == mem4) .or. &
                (MolData(molType)%torsion(iTors)%mem2 == mem4) .or. &
                (MolData(molType)%torsion(iTors)%mem3 == mem4) .or. &
                (MolData(molType)%torsion(iTors)%mem4 == mem4)) then
              torsType = MolData(molType)%torsion(iTors)%TorsType
              return
            endif
          endif
        endif
      endif
    enddo

  end subroutine
!=======================================================================
end module
!==========================================================================================
