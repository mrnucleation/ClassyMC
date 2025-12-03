!==========================================================================================
! CBMC style regrowth algorithm for Branched molecules. Much more general purpose than
! the linear version, but
!==========================================================================================
module MolCon_BranchCBMC
  use CoordinateTypes, only: Perturbation, Addition, Displacement, SingleMol, Deletion
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision

  type, public :: RosenTrial

  end type

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
    integer,private :: nAtoms 
    integer,private :: firstAtom = 1 
    integer, private :: rosenNeighList = 1
    integer,private :: nRosenTrials = 1

    !RosenProb => Weight array for chosing which trial position to use
    real(dp), private, allocatable :: GenProb(:) 
    real(dp), private, allocatable :: RosenProb(:) 
    real(dp), private, allocatable :: tempcoords(:, :) 
    real(dp), private, allocatable :: newconfig(:, :) 


    !trialcoords => Temporary storage bucket for trial atom positions
    integer, private :: nTrialAtoms
    real(dp), private, allocatable :: trialcoords(:, :, :) 
    real(dp), private, allocatable :: trialGenProb(:)

    !nGrown => Counter which tracks how many atoms are already regrown. 
    !grown => Array which specifies if an atom has been added to the system or not
    !scratchschedule => !The order in which a molecule will be grown if growing from scratch.
    integer, private :: nGrown
    logical, private, allocatable :: grown(:) 
!    integer, private, allocatable :: patharray(:) 
!    integer, private, allocatable :: pathposition(:) 
    integer, private, allocatable :: tempList(:,:), tempNNei(:)
    integer, private, allocatable :: atomtypes(:)
    real(dp), private, allocatable :: posN(:,:)

    integer, private :: nNewBonds, nNewAngles, nNewTorsions
    integer, private, allocatable :: newBonds(:)
    integer, private, allocatable :: newBondType(:)
    integer, private, allocatable :: newAngles(:)
    integer, private, allocatable :: newAngType(:)
    integer, private, allocatable :: newTorsions(:)
    integer, private, allocatable :: newTorsType(:)

    !freq => The number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds  = 0 -> Isolated Atom
    ! nBonds  = 1 -> Terminal Atom located on the end of a chain
    ! nBonds  = 2 -> Linker Atom located in the middle of a linear chain 
    ! nBonds >= 3 -> Branch Atom has two or more potential paths one can wander down.
    integer, private, allocatable :: freq(:) 
    integer, private, allocatable :: branchnext(:,:)
    contains
!      procedure, public, pass :: Constructor => BranchCBMC_Constructor
      procedure, public, pass :: Prologue => BranchCBMC_Prologue
      procedure, public, pass :: CreateSchedule => BranchCBMC_CreateSchedule
      procedure, public, pass :: GenerateConfig => BranchCBMC_GenerateConfig
      procedure, public, pass :: ReverseConfig => BranchCBMC_ReverseConfig
      procedure, public, pass :: RosenBluth => BranchCBMC_RosenBluth
      procedure, public, pass :: GetPath => BranchCBMC_GetPath
      procedure, public, pass :: GasConfig => BranchCBMC_GasConfig
      procedure, public, pass :: FindAtomsFromPath => BranchCBMC_FindAtomsFromPath
      procedure, public, pass :: ProcessIO => BranchCBMC_ProcessIO
!      procedure, public, pass :: GetNInsertPoints
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine BranchCBMC_Prologue(self)
    use Common_MolInfo, only: MolData
    use MolSearch, only: FindBond
!    use ParallelVar, only: nout
    implicit none
    class(BranchCBMC), intent(inout) :: self
!    integer, intent(in) :: molType
    integer :: iBond, iAtom
    integer :: atm1, atm2
    integer :: iError = 0
    integer :: iSchedule, nBranch
    integer :: nextAtm, prevAtm, curAtm, iPath
    integer :: leftdist, rightdist, nfirst
    integer :: largestVal

    !If there's no terminal atoms the molecule is likely cyclic
    if( all(self%freq /= 1) ) then
      iError = -1
    endif

    !Since this module is for linear molecules, if a branch is detected or
    !the atom is simply an isolated atom throw an error since it is not compatible with.
    !this module
    if(iError /= 0) then
      write(0,*) "Linear CBMC has been invoked on a molecule that is not linear. Stopping Classy."
      error stop
    endif

    if(.not. allocated(self%RosenProb)) then
      allocate(self%GenProb(1:self%nRosenTrials))
      allocate(self%RosenProb(1:self%nRosenTrials))
      allocate(self%tempcoords(1:3, 1:self%nRosenTrials))
    endif


    !Create the path array which is used to determine regrowth order.
    if(.not. allocated(self%patharray) ) then
      allocate( self%newconfig(1:3, 1:self%nAtoms) )
      allocate( self%trialcoords(1:3, 1:self%nAtoms, 1:self%nRosenTrials) )
      allocate( self%trialGenProb(1:self%nRosenTrials) )
      allocate( self%grown(1:self%nAtoms) )  
      allocate( self%patharray(1:self%nAtoms) )  
      allocate( self%pathposition(1:self%nAtoms) )  
      allocate( self%nBranches(1:self%nAtoms) )
      allocate( self%branchnext(1:self%nAtoms, 1:self%nAtoms) )  
      allocate( self%schedule(1:self%nAtoms) )  
      allocate( self%scratchschedule(1:self%nAtoms) )  
    endif

    largestVal = 0
    do iAtom = 1, self%nAtoms
      largestVal = max(largestVal, MolData(self%molType)%nAtmBonds(iAtom))
    enddo
    if(.not. allocated(self%newBonds)) allocate(self%newBonds(1:largestVal))
    if(.not. allocated(self%newBondTypes)) allocate(self%newBondTypes(1:largestVal))

    largestVal = 0
    do iAtom = 1, self%nAtoms
      largestVal = max(largestVal, MolData(self%molType)%nAtmAngles(iAtom))
    enddo
    if(.not. allocated(self%newAngles)) allocate(self%newAngles(1:largestVal))
    if(.not. allocated(self%newAngleTypes)) allocate(self%newAngleTypes(1:largestVal))

    largestVal = 0
    do iAtom = 1, self%nAtoms
      largestVal = max(largestVal, MolData(self%molType)%nAtmTorsions(iAtom))
    enddo
    if(.not. allocated(self%newTorsions)) allocate(self%newTorsions(1:largestVal))
    if(.not. allocated(self%newTorsionTypes)) allocate(self%newTorsionTypes(1:largestVal))

    !Count the number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds  = 0 >> Isolated Atom
    ! nBonds  = 1 >> Terminal Atom located on the end of a chain
    ! nBonds  = 2 >> Linker Atom. Located in the middle of a linear chain
    ! nBonds >= 3 >> Branch Atom. Has two or more potential paths one can wander down.
    self%nAtoms = MolData(self%molType)%nAtoms
    self%include15 = .false.
    self%nBranches = 0
    self%branchnext = 0
    if(self%nAtoms > 1) then
      allocate( self%freq(1:self%nAtoms) )  
      self%freq = 0
      do iBond = 1, size(MolData(self%molType)%bond)
        atm1 = MolData(self%molType)%bond(iBond)%mem1
        self%freq(atm1) = self%freq(atm1) + 1
        atm2 = MolData(self%molType)%bond(iBond)%mem2
        self%freq(atm2) = self%freq(atm2) + 1
         !Branchnext is basically a neighborlist
         !keeping track of which atoms are directly connected.
         !Used to determine which atoms to grow next
        nBranch = self%nBranches(atm1) + 1
        self%nBranches(atm1) = nBranch
        self%branchnext(nBranch, atm1) = atm2

        nBranch = self%nBranches(atm2) + 1
        self%nBranches(atm2) = nBranch
        self%branchnext(nBranch, atm2) = atm1
      enddo
    else 
      iError = -1
    endif


  end subroutine
!==========================================================================================
  subroutine BranchCBMC_GenerateConfig(self, trialBox, disp, probconstruct, accept, insPoint, insProb)
    use Common_MolInfo, only: BondData, AngleData, TorsionData
    use MolSearch, only: FindBond, FindAngle, FindTorsion
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone, Generate_UnitTorsion, ListRNG
    use ForcefieldData, only: ECalcArray
    use SearchSort, only: QSort

    implicit none
    class(BranchCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:,:)
    real(dp), intent(in), optional :: insProb(:)
    real(dp), intent(out) :: probconstruct 
    logical, intent(out) :: accept

    integer :: iNei
    integer :: nNext
    integer :: nextlist(1:self%nAtoms)
    integer :: nSchedule
    integer :: schedule(1:self%nAtoms)
    integer :: slice(1:2)
    integer :: curAtm
    real(dp), pointer :: atoms(:,:) => null()


    probconstruct = 1E0_dp
    nNext = 0
    select type(disp)
      class is(Displacement)
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
        slice(1) = molStart
        slice(2) = molEnd
        call trialbox%GetCoordinates(atoms, slice=slice)
        self%newconfig(1:3, 1:self%nAtoms) = atoms(1:3, 1:self%nAtoms)

      class is(Addition)
        self%grown = .false.
        self%nGrown = 0
        molindx = disp(1)%molindx
        call trialbox%GetMolData(molindx, molStart=molStart)
        do iDisp = 1, size(disp)
          atm1 = disp(iDisp)%atmindx - molStart + 1
          dispsubindx(iDisp) = atm1
          atmdispindx(atm1) = iDisp
        enddo

      class default
        write(0,*) "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
        error stop 
    end select

    if(present(insPoint)) then
      if(size(insPoint,2) /= self%nRosenTrials) then
        write(0,*) "ERROR! Branch CBMC Regrowth received a different number of insertion points"
        write(0,*) "than it was expecting!"
        error stop
      endif


    endif

    !We
    do while(nNext > 0)
      curAtm = nextlist(nNext)
      nSchedule = 0
       !Figure out which neighboring atoms have already been
       !grown and which ones haven't. 
      do iNei = 1, self%nBranches(curAtm)
        if(.not. self%grown(self%branchnext(iNei, curAtm)) ) then
          nSchedule = nSchedule + 1
          schedule(nSchedule) = self%branchnext(iNei, curAtm)
        endif
      enddo

      call QSort( schedule(1:nSchedule) )
      call self%CreateTrials( schedule(1:nSchedule) ) 
      call self%SelectTrial()

      nNext = nNext - 1
    enddo


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

  end subroutine
!======================================================================
!==========================================================================================
  subroutine BranchCBMC_GasConfig(self, probGas)

!  Routine for simulating the probability of an isolated molecule in the gas phase
!  Used primarily for swap moves with an implict gas box.
    implicit none
    class(BranchCBMC), intent(inout) :: self
    real(dp), intent(out) :: probGas

    probGas = 1E0_dp
    if(.not. self%include15) then
      probGas = 1E0_dp/real(self%nRosenTrials, dp)**(self%nAtoms)
    endif
!    write(*,*) probGas, self%nAtoms, self%nRosenTrials
    

  end subroutine
!=======================================================================
  subroutine BranchCBMC_CreateTrials(self, curAtm, nTrials, schedule, trialbox)
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData
    use SearchSort, only: BinarySearch
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(in) :: curAtm
    integer, intent(in) :: nTrials
    real(dp), intent(in) :: schedule(:)
    class(SimBox), intent(inout) :: trialBox

    integer :: iTrial
    integer :: atm1, atm2, atm3, atm4
    integer :: iBond, iAngle, iTorsion
    integer :: nBonds, nAngles, nTorsions
    integer :: curbond, curangle, curtorsion
    integer :: originBond = -1
    real(dp) :: originAtom(1:3)
    real(dp) :: v1(1:3)

    nBonds = MolData(self%molType)%nAtmBonds(curAtm)
    nAngles = MolData(self%molType)%nAtmAngles(curAtm)
    nTorsions = MolData(self%molType)%nAtmTorsions(curAtm)

    self%newBonds = 0
    nNewBonds = 0
    do iBond = 1, nBonds
      curbond = MolData(self%molType)%atmBonds(iBond, curAtm)
      atm1 = MolData(self%molType)%bond(curbond)%mem1
      atm2 = MolData(self%molType)%bond(curbond)%mem2
      if( 
      if( self%grown(atm1) .and.  self%grown(atm2) ) then
        if(originBond > 0) stop "ERROR! BranchCBMC has been used for a molecule that's potentially cyclic"
        originBond = curbond
      else
        nNewBonds = nNewBonds + 1
        self%newBonds(nNewBonds) = curbond
      endif
    enddo

    self%newAngles = 0
    nNewAngles = 0
    do iAngle = 1, nAngles
      curangle = MolData(self%molType)%atmAngles(iAngle, curAtm)
      atm1 = MolData(self%molType)%angle(curangle)%mem1
      atm2 = MolData(self%molType)%angle(curangle)%mem2
      atm3 = MolData(self%molType)%angle(curangle)%mem3
      if( .not. (self%grown(atm1) .and.  self%grown(atm2) .and. self%grown(atm3) ) ) then
        nNewAngles = nNewAngles + 1
        self%newAngles(nNewAngles) = curangle
      endif
    enddo

    self%newTorsions = 0
    nNewTorsions = 0
    do iTorsion = 1, nTorsions
      curtorsion = MolData(self%molType)%atmTorsions(iTorsion, curAtm)
      atm1 = MolData(self%molType)%torsion(curtorsion)%mem1
      atm2 = MolData(self%molType)%torsion(curtorsion)%mem2
      atm3 = MolData(self%molType)%torsion(curtorsion)%mem3
      atm4 = MolData(self%molType)%torsion(curtorsion)%mem4
      if( .not. (self%grown(atm1) .and.  self%grown(atm2) .and. self%grown(atm3).and. self%grown(atm4) ) ) then
        nNewTorsions = nNewTorsions + 1
        self%newTorsions(nNewTorsions) = curtorsion
      endif
    enddo

     !If we have an origin bond (IE atom was generated by previously growing a bond)
     !We use that as our reference vector
    if(originBond > 0) then
      atm1 = MolData(self%molType)%bond(originBond)%mem1
      atm2 = MolData(self%molType)%bond(originBond)%mem2
      if(atm1 == curAtm) then
        v1(1:3) = self%newConfig(1:3, atm2) - self%newConfig(1:3, atm1)
        originAtom(1:3) = self%newConfig(1:3, atm1)
      else
        v1(1:3) = self%newConfig(1:3, atm1) - self%newConfig(1:3, atm2)
        originAtom(1:3) = self%newConfig(1:3, atm2)
      endif
      call trialBox%Boundary(v1(1), v1(2), v1(3))
    else
      v1 = 0E0_dp
    endif

    do iTrial = 1, nTrials
       !If there's no origin we generate a random orientation to start from.

      select case(nNewBonds)
        case(1) !Linear Single Branch
          call self%OneBranch(v1, curatm, iTrial, originAtom, originBond, trialBox, probGen)
        case(2) !Two Branch
          call self%TwoBranch
        case(3) !Three Branch
          call self%ThreeBranch
!        case(4) !Four Branch
        case default
          write(0,*) nNewBonds
          stop "ERROR! Branch CBMC has not been set up for this number of branches yet!"
      end select
    enddo

  end subroutine 
!=======================================================================
  subroutine BranchCBMC_OneBranch(self, v1, curatm, iTrial, originAtom, originBond, trialBox, probGen)
     !To create a new single branch we need to generate three variables
     !The bond distance (r), the bond angle(theta), and torsion angle(phi)
     !If there's no new bond angles or torsion angles, theta and phi will simply
     !be assigned uniformly.
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData
    use Math_Coordinates, only: GetPerpendicularVector
    implicit none
    class(BranchCBMC), intent(inout) :: self
    real(dp), intent(in) :: v1(1:3)
    integer, intent(in) :: curatm, originBond, iTrial
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in) :: originAtom(1:3)
    real(dp), intent(out) :: probGen
    logical :: propertors(1:self%nNewTorsions)
    integer :: bondIndx, bondType
    integer :: angleIndx, angleType
    integer :: torsionIndx, torsionType
    integer :: iTors
    integer :: atm1, atm2, atm3, atm4
    integer :: tors
    real(dp) :: atompos(1:3, 1:4)
    real(dp) :: r, theta, phi
    real(dp) :: prob_r, prob_theta, prob_phi
    real(dp) :: prevphi(1:self%nNewTorsions)
    real(dp) :: phi_tors
    real(dp) :: v2(1:3)
    real(dp) :: vpervec(1:3)
    real(dp) :: E_Tors

    bondIndx = self%newBonds(1)
    bondType = MolData(self%molType)%bond(bondIndx)%bondType
    call BondData(bondType) % bondFF % GenerateDist(trialBox%beta, r, prob_r)



    prob_theta = 1E0_dp
    prob_phi = 1E0_dp
    if( (self%nNewAngles < 1) .and. (self%nNewTorsions < 1) ) then
      call Generate_UnitSphere(v2(1), v2(2), v2(3))
    else if(self%nNewTorsions < 1) then
       ! For a single branch in a non-cyclic molecule, there can only be one
       ! bending angle for this branch therefore 
      angleIndx = self%newAngle(1)
      angleType = MolData(self%molType)%angle(angleIndx)%angleType
      call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, theta, prob_theta)

       ! For the torsional component however, we can easily have multiple torsion angles
       ! thus we need to generate according to that. We may have several
       ! rejections before a good angle is found so for speed, we need to get a common reference point
       ! so we don't have to recompute the torsion energies from the raw coordinates every single time.
       ! We only include proper torsion in this stage.
      prevphi = 0E0_dp
      propertors = .true.
      call GetPerpendicularVector(v1, vpervec, phi=0E0_dp)
      do iTors = 1, self%nNewTorsion
        torsionIndx = self%newTorsion(iTors)
        torsionType = self%newTorsionType(iTors)
        atm1 = MolData(self%molType)%torsion(torsionIndx)%mem1
        atm4 = MolData(self%molType)%torsion(torsionIndx)%mem4
        if(curatm == atm1) then
          atm1 = MolData(self%molType)%torsion(torsionIndx)%mem4
          atm2 = MolData(self%molType)%torsion(torsionIndx)%mem3
          atm3 = MolData(self%molType)%torsion(torsionIndx)%mem2
          atm4 = MolData(self%molType)%torsion(torsionIndx)%mem1
        else if(curatm == atm4)
          atm1 = MolData(self%molType)%torsion(torsionIndx)%mem1
          atm2 = MolData(self%molType)%torsion(torsionIndx)%mem2
          atm3 = MolData(self%molType)%torsion(torsionIndx)%mem3
          atm4 = MolData(self%molType)%torsion(torsionIndx)%mem4
        else
          !If the newly grown atom is
          properTors(iTors) = .false.
          cycle
        endif
        atompos(1:3, 1) = self%newconfig(1:3, atm1)
        atompos(1:3, 2) = self%newconfig(1:3, atm2)
        atompos(1:3, 3) = self%newconfig(1:3, atm3)
        atompos(1:3, 4) = self%newconfig(1:3, atm3) + vpervec(1:3) 
        prevphi(iTors) = TorsionData(torsionType) % torsionFF % ComputeTors(trialBox, atompos)
      enddo

      !Now that we have the reference points, we can simply compute all torsion angles
      !by adding the offset from the reference angle. We generate an angle
      !that
      do
        E_Tors = 0E0_dp
        phi = two_pi*grnd()
        do iTors = 1, self%nNewTorsion
          if( .not. properTors(iTors) ) cycle
          torsionIndx = self%newTorsion(iTors)
          torsionType = self%newTorsionType(iTors)
          phi_tors = prevphi(iTors) + phi
          E_Tors = E_Tors + TorsionData(torsionType) % torsionFF % EFunc(phi_tors)
        enddo
        prob_phi = exp(-trialBox%beta*E_Tors)
        if( prob_phi > grnd() ) exit
      enddo
      call Generate_RelativeVector(v1, r2, theta, phi, v2)

    endif
     !Now that we have the vectors relative to the previous Atom, we get the
     !absolute position of the new atom and store the position in the trial
     !array
    v2(1:3) = originAtom(1:3) + v2(1:3)
    call trialBox%Boundary( v2(1), v2(2), v2(3) )
    self%trialcoords(1:3, curatm, iTrial) = v2(1:3)
    probGen = prob_r * prob_theta * prob_phi
    self%trialProbGen(iTrial) = probGen

  end subroutine 
!=======================================================================
  subroutine BranchCBMC_RosenBluth(self, trialBox, atmsubindx, iDisp, disp)
    ! Computes the Rosenbluth Weight of each trial position.  This is used
    ! to select which trial position should be used to regrow.
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

    logical :: accept
    logical :: overlap(1:self%nRosenTrials)
    integer :: molIndx, molType, atmIndx, molStart
    integer :: iRosen, jAtom, jNei, neiSize
    integer :: atmtype1, atmNeiIndx
    real(dp) :: E_Atom, E_Min, norm
    real(dp) :: pos1(1:3)
    type(Displacement) :: tempdisp(1:1)


    if(self%nRosenTrials == 1) then
      self%RosenProb(1) = 1E0_dp
      norm = 1E0_dp
      return
    endif

!    self%RosenProb(1:self%nRosenTrials) = 1E0_dp
!    norm = sum(self%RosenProb(1:self%nRosenTrials))
!    return

    if( (iDisp < 1) .or. (iDisp > size(disp)) ) then
      write(0,*) "Invalid Perturbation Index!"
      write(0,*) iDisp, size(disp)
      error stop
    endif

    select type(trialbox)
      class is(SimpleBox)
        call trialbox%GetCoordinates(atoms)
        call trialbox%GetEFunc(EFunc)
        call trialBox%GetNeighborList(self%rosenNeighList, neighlist, nNeigh)
    end select

    if(.not. allocated(self%tempList)) then
      neiSize = size(neighlist, 1)
      allocate( self%tempList(1:neiSize, 1:1) )
      allocate( self%tempNNei(1:neiSize)      )
      allocate( self%atomtypes(1:neiSize)     )
      allocate( self%posN(1:3, 1:neiSize)     )
    endif

    select type(disp)
      class is(SingleMol)
        tempdisp(1)%molType = disp(iDisp)%molType
        tempdisp(1)%molindx = disp(iDisp)%molindx
!        tempdisp(1)%atmindx = disp(iDisp)%atmindx
    end select

    call trialBox%GetMolData(tempdisp(1)%molindx, molStart=molStart)

    overlap = .false.

    atmIndx = atmSubIndx + molStart - 1
    tempdisp(1)%atmindx = atmIndx
    atmNeiIndx = atmIndx
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
            call trialbox%GetNewNeighborList(self%rosenNeighList, 1, self%tempList, self%tempNNei, tempdisp(1))
          end select
          neighlist => self%tempList
          nNeigh => self%tempNNei 
          atmNeiIndx = 1
      end select
      pos1(1:3) = self%tempcoords(1:3, iRosen)
      accept = .true.
      if(nNeigh(atmNeiIndx) > 0) then
        do jNei = 1, nNeigh(atmNeiIndx)
          jAtom = neighlist(jNei, atmNeiIndx)
          call trialBox%GetAtomData(jAtom, atomtype=self%atomtypes(jNei))
          self%posN(1:3, jNei) = atoms(1:3, jAtom)
        enddo
        !Add 1-5 terms later.
        call EFunc % Method % ManyBody(trialbox,& 
                           atmtype1,& 
                           pos1,& 
                           self%atomtypes(1:nNeigh(atmNeiIndx)),& 
                           self%posN(1:3,1:nNeigh(atmNeiIndx)),&
                           E_Atom,&
                           accept)
      else
        E_Atom = 0E0_dp
      endif
      if(accept) then
          if(E_Atom < E_Min) then
            E_Min = E_Atom
          endif
          self%RosenProb(iRosen) = E_Atom
      else
          overlap(iRosen) = .true.
          self%RosenProb(iRosen) = 0E0_dp
      endif
    enddo

    !We now compute the weight of the Rosenbluth trial.  P(E_i) = exp(-E_i/kt)/N. 
    !All trials are normalized by E_Max to prevent floating point overflow, but the
    !extra term cancels out in probability leaving the result unchanged. 
    norm = 0E0_dp
    do iRosen = 1, self%nRosenTrials
      if( overlap(iRosen) ) cycle

      self%RosenProb(iRosen) = (self%RosenProb(iRosen)-E_Min)*trialbox%beta
      self%RosenProb(iRosen) = exp(-self%RosenProb(iRosen))
      if(IsNan(self%RosenProb(iRosen)) .or. IsInf(self%RosenProb(iRosen))) then
        write(0,*) "Invalid Weight Problem found in CBMC module!"
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
    read(command,*) nRosen
    self%inspoints = nRosen
    self%nRosenTrials = nRosen

  end subroutine
!==========================================================================================
  function BranchCBMC_InSchedule(self, schedule, atmcur, atm1) result(
    implicit none
    class(BranchCBMC), intent(in) :: self
    integer, intent(in) :: schedule(:)
    logical :: 

    


  end function
!=======================================================================
end module
!==========================================================================================
