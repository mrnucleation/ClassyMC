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

    !GenProb => Weight array for chosing which trial position to use
    !RosenProb => Weight array for chosing which trial position to use
    !tempcoords => Temporary storage bucket for trial atom positions
    real(dp), private, allocatable :: GenProb(:) 
    real(dp), private, allocatable :: RosenProb(:) 
    real(dp), private, allocatable :: tempcoords(:, :) 
    real(dp), private, allocatable :: newconfig(:, :) 

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
    integer, private, allocatable :: newBonds(:)
    integer, private, allocatable :: newAngles(:)
    integer, private, allocatable :: newTorsions(:)

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

    largestVal = 0
    do iAtom = 1, self%nAtoms
      largestVal = max(largestVal, MolData(self%molType)%nAtmAngles(iAtom))
    enddo
    if(.not. allocated(self%newAngles)) allocate(self%newAngles(1:largestVal))

    largestVal = 0
    do iAtom = 1, self%nAtoms
      largestVal = max(largestVal, MolData(self%molType)%nAtmTorsions(iAtom))
    enddo
    if(.not. allocated(self%newTorsions)) allocate(self%newTorsions(1:largestVal))

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

      call self%CreateTrials( schedule(1:nSchedule) )

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


    integer :: dispsubindx(1:self%nAtoms), atmdispindx(1:self%nAtoms)
    integer :: bondType, angleType, torsType
    integer :: molindx, molStart, molEnd
    integer :: atm1, atm2,atm3,atm4, iDisp, iRosen
    integer :: lastGrown
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r
    real(dp) :: norm
    real(dp) :: bend_angle,tors_angle
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
        error stop "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
    end select

    if(present(insPoint) .neqv. present(insProb)) then
      write(0,*) "insPoint and insProb must both be present or abscent!"
      error stop
    endif

    if(present(insPoint)) then
      if(size(insPoint,2) /= self%nRosenTrials) then
        write(0,*) "ERROR! Linear CBMC Regrowth received a different number of insertion points"
        write(0,*) "than it was expecting!"
        error stop
      endif
    endif


    do while(self%nGrown < self%nAtoms)
      if(self%nGrown == 0) then
          Atm1 = self%schedule(1)
          if(present(inspoint) .and. present(insProb)) then
            do iRosen = 1, self%nRosenTrials-1
              self%tempcoords(1:3, iRosen) = insPoint(1:3,iRosen)
              self%GenProb(iRosen) = insProb(iRosen)
            enddo
            self%tempcoords(1:3, self%nRosenTrials) = self%newconfig(1:3, Atm1)
            self%GenProb(self%nRosenTrials) = insProb(self%nRosenTrials)
          else
            error stop "Full Regrowth has been requsted without passing in insertion points"
          endif 
          lastGrown = Atm1

      elseif(self%nGrown == 1) then
        Atm1 = self%schedule(1)
        Atm2 = self%schedule(2)
        if(.not. self%grown(Atm1) ) then
          write(0,*)  Atm1, Atm2
          write(0,*) self%grown(1:self%nAtoms)
          error stop "Unexpected gap in Linear CBMC!"
        endif
        call FindBond(self%molType, Atm1, Atm2, bondType)
        do iRosen = 1, self%nRosenTrials-1
          call BondData(bondType) % bondFF % GenerateDist(trialBox%beta, r, prob_r)
          probgen = prob_r
          call Generate_UnitSphere(dx, dy, dz)
          self%tempcoords(1, iRosen) = r*dx + self%newconfig(1, Atm1)
          self%tempcoords(2, iRosen) = r*dy + self%newconfig(2, Atm1)
          self%tempcoords(3, iRosen) = r*dz + self%newconfig(3, Atm1)
          call trialBox%Boundary(self%tempcoords(1, iRosen),&
                                 self%tempcoords(2, iRosen),&
                                 self%tempcoords(3, iRosen))
          self%GenProb(iRosen) = prob_r
        enddo
        !Compute the generation probability of the old position

        oldpos(1:3, 1) = atoms(1:3, Atm1)
        oldpos(1:3, 2) = atoms(1:3, Atm2)   

        call BondData(bondType) % bondFF % GenerateReverseDist(trialbox, &
                                                     oldpos(1:3,1:2), &
                                                     prob_r)
!        write(*,*) "2", prob_r
        self%GenProb(self%nRosenTrials) = prob_r
        self%tempcoords(1:3, self%nRosenTrials) = self%newconfig(1:3, Atm2)
        lastGrown = Atm2

      elseif(self%nGrown == 2) then
        Atm3 = self%schedule(3) 
        Atm2 = self%schedule(2) 
        if(self%freq(Atm2) == 2) then
            Atm1 = self%schedule(1) 
            Atm2 = self%schedule(2) 
        else
            Atm1 = self%schedule(2) 
            Atm2 = self%schedule(1) 
        endif

        v1(1:3) = self%newconfig(1:3, Atm1) - self%newconfig(1:3, Atm2)
        call trialBox%Boundary(v1(1), v1(2), v1(3))

        call FindBond(self%molType, Atm2, Atm3, bondType)
        call FindAngle(self%molType, Atm1, Atm2, Atm3, angleType)
        do iRosen = 1, self%nRosenTrials-1
          call BondData(bondType) % bondFF % GenerateDist(trialBox%beta,r, prob_r)
          call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, bend_angle, prob_ang)
          probgen = prob_r * prob_ang
          call Generate_UnitCone(v1, r, bend_angle, v2)
          self%tempcoords(1:3, iRosen) = v2(1:3) + self%newconfig(1:3, Atm2)
          call trialBox%Boundary(self%tempcoords(1, iRosen),&
                                 self%tempcoords(2, iRosen),&
                                 self%tempcoords(3, iRosen))

          self%GenProb(iRosen) = probgen
        enddo
        oldpos(1:3, 1) = atoms(1:3, Atm1)
        oldpos(1:3, 2) = atoms(1:3, Atm2)
        oldpos(1:3, 3) = atoms(1:3, Atm3)
        call BondData(bondType) % bondFF % GenerateReverseDist(trialbox, &
                                                     oldpos(1:3,2:3), &
                                                     prob_r)


        call AngleData(angleType) % angleFF % GenerateReverseDist(trialbox, &
                                                     oldpos(1:3,1:3), &
                                                     prob_ang)

!        write(*,*) "3", prob_r, prob_ang
        probgen = prob_r * prob_ang
        self%GenProb(self%nRosenTrials) = probgen
        self%tempcoords(1:3, self%nRosenTrials) = self%newconfig(1:3, Atm3)
        lastGrown = Atm3

      elseif(self%nGrown < self%nAtoms) then
            Atm4 = self%schedule(self%nGrown+1) 
            call self%FindAtomsFromPath(Atm4, Atm1, Atm2, Atm3)

            call FindBond(self%molType, Atm3, Atm4, bondType)
            call FindAngle(self%molType, Atm2, Atm3, Atm4, angleType)
            call FindTorsion(self%molType, Atm1, Atm2, Atm3, Atm4, torsType)

            v1(1:3) = self%newconfig(1:3, Atm1) - self%newconfig(1:3, Atm3)
            call trialBox%Boundary(v1(1), v1(2), v1(3))
            v2(1:3) = self%newconfig(1:3, Atm2) - self%newconfig(1:3, Atm3)
            call trialBox%Boundary(v2(1), v2(2), v2(3))

            do iRosen = 1, self%nRosenTrials-1
              call BondData(bondType) % bondFF % GenerateDist(trialBox%beta, r, prob_r)
              call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, bend_angle, prob_ang)
              call TorsionData(torsType) % torsionFF % GenerateDist(trialBox%beta, tors_angle, prob_tors)
              probgen = prob_r*prob_ang*prob_tors
              call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
              self%tempcoords(1:3, iRosen)  = v3(1:3) + self%newconfig(1:3, Atm3)
              call trialBox%Boundary(self%tempcoords(1, iRosen),&
                                     self%tempcoords(2, iRosen),&
                                     self%tempcoords(3, iRosen))

              self%GenProb(iRosen) = probgen
            enddo
            oldpos(1:3, 1) = atoms(1:3, Atm1)
            oldpos(1:3, 2) = atoms(1:3, Atm2)
            oldpos(1:3, 3) = atoms(1:3, Atm3)
            oldpos(1:3, 4) = atoms(1:3, Atm4)

            call BondData(bondType) % bondFF % GenerateReverseDist(trialbox, &
                                                     oldpos(1:3,3:4), &
                                                     prob_r)


            call AngleData(angleType) % angleFF % GenerateReverseDist(trialbox, &
                                                     oldpos(1:3,2:4), &
                                                     prob_ang)

            call TorsionData(torsType) % torsionFF % GenerateReverseDist(trialbox, &
                                                     oldpos(1:3,1:4), &
                                                     prob_tors)
            probgen = prob_r * prob_ang * prob_tors
            self%GenProb(self%nRosenTrials) = probgen
            self%tempcoords(1:3, self%nRosenTrials) = self%newconfig(1:3, Atm4)
            lastGrown = Atm4
        else
            error stop "nGrown is some invalid number."
        endif


        iDisp = atmdispindx(lastGrown)
        call self%RosenBluth(trialBox, lastGrown, iDisp, disp)
        norm = 0E0_dp
        !For this method we leave off the prob term since we can generate the angles and such
        !directly. 
        do iRosen = 1, self%nRosenTrials
          norm = norm + self%RosenProb(iRosen)
        enddo
        if(norm < 1E-10_dp) then
          accept = .false.
          return
        endif
        probconstruct = probconstruct * self%GenProb(self%nRosenTrials) * self%RosenProb(self%nRosenTrials)/norm

!        write(*,*) self%nRosenTrials, self%GenProb(self%nRosenTrials), self%RosenProb(1:self%nRosenTrials)/norm
        self%grown(lastGrown) = .true.
        self%nGrown = self%nGrown + 1
     enddo

  end subroutine
!======================================================================
!  Routine for simulating the probability of an isolated molecule in the gas phase
!  Used primarily for swap moves with an implict gas box.
!==========================================================================================
  subroutine BranchCBMC_GasConfig(self, probGas)
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
  subroutine BranchCBMC_CreateTrials(self, curAtm, nTrials, schedule)
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData
    implicit none
    class(BranchCBMC), intent(inout) :: self
    integer, intent(in) :: curAtm
    integer, intent(in) :: nTrials
    real(dp), intent(in) :: schedule(:)

    integer :: iTrial
    integer :: atm1, atm2, atm3, atm4
    integer :: iBond, iAngle, iTorsion
    integer :: nBonds, nAngles, nTorsions
    integer :: nNewBonds, nNewAngles, nNewTorsions
    integer :: curbond, curangle, curtorsion
    integer :: originBond = -1
    real(dp) :: v1(1:3)

    nBonds = MolData(self%molType)%nAtmBonds(curAtm)
    nAngles = MolData(self%molType)%nAtmAngles(curAtm)
    nTorsions = MolData(self%molType)%nAtmTorsions(curAtm)
    nNewBonds = 0
    nNewAngles = 0
    nNewTorsions = 0

    self%newBonds = 0
    do iBond = 1, nBonds
      curbond = MolData(self%molType)%atmBonds(iBond, curAtm)
      atm1 = MolData(self%molType)%bond(curbond)%mem1
      atm2 = MolData(self%molType)%bond(curbond)%mem2
      if( self%grown(atm1) .and.  self%grown(atm2) ) then
        if(originBond > 0) stop "ERROR! Branch has been used for a molecule that's potentially cyclic"
        originBond = curbond
      else
        nNewBonds = nNewBonds + 1
        self%newBonds(nNewBonds) = curbond
      endif
    enddo

    self%newAngles = 0
    do iAngle = 1, nAngles
      curangle = MolData(self%molType)%atmAngles(iAngle, curAtm)
      atm1 = MolData(self%molType)%angle(curangle)%mem1
      atm2 = MolData(self%molType)%angle(curangle)%mem2
      atm3 = MolData(self%molType)%angle(curangle)%mem3
      if( .not. (self%grown(atm1) .and.  self%grown(atm2) .and. self%grown(atm3) ) then
        nNewAngles = nNewAngles + 1
        self%newAngles(nNewAngles) = curangle
      endif
    enddo

    self%newTorsions = 0
    do iTorsion = 1, nTorsions
      curtorsion = MolData(self%molType)%atmTorsions(iTorsion, curAtm)
      atm1 = MolData(self%molType)%torsion(curtorsion)%mem1
      atm2 = MolData(self%molType)%torsion(curtorsion)%mem2
      atm3 = MolData(self%molType)%torsion(curtorsion)%mem3
      atm4 = MolData(self%molType)%torsion(curtorsion)%mem4
      if( .not. (self%grown(atm1) .and.  self%grown(atm2) .and. self%grown(atm3).and. self%grown(atm4) ) then
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
      else
        v1(1:3) = self%newConfig(1:3, atm1) - self%newConfig(1:3, atm2)
      endif
    endif

    do iTrial = 1, nTrials
       !If there's no origin we generate a random orientation to start from.
      if(originBond < 1) then
        call Generate_UnitSphere(v1(1),v1(2),v1(3))
      endif
      select case(nNewBonds)
        case(1) !Linear Single Branch
        case(2) !Two Branch
        case(3) !Three Branch
        case(4) !Four Branch
        case default
          write(0,*) nNewBonds
          stop "ERROR! Branch CBMC has not been set up for this number of branches yet!"
      end select
    enddo

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
!=======================================================================
end module
!==========================================================================================
