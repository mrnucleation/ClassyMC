!==========================================================================================
! Simple Regrowth Object
!==========================================================================================
module MolCon_LinearCBMC
  use CoordinateTypes, only: Perturbation, Addition
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: LinearCBMC
!    integer :: inspoints = 1 !
    integer :: firstAtom = 1 !Index of the first atom to be inserted into the system
    integer :: nRosenTrials = 8 !Number of Rosenbluth trials to use.  This is the number of potential positions
                                !that are generated for each new atom as it is being regrown.  One of these positions
                                !will be selected for the final atom position.  

    real(dp), allocatable :: GenProb(:) !Weight array for chosing which trial position to use
    real(dp), allocatable :: RosenProb(:) !Weight array for chosing which trial position to use
    real(dp), allocatable :: tempcoords(:, :) !Temporary storage bucket for trial atom positions

    logical, allocatable :: grown(:) !Array which specifies if an atom has been added to the system or not
    integer, allocatable :: patharray(:) !Array which lays out the topology of the molecule 
    integer, allocatable :: pathposition(:) !Array which takes an atom's ID and gives the position in the patharray
    integer, allocatable :: schedule(:) !The order in which a molecule will be grown

    integer, allocatable :: freq(:) !The number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds  = 0 >> Isolated Atom
    ! nBonds  = 1 >> Terminal Atom located on the end of a chain
    ! nBonds  = 2 >> Linker Atom located in the middle of a linear chain
    ! nBonds >= 3 >> Branch Atom has two or more potential paths one can wander down.

    contains
!      procedure, public, pass :: Constructor => LinearCBMC_Constructor
      procedure, public, pass :: Prologue => LinearCBMC_Prologue
      procedure, public, pass :: GenerateConfig => LinearCBMC_GenerateConfig
      procedure, public, pass :: ReverseConfig => LinearCBMC_ReverseConfig
      procedure, public, pass :: FindAtomsFromPath => LinearCBMC_FindAtomsFromPath
!      procedure, public, pass :: GetNInsertPoints
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine LinearCBMC_Prologue(self)
    use Common_MolInfo, only: MolData, BondData, nMolTypes
    use MolSearch, only: FindBond
    use ParallelVar, only: nout
    implicit none
    class(LinearCBMC), intent(inout) :: self
!    integer, intent(in) :: molType
    integer :: iType, iBond, iAtom, curMax, nAtoms
    integer :: atm1, atm2
    integer :: iError = 0
    integer :: iSchedule
    integer :: nextAtm, prevAtm, curAtm, iPath
    integer ::  leftdist, rightdist, nleft, nright, nfirst

    !Count the number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds = 0 >> Isolated Atom
    ! nBonds = 1 >> Terminal Atom located on the end of a chain
    ! nBonds = 2 >> Linker Atom.  Located in the middle of a linear chain
    ! nBonds >= 3 >> Branch Atom.  Has two or more potential paths one can wander down.
    nAtoms = MolData(self%molType)%nAtoms
    if(MolData(self%molType)%nAtoms > 1) then
      allocate( self%freq(1:nAtoms) )  
      self%freq = 0
      do iBond = 1, size(MolData(self%molType)%bond)
        atm1 = MolData(self%molType)%bond(iBond)%mem1
        self%freq(atm1) = self%freq(atm1) + 1
        atm2 = MolData(self%molType)%bond(iBond)%mem2
        self%freq(atm2) = self%freq(atm2) + 1
      enddo
    else 
      iError = -1
    endif


    !If any atom has more than 3 bonds it's a branched molecule and not suited for this algorithm.
    if( any(self%freq > 2) ) then
      iError = -1
    endif

    !If there's no terminal atoms the molecule is likely cyclic
    if( all(self%freq /= 1) ) then
      iError = -1
    endif


    !Since this module is for linear molecules, if a branch is detected or
    !the atom is simply an isolated atom throw an error since it is not compatible with.
    !this module
    if(iError /= 0) then
      write(0,*) "Linear CBMC has been invoked on a molecule that is not linear. Stopping Classy."
      stop
    endif

    if(.not. allocated(self%RosenProb)) then
      allocate(self%RosenProb(1:self%nRosenTrials))
      allocate(self%tempcoords(1:3, 1:self%nRosenTrials))
    endif


    !Create the path array which is used to determine regrowth order.
    if(.not. allocated(self%patharray) ) then
      allocate( self%grown(1:nAtoms) )  
      allocate( self%patharray(1:nAtoms) )  
      allocate( self%pathposition(1:nAtoms) )  
      allocate( self%schedule(1:nAtoms) )  
    endif
    self%patharray = 0
    self%pathposition = 0
    !Pick an terminal atom to start building the patharray
    prevAtm = 0
    do iAtom = 1,nAtoms
      if(freq(iAtom) == 1) then
        curatm = iAtom
        exit
      endif
    enddo

    self%patharray(1) = curatm
    self%pathposition(curatm) = 1
    iPath = 1
    ! Begin the building process by traversing from bond to bond.
    do iAtom = 1, nAtoms-1
      nextatm = 0
      do iBond = 1, size(MolData(self%molType)%bond)
        atm1 = MolData(self%molType)%bond(iBond)%mem1
        atm2 = MolData(self%molType)%bond(iBond)%mem2
        if(atm1 == curatm) then
          if(atm2 /= prevatm) then
            nextatm = atm2
            exit
          endif
        endif
        if(atm2 == curatm) then
          if(atm1 /= prevatm) then
            nextatm = atm1
            exit
          endif
        endif
      enddo
!      write(*,*) prevatm, curatm, nextatm
      if(nextatm /= 0) then
         iPath = iPath + 1
         self%patharray(iPath) = nextatm
         self%pathposition(nextatm) = iPath
         prevatm = curatm
         curatm = nextatm
      else
        write(0,*) "Linear CBMC: DEAD END ERROR! There is a problem in the molecule definition!"
        write(0,*) "Ensure your molecule's bonds are properly connected."
        error stop
      endif
    enddo


    if(any(self%patharray == 0)) then
      write(0,*) "LINEAR CBMC: ERROR! There are atoms which are unaccounted for by the path building algorithm."
      write(0,*) "Ensure your molecule's bonds are properly connected."
      error stop "User input error!"
    endif


!    write(*,*) self%patharray
!    write(*,*) self%pathposition
    ! We now construct the scheduling array which is responsible for determing the order
    ! that the atoms are inserted into the system.  To do this first we determine
    ! after the first atom is inserted into the system
    nfirst = self%pathposition(self%firstAtom)
    leftdist = self%patharray(nfirst) - self%patharray(1)
    rightdist =  self%patharray(nAtoms) - self%patharray(nfirst)
    self%schedule(1) = self%firstAtom
    iSchedule = 1
    if(leftdist > rightdist) then
       do iAtom = nfirst+1, natoms
         iSchedule = iSchedule + 1
         self%schedule(iSchedule) = self%patharray(iAtom)
       enddo
!       do iAtom = 1, nfirst-1, -1
       do iAtom = nfirst-1, 1, -1
         iSchedule = iSchedule + 1
         self%schedule(iSchedule) = self%patharray(iAtom)
       enddo
    else
       do iAtom = nfirst-1, 1, -1
         iSchedule = iSchedule + 1
         self%schedule(iSchedule) = self%patharray(iAtom)
       enddo
       do iAtom = nfirst+1, natoms
         iSchedule = iSchedule + 1
         self%schedule(iSchedule) = self%patharray(iAtom)
       enddo
    endif




  end subroutine
!==========================================================================================
  subroutine LinearCBMC_GenerateConfig(self, trialBox, disp, probconstruct, insPoint, templist, tempNNei)
    use Common_MolInfo, only: MolData, BondData, AngleData, nMolTypes
    use MolSearch, only: FindBond, FindAngle
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone
    implicit none
    class(LinearCBMC), intent(inout) :: self
    class(Addition), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)
    real(dp), intent(in), optional :: insPoint(:)
    real(dp), intent(out) :: probconstruct
    integer, intent(in), optional :: templist(:,:), tempNNei(:)

    integer :: bondType, angleType, molType
    integer :: iRosen
    integer :: atm1, atm2,atm3,atm4, iDisp, iAtom
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r
    real(dp) :: r1, r2
    real(dp) :: prob_r, prob_ang, prob_tors, probgen
    real(dp) :: ang1, ang2

    probconstruct = 1E0_dp
    select type(disp)
      class is(Addition)
        molType = disp(1)%molType
      class default
        stop "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
    end select

    if(size(insPoint) /= self%nRosenTrials) then
      write(0,*) "ERROR! Linear CBMC Regrowth received a different number of insertion points"
      write(0,*) "than it was expecting!"
      error stop
    endif


    self%grown = .false.
    nRegrown = 0


    !Weight Insertion Points and chose one of them.
    Atm1 = self%schedule(1)
    do iRosen = 1, self%nRosenTrials

    enddo
    nRegrown = nRegrown + 1
    self%grown(Atm1) = .true.

    if(natoms == 1) then
      return
    endif


    !Starting from the first inserted atom, begin growing the second atom by
    !
    Atm2 = self%schedule(2)
    call FindBond(nType, Atm1, Atm2, bondType)
    do iRosen = 1, self%nRosenTrials
      call BondData(bondType) % bondFF % GenerateDist(trialBox%beta, r2, prob_r)
      call Generate_UnitSphere(dx, dy, dz)
      self%temppos(1, iRosen) = r2*dx + disp(Atm1)%x_new 
      self%temppos(2, iRosen) = r2*dy + disp(Atm1)%y_new
      self%temppos(3, iRosen) = r2*dz + disp(Atm1)%z_new
    enddo
    self%grown(Atm2) = .true.
    nRegrown = nRegrown + 1

    if(natoms == 2) then
      return
    endif

!    For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!    we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!    then the bending angle must use the 1st atom as the central atom for the angle generation.  
!    However if it was a link in the chain then the 2nd atom regrown should be used as the central atom.
    Atm3 = self%schedule(3) 
    if(self%freq(Atm2) == 2) then
        Atm1 = self%schedule(1) 
        Atm2 = self%schedule(2) 
    else
        Atm1 = self%schedule(2) 
        Atm2 = self%schedule(1) 
    endif
    v1(1) = disp(Atm1)%x_new - disp(Atm2)%x_new
    v1(2) = disp(Atm1)%y_new - disp(Atm2)%y_new
    v1(3) = disp(Atm1)%z_new - disp(Atm2)%z_new
    call FindBond(nType, Atm2, Atm3, bondType)
    call FindAngle(nType, Atm1, Atm2, Atm3, angleType)
    do iRosen = 1, nRosenTrials(nType)
      call BondData(bondType) % bondFF % GenerateDist(trialBox%beta,r2, prob_r)
      call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, bond_ang, prob_ang)
      probgen = prob_r * prob_ang
      call Generate_UnitCone(v1, r2, bond_ang, v2)
      self%tempcoords(1, iRosen) = v2(1) + disp(Atm2)%x_new
      self%tempcoords(2, iRosen) = v2(2) + disp(Atm2)%y_new
      self%tempcoords(3, iRosen) = v2(3) + disp(Atm2)%z_new
    enddo
    nRegrown = nRegrown + 1
    self%grown(Atm3) = .true.


    if(natoms == 3) then
      return
    endif

    do while(nRegrown < nAtoms)
        nRegrown = nRegrown + 1
        atm4 = regrowOrder(nRegrown) 
        call self%FindAtomsFromPath(Atm4, Atm1, Atm2, Atm3)
        call FindBond(Atm3, Atm4, bondType)
        call FindAngle(Atm2, Atm3, Atm4, bendType)
        call FindTorsion(Atm1, Atm2, Atm3, Atm4, torsType)
        v1(1) = disp(Atm1)%x_new - disp(Atm3)%x_new
        v1(2) = disp(Atm1)%y_new - disp(Atm3)%y_new
        v1(3) = disp(Atm1)%z_new - disp(Atm3)%z_new

        v2(1) = disp(Atm2)%x_new - disp(Atm3)%x_new
        v2(2) = disp(Atm2)%y_new - disp(Atm3)%y_new
        v2(3) = disp(Atm2)%z_new - disp(Atm3)%z_new
        overlap = .false.
        E_Trial = 0E0_dp
        do iRosen = 1, nRosenTrials(nType)
          call BondData(bondType) % bondFF % GenerateDist(trialBox%beta,r2, prob)
          call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, bond_ang, prob)
          call TorsionData(torsType) % torsFF % GenerateDist(trialBox%beta, tors_ang, prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)

          self%tempcoords(1, iRosen)  = v3(1) + disp(Atm3)%x_new 
          self%tempcoords(2, iRosen)  = v3(2) + disp(Atm3)%y_new 
          self%tempcoords(3, iRosen)  = v3(3) + disp(Atm3)%z_new 
        enddo

        self%grown(Atm4) = .true.
        disp(Atm4)%x_new =  self%tempcoords(1, nSel) 
        disp(Atm4)%y_new =  self%tempcoords(2, nSel) 
        disp(Atm4)%z_new =  self%tempcoords(3, nSel) 
    enddo

    probconstruct = probconstruct*ProbRosen(nSel)*dble(self%nRosenTrials)/rosenNorm

  end subroutine
!======================================================================================
  subroutine LinearCBMC_ReverseConfig(self, trialBox, probconstruct, accept)
    implicit none
    class(LinearCBMC), intent(inout) :: self
!    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct 
    logical, intent(out) :: accept

      integer :: i, iRosen, iAtom, nSel, nIndx, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: atom1_Pos
      integer :: bondType, bendType, torsType, cnt
      logical :: overlap
      logical :: isIncluded(1:maxMol)
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial)
      real(dp) :: grnd,rotang
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob
      real(dp) :: k_bend, ang_eq, bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos
      type(SimpleAtomCoords) :: v1, v2, v3
      
      ProbRosen = 0d0      
      E_Trial = 0d0      
      nTargetMol = subIndxList(nTarget)
      nIndx = molArray(nType)%mol(nMol)%indx
      regrown = .false.

!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      E_Trial = 0d0
!      E_Complete = 0d0
      probconstruct = 1d0
!      call Rosen_BoltzWeight_Atom_Old(nType, nMol, 1, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
          !ComputeWeight
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      probconstruct = probconstruct*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
      regrown(1) = .true.


!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      Atm1 = 1
      Atm2 = regrowOrder(nType, 2) 
      call FindBond(nType,Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
!      call Rosen_BoltzWeight_Atom_Old(nType, nMol, Atm2, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        call GenerateBondLength(r, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r*dx + molArray(nType)%mol(nMol)%x(Atm1)
        trialPos%y = r*dy + molArray(nType)%mol(nMol)%y(Atm1)
        trialPos%z = r*dz + molArray(nType)%mol(nMol)%z(Atm1)
        call Rosen_BoltzWeight_Atom_New(nType, Atm2, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      probconstruct = probconstruct*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
      regrown(Atm2) = .true.

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif

      v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm2)
      v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm2)
      v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm2)
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!      k_bend = bendData(bendType)%k_eq
!      ang_eq = bendData(bendType)%ang_eq

      call Rosen_BoltzWeight_Atom_Old(nType, nMol, Atm3, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        call GenerateBondLength(r, k_bond, r_eq, Prob)
!        call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
        call GenerateBendAngle(bend_angle, bendType, Prob)
        call Generate_UnitCone(v1, r, bend_angle, v2)
        trialPos%x = v2%x + molArray(nType)%mol(nMol)%x(Atm2)
        trialPos%y = v2%y + molArray(nType)%mol(nMol)%y(Atm2)
        trialPos%z = v2%z + molArray(nType)%mol(nMol)%z(Atm2)
        call Rosen_BoltzWeight_Atom_New(nType, Atm3, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo

      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      probconstruct = probconstruct*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
      regrown(Atm3) = .true.

!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm3)
        v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm3)
        v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm3)

        v2%x = molArray(nType)%mol(nMol)%x(Atm2) - molArray(nType)%mol(nMol)%x(Atm3)
        v2%y = molArray(nType)%mol(nMol)%y(Atm2) - molArray(nType)%mol(nMol)%y(Atm3)
        v2%z = molArray(nType)%mol(nMol)%z(Atm2) - molArray(nType)%mol(nMol)%z(Atm3)
        overlap = .false.
        call Rosen_BoltzWeight_Atom_Old(nType, nMol, Atm4, isIncluded,  E_Trial(1))
        do iRosen = 2, nRosenTrials(nType)
          call GenerateBondLength(r, k_bond, r_eq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos%x = v3%x + molArray(nType)%mol(nMol)%x(Atm3)
          trialPos%y = v3%y + molArray(nType)%mol(nMol)%y(Atm3)
          trialPos%z = v3%z + molArray(nType)%mol(nMol)%z(Atm3)
          call Rosen_BoltzWeight_Atom_New(nType, Atm4, trialPos, isIncluded,  E_Trial(iRosen), overlap)
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        probconstruct = probconstruct*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
        regrown(Atm4) = .true.
      enddo

    accept = .true.
    probconstruct = 1E0_dp
  end subroutine
!======================================================================
!  Routine for simulating the probability of an isolated molecule in the gas phase
!  
  subroutine LinearCBMC_GasConfig(self,  probGas)
    implicit none
    class(LinearCBMC), intent(inout) :: self
    real(dp), intent(out) :: probGas


    probGas = 1E0_dp
  end subroutine
!=======================================================================
  subroutine LinearCBMC_FindAtomsFromPath(self, Atm4, Atm1, Atm2, Atm3)
    use Forcefield
    use SimParameters
    use Coords
    use CBMC_Variables
    implicit none
    class(LinearCBMC), intent(inout) :: self
    integer, intent(in) :: Atm4
    integer, intent(out) ::  Atm1, Atm2, Atm3

    integer :: i, atm4Pos
    integer :: atm4plus1, atm4minus1

    atm4Pos = self%pathposition(atm4)

    atm4plus1 = atm4Pos + 1
    atm4minus1 = atm4Pos - 1
    if( atm4plus1 > natoms ) then
        Atm3 = self%patharray(atm4Pos+1)
        Atm2 = self%patharray(atm4Pos+2)
        Atm1 = self%patharray(atm4Pos+3)
        return
    else if( atm4minus1 < 1 ) then
        Atm3 = self%patharray(atm4Pos-1)
        Atm2 = self%patharray(atm4Pos-2)
        Atm1 = self%patharray(atm4Pos-3)
        return
    endif

    if( self%grown(atm4plus1) ) then
        Atm3 = self%patharray(atm4Pos+1)
        Atm2 = self%patharray(atm4Pos+2)
        Atm1 = self%patharray(atm4Pos+3)
    else if( self%grown(atmminus1) ) then
        Atm3 = self%patharray(atm4Pos-1)
        Atm2 = self%patharray(atm4Pos-2)
        Atm1 = self%patharray(atm4Pos-3)
    else
        write(0,*) "ERROR in FindAtomsFromPath subroutine! The atom is being regrown before any neighboring atoms have! "
        write(0,*) "Atom ID:", Atm4
        error stop
    endif

    if( (.not. self%grown(Atm3)) .or.
        (.not. self%grown(Atm2)) .or.
        (.not. self%grown(Atm1)) ) then
        write(0,*) "ERROR in FindAtomsFromPath subroutine! The atom used for torsional selection has not be grown yet!"
        write(0,*) "Atom ID:", Atm4
        error stop
    endif
  end subroutine 
!==========================================================================================
end module
!==========================================================================================
