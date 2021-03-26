!==========================================================================================
! Simple Regrowth Object
!==========================================================================================
module MolCon_LinearCBMC
  use CoordinateTypes, only: Perturbation, Addition, Displacement
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: LinearCBMC
!    inspoints => Number of insertion points to request from the move
!                 calling this algorithm. Inherited from parent class.
!    integer :: inspoints = 1 

    !nAtoms => Number of Atoms in this molecule.
    !firstAtom => Relative Atom Index of the first atom to be inserted into the system
    !nRosenTrial => Number of Rosenbluth trials to use.  
    !                This is the number of potential positions
    !                that are generated for each new atom as it is being regrown.  
    !                One of these positions will be selected for the final atom position.  
    integer,private :: nAtoms 
    integer,private :: firstAtom = 1 
    integer, private :: rosenNeighList = 1
    integer,private :: nRosenTrials = 8 

    !GenProb => Weight array for chosing which trial position to use
    !RosenProb => Weight array for chosing which trial position to use
    !tempcoords => Temporary storage bucket for trial atom positions
    real(dp), private, allocatable :: GenProb(:) 
    real(dp), private, allocatable :: RosenProb(:) 
    real(dp), private, allocatable :: tempcoords(:, :) 
    real(dp), private, allocatable :: newconfig(:, :) 

    !nGrown => Counter which tracks how many atoms are already regrown. 
    !grown => Array which specifies if an atom has been added to the system or not
    !patharray => Array which lays out the topology of the molecule 
    !pathposition => Array which takes an atom's ID and gives the position in the patharray
    !schedule => !The order in which atoms in a molecule will be grown.
    !scratchschedule => !The order in which a molecule will be grown if growing from scratch.
    integer, private :: nGrown
    logical, private, allocatable :: grown(:) 
    integer, private, allocatable :: patharray(:) 
    integer, private, allocatable :: pathposition(:) 
    integer, private, allocatable :: schedule(:) 
    integer, private, allocatable :: scratchschedule(:) 

    integer, private, allocatable :: freq(:) 
    !freq => The number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds  = 0 -> Isolated Atom
    ! nBonds  = 1 -> Terminal Atom located on the end of a chain
    ! nBonds  = 2 -> Linker Atom located in the middle of a linear chain 
    ! nBonds >= 3 -> Branch Atom has two or more potential paths one can wander down.

    contains
!      procedure, public, pass :: Constructor => LinearCBMC_Constructor
      procedure, public, pass :: Prologue => LinearCBMC_Prologue
      procedure, public, pass :: CreateSchedule => LinearCBMC_CreateSchedule
      procedure, public, pass :: GenerateConfig => LinearCBMC_GenerateConfig
      procedure, public, pass :: ReverseConfig => LinearCBMC_ReverseConfig
      procedure, public, pass :: RosenBluth => LinearCBMC_RosenBluth
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
    integer :: iType, iBond, iAtom, curMax
    integer :: atm1, atm2
    integer :: iError = 0
    integer :: iSchedule
    integer :: nextAtm, prevAtm, curAtm, iPath
    integer :: leftdist, rightdist, nleft, nright, nfirst

    !Count the number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds  = 0 >> Isolated Atom
    ! nBonds  = 1 >> Terminal Atom located on the end of a chain
    ! nBonds  = 2 >> Linker Atom. Located in the middle of a linear chain
    ! nBonds >= 3 >> Branch Atom. Has two or more potential paths one can wander down.
    self%nAtoms = MolData(self%molType)%nAtoms
    if(self%nAtoms > 1) then
      allocate( self%freq(1:self%nAtoms) )  
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
      error stop
    endif

    if(.not. allocated(self%RosenProb)) then
      allocate(self%GenProb(1:self%nRosenTrials))
      allocate(self%RosenProb(1:self%nRosenTrials))
      allocate(self%tempcoords(1:3, 1:self%nRosenTrials))
    endif


    !Create the path array which is used to determine regrowth order.
    if(.not. allocated(self%patharray) ) then
      allocate( self%grown(1:self%nAtoms) )  
      allocate( self%patharray(1:self%nAtoms) )  
      allocate( self%pathposition(1:self%nAtoms) )  
      allocate( self%schedule(1:self%nAtoms) )  
      allocate( self%scratchschedule(1:self%nAtoms) )  
    endif
    self%patharray = 0
    self%pathposition = 0
    !Pick an terminal atom to start building the patharray
    prevAtm = 0
    do iAtom = 1,self%nAtoms
      if(self%freq(iAtom) == 1) then
        curatm = iAtom
        exit
      endif
    enddo

    self%patharray(1) = curatm
    self%pathposition(curatm) = 1
    iPath = 1
    ! Begin the building process by traversing from bond to bond.
    do iAtom = 1, self%nAtoms-1
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
    rightdist =  self%patharray(self%nAtoms) - self%patharray(nfirst)
    self%scratchschedule(1) = self%firstAtom
    iSchedule = 1
    if(leftdist > rightdist) then
       do iAtom = nfirst+1, self%nAtoms
         iSchedule = iSchedule + 1
         self%scratchschedule(iSchedule) = self%patharray(iAtom)
       enddo
       do iAtom = nfirst-1, 1, -1
         iSchedule = iSchedule + 1
         self%scratchschedule(iSchedule) = self%patharray(iAtom)
       enddo
    else
       do iAtom = nfirst-1, 1, -1
         iSchedule = iSchedule + 1
         self%scratchschedule(iSchedule) = self%patharray(iAtom)
       enddo
       do iAtom = nfirst+1, self%nAtoms
         iSchedule = iSchedule + 1
         self%scratchschedule(iSchedule) = self%patharray(iAtom)
       enddo
    endif




  end subroutine
!==========================================================================================
  subroutine LinearCBMC_GenerateConfig(self, trialBox, disp, probconstruct, insPoint, insProb)
    use Common_MolInfo, only: MolData, BondData, AngleData, TorsionData, nMolTypes
    use MolSearch, only: FindBond, FindAngle, FindTorsion
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone, Generate_UnitTorsion
    use ForcefieldData, only: ECalcArray

    implicit none
    class(LinearCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)
    real(dp), intent(in), optional :: insProb(:)
    real(dp), intent(out) :: probconstruct 

    class(ECalcArray), pointer :: EFunc => null()
    integer :: bondType, angleType, torsType, molType
    integer :: molindx, molStart
    integer :: atm1, atm2,atm3,atm4, iDisp, iRosen
    integer :: lastGrown, nSel
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r
    real(dp) :: r1, r2
    real(dp) :: bend_angle,tors_angle
    real(dp) :: prob_r, prob_ang, prob_tors, probgen



    call trialbox%GetEFunc(EFunc)

    probconstruct = 1E0_dp
    select type(disp)
      class is(Displacement)
        self%grown = .true.
        self%nGrown = self%nAtoms
        molindx = disp(1)%molindx
        call trialbox%GetMolData(molindx, molStart=molStart)
        do iDisp = 1, size(disp)
          atm1 = disp(iDisp)%atmindx - molStart + 1
          self%grown(atm1) = .false.
          self%nGrown = self%nGrown - 1
        enddo

      class is(Addition)
        self%grown = .false.
        self%nGrown = 0
        self%schedule = self%scratchschedule

      class default
        error stop "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
    end select

    call self%CreateSchedule
    if(present(insPont)) then
      if(size(insPoint) /= self%nRosenTrials) then
        write(0,*) "ERROR! Linear CBMC Regrowth received a different number of insertion points"
        write(0,*) "than it was expecting!"
        error stop
      endif
    endif


    do while(self%nGrown < self%nAtoms)
        !Weight Insertion Points and chose one of them.
      if(self%nGrown == 0) then
          Atm1 = self%schedule(1)
          if(present(inspoint) .and. present(insProb)) then
            do iRosen = 1, self%nRosenTrials
              self%tempcoords(1, iRosen) = insPoint(3*iRosen-2)
              self%tempcoords(2, iRosen) = insPoint(3*iRosen-1)
              self%tempcoords(3, iRosen) = insPoint(3*iRosen-0)
              self%GenProb(iRosen) = insProb(iRosen)
            enddo
          else
            error stop "Full Regrowth has been requsted without passing in insertion points"
          endif 
          lastGrown = Atm1

      elseif(self%nGrown == 1) then
        !Starting from the first inserted atom, begin growing the second atom by
        Atm1 = self%schedule(1)
        Atm2 = self%schedule(2)
        call FindBond(self%molType, Atm1, Atm2, bondType)
        do iRosen = 1, self%nRosenTrials
          call BondData(bondType) % bondFF % GenerateDist(trialBox%beta, r2, prob_r)
          probgen = prob_r
          call Generate_UnitSphere(dx, dy, dz)
          self%tempcoords(1, iRosen) = r2*dx + self%newconfig(1, Atm1)
          self%tempcoords(2, iRosen) = r2*dy + self%newconfig(2, Atm1)
          self%tempcoords(3, iRosen) = r2*dz + self%newconfig(3, Atm1)
          self%GenProb(iRosen) = prob_r
        enddo
        lastGrown = Atm2

      elseif(self%nGrown == 2) then
    !    For the third atom we must begin to include the bending angle 
    !    in choosing its trial positions. The first thing
    !    we must consider is if the second regrown atom was a terminal atom or a linker atom. 
    !    If it was a terminal atom  then the bending angle must use the 
    !    1st atom as the central atom for the angle generation.  
    !    However if it was a link in the chain then the 2nd atom 
    !    regrown should be used as the central atom.
        Atm3 = self%schedule(3) 
        if(self%freq(Atm2) == 2) then
            Atm1 = self%schedule(1) 
            Atm2 = self%schedule(2) 
        else
            Atm1 = self%schedule(2) 
            Atm2 = self%schedule(1) 
        endif
        v1(1) = self%newconfig(1, Atm1) - self%newconfig(1, Atm2)
        v1(2) = self%newconfig(2, Atm1) - self%newconfig(2, Atm2)
        v1(3) = self%newconfig(3, Atm1) - self%newconfig(3, Atm2)
        call FindBond(self%molType, Atm2, Atm3, bondType)
        call FindAngle(self%molType, Atm1, Atm2, Atm3, angleType)
        do iRosen = 1, self%nRosenTrials
          call BondData(bondType) % bondFF % GenerateDist(trialBox%beta,r2, prob_r)
          call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, bend_angle, prob_ang)
          probgen = prob_r * prob_ang
          call Generate_UnitCone(v1, r2, bend_angle, v2)
          self%tempcoords(1, iRosen) = v2(1) + self%newconfig(1, Atm2)
          self%tempcoords(2, iRosen) = v2(2) + self%newconfig(2, Atm2)
          self%tempcoords(3, iRosen) = v2(3) + self%newconfig(3, Atm2)
          self%GenProb(iRosen) = probgen
        enddo
        lastGrown = Atm3

      elseif(self%nGrown < self%nAtoms) then
            Atm4 = self%schedule(self%nGrown+1) 
            call self%FindAtomsFromPath(Atm4, Atm1, Atm2, Atm3)
            call FindBond(self%molType, Atm3, Atm4, bondType)
            call FindAngle(self%molType, Atm2, Atm3, Atm4, angleType)
            call FindTorsion(self%molType, Atm1, Atm2, Atm3, Atm4, torsType)
            v1(1) = self%newconfig(1, Atm1) - self%newconfig(1, Atm3)
            v1(2) = self%newconfig(2, Atm1) - self%newconfig(2, Atm3)
            v1(3) = self%newconfig(3, Atm1) - self%newconfig(3, Atm3)


            v1(1) = self%newconfig(1, Atm2) - self%newconfig(1, Atm3)
            v1(2) = self%newconfig(2, Atm2) - self%newconfig(2, Atm3)
            v1(3) = self%newconfig(2, Atm2) - self%newconfig(3, Atm3)
            do iRosen = 1, self%nRosenTrials
              call BondData(bondType) % bondFF % GenerateDist(trialBox%beta, r2, prob_r)
              call AngleData(angleType) % angleFF % GenerateDist(trialBox%beta, bend_angle, prob_ang)
              call TorsionData(torsType) % torsionFF % GenerateDist(trialBox%beta, tors_angle, prob_tors)
              probgen = prob_r*prob_ang*prob_tors
              call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
              self%tempcoords(1, iRosen)  = v3(1) + self%newconfig(1, Atm3)
              self%tempcoords(2, iRosen)  = v3(2) + self%newconfig(2, Atm3)
              self%tempcoords(3, iRosen)  = v3(3) + self%newconfig(3, Atm3)
              self%GenProb(iRosen) = probgen
            enddo
            lastGrown = Atm4
        else
            error stop "nGrown is some invalid number."
        endif

        self%grown(lastGrown) = .true.
        self%nGrown = self%nGrown + 1


     enddo

     do iDisp = 1, size(disp)
      select type(disp)
       class is(Displacement)
          atmindx = disp(iDisp)%atmindx - molStart + 1
          disp(iDisp)%x_new = self%tempcoords(1, atmindx) 
          disp(iDisp)%y_new = self%tempcoords(2, atmindx) 
          disp(iDisp)%z_new = self%tempcoords(3, atmindx) 
       class is(Addition)
          atmindx = disp(iDisp)%atmindx - molStart + 1
          disp(iDisp)%x_new = self%tempcoords(1, atmindx) 
          disp(iDisp)%y_new = self%tempcoords(2, atmindx) 
          disp(iDisp)%z_new = self%tempcoords(3, atmindx) 
       end select
    enddo

!    probconstruct = probconstruct*ProbRosen(nSel)*dble(self%nRosenTrials)/rosenNorm

  end subroutine
!======================================================================================
  subroutine LinearCBMC_ReverseConfig(self, disp, trialBox, probconstruct, insPoint, insProb, accept)
    implicit none
    class(LinearCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct 

    real(dp), intent(in), optional :: insPoint(:)
    real(dp), intent(in), optional :: insProb(:)
    logical, intent(out) :: accept

  end subroutine
!======================================================================
!  Routine for simulating the probability of an isolated molecule in the gas phase
!  Used primarily for swap moves with an implict gas box.
  subroutine LinearCBMC_GasConfig(self,  probGas)
    implicit none
    class(LinearCBMC), intent(inout) :: self
    real(dp), intent(out) :: probGas


    probGas = 1E0_dp
  end subroutine
!======================================================================
!  Routine for simulating the probability of an isolated molecule in the gas phase
!  Used primarily for swap moves with an implict gas box.
  subroutine LinearCBMC_CreateSchedule(self)
    implicit none
    class(LinearCBMC), intent(inout) :: self
    integer :: iPath, iAtom, iSchedule, lastatom
    logical :: lastgrown
    logical :: reverse
    integer :: breakpoint


    lastgrown = self%grown(self%patharray(1))
    lastatom = self%patharray(1)
    !First we have to figure out what's still left to regrow. We find the point
    !where the first ungrown atom is in the path.
    do iPath = 2, self%nAtoms
      iAtom = self%patharray(iPath)
      if(self%grown(iAtom) .neqv. lastgrown) then
        if(.not. self%grown(iAtom)) then
          reverse = .false.
          breakpoint = iAtom
        else
          reverse = .true.
          breakpoint = lastatom
        endif
        exit
      endif
      lastatom = iAtom
    enddo


    if(.not. reverse) then

      do iPath = 1, self%nAtoms
        self%schedule(iPath) = self%patharray(iPath)
      enddo 
    else
      do iPath = 1, self%nAtoms
        self%schedule(iPath) = self%patharray(self%nAtoms-iPath+1)
      enddo
    endif


  end subroutine
!=======================================================================
  subroutine LinearCBMC_FindAtomsFromPath(self, Atm4, Atm1, Atm2, Atm3)
    ! Figures out which atoms
    implicit none
    class(LinearCBMC), intent(inout) :: self
    integer, intent(in) :: Atm4
    integer, intent(out) ::  Atm1, Atm2, Atm3

    integer :: i, atm4Pos
    integer :: atm4plus1, atm4minus1

    atm4Pos = self%pathposition(atm4)

    atm4plus1 = atm4Pos + 1
    atm4minus1 = atm4Pos - 1
    if( atm4plus1 > self%nAtoms ) then
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
    else if( self%grown(atm4minus1) ) then
        Atm3 = self%patharray(atm4Pos-1)
        Atm2 = self%patharray(atm4Pos-2)
        Atm1 = self%patharray(atm4Pos-3)
    else
        write(0,*) "ERROR in FindAtomsFromPath subroutine! The atom is being regrown before any neighboring atoms have! "
        write(0,*) "Atom ID:", Atm4
        error stop
    endif

    if( (.not. self%grown(Atm3)) .or. &
        (.not. self%grown(Atm2)) .or. &
        (.not. self%grown(Atm1)) ) then
        write(0,*) "ERROR in FindAtomsFromPath subroutine! The atom used for torsional selection has not be grown yet!"
        write(0,*) "Atom ID:", Atm4
        error stop
    endif
  end subroutine 
!=======================================================================
  subroutine LinearCBMC_RosenBluth(self, trialBox, disp)
    ! Figures out which atoms
    implicit none
    class(LinearCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    type(Displacement) :: tempdisp(1:1)

    integer :: tempList(1:100,1:self%nAtoms), tempNNei(1:self%nAtoms)

    tempdisp


    call trialbox%GetNewList(self%rosenNeighList, 1, tempList, tempNNei, disp)



  end subroutine 
!=======================================================================
end module
!==========================================================================================
