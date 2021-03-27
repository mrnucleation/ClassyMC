!========================================================
! Move which cuts and regrows a molecule in conjunction with
! the CBMC MolCon Regrowth Functions.  
!========================================================
module MCMove_CBMC
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: CBMC
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

    real(dp), private, allocatable :: boxatmps(:)
    real(dp), private , allocatable:: boxaccpt(:)
    logical, private :: verbose = .true.
    logical, private :: proportional = .true.

    logical, private, allocatable :: validtype(:)
    integer, private, allocatable :: patharrays(:, :)

    type(Displacement), allocatable :: disp(:)

    !Rejection Counters
    integer, private :: ovlaprej = 0 
    integer, private :: constrainrej = 0 
    integer, private :: detailedrej = 0 

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => CBMC_Constructor
!      procedure, pass :: GeneratePosition => CBMC_GeneratePosition
      procedure, pass :: FullMove => CBMC_FullMove
!      procedure, pass :: Maintenance => CBMC_Maintenance
      procedure, pass :: Prologue => CBMC_Prologue
      procedure, pass :: Update => CBMC_Update
      procedure, pass :: Epilogue => CBMC_Epilogue
      procedure, pass :: ProcessIO => CBMC_ProcessIO
  end type
!========================================================
 contains
!========================================================
  subroutine CBMC_Constructor(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes
    use MolCon_LinearCBMC, only: LinearCBMC
    use Template_MolConstructor, only: MolConstructor
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self
    integer :: iType, maxAtoms, nAtoms, nBoxes

    nBoxes = size(boxArray)
    if(.not. allocated(self%boxProb)) then
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    allocate( self%validtype(1:nMolTypes) )

    allocate( self%boxatmps(1:nBoxes) )
    allocate( self%boxaccpt(1:nBoxes) )
    self%boxatmps = 1e-50_dp
    self%boxaccpt = 0E0_dp

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo


    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(1000, maxAtoms) )
    allocate( self%disp(1:maxAtoms) )
    allocate( self%patharrays(1:maxAtoms, 1:nMolTypes) )
    self%patharrays = 0
    self%validtype = .true.

    do iType = 1, nMolTypes
      select type(molcon => MolData(iType)%molConstruct)
        class is(LinearCBMC)
          nAtoms = MolData(iType)%nAtoms 
          call molcon % GetPath( self%patharrays(1:nAtoms, iType) )

        class default
          self%validtype(iType) = .false.
          error stop
      end select
    enddo

    write(nout,*) "CBMC Valid Molecule Types: ", self%validtype(1:nMolTypes)
    if( all(.not. self%validtype) ) then
      stop "CBMC Moves have been specified on molecules which do not use a CBMC regrowth style"
    endif


  end subroutine
!===============================================
  subroutine CBMC_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(CBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    logical :: reverse
    integer :: boxID, iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    integer :: nRegrow, cutpoint
    integer :: lowIndx, highIndx, iDisp
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, E_Inter, E_Intra, biasE
    real(dp) :: Prob, ProbSub

    boxID = trialBox % boxID
    self % atmps = self % atmps + 1E0_dp
    self % boxatmps(boxID) = self % boxatmps(boxID) + 1E0_dp
    accept = .true.

    !Propose move
    nMove = self%UniformMoleculeSelect(trialBox, restrict=self%validtype(1:nMolTypes))
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)
    nAtoms = MolData(molType)%nAtoms

    cutpoint = floor( nAtoms * grnd() + 1E0_dp )
    if( grnd() < 0.5E0_dp) then
      reverse = .false.
    else
      reverse = .true.
    endif

 
    if(reverse) then
      if(cutpoint /= 1) then
        lowIndx = cutpoint
      else
        lowIndx = 2
      endif
      highIndx = nAtoms
    else
      lowIndx = 1
      if(cutpoint /= nAtoms) then
        highIndx = cutpoint
      else
        highIndx = nAtoms - 1
      endif     
    endif

    iDisp = 0
    nRegrow = 0
    do iAtom = lowIndx, highIndx
      iDisp = iDisp + 1
      atomIndx = molStart + iAtom - 1
      nRegrow = nRegrow + 1

      self%disp(iDisp)%molType = molType
      self%disp(iDisp)%molIndx = nMove
      self%disp(iDisp)%atmIndx = self%patharrays(iAtom, molType)

      self%disp(iDisp)%x_new = 0E0_dp
      self%disp(iDisp)%y_new = 0E0_dp
      self%disp(iDisp)%z_new = 0E0_dp

      self%disp(iDisp)%newlist = .false.
      self%disp(iDisp)%listIndex = iDisp
    enddo

    write(*,*) nRegrow
    call MolData(molType) % molConstruct % GenerateConfig(trialBox, self%disp(1:nRegrow), ProbSub)
    Prob = 1E0_dp/ProbSub
    write(*,*) "Prob:", prob
    iDisp = 0
    do iAtom = lowIndx, highIndx
      iDisp = iDisp + 1
      write(*,*) self%disp(iDisp)%x_new, self%disp(iDisp)%y_new, self%disp(iDisp)%z_new 
    enddo


    !If the particle moved a large distance get a temporary neighborlist
!    if(any([dx,dy,dz] > neighSkin)) then
!      call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%disp(1))
!      self%disp(1)%newlist = .true.
!    else

!    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nRegrow) )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif

    !Energy Calculation
    call trialBox%ComputeEnergyDelta(self%disp(1:nRegrow),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)
    write(*,*) "E_Diff", E_Diff
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:nRegrow), E_Diff )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif



    call MolData(molType) % molConstruct % ReverseConfig(self%disp(1:nRegrow), trialBox, ProbSub, accept)

    Prob = Prob * ProbSub
    write(*,*) Prob

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:nRegrow), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % boxaccpt(boxID) = self % boxaccpt(boxID) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    else
      self%detailedrej = self%detailedrej + 1
!      write(*,*) E_Diff, trialBox%beta, Prob
    endif

  end subroutine
!=========================================================================
  subroutine CBMC_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      

  end subroutine
!=========================================================================
  subroutine CBMC_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "CBMC Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "CBMC Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "CBMC Acceptance Rate: ", accptRate


    if(self%verbose) then
      write(nout, "(1x,A,I15)") "CBMC, Rejections due to overlap:", self%ovlaprej
      write(nout, "(1x,A,I15)") "CBMC, Rejections due to constraint:", self%constrainrej
      write(nout, "(1x,A,I15)") "CBMC, Rejections due to detailed balance:", self%detailedrej
    endif

  end subroutine
!=========================================================================
  subroutine CBMC_Update(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self
    integer :: iBox
    real(dp) :: norm

      if(self%proportional) then
        norm = 0E0_dp
        do iBox = 1, size(BoxArray)
          norm = norm + real(BoxArray(ibox)%box%nMolTotal,dp)
        enddo
        do iBox = 1, size(BoxArray)
          self%boxprob(iBox) = real(BoxArray(ibox)%box%nMolTotal,dp)/norm
        enddo
      endif
  end subroutine
!=========================================================================
  subroutine CBMC_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(CBMC), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("proportional")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%proportional = logicVal


      case("updatefreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%maintFreq = intVal

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine
!========================================================
end module
!========================================================
