!========================================================
module MCMove_MolTranslation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: MolTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    logical :: tuneMax = .true.
    real(dp) :: limit = 3.00E0_dp
    real(dp) :: targAccpt = 50E0_dp
    real(dp) :: max_dist = 0.05E0_dp
    type(Displacement), allocatable :: disp(:)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => MolTrans_Constructor
!      procedure, pass :: GeneratePosition => MolTrans_GeneratePosition
      procedure, pass :: FullMove => MolTrans_FullMove
      procedure, pass :: Maintenance => MolTrans_Maintenance
      procedure, pass :: Prologue => MolTrans_Prologue
      procedure, pass :: Epilogue => MolTrans_Epilogue
      procedure, pass :: ProcessIO => MolTrans_ProcessIO
  end type
!========================================================
 contains
!========================================================
  subroutine MolTrans_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    use BoxData, only: BoxArray
    implicit none
    class(MolTranslate), intent(inout) :: self
    integer :: iType, maxAtoms, nBoxes

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo


    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(1000, maxAtoms) )
    allocate( self%disp(1:maxAtoms) )
  end subroutine
!========================================================
!  subroutine MolTrans_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(MolTranslate), intent(in) :: self
!    type(Displacement), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine MolTrans_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(MolTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE
    real(dp), parameter :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp )
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)

    dx = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dy = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
    dz = self % max_dist * (2E0_dp * grnd() - 1E0_dp)
 
    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1

      self%disp(iAtom)%molType = molType
      self%disp(iAtom)%molIndx = nMove
      self%disp(iAtom)%atmIndx = atomIndx

      self%disp(iAtom)%x_new = trialBox%atoms(1, atomIndx) + dx
      self%disp(iAtom)%y_new = trialBox%atoms(2, atomIndx) + dy
      self%disp(iAtom)%z_new = trialBox%atoms(3, atomIndx) + dz

      self%disp(iAtom)%newlist = .false.
      self%disp(iAtom)%listIndex = iAtom
    enddo
!    write(*,*) dx, dy, dz

    !If the particle moved a large distance get a temporary neighborlist
!    if(any([dx,dy,dz] > neighSkin)) then
!      call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%disp(1))
!      self%disp(1)%newlist = .true.
!    else

!    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:nAtoms), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:nAtoms), E_Diff )
    if(.not. accept) then
      return
    endif




    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:nAtoms), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine MolTrans_Maintenance(self)
    implicit none
    class(MolTranslate), intent(inout) :: self
!    real(dp), parameter :: limit = 3.0E0_dp
      
    if(self%tuneMax) then
      if(self%atmps .lt. 0.5E0_dp) then
        return
      endif

      if(self%GetAcceptRate() > self%targAccpt) then
        if(self%max_dist*1.01E0_dp .lt. self%limit) then
          self%max_dist = self%max_dist * 1.01E0_dp
        else 
          self%max_dist = self%limit       
        endif
      else
        self%max_dist = self%max_dist * 0.99E0_dp
      endif
    endif

  end subroutine
!=========================================================================
  subroutine MolTrans_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(MolTranslate), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      

    write(nout,"(1x,A,F15.8)") "(Molecule Translate) Maximum Displacement: ", self%max_dist

  end subroutine
!=========================================================================
  subroutine MolTrans_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(MolTranslate), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Molecule Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Molecule Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Molecule Translation Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,F15.8)") "Final Maximum Displacement: ", self%max_dist
    endif
 

  end subroutine
!=========================================================================
  subroutine MolTrans_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(MolTranslate), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("tunemax")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%tunemax = logicVal

      case("dynamiclimit")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%limit = realVal

      case("dynamictarget")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%targAccpt = realVal

      case("maxdisplace")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%max_dist = realVal


      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine
!========================================================
end module
!========================================================
