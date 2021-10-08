!========================================================
! Standard Molecular Translation move with the ability to restrict the
! degrees of motion by turning off movement in the x,y, or z directions. 
! This can be used to constraint to the system to 2D or 1D gemoetries. 
!========================================================
module MCMove_PlaneTranslation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: PlaneTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

    real(dp), allocatable :: boxatmps(:)
    real(dp) , allocatable:: boxaccpt(:)
    logical :: proportional = .true.
    logical :: tuneMax = .true.
    logical :: xdir = .true.
    logical :: ydir = .true.
    logical :: zdir = .true.
    real(dp) :: limit = 3.00E0_dp
    real(dp) :: targAccpt = 50E0_dp
    real(dp) :: max_dist = 0.05E0_dp
    real(dp), allocatable :: boxlimit(:)
    real(dp), allocatable :: boxtargAccpt(:)
    real(dp), allocatable :: boxmax_dist(:)

    type(Displacement), allocatable :: disp(:)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => PlaneTranslate_Constructor
!      procedure, pass :: GeneratePosition => PlaneTranslate_GeneratePosition
      procedure, pass :: FullMove => PlaneTranslate_FullMove
      procedure, pass :: Maintenance => PlaneTranslate_Maintenance
      procedure, pass :: Prologue => PlaneTranslate_Prologue
      procedure, pass :: Update => PlaneTranslate_Update
      procedure, pass :: Epilogue => PlaneTranslate_Epilogue
      procedure, pass :: ProcessIO => PlaneTranslate_ProcessIO
  end type
!========================================================
 contains
!========================================================
  subroutine PlaneTranslate_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    use BoxData, only: BoxArray
    implicit none
    class(PlaneTranslate), intent(inout) :: self
    integer :: iType, maxAtoms, nBoxes

    nBoxes = size(boxArray)
    if(.not. allocated(self%boxProb)) then
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    allocate( self%boxatmps(1:nBoxes) )
    self%boxatmps = 1e-50_dp
    allocate( self%boxaccpt(1:nBoxes) )
    self%boxaccpt = 0E0_dp

    allocate( self%boxLimit(1:nBoxes) )
    self%boxLimit = self%limit
    allocate( self%boxmax_dist(1:nBoxes) )
    self%boxmax_dist = self%max_dist
    allocate( self%boxtargAccpt(1:nBoxes) )
    self%boxtargAccpt = self%targAccpt

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
!  subroutine PlaneTranslate_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(PlaneTranslate), intent(in) :: self
!    type(Displacement), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine PlaneTranslate_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use RandomGen, only: grnd
    implicit none
    class(PlaneTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: boxID, iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, E_Inter, E_Intra, biasE
    real(dp), parameter :: Prob = 1E0_dp

    boxID = trialBox % boxID
    self % atmps = self % atmps + 1E0_dp
    self % boxatmps(boxID) = self % boxatmps(boxID) + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp )
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)


    dx = 0E0_dp
    dy = 0E0_dp
    dz = 0E0_dp

    if(self%xdir) then
      dx = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    endif
    if(self%ydir) then
      dy = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    endif
    if(self%zdir) then
      dz = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    endif
 
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
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:nAtoms), self%tempList, self%tempNNei, E_Diff, accept)
    call trialBox%ComputeEnergyDelta(self%disp(1:nAtoms),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.false.)

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
      self % boxaccpt(boxID) = self % boxaccpt(boxID) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine PlaneTranslate_Maintenance(self)
    implicit none
    class(PlaneTranslate), intent(inout) :: self
    integer :: iBox
    real(dp) :: accRate
!    real(dp), parameter :: lowerlimit = 0.1E0_dp
      
    if(self%tuneMax) then
      do iBox = 1, size(self%boxatmps)
        if(self%boxatmps(iBox) < 0.5E0_dp) then
          cycle
        endif
        accRate = 1e2_dp*self%boxaccpt(iBox)/self%boxatmps(iBox)

        if(accRate > self%boxtargAccpt(iBox)) then
          if(self%boxmax_dist(iBox)*1.01E0_dp < self%boxlimit(iBox)) then
            self%boxmax_dist(iBox) = self%boxmax_dist(iBox) * 1.01E0_dp
          else 
            self%boxmax_dist(iBox) = self%boxlimit(iBox)
          endif
        else
          self%boxmax_dist(iBox) = self%boxmax_dist(iBox) * 0.99E0_dp
        endif
      enddo
    endif

  end subroutine
!=========================================================================
  subroutine PlaneTranslate_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(PlaneTranslate), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      

    write(nout,"(1x,A,F15.8)") "(Plane Translate) Maximum Displacement: ", self%max_dist

  end subroutine
!=========================================================================
  subroutine PlaneTranslate_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(PlaneTranslate), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Plane Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Plane Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Plane Translation Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,100F15.8)") "Final Maximum Displacement: ", self%boxmax_dist(1:)
    endif
 

  end subroutine
!=========================================================================
  subroutine PlaneTranslate_Update(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(PlaneTranslate), intent(inout) :: self
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
  subroutine PlaneTranslate_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(PlaneTranslate), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("directions")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%xDir = logicVal

        call GetXCommand(line, command, 6, lineStat)
        read(command, *) logicVal
        self%yDir = logicVal

        call GetXCommand(line, command, 7, lineStat)
        read(command, *) logicVal
        self%zDir = logicVal

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
