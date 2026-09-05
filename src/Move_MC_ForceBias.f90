!=========================================================================
! Force-bias / smart Monte Carlo (Rossky, Doll, Friedman 1978).
! Collective same-N move using NewState_IsoMol:
!   r' = r + lambda * F(r) + sigma * Gaussian()
! Reverse kernel uses F(r') on the trial coordinates.
! If lambda is not set, lambda = 0.5 * sigma^2 * beta (theoretically consistent).
! Cluster constraints (distancecriteria / allmoldistancecriteria) see the
! full trial coordinate array and rebuild the cluster graph from it.
!=========================================================================
module MCMove_ForceBias
  use CoordinateTypes, only: NewState_IsoMol
  use SimpleSimBox, only: SimpleBox
  use VarPrecision
  use MoveClassDef

  type, public, extends(MCMove) :: ForceBias
    logical :: verbose = .true.
    logical :: autoLambda = .true.
    integer :: molType = 0
    real(dp) :: sigma = 0.10E0_dp
    real(dp) :: lambda = 0E0_dp

    integer :: ovlaprej = 0
    integer :: constrainrej = 0
    integer :: detailedrej = 0

    type(NewState_IsoMol) :: disp(1:1)
    real(dp), allocatable :: forces_old(:,:)
    real(dp), allocatable :: forces_new(:,:)
    contains
      procedure, pass :: Constructor => ForceBias_Constructor
      procedure, pass :: FullMove => ForceBias_FullMove
      procedure, pass :: Prologue => ForceBias_Prologue
      procedure, pass :: Epilogue => ForceBias_Epilogue
      procedure, pass :: ProcessIO => ForceBias_ProcessIO
  end type

 contains
!========================================================
  subroutine ForceBias_Constructor(self)
    use BoxData, only: BoxArray
    implicit none
    class(ForceBias), intent(inout) :: self
    integer :: iBox, nBoxes, nMax

    nBoxes = size(BoxArray)
    if(.not. allocated(self%boxProb)) then
      allocate(self%boxProb(1:nBoxes))
      self%boxProb = 1E0_dp / real(nBoxes, dp)
    endif

    nMax = 0
    do iBox = 1, nBoxes
      nMax = max(nMax, BoxArray(iBox)%box%nMaxAtoms)
    enddo

    if(.not. allocated(self%disp(1)%newAtoms)) then
      allocate(self%disp(1)%newAtoms(1:3, 1:nMax))
    endif
    if(.not. allocated(self%forces_old)) then
      allocate(self%forces_old(1:3, 1:nMax))
    endif
    if(.not. allocated(self%forces_new)) then
      allocate(self%forces_new(1:3, 1:nMax))
    endif
    self%disp(1)%newAtoms = 0E0_dp
    self%forces_old = 0E0_dp
    self%forces_new = 0E0_dp

    call self%CreateTempArray(1)
  end subroutine
!=========================================================================
  subroutine ForceBias_FullMove(self, trialBox, accept)
    use CommonSampling, only: sampling
    use RandomGen, only: Gaussian
    implicit none
    class(ForceBias), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: iAtom, nMax, nMoved
    real(dp) :: lambda, inv2sig2
    real(dp) :: dx, dy, dz
    real(dp) :: ssq_fwd, ssq_rev, logProb, extraTerms
    real(dp) :: E_Diff, E_Inter, E_Intra
    logical :: moveAtom

    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    if(.not. allocated(self%disp(1)%newAtoms)) then
      call self%Constructor
    endif

    call self%LoadBoxInfo(trialBox, self%disp)

    nMax = min(trialBox%nMaxAtoms, size(self%disp(1)%newAtoms, 2))
    if(self%autoLambda) then
      lambda = 0.5E0_dp * self%sigma * self%sigma * trialBox%beta
    else
      lambda = self%lambda
    endif

    if(trialBox%EFunc%Method%HasComputeForces()) then
      call trialBox%EFunc%Method%ComputeForces(trialBox, self%forces_old, trialBox%atoms)
    else
      call trialBox%EFunc%Method%ComputeForcesFiniteDiff(trialBox, self%forces_old, trialBox%atoms)
    endif

    self%disp(1)%newAtoms = 0E0_dp
    nMoved = 0
    do iAtom = 1, nMax
      if(.not. trialBox%IsActive(iAtom)) then
        cycle
      endif
      self%disp(1)%newAtoms(1:3, iAtom) = trialBox%atoms(1:3, iAtom)
      moveAtom = .true.
      if(self%molType > 0) then
        if(trialBox%MolType(iAtom) /= self%molType) then
          moveAtom = .false.
        endif
      endif
      if(.not. moveAtom) then
        cycle
      endif
      nMoved = nMoved + 1
      self%disp(1)%newAtoms(1, iAtom) = trialBox%atoms(1, iAtom) + &
           lambda * self%forces_old(1, iAtom) + self%sigma * Gaussian()
      self%disp(1)%newAtoms(2, iAtom) = trialBox%atoms(2, iAtom) + &
           lambda * self%forces_old(2, iAtom) + self%sigma * Gaussian()
      self%disp(1)%newAtoms(3, iAtom) = trialBox%atoms(3, iAtom) + &
           lambda * self%forces_old(3, iAtom) + self%sigma * Gaussian()
    enddo
    self%disp(1)%nMoved = nMoved

    if(nMoved < 1) then
      accept = .false.
      return
    endif

    accept = trialBox % CheckConstraint(self%disp(1:1))
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif

    call trialBox%ComputeEnergyDelta(self%disp(1:1), &
                                     self%tempList, &
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      return
    endif

    accept = trialBox % CheckPostEnergy(self%disp(1:1), E_Diff)
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif

    if(trialBox%EFunc%Method%HasComputeForces()) then
      call trialBox%EFunc%Method%ComputeForces(trialBox, self%forces_new, self%disp(1)%newAtoms)
    else
      call trialBox%EFunc%Method%ComputeForcesFiniteDiff(trialBox, self%forces_new, self%disp(1)%newAtoms)
    endif

    ssq_fwd = 0E0_dp
    ssq_rev = 0E0_dp
    do iAtom = 1, nMax
      if(.not. trialBox%IsActive(iAtom)) cycle
      if(self%molType > 0) then
        if(trialBox%MolType(iAtom) /= self%molType) cycle
      endif
      dx = self%disp(1)%newAtoms(1, iAtom) - trialBox%atoms(1, iAtom) - lambda * self%forces_old(1, iAtom)
      dy = self%disp(1)%newAtoms(2, iAtom) - trialBox%atoms(2, iAtom) - lambda * self%forces_old(2, iAtom)
      dz = self%disp(1)%newAtoms(3, iAtom) - trialBox%atoms(3, iAtom) - lambda * self%forces_old(3, iAtom)
      ssq_fwd = ssq_fwd + dx*dx + dy*dy + dz*dz

      dx = trialBox%atoms(1, iAtom) - self%disp(1)%newAtoms(1, iAtom) - lambda * self%forces_new(1, iAtom)
      dy = trialBox%atoms(2, iAtom) - self%disp(1)%newAtoms(2, iAtom) - lambda * self%forces_new(2, iAtom)
      dz = trialBox%atoms(3, iAtom) - self%disp(1)%newAtoms(3, iAtom) - lambda * self%forces_new(3, iAtom)
      ssq_rev = ssq_rev + dx*dx + dy*dy + dz*dz
    enddo

    inv2sig2 = 1.0E0_dp / (2.0E0_dp * self%sigma * self%sigma)
    logProb = (-ssq_rev + ssq_fwd) * inv2sig2

    extraTerms = sampling % GetExtraTerms(self%disp(1:1), trialBox)
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:1), logProb=logProb, extraIn=extraTerms)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList, self%tempNNei)
    else
      self%detailedrej = self%detailedrej + 1
    endif

  end subroutine
!=========================================================================
  subroutine ForceBias_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(ForceBias), intent(inout) :: self

    if(.not. allocated(self%disp(1)%newAtoms)) then
      call self%Constructor
    endif

    if(self%sigma <= 0E0_dp) then
      write(0,*) "ERROR! Force-bias move requires sigma > 0"
      error stop
    endif

    write(nout,"(1x,A,F15.8)") "(Force-Bias) Noise sigma: ", self%sigma
    if(self%autoLambda) then
      write(nout,"(1x,A)") "(Force-Bias) Lambda: 0.5 * sigma^2 * beta (per box)"
    else
      write(nout,"(1x,A,F15.8)") "(Force-Bias) Lambda: ", self%lambda
    endif
    if(self%molType > 0) then
      write(nout,"(1x,A,I8)") "(Force-Bias) Restricted to molecule type: ", self%molType
    else
      write(nout,"(1x,A)") "(Force-Bias) All active atoms are displaced"
    endif

  end subroutine
!=========================================================================
  subroutine ForceBias_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(ForceBias), intent(inout) :: self
    real(dp) :: accptRate

    write(nout,"(1x,A,I15)") "Force-Bias Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Force-Bias Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Force-Bias Acceptance Rate: ", accptRate

    if(self%verbose) then
      write(nout, "(1x,A,I15)") "Force-Bias, Rejections due to overlap:", self%ovlaprej
      write(nout, "(1x,A,I15)") "Force-Bias, Rejections due to constraint:", self%constrainrej
      write(nout, "(1x,A,I15)") "Force-Bias, Rejections due to detailed balance:", self%detailedrej
    endif

  end subroutine
!=========================================================================
  subroutine ForceBias_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(ForceBias), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("sigma")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%sigma = realVal

      case("lambda")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%lambda = realVal
        self%autoLambda = .false.

      case("moltype")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%molType = intVal

      case("verbose")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%verbose = logicVal

      case default
        lineStat = -1
        return
    end select
    lineStat = 0

  end subroutine
!=========================================================================
end module
!=========================================================================
