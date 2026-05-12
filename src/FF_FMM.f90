!================================================================================
! Fast Multipole Method (FMM) Module for Monte Carlo Simulations
!
! This module implements the Fast Multipole Method for computing long-range
! electrostatic interactions in periodic systems. The implementation is
! specifically optimized for Monte Carlo simulations with efficient O(1)
! per-move updates using incremental multipole updates.
!
! The FMM achieves O(N) scaling by:
!   1. Organizing particles into a hierarchical octree
!   2. Representing far-field interactions via multipole expansions
!   3. Using efficient M2M, M2L, L2L translation operators
!   4. Computing near-field interactions directly
!
! For periodic systems, uses lattice summation at the root level.
!
! References:
!   - Greengard & Rokhlin (1987) - A Fast Algorithm for Particle Simulations
!   - Cheng, Greengard, Rokhlin (1999) - Fast Adaptive Multipole Algorithm
!   - Kudin & Scuseria (1998) - Periodic FMM for crystalline systems
!================================================================================
module FF_FMM
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use ClassyConstants, only: coulombConst, pi, two_pi
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  use SphericalHarmonics

  implicit none

  !-----------------------------------------------------------------------------
  ! Octree node type
  !-----------------------------------------------------------------------------
  type :: OctreeNode
    integer :: level                     ! Tree level (0 = root)
    integer :: nodeIndex                 ! Linear index in nodes array
    real(dp) :: center(3)                ! Cell center coordinates
    real(dp) :: halfWidth                ! Half of cell side length
    integer :: parent                    ! Parent node index (-1 for root)
    integer :: children(8)               ! Child node indices (-1 if no child)
    integer :: nParticles                ! Number of particles in this cell
    integer, allocatable :: particles(:) ! Particle indices (only for leaves)
    logical :: isLeaf                    ! True if this is a leaf node
    
    ! Multipole and local expansion coefficients
    complex(dp), allocatable :: M(:)     ! Multipole moments M_l^m
    complex(dp), allocatable :: L(:)     ! Local expansion coefficients L_l^m
    complex(dp), allocatable :: M_old(:) ! Backup for move rejection
    complex(dp), allocatable :: L_old(:) ! Backup for move rejection
    
    ! Interaction lists
    integer, allocatable :: nearList(:)      ! Near-field interaction list
    integer, allocatable :: farList(:)       ! Far-field (M2L) interaction list
    integer :: nNear                         ! Number of near-field cells
    integer :: nFar                          ! Number of far-field cells
  end type OctreeNode

  !-----------------------------------------------------------------------------
  ! Main FMM force field type
  !-----------------------------------------------------------------------------
  type, extends(forcefield) :: Pair_FMM
    ! FMM parameters
    integer :: expansionOrder = 6        ! Multipole expansion order (p)
    integer :: maxLevel = 4              ! Maximum tree depth
    integer :: nCrit = 40                ! Max particles per leaf
    logical :: isPeriodic = .true.       ! Periodic boundary conditions
    
    ! Precision control
    real(dp) :: fmmPrecision = 1.0E-6_dp ! Target precision
    
    ! Charges
    real(dp), allocatable :: q(:)        ! Charge per atom type
    real(dp), allocatable :: qTable(:,:) ! q_i * q_j * coulombConst
    
    ! Octree structure
    type(OctreeNode), allocatable :: nodes(:)  ! All octree nodes
    integer :: nNodes                          ! Current number of nodes
    integer :: maxNodes                        ! Maximum nodes allocated
    integer :: rootNode                        ! Index of root node
    integer, allocatable :: leafNodes(:)       ! Indices of all leaf nodes
    integer :: nLeaves                         ! Number of leaf nodes
    
    ! Particle to node mapping
    integer, allocatable :: particleToLeaf(:)  ! Which leaf contains each particle
    integer, allocatable :: particleToLeaf_old(:) ! Backup for rejection
    
    ! Box dimensions cache
    real(dp) :: boxCenter(3)
    real(dp) :: boxHalfWidth
    real(dp) :: cachedLx, cachedLy, cachedLz
    
    ! Expansion size
    integer :: nExpansion                ! (p+1)^2 coefficients
    
    ! Energy storage
    real(dp) :: E_near_total             ! Near-field energy
    real(dp) :: E_far_total              ! Far-field energy
    real(dp) :: E_self_total             ! Self-energy correction
    real(dp) :: E_periodic_total         ! Periodic correction
    
    ! Minimum distance for overlap detection
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    
    ! Flags
    logical :: initialized = .false.
    logical :: treeBuilt = .false.
    logical :: useMultipole = .false.
    
    ! MC move tracking
    integer :: lastMovedParticle
    integer :: lastOldLeaf
    integer :: lastNewLeaf
    
    ! Pending move tracking for incremental M updates (O(N) far-field)
    logical :: hasPendingMove = .false.
    integer :: pendingAtomIdx = 0
    real(dp) :: pendingOldX = 0.0_dp
    real(dp) :: pendingOldY = 0.0_dp
    real(dp) :: pendingOldZ = 0.0_dp
    integer :: nSavedLeaves = 0
    integer :: savedLeaves(20) = 0
    
  contains
    procedure, pass :: Constructor => Constructor_FMM
    procedure, pass :: DetailedECalc => Detailed_FMM
    procedure, pass :: DiffECalc => DiffECalc_FMM
    procedure, pass :: ShiftECalc_Single => Shift_FMM_Single
    procedure, pass :: NewECalc => New_FMM
    procedure, pass :: OldECalc => Old_FMM
    procedure, pass :: OrthoVolECalc => OrthoVol_FMM
    procedure, pass :: ScaleOctreeForVolume => ScaleOctree_FMM  ! Uniform scale centers/sizes for IsoVol (tree bins)
    procedure, pass :: ManyBody => ManyBody_FMM
    procedure, pass :: DiagnosticMultipoleCompare => DiagnosticMultipoleCompare_FMM
    procedure, pass :: SingleCellM2LTest => SingleCellM2LTest_FMM
    procedure, pass :: ProcessIO => ProcessIO_FMM
    procedure, pass :: Prologue => Prologue_FMM
    procedure, pass :: GetCutOff => GetCutOff_FMM
    procedure, pass :: Epilogue => Epilogue_FMM
    
    ! FMM-specific procedures
    procedure, pass :: BuildOctree => BuildOctree_FMM
    procedure, pass :: BuildInteractionLists => BuildInteractionLists_FMM
    procedure, pass :: ComputeMultipoles => ComputeMultipoles_FMM
    procedure, pass :: ComputeLocalExpansions => ComputeLocalExpansions_FMM
    procedure, pass :: ComputeNearField => ComputeNearField_FMM
    procedure, pass :: ComputeNearFieldOctree => ComputeNearField_Octree_FMM
    procedure, pass :: ComputeFarField => ComputeFarField_FMM
    procedure, pass :: ComputeFarFieldMultipole => ComputeFarFieldMultipole_FMM
    procedure, pass :: ComputeSelfEnergy => ComputeSelfEnergy_FMM
    procedure, pass :: ComputePeriodicCorrection => ComputePeriodicCorrection_FMM
    procedure, pass :: UpdateMultipoleForMove => UpdateMultipoleForMove_FMM
    procedure, pass :: AcceptUpdate => AcceptUpdate_FMM
    procedure, pass :: RejectUpdate => RejectUpdate_FMM
    procedure, pass :: FindLeafForPosition => FindLeafForPosition_FMM
    procedure, pass :: GetParticleEnergy => GetParticleEnergy_FMM
    procedure, pass :: ResolvePendingMove => ResolvePendingMove_FMM
    procedure, pass :: EvalFarFieldPotential => EvalFarFieldPotential_FMM
    procedure, pass :: TentativeMultipoleUpdate => TentativeMultipoleUpdate_FMM
    procedure, pass :: VerifyMultipoles => VerifyMultipoles_FMM
    
  end type Pair_FMM

contains

!=============================================================================
subroutine Constructor_FMM(self)
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_FMM), intent(inout) :: self
  integer :: AllocateStat

  allocate(self%q(1:nAtomTypes), stat=AllocateStat)
  allocate(self%qTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  allocate(self%rMin(1:nAtomTypes), stat=AllocateStat)
  allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)

  self%q = 0.0_dp
  self%qTable = 0.0_dp
  self%rMin = 0.5_dp
  self%rMinTable = 0.25_dp  ! Stored as r^2
  
  ! Default FMM parameters
  self%expansionOrder = 6
  self%maxLevel = 4
  self%nCrit = 40
  self%isPeriodic = .true.
  self%rCut = 10.0_dp   ! Default near-field cutoff (should cover ~3 leaf cells)
  self%rCutSq = self%rCut**2
  
  self%initialized = .false.
  self%treeBuilt = .false.
  self%nNodes = 0
  self%nLeaves = 0
  
  ! Compute expansion size
  self%nExpansion = SH_GetExpansionSize(self%expansionOrder)
  
  ! Initialize spherical harmonics module
  call SH_Init(self%expansionOrder * 2)  ! Need 2p for M2L

  if (AllocateStat /= 0) error stop "Allocation error in FMM Constructor"

end subroutine
!=============================================================================
subroutine Detailed_FMM(self, curbox, E_T, accept)
  use ParallelVar, only: nout
  use Common_MolInfo, only: nMolTypes
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(inout) :: E_T
  logical, intent(out) :: accept
  
  real(dp) :: E_coulomb, E_periodic
  real(dp) :: dipole(3), dipole_sq
  real(dp) :: Lx, Ly, Lz, volume
  real(dp) :: qi, E_near, E_far
  integer :: iAtom, atmType1
  real(dp), pointer :: atoms(:,:) => null()

  accept = .true.
  if (.not. self%isSubPair) curbox%ETable = 0.0_dp
  
  call curbox%GetCoordinates(atoms)
  
  call self%BuildOctree(curbox)
  call self%BuildInteractionLists()
  call self%ComputeMultipoles(curbox)
  call self%ComputeLocalExpansions()
  call self%ComputeNearFieldOctree(curbox, E_near, accept)
  if (.not. accept) return

  call self%ComputeFarField(curbox, E_far)
  write(nout,'(A,I6)') " FMM nLeaves=", self%nLeaves
  E_coulomb = E_near + E_far
  
  call curbox%GetCoordinates(atoms)
  
  ! Get box dimensions and cache them
  select type(curbox)
    class is(OrthoBox)
      Lx = curbox%boxLx
      Ly = curbox%boxLy
      Lz = curbox%boxLz
    class is(CubeBox)
      Lx = curbox%boxL
      Ly = curbox%boxL
      Lz = curbox%boxL
    class default
      Lx = 30.0_dp
      Ly = 30.0_dp
      Lz = 30.0_dp
  end select
  
  ! Cache box dimensions for use in shift calculations
  self%cachedLx = Lx
  self%cachedLy = Ly
  self%cachedLz = Lz
  volume = Lx * Ly * Lz
  
  ! Compute periodic correction if periodic boundaries
  E_periodic = 0.0_dp
  if (self%isPeriodic) then
    ! Compute total dipole moment
    dipole = 0.0_dp
    do iAtom = 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(iAtom)) cycle
      atmType1 = curbox%AtomType(iAtom)
      qi = self%q(atmType1)
      dipole(1) = dipole(1) + qi * atoms(1, iAtom)
      dipole(2) = dipole(2) + qi * atoms(2, iAtom)
      dipole(3) = dipole(3) + qi * atoms(3, iAtom)
    enddo
    
    dipole_sq = dipole(1)**2 + dipole(2)**2 + dipole(3)**2
    
    ! Surface dipole correction (tin-foil boundary conditions)
    E_periodic = 2.0_dp * pi * coulombConst * dipole_sq / (3.0_dp * volume)
  endif
  
  self%E_near_total = E_near
  self%E_far_total = E_far
  self%E_self_total = 0.0_dp
  self%E_periodic_total = E_periodic
  
  write(nout,*) "FMM Near-Field Energy:", E_near
  write(nout,*) "FMM Far-Field Energy:", E_far
  write(nout,*) "FMM Periodic Correction:", E_periodic
  write(nout,*) "FMM Total Energy:", E_coulomb + E_periodic
  
  E_T = E_coulomb + E_periodic

end subroutine
!=============================================================================
subroutine DiffECalc_FMM(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(simBox), intent(inout) :: curbox
  class(Perturbation), intent(inout), target :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  accept = .true.
  if (.not. self%isSubPair) curbox%dETable = 0.0_dp
  E_Diff = 0.0_dp

  select type(disp)
    class is(Displacement)
      call self%ShiftECalc_Single(curbox, disp, E_Diff, accept)
      
    class is(Addition)
      call self%NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)
      
    class is(Deletion)
      call self%OldECalc(curbox, disp, E_Diff)
      
    class is(OrthoVolChange)
      call self%OrthoVolECalc(curbox, disp, E_Diff, accept)
      
    class default
      write(0,*) "Unknown Perturbation Type in FMM DiffECalc"
      error stop
  end select

end subroutine
!=============================================================================
subroutine Shift_FMM_Single(self, curbox, disp, E_Diff, accept)
  ! Computes energy change when atoms are displaced
  ! Uses neighbor list for near-field (short-range) interactions
  ! Also includes the periodic dipole correction change
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jAtom, dispLen, jNei
  integer :: atmType1, atmType2
  integer :: leafIdx, newLeaf, nodeJ, j, k   ! For tree lists in pair delta
  real(dp) :: qi, qj, rmin_ij
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()
  
  ! For periodic correction
  real(dp) :: dipole_old(3), dipole_new(3)
  real(dp) :: dipole_sq_old, dipole_sq_new
  real(dp) :: volume, Lx, Ly, Lz
  real(dp) :: E_periodic_diff
  integer :: kAtom, kType
  real(dp) :: xnew, ynew, znew
  
  call curbox%GetCoordinates(atoms)

  dispLen = size(disp)
  E_Diff = 0.0_dp
  accept = .true.

  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)
    qi = self%q(atmType1)
    
    if (abs(qi) < 1.0E-15_dp) cycle

    leafIdx = self%particleToLeaf(iAtom)
    newLeaf = leafIdx
    if (self%treeBuilt) then
      newLeaf = self%FindLeafForPosition(disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new)
    endif

    if (leafIdx <= 0 .or. newLeaf <= 0 .or. newLeaf /= leafIdx) then
      do jAtom = 1, curbox%nMaxAtoms
        if (jAtom == iAtom) cycle
        if (.not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
        call ComputePairDelta_FMM(self, curbox, atoms, iAtom, qi, jAtom, disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new, E_Diff, accept)
        if (.not. accept) return
      enddo
    else
      ! 1. Own leaf particles (skip self/intra)
      if (allocated(self%nodes(leafIdx)%particles)) then
        do j = 1, size(self%nodes(leafIdx)%particles)
          jAtom = self%nodes(leafIdx)%particles(j)
          if (jAtom == iAtom .or. .not. curbox%IsActive(jAtom) .or. curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
          call ComputePairDelta_FMM(self, curbox, atoms, iAtom, qi, jAtom, disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new, E_Diff, accept)
          if (.not. accept) return
        enddo
      endif

      ! 2. Near cells
      do k = 1, self%nodes(leafIdx)%nNear
        nodeJ = self%nodes(leafIdx)%nearList(k)
        if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
        do j = 1, size(self%nodes(nodeJ)%particles)
          jAtom = self%nodes(nodeJ)%particles(j)
          if (.not. curbox%IsActive(jAtom) .or. curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
          call ComputePairDelta_FMM(self, curbox, atoms, iAtom, qi, jAtom, disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new, E_Diff, accept)
          if (.not. accept) return
        enddo
      enddo

      ! 3. Far-field: direct pair sums over tree far-field cells
      do k = 1, self%nodes(leafIdx)%nFar
        nodeJ = self%nodes(leafIdx)%farList(k)
        if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
        do j = 1, size(self%nodes(nodeJ)%particles)
          jAtom = self%nodes(nodeJ)%particles(j)
          if (.not. curbox%IsActive(jAtom) .or. curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
          call ComputePairDelta_FMM(self, curbox, atoms, iAtom, qi, jAtom, disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new, E_Diff, accept)
          if (.not. accept) return
        enddo
      enddo
    endif
  enddo
  
  ! Compute periodic correction change if periodic boundaries
  if (self%isPeriodic) then
    ! Get box dimensions
    select type(curbox)
      class is(OrthoBox)
        Lx = curbox%boxLx
        Ly = curbox%boxLy
        Lz = curbox%boxLz
      class is(CubeBox)
        Lx = curbox%boxL
        Ly = curbox%boxL
        Lz = curbox%boxL
      class default
        Lx = self%cachedLx
        Ly = self%cachedLy
        Lz = self%cachedLz
    end select
    volume = Lx * Ly * Lz
    
    ! Compute old and new total dipole moments
    dipole_old = 0.0_dp
    dipole_new = 0.0_dp
    
    do kAtom = 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(kAtom)) cycle
      kType = curbox%AtomType(kAtom)
      qj = self%q(kType)
      
      if (abs(qj) < 1.0E-15_dp) cycle
      
      ! Check if this atom is being displaced
      do iDisp = 1, dispLen
        if (disp(iDisp)%atmIndx == kAtom) then
          ! Use old and new positions
          dipole_old(1) = dipole_old(1) + qj * atoms(1, kAtom)
          dipole_old(2) = dipole_old(2) + qj * atoms(2, kAtom)
          dipole_old(3) = dipole_old(3) + qj * atoms(3, kAtom)

          ! IMPORTANT: proposed coordinates may be outside the primary cell.
          ! The box wraps positions on accept (UpdatePosition), so we must wrap
          ! here too to keep the dipole correction consistent with stored coords.
          xnew = disp(iDisp)%x_new
          ynew = disp(iDisp)%y_new
          znew = disp(iDisp)%z_new
          call curbox%Boundary(xnew, ynew, znew)

          dipole_new(1) = dipole_new(1) + qj * xnew
          dipole_new(2) = dipole_new(2) + qj * ynew
          dipole_new(3) = dipole_new(3) + qj * znew
          goto 100
        endif
      enddo
      
      ! Not displaced - same contribution to old and new
      dipole_old(1) = dipole_old(1) + qj * atoms(1, kAtom)
      dipole_old(2) = dipole_old(2) + qj * atoms(2, kAtom)
      dipole_old(3) = dipole_old(3) + qj * atoms(3, kAtom)
      
      dipole_new(1) = dipole_new(1) + qj * atoms(1, kAtom)
      dipole_new(2) = dipole_new(2) + qj * atoms(2, kAtom)
      dipole_new(3) = dipole_new(3) + qj * atoms(3, kAtom)
      
      100 continue
    enddo
    
    dipole_sq_old = dipole_old(1)**2 + dipole_old(2)**2 + dipole_old(3)**2
    dipole_sq_new = dipole_new(1)**2 + dipole_new(2)**2 + dipole_new(3)**2
    
    ! Periodic correction: E = 2*pi/(3*V) * |dipole|^2 * coulombConst
    E_periodic_diff = 2.0_dp * pi * coulombConst * (dipole_sq_new - dipole_sq_old) / (3.0_dp * volume)
    
    E_Diff = E_Diff + E_periodic_diff
  endif
  
  ! Debug output (can be removed after debugging)
  !write(*,*) "SHIFT DEBUG: E_coulomb_diff=", E_new - E_old, " E_periodic_diff=", E_periodic_diff

end subroutine
!=============================================================================
subroutine New_FMM(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  ! Energy change for particle addition
  ! Uses octree structure for near-field and far-field (long-range) interactions
  ! This is the "new position" half of the Shift calculation
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Addition), intent(in) :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jAtom, j, k
  integer :: atmType1, atmType2
  integer :: leafIdx, nodeJ
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_pair, qi, qj
  real(dp), pointer :: atoms(:,:) => null()
  
  ! For periodic correction
  real(dp) :: dipole_old(3), dipole_new(3)
  real(dp) :: dipole_sq_old, dipole_sq_new
  real(dp) :: volume, Lx, Ly, Lz
  real(dp) :: E_periodic_diff
  integer :: kAtom, kType
  real(dp) :: xnew, ynew, znew

  call curbox%GetCoordinates(atoms)

  E_Diff = 0.0_dp
  accept = .true.

  do iDisp = 1, size(disp)
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)
    qi = self%q(atmType1)
    
    if (abs(qi) < 1.0E-15_dp) cycle

    ! Find the leaf that would contain the new position
    leafIdx = self%FindLeafForPosition(disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new)

    ! 1. Own leaf particles
    if (allocated(self%nodes(leafIdx)%particles)) then
      do j = 1, size(self%nodes(leafIdx)%particles)
        jAtom = self%nodes(leafIdx)%particles(j)
        if (jAtom == iAtom .or. .not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle

        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle

        rx = disp(iDisp)%x_new - atoms(1, jAtom)
        ry = disp(iDisp)%y_new - atoms(2, jAtom)
        rz = disp(iDisp)%z_new - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz

        if (rsq < self%rMinTable(atmType1, atmType2)) then
          accept = .false.
          return
        endif

        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_Diff = E_Diff + E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair
      enddo
    endif

    ! 2. Near cells
    do k = 1, self%nodes(leafIdx)%nNear
      nodeJ = self%nodes(leafIdx)%nearList(k)
      if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
      do j = 1, size(self%nodes(nodeJ)%particles)
        jAtom = self%nodes(nodeJ)%particles(j)
        if (.not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle

        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle

        rx = disp(iDisp)%x_new - atoms(1, jAtom)
        ry = disp(iDisp)%y_new - atoms(2, jAtom)
        rz = disp(iDisp)%z_new - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz

        if (rsq < self%rMinTable(atmType1, atmType2)) then
          accept = .false.
          return
        endif

        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_Diff = E_Diff + E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair
      enddo
    enddo

    ! 3. Far cells (long-range interactions)
    do k = 1, self%nodes(leafIdx)%nFar
      nodeJ = self%nodes(leafIdx)%farList(k)
      if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
      do j = 1, size(self%nodes(nodeJ)%particles)
        jAtom = self%nodes(nodeJ)%particles(j)
        if (.not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle

        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle

        rx = disp(iDisp)%x_new - atoms(1, jAtom)
        ry = disp(iDisp)%y_new - atoms(2, jAtom)
        rz = disp(iDisp)%z_new - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz

        if (rsq < self%rMinTable(atmType1, atmType2)) then
          accept = .false.
          return
        endif

        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_Diff = E_Diff + E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair
      enddo
    enddo
  enddo
  
  ! Compute periodic correction change if periodic boundaries
  if (self%isPeriodic) then
    ! Get box dimensions
    select type(curbox)
      class is(OrthoBox)
        Lx = curbox%boxLx
        Ly = curbox%boxLy
        Lz = curbox%boxLz
      class is(CubeBox)
        Lx = curbox%boxL
        Ly = curbox%boxL
        Lz = curbox%boxL
      class default
        Lx = self%cachedLx
        Ly = self%cachedLy
        Lz = self%cachedLz
    end select
    volume = Lx * Ly * Lz
    
    ! Compute old dipole (without new particles)
    dipole_old = 0.0_dp
    do kAtom = 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(kAtom)) cycle
      kType = curbox%AtomType(kAtom)
      qj = self%q(kType)
      if (abs(qj) < 1.0E-15_dp) cycle
      
      ! Skip atoms being added
      do iDisp = 1, size(disp)
        if (disp(iDisp)%atmIndx == kAtom) goto 200
      enddo
      
      dipole_old(1) = dipole_old(1) + qj * atoms(1, kAtom)
      dipole_old(2) = dipole_old(2) + qj * atoms(2, kAtom)
      dipole_old(3) = dipole_old(3) + qj * atoms(3, kAtom)
      200 continue
    enddo
    
    ! New dipole = old + contribution from new particles
    dipole_new = dipole_old
    do iDisp = 1, size(disp)
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox%AtomType(iAtom)
      qi = self%q(atmType1)
      
      ! Proposed insertion coordinates may not be wrapped yet; match UpdatePosition.
      xnew = disp(iDisp)%x_new
      ynew = disp(iDisp)%y_new
      znew = disp(iDisp)%z_new
      call curbox%Boundary(xnew, ynew, znew)

      dipole_new(1) = dipole_new(1) + qi * xnew
      dipole_new(2) = dipole_new(2) + qi * ynew
      dipole_new(3) = dipole_new(3) + qi * znew
    enddo
    
    dipole_sq_old = dipole_old(1)**2 + dipole_old(2)**2 + dipole_old(3)**2
    dipole_sq_new = dipole_new(1)**2 + dipole_new(2)**2 + dipole_new(3)**2
    
    E_periodic_diff = 2.0_dp * pi * coulombConst * (dipole_sq_new - dipole_sq_old) / (3.0_dp * volume)
    E_Diff = E_Diff + E_periodic_diff
  endif

end subroutine
!=============================================================================
subroutine Old_FMM(self, curbox, disp, E_Diff)
  ! Energy change for particle deletion
  ! Uses octree structure for near-field and far-field (long-range) interactions
  ! This is the "old position" half of the Shift calculation (negated)
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Deletion), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff

  integer :: iAtom, jAtom, j, k
  integer :: atmType1, atmType2
  integer :: molEnd, molStart
  integer :: leafIdx, nodeJ
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_pair, qi, qj
  real(dp), pointer :: atoms(:,:) => null()
  
  ! For periodic correction
  real(dp) :: dipole_old(3), dipole_new(3)
  real(dp) :: dipole_sq_old, dipole_sq_new
  real(dp) :: volume, Lx, Ly, Lz
  real(dp) :: E_periodic_diff
  integer :: kAtom, kType

  call curbox%GetCoordinates(atoms)

  E_Diff = 0.0_dp

  call curBox%GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

  ! Compute energy of deleted atoms with all other atoms using octree
  do iAtom = molStart, molEnd
    atmType1 = curbox%AtomType(iAtom)
    qi = self%q(atmType1)
    
    if (abs(qi) < 1.0E-15_dp) cycle

    ! Get leaf for this atom from tree mapping
    leafIdx = self%particleToLeaf(iAtom)

    ! 1. Own leaf particles
    if (allocated(self%nodes(leafIdx)%particles)) then
      do j = 1, size(self%nodes(leafIdx)%particles)
        jAtom = self%nodes(leafIdx)%particles(j)
        if (jAtom == iAtom .or. .not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle

        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle

        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        if (self%isPeriodic) then
          call curbox%Boundary(rx, ry, rz)
        endif
        rsq = rx*rx + ry*ry + rz*rz

        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_Diff = E_Diff - E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair
      enddo
    endif

    ! 2. Near cells
    do k = 1, self%nodes(leafIdx)%nNear
      nodeJ = self%nodes(leafIdx)%nearList(k)
      if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
      do j = 1, size(self%nodes(nodeJ)%particles)
        jAtom = self%nodes(nodeJ)%particles(j)
        if (.not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle

        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle

        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        if (self%isPeriodic) then
          call curbox%Boundary(rx, ry, rz)
        endif
        rsq = rx*rx + ry*ry + rz*rz

        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_Diff = E_Diff - E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair
      enddo
    enddo

    ! 3. Far cells (long-range interactions)
    do k = 1, self%nodes(leafIdx)%nFar
      nodeJ = self%nodes(leafIdx)%farList(k)
      if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
      do j = 1, size(self%nodes(nodeJ)%particles)
        jAtom = self%nodes(nodeJ)%particles(j)
        if (.not. curbox%IsActive(jAtom)) cycle
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle

        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle

        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        if (self%isPeriodic) then
          call curbox%Boundary(rx, ry, rz)
        endif
        rsq = rx*rx + ry*ry + rz*rz

        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_Diff = E_Diff - E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair
      enddo
    enddo
  enddo
  
  ! Compute periodic correction change if periodic boundaries
  if (self%isPeriodic) then
    ! Get box dimensions
    select type(curbox)
      class is(OrthoBox)
        Lx = curbox%boxLx
        Ly = curbox%boxLy
        Lz = curbox%boxLz
      class is(CubeBox)
        Lx = curbox%boxL
        Ly = curbox%boxL
        Lz = curbox%boxL
      class default
        Lx = self%cachedLx
        Ly = self%cachedLy
        Lz = self%cachedLz
    end select
    volume = Lx * Ly * Lz
    
    ! Compute old dipole (with particles being deleted)
    dipole_old = 0.0_dp
    do kAtom = 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(kAtom)) cycle
      kType = curbox%AtomType(kAtom)
      qj = self%q(kType)
      if (abs(qj) < 1.0E-15_dp) cycle
      
      dipole_old(1) = dipole_old(1) + qj * atoms(1, kAtom)
      dipole_old(2) = dipole_old(2) + qj * atoms(2, kAtom)
      dipole_old(3) = dipole_old(3) + qj * atoms(3, kAtom)
    enddo
    
    ! New dipole = old - contribution from deleted particles
    dipole_new = dipole_old
    do iAtom = molStart, molEnd
      atmType1 = curbox%AtomType(iAtom)
      qi = self%q(atmType1)
      
      dipole_new(1) = dipole_new(1) - qi * atoms(1, iAtom)
      dipole_new(2) = dipole_new(2) - qi * atoms(2, iAtom)
      dipole_new(3) = dipole_new(3) - qi * atoms(3, iAtom)
    enddo
    
    dipole_sq_old = dipole_old(1)**2 + dipole_old(2)**2 + dipole_old(3)**2
    dipole_sq_new = dipole_new(1)**2 + dipole_new(2)**2 + dipole_new(3)**2
    
    E_periodic_diff = 2.0_dp * pi * coulombConst * (dipole_sq_new - dipole_sq_old) / (3.0_dp * volume)
    E_Diff = E_Diff + E_periodic_diff
  endif

end subroutine
!=============================================================================
subroutine OrthoVol_FMM(self, curbox, disp, E_Diff, accept)
  ! Energy difference for volume change (NPT ensemble)
  ! For FMM, volume changes require:
  !   1. Particle positions scale with their center of mass
  !   2. All pairwise distances change
  !   3. Octree structure would need rebuilding
  !   4. Periodic correction changes due to volume and dipole changes
  !
  ! Strategy: Compute full energy with scaled positions, subtract old energy
  use Common_MolInfo, only: nMolTypes
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(OrthoVolChange), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iAtom, jAtom
  integer :: iType, typeStart, typeEnd
  integer :: atmType1, atmType2
  integer :: molIndx1, molIndx2
  real(dp) :: iTemp(1:3), jTemp(1:3)
  real(dp) :: dxi, dyi, dzi
  real(dp) :: dxj, dyj, dzj
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: xnew, ynew, znew
  real(dp) :: qi, qj
  real(dp) :: rmin_ij
  real(dp) :: E_Pair, E_New
  
  ! For periodic correction
  real(dp) :: dipole_old(3), dipole_new(3)
  real(dp) :: dipole_sq_old, dipole_sq_new
  real(dp) :: volume_old, volume_new
  real(dp) :: Lx, Ly, Lz
  real(dp) :: E_periodic_old, E_periodic_new

  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)
  
  accept = .true.
  E_New = 0.0_dp
  E_Diff = 0.0_dp
  if (.not. self%isSubPair) curbox%dETable = 0.0_dp

  ! Scale tree boxes uniformly (centers/sizes * scale) for volume change consistency.
  ! (Commented: trial scale would require inverse on reject; Detailed rebuild suffices post-accept for test.)
  ! call self%ScaleOctreeForVolume(disp(1)%xScale, disp(1)%yScale, disp(1)%zScale)

  ! Compute new energy with scaled positions
  ! For volume moves, molecules are translated by scaling their center of mass
  ! Atoms within a molecule maintain their relative positions to the COM
  do iType = 1, nMolTypes
    call curbox%GetTypeAtoms(iType, typeStart, typeEnd)
    if (typeStart < 1) cycle

    do iAtom = typeStart, typeEnd
      if (.not. curbox%IsActive(iAtom)) cycle

      atmType1 = curbox%AtomType(iAtom)
      qi = self%q(atmType1)
      
      if (abs(qi) < 1.0E-15_dp) cycle
      
      molIndx1 = curbox%MolIndx(iAtom)
      
      ! Compute displacement of this atom due to COM scaling
      dxi = curbox%centerMass(1, molIndx1) * (disp(1)%xScale - 1.0_dp)
      dyi = curbox%centerMass(2, molIndx1) * (disp(1)%yScale - 1.0_dp)
      dzi = curbox%centerMass(3, molIndx1) * (disp(1)%zScale - 1.0_dp)
      
      ! Get atom position relative to COM, apply boundary, add back COM
      iTemp(1:3) = atoms(1:3, iAtom) - curbox%centerMass(1:3, molIndx1)
    if (self%isPeriodic) then
      call curbox%Boundary(iTemp(1), iTemp(2), iTemp(3))
    endif
      iTemp(1:3) = iTemp(1:3) + curbox%centerMass(1:3, molIndx1)
      
      ! Loop over higher-indexed atoms (to count each pair once)
      do jAtom = iAtom + 1, curbox%nMaxAtoms
        if (.not. curbox%IsActive(jAtom)) cycle
        
        ! Skip intramolecular pairs
        molIndx2 = curbox%MolIndx(jAtom)
        if (molIndx2 == molIndx1) cycle
        
        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        
        if (abs(qj) < 1.0E-15_dp) cycle
        
        ! Compute displacement of j atom due to COM scaling
        dxj = curbox%centerMass(1, molIndx2) * (disp(1)%xScale - 1.0_dp)
        dyj = curbox%centerMass(2, molIndx2) * (disp(1)%yScale - 1.0_dp)
        dzj = curbox%centerMass(3, molIndx2) * (disp(1)%zScale - 1.0_dp)
        
        ! Get j atom position relative to COM, apply boundary, add back COM
        jTemp(1:3) = atoms(1:3, jAtom) - curbox%centerMass(1:3, molIndx2)
        if (self%isPeriodic) then
          call curbox%Boundary(jTemp(1), jTemp(2), jTemp(3))
        endif
        jTemp(1:3) = jTemp(1:3) + curbox%centerMass(1:3, molIndx2)
        
        ! Compute distance with scaled positions
        rx = iTemp(1) + dxi - jTemp(1) - dxj
        ry = iTemp(2) + dyi - jTemp(2) - dyj
        rz = iTemp(3) + dzi - jTemp(3) - dzj
        
        ! Apply boundary with new box dimensions
        if (self%isPeriodic) then
          call curbox%BoundaryNew(rx, ry, rz, disp)
        endif
        rsq = rx*rx + ry*ry + rz*rz
        
        ! Check for overlap
        rmin_ij = self%rMinTable(atmType1, atmType2)
        if (rsq < rmin_ij) then
          accept = .false.
          return
        endif
        
        ! Compute Coulomb energy (full direct, no rcut - includes far field for consistency with DetailedECalc)
        ! (rmin check already done above; tree not used here as volume changes all positions)
        r = sqrt(rsq)
        E_Pair = qi * qj * coulombConst / r
        E_New = E_New + E_Pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_Pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_Pair
      enddo
    enddo
  enddo

  ! Add periodic correction difference if periodic boundaries
  if (self%isPeriodic) then
    ! Get old box dimensions
    select type(curbox)
      class is(OrthoBox)
        Lx = curbox%boxLx
        Ly = curbox%boxLy
        Lz = curbox%boxLz
      class is(CubeBox)
        Lx = curbox%boxL
        Ly = curbox%boxL
        Lz = curbox%boxL
      class default
        Lx = self%cachedLx
        Ly = self%cachedLy
        Lz = self%cachedLz
    end select
    volume_old = Lx * Ly * Lz
    volume_new = Lx * disp(1)%xScale * Ly * disp(1)%yScale * Lz * disp(1)%zScale
    
    ! Compute old and new dipole moments
    dipole_old = 0.0_dp
    dipole_new = 0.0_dp
    
    do iAtom = 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(iAtom)) cycle
      atmType1 = curbox%AtomType(iAtom)
      qi = self%q(atmType1)
      
      if (abs(qi) < 1.0E-15_dp) cycle
      
      molIndx1 = curbox%MolIndx(iAtom)
      
      ! Old dipole contribution
      dipole_old(1) = dipole_old(1) + qi * atoms(1, iAtom)
      dipole_old(2) = dipole_old(2) + qi * atoms(2, iAtom)
      dipole_old(3) = dipole_old(3) + qi * atoms(3, iAtom)
      
      ! New position: match UpdatePosition (wrap relative to COM, then scale COM)
      dxi = curbox%centerMass(1, molIndx1) * (disp(1)%xScale - 1.0_dp)
      dyi = curbox%centerMass(2, molIndx1) * (disp(1)%yScale - 1.0_dp)
      dzi = curbox%centerMass(3, molIndx1) * (disp(1)%zScale - 1.0_dp)

      iTemp(1:3) = atoms(1:3, iAtom) - curbox%centerMass(1:3, molIndx1)
      if (self%isPeriodic) then
        call curbox%Boundary(iTemp(1), iTemp(2), iTemp(3))
      endif
      iTemp(1:3) = iTemp(1:3) + curbox%centerMass(1:3, molIndx1)

      xnew = iTemp(1) + dxi
      ynew = iTemp(2) + dyi
      znew = iTemp(3) + dzi
      if (self%isPeriodic) then
        call curbox%BoundaryNew(xnew, ynew, znew, disp)
      endif

      dipole_new(1) = dipole_new(1) + qi * xnew
      dipole_new(2) = dipole_new(2) + qi * ynew
      dipole_new(3) = dipole_new(3) + qi * znew
    enddo
    
    dipole_sq_old = dipole_old(1)**2 + dipole_old(2)**2 + dipole_old(3)**2
    dipole_sq_new = dipole_new(1)**2 + dipole_new(2)**2 + dipole_new(3)**2
    
    ! Periodic correction: E = 2*pi/(3*V) * |dipole|^2 * coulombConst
    E_periodic_old = 2.0_dp * pi * coulombConst * dipole_sq_old / (3.0_dp * volume_old)
    E_periodic_new = 2.0_dp * pi * coulombConst * dipole_sq_new / (3.0_dp * volume_new)
    
    E_New = E_New + E_periodic_new
  endif

  ! Only subtract E_Inter if this is the main forcefield (not a sub-pair in a hybrid)
  ! For sub-pairs, the hybrid handles the E_Inter subtraction
  if (.not. self%isSubPair) then
    E_Diff = E_New - curbox%E_Inter
    curbox%dETable = curbox%dETable - curbox%ETable
  else
    E_Diff = E_New
  endif

end subroutine
!=============================================================================
subroutine ScaleOctree_FMM(self, scaleX, scaleY, scaleZ)
  ! Uniformly scale octree node centers and sizes by volume scale factors (for IsoVol).
  ! Interaction lists remain valid since relative separations preserved.
  ! (Used to keep tree consistent with scaled positions after volume accept; rebuild in Detailed also works.)
  ! Note: If multipole active, M/L would require recompute here.
  implicit none
  class(Pair_FMM), intent(inout) :: self
  real(dp), intent(in) :: scaleX, scaleY, scaleZ
  integer :: iNode

  do iNode = 1, self%nNodes
    self%nodes(iNode)%center(1) = self%nodes(iNode)%center(1) * scaleX
    self%nodes(iNode)%center(2) = self%nodes(iNode)%center(2) * scaleY
    self%nodes(iNode)%center(3) = self%nodes(iNode)%center(3) * scaleZ
    self%nodes(iNode)%halfWidth = self%nodes(iNode)%halfWidth * scaleX  ! Isotropic for IsoVol; lists valid
  enddo
end subroutine
!=============================================================================
subroutine ManyBody_FMM(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept)
  ! For trial insertion moves (CBMC, etc.)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(simBox), intent(inout) :: curbox
  integer, intent(in) :: atmtype1
  integer, intent(in) :: atmtypes(:)
  real(dp), intent(in) :: pos1(:)
  real(dp), intent(in) :: posN(:,:)
  logical, intent(out) :: accept
  real(dp), intent(out) :: E_Many

  integer :: jAtom, atmType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_pair, qi, qj

  accept = .true.
  E_Many = 0.0_dp
  
  qi = self%q(atmtype1)
  if (abs(qi) < 1.0E-15_dp) return
  
  ! Only near-field contribution for quick trial evaluation
  do jAtom = 1, size(posN, 2)
    atmType2 = atmtypes(jAtom)
    qj = self%q(atmType2)
    
    rx = pos1(1) - posN(1, jAtom)
    ry = pos1(2) - posN(2, jAtom)
    rz = pos1(3) - posN(3, jAtom)
    if (self%isPeriodic) then
      call curbox%Boundary(rx, ry, rz)
    endif
    rsq = rx*rx + ry*ry + rz*rz
    
    if (rsq < self%rMinTable(atmtype1, atmType2)) then
      accept = .false.
      return
    endif
    
    if (rsq < self%rCutSq .and. rsq > 1.0E-10_dp) then
      r = sqrt(rsq)
      E_pair = qi * qj * coulombConst / r
      E_Many = E_Many + E_pair
    endif
  enddo

end subroutine
!=============================================================================
subroutine BuildOctree_FMM(self, curbox)
  ! Build the octree from particle positions
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  real(dp) :: Lx, Ly, Lz
  real(dp), pointer :: atoms(:,:) => null()
  integer :: iAtom, nActive
  integer :: maxNodesNeeded
  
  call curbox%GetCoordinates(atoms)
  
  ! Get box dimensions
  select type(curbox)
    class is(OrthoBox)
      Lx = curbox%boxLx
      Ly = curbox%boxLy
      Lz = curbox%boxLz
    class is(CubeBox)
      Lx = curbox%boxL
      Ly = curbox%boxL
      Lz = curbox%boxL
    class default
      Lx = 30.0_dp
      Ly = 30.0_dp
      Lz = 30.0_dp
  end select
  
  ! Cache box dimensions
  self%cachedLx = Lx
  self%cachedLy = Ly
  self%cachedLz = Lz
  
  ! Box center and half-width
  !self%boxCenter = [Lx/2.0_dp, Ly/2.0_dp, Lz/2.0_dp]
  self%boxCenter = [0.0_dp, 0.0_dp, 0.0_dp]
  self%boxHalfWidth = max(Lx, Ly, Lz) / 2.0_dp
  
  ! Count active atoms
  nActive = 0
  do iAtom = 1, curbox%nMaxAtoms
    if (curbox%IsActive(iAtom)) nActive = nActive + 1
  enddo
  
  ! Estimate maximum nodes needed (8^maxLevel for full tree)
  maxNodesNeeded = min(8**self%maxLevel, 8 * nActive / self%nCrit + 100)
  self%maxNodes = maxNodesNeeded
  
  ! Allocate or reallocate nodes array
  if (allocated(self%nodes)) then
    if (size(self%nodes) < maxNodesNeeded) then
      deallocate(self%nodes)
      allocate(self%nodes(maxNodesNeeded))
    endif
  else
    allocate(self%nodes(maxNodesNeeded))
  endif
  
  ! Allocate particle to leaf mapping
  if (allocated(self%particleToLeaf)) deallocate(self%particleToLeaf)
  if (allocated(self%particleToLeaf_old)) deallocate(self%particleToLeaf_old)
  allocate(self%particleToLeaf(curbox%nMaxAtoms))
  allocate(self%particleToLeaf_old(curbox%nMaxAtoms))
  self%particleToLeaf = -1
  
  ! Allocate leaf nodes array
  if (allocated(self%leafNodes)) deallocate(self%leafNodes)
  allocate(self%leafNodes(maxNodesNeeded))
  
  ! Create root node
  self%nNodes = 1
  self%rootNode = 1
  call InitNode(self%nodes(1), 1, 0, self%boxCenter, self%boxHalfWidth, &
                -1, self%nExpansion)
  
  ! Insert all particles
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    call InsertParticle(self, iAtom, atoms(1:3, iAtom), 1)
  enddo
  
  ! Collect leaf nodes
  self%nLeaves = 0
  call CollectLeaves(self, 1)
  
  self%treeBuilt = .true.

contains

  subroutine InitNode(node, idx, level, center, halfWidth, parent, nExp)
    type(OctreeNode), intent(inout) :: node
    integer, intent(in) :: idx, level, parent, nExp
    real(dp), intent(in) :: center(3), halfWidth
    
    node%level = level
    node%nodeIndex = idx
    node%center = center
    node%halfWidth = halfWidth
    node%parent = parent
    node%children = -1
    node%nParticles = 0
    node%isLeaf = .true.
    node%nNear = 0
    node%nFar = 0
    
    if (allocated(node%particles)) deallocate(node%particles)
    if (allocated(node%M)) deallocate(node%M)
    if (allocated(node%L)) deallocate(node%L)
    if (allocated(node%M_old)) deallocate(node%M_old)
    if (allocated(node%L_old)) deallocate(node%L_old)
    if (allocated(node%nearList)) deallocate(node%nearList)
    if (allocated(node%farList)) deallocate(node%farList)
    
    allocate(node%M(nExp))
    allocate(node%L(nExp))
    allocate(node%M_old(nExp))
    allocate(node%L_old(nExp))
    node%M = (0.0_dp, 0.0_dp)
    node%L = (0.0_dp, 0.0_dp)
    
  end subroutine
  
  recursive subroutine InsertParticle(self, iAtom, pos, nodeIdx)
    class(Pair_FMM), intent(inout) :: self
    integer, intent(in) :: iAtom, nodeIdx
    real(dp), intent(in) :: pos(3)
    
    integer :: childIdx, octant
    real(dp) :: childCenter(3), childHalf
    
    self%nodes(nodeIdx)%nParticles = self%nodes(nodeIdx)%nParticles + 1
    
    if (self%nodes(nodeIdx)%isLeaf) then
      ! If this is a leaf, either add particle or subdivide
      if (self%nodes(nodeIdx)%nParticles <= self%nCrit .or. &
          self%nodes(nodeIdx)%level >= self%maxLevel) then
        ! Add to leaf
        call AddParticleToLeaf(self%nodes(nodeIdx), iAtom)
        self%particleToLeaf(iAtom) = nodeIdx
      else
        ! Subdivide and redistribute existing particles
        call SubdivideNode(self, nodeIdx)
        ! Now insert this particle - node is no longer a leaf, so it will go to else branch
        octant = GetOctant(pos, self%nodes(nodeIdx)%center)
        ! Create child if it doesn't exist (may not have been created during subdivision)
        if (self%nodes(nodeIdx)%children(octant) < 0) then
          self%nNodes = self%nNodes + 1
          childIdx = self%nNodes
          self%nodes(nodeIdx)%children(octant) = childIdx
          childHalf = self%nodes(nodeIdx)%halfWidth / 2.0_dp
          childCenter = GetChildCenter(self%nodes(nodeIdx)%center, octant, childHalf)
          call InitNode(self%nodes(childIdx), childIdx, self%nodes(nodeIdx)%level + 1, &
                        childCenter, childHalf, nodeIdx, self%nExpansion)
        endif
        call InsertParticle(self, iAtom, pos, self%nodes(nodeIdx)%children(octant))
      endif
    else
      ! Insert into appropriate child
      octant = GetOctant(pos, self%nodes(nodeIdx)%center)
      if (self%nodes(nodeIdx)%children(octant) < 0) then
        ! Create child
        self%nNodes = self%nNodes + 1
        childIdx = self%nNodes
        self%nodes(nodeIdx)%children(octant) = childIdx
        childHalf = self%nodes(nodeIdx)%halfWidth / 2.0_dp
        childCenter = GetChildCenter(self%nodes(nodeIdx)%center, octant, childHalf)
        call InitNode(self%nodes(childIdx), childIdx, self%nodes(nodeIdx)%level + 1, &
                      childCenter, childHalf, nodeIdx, self%nExpansion)
      endif
      call InsertParticle(self, iAtom, pos, self%nodes(nodeIdx)%children(octant))
    endif
    
  end subroutine
  
  subroutine SubdivideNode(self, nodeIdx)
    class(Pair_FMM), intent(inout) :: self
    integer, intent(in) :: nodeIdx
    
    integer :: i, iAtom, octant, childIdx, nPart
    real(dp) :: childHalf, childCenter(3), pos(3)
    real(dp), pointer :: atoms(:,:)
    integer, allocatable :: tempParticles(:)
    
    call curbox%GetCoordinates(atoms)
    self%nodes(nodeIdx)%isLeaf = .false.
    childHalf = self%nodes(nodeIdx)%halfWidth / 2.0_dp
    
    ! Redistribute existing particles to children
    if (allocated(self%nodes(nodeIdx)%particles)) then
      nPart = size(self%nodes(nodeIdx)%particles)
      allocate(tempParticles(nPart))
      tempParticles = self%nodes(nodeIdx)%particles
      deallocate(self%nodes(nodeIdx)%particles)
      
      do i = 1, nPart
        iAtom = tempParticles(i)
        pos = atoms(1:3, iAtom)
        octant = GetOctant(pos, self%nodes(nodeIdx)%center)
        
        if (self%nodes(nodeIdx)%children(octant) < 0) then
          self%nNodes = self%nNodes + 1
          childIdx = self%nNodes
          self%nodes(nodeIdx)%children(octant) = childIdx
          childCenter = GetChildCenter(self%nodes(nodeIdx)%center, octant, childHalf)
          call InitNode(self%nodes(childIdx), childIdx, self%nodes(nodeIdx)%level + 1, &
                        childCenter, childHalf, nodeIdx, self%nExpansion)
        endif
        
        childIdx = self%nodes(nodeIdx)%children(octant)
        self%nodes(childIdx)%nParticles = self%nodes(childIdx)%nParticles + 1
        call AddParticleToLeaf(self%nodes(childIdx), iAtom)
        self%particleToLeaf(iAtom) = childIdx
      enddo
      
      deallocate(tempParticles)
    endif
    self%nodes(nodeIdx)%nParticles = 0
    
  end subroutine
  
  subroutine AddParticleToLeaf(node, iAtom)
    type(OctreeNode), intent(inout) :: node
    integer, intent(in) :: iAtom
    
    integer, allocatable :: temp(:)
    integer :: n
    
    if (.not. allocated(node%particles)) then
      allocate(node%particles(1))
      node%particles(1) = iAtom
    else
      n = size(node%particles)
      allocate(temp(n+1))
      temp(1:n) = node%particles
      temp(n+1) = iAtom
      call move_alloc(temp, node%particles)
    endif
    
  end subroutine
  
  pure function GetOctant(pos, center) result(octant)
    real(dp), intent(in) :: pos(3), center(3)
    integer :: octant
    
    octant = 1
    if (pos(1) >= center(1)) octant = octant + 1
    if (pos(2) >= center(2)) octant = octant + 2
    if (pos(3) >= center(3)) octant = octant + 4
    
  end function
  
  pure function GetChildCenter(parentCenter, octant, childHalf) result(center)
    real(dp), intent(in) :: parentCenter(3), childHalf
    integer, intent(in) :: octant
    real(dp) :: center(3)
    
    center = parentCenter
    if (mod(octant-1, 2) == 0) then
      center(1) = center(1) - childHalf
    else
      center(1) = center(1) + childHalf
    endif
    if (mod((octant-1)/2, 2) == 0) then
      center(2) = center(2) - childHalf
    else
      center(2) = center(2) + childHalf
    endif
    if (mod((octant-1)/4, 2) == 0) then
      center(3) = center(3) - childHalf
    else
      center(3) = center(3) + childHalf
    endif
    
  end function
  
  recursive subroutine CollectLeaves(self, nodeIdx)
    class(Pair_FMM), intent(inout) :: self
    integer, intent(in) :: nodeIdx
    
    integer :: i, childIdx
    
    if (self%nodes(nodeIdx)%isLeaf) then
      self%nLeaves = self%nLeaves + 1
      self%leafNodes(self%nLeaves) = nodeIdx
    else
      do i = 1, 8
        childIdx = self%nodes(nodeIdx)%children(i)
        if (childIdx > 0) call CollectLeaves(self, childIdx)
      enddo
    endif
    
  end subroutine
  
end subroutine
!=============================================================================
subroutine BuildInteractionLists_FMM(self)
  ! Build near-field and far-field interaction lists for all leaves
  implicit none
  class(Pair_FMM), intent(inout) :: self

  integer :: iLeaf, jLeaf, nodeI, nodeJ
  real(dp) :: dx, dy, dz, dist
  real(dp) :: threshold
  integer :: nearCount, farCount
  integer, allocatable :: tempNear(:), tempFar(:)
  
  allocate(tempNear(self%nLeaves))
  allocate(tempFar(self%nLeaves))
  
  ! For each leaf, determine which other leaves are near or far
  do iLeaf = 1, self%nLeaves
    nodeI = self%leafNodes(iLeaf)
    
    nearCount = 0
    farCount = 0
    
    do jLeaf = 1, self%nLeaves
      if (iLeaf == jLeaf) cycle
      
      nodeJ = self%leafNodes(jLeaf)
      
      ! Well-separation criterion must account for BOTH cell sizes
      threshold = 3.0_dp * (self%nodes(nodeI)%halfWidth + self%nodes(nodeJ)%halfWidth)
      
      dx = self%nodes(nodeI)%center(1) - self%nodes(nodeJ)%center(1)
      dy = self%nodes(nodeI)%center(2) - self%nodes(nodeJ)%center(2)
      dz = self%nodes(nodeI)%center(3) - self%nodes(nodeJ)%center(3)
      
      ! Apply periodic boundaries
      if (self%isPeriodic) then
        if (dx > self%cachedLx/2.0_dp) dx = dx - self%cachedLx
        if (dx < -self%cachedLx/2.0_dp) dx = dx + self%cachedLx
        if (dy > self%cachedLy/2.0_dp) dy = dy - self%cachedLy
        if (dy < -self%cachedLy/2.0_dp) dy = dy + self%cachedLy
        if (dz > self%cachedLz/2.0_dp) dz = dz - self%cachedLz
        if (dz < -self%cachedLz/2.0_dp) dz = dz + self%cachedLz
      endif
      
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      
      if (dist < threshold) then
        nearCount = nearCount + 1
        tempNear(nearCount) = nodeJ
      else
        ! For periodic systems, ensure the expansion displacement (R - r) 
        ! stays within the minimum image range for ALL source-target pairs.
        ! max|R_a - r_a| = |d_a| + hw_source + hw_target per axis.
        block
          logical :: periodic_ok
          real(dp) :: hwI, hwJ
          periodic_ok = .true.
          if (self%isPeriodic) then
            hwI = self%nodes(nodeI)%halfWidth
            hwJ = self%nodes(nodeJ)%halfWidth
            if (abs(dx) + hwI + hwJ >= self%cachedLx/2.0_dp) periodic_ok = .false.
            if (abs(dy) + hwI + hwJ >= self%cachedLy/2.0_dp) periodic_ok = .false.
            if (abs(dz) + hwI + hwJ >= self%cachedLz/2.0_dp) periodic_ok = .false.
          endif
          if (periodic_ok) then
            farCount = farCount + 1
            tempFar(farCount) = nodeJ
          else
            nearCount = nearCount + 1
            tempNear(nearCount) = nodeJ
          endif
        end block
      endif
    enddo
    
    ! Store interaction lists
    self%nodes(nodeI)%nNear = nearCount
    self%nodes(nodeI)%nFar = farCount
    
    if (allocated(self%nodes(nodeI)%nearList)) deallocate(self%nodes(nodeI)%nearList)
    if (allocated(self%nodes(nodeI)%farList)) deallocate(self%nodes(nodeI)%farList)
    
    if (nearCount > 0) then
      allocate(self%nodes(nodeI)%nearList(nearCount))
      self%nodes(nodeI)%nearList = tempNear(1:nearCount)
    endif
    
    if (farCount > 0) then
      allocate(self%nodes(nodeI)%farList(farCount))
      self%nodes(nodeI)%farList = tempFar(1:farCount)
    endif
  enddo
  
  deallocate(tempNear, tempFar)

end subroutine
!=============================================================================
subroutine ComputeMultipoles_FMM(self, curbox)
  ! Compute multipole moments for all leaves and propagate up (M2M)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  integer :: iLeaf, nodeIdx, iAtom, j, l, m
  integer :: atmType
  real(dp) :: qi, dx, dy, dz
  real(dp), pointer :: atoms(:,:) => null()
  complex(dp) :: Rlm(self%nExpansion)
  
  call curbox%GetCoordinates(atoms)
  
  ! Clear all multipoles
  do nodeIdx = 1, self%nNodes
    self%nodes(nodeIdx)%M = (0.0_dp, 0.0_dp)
    self%nodes(nodeIdx)%L = (0.0_dp, 0.0_dp)
  enddo
  
  ! Compute multipole moments for each leaf (P2M)
  do iLeaf = 1, self%nLeaves
    nodeIdx = self%leafNodes(iLeaf)
    
    if (.not. allocated(self%nodes(nodeIdx)%particles)) cycle
    
    do j = 1, size(self%nodes(nodeIdx)%particles)
      iAtom = self%nodes(nodeIdx)%particles(j)
      atmType = curbox%AtomType(iAtom)
      qi = self%q(atmType)
      
      if (abs(qi) < 1.0E-15_dp) cycle
      
      ! Position relative to cell center
      dx = atoms(1, iAtom) - self%nodes(nodeIdx)%center(1)
      dy = atoms(2, iAtom) - self%nodes(nodeIdx)%center(2)
      dz = atoms(3, iAtom) - self%nodes(nodeIdx)%center(3)
      
      ! Compute regular solid harmonics at particle position relative to cell center
      call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm)
      
      ! P2M: M_l^m = sum_i q_i * R_l^{-m}(y_i - c)
      do l = 0, self%expansionOrder
        do m = -l, l
          self%nodes(nodeIdx)%M(SH_Index(l,m)) = self%nodes(nodeIdx)%M(SH_Index(l,m)) &
            + qi * Rlm(SH_Index(l,-m))
        enddo
      enddo
    enddo
  enddo
  
  ! M2M: propagate multipoles up the tree
  call M2M_Upward(self, self%rootNode)

contains

  recursive subroutine M2M_Upward(self, nodeIdx)
    class(Pair_FMM), intent(inout) :: self
    integer, intent(in) :: nodeIdx
    
    integer :: i, childIdx
    real(dp) :: dx, dy, dz
    
    if (self%nodes(nodeIdx)%isLeaf) return
    
    ! First, process all children
    do i = 1, 8
      childIdx = self%nodes(nodeIdx)%children(i)
      if (childIdx > 0) call M2M_Upward(self, childIdx)
    enddo
    
    ! Then, translate child multipoles to this node
    do i = 1, 8
      childIdx = self%nodes(nodeIdx)%children(i)
      if (childIdx <= 0) cycle
      
      ! Translation vector from child to parent center
      dx = self%nodes(childIdx)%center(1) - self%nodes(nodeIdx)%center(1)
      dy = self%nodes(childIdx)%center(2) - self%nodes(nodeIdx)%center(2)
      dz = self%nodes(childIdx)%center(3) - self%nodes(nodeIdx)%center(3)
      
      call SH_M2M_Coeff(dx, dy, dz, self%expansionOrder, &
                        self%nodes(childIdx)%M, self%nodes(nodeIdx)%M)
    enddo
    
  end subroutine
  
end subroutine
!=============================================================================
subroutine ComputeLocalExpansions_FMM(self)
  ! Compute local expansions using M2L at leaf level and propagate down via L2L.
  !
  ! Full FMM pipeline:
  !   1. M2L: Convert far-field multipoles to local expansions at each leaf
  !   2. L2L: Propagate local expansions from parent to children (downward pass)
  !
  ! Currently uses flat leaf-only interaction lists. The M2L is done at leaf level
  ! and L2L propagates any parent-level local contributions down to children.
  implicit none
  class(Pair_FMM), intent(inout) :: self

  integer :: iLeaf, nodeI, k, nodeJ
  real(dp) :: dx, dy, dz
  
  ! Clear all local expansions before computing
  do nodeI = 1, self%nNodes
    self%nodes(nodeI)%L = (0.0_dp, 0.0_dp)
  enddo
  
  ! M2L: for each leaf, convert far-field multipoles to local
  do iLeaf = 1, self%nLeaves
    nodeI = self%leafNodes(iLeaf)
    
    do k = 1, self%nodes(nodeI)%nFar
      nodeJ = self%nodes(nodeI)%farList(k)
      
      ! Translation vector from source to target (d = target - source)
      dx = self%nodes(nodeI)%center(1) - self%nodes(nodeJ)%center(1)
      dy = self%nodes(nodeI)%center(2) - self%nodes(nodeJ)%center(2)
      dz = self%nodes(nodeI)%center(3) - self%nodes(nodeJ)%center(3)
      
      ! Apply periodic boundaries
      if (self%isPeriodic) then
        if (dx > self%cachedLx/2.0_dp) dx = dx - self%cachedLx
        if (dx < -self%cachedLx/2.0_dp) dx = dx + self%cachedLx
        if (dy > self%cachedLy/2.0_dp) dy = dy - self%cachedLy
        if (dy < -self%cachedLy/2.0_dp) dy = dy + self%cachedLy
        if (dz > self%cachedLz/2.0_dp) dz = dz - self%cachedLz
        if (dz < -self%cachedLz/2.0_dp) dz = dz + self%cachedLz
      endif
      
      call SH_M2L_Coeff(dx, dy, dz, self%expansionOrder, &
                        self%nodes(nodeJ)%M, self%nodes(nodeI)%L)
    enddo
  enddo
  
  ! L2L downward pass: propagate local expansions from parent to children
  ! This is needed when M2L is done at higher levels (hierarchical FMM).
  ! For the current flat leaf-only lists, M2L is already at leaf level,
  ! but we still do L2L in case any parent nodes accumulated L contributions.
  call L2L_Downward(self, self%rootNode)

contains

  recursive subroutine L2L_Downward(self, nodeIdx)
    class(Pair_FMM), intent(inout) :: self
    integer, intent(in) :: nodeIdx
    
    integer :: i, childIdx
    real(dp) :: dx, dy, dz
    
    if (self%nodes(nodeIdx)%isLeaf) return
    
    ! Translate this node's local expansion to each child
    do i = 1, 8
      childIdx = self%nodes(nodeIdx)%children(i)
      if (childIdx <= 0) cycle
      
      ! Translation vector: tau = child_center - parent_center
      dx = self%nodes(childIdx)%center(1) - self%nodes(nodeIdx)%center(1)
      dy = self%nodes(childIdx)%center(2) - self%nodes(nodeIdx)%center(2)
      dz = self%nodes(childIdx)%center(3) - self%nodes(nodeIdx)%center(3)
      
      call SH_L2L_Coeff(dx, dy, dz, self%expansionOrder, &
                        self%nodes(nodeIdx)%L, self%nodes(childIdx)%L)
    enddo
    
    ! Then recurse into children
    do i = 1, 8
      childIdx = self%nodes(nodeIdx)%children(i)
      if (childIdx > 0) call L2L_Downward(self, childIdx)
    enddo
    
  end subroutine

end subroutine
!=============================================================================
subroutine ComputeNearField_FMM(self, curbox, E_near, accept)
  ! Compute near-field energy by looping over all atom pairs
  ! (like Ewald does for detailed calculation)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_near
  logical, intent(out) :: accept

  integer :: iAtom, jAtom
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: qi, qj, E_pair
  real(dp), pointer :: atoms(:,:) => null()
  
  call curbox%GetCoordinates(atoms)
  
  E_near = 0.0_dp
  accept = .true.
  
  ! Loop over all unique atom pairs (like Ewald's detailed calculation)
  do iAtom = 1, curbox%nMaxAtoms - 1
    if (.not. curbox%IsActive(iAtom)) cycle
    
    atmType1 = curbox%AtomType(iAtom)
    qi = self%q(atmType1)
    if (abs(qi) < 1.0E-15_dp) cycle
    
    do jAtom = iAtom + 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(jAtom)) cycle
      
      ! Skip intramolecular pairs
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      
      atmType2 = curbox%AtomType(jAtom)
      qj = self%q(atmType2)
      if (abs(qj) < 1.0E-15_dp) cycle
      
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      if (self%isPeriodic) then
        call curbox%Boundary(rx, ry, rz)
      endif
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rMinTable(atmType1, atmType2)) then
        accept = .false.
        return
      endif
      
      ! Compute Coulomb interaction (no cutoff for full electrostatics)
      r = sqrt(rsq)
      E_pair = qi * qj * coulombConst / r
      E_near = E_near + E_pair
      curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_pair
      curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_pair
    enddo
  enddo

end subroutine
!=============================================================================
subroutine ComputeNearField_Octree_FMM(self, curbox, E_near, accept)
  ! Compute near-field energy using octree structure
  ! Near-field = direct Coulomb sum for particles in:
  !   1. Same leaf cell (self-interactions)
  !   2. Neighboring leaf cells (nearList)
  ! No distance cutoff is used - all pairs in near cells are computed exactly
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_near
  logical, intent(out) :: accept

  integer :: iLeaf, nodeI, k, nodeJ
  integer :: i, j, iAtom, jAtom
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: qi, qj, E_pair
  real(dp), pointer :: atoms(:,:) => null()
  
  call curbox%GetCoordinates(atoms)
  
  E_near = 0.0_dp
  accept = .true.
  
  ! Self-interactions within each leaf (all pairs in same cell)
  do iLeaf = 1, self%nLeaves
    nodeI = self%leafNodes(iLeaf)
    
    if (.not. allocated(self%nodes(nodeI)%particles)) cycle
    if (size(self%nodes(nodeI)%particles) < 2) cycle
    
    do i = 1, size(self%nodes(nodeI)%particles) - 1
      iAtom = self%nodes(nodeI)%particles(i)
      atmType1 = curbox%AtomType(iAtom)
      qi = self%q(atmType1)
      
      if (abs(qi) < 1.0E-15_dp) cycle
      
      do j = i + 1, size(self%nodes(nodeI)%particles)
        jAtom = self%nodes(nodeI)%particles(j)
        
        ! Skip intramolecular pairs
        if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
        
        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        
        if (abs(qj) < 1.0E-15_dp) cycle
        
        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        if (self%isPeriodic) then
          call curbox%Boundary(rx, ry, rz)
        endif
        rsq = rx*rx + ry*ry + rz*rz
        
        if (rsq < self%rMinTable(atmType1, atmType2)) then
          accept = .false.
          return
        endif
        
        ! Direct Coulomb - no cutoff within near-field cells
        r = sqrt(rsq)
        E_pair = qi * qj * coulombConst / r
        E_near = E_near + E_pair
        curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_pair
        curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_pair
      enddo
    enddo
  enddo
  
  ! Interactions between neighboring leaves (cells in nearList)
  do iLeaf = 1, self%nLeaves
    nodeI = self%leafNodes(iLeaf)
    
    if (.not. allocated(self%nodes(nodeI)%particles)) cycle
    
    do k = 1, self%nodes(nodeI)%nNear
      nodeJ = self%nodes(nodeI)%nearList(k)
      
      ! Only count each pair once (nodeI < nodeJ)
      if (nodeI >= nodeJ) cycle
      
      if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
      
      do i = 1, size(self%nodes(nodeI)%particles)
        iAtom = self%nodes(nodeI)%particles(i)
        atmType1 = curbox%AtomType(iAtom)
        qi = self%q(atmType1)
        
        if (abs(qi) < 1.0E-15_dp) cycle
        
        do j = 1, size(self%nodes(nodeJ)%particles)
          jAtom = self%nodes(nodeJ)%particles(j)
          
          if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
          
          atmType2 = curbox%AtomType(jAtom)
          qj = self%q(atmType2)
          
          if (abs(qj) < 1.0E-15_dp) cycle
          
          rx = atoms(1, iAtom) - atoms(1, jAtom)
          ry = atoms(2, iAtom) - atoms(2, jAtom)
          rz = atoms(3, iAtom) - atoms(3, jAtom)
          if (self%isPeriodic) then
            call curbox%Boundary(rx, ry, rz)
          endif
          rsq = rx*rx + ry*ry + rz*rz
          
          if (rsq < self%rMinTable(atmType1, atmType2)) then
            accept = .false.
            return
          endif
          
          ! Direct Coulomb - no cutoff within near-field cells
          r = sqrt(rsq)
          E_pair = qi * qj * coulombConst / r
          E_near = E_near + E_pair
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_pair
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_pair
        enddo
      enddo
    enddo
  enddo

end subroutine
!=============================================================================
subroutine ComputeFarField_FMM(self, curbox, E_far)
  ! Compute far-field energy using direct Coulomb summation over far cells
  !
  ! This computes exact interactions between particles in cells that are in
  ! each other's far-field lists. This is O(N^2) in the worst case but ensures
  ! correctness while we debug the multipole approximation.
  !
  ! The multipole-based far-field (ComputeFarFieldMultipole) can be used once
  ! the translation operators are verified to be correct.
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_far

  integer :: iLeaf, nodeI, k, nodeJ
  integer :: i, j, iAtom, jAtom
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: qi, qj, E_pair
  real(dp), pointer :: atoms(:,:) => null()
  
  call curbox%GetCoordinates(atoms)
  
  E_far = 0.0_dp
  
  ! Loop over all leaves and their far-field cells
  do iLeaf = 1, self%nLeaves
    nodeI = self%leafNodes(iLeaf)
    
    if (.not. allocated(self%nodes(nodeI)%particles)) cycle
    
    do k = 1, self%nodes(nodeI)%nFar
      nodeJ = self%nodes(nodeI)%farList(k)
      
      ! Only count each cell pair once (nodeI < nodeJ)
      if (nodeI >= nodeJ) cycle
      
      if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
      
      ! Compute all pairwise interactions between particles in these cells
      do i = 1, size(self%nodes(nodeI)%particles)
        iAtom = self%nodes(nodeI)%particles(i)
        atmType1 = curbox%AtomType(iAtom)
        qi = self%q(atmType1)
        
        if (abs(qi) < 1.0E-15_dp) cycle
        
        do j = 1, size(self%nodes(nodeJ)%particles)
          jAtom = self%nodes(nodeJ)%particles(j)
          
          ! Skip intramolecular pairs
          if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
          
          atmType2 = curbox%AtomType(jAtom)
          qj = self%q(atmType2)
          
          if (abs(qj) < 1.0E-15_dp) cycle
          
          rx = atoms(1, iAtom) - atoms(1, jAtom)
          ry = atoms(2, iAtom) - atoms(2, jAtom)
          rz = atoms(3, iAtom) - atoms(3, jAtom)
          if (self%isPeriodic) then
            call curbox%Boundary(rx, ry, rz)
          endif
          rsq = rx*rx + ry*ry + rz*rz
          
          ! Direct Coulomb interaction
          r = sqrt(rsq)
          E_pair = qi * qj * coulombConst / r
          E_far = E_far + E_pair
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_pair
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_pair
        enddo
      enddo
    enddo
  enddo

end subroutine
!=============================================================================
subroutine ComputeFarFieldMultipole_FMM(self, curbox, E_far)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_far

  integer :: iLeaf, nodeIdx, j, iAtom
  integer :: atmType
  real(dp) :: qi, dx, dy, dz
  real(dp), pointer :: atoms(:,:) => null()
  complex(dp) :: Rlm(self%nExpansion)
  complex(dp) :: phi
  integer :: l, m, idx
  
  integer :: molIdx, molStart, molEnd, iA, iB
  integer :: leafA, leafB, k, atmTypeA, atmTypeB
  real(dp) :: qA, qB, rx, ry, rz, r, E_pair, E_intra_corr
  logical :: isNear
  
  call curbox%GetCoordinates(atoms)
  
  E_far = 0.0_dp
  
  do iLeaf = 1, self%nLeaves
    nodeIdx = self%leafNodes(iLeaf)
    
    if (.not. allocated(self%nodes(nodeIdx)%particles)) cycle
    
    do j = 1, size(self%nodes(nodeIdx)%particles)
      iAtom = self%nodes(nodeIdx)%particles(j)
      atmType = curbox%AtomType(iAtom)
      qi = self%q(atmType)
      
      if (abs(qi) < 1.0E-15_dp) cycle
      
      dx = atoms(1, iAtom) - self%nodes(nodeIdx)%center(1)
      dy = atoms(2, iAtom) - self%nodes(nodeIdx)%center(2)
      dz = atoms(3, iAtom) - self%nodes(nodeIdx)%center(3)
      
      call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm)
      
      phi = (0.0_dp, 0.0_dp)
      do l = 0, self%expansionOrder
        do m = -l, l
          idx = SH_Index(l, m)
          phi = phi + self%nodes(nodeIdx)%L(idx) * Rlm(idx)
        enddo
      enddo
      
      E_far = E_far + 0.5_dp * qi * real(phi, dp) * coulombConst
    enddo
  enddo

  ! Intramolecular far-field correction: the local expansions include
  ! contributions from ALL far-field charges, including same-molecule atoms
  ! whose leaf cells happen to be in the far-field.  Subtract those.
  E_intra_corr = 0.0_dp
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    molIdx = curbox%MolIndx(iAtom)
    molStart = curbox%MolStartIndx(molIdx)
    if (iAtom /= molStart) cycle
    molEnd = curbox%MolEndIndx(molIdx)
    if (molEnd == molStart) cycle

    do iA = molStart, molEnd - 1
      if (.not. curbox%IsActive(iA)) cycle
      leafA = self%particleToLeaf(iA)
      if (leafA <= 0) cycle
      atmTypeA = curbox%AtomType(iA)
      qA = self%q(atmTypeA)
      if (abs(qA) < 1.0E-15_dp) cycle

      do iB = iA + 1, molEnd
        if (.not. curbox%IsActive(iB)) cycle
        leafB = self%particleToLeaf(iB)
        if (leafB <= 0) cycle
        if (leafA == leafB) cycle

        atmTypeB = curbox%AtomType(iB)
        qB = self%q(atmTypeB)
        if (abs(qB) < 1.0E-15_dp) cycle

        ! Check A→B: is leafB in leafA's far-field?
        isNear = .false.
        do k = 1, self%nodes(leafA)%nNear
          if (self%nodes(leafA)%nearList(k) == leafB) then
            isNear = .true.
            exit
          endif
        enddo
        if (isNear) then
          idx = 0
        else
          idx = 1
        endif

        ! Check B→A: is leafA in leafB's far-field?
        isNear = .false.
        do k = 1, self%nodes(leafB)%nNear
          if (self%nodes(leafB)%nearList(k) == leafA) then
            isNear = .true.
            exit
          endif
        enddo
        if (.not. isNear) idx = idx + 1

        if (idx > 0) then
          rx = atoms(1, iA) - atoms(1, iB)
          ry = atoms(2, iA) - atoms(2, iB)
          rz = atoms(3, iA) - atoms(3, iB)
          if (self%isPeriodic) call curbox%Boundary(rx, ry, rz)
          r = sqrt(rx*rx + ry*ry + rz*rz)

          E_pair = qA * qB * coulombConst / r
          E_intra_corr = E_intra_corr + 0.5_dp * real(idx, dp) * E_pair
        endif
      enddo
    enddo
  enddo

  block
    use ParallelVar, only: nout
    write(nout,*) " E_far_raw=", E_far, " E_intra_corr=", E_intra_corr
  end block
  E_far = E_far - E_intra_corr

end subroutine
!=============================================================================
subroutine DiagnosticMultipoleCompare_FMM(self, curbox, E_far_direct)
  use ParallelVar, only: nout
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(in) :: E_far_direct

  integer :: iLeaf, nodeIdx, j, iAtom, jAtom
  integer :: atmType, atmType2, k, nodeJ, ii
  real(dp) :: qi, qj, dx, dy, dz, rx, ry, rz, rsq, r
  real(dp), pointer :: atoms(:,:) => null()
  complex(dp) :: Rlm(self%nExpansion)
  complex(dp) :: phi
  integer :: l, m, idx
  real(dp) :: E_far_multi, phi_multi, phi_direct
  integer :: nSampled
  real(dp) :: maxErr, avgErr, relErr

  call curbox%GetCoordinates(atoms)

  E_far_multi = 0.0_dp
  maxErr = 0.0_dp
  avgErr = 0.0_dp
  nSampled = 0

  do iLeaf = 1, self%nLeaves
    nodeIdx = self%leafNodes(iLeaf)
    if (.not. allocated(self%nodes(nodeIdx)%particles)) cycle

    do j = 1, size(self%nodes(nodeIdx)%particles)
      iAtom = self%nodes(nodeIdx)%particles(j)
      atmType = curbox%AtomType(iAtom)
      qi = self%q(atmType)
      if (abs(qi) < 1.0E-15_dp) cycle

      dx = atoms(1, iAtom) - self%nodes(nodeIdx)%center(1)
      dy = atoms(2, iAtom) - self%nodes(nodeIdx)%center(2)
      dz = atoms(3, iAtom) - self%nodes(nodeIdx)%center(3)
      call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm)

      phi = (0.0_dp, 0.0_dp)
      do l = 0, self%expansionOrder
        do m = -l, l
          idx = SH_Index(l, m)
          phi = phi + self%nodes(nodeIdx)%L(idx) * Rlm(idx)
        enddo
      enddo
      phi_multi = real(phi, dp)

      phi_direct = 0.0_dp
      do k = 1, self%nodes(nodeIdx)%nFar
        nodeJ = self%nodes(nodeIdx)%farList(k)
        if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
        do ii = 1, size(self%nodes(nodeJ)%particles)
          jAtom = self%nodes(nodeJ)%particles(ii)
          atmType2 = curbox%AtomType(jAtom)
          qj = self%q(atmType2)
          if (abs(qj) < 1.0E-15_dp) cycle
          rx = atoms(1, iAtom) - atoms(1, jAtom)
          ry = atoms(2, iAtom) - atoms(2, jAtom)
          rz = atoms(3, iAtom) - atoms(3, jAtom)
          if (self%isPeriodic) call curbox%Boundary(rx, ry, rz)
          rsq = rx*rx + ry*ry + rz*rz
          r = sqrt(rsq)
          phi_direct = phi_direct + qj / r
        enddo
      enddo

      E_far_multi = E_far_multi + 0.5_dp * qi * phi_multi * coulombConst

      if (abs(phi_direct) > 1.0E-15_dp) then
        relErr = abs(phi_multi - phi_direct) / abs(phi_direct)
        if (relErr > maxErr) maxErr = relErr
        avgErr = avgErr + relErr
        nSampled = nSampled + 1
      endif

      if (nSampled <= 5) then
        ! Count intramolecular atoms in far-field cells
        block
          integer :: nIntraFar, molI, jj, jjAtom
          real(dp) :: phi_intra
          nIntraFar = 0
          phi_intra = 0.0_dp
          molI = curbox%MolIndx(iAtom)
          do k = 1, self%nodes(nodeIdx)%nFar
            nodeJ = self%nodes(nodeIdx)%farList(k)
            if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
            do jj = 1, size(self%nodes(nodeJ)%particles)
              jjAtom = self%nodes(nodeJ)%particles(jj)
              if (curbox%MolIndx(jjAtom) == molI) then
                nIntraFar = nIntraFar + 1
                rx = atoms(1,iAtom)-atoms(1,jjAtom)
                ry = atoms(2,iAtom)-atoms(2,jjAtom)
                rz = atoms(3,iAtom)-atoms(3,jjAtom)
                if (self%isPeriodic) call curbox%Boundary(rx,ry,rz)
                phi_intra = phi_intra + self%q(curbox%AtomType(jjAtom)) / sqrt(rx*rx+ry*ry+rz*rz)
              endif
            enddo
          enddo
          write(nout, '(A,I6,A,ES14.6,A,ES14.6,A,F8.4,A,I3,A,ES12.4,A,F8.4)') &
            "  Atom", iAtom, " multi=", phi_multi, " direct=", phi_direct, &
            " err=", merge(abs(phi_multi-phi_direct)/abs(phi_direct), 0.0_dp, abs(phi_direct)>1e-15_dp)*100.0_dp, &
            "% intraFar=", nIntraFar, " phiIntra=", phi_intra, " hw=", self%nodes(nodeIdx)%halfWidth
        end block
      endif
    enddo
  enddo

  if (nSampled > 0) avgErr = avgErr / real(nSampled, dp)

  write(nout,*) "=== Multipole vs Direct Far-Field Diagnostic ==="
  write(nout,*) "  E_far (multipole):", E_far_multi
  write(nout,*) "  E_far (direct):   ", E_far_direct
  if (abs(E_far_direct) > 1.0E-15_dp) then
    write(nout,*) "  Rel. error:       ", abs(E_far_multi - E_far_direct)/abs(E_far_direct)*100.0_dp, "%"
  endif
  write(nout,*) "  Per-particle max err:", maxErr*100.0_dp, "%"
  write(nout,*) "  Per-particle avg err:", avgErr*100.0_dp, "%"
  write(nout,*) "  Particles sampled:  ", nSampled

  ! Check if error correlates with whether atom is in the same cell as other atoms
  ! in its molecule
  block
    integer :: nSameCell, nDiffCell
    integer :: aIdx, molI2, molStart, molEnd, bIdx, leafA, leafB
    real(dp) :: errSame, errDiff, phiM, phiD, re2
    nSameCell = 0; nDiffCell = 0; errSame = 0.0_dp; errDiff = 0.0_dp
    do iLeaf = 1, self%nLeaves
      nodeIdx = self%leafNodes(iLeaf)
      if (.not. allocated(self%nodes(nodeIdx)%particles)) cycle
      do j = 1, size(self%nodes(nodeIdx)%particles)
        aIdx = self%nodes(nodeIdx)%particles(j)
        atmType = curbox%AtomType(aIdx)
        qi = self%q(atmType)
        if (abs(qi) < 1.0E-15_dp) cycle
        leafA = nodeIdx
        molI2 = curbox%MolIndx(aIdx)
        molStart = curbox%MolStartIndx(molI2)
        molEnd = curbox%MolEndIndx(molI2)
        ! Check if all atoms in molecule are in this same leaf
        block
          logical :: allSame
          integer :: bb
          allSame = .true.
          do bb = molStart, molEnd
            if (bb == aIdx) cycle
            if (.not. curbox%IsActive(bb)) cycle
            if (self%particleToLeaf(bb) /= leafA) then
              allSame = .false.
              exit
            endif
          enddo
          ! Compute phi_multi and phi_direct for this atom
          dx = atoms(1,aIdx) - self%nodes(nodeIdx)%center(1)
          dy = atoms(2,aIdx) - self%nodes(nodeIdx)%center(2)
          dz = atoms(3,aIdx) - self%nodes(nodeIdx)%center(3)
          call SH_ComputeSolidHarmonic(dx,dy,dz,self%expansionOrder,Rlm)
          phi = (0.0_dp,0.0_dp)
          do l = 0, self%expansionOrder
            do m = -l, l
              idx = SH_Index(l,m)
              phi = phi + self%nodes(nodeIdx)%L(idx) * Rlm(idx)
            enddo
          enddo
          phiM = real(phi, dp)
          phiD = 0.0_dp
          do k = 1, self%nodes(nodeIdx)%nFar
            nodeJ = self%nodes(nodeIdx)%farList(k)
            if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
            do ii = 1, size(self%nodes(nodeJ)%particles)
              bIdx = self%nodes(nodeJ)%particles(ii)
              atmType2 = curbox%AtomType(bIdx)
              qj = self%q(atmType2)
              if (abs(qj) < 1.0E-15_dp) cycle
              rx = atoms(1,aIdx)-atoms(1,bIdx)
              ry = atoms(2,aIdx)-atoms(2,bIdx)
              rz = atoms(3,aIdx)-atoms(3,bIdx)
              if (self%isPeriodic) call curbox%Boundary(rx,ry,rz)
              rsq = rx*rx+ry*ry+rz*rz
              r = sqrt(rsq)
              phiD = phiD + qj/r
            enddo
          enddo
          if (abs(phiD) > 1e-15_dp) then
            re2 = abs(phiM - phiD)/abs(phiD)
            if (allSame) then
              nSameCell = nSameCell + 1
              errSame = errSame + re2
            else
              nDiffCell = nDiffCell + 1
              errDiff = errDiff + re2
            endif
          endif
        end block
      enddo
    enddo
    if (nSameCell > 0) errSame = errSame / real(nSameCell, dp)
    if (nDiffCell > 0) errDiff = errDiff / real(nDiffCell, dp)
    write(nout,*) "  Atoms w/ all mol-mates in same cell:", nSameCell, " avgErr:", errSame*100.0_dp, "%"
    write(nout,*) "  Atoms w/ mol-mates in diff cell:   ", nDiffCell, " avgErr:", errDiff*100.0_dp, "%"
  end block
  write(nout,*) "================================================="

end subroutine
!=============================================================================
subroutine SingleCellM2LTest_FMM(self, curbox)
  use ParallelVar, only: nout
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  integer :: tgtLeaf, tgtNode, srcNode, p
  integer :: j, k, ii, iAtom, jAtom, atmType, atmType2, nSrc
  real(dp) :: qi, qj, dx, dy, dz, rx, ry, rz, rsq, r
  real(dp), pointer :: atoms(:,:) => null()
  complex(dp), allocatable :: M_src(:), L_single(:), Rlm(:), Slm(:)
  complex(dp) :: phi
  real(dp) :: phi_m2l, phi_exact
  integer :: l, m, idx

  call curbox%GetCoordinates(atoms)
  p = self%expansionOrder

  allocate(M_src((p+1)*(p+1)), L_single((p+1)*(p+1)), Rlm((p+1)*(p+1)))

  ! Pick first leaf with particles and a far-field source
  tgtLeaf = 0
  do ii = 1, self%nLeaves
    tgtNode = self%leafNodes(ii)
    if (.not. allocated(self%nodes(tgtNode)%particles)) cycle
    if (self%nodes(tgtNode)%nFar < 1) cycle
    tgtLeaf = ii
    exit
  enddo
  if (tgtLeaf == 0) then
    write(nout,*) "No suitable target cell for M2L test"
    deallocate(M_src, L_single, Rlm)
    return
  endif
  tgtNode = self%leafNodes(tgtLeaf)
  srcNode = self%nodes(tgtNode)%farList(1)

  ! Count source particles
  nSrc = 0
  if (allocated(self%nodes(srcNode)%particles)) nSrc = size(self%nodes(srcNode)%particles)

  write(nout,*) "=== Single Cell M2L Test ==="
  write(nout,*) "  Target leaf:", tgtLeaf, " node:", tgtNode, " nPart:", size(self%nodes(tgtNode)%particles)
  write(nout,*) "  Source node:", srcNode, " nPart:", nSrc
  write(nout,*) "  Target center:", self%nodes(tgtNode)%center
  write(nout,*) "  Source center:", self%nodes(srcNode)%center

  ! Compute source multipole manually (P2M for source cell)
  M_src = (0.0_dp, 0.0_dp)
  if (nSrc > 0) then
    do j = 1, nSrc
      jAtom = self%nodes(srcNode)%particles(j)
      atmType2 = curbox%AtomType(jAtom)
      qj = self%q(atmType2)
      if (abs(qj) < 1.0E-15_dp) cycle
      dx = atoms(1,jAtom) - self%nodes(srcNode)%center(1)
      dy = atoms(2,jAtom) - self%nodes(srcNode)%center(2)
      dz = atoms(3,jAtom) - self%nodes(srcNode)%center(3)
      call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
      do l = 0, p
        do m = -l, l
          M_src(SH_Index(l,m)) = M_src(SH_Index(l,m)) + qj * Rlm(SH_Index(l,-m))
        enddo
      enddo
    enddo
  endif
  write(nout,*) "  M_src(0,0) =", M_src(SH_Index(0,0)), " (should be ~0 if neutral)"

  ! M2L from this single source
  L_single = (0.0_dp, 0.0_dp)
  dx = self%nodes(tgtNode)%center(1) - self%nodes(srcNode)%center(1)
  dy = self%nodes(tgtNode)%center(2) - self%nodes(srcNode)%center(2)
  dz = self%nodes(tgtNode)%center(3) - self%nodes(srcNode)%center(3)
  if (self%isPeriodic) then
    if (dx > self%cachedLx/2.0_dp) dx = dx - self%cachedLx
    if (dx < -self%cachedLx/2.0_dp) dx = dx + self%cachedLx
    if (dy > self%cachedLy/2.0_dp) dy = dy - self%cachedLy
    if (dy < -self%cachedLy/2.0_dp) dy = dy + self%cachedLy
    if (dz > self%cachedLz/2.0_dp) dz = dz - self%cachedLz
    if (dz < -self%cachedLz/2.0_dp) dz = dz + self%cachedLz
  endif
  write(nout,*) "  Translation d:", dx, dy, dz, " |d|:", sqrt(dx*dx+dy*dy+dz*dz)
  call SH_M2L_Coeff(dx, dy, dz, p, M_src, L_single)

  ! Evaluate L2P at each target particle and compare with direct sum from source
  do j = 1, size(self%nodes(tgtNode)%particles)
    iAtom = self%nodes(tgtNode)%particles(j)
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    if (abs(qi) < 1.0E-15_dp) cycle

    dx = atoms(1,iAtom) - self%nodes(tgtNode)%center(1)
    dy = atoms(2,iAtom) - self%nodes(tgtNode)%center(2)
    dz = atoms(3,iAtom) - self%nodes(tgtNode)%center(3)
    call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)

    phi = (0.0_dp, 0.0_dp)
    do l = 0, p
      do m = -l, l
        idx = SH_Index(l, m)
        phi = phi + L_single(idx) * Rlm(idx)
      enddo
    enddo
    phi_m2l = real(phi, dp)

    phi_exact = 0.0_dp
    if (nSrc > 0) then
      do k = 1, nSrc
        jAtom = self%nodes(srcNode)%particles(k)
        atmType2 = curbox%AtomType(jAtom)
        qj = self%q(atmType2)
        if (abs(qj) < 1.0E-15_dp) cycle
        rx = atoms(1,iAtom) - atoms(1,jAtom)
        ry = atoms(2,iAtom) - atoms(2,jAtom)
        rz = atoms(3,iAtom) - atoms(3,jAtom)
        if (self%isPeriodic) call curbox%Boundary(rx,ry,rz)
        rsq = rx*rx+ry*ry+rz*rz
        r = sqrt(rsq)
        phi_exact = phi_exact + qj/r
      enddo
    endif

    ! Also compute L' direct: L'_l^m = sum_i q_i S_l^m(y_i - c_T)
    block
      complex(dp), allocatable :: L_direct(:), Slm_src(:)
      real(dp) :: phi_direct_L, sx, sy, sz
      complex(dp) :: phi_d
      allocate(L_direct((p+1)*(p+1)), Slm_src((p+1)*(p+1)))
      L_direct = (0.0_dp, 0.0_dp)
      if (nSrc > 0) then
        do k = 1, nSrc
          jAtom = self%nodes(srcNode)%particles(k)
          atmType2 = curbox%AtomType(jAtom)
          qj = self%q(atmType2)
          if (abs(qj) < 1.0E-15_dp) cycle
          sx = atoms(1,jAtom) - self%nodes(tgtNode)%center(1)
          sy = atoms(2,jAtom) - self%nodes(tgtNode)%center(2)
          sz = atoms(3,jAtom) - self%nodes(tgtNode)%center(3)
          if (self%isPeriodic) call curbox%Boundary(sx,sy,sz)
          call SH_ComputeIrregularSolid(sx, sy, sz, p, Slm_src)
          L_direct = L_direct + qj * Slm_src
        enddo
      endif
      ! Evaluate: phi = sum R'(l,-m) * L_direct(l,m)
      phi_d = (0.0_dp, 0.0_dp)
      do l = 0, p
        do m = -l, l
          phi_d = phi_d + Rlm(SH_Index(l,-m)) * L_direct(SH_Index(l,m))
        enddo
      enddo
      phi_direct_L = real(phi_d, dp)
      write(nout,'(A,I6,A,ES14.6,A,ES14.6,A,ES14.6,A,F8.4,A,F8.4,A)') &
        "  Atom", iAtom, " m2l=", phi_m2l, " directL=", phi_direct_L, &
        " exact=", phi_exact, " m2l_err=", &
        merge(abs(phi_m2l-phi_exact)/abs(phi_exact), 0.0_dp, abs(phi_exact)>1e-15_dp)*100.0_dp, &
        "% dL_err=", &
        merge(abs(phi_direct_L-phi_exact)/abs(phi_exact), 0.0_dp, abs(phi_exact)>1e-15_dp)*100.0_dp, "%"

      ! Print L coefficient comparison for first atom
      if (j == 1) then
        write(nout,*) "  --- L coefficients: M2L vs Direct (L'(-k)) ---"
        do l = 0, min(3, p)
          do m = -l, l
            idx = SH_Index(l, m)
            write(nout,'(A,I2,A,I3,A,2ES12.4,A,2ES12.4)') &
              "    L(", l, ",", m, ") M2L=", real(L_single(idx)), aimag(L_single(idx)), &
              "  direct(-k)=", real(L_direct(SH_Index(l,-m))), aimag(L_direct(SH_Index(l,-m)))
          enddo
        enddo
      endif

      deallocate(L_direct, Slm_src)
    end block
  enddo

  ! Now compare accumulated L from pipeline vs manual accumulation
  write(nout,*) "--- Accumulated L comparison ---"
  block
    complex(dp), allocatable :: L_manual(:), Slm_tmp(:)
    integer :: kk, sNode
    real(dp) :: ddx, ddy, ddz

    allocate(L_manual((p+1)*(p+1)), Slm_tmp((2*p+1)*(2*p+1)))
    L_manual = (0.0_dp, 0.0_dp)

    do kk = 1, self%nodes(tgtNode)%nFar
      sNode = self%nodes(tgtNode)%farList(kk)
      ddx = self%nodes(tgtNode)%center(1) - self%nodes(sNode)%center(1)
      ddy = self%nodes(tgtNode)%center(2) - self%nodes(sNode)%center(2)
      ddz = self%nodes(tgtNode)%center(3) - self%nodes(sNode)%center(3)
      if (self%isPeriodic) then
        if (ddx > self%cachedLx/2.0_dp) ddx = ddx - self%cachedLx
        if (ddx < -self%cachedLx/2.0_dp) ddx = ddx + self%cachedLx
        if (ddy > self%cachedLy/2.0_dp) ddy = ddy - self%cachedLy
        if (ddy < -self%cachedLy/2.0_dp) ddy = ddy + self%cachedLy
        if (ddz > self%cachedLz/2.0_dp) ddz = ddz - self%cachedLz
        if (ddz < -self%cachedLz/2.0_dp) ddz = ddz + self%cachedLz
      endif
      call SH_M2L_Coeff(ddx, ddy, ddz, p, self%nodes(sNode)%M, L_manual)
    enddo

    write(nout,*) "  nFar:", self%nodes(tgtNode)%nFar
    write(nout,*) "  Pipeline L(0,0):", self%nodes(tgtNode)%L(SH_Index(0,0))
    write(nout,*) "  Manual   L(0,0):", L_manual(SH_Index(0,0))
    write(nout,*) "  Pipeline L(1,0):", self%nodes(tgtNode)%L(SH_Index(1,0))
    write(nout,*) "  Manual   L(1,0):", L_manual(SH_Index(1,0))
    write(nout,*) "  Pipeline L(1,1):", self%nodes(tgtNode)%L(SH_Index(1,1))
    write(nout,*) "  Manual   L(1,1):", L_manual(SH_Index(1,1))

    ! L2P from manual L vs from pipeline L
    do j = 1, min(3, size(self%nodes(tgtNode)%particles))
      iAtom = self%nodes(tgtNode)%particles(j)
      atmType = curbox%AtomType(iAtom)
      qi = self%q(atmType)
      if (abs(qi) < 1.0E-15_dp) cycle

      dx = atoms(1,iAtom) - self%nodes(tgtNode)%center(1)
      dy = atoms(2,iAtom) - self%nodes(tgtNode)%center(2)
      dz = atoms(3,iAtom) - self%nodes(tgtNode)%center(3)
      call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)

      phi = (0.0_dp, 0.0_dp)
      do l = 0, p
        do m = -l, l
          idx = SH_Index(l, m)
          phi = phi + self%nodes(tgtNode)%L(idx) * Rlm(idx)
        enddo
      enddo
      phi_m2l = real(phi, dp)

      phi = (0.0_dp, 0.0_dp)
      do l = 0, p
        do m = -l, l
          idx = SH_Index(l, m)
          phi = phi + L_manual(idx) * Rlm(idx)
        enddo
      enddo

      phi_exact = 0.0_dp
      do kk = 1, self%nodes(tgtNode)%nFar
        sNode = self%nodes(tgtNode)%farList(kk)
        if (.not. allocated(self%nodes(sNode)%particles)) cycle
        do ii = 1, size(self%nodes(sNode)%particles)
          jAtom = self%nodes(sNode)%particles(ii)
          atmType2 = curbox%AtomType(jAtom)
          qj = self%q(atmType2)
          if (abs(qj) < 1.0E-15_dp) cycle
          rx = atoms(1,iAtom)-atoms(1,jAtom)
          ry = atoms(2,iAtom)-atoms(2,jAtom)
          rz = atoms(3,iAtom)-atoms(3,jAtom)
          if (self%isPeriodic) call curbox%Boundary(rx,ry,rz)
          rsq = rx*rx+ry*ry+rz*rz
          r = sqrt(rsq)
          phi_exact = phi_exact + qj/r
        enddo
      enddo

      write(nout,'(A,I6,A,ES14.6,A,ES14.6,A,ES14.6)') &
        "  Atom", iAtom, " pipeline=", phi_m2l, " manual=", real(phi,dp), " exact=", phi_exact
    enddo

    deallocate(L_manual, Slm_tmp)
  end block
  write(nout,*) "============================"

  deallocate(M_src, L_single, Rlm)
end subroutine
!=============================================================================
subroutine ComputeSelfEnergy_FMM(self, curbox, E_self)
  ! Self-energy correction for FMM
  !
  ! Unlike Ewald summation, FMM does not use a screening function (erfc),
  ! so there is no explicit self-interaction to remove from the real-space sum.
  !
  ! In FMM, self-interactions are naturally avoided:
  !   - Near-field: we explicitly skip i=j pairs
  !   - Far-field: a particle's own cell is in the near-field, not far-field
  !
  ! Therefore, no self-energy correction is needed for the standard FMM.
  ! The periodic correction (dipole term) handles the periodic images separately.
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_self

  ! No self-energy correction for standard FMM (unlike Ewald)
  E_self = 0.0_dp

end subroutine
!=============================================================================
subroutine ComputePeriodicCorrection_FMM(self, curbox, E_periodic)
  ! Compute periodic correction using lattice sum at root level
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_periodic

  real(dp) :: Lx, Ly, Lz, volume
  real(dp) :: dipole(3), dipole_sq
  integer :: iAtom, atmType
  real(dp) :: qi
  real(dp), pointer :: atoms(:,:) => null()
  
  call curbox%GetCoordinates(atoms)
  
  ! Get box dimensions
  select type(curbox)
    class is(OrthoBox)
      Lx = curbox%boxLx
      Ly = curbox%boxLy
      Lz = curbox%boxLz
    class is(CubeBox)
      Lx = curbox%boxL
      Ly = curbox%boxL
      Lz = curbox%boxL
    class default
      E_periodic = 0.0_dp
      return
  end select
  
  volume = Lx * Ly * Lz
  
  ! Compute total dipole moment
  dipole = 0.0_dp
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    dipole(1) = dipole(1) + qi * atoms(1, iAtom)
    dipole(2) = dipole(2) + qi * atoms(2, iAtom)
    dipole(3) = dipole(3) + qi * atoms(3, iAtom)
  enddo
  
  dipole_sq = dipole(1)**2 + dipole(2)**2 + dipole(3)**2
  
  ! Surface dipole correction (tin-foil boundary conditions)
  E_periodic = 2.0_dp * pi * coulombConst * dipole_sq / (3.0_dp * volume)

end subroutine
!=============================================================================
subroutine UpdateMultipoleForMove_FMM(self, curbox, iAtom, x_old, y_old, z_old, &
                                       x_new, y_new, z_new)
  ! Incrementally update multipoles when a particle moves
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  integer, intent(in) :: iAtom
  real(dp), intent(in) :: x_old, y_old, z_old
  real(dp), intent(in) :: x_new, y_new, z_new

  integer :: oldLeaf, newLeaf, atmType, l, m
  real(dp) :: qi, dx, dy, dz
  complex(dp) :: Rlm_old(self%nExpansion), Rlm_new(self%nExpansion)
  
  atmType = curbox%AtomType(iAtom)
  qi = self%q(atmType)
  
  if (abs(qi) < 1.0E-15_dp) return
  
  oldLeaf = self%particleToLeaf(iAtom)
  newLeaf = self%FindLeafForPosition(x_new, y_new, z_new)
  
  ! Save for potential rejection
  self%lastOldLeaf = oldLeaf
  self%lastNewLeaf = newLeaf
  self%particleToLeaf_old = self%particleToLeaf
  
  ! Save old multipoles for nodes that will be affected
  self%nodes(oldLeaf)%M_old = self%nodes(oldLeaf)%M
  if (newLeaf /= oldLeaf) then
    self%nodes(newLeaf)%M_old = self%nodes(newLeaf)%M
  endif
  
  ! Remove old contribution from old leaf
  dx = x_old - self%nodes(oldLeaf)%center(1)
  dy = y_old - self%nodes(oldLeaf)%center(2)
  dz = z_old - self%nodes(oldLeaf)%center(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm_old)
  do l = 0, self%expansionOrder
    do m = -l, l
      self%nodes(oldLeaf)%M(SH_Index(l,m)) = self%nodes(oldLeaf)%M(SH_Index(l,m)) &
        - qi * Rlm_old(SH_Index(l,-m))
    enddo
  enddo
  
  ! Add new contribution to new leaf
  dx = x_new - self%nodes(newLeaf)%center(1)
  dy = y_new - self%nodes(newLeaf)%center(2)
  dz = z_new - self%nodes(newLeaf)%center(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm_new)
  do l = 0, self%expansionOrder
    do m = -l, l
      self%nodes(newLeaf)%M(SH_Index(l,m)) = self%nodes(newLeaf)%M(SH_Index(l,m)) &
        + qi * Rlm_new(SH_Index(l,-m))
    enddo
  enddo
  
  ! Update particle to leaf mapping
  self%particleToLeaf(iAtom) = newLeaf

  ! If crossed to new leaf, update particles lists (for near direct sums)
  if (newLeaf /= oldLeaf) then
    call MoveParticleBetweenLeaves(self%nodes(oldLeaf), self%nodes(newLeaf), iAtom)
  endif

end subroutine
!=============================================================================
subroutine AcceptUpdate_FMM(self)
  ! Called when a move is accepted
  implicit none
  class(Pair_FMM), intent(inout) :: self

  ! Nothing to do - multipoles already updated

end subroutine
!=============================================================================
subroutine RejectUpdate_FMM(self)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  integer :: iS, nodeIdx

  ! Restore multipoles for all saved leaves
  do iS = 1, self%nSavedLeaves
    nodeIdx = self%savedLeaves(iS)
    if (nodeIdx > 0 .and. nodeIdx <= self%nNodes) then
      self%nodes(nodeIdx)%M = self%nodes(nodeIdx)%M_old
    endif
  enddo
  self%nSavedLeaves = 0
  
  ! Restore particle mapping
  self%particleToLeaf = self%particleToLeaf_old

end subroutine
!=============================================================================
subroutine ResolvePendingMove_FMM(self, curbox)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), pointer :: atoms(:,:) => null()
  
  if (.not. self%hasPendingMove) return
  if (.not. self%useMultipole .or. .not. self%treeBuilt) then
    self%hasPendingMove = .false.
    return
  endif
  
  call curbox%GetCoordinates(atoms)
  
  if (atoms(1, self%pendingAtomIdx) /= self%pendingOldX .or. &
      atoms(2, self%pendingAtomIdx) /= self%pendingOldY .or. &
      atoms(3, self%pendingAtomIdx) /= self%pendingOldZ) then
    ! Accepted: M is already tentatively updated, keep it
    call self%AcceptUpdate()
  else
    ! Rejected: revert M to pre-move state
    call self%RejectUpdate()
  endif
  
  self%hasPendingMove = .false.
end subroutine
!=============================================================================
function EvalFarFieldPotential_FMM(self, leafIdx, x, y, z) result(phi)
  ! Recompute L on-the-fly from current M of all far-field cells,
  ! then evaluate L2P at the given position. O(n_far * p^4).
  implicit none
  class(Pair_FMM), intent(in) :: self
  integer, intent(in) :: leafIdx
  real(dp), intent(in) :: x, y, z
  real(dp) :: phi
  
  complex(dp) :: L_temp(self%nExpansion), Rlm(self%nExpansion)
  complex(dp) :: phi_c
  integer :: k, nodeJ, l, m, idx
  real(dp) :: dx, dy, dz
  
  L_temp = (0.0_dp, 0.0_dp)
  do k = 1, self%nodes(leafIdx)%nFar
    nodeJ = self%nodes(leafIdx)%farList(k)
    dx = self%nodes(leafIdx)%center(1) - self%nodes(nodeJ)%center(1)
    dy = self%nodes(leafIdx)%center(2) - self%nodes(nodeJ)%center(2)
    dz = self%nodes(leafIdx)%center(3) - self%nodes(nodeJ)%center(3)
    if (self%isPeriodic) then
      if (dx >  self%cachedLx*0.5_dp) dx = dx - self%cachedLx
      if (dx < -self%cachedLx*0.5_dp) dx = dx + self%cachedLx
      if (dy >  self%cachedLy*0.5_dp) dy = dy - self%cachedLy
      if (dy < -self%cachedLy*0.5_dp) dy = dy + self%cachedLy
      if (dz >  self%cachedLz*0.5_dp) dz = dz - self%cachedLz
      if (dz < -self%cachedLz*0.5_dp) dz = dz + self%cachedLz
    endif
    call SH_M2L_Coeff(dx, dy, dz, self%expansionOrder, self%nodes(nodeJ)%M, L_temp)
  enddo
  
  dx = x - self%nodes(leafIdx)%center(1)
  dy = y - self%nodes(leafIdx)%center(2)
  dz = z - self%nodes(leafIdx)%center(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm)
  
  phi_c = (0.0_dp, 0.0_dp)
  do l = 0, self%expansionOrder
    do m = -l, l
      idx = SH_Index(l, m)
      phi_c = phi_c + L_temp(idx) * Rlm(idx)
    enddo
  enddo
  
  phi = real(phi_c, dp)
end function
!=============================================================================
subroutine TentativeMultipoleUpdate_FMM(self, curbox, disp, dispLen, atoms)
  ! After a displacement DiffECalc, tentatively update M for all displaced atoms.
  ! Saves M_old collectively so that RejectUpdate can revert all at once.
  ! IMPORTANT: new positions must be wrapped to the primary cell (matching
  ! what UpdatePosition does) before computing solid harmonics.
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  integer, intent(in) :: dispLen
  real(dp), pointer, intent(in) :: atoms(:,:)
  
  integer :: iDisp, iAtom, oldLeaf, newLeaf, atmType, l, m, iS
  real(dp) :: qi, dx, dy, dz
  real(dp) :: xw, yw, zw
  complex(dp) :: Rlm_tmp(self%nExpansion)
  logical :: alreadySaved
  
  self%particleToLeaf_old = self%particleToLeaf
  self%nSavedLeaves = 0
  
  ! First pass: collect and save M_old for all affected leaves
  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    if (abs(qi) < 1.0E-15_dp) cycle
    
    oldLeaf = self%particleToLeaf(iAtom)
    
    ! Wrap new position to primary cell (same as UpdatePosition does)
    xw = disp(iDisp)%x_new
    yw = disp(iDisp)%y_new
    zw = disp(iDisp)%z_new
    call curbox%Boundary(xw, yw, zw)
    
    newLeaf = self%FindLeafForPosition(xw, yw, zw)
    if (oldLeaf <= 0) cycle
    if (newLeaf <= 0) newLeaf = oldLeaf
    
    alreadySaved = .false.
    do iS = 1, self%nSavedLeaves
      if (self%savedLeaves(iS) == oldLeaf) then
        alreadySaved = .true.
        exit
      endif
    enddo
    if (.not. alreadySaved .and. self%nSavedLeaves < 20) then
      self%nSavedLeaves = self%nSavedLeaves + 1
      self%savedLeaves(self%nSavedLeaves) = oldLeaf
      self%nodes(oldLeaf)%M_old = self%nodes(oldLeaf)%M
    endif
    
    if (newLeaf /= oldLeaf) then
      alreadySaved = .false.
      do iS = 1, self%nSavedLeaves
        if (self%savedLeaves(iS) == newLeaf) then
          alreadySaved = .true.
          exit
        endif
      enddo
      if (.not. alreadySaved .and. self%nSavedLeaves < 20) then
        self%nSavedLeaves = self%nSavedLeaves + 1
        self%savedLeaves(self%nSavedLeaves) = newLeaf
        self%nodes(newLeaf)%M_old = self%nodes(newLeaf)%M
      endif
    endif
  enddo
  
  ! Second pass: update M for each displaced atom (using wrapped positions)
  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    if (abs(qi) < 1.0E-15_dp) cycle
    
    oldLeaf = self%particleToLeaf(iAtom)
    
    xw = disp(iDisp)%x_new
    yw = disp(iDisp)%y_new
    zw = disp(iDisp)%z_new
    call curbox%Boundary(xw, yw, zw)
    
    newLeaf = self%FindLeafForPosition(xw, yw, zw)
    if (oldLeaf <= 0) cycle
    if (newLeaf <= 0) newLeaf = oldLeaf
    
    ! Remove old contribution from old leaf
    dx = atoms(1, iAtom) - self%nodes(oldLeaf)%center(1)
    dy = atoms(2, iAtom) - self%nodes(oldLeaf)%center(2)
    dz = atoms(3, iAtom) - self%nodes(oldLeaf)%center(3)
    call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm_tmp)
    do l = 0, self%expansionOrder
      do m = -l, l
        self%nodes(oldLeaf)%M(SH_Index(l,m)) = self%nodes(oldLeaf)%M(SH_Index(l,m)) &
          - qi * Rlm_tmp(SH_Index(l,-m))
      enddo
    enddo
    
    ! Add new contribution to new leaf (using WRAPPED position)
    dx = xw - self%nodes(newLeaf)%center(1)
    dy = yw - self%nodes(newLeaf)%center(2)
    dz = zw - self%nodes(newLeaf)%center(3)
    call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm_tmp)
    do l = 0, self%expansionOrder
      do m = -l, l
        self%nodes(newLeaf)%M(SH_Index(l,m)) = self%nodes(newLeaf)%M(SH_Index(l,m)) &
          + qi * Rlm_tmp(SH_Index(l,-m))
      enddo
    enddo
    
    self%particleToLeaf(iAtom) = newLeaf
    
    if (newLeaf /= oldLeaf) then
      call MoveParticleBetweenLeaves(self%nodes(oldLeaf), self%nodes(newLeaf), iAtom)
    endif
  enddo
  
  self%hasPendingMove = .true.
  self%pendingAtomIdx = disp(1)%atmIndx
  self%pendingOldX = atoms(1, disp(1)%atmIndx)
  self%pendingOldY = atoms(2, disp(1)%atmIndx)
  self%pendingOldZ = atoms(3, disp(1)%atmIndx)
end subroutine
!=============================================================================
subroutine MoveParticleBetweenLeaves(oldNode, newNode, iAtom)
  ! Helper to move particle between leaf lists when crossing cells
  ! (minimal impl for near-field direct sums)
  implicit none
  type(OctreeNode), intent(inout) :: oldNode, newNode
  integer, intent(in) :: iAtom
  integer :: n, i, pos
  integer, allocatable :: temp(:)

  ! Remove from old leaf list
  if (allocated(oldNode%particles)) then
    n = size(oldNode%particles)
    pos = 0
    do i = 1, n
      if (oldNode%particles(i) == iAtom) pos = i
    enddo
    if (pos > 0) then
      if (n > 1) then
        allocate(temp(n-1))
        temp(1:pos-1) = oldNode%particles(1:pos-1)
        temp(pos:n-1) = oldNode%particles(pos+1:n)
        call move_alloc(temp, oldNode%particles)
      else
        deallocate(oldNode%particles)
      endif
    endif
  endif

  ! Add to new leaf list (duplicate simple add)
  if (.not. allocated(newNode%particles)) then
    allocate(newNode%particles(1))
    newNode%particles(1) = iAtom
  else
    n = size(newNode%particles)
    allocate(temp(n+1))
    temp(1:n) = newNode%particles
    temp(n+1) = iAtom
    call move_alloc(temp, newNode%particles)
  endif

end subroutine
!=============================================================================
subroutine VerifyMultipoles_FMM(self, curbox)
  use ParallelVar, only: nout
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  
  integer :: iLeaf, nodeIdx, iAtom, atmType, l, m
  real(dp) :: qi, dx, dy, dz, maxErr
  real(dp), pointer :: atoms(:,:) => null()
  complex(dp) :: Rlm(self%nExpansion)
  complex(dp), allocatable :: M_scratch(:,:)
  integer, save :: callCount = 0
  integer :: leafForAtom, nBadLeaves
  
  callCount = callCount + 1
  if (callCount > 10) return
  
  call curbox%GetCoordinates(atoms)
  allocate(M_scratch(self%nExpansion, self%nNodes))
  M_scratch = (0.0_dp, 0.0_dp)
  
  ! Recompute M from scratch using particleToLeaf mapping (canonical source)
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    leafForAtom = self%particleToLeaf(iAtom)
    if (leafForAtom <= 0 .or. leafForAtom > self%nNodes) cycle
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    if (abs(qi) < 1.0E-15_dp) cycle
    dx = atoms(1, iAtom) - self%nodes(leafForAtom)%center(1)
    dy = atoms(2, iAtom) - self%nodes(leafForAtom)%center(2)
    dz = atoms(3, iAtom) - self%nodes(leafForAtom)%center(3)
    call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm)
    do l = 0, self%expansionOrder
      do m = -l, l
        M_scratch(SH_Index(l,m), leafForAtom) = M_scratch(SH_Index(l,m), leafForAtom) &
          + qi * Rlm(SH_Index(l,-m))
      enddo
    enddo
  enddo
  
  maxErr = 0.0_dp
  nBadLeaves = 0
  do iLeaf = 1, self%nLeaves
    nodeIdx = self%leafNodes(iLeaf)
    do l = 0, self%expansionOrder
      do m = -l, l
        maxErr = max(maxErr, abs(self%nodes(nodeIdx)%M(SH_Index(l,m)) - M_scratch(SH_Index(l,m), nodeIdx)))
      enddo
    enddo
    if (abs(self%nodes(nodeIdx)%M(1) - M_scratch(1, nodeIdx)) > 1.0E-8_dp) then
      nBadLeaves = nBadLeaves + 1
      if (nBadLeaves <= 3) then
        write(nout,'(A,I4,A,ES14.6,A,ES14.6)') " Bad leaf=", nodeIdx, &
          " M(1)stored=", real(self%nodes(nodeIdx)%M(1)), " M(1)scratch=", real(M_scratch(1, nodeIdx))
      endif
    endif
  enddo
  
  write(nout,'(A,I4,A,ES12.4,A,I4)') " M_VERIFY call=", callCount, &
    " maxErr=", maxErr, " badLeaves=", nBadLeaves
  
  deallocate(M_scratch)
end subroutine
!=============================================================================
function FindLeafForPosition_FMM(self, x, y, z) result(leafIdx)
  ! Find which leaf node contains the given position
  implicit none
  class(Pair_FMM), intent(in) :: self
  real(dp), intent(in) :: x, y, z
  integer :: leafIdx

  integer :: nodeIdx, octant
  real(dp) :: pos(3)
  
  pos = [x, y, z]
  nodeIdx = self%rootNode
  
  do while (.not. self%nodes(nodeIdx)%isLeaf)
    octant = GetOctant(pos, self%nodes(nodeIdx)%center)
    if (self%nodes(nodeIdx)%children(octant) > 0) then
      nodeIdx = self%nodes(nodeIdx)%children(octant)
    else
      exit
    endif
  enddo
  
  leafIdx = nodeIdx

contains

  pure function GetOctant(pos, center) result(octant)
    real(dp), intent(in) :: pos(3), center(3)
    integer :: octant
    
    octant = 1
    if (pos(1) >= center(1)) octant = octant + 1
    if (pos(2) >= center(2)) octant = octant + 2
    if (pos(3) >= center(3)) octant = octant + 4
    
  end function

end function
!=============================================================================
subroutine GetParticleEnergy_FMM(self, curbox, iAtom, x, y, z, E, accept)
  ! Get electrostatic energy of a single particle at position (x,y,z)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  integer, intent(in) :: iAtom
  real(dp), intent(in) :: x, y, z
  real(dp), intent(out) :: E
  logical, intent(out) :: accept

  integer :: jAtom, atmType1, atmType2
  integer :: leafIdx, k, j, nodeJ
  real(dp) :: qi, qj, rx, ry, rz, rsq, r
  real(dp) :: dx, dy, dz
  real(dp), pointer :: atoms(:,:) => null()
  complex(dp) :: Rlm(self%nExpansion)
  complex(dp) :: phi
  integer :: l, m, idx
  
  call curbox%GetCoordinates(atoms)
  
  E = 0.0_dp
  accept = .true.
  
  atmType1 = curbox%AtomType(iAtom)
  qi = self%q(atmType1)
  
  if (abs(qi) < 1.0E-15_dp) return
  
  ! Find leaf for this position
  leafIdx = self%FindLeafForPosition(x, y, z)
  
  ! Near-field: direct sum to own leaf + nearList cells (full pair energy)
  ! Own leaf
  if (allocated(self%nodes(leafIdx)%particles)) then
    do j = 1, size(self%nodes(leafIdx)%particles)
      jAtom = self%nodes(leafIdx)%particles(j)
      if (jAtom == iAtom) cycle
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      atmType2 = curbox%AtomType(jAtom)
      qj = self%q(atmType2)
      rx = x - atoms(1, jAtom)
      ry = y - atoms(2, jAtom)
      rz = z - atoms(3, jAtom)
      if (self%isPeriodic) then
        call curbox%Boundary(rx, ry, rz)
      endif
      rsq = rx*rx + ry*ry + rz*rz
      if (rsq < self%rMinTable(atmType1, atmType2)) then
        accept = .false.
        return
      endif
      r = sqrt(rsq)
      E = E + qi * qj * coulombConst / r
    enddo
  endif
  ! Near cells
  do k = 1, self%nodes(leafIdx)%nNear
    nodeJ = self%nodes(leafIdx)%nearList(k)
    if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
    do j = 1, size(self%nodes(nodeJ)%particles)
      jAtom = self%nodes(nodeJ)%particles(j)
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      atmType2 = curbox%AtomType(jAtom)
      qj = self%q(atmType2)
      rx = x - atoms(1, jAtom)
      ry = y - atoms(2, jAtom)
      rz = z - atoms(3, jAtom)
      if (self%isPeriodic) then
        call curbox%Boundary(rx, ry, rz)
      endif
      rsq = rx*rx + ry*ry + rz*rz
      if (rsq < self%rMinTable(atmType1, atmType2)) then
        accept = .false.
        return
      endif
      r = sqrt(rsq)
      E = E + qi * qj * coulombConst / r
    enddo
  enddo
  
  ! Far-field: direct sum over farList cells (full, consistent with Detailed/ComputeFarField direct)
  do k = 1, self%nodes(leafIdx)%nFar
    nodeJ = self%nodes(leafIdx)%farList(k)
    if (.not. allocated(self%nodes(nodeJ)%particles)) cycle
    do j = 1, size(self%nodes(nodeJ)%particles)
      jAtom = self%nodes(nodeJ)%particles(j)
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      atmType2 = curbox%AtomType(jAtom)
      qj = self%q(atmType2)
      rx = x - atoms(1, jAtom)
      ry = y - atoms(2, jAtom)
      rz = z - atoms(3, jAtom)
      if (self%isPeriodic) then
        call curbox%Boundary(rx, ry, rz)
      endif
      rsq = rx*rx + ry*ry + rz*rz
      if (rsq < self%rMinTable(atmType1, atmType2)) then
        accept = .false.
        return
      endif
      r = sqrt(rsq)
      E = E + qi * qj * coulombConst / r
    enddo
  enddo

end subroutine
!=============================================================================
subroutine ComputePairDelta_FMM(self, curbox, atoms, iAtom, qi, jAtom, xnew, ynew, znew, E_Diff, accept)
  ! Helper for Shift: compute E_pair delta (old vs new pos for i-j), add to E_Diff and both dETable
  ! Ensures table consistency (individuals + sum) and no drift in tablecheck.
  implicit none
  class(Pair_FMM), intent(in) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), pointer :: atoms(:,:)
  integer, intent(in) :: iAtom, jAtom
  real(dp), intent(in) :: qi, xnew, ynew, znew
  real(dp), intent(inout) :: E_Diff
  logical, intent(inout) :: accept
  integer :: atmType2
  real(dp) :: qj, rx, ry, rz, rsq, r, E_old_pair, E_new_pair, delta

  atmType2 = curbox%AtomType(jAtom)
  qj = self%q(atmType2)
  if (abs(qj) < 1.0E-15_dp) return

  ! Old pair (current positions)
  rx = atoms(1, iAtom) - atoms(1, jAtom)
  ry = atoms(2, iAtom) - atoms(2, jAtom)
  rz = atoms(3, iAtom) - atoms(3, jAtom)
  if (self%isPeriodic) then
    call curbox%Boundary(rx, ry, rz)
  endif
  rsq = rx*rx + ry*ry + rz*rz
  r = sqrt(rsq)
  E_old_pair = qi * qj * coulombConst / r

  ! New pair at trial pos
  rx = xnew - atoms(1, jAtom)
  ry = ynew - atoms(2, jAtom)
  rz = znew - atoms(3, jAtom)
  if (self%isPeriodic) then
    call curbox%Boundary(rx, ry, rz)
  endif
  rsq = rx*rx + ry*ry + rz*rz
  if (rsq < self%rMinTable(curbox%AtomType(iAtom), atmType2)) then
    accept = .false.
    return
  endif
  r = sqrt(rsq)
  E_new_pair = qi * qj * coulombConst / r

  delta = E_new_pair - E_old_pair
  E_Diff = E_Diff + delta
  curbox%dETable(iAtom) = curbox%dETable(iAtom) + delta
  curbox%dETable(jAtom) = curbox%dETable(jAtom) + delta
end subroutine
!=============================================================================
subroutine ProcessIO_FMM(self, line)
  use Common_MolInfo, only: nAtomTypes
  use Input_Format, only: CountCommands, GetXCommand
  use Input_Format, only: maxLineLen
  use Units, only: inEngUnit, inLenUnit
  implicit none
  class(Pair_FMM), intent(inout) :: self
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: command
  integer :: jType, lineStat, nPar
  integer :: type1, intVal
  real(dp) :: rCut, q, rMin, realVal

  call GetXCommand(line, command, 1, lineStat)
  call CountCommands(line, nPar)

  select case(trim(adjustl(command)))
    case("rcut")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) rCut
      self%rCut = rCut * inLenUnit
      self%rCutSq = self%rCut**2

    case("p", "order", "expansion_order")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) intVal
      self%expansionOrder = intVal
      self%nExpansion = SH_GetExpansionSize(intVal)
      call SH_Init(intVal * 2)

    case("maxlevel", "max_level", "depth")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) intVal
      self%maxLevel = intVal

    case("ncrit", "n_crit")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) intVal
      self%nCrit = intVal

    case("periodic")
      call GetXCommand(line, command, 2, lineStat)
      select case(trim(adjustl(command)))
        case("true", "yes", "1", ".true.")
          self%isPeriodic = .true.
        case("false", "no", "0", ".false.")
          self%isPeriodic = .false.
      end select

    case("usemultipole", "use_multipole", "multipole")
      call GetXCommand(line, command, 2, lineStat)
      select case(trim(adjustl(command)))
        case("true", "yes", "1", ".true.")
          self%useMultipole = .true.
        case("false", "no", "0", ".false.")
          self%useMultipole = .false.
      end select

    case("precision")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) self%fmmPrecision

    case("charge", "q")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) q
      
      self%q(type1) = q
      do jType = 1, nAtomTypes
        self%qTable(type1, jType) = q * self%q(jType) * coulombConst
        self%qTable(jType, type1) = q * self%q(jType) * coulombConst
      enddo

    case("rmin")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) rMin
      
      rMin = rMin * inLenUnit
      self%rMin(type1) = rMin
      do jType = 1, nAtomTypes
        self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
        self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
      enddo

    case default
      ! Try to parse as: type q rMin
      if (nPar >= 3) then
        read(line, *, iostat=lineStat) type1, q, rMin
        if (lineStat == 0) then
          rMin = rMin * inLenUnit
          self%q(type1) = q
          self%rMin(type1) = rMin
          
          do jType = 1, nAtomTypes
            self%qTable(type1, jType) = q * self%q(jType) * coulombConst
            self%qTable(jType, type1) = q * self%q(jType) * coulombConst
            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        endif
      endif
  end select

end subroutine
!=============================================================================
subroutine Prologue_FMM(self)
  use ParallelVar, only: nout
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_FMM), intent(inout) :: self
  integer :: i

  write(nout, *) "============================================"
  write(nout, *) "Fast Multipole Method (FMM) Parameters:"
  write(nout, *) "  Expansion order (p):", self%expansionOrder
  write(nout, *) "  Maximum tree depth:", self%maxLevel
  write(nout, *) "  Particles per leaf (nCrit):", self%nCrit
  write(nout, *) "  Near-field cutoff:", self%rCut
  write(nout, *) "  Periodic boundaries:", self%isPeriodic
  write(nout, *) "  Use multipole far-field:", self%useMultipole
  write(nout, *) "  Target precision:", self%fmmPrecision
  write(nout, *) ""
  write(nout, *) "Charges per atom type:"
  do i = 1, nAtomTypes
    write(nout, *) "  Type", i, ":", self%q(i)
  enddo
  write(nout, *) "============================================"

end subroutine
!=============================================================================
subroutine Epilogue_FMM(self)
  implicit none
  class(Pair_FMM), intent(inout) :: self

  ! Cleanup spherical harmonics
  call SH_Cleanup()

end subroutine
!=============================================================================
function GetCutOff_FMM(self) result(rCut)
  implicit none
  class(Pair_FMM), intent(inout) :: self
  real(dp) :: rCut

  rCut = self%rCut

end function
!=============================================================================
subroutine VerifyAdditionTheorem(p)
  use ParallelVar, only: nout
  use SphericalHarmonics, only: SH_ComputeSolidHarmonic, SH_ComputeIrregularSolid, SH_Index
  use VarPrecision
  implicit none
  integer, intent(in) :: p

  real(dp) :: a(3), b(3), ab(3), r_ab
  real(dp) :: exact, approx_noSign, approx_negm
  complex(dp), allocatable :: Rlm(:), Slm(:)
  integer :: l, m, idx_lm, idx_lnm, nR, nS
  complex(dp) :: sum1, sum2

  nR = (p+1)*(p+1)
  nS = (2*p+1)*(2*p+1)
  allocate(Rlm(nR), Slm(nS))

  a = [8.0_dp, 5.0_dp, 3.0_dp]
  b = [1.0_dp, 0.5_dp, -0.3_dp]
  ab = a - b
  r_ab = sqrt(ab(1)**2 + ab(2)**2 + ab(3)**2)
  exact = 1.0_dp / r_ab

  call SH_ComputeSolidHarmonic(b(1), b(2), b(3), p, Rlm)
  call SH_ComputeIrregularSolid(a(1), a(2), a(3), 2*p, Slm)

  ! Test: sum_{l,m} R_l^{-m}(b) S_l^m(a) should = 1/|a-b|
  sum1 = (0.0_dp, 0.0_dp)
  ! Test: sum_{l,m} (-1)^m R_l^{-m}(b) S_l^m(a)
  sum2 = (0.0_dp, 0.0_dp)
  do l = 0, p
    do m = -l, l
      idx_lm = SH_Index(l, m)
      idx_lnm = SH_Index(l, -m)
      sum1 = sum1 + Rlm(idx_lnm) * Slm(idx_lm)
      if (mod(abs(m), 2) == 0) then
        sum2 = sum2 + Rlm(idx_lnm) * Slm(idx_lm)
      else
        sum2 = sum2 - Rlm(idx_lnm) * Slm(idx_lm)
      endif
    enddo
  enddo

  write(nout,*) "=== Addition Theorem Verification ==="
  write(nout,*) "  Exact 1/|a-b|:             ", exact
  write(nout,*) "  sum R(-m) S(m) (no sign):  ", real(sum1, dp), " err:", abs(real(sum1,dp)-exact)/exact*100, "%"
  write(nout,*) "  sum (-1)^m R(-m) S(m):     ", real(sum2, dp), " err:", abs(real(sum2,dp)-exact)/exact*100, "%"

  ! Also test: sum R_l^m(b) S_l^m(a) (same sign for both)
  sum1 = (0.0_dp, 0.0_dp)
  do l = 0, p
    do m = -l, l
      idx_lm = SH_Index(l, m)
      sum1 = sum1 + Rlm(idx_lm) * Slm(idx_lm)
    enddo
  enddo
  write(nout,*) "  sum R(m) S(m):             ", real(sum1, dp), " err:", abs(real(sum1,dp)-exact)/exact*100, "%"

  ! Test with R(-m) (should equal kernel for correct convention)
  sum1 = (0.0_dp, 0.0_dp)
  do l = 0, p
    do m = -l, l
      idx_lm = SH_Index(l, m)
      idx_lnm = SH_Index(l, -m)
      sum1 = sum1 + Rlm(idx_lnm) * Slm(idx_lm)
    enddo
  enddo
  write(nout,*) "  sum R(-m) S(m) (dup check):", real(sum1, dp), " err:", abs(real(sum1,dp)-exact)/exact*100, "%"

  ! Verify the S addition theorem: S_n^m(a+b) = sum_{j,k} C * R_j^k(a) * S_{n+j}^{m-k}(b)
  ! Test with S_0^0(delta+d) = 1/|delta+d|
  block
    real(dp) :: delta(3), d(3), apb(3), r_apb
    real(dp) :: exact_S00
    complex(dp), allocatable :: R_delta(:), S_d(:), S_apb(:)
    complex(dp) :: sum_nj, sum_no_sign
    integer :: j2, k2, idx_jk2, njp2, mmk2, idx_s2

    delta = [0.3_dp, -0.2_dp, 0.1_dp]
    d = [8.0_dp, 5.0_dp, 3.0_dp]
    apb = delta + d
    r_apb = sqrt(apb(1)**2 + apb(2)**2 + apb(3)**2)
    exact_S00 = 1.0_dp / r_apb

    allocate(R_delta(nR), S_d(nS), S_apb(nS))
    call SH_ComputeSolidHarmonic(delta(1), delta(2), delta(3), p, R_delta)
    call SH_ComputeIrregularSolid(d(1), d(2), d(3), 2*p, S_d)
    call SH_ComputeIrregularSolid(apb(1), apb(2), apb(3), 2*p, S_apb)

    ! Test S_0^0(delta+d) = sum_{j,k} C * R_j^k(delta) * S_j^{-k}(d)
    ! Try C = (-1)^j
    sum_nj = (0.0_dp, 0.0_dp)
    sum_no_sign = (0.0_dp, 0.0_dp)
    do j2 = 0, p
      do k2 = -j2, j2
        idx_jk2 = SH_Index(j2, k2)
        mmk2 = -k2
        if (abs(mmk2) > j2) cycle
        idx_s2 = SH_Index(j2, mmk2)
        if (mod(j2, 2) == 0) then
          sum_nj = sum_nj + R_delta(idx_jk2) * S_d(idx_s2)
        else
          sum_nj = sum_nj - R_delta(idx_jk2) * S_d(idx_s2)
        endif
        sum_no_sign = sum_no_sign + R_delta(idx_jk2) * S_d(idx_s2)
      enddo
    enddo

    write(nout,*) "--- S addition theorem: S_0^0(delta+d) ---"
    write(nout,*) "  Exact S_0^0:                ", exact_S00
    write(nout,*) "  Direct S_0^0(delta+d):      ", real(S_apb(SH_Index(0,0)), dp)
    write(nout,*) "  sum (-1)^j R S (n=0,m=0):   ", real(sum_nj, dp), " err:", &
      abs(real(sum_nj,dp) - exact_S00)/exact_S00*100, "%"
    write(nout,*) "  sum R S no sign (n=0,m=0):  ", real(sum_no_sign, dp), " err:", &
      abs(real(sum_no_sign,dp) - exact_S00)/exact_S00*100, "%"

    ! Verify R addition theorem (should be EXACT since R is polynomial)
    block
      complex(dp), allocatable :: R_a(:), R_b(:), R_apb(:)
      complex(dp) :: sum_R, exact_R
      integer :: n_t, m_t, j_t, k_t, nj_t, mk_t
      real(dp) :: max_err_R
      integer :: worst_n, worst_m

      allocate(R_a(nR), R_b(nR), R_apb(nR))
      call SH_ComputeSolidHarmonic(delta(1), delta(2), delta(3), p, R_a)
      call SH_ComputeSolidHarmonic(d(1), d(2), d(3), p, R_b)
      call SH_ComputeSolidHarmonic(apb(1), apb(2), apb(3), p, R_apb)

      write(nout,*) "--- R addition theorem test (should be exact) ---"
      max_err_R = 0.0_dp
      worst_n = 0; worst_m = 0
      do n_t = 0, p
        do m_t = -n_t, n_t
          exact_R = R_apb(SH_Index(n_t, m_t))
          sum_R = (0.0_dp, 0.0_dp)
          do j_t = 0, n_t
            do k_t = -j_t, j_t
              nj_t = n_t - j_t
              mk_t = m_t - k_t
              if (abs(mk_t) > nj_t) cycle
              sum_R = sum_R + R_a(SH_Index(j_t, k_t)) * R_b(SH_Index(nj_t, mk_t))
            enddo
          enddo
          if (abs(exact_R) > 1e-15_dp) then
            if (abs(sum_R - exact_R)/abs(exact_R)*100 > max_err_R) then
              max_err_R = abs(sum_R - exact_R)/abs(exact_R)*100
              worst_n = n_t; worst_m = m_t
            endif
          endif
          if (n_t <= 2) then
            write(nout,'(A,I2,A,I3,A,2ES12.4,A,2ES12.4,A,ES10.2,A)') &
              "  R(", n_t, ",", m_t, ") exact=", real(exact_R), aimag(exact_R), &
              "  add=", real(sum_R), aimag(sum_R), &
              "  err=", merge(abs(sum_R-exact_R)/abs(exact_R)*100, 0.0_dp, abs(exact_R)>1e-15_dp), "%"
          endif
        enddo
      enddo
      write(nout,'(A,ES10.2,A,I2,A,I3,A)') "  Max R add.thm err: ", max_err_R, &
        "% at R(", worst_n, ",", worst_m, ")"
      deallocate(R_a, R_b, R_apb)
    end block

    ! Test S addition theorem for various (n,m) including m≠0
    block
      complex(dp) :: sum_S, exact_S_val
      integer :: n_s, m_s, j_s, k_s, njp_s, mmk_s
      write(nout,*) "--- S addition theorem for S(n,m) ---"
      do n_s = 0, 2
        do m_s = -n_s, n_s
          exact_S_val = S_apb(SH_Index(n_s, m_s))
          sum_S = (0.0_dp, 0.0_dp)
          do j_s = 0, p
            do k_s = -j_s, j_s
              njp_s = n_s + j_s
              mmk_s = m_s - k_s
              if (abs(mmk_s) > njp_s .or. njp_s > 2*p) cycle
              if (mod(j_s, 2) == 0) then
                sum_S = sum_S + R_delta(SH_Index(j_s, k_s)) * S_d(SH_Index(njp_s, mmk_s))
              else
                sum_S = sum_S - R_delta(SH_Index(j_s, k_s)) * S_d(SH_Index(njp_s, mmk_s))
              endif
            enddo
          enddo
          write(nout,'(A,I2,A,I3,A,2ES14.6,A,2ES14.6,A,ES10.2,A)') &
            "  S(", n_s, ",", m_s, ") exact=", real(exact_S_val), aimag(exact_S_val), &
            "  add=", real(sum_S), aimag(sum_S), &
            "  err=", abs(sum_S - exact_S_val)/max(abs(exact_S_val), 1d-30)*100, "%"
        enddo
      enddo
    end block

    ! Direct S values at delta+d
    block
      real(dp) :: r_pt(3), y_pt(3), dist, exact_val
      complex(dp), allocatable :: R_r(:), S_y(:)
      complex(dp) :: sum_direct
      integer :: lt, mt

      allocate(R_r(nR), S_y(nS))
      r_pt = apb  ! delta + d
      y_pt = [0.0_dp, 0.0_dp, 0.0_dp]

      dist = sqrt(sum((r_pt - y_pt)**2))
      exact_val = 1.0_dp / dist

      call SH_ComputeSolidHarmonic(y_pt(1), y_pt(2), y_pt(3), p, R_r)
      call SH_ComputeIrregularSolid(r_pt(1), r_pt(2), r_pt(3), 2*p, S_y)

      write(nout,*) "--- Direct S values at delta+d ---"
      do lt = 0, min(2, p)
        do mt = -lt, lt
          write(nout,'(A,I2,A,I3,A,2ES14.6)') "  S(", lt, ",", mt, ") = ", &
            real(S_y(SH_Index(lt,mt))), aimag(S_y(SH_Index(lt,mt)))
        enddo
      enddo
      deallocate(R_r, S_y)
    end block

    deallocate(R_delta, S_d, S_apb)
  end block
  write(nout,*) "======================================="

  deallocate(Rlm, Slm)
end subroutine
!=============================================================================
subroutine TestM2L_SinglePair(self, curbox)
  use ParallelVar, only: nout
  implicit none
  class(Pair_FMM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  integer :: p
  real(dp) :: q_test, ys(3), rs(3), ds(3), exact_phi, phi_val
  complex(dp), allocatable :: M_test(:), L_test(:), Rlm(:), Slm(:)
  complex(dp) :: phi
  integer :: l, m, n, m_n, j, k, idx, njp, mmk, nExp
  integer :: idx_jk, idx_nm, idx_s
  real(dp) :: dx, dy, dz, sign_factor, r_exact

  p = self%expansionOrder
  nExp = (p+1)*(p+1)
  allocate(M_test(nExp), L_test(nExp), Rlm(nExp), Slm((2*p+1)*(2*p+1)))

  q_test = 1.0_dp
  ys = [1.0_dp, 0.5_dp, -0.3_dp]
  ds = [8.0_dp, 5.0_dp, 3.0_dp]
  rs = ds + [0.2_dp, -0.1_dp, 0.15_dp]

  ! Exact potential at rs from charge q_test at ys
  dx = rs(1) - ys(1)
  dy = rs(2) - ys(2)
  dz = rs(3) - ys(3)
  r_exact = sqrt(dx*dx + dy*dy + dz*dz)
  exact_phi = q_test / r_exact

  ! P2M: source center at origin
  M_test = (0.0_dp, 0.0_dp)
  call SH_ComputeSolidHarmonic(ys(1), ys(2), ys(3), p, Rlm)
  do l = 0, p
    do m = -l, l
      M_test(SH_Index(l,m)) = q_test * Rlm(SH_Index(l,-m))
    enddo
  enddo

  ! Try different sign conventions for M2L
  ! Convention A: (-1)^{m+k}  (current)
  L_test = (0.0_dp, 0.0_dp)
  call SH_ComputeIrregularSolid(ds(1), ds(2), ds(3), 2*p, Slm)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(abs(m_n + k), 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "=== M2L SIGN TEST ==="
  write(nout,*) "  Exact phi:", exact_phi
  write(nout,*) "  (-1)^{m+k}:", real(phi, dp)

  ! Convention B: (-1)^{j-k}  (original code)
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(abs(j - k), 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^{j-k}:", real(phi, dp)

  ! Convention C: (-1)^n
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(n, 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^n:", real(phi, dp)

  ! Convention D: (-1)^{n+j-k}
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(abs(n + j - k), 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^{n+j-k}:", real(phi, dp)

  ! Convention E: (-1)^j (correct for CS-phase Legendre convention)
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(j, 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^j:", real(phi, dp)

  ! Convention E2: (-1)^{j-k} with conj(S) i.e. index k-m instead of m-k
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = k - m_n
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(abs(j - k), 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^{j-k} conj(S):", real(phi, dp)

  ! Convention F: (-1)^n with conj(S)
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = k - m_n
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(n, 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^n conj(S):", real(phi, dp)

  ! Convention G: (-1)^{m+j-k} with original S (m-k)
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(abs(m_n + j - k), 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^{m+j-k}:", real(phi, dp)

  ! Convention H: (-1)^n with original S (m-k), include (-1)^m in front
  L_test = (0.0_dp, 0.0_dp)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = 0, p
        do m_n = -n, n
          njp = n + j; mmk = m_n - k
          if (abs(mmk) > njp .or. njp > 2*p) cycle
          if (mod(abs(n + m_n), 2) == 0) then; sign_factor = 1.0_dp
          else; sign_factor = -1.0_dp; endif
          idx_nm = SH_Index(n, m_n); idx_s = SH_Index(njp, mmk)
          L_test(idx_jk) = L_test(idx_jk) + sign_factor * M_test(idx_nm) * Slm(idx_s)
        enddo
      enddo
    enddo
  enddo
  dx = rs(1) - ds(1); dy = rs(2) - ds(2); dz = rs(3) - ds(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  phi = (0.0_dp, 0.0_dp)
  do l = 0, p; do m = -l, l
    idx = SH_Index(l, m); phi = phi + L_test(idx) * Rlm(idx)
  enddo; enddo
  write(nout,*) "  (-1)^{n+m}:", real(phi, dp)
  write(nout,*) "======================"

  deallocate(M_test, L_test, Rlm, Slm)
end subroutine
!=============================================================================
end module FF_FMM
!=============================================================================
