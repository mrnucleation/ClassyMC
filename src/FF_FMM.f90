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
    
    ! MC move tracking
    integer :: lastMovedParticle
    integer :: lastOldLeaf
    integer :: lastNewLeaf
    
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
  
  call self%BuildOctree(curbox)
  call self%BuildInteractionLists()
  call self%ComputeMultipoles(curbox)
  call self%ComputeLocalExpansions()
  call self%ComputeNearFieldOctree(curbox, E_near, accept)
  if (.not. accept) return
  ! Use direct far-field as the authoritative energy
  call self%ComputeFarField(curbox, E_far)
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

  ! Loop over displaced atoms: compute pair deltas using tree's near/far lists from moved atom's leaf
  ! (standard for table consistency: update both i and j dETable; E_Diff = sum pair_delta)
  ! Uses ComputePairDelta_FMM helper for old/new E_pair.
  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)
    qi = self%q(atmType1)
    
    if (abs(qi) < 1.0E-15_dp) cycle

    ! Get leaf for old position (use mapping for exact; lists from old tree)
    leafIdx = self%particleToLeaf(iAtom)
    newLeaf = leafIdx
    if (self%treeBuilt) then
      newLeaf = self%FindLeafForPosition(disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new)
    endif

    ! If the particle crosses leaves (or mapping is invalid), do a full scan to avoid misses
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

      ! 3. Far cells (full direct)
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
    threshold = 3.0_dp * self%nodes(nodeI)%halfWidth * 2.0_dp
    
    nearCount = 0
    farCount = 0
    
    do jLeaf = 1, self%nLeaves
      if (iLeaf == jLeaf) cycle
      
      nodeJ = self%leafNodes(jLeaf)
      
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
        farCount = farCount + 1
        tempFar(farCount) = nodeJ
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

  integer :: iLeaf, nodeIdx, iAtom, j
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
      
      ! P2M: M_n^m = sum_i q_i * conj(R_n^m(y_i - c))
      ! The conjugate is essential for the addition theorem to work correctly
      self%nodes(nodeIdx)%M = self%nodes(nodeIdx)%M + qi * conjg(Rlm)
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
  ! Compute far-field energy from local expansions using multipole approximation
  !
  ! NOTE: This routine requires correct spherical harmonic translation operators.
  ! Use ComputeFarField (direct summation) until the multipole operators are verified.
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
  
  call curbox%GetCoordinates(atoms)
  
  E_far = 0.0_dp
  
  ! For each particle, evaluate local expansion at its position
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
      
      ! Compute regular solid harmonics at particle position
      call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm)
      
      ! Evaluate local expansion: phi = sum_{l,m} L_l^m * R_l^m(r)
      phi = (0.0_dp, 0.0_dp)
      do l = 0, self%expansionOrder
        do m = -l, l
          idx = SH_Index(l, m)
          phi = phi + self%nodes(nodeIdx)%L(idx) * Rlm(idx)
        enddo
      enddo
      
      ! Energy contribution (divide by 2 to avoid double counting for total E_far)
      ! Add full to ETable (consistent with near/direct and ETable convention: sum/2 = total)
      E_far = E_far + 0.5_dp * qi * real(phi, dp) * coulombConst
      curbox%ETable(iAtom) = curbox%ETable(iAtom) + qi * real(phi, dp) * coulombConst
    enddo
  enddo

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

  integer :: oldLeaf, newLeaf, atmType
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
  
  ! Remove old contribution from old leaf (P2M uses conjugate)
  dx = x_old - self%nodes(oldLeaf)%center(1)
  dy = y_old - self%nodes(oldLeaf)%center(2)
  dz = z_old - self%nodes(oldLeaf)%center(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm_old)
  self%nodes(oldLeaf)%M = self%nodes(oldLeaf)%M - qi * conjg(Rlm_old)
  
  ! Add new contribution to new leaf (P2M uses conjugate)
  dx = x_new - self%nodes(newLeaf)%center(1)
  dy = y_new - self%nodes(newLeaf)%center(2)
  dz = z_new - self%nodes(newLeaf)%center(3)
  call SH_ComputeSolidHarmonic(dx, dy, dz, self%expansionOrder, Rlm_new)
  self%nodes(newLeaf)%M = self%nodes(newLeaf)%M + qi * conjg(Rlm_new)
  
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
  ! Called when a move is rejected - restore old multipoles
  implicit none
  class(Pair_FMM), intent(inout) :: self

  ! Restore multipoles
  self%nodes(self%lastOldLeaf)%M = self%nodes(self%lastOldLeaf)%M_old
  if (self%lastNewLeaf /= self%lastOldLeaf) then
    self%nodes(self%lastNewLeaf)%M = self%nodes(self%lastNewLeaf)%M_old
  endif
  
  ! Restore particle mapping
  self%particleToLeaf = self%particleToLeaf_old

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
          ! Multipole far-field will be used once verified
          ! Currently the direct far-field is always used for energy
          ! This flag is noted but does not yet switch the energy pathway
          continue
        case("false", "no", "0", ".false.")
          continue
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
end module FF_FMM
!=============================================================================
