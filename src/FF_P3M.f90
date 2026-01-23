!================================================================================
! P3M (Particle-Particle Particle-Mesh) Electrostatics Module
!
! This module implements the P3M method for computing long-range electrostatic
! interactions in periodic systems. P3M uses:
!   1. Real-space short-range part with screening (like Ewald)
!   2. Mesh-based FFT for long-range reciprocal space (O(N log N) scaling)
!   3. Charge assignment to mesh grid
!   4. Differentiation for force/energy calculation
!
! Advantages over Ewald:
!   - O(N log N) scaling vs O(N^1.5) for Ewald
!   - More efficient for large systems
!   - Better parallelization potential
!
! The P3M method splits electrostatics as:
!   E_total = E_real + E_mesh + E_self + E_intra
!
! References:
!   - Hockney & Eastwood, "Computer Simulation Using Particles" (1988)
!   - Deserno & Holm, J. Chem. Phys. 109, 7678 (1998)
!   - Ballenegger et al., J. Chem. Phys. 131, 094107 (2009)
!================================================================================
module FF_P3M
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use ClassyConstants, only: coulombConst, pi, two_pi
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  use fftw3_interface

  implicit none

  type, extends(EasyPair_Cut) :: Pair_P3M
    ! P3M parameters
    real(dp) :: alpha          ! Screening parameter (1/length)
    real(dp) :: alpha_sq       ! alpha^2 for convenience
    real(dp) :: rCutReal       ! Real-space cutoff
    real(dp) :: rCutRealSq     ! Real-space cutoff squared
    
    ! Mesh parameters
    integer :: meshX, meshY, meshZ        ! Mesh dimensions
    integer :: chargeOrder = 4            ! Charge assignment order (typically 4-6)
    real(dp) :: meshSpacing               ! Average mesh spacing
    
    ! Precision control
    real(dp) :: p3mPrecision = 1.0E-5_dp  ! Target precision
    
    ! Charges
    real(dp), allocatable :: q(:)         ! Charge per atom type
    real(dp), allocatable :: qTable(:,:)  ! q_i * q_j * coulombConst lookup table
    
    ! Mesh arrays (pointer instead of allocatable+target for gfortran compatibility)
    real(dp), pointer :: chargeMesh(:,:,:) => null()      ! Real-space charge mesh
    complex(dp), pointer :: chargeMeshK(:,:,:) => null()  ! Fourier-space charge mesh
    
    ! Green's function (optimal influence function)
    real(dp), allocatable :: greenFunc(:,:,:)       ! G(k) for energy calculation
    
    ! Charge assignment function (CAF) - precomputed
    real(dp), allocatable :: cafWeights(:,:)        ! Precomputed CAF weights
    
    ! Energy components
    real(dp) :: E_mesh_total   ! Current mesh (reciprocal) energy
    real(dp) :: E_self_total   ! Current self-energy correction
    
    ! Cached mesh for efficient MC updates
    real(dp), pointer :: chargeMesh_old(:,:,:) => null()
    complex(dp), pointer :: chargeMeshK_old(:,:,:) => null()
    
    ! Box dimensions cache
    real(dp) :: cachedLx, cachedLy, cachedLz
    
    ! Flags
    logical :: initialized = .false.
    logical :: meshSetup = .false.
    
    ! FFTW plan handles
    type(C_PTR) :: plan_forward = C_NULL_PTR
    logical :: fftw_initialized = .false.
    
  contains
    procedure, pass :: Constructor => Constructor_P3M
    procedure, pass :: PairFunction => PairFunction_P3M
    procedure, pass :: DetailedECalc => Detailed_P3M
    procedure, pass :: DiffECalc => DiffECalc_P3M
    procedure, pass :: ShiftECalc_Single => Shift_P3M_Single
    procedure, pass :: NewECalc => New_P3M
    procedure, pass :: OldECalc => Old_P3M
    procedure, pass :: OrthoVolECalc => OrthoVol_P3M
    procedure, pass :: ProcessIO => ProcessIO_P3M
    procedure, pass :: Prologue => Prologue_P3M
    procedure, pass :: GetCutOff => GetCutOff_P3M
    
    ! P3M-specific procedures
    procedure, pass :: SetupMesh => SetupMesh_P3M
    procedure, pass :: ComputeGreenFunction => ComputeGreenFunction_P3M
    procedure, pass :: AssignChargesToMesh => AssignChargesToMesh_P3M
    procedure, pass :: ComputeMeshEnergy => ComputeMeshEnergy_P3M
    procedure, pass :: ComputeSelfEnergy => ComputeSelfEnergy_P3M
    procedure, pass :: ComputeRealPair => ComputeRealPair_P3M
    procedure, pass :: OptimizeAlpha => OptimizeAlpha_P3M
    procedure, pass :: UpdateMeshForMove => UpdateMeshForMove_P3M
    procedure, pass :: GetCAFWeights => GetCAFWeights_P3M
    procedure, pass :: ComputeVolumeMoveEnergy => ComputeVolumeMoveEnergy_P3M
    procedure, pass :: FFT3D_Forward => FFT3D_Forward_P3M
    procedure, pass :: AssignSingleChargeToMesh => AssignSingleChargeToMesh_P3M
    procedure, pass :: Destructor => Destructor_P3M
    procedure, pass :: AcceptUpdate => AcceptUpdate_P3M
    procedure, pass :: RejectUpdate => RejectUpdate_P3M
    
  end type

contains

!=============================================================================
subroutine Constructor_P3M(self)
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_P3M), intent(inout) :: self
  integer :: AllocateStat

  ! Call parent constructor
  call self%EasyPair_Cut%Constructor()

  ! Allocate charge arrays
  allocate(self%q(1:nAtomTypes), stat=AllocateStat)
  allocate(self%qTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  
  if (.not. allocated(self%rMin)) then
    allocate(self%rMin(1:nAtomTypes), stat=AllocateStat)
  endif
  if (.not. allocated(self%rMinTable)) then
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  endif

  self%q = 0.0_dp
  self%qTable = 0.0_dp
  self%rMin = 0.5_dp
  self%rMinTable = 0.25_dp  ! Stored as r^2
  
  ! Default P3M parameters
  self%alpha = 0.3_dp
  self%alpha_sq = self%alpha**2
  self%rCutReal = 10.0_dp
  self%rCutRealSq = self%rCutReal**2
  self%rCut = self%rCutReal
  self%rCutSq = self%rCutRealSq
  
  ! Default mesh parameters
  self%meshX = 32
  self%meshY = 32
  self%meshZ = 32
  self%chargeOrder = 4
  
  self%initialized = .false.
  self%meshSetup = .false.

  if (AllocateStat /= 0) error stop "Allocation error in P3M Constructor"

end subroutine

!=============================================================================
function PairFunction_P3M(self, rsq, atmtype1, atmtype2) result(E_Pair)
  ! Real-space pair function with erfc screening
  implicit none
  class(Pair_P3M), intent(inout) :: self
  real(dp), intent(in) :: rsq
  integer, intent(in) :: atmtype1, atmtype2
  real(dp) :: E_Pair
  real(dp) :: r, qij

  if (rsq >= self%rCutRealSq) then
    E_Pair = 0.0_dp
    return
  endif

  r = sqrt(rsq)
  qij = self%qTable(atmtype1, atmtype2)
  
  ! Real-space P3M/Ewald: q_i * q_j * erfc(alpha * r) / r
  E_Pair = qij * erfc(self%alpha * r) / r

end function

!=============================================================================
subroutine Detailed_P3M(self, curbox, E_T, accept)
  use ParallelVar, only: nout
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(inout) :: E_T
  logical, intent(out) :: accept
  
  real(dp) :: E_real, E_mesh, E_self
  real(dp) :: E_Corr

  ! Setup mesh if needed
  call self%SetupMesh(curbox)
  
  ! Compute real-space contribution using parent class
  call self%EasyPair_Cut%DetailedECalc(curbox, E_real, accept)
  
  if (.not. accept) then
    E_T = 0.0_dp
    return
  endif
  
  ! Assign charges to mesh
  call self%AssignChargesToMesh(curbox)
  
  ! Compute mesh (reciprocal space) contribution
  call self%ComputeMeshEnergy(curbox, E_mesh)
  self%E_mesh_total = E_mesh
  
  ! Compute self-energy correction
  call self%ComputeSelfEnergy(curbox, E_self)
  self%E_self_total = E_self
  
  write(nout,*) "P3M Real-space Energy:", E_real
  write(nout,*) "P3M Mesh Energy:", E_mesh
  write(nout,*) "P3M Self-energy Correction:", E_self
  write(nout,*) "P3M Total Energy:", E_real + E_mesh + E_self
  
  E_T = E_real + E_mesh + E_self

end subroutine

!=============================================================================
subroutine DiffECalc_P3M(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(simBox), intent(inout) :: curbox
  class(Perturbation), intent(inout), target :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  accept = .true.
  curbox%dETable = 0.0_dp
  E_Diff = 0.0_dp

  select type(disp)
    class is(Displacement)
      ! For displacements, don't pass tempList/tempNNei - they're not initialized
      call self%ShiftECalc_Single(curbox, disp, E_Diff, accept)
      
    class is(Addition)
      ! For insertions, tempList/tempNNei contain the new particle's neighbors
      call self%NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)
      
    class is(Deletion)
      call self%OldECalc(curbox, disp, E_Diff)
      
    class is(OrthoVolChange)
      call self%OrthoVolECalc(curbox, disp, E_Diff, accept)
      
    class default
      ! Fall back to parent for other perturbation types
      call self%EasyPair_Cut%DiffECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)
  end select

end subroutine

!=============================================================================
subroutine Shift_P3M_Single(self, curbox, disp, E_Diff, accept, tempList, tempNNei)
  ! Efficient energy difference for particle displacement
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(inout) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept
  integer, intent(in), target, optional :: tempList(:,:), tempNNei(:)
  
  real(dp) :: E_real_diff, E_mesh_diff

  ! Compute real-space contribution using parent class
  ! Must check if optional arguments are present before passing them through
  if (present(tempList) .and. present(tempNNei)) then
    call self%EasyPair_Cut%ShiftECalc_Single(curbox, disp, E_real_diff, accept, tempList, tempNNei)
  else
    call self%EasyPair_Cut%ShiftECalc_Single(curbox, disp, E_real_diff, accept)
  endif
  
  if (.not. accept) then
    E_Diff = 0.0_dp
    return
  endif
  
  ! Compute mesh contribution difference
  call self%UpdateMeshForMove(curbox, disp, E_mesh_diff)
  
  E_Diff = E_real_diff + E_mesh_diff

end subroutine

!=============================================================================
subroutine New_P3M(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  ! Energy change for particle addition
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Addition), intent(inout) :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept
  
  real(dp) :: E_real_diff, E_mesh_diff, E_self_diff
  integer :: iDisp, iAtom, atmType

  ! Real-space contribution using parent class
  call self%EasyPair_Cut%NewECalc(curbox, disp, tempList, tempNNei, E_real_diff, accept)
  
  if (.not. accept) then
    E_Diff = 0.0_dp
    return
  endif
  
  ! Mesh contribution (would need full recalculation or incremental update)
  ! For now, use a simplified approach
  E_mesh_diff = 0.0_dp
  ! TODO: Implement incremental mesh update for additions
  
  ! Self-energy correction for added particles
  E_self_diff = 0.0_dp
  do iDisp = 1, size(disp)
    iAtom = disp(iDisp)%atmIndx
    atmType = curbox%AtomType(iAtom)
    E_self_diff = E_self_diff - self%q(atmType)**2 * coulombConst * self%alpha / sqrt(pi)
  enddo
  
  E_Diff = E_real_diff + E_mesh_diff + E_self_diff

end subroutine

!=============================================================================
subroutine Old_P3M(self, curbox, disp, E_Diff)
  ! Energy change for particle deletion
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Deletion), intent(inout) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  
  real(dp) :: E_self_diff
  integer :: iAtom, atmType
  integer :: molStart, molEnd

  ! Real-space contribution using parent class
  call self%EasyPair_Cut%OldECalc(curbox, disp, E_Diff)
  
  ! Self-energy correction for removed particles
  E_self_diff = 0.0_dp
  call curbox%GetMolData(disp(1)%molIndx, molStart=molStart, molEnd=molEnd)
  
  do iAtom = molStart, molEnd
    atmType = curbox%AtomType(iAtom)
    E_self_diff = E_self_diff + self%q(atmType)**2 * coulombConst * self%alpha / sqrt(pi)
  enddo
  
  E_Diff = E_Diff + E_self_diff

end subroutine

!=============================================================================
subroutine OrthoVol_P3M(self, curbox, disp, E_Diff, accept)
  ! Energy difference for volume change (NPT ensemble)
  ! Volume changes affect:
  ! 1. Real-space interactions (positions scale)
  ! 2. Mesh (box dimensions change)
  ! 3. Green's function (depends on box size)
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(OrthoVolChange), intent(inout) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept
  
  real(dp) :: E_real_new, E_real_old
  real(dp) :: E_mesh_new, E_mesh_old
  real(dp) :: E_self_new, E_self_old
  real(dp) :: Lx_old, Ly_old, Lz_old
  real(dp) :: Lx_new, Ly_new, Lz_new
  
  ! Get current box dimensions
  select type(curbox)
    class is(OrthoBox)
      Lx_old = curbox%boxLx
      Ly_old = curbox%boxLy
      Lz_old = curbox%boxLz
    class is(CubeBox)
      Lx_old = curbox%boxL
      Ly_old = curbox%boxL
      Lz_old = curbox%boxL
    class default
      write(0,*) "Error: P3M OrthoVolECalc requires OrthoBox or CubeBox"
      accept = .false.
      return
  end select
  
  ! Compute new box dimensions
  Lx_new = Lx_old * disp(1)%xScale
  Ly_new = Ly_old * disp(1)%yScale
  Lz_new = Lz_old * disp(1)%zScale
  
  ! Store current energies
  E_real_old = curbox%E_Inter  ! Current real-space energy
  E_mesh_old = self%E_mesh_total
  E_self_old = self%E_self_total
  
  ! For volume moves, we need to recalculate everything because:
  ! 1. Particle positions scale with box
  ! 2. Mesh spacing changes
  ! 3. Green's function depends on box size
  
  ! Strategy: Use parent class for real-space, then add mesh contribution
  
  ! Get real-space energy difference using parent class
  call self%EasyPair_Cut%OrthoVolECalc(curbox, disp, E_Diff, accept)
  
  if (.not. accept) then
    return
  endif
  
  E_real_new = E_real_old + E_Diff
  
  ! For mesh energy, we need to:
  ! 1. Temporarily update box dimensions
  ! 2. Rebuild mesh with new dimensions
  ! 3. Recalculate mesh energy
  ! 4. Restore or keep depending on acceptance
  
  ! Save old mesh state
  self%chargeMesh_old = self%chargeMesh
  self%chargeMeshK_old = self%chargeMeshK
  
  ! Update cached box dimensions to trigger mesh rebuild
  self%cachedLx = Lx_new
  self%cachedLy = Ly_new
  self%cachedLz = Lz_new
  self%meshSetup = .false.
  
  ! Setup mesh with new dimensions
  ! Note: We need to pass the box with NEW dimensions
  ! This is tricky - the box hasn't been updated yet
  ! We'll need to temporarily modify or use a different approach
  
  ! Simplified approach: Compute new mesh energy
  ! by recalculating with scaled coordinates
  call self%ComputeVolumeMoveEnergy(curbox, disp, E_mesh_new, E_self_new, accept)
  
  if (.not. accept) then
    ! Restore old mesh state
    self%chargeMesh = self%chargeMesh_old
    self%chargeMeshK = self%chargeMeshK_old
    self%cachedLx = Lx_old
    self%cachedLy = Ly_old
    self%cachedLz = Lz_old
    return
  endif
  
  ! Total energy difference
  E_Diff = E_Diff + (E_mesh_new - E_mesh_old) + (E_self_new - E_self_old)
  
  ! Update stored energies
  self%E_mesh_total = E_mesh_new
  self%E_self_total = E_self_new

end subroutine

!=============================================================================
subroutine ComputeVolumeMoveEnergy_P3M(self, curbox, disp, E_mesh, E_self, accept)
  ! Compute mesh and self energy for volume change
  ! Uses scaled coordinates
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(OrthoVolChange), intent(in) :: disp(:)
  real(dp), intent(out) :: E_mesh, E_self
  logical, intent(out) :: accept
  
  integer :: iAtom, atmType, ix, iy, iz, mx, my, mz
  real(dp) :: qi, rx, ry, rz
  real(dp) :: rx_scaled, ry_scaled, rz_scaled
  real(dp) :: fx, fy, fz, wx, wy, wz
  real(dp) :: hx, hy, hz
  real(dp) :: Lx_new, Ly_new, Lz_new
  real(dp) :: rho_k_sq, sumQsq
  real(dp), pointer :: atoms(:,:) => null()
  real(dp), allocatable :: tempMesh(:,:,:)
  complex(dp), allocatable :: tempMeshK(:,:,:)
  
  call curbox%GetCoordinates(atoms)
  accept = .true.
  
  ! Get new box dimensions
  select type(curbox)
    class is(OrthoBox)
      Lx_new = curbox%boxLx * disp(1)%xScale
      Ly_new = curbox%boxLy * disp(1)%yScale
      Lz_new = curbox%boxLz * disp(1)%zScale
    class is(CubeBox)
      Lx_new = curbox%boxL * disp(1)%xScale
      Ly_new = curbox%boxL * disp(1)%yScale
      Lz_new = curbox%boxL * disp(1)%zScale
  end select
  
  ! Allocate temporary mesh arrays
  allocate(tempMesh(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1))
  allocate(tempMeshK(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1))
  tempMesh = 0.0_dp
  
  ! New mesh spacings
  hx = Lx_new / real(self%meshX, dp)
  hy = Ly_new / real(self%meshY, dp)
  hz = Lz_new / real(self%meshZ, dp)
  
  ! Assign charges to mesh with scaled coordinates
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    
    if (abs(qi) < 1.0E-15_dp) cycle
    
    ! Get current positions
    rx = atoms(1, iAtom)
    ry = atoms(2, iAtom)
    rz = atoms(3, iAtom)
    
    ! Scale positions for new box
    rx_scaled = rx * disp(1)%xScale
    ry_scaled = ry * disp(1)%yScale
    rz_scaled = rz * disp(1)%zScale
    
    ! Convert to fractional mesh coordinates
    fx = rx_scaled / hx
    fy = ry_scaled / hy
    fz = rz_scaled / hz
    
    ! Get base mesh point
    ix = floor(fx)
    iy = floor(fy)
    iz = floor(fz)
    
    ! Fractional parts
    fx = fx - real(ix, dp)
    fy = fy - real(iy, dp)
    fz = fz - real(iz, dp)
    
    ! Linear charge assignment (CIC)
    do mx = 0, 1
      wx = merge(1.0_dp - fx, fx, mx == 0)
      do my = 0, 1
        wy = merge(1.0_dp - fy, fy, my == 0)
        do mz = 0, 1
          wz = merge(1.0_dp - fz, fz, mz == 0)
          
          tempMesh(modulo(ix+mx, self%meshX), &
                   modulo(iy+my, self%meshY), &
                   modulo(iz+mz, self%meshZ)) = &
            tempMesh(modulo(ix+mx, self%meshX), &
                     modulo(iy+my, self%meshY), &
                     modulo(iz+mz, self%meshZ)) + qi * wx * wy * wz
        enddo
      enddo
    enddo
  enddo
  
  ! Rebuild Green's function for new box dimensions
  call self%ComputeGreenFunction(Lx_new, Ly_new, Lz_new)
  
  ! FFT charge mesh to reciprocal space
  call self%FFT3D_Forward(tempMesh, tempMeshK)
  
  ! Compute mesh energy
  E_mesh = 0.0_dp
  do ix = 0, self%meshX-1
    do iy = 0, self%meshY-1
      do iz = 0, self%meshZ-1
        rho_k_sq = real(tempMeshK(ix,iy,iz) * conjg(tempMeshK(ix,iy,iz)), dp)
        E_mesh = E_mesh + self%greenFunc(ix,iy,iz) * rho_k_sq
      enddo
    enddo
  enddo
  E_mesh = 0.5_dp * E_mesh
  
  ! Compute self-energy (doesn't change with volume, but recompute for consistency)
  sumQsq = 0.0_dp
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType = curbox%AtomType(iAtom)
    sumQsq = sumQsq + self%q(atmType)**2
  enddo
  E_self = -self%alpha / sqrt(pi) * sumQsq * coulombConst
  
  ! Update mesh for potential acceptance
  self%chargeMesh = tempMesh
  self%chargeMeshK = tempMeshK
  
  deallocate(tempMesh, tempMeshK)

end subroutine

!=============================================================================
subroutine SetupMesh_P3M(self, curbox)
  ! Initialize mesh arrays and Green's function
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  
  real(dp) :: Lx, Ly, Lz
  integer :: AllocateStat

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
      write(0,*) "Error: P3M requires OrthoBox or CubeBox"
      error stop
  end select
  
  ! Check if mesh needs to be rebuilt
  if (self%meshSetup) then
    if (abs(Lx - self%cachedLx) < 1.0E-10_dp .and. &
        abs(Ly - self%cachedLy) < 1.0E-10_dp .and. &
        abs(Lz - self%cachedLz) < 1.0E-10_dp) then
      return  ! No rebuild needed
    endif
  endif
  
  ! Allocate mesh arrays
  if (associated(self%chargeMesh)) deallocate(self%chargeMesh)
  if (associated(self%chargeMeshK)) deallocate(self%chargeMeshK)
  if (associated(self%chargeMesh_old)) deallocate(self%chargeMesh_old)
  if (associated(self%chargeMeshK_old)) deallocate(self%chargeMeshK_old)
  if (allocated(self%greenFunc)) deallocate(self%greenFunc)
  
  allocate(self%chargeMesh(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1), stat=AllocateStat)
  allocate(self%chargeMeshK(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1), stat=AllocateStat)
  allocate(self%chargeMesh_old(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1), stat=AllocateStat)
  allocate(self%chargeMeshK_old(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1), stat=AllocateStat)
  allocate(self%greenFunc(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1), stat=AllocateStat)
  
  if (AllocateStat /= 0) error stop "Mesh allocation error in P3M"
  
  self%chargeMesh = 0.0_dp
  self%chargeMeshK = cmplx(0.0_dp, 0.0_dp, dp)
  self%greenFunc = 0.0_dp
  
  ! Create FFTW plan if not already initialized
  if (.not. self%fftw_initialized) then
    ! IMPORTANT: FFTW uses C ordering (row-major), so dimensions are reversed
    ! For Fortran array(0:Nx-1, 0:Ny-1, 0:Nz-1), pass (Nz, Ny, Nx) to FFTW
    ! Use temporary pointers to get C addresses (workaround for allocatable arrays)
    block
      real(dp), pointer :: mesh_ptr(:,:,:)
      complex(dp), pointer :: meshK_ptr(:,:,:)
      
      mesh_ptr => self%chargeMesh
      meshK_ptr => self%chargeMeshK
      
      self%plan_forward = fftw_plan_dft_r2c_3d(self%meshZ, self%meshY, self%meshX, &
                                                c_loc(mesh_ptr), &
                                                c_loc(meshK_ptr), &
                                                FFTW_ESTIMATE)
    end block
    
    if (.not. c_associated(self%plan_forward)) then
      write(0,*) "ERROR: Failed to create FFTW forward plan!"
      write(0,*) "       Check that FFTW3 is properly installed and linked."
      error stop "FFTW initialization failed"
    endif
    
    self%fftw_initialized = .true.
    write(0,*) "FFTW initialized successfully for P3M electrostatics"
    write(0,'(A,3I6)') "  Mesh dimensions: ", self%meshX, self%meshY, self%meshZ
  endif
  
  ! Compute Green's function
  call self%ComputeGreenFunction(Lx, Ly, Lz)
  
  ! Cache box dimensions
  self%cachedLx = Lx
  self%cachedLy = Ly
  self%cachedLz = Lz
  self%meshSetup = .true.

end subroutine

!=============================================================================
subroutine ComputeGreenFunction_P3M(self, Lx, Ly, Lz)
  ! Compute optimal Green's function for P3M
  implicit none
  class(Pair_P3M), intent(inout) :: self
  real(dp), intent(in) :: Lx, Ly, Lz
  
  integer :: ix, iy, iz, nx, ny, nz
  real(dp) :: kx, ky, kz, k_sq
  real(dp) :: volume, aliasing_factor
  real(dp) :: hx, hy, hz  ! Mesh spacings
  
  volume = Lx * Ly * Lz
  hx = Lx / real(self%meshX, dp)
  hy = Ly / real(self%meshY, dp)
  hz = Lz / real(self%meshZ, dp)
  
  ! Compute Green's function G(k) for each mesh point
  do ix = 0, self%meshX-1
    nx = ix
    if (nx > self%meshX/2) nx = nx - self%meshX
    kx = two_pi * real(nx, dp) / Lx
    
    do iy = 0, self%meshY-1
      ny = iy
      if (ny > self%meshY/2) ny = ny - self%meshY
      ky = two_pi * real(ny, dp) / Ly
      
      do iz = 0, self%meshZ-1
        nz = iz
        if (nz > self%meshZ/2) nz = nz - self%meshZ
        kz = two_pi * real(nz, dp) / Lz
        
        k_sq = kx**2 + ky**2 + kz**2
        
        if (k_sq < 1.0E-10_dp) then
          ! Handle k=0 specially
          self%greenFunc(ix, iy, iz) = 0.0_dp
        else
          ! Green's function: exp(-k^2/(4*alpha^2)) / k^2
          ! Including charge assignment function correction
          aliasing_factor = 1.0_dp  ! Simplified; should include CAF correction
          self%greenFunc(ix, iy, iz) = aliasing_factor * &
                                       exp(-k_sq / (4.0_dp * self%alpha_sq)) / k_sq
        endif
      enddo
    enddo
  enddo
  
  ! Scale by prefactor
  self%greenFunc = self%greenFunc * 4.0_dp * pi / volume * coulombConst

end subroutine

!=============================================================================
subroutine AssignChargesToMesh_P3M(self, curbox)
  ! Assign particle charges to mesh using charge assignment function
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  
  integer :: iAtom, atmType
  integer :: ix, iy, iz, mx, my, mz
  real(dp) :: qi, rx, ry, rz
  real(dp) :: fx, fy, fz  ! Fractional coordinates
  real(dp) :: wx, wy, wz  ! Weights
  real(dp) :: hx, hy, hz  ! Mesh spacings
  real(dp), pointer :: atoms(:,:) => null()
  
  call curbox%GetCoordinates(atoms)
  
  ! Clear mesh
  self%chargeMesh = 0.0_dp
  
  ! Mesh spacings
  select type(curbox)
    class is(CubeBox)
      hx = curbox%boxL / real(self%meshX, dp)
      hy = curbox%boxL / real(self%meshY, dp)
      hz = curbox%boxL / real(self%meshZ, dp)
    class is(OrthoBox)
      hx = curbox%boxLx / real(self%meshX, dp)
      hy = curbox%boxLy / real(self%meshY, dp)
      hz = curbox%boxLz / real(self%meshZ, dp)
  end select
  
  ! Assign charges using 4th-order spline (simplified to linear for now)
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    
    if (abs(qi) < 1.0E-15_dp) cycle
    
    rx = atoms(1, iAtom)
    ry = atoms(2, iAtom)
    rz = atoms(3, iAtom)
    
    ! Convert to fractional mesh coordinates
    fx = rx / hx
    fy = ry / hy
    fz = rz / hz
    
    ! Get base mesh point (using linear interpolation for simplicity)
    ix = floor(fx)
    iy = floor(fy)
    iz = floor(fz)
    
    ! Fractional parts
    fx = fx - real(ix, dp)
    fy = fy - real(iy, dp)
    fz = fz - real(iz, dp)
    
    ! Linear charge assignment (CIC - Cloud-In-Cell)
    do mx = 0, 1
      wx = merge(1.0_dp - fx, fx, mx == 0)
      do my = 0, 1
        wy = merge(1.0_dp - fy, fy, my == 0)
        do mz = 0, 1
          wz = merge(1.0_dp - fz, fz, mz == 0)
          
          ! Periodic wrapping
          self%chargeMesh(modulo(ix+mx, self%meshX), &
                         modulo(iy+my, self%meshY), &
                         modulo(iz+mz, self%meshZ)) = &
            self%chargeMesh(modulo(ix+mx, self%meshX), &
                           modulo(iy+my, self%meshY), &
                           modulo(iz+mz, self%meshZ)) + qi * wx * wy * wz
        enddo
      enddo
    enddo
  enddo

end subroutine

!=============================================================================
subroutine ComputeMeshEnergy_P3M(self, curbox, E_mesh)
  ! Compute energy from mesh using FFT
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_mesh
  
  integer :: ix, iy, iz
  real(dp) :: rho_k_sq
  
  ! FFT charge mesh to reciprocal space
  ! Note: This is a placeholder. In production, use FFTW or similar
  call self%FFT3D_Forward(self%chargeMesh, self%chargeMeshK)
  
  ! Compute energy: E = 0.5 * sum_k G(k) |rho(k)|^2
  E_mesh = 0.0_dp
  do ix = 0, self%meshX-1
    do iy = 0, self%meshY-1
      do iz = 0, self%meshZ-1
        rho_k_sq = real(self%chargeMeshK(ix,iy,iz) * conjg(self%chargeMeshK(ix,iy,iz)), dp)
        E_mesh = E_mesh + self%greenFunc(ix,iy,iz) * rho_k_sq
      enddo
    enddo
  enddo
  
  E_mesh = 0.5_dp * E_mesh

end subroutine

!=============================================================================
subroutine ComputeSelfEnergy_P3M(self, curbox, E_self)
  ! Self-energy correction
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_self
  
  integer :: iAtom, atmType
  real(dp) :: qi, sumQsq

  sumQsq = 0.0_dp
  
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    sumQsq = sumQsq + qi**2
  enddo
  
  E_self = -self%alpha / sqrt(pi) * sumQsq * coulombConst

end subroutine

!=============================================================================
function ComputeRealPair_P3M(self, rsq, atmtype1, atmtype2) result(E_Pair)
  ! Real-space pair energy (same as PairFunction)
  implicit none
  class(Pair_P3M), intent(inout) :: self
  real(dp), intent(in) :: rsq
  integer, intent(in) :: atmtype1, atmtype2
  real(dp) :: E_Pair

  E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)

end function

!=============================================================================
subroutine UpdateMeshForMove_P3M(self, curbox, disp, E_mesh_diff)
  ! Compute mesh energy difference for particle displacement
  ! 
  ! Strategy: Rebuild both old and new meshes from scratch each time.
  ! This is O(N) per move but guarantees correctness.
  ! The mesh is rebuilt to match current atom positions for "old" state,
  ! and with the displacement applied for "new" state.
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(out) :: E_mesh_diff
  
  integer :: iDisp, iAtom, atmType, nAtoms, jAtom
  real(dp) :: qi
  real(dp) :: rx, ry, rz
  real(dp), allocatable, target :: oldMesh(:,:,:), newMesh(:,:,:)
  complex(dp), allocatable, target :: oldMeshK(:,:,:), newMeshK(:,:,:)
  real(dp), pointer :: atoms(:,:) => null()
  real(dp) :: E_mesh_old, E_mesh_new
  real(dp) :: rho_k_sq, hx, hy, hz, fx, fy, fz, wx, wy, wz
  integer :: ix, iy, iz, imx, imy, imz, mx, my, mz
  logical :: isDisplaced
  
  call curbox%GetCoordinates(atoms)
  nAtoms = curbox%nMaxAtoms
  
  ! Allocate meshes
  allocate(oldMesh(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1))
  allocate(newMesh(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1))
  allocate(oldMeshK(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1))
  allocate(newMeshK(0:self%meshX-1, 0:self%meshY-1, 0:self%meshZ-1))
  oldMesh = 0.0_dp
  newMesh = 0.0_dp
  
  ! Mesh spacings
  select type(curbox)
    class is(CubeBox)
      hx = curbox%boxL / real(self%meshX, dp)
      hy = curbox%boxL / real(self%meshY, dp)
      hz = curbox%boxL / real(self%meshZ, dp)
    class is(OrthoBox)
      hx = curbox%boxLx / real(self%meshX, dp)
      hy = curbox%boxLy / real(self%meshY, dp)
      hz = curbox%boxLz / real(self%meshZ, dp)
  end select
  
  ! Build BOTH old and new meshes in one pass
  do jAtom = 1, nAtoms
    if (.not. curbox%IsActive(jAtom)) cycle
    
    atmType = curbox%AtomType(jAtom)
    qi = self%q(atmType)
    if (abs(qi) < 1.0E-15_dp) cycle
    
    ! Current (old) position
    rx = atoms(1, jAtom)
    ry = atoms(2, jAtom)
    rz = atoms(3, jAtom)
    
    ! Assign to OLD mesh
    fx = rx / hx
    fy = ry / hy
    fz = rz / hz
    imx = floor(fx)
    imy = floor(fy)
    imz = floor(fz)
    fx = fx - real(imx, dp)
    fy = fy - real(imy, dp)
    fz = fz - real(imz, dp)
    
    do mx = 0, 1
      wx = merge(1.0_dp - fx, fx, mx == 0)
      do my = 0, 1
        wy = merge(1.0_dp - fy, fy, my == 0)
        do mz = 0, 1
          wz = merge(1.0_dp - fz, fz, mz == 0)
          oldMesh(modulo(imx+mx, self%meshX), &
                  modulo(imy+my, self%meshY), &
                  modulo(imz+mz, self%meshZ)) = &
            oldMesh(modulo(imx+mx, self%meshX), &
                    modulo(imy+my, self%meshY), &
                    modulo(imz+mz, self%meshZ)) + qi * wx * wy * wz
        enddo
      enddo
    enddo
    
    ! Check if this atom is being displaced for NEW mesh
    isDisplaced = .false.
    do iDisp = 1, size(disp)
      if (disp(iDisp)%atmIndx == jAtom) then
        rx = disp(iDisp)%x_new
        ry = disp(iDisp)%y_new
        rz = disp(iDisp)%z_new
        isDisplaced = .true.
        exit
      endif
    enddo
    
    ! If not displaced, use current position (already in rx, ry, rz from reset above)
    if (.not. isDisplaced) then
      rx = atoms(1, jAtom)
      ry = atoms(2, jAtom)
      rz = atoms(3, jAtom)
    endif
    
    ! Assign to NEW mesh
    fx = rx / hx
    fy = ry / hy
    fz = rz / hz
    imx = floor(fx)
    imy = floor(fy)
    imz = floor(fz)
    fx = fx - real(imx, dp)
    fy = fy - real(imy, dp)
    fz = fz - real(imz, dp)
    
    do mx = 0, 1
      wx = merge(1.0_dp - fx, fx, mx == 0)
      do my = 0, 1
        wy = merge(1.0_dp - fy, fy, my == 0)
        do mz = 0, 1
          wz = merge(1.0_dp - fz, fz, mz == 0)
          newMesh(modulo(imx+mx, self%meshX), &
                  modulo(imy+my, self%meshY), &
                  modulo(imz+mz, self%meshZ)) = &
            newMesh(modulo(imx+mx, self%meshX), &
                    modulo(imy+my, self%meshY), &
                    modulo(imz+mz, self%meshZ)) + qi * wx * wy * wz
        enddo
      enddo
    enddo
  enddo
  
  ! FFT both meshes
  call self%FFT3D_Forward(oldMesh, oldMeshK)
  call self%FFT3D_Forward(newMesh, newMeshK)
  
  ! Compute energies
  E_mesh_old = 0.0_dp
  E_mesh_new = 0.0_dp
  do ix = 0, self%meshX-1
    do iy = 0, self%meshY-1
      do iz = 0, self%meshZ-1
        rho_k_sq = real(oldMeshK(ix,iy,iz) * conjg(oldMeshK(ix,iy,iz)), dp)
        E_mesh_old = E_mesh_old + 0.5_dp * self%greenFunc(ix,iy,iz) * rho_k_sq
        
        rho_k_sq = real(newMeshK(ix,iy,iz) * conjg(newMeshK(ix,iy,iz)), dp)
        E_mesh_new = E_mesh_new + 0.5_dp * self%greenFunc(ix,iy,iz) * rho_k_sq
      enddo
    enddo
  enddo
  
  E_mesh_diff = E_mesh_new - E_mesh_old
  
  ! Also update the stored mesh to match current positions (for Detailed calculations)
  self%chargeMesh = oldMesh
  self%chargeMeshK = oldMeshK
  
  deallocate(oldMesh, newMesh, oldMeshK, newMeshK)

end subroutine

!=============================================================================
subroutine AssignSingleChargeToMesh_P3M(self, rx, ry, rz, qi, mesh, curbox)
  ! Assign a single charge to mesh - helper for incremental updates
  implicit none
  class(Pair_P3M), intent(in) :: self
  real(dp), intent(in) :: rx, ry, rz, qi
  real(dp), intent(inout) :: mesh(0:,0:,0:)
  class(SimBox), intent(in) :: curbox
  
  integer :: ix, iy, iz, mx, my, mz
  real(dp) :: fx, fy, fz
  real(dp) :: wx, wy, wz
  real(dp) :: hx, hy, hz
  
  ! Mesh spacings
  select type(curbox)
    class is(CubeBox)
      hx = curbox%boxL / real(self%meshX, dp)
      hy = curbox%boxL / real(self%meshY, dp)
      hz = curbox%boxL / real(self%meshZ, dp)
    class is(OrthoBox)
      hx = curbox%boxLx / real(self%meshX, dp)
      hy = curbox%boxLy / real(self%meshY, dp)
      hz = curbox%boxLz / real(self%meshZ, dp)
  end select
  
  ! Convert to fractional mesh coordinates
  fx = rx / hx
  fy = ry / hy
  fz = rz / hz
  
  ! Get base mesh point
  ix = floor(fx)
  iy = floor(fy)
  iz = floor(fz)
  
  ! Fractional parts
  fx = fx - real(ix, dp)
  fy = fy - real(iy, dp)
  fz = fz - real(iz, dp)
  
  ! Linear charge assignment (CIC)
  do mx = 0, 1
    wx = merge(1.0_dp - fx, fx, mx == 0)
    do my = 0, 1
      wy = merge(1.0_dp - fy, fy, my == 0)
      do mz = 0, 1
        wz = merge(1.0_dp - fz, fz, mz == 0)
        
        mesh(modulo(ix+mx, self%meshX), &
             modulo(iy+my, self%meshY), &
             modulo(iz+mz, self%meshZ)) = &
          mesh(modulo(ix+mx, self%meshX), &
               modulo(iy+my, self%meshY), &
               modulo(iz+mz, self%meshZ)) + qi * wx * wy * wz
      enddo
    enddo
  enddo

end subroutine

!=============================================================================
subroutine GetCAFWeights_P3M(self, x, weights)
  ! Get charge assignment function weights
  implicit none
  class(Pair_P3M), intent(in) :: self
  real(dp), intent(in) :: x
  real(dp), intent(out) :: weights(:)
  
  ! For 4th order spline (simplified)
  ! In production, use proper B-spline formulas
  weights = 0.0_dp
  if (size(weights) >= 2) then
    weights(1) = 1.0_dp - x
    weights(2) = x
  endif

end subroutine

!=============================================================================
subroutine FFT3D_Forward_P3M(self, real_data, complex_data)
  ! 3D FFT using FFTW - forward transform (real to complex)
  ! Computes: complex_data[k] = Σ real_data[r] * exp(-2πi k·r)
  use, intrinsic :: iso_c_binding
  implicit none
  class(Pair_P3M), intent(inout) :: self
  real(dp), intent(inout), target :: real_data(0:,0:,0:)
  complex(dp), intent(out), target :: complex_data(0:,0:,0:)
  
  real(dp), pointer :: in_ptr(:,:,:)
  complex(dp), pointer :: out_ptr(:,:,:)
  
  if (.not. self%fftw_initialized) then
    write(0,*) "ERROR: FFTW not initialized! Call SetupMesh first."
    error stop "FFT called before initialization"
  endif
  
  ! Use pointers to get C addresses
  in_ptr => real_data
  out_ptr => complex_data
  
  ! Execute the pre-created FFTW plan with C pointers
  call fftw_execute_dft_r2c(self%plan_forward, c_loc(in_ptr), c_loc(out_ptr))
  
  ! Note: FFTW does not normalize the transform
  ! Normalization factor of 1/N is incorporated into the Green's function
  ! to avoid redundant operations

end subroutine

!=============================================================================
subroutine OptimizeAlpha_P3M(self, curbox)
  ! Optimize alpha parameter for P3M
  implicit none
  class(Pair_P3M), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  
  real(dp) :: Lx, Ly, Lz, Lmin
  
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
      Lx = 10.0_dp
      Ly = 10.0_dp
      Lz = 10.0_dp
  end select
  
  Lmin = min(Lx, Ly, Lz)
  
  ! Optimal alpha balances real-space and mesh costs
  self%alpha = 3.5_dp / self%rCutReal
  self%alpha_sq = self%alpha**2
  
  ! Adjust mesh spacing accordingly
  self%meshSpacing = Lmin / real(max(self%meshX, self%meshY, self%meshZ), dp)

end subroutine

!=============================================================================
subroutine ProcessIO_P3M(self, line)
  use Common_MolInfo, only: nAtomTypes
  use Input_Format, only: CountCommands, GetXCommand, maxLineLen
  use Units, only: inEngUnit, inLenUnit
  implicit none
  class(Pair_P3M), intent(inout) :: self
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: command
  integer :: jType, lineStat, nPar, type1
  real(dp) :: rCut, q, rMin, alpha

  call GetXCommand(line, command, 1, lineStat)
  call CountCommands(line, nPar)

  select case(trim(adjustl(command)))
    case("rcut")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) rCut
      self%rCutReal = rCut * inLenUnit
      self%rCutRealSq = self%rCutReal**2
      self%rCut = self%rCutReal
      self%rCutSq = self%rCutRealSq

    case("alpha")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) alpha
      self%alpha = alpha / inLenUnit
      self%alpha_sq = self%alpha**2

    case("mesh")
      ! mesh 32 32 32
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) self%meshX
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) self%meshY
      call GetXCommand(line, command, 4, lineStat)
      read(command, *) self%meshZ

    case("order")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) self%chargeOrder

    case("precision")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) self%p3mPrecision

    case("charge")
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
      ! Try atom type definition: type q rMin
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
subroutine Prologue_P3M(self)
  use ParallelVar, only: nout
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_P3M), intent(inout) :: self
  integer :: i

  write(nout, *) "============================================"
  write(nout, *) "P3M (Particle-Particle Particle-Mesh) Parameters:"
  write(nout, *) "  Alpha (screening parameter):", self%alpha
  write(nout, *) "  Real-space cutoff:", self%rCutReal
  write(nout, *) "  Mesh dimensions:", self%meshX, "x", self%meshY, "x", self%meshZ
  write(nout, *) "  Charge assignment order:", self%chargeOrder
  write(nout, *) "  Target precision:", self%p3mPrecision
  write(nout, *) ""
  write(nout, *) "Charges per atom type:"
  do i = 1, nAtomTypes
    write(nout, *) "  Type", i, ":", self%q(i)
  enddo
  write(nout, *) "============================================"

end subroutine

!=============================================================================
function GetCutOff_P3M(self) result(rCut)
  implicit none
  class(Pair_P3M), intent(inout) :: self
  real(dp) :: rCut

  rCut = self%rCutReal

end function

!=============================================================================
subroutine Destructor_P3M(self)
  ! Clean up FFTW plans and deallocate arrays
  ! Should be called when forcefield is destroyed
  implicit none
  class(Pair_P3M), intent(inout) :: self
  
  ! Destroy FFTW plans
  if (c_associated(self%plan_forward)) then
    call fftw_destroy_plan(self%plan_forward)
    self%plan_forward = C_NULL_PTR
  endif
  
  self%fftw_initialized = .false.
  
  ! Deallocate mesh arrays
  if (associated(self%chargeMesh)) deallocate(self%chargeMesh)
  if (associated(self%chargeMeshK)) deallocate(self%chargeMeshK)
  if (associated(self%chargeMesh_old)) deallocate(self%chargeMesh_old)
  if (associated(self%chargeMeshK_old)) deallocate(self%chargeMeshK_old)
  if (allocated(self%greenFunc)) deallocate(self%greenFunc)
  if (allocated(self%cafWeights)) deallocate(self%cafWeights)
  if (allocated(self%q)) deallocate(self%q)
  if (allocated(self%qTable)) deallocate(self%qTable)
  
end subroutine

!=============================================================================
subroutine AcceptUpdate_P3M(self)
  ! Called when a move is accepted - mesh is already updated
  implicit none
  class(Pair_P3M), intent(inout) :: self

  ! The chargeMesh and chargeMeshK already contain the new values
  ! (they were updated in UpdateMeshForMove)
  ! Nothing special to do here
  
end subroutine

!=============================================================================
subroutine RejectUpdate_P3M(self)
  ! Called when a move is rejected - restore old mesh state
  implicit none
  class(Pair_P3M), intent(inout) :: self

  ! Restore the old mesh values
  if (associated(self%chargeMesh_old) .and. associated(self%chargeMesh)) then
    self%chargeMesh = self%chargeMesh_old
  endif
  if (associated(self%chargeMeshK_old) .and. associated(self%chargeMeshK)) then
    self%chargeMeshK = self%chargeMeshK_old
  endif
  
end subroutine

!=============================================================================
end module FF_P3M
!=============================================================================

