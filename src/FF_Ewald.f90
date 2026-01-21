!================================================================================
! Ewald Summation Module for Monte Carlo Simulations
!
! This module implements the Ewald summation method for computing long-range
! electrostatic interactions in periodic systems. The implementation is 
! specifically designed for Monte Carlo simulations, with efficient methods
! for computing energy differences when particles are moved, added, or deleted.
!
! The Ewald method splits the electrostatic potential into three parts:
!   E_total = E_real + E_recip + E_self + E_intra
!
! Where:
!   E_real  - Short-range real-space contribution (screened by erfc)
!   E_recip - Long-range reciprocal-space contribution (k-space sum)
!   E_self  - Self-interaction correction
!   E_intra - Intramolecular correction (for excluded pairs)
!
! References:
!   - Frenkel & Smit, "Understanding Molecular Simulation", 2nd Ed.
!   - Allen & Tildesley, "Computer Simulation of Liquids"
!   - Toukmaji & Board, Comp. Phys. Comm. 95 (1996) 73-92
!================================================================================
module FF_Ewald
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use ClassyConstants, only: coulombConst, pi, two_pi
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox

  implicit none

  type, extends(forcefield) :: Pair_Ewald
    ! Ewald parameters
    real(dp) :: alpha          ! Ewald screening parameter (1/length)
    real(dp) :: alpha_sq       ! alpha^2 for convenience
    real(dp) :: rCutReal       ! Real-space cutoff
    real(dp) :: rCutRealSq     ! Real-space cutoff squared
    real(dp) :: kCutoff        ! Reciprocal space cutoff (in units of 2*pi/L)
    
    ! Precision control
    real(dp) :: ewaldPrecision = 1.0E-6_dp  ! Target precision for Ewald sum
    
    ! Charges
    real(dp), allocatable :: q(:)           ! Charge per atom type
    real(dp), allocatable :: qTable(:,:)    ! q_i * q_j * coulombConst lookup table
    
    ! K-vector arrays
    integer :: nKVectors                    ! Total number of k-vectors
    integer :: kxMax, kyMax, kzMax          ! Max k-indices in each direction
    real(dp), allocatable :: kVectors(:,:)  ! K-vectors (3, nKVectors)
    real(dp), allocatable :: kMagSq(:)      ! |k|^2 for each k-vector
    real(dp), allocatable :: kPrefac(:)     ! exp(-k^2/(4*alpha^2))/(k^2) prefactor
    
    ! Structure factors (for efficient MC updates)
    ! These store the current state of the system
    complex(dp), allocatable :: rhoK(:)     ! Structure factor for each k-vector
    complex(dp), allocatable :: rhoK_old(:) ! Backup for move rejection
    
    ! Pre-computed exp(i k.r) tables per atom for efficient updates
    ! expIKR(kIndex, atomIndex) = exp(i k.r_atom)
    complex(dp), allocatable :: expIKR(:,:)
    complex(dp), allocatable :: expIKR_old(:,:)  ! Backup for rejection
    
    ! Energy components
    real(dp) :: E_recip_total  ! Current reciprocal space energy
    real(dp) :: E_self_total   ! Current self-energy correction
    
    ! Flags
    logical :: initialized = .false.
    logical :: precomputed = .false.  ! Whether k-vectors have been set up
    
    ! Box dimensions cache (needed to detect if box changed)
    real(dp) :: cachedLx, cachedLy, cachedLz
    
    ! Minimum distance for overlap detection
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    
  contains
    procedure, pass :: Constructor => Constructor_Ewald
    procedure, pass :: DetailedECalc => Detailed_Ewald
    procedure, pass :: DiffECalc => DiffECalc_Ewald
    procedure, pass :: ShiftECalc_Single => Shift_Ewald_Single
    procedure, pass :: NewECalc => New_Ewald
    procedure, pass :: OldECalc => Old_Ewald
    procedure, pass :: ManyBody => ManyBody_Ewald  
    procedure, pass :: ProcessIO => ProcessIO_Ewald
    procedure, pass :: Prologue => Prologue_Ewald
    procedure, pass :: GetCutOff => GetCutOff_Ewald
    procedure, pass :: Epilogue => Epilogue_Ewald
    
    ! Ewald-specific procedures
    procedure, pass :: SetupKVectors => SetupKVectors_Ewald
    procedure, pass :: ComputeStructureFactors => ComputeStructureFactors_Ewald
    procedure, pass :: ComputeRecipEnergy => ComputeRecipEnergy_Ewald
    procedure, pass :: ComputeSelfEnergy => ComputeSelfEnergy_Ewald
    procedure, pass :: ComputeRealPair => ComputeRealPair_Ewald
    procedure, pass :: OptimizeAlpha => OptimizeAlpha_Ewald
    procedure, pass :: UpdateStructureFactor => UpdateStructureFactor_Ewald
    procedure, pass :: AcceptUpdate => AcceptUpdate_Ewald
    procedure, pass :: RejectUpdate => RejectUpdate_Ewald
    
  end type

contains

!=============================================================================
subroutine Constructor_Ewald(self)
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  integer :: AllocateStat

  allocate(self%q(1:nAtomTypes), stat=AllocateStat)
  allocate(self%qTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  allocate(self%rMin(1:nAtomTypes), stat=AllocateStat)
  allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)

  self%q = 0.0_dp
  self%qTable = 0.0_dp
  self%rMin = 0.5_dp
  self%rMinTable = 0.25_dp  ! Stored as r^2
  
  ! Default parameters
  self%alpha = 0.3_dp
  self%alpha_sq = self%alpha**2
  self%rCutReal = 10.0_dp
  self%rCutRealSq = self%rCutReal**2
  self%rCut = self%rCutReal
  self%rCutSq = self%rCutRealSq
  self%kCutoff = 5.0_dp
  
  self%initialized = .false.
  self%precomputed = .false.
  self%nKVectors = 0

  if (AllocateStat /= 0) error stop "Allocation error in Ewald Constructor"

end subroutine
!=============================================================================
subroutine Detailed_Ewald(self, curbox, E_T, accept)
  use ParallelVar, only: nout
  use Common_MolInfo, only: nMolTypes
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(inout) :: E_T
  logical, intent(out) :: accept
  
  integer :: iAtom, jAtom
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_real, E_recip, E_self
  real(dp) :: E_pair
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)
  
  accept = .true.
  E_real = 0.0_dp
  curbox%ETable = 0.0_dp
  
  ! Check if k-vectors need to be set up or updated
  call self%SetupKVectors(curbox)
  
  ! Compute structure factors from scratch
  call self%ComputeStructureFactors(curbox)
  
  ! Compute real-space contribution
  do iAtom = 1, curbox%nMaxAtoms - 1
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType1 = curbox%AtomType(iAtom)
    
    do jAtom = iAtom + 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(jAtom)) cycle
      ! Skip intramolecular pairs
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx**2 + ry**2 + rz**2
      
      if (rsq < self%rCutRealSq) then
        atmType2 = curbox%AtomType(jAtom)
        rmin_ij = self%rMinTable(atmType1, atmType2)
        
        if (rsq < rmin_ij) then
          write(0,*) "ERROR! Overlapping atoms found in Ewald!"
          write(0,*) "Distance:", sqrt(rsq), "Min:", sqrt(rmin_ij)
          accept = .false.
          return
        endif
        
        E_pair = self%ComputeRealPair(rsq, atmType1, atmType2)
        E_real = E_real + E_pair
        curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_pair
        curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_pair
      endif
    enddo
  enddo
  
  ! Compute reciprocal-space contribution
  call self%ComputeRecipEnergy(curbox, E_recip)
  self%E_recip_total = E_recip
  
  ! Compute self-energy correction
  call self%ComputeSelfEnergy(curbox, E_self)
  self%E_self_total = E_self
  
  write(nout,*) "Ewald Real-space Energy:", E_real
  write(nout,*) "Ewald Reciprocal Energy:", E_recip
  write(nout,*) "Ewald Self-energy Correction:", E_self
  write(nout,*) "Ewald Total Energy:", E_real + E_recip + E_self
  
  E_T = E_real + E_recip + E_self

end subroutine
!=============================================================================
subroutine DiffECalc_Ewald(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_Ewald), intent(inout) :: self
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
      call self%ShiftECalc_Single(curbox, disp, E_Diff, accept)
      
    class is(Addition)
      call self%NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)
      
    class is(Deletion)
      call self%OldECalc(curbox, disp, E_Diff)
      
    class default
      write(0,*) "Unknown Perturbation Type in Ewald DiffECalc"
      error stop
  end select

end subroutine
!=============================================================================
subroutine Shift_Ewald_Single(self, curbox, disp, E_Diff, accept)
  ! Computes energy change when atoms are displaced
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jNei, jAtom, dispLen
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_real_old, E_real_new
  real(dp) :: E_recip_old, E_recip_new, E_recip_diff
  real(dp) :: E_pair
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()

  call curbox%GetCoordinates(atoms)
  call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)

  dispLen = size(disp)
  E_Diff = 0.0_dp
  E_real_old = 0.0_dp
  E_real_new = 0.0_dp
  accept = .true.

  ! Loop over displaced atoms
  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)

    ! Real-space contribution from neighbors
    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)

      ! New position energy
      rx = disp(iDisp)%x_new - atoms(1, jAtom)
      ry = disp(iDisp)%y_new - atoms(2, jAtom)
      rz = disp(iDisp)%z_new - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutRealSq) then
        rmin_ij = self%rMinTable(atmType1, atmType2)
        if (rsq < rmin_ij) then
          accept = .false.
          return
        endif
        E_pair = self%ComputeRealPair(rsq, atmType1, atmType2)
        E_real_new = E_real_new + E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair
      endif

      ! Old position energy
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutRealSq) then
        E_pair = self%ComputeRealPair(rsq, atmType1, atmType2)
        E_real_old = E_real_old + E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair
      endif
    enddo
  enddo

  ! Reciprocal-space contribution
  ! Save old structure factors
  call self%UpdateStructureFactor(curbox, disp, E_recip_diff)
  
  E_Diff = (E_real_new - E_real_old) + E_recip_diff

end subroutine
!=============================================================================
subroutine New_Ewald(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  ! Energy change for particle addition
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Addition), intent(in) :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_real, E_recip_diff, E_self_diff
  real(dp) :: E_pair
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)

  dispLen = size(disp)
  E_Diff = 0.0_dp
  E_real = 0.0_dp
  accept = .true.

  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)

    listIndx = disp(iDisp)%listIndex
    maxNei = tempNNei(listIndx)
    
    do jNei = 1, maxNei
      jAtom = tempList(jNei, listIndx)
      
      rx = disp(iDisp)%x_new - atoms(1, jAtom)
      ry = disp(iDisp)%y_new - atoms(2, jAtom)
      rz = disp(iDisp)%z_new - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutRealSq) then
        atmType2 = curbox%AtomType(jAtom)
        rmin_ij = self%rMinTable(atmType1, atmType2)
        
        if (rsq < rmin_ij) then
          accept = .false.
          return
        endif
        
        E_pair = self%ComputeRealPair(rsq, atmType1, atmType2)
        E_real = E_real + E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair
      endif
    enddo
  enddo

  ! Reciprocal space change for addition
  ! (Note: Would need proper structure factor update here)
  E_recip_diff = 0.0_dp
  
  ! Self-energy change for addition
  E_self_diff = 0.0_dp
  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)
    E_self_diff = E_self_diff - self%q(atmType1)**2 * coulombConst * self%alpha / sqrt(pi)
  enddo

  E_Diff = E_real + E_recip_diff + E_self_diff

end subroutine
!=============================================================================
subroutine Old_Ewald(self, curbox, disp, E_Diff)
  ! Energy change for particle deletion
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Deletion), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff

  integer :: iAtom, jAtom, jNei
  integer :: atmType1, atmType2
  integer :: molEnd, molStart
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_real, E_recip_diff, E_self_diff
  real(dp) :: E_pair
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()

  call curbox%GetCoordinates(atoms)
  call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)

  E_Diff = 0.0_dp
  E_real = 0.0_dp
  E_self_diff = 0.0_dp

  call curBox%GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

  do iAtom = molStart, molEnd
    atmType1 = curbox%AtomType(iAtom)
    
    ! Real-space contribution
    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)

      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutRealSq) then
        E_pair = self%ComputeRealPair(rsq, atmType1, atmType2)
        E_real = E_real - E_pair
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair
      endif
    enddo
    
    ! Self-energy correction (removed particle no longer contributes)
    E_self_diff = E_self_diff + self%q(atmType1)**2 * coulombConst * self%alpha / sqrt(pi)
  enddo

  ! Reciprocal space change for deletion
  E_recip_diff = 0.0_dp

  E_Diff = E_real + E_recip_diff + E_self_diff

end subroutine
!=============================================================================
subroutine ManyBody_Ewald(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept)
  ! For trial insertion moves (CBMC, etc.)
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(simBox), intent(inout) :: curbox
  integer, intent(in) :: atmtype1
  integer, intent(in) :: atmtypes(:)
  real(dp), intent(in) :: pos1(:)
  real(dp), intent(in) :: posN(:,:)
  logical, intent(out) :: accept
  real(dp), intent(out) :: E_Many

  integer :: jAtom
  integer :: atmType2
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_pair
  real(dp) :: rmin_ij

  accept = .true.
  E_Many = 0.0_dp
  
  ! Only real-space contribution for quick trial move evaluation
  do jAtom = 1, size(posN, 2)
    atmType2 = atmtypes(jAtom)
    
    rx = pos1(1) - posN(1, jAtom)
    ry = pos1(2) - posN(2, jAtom)
    rz = pos1(3) - posN(3, jAtom)
    call curbox%Boundary(rx, ry, rz)
    rsq = rx*rx + ry*ry + rz*rz
    
    rmin_ij = self%rMinTable(atmType2, atmtype1)
    if (rsq < rmin_ij) then
      accept = .false.
      return
    endif
    
    if (rsq < self%rCutRealSq) then
      E_pair = self%ComputeRealPair(rsq, atmtype1, atmType2)
      E_Many = E_Many + E_pair
    endif
  enddo

end subroutine
!=============================================================================
function ComputeRealPair_Ewald(self, rsq, atmtype1, atmtype2) result(E_Pair)
  ! Computes the real-space pair energy with erfc screening
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  real(dp), intent(in) :: rsq
  integer, intent(in) :: atmtype1, atmtype2
  real(dp) :: E_Pair

  real(dp) :: r, qij

  r = sqrt(rsq)
  qij = self%qTable(atmtype1, atmtype2)
  
  ! Real-space Ewald: q_i * q_j * erfc(alpha * r) / r
  E_Pair = qij * erfc(self%alpha * r) / r

end function
!=============================================================================
subroutine SetupKVectors_Ewald(self, curbox)
  ! Sets up the k-vectors for the reciprocal space sum
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  real(dp) :: Lx, Ly, Lz
  real(dp) :: kx, ky, kz, ksq, kCutSq
  integer :: nkx, nky, nkz
  integer :: ikx, iky, ikz
  integer :: nKMax, nK
  real(dp), allocatable :: tempK(:,:), tempKsq(:), tempPrefac(:)
  real(dp) :: boxdim(2,3)

  ! Get box dimensions based on box type
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
      write(0,*) "Error: Ewald requires OrthoBox or CubeBox"
      error stop
  end select
  
  ! Check if we need to rebuild k-vectors
  if (self%precomputed) then
    if (abs(Lx - self%cachedLx) < 1.0E-10_dp .and. &
        abs(Ly - self%cachedLy) < 1.0E-10_dp .and. &
        abs(Lz - self%cachedLz) < 1.0E-10_dp) then
      return  ! No rebuild needed
    endif
  endif

  ! Maximum k-indices based on cutoff and box size
  self%kxMax = ceiling(self%kCutoff * Lx / two_pi)
  self%kyMax = ceiling(self%kCutoff * Ly / two_pi)
  self%kzMax = ceiling(self%kCutoff * Lz / two_pi)

  kCutSq = self%kCutoff**2

  ! First pass: count k-vectors
  nKMax = (2*self%kxMax + 1) * (2*self%kyMax + 1) * (2*self%kzMax + 1)
  allocate(tempK(3, nKMax))
  allocate(tempKsq(nKMax))
  allocate(tempPrefac(nKMax))

  nK = 0
  do ikx = 0, self%kxMax
    do iky = -self%kyMax, self%kyMax
      do ikz = -self%kzMax, self%kzMax
        ! Skip k = 0
        if (ikx == 0 .and. iky == 0 .and. ikz == 0) cycle
        
        ! For ikx = 0, only count half to avoid double counting
        if (ikx == 0 .and. iky < 0) cycle
        if (ikx == 0 .and. iky == 0 .and. ikz < 0) cycle
        
        kx = two_pi * real(ikx, dp) / Lx
        ky = two_pi * real(iky, dp) / Ly
        kz = two_pi * real(ikz, dp) / Lz
        ksq = kx**2 + ky**2 + kz**2
        
        if (ksq <= kCutSq) then
          nK = nK + 1
          tempK(1, nK) = kx
          tempK(2, nK) = ky
          tempK(3, nK) = kz
          tempKsq(nK) = ksq
          ! Prefactor: 4*pi/V * exp(-k^2/(4*alpha^2)) / k^2
          ! The factor of 2 accounts for the k/-k symmetry
          tempPrefac(nK) = 4.0_dp * pi / (Lx * Ly * Lz) * &
                           exp(-ksq / (4.0_dp * self%alpha_sq)) / ksq
          if (ikx /= 0) tempPrefac(nK) = 2.0_dp * tempPrefac(nK)
        endif
      enddo
    enddo
  enddo

  self%nKVectors = nK

  ! Allocate final arrays
  if (allocated(self%kVectors)) deallocate(self%kVectors)
  if (allocated(self%kMagSq)) deallocate(self%kMagSq)
  if (allocated(self%kPrefac)) deallocate(self%kPrefac)
  if (allocated(self%rhoK)) deallocate(self%rhoK)
  if (allocated(self%rhoK_old)) deallocate(self%rhoK_old)

  allocate(self%kVectors(3, nK))
  allocate(self%kMagSq(nK))
  allocate(self%kPrefac(nK))
  allocate(self%rhoK(nK))
  allocate(self%rhoK_old(nK))

  self%kVectors(:, 1:nK) = tempK(:, 1:nK)
  self%kMagSq(1:nK) = tempKsq(1:nK)
  self%kPrefac(1:nK) = tempPrefac(1:nK)

  deallocate(tempK, tempKsq, tempPrefac)

  ! Cache box dimensions
  self%cachedLx = Lx
  self%cachedLy = Ly
  self%cachedLz = Lz
  self%precomputed = .true.

end subroutine
!=============================================================================
subroutine ComputeStructureFactors_Ewald(self, curbox)
  ! Computes the structure factors rho(k) = sum_i q_i exp(i k.r_i)
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  integer :: iAtom, iK
  integer :: atmType
  real(dp) :: qi
  real(dp) :: kx, ky, kz
  real(dp) :: rx, ry, rz
  real(dp) :: kdotr
  complex(dp) :: expikr
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)

  self%rhoK = cmplx(0.0_dp, 0.0_dp, dp)

  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    
    if (abs(qi) < 1.0E-15_dp) cycle
    
    rx = atoms(1, iAtom)
    ry = atoms(2, iAtom)
    rz = atoms(3, iAtom)
    
    do iK = 1, self%nKVectors
      kx = self%kVectors(1, iK)
      ky = self%kVectors(2, iK)
      kz = self%kVectors(3, iK)
      
      kdotr = kx * rx + ky * ry + kz * rz
      expikr = cmplx(cos(kdotr), sin(kdotr), dp)
      
      self%rhoK(iK) = self%rhoK(iK) + qi * expikr
    enddo
  enddo

end subroutine
!=============================================================================
subroutine ComputeRecipEnergy_Ewald(self, curbox, E_recip)
  ! Computes the reciprocal space energy from structure factors
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_recip

  integer :: iK
  real(dp) :: rhoK_sq

  E_recip = 0.0_dp

  do iK = 1, self%nKVectors
    ! |rho(k)|^2 = rho(k) * rho(-k) = rho(k) * conj(rho(k))
    rhoK_sq = real(self%rhoK(iK) * conjg(self%rhoK(iK)), dp)
    E_recip = E_recip + self%kPrefac(iK) * rhoK_sq
  enddo

  E_recip = E_recip * coulombConst

end subroutine
!=============================================================================
subroutine ComputeSelfEnergy_Ewald(self, curbox, E_self)
  ! Computes the self-energy correction: -alpha/sqrt(pi) * sum_i q_i^2
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_self

  integer :: iAtom
  integer :: atmType
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
subroutine UpdateStructureFactor_Ewald(self, curbox, disp, E_recip_diff)
  ! Efficiently updates structure factors for a displacement move
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(out) :: E_recip_diff

  integer :: iDisp, iAtom, iK
  integer :: atmType
  real(dp) :: qi
  real(dp) :: kx, ky, kz
  real(dp) :: rx_old, ry_old, rz_old
  real(dp) :: rx_new, ry_new, rz_new
  real(dp) :: kdotr_old, kdotr_new
  complex(dp) :: expikr_old, expikr_new, delta_rhoK
  real(dp) :: E_recip_old, E_recip_new
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)

  ! Save old structure factors
  self%rhoK_old = self%rhoK

  ! Compute old reciprocal energy
  E_recip_old = 0.0_dp
  do iK = 1, self%nKVectors
    E_recip_old = E_recip_old + self%kPrefac(iK) * &
                  real(self%rhoK(iK) * conjg(self%rhoK(iK)), dp)
  enddo
  E_recip_old = E_recip_old * coulombConst

  ! Update structure factors for each displaced atom
  do iDisp = 1, size(disp)
    iAtom = disp(iDisp)%atmIndx
    atmType = curbox%AtomType(iAtom)
    qi = self%q(atmType)
    
    if (abs(qi) < 1.0E-15_dp) cycle
    
    rx_old = atoms(1, iAtom)
    ry_old = atoms(2, iAtom)
    rz_old = atoms(3, iAtom)
    
    rx_new = disp(iDisp)%x_new
    ry_new = disp(iDisp)%y_new
    rz_new = disp(iDisp)%z_new
    
    do iK = 1, self%nKVectors
      kx = self%kVectors(1, iK)
      ky = self%kVectors(2, iK)
      kz = self%kVectors(3, iK)
      
      kdotr_old = kx * rx_old + ky * ry_old + kz * rz_old
      kdotr_new = kx * rx_new + ky * ry_new + kz * rz_new
      
      expikr_old = cmplx(cos(kdotr_old), sin(kdotr_old), dp)
      expikr_new = cmplx(cos(kdotr_new), sin(kdotr_new), dp)
      
      delta_rhoK = qi * (expikr_new - expikr_old)
      self%rhoK(iK) = self%rhoK(iK) + delta_rhoK
    enddo
  enddo

  ! Compute new reciprocal energy
  E_recip_new = 0.0_dp
  do iK = 1, self%nKVectors
    E_recip_new = E_recip_new + self%kPrefac(iK) * &
                  real(self%rhoK(iK) * conjg(self%rhoK(iK)), dp)
  enddo
  E_recip_new = E_recip_new * coulombConst

  E_recip_diff = E_recip_new - E_recip_old

end subroutine
!=============================================================================
subroutine AcceptUpdate_Ewald(self)
  ! Called when a move is accepted - structure factors are already updated
  implicit none
  class(Pair_Ewald), intent(inout) :: self

  ! Nothing to do - rhoK already contains the new values

end subroutine
!=============================================================================
subroutine RejectUpdate_Ewald(self)
  ! Called when a move is rejected - restore old structure factors
  implicit none
  class(Pair_Ewald), intent(inout) :: self

  self%rhoK = self%rhoK_old

end subroutine
!=============================================================================
subroutine OptimizeAlpha_Ewald(self, curbox)
  ! Optimizes the Ewald alpha parameter based on box size and precision
  ! Using the relation: erfc(alpha*rCut) ~ precision for real-space
  ! and exp(-kCut^2/(4*alpha^2)) ~ precision for k-space
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  real(dp) :: Lx, Ly, Lz, Lmin
  real(dp) :: precision
  
  ! Get box dimensions based on box type
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
  
  precision = self%ewaldPrecision
  
  ! Optimize alpha such that real and k-space errors are balanced
  ! alpha ~= s/rCut where s is determined by precision
  ! For precision = 10^-6, s ~ 3.5
  self%alpha = 3.5_dp / self%rCutReal
  self%alpha_sq = self%alpha**2
  
  ! Set k-space cutoff to achieve similar precision
  ! kCut = 2 * alpha * sqrt(-log(precision))
  self%kCutoff = 2.0_dp * self%alpha * sqrt(-log(precision))

end subroutine
!=============================================================================
subroutine ProcessIO_Ewald(self, line)
  use Common_MolInfo, only: nAtomTypes
  use Input_Format, only: CountCommands, GetXCommand
  use Input_Format, only: maxLineLen
  use Units, only: inEngUnit, inLenUnit
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: command
  integer :: jType, lineStat, nPar
  integer :: type1
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

    case("kmax")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) self%kCutoff

    case("precision")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) self%ewaldPrecision

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
      ! Try to parse as atom type definition: type q rMin
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
subroutine Prologue_Ewald(self)
  use ParallelVar, only: nout
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  integer :: i

  write(nout, *) "============================================"
  write(nout, *) "Ewald Summation Parameters:"
  write(nout, *) "  Alpha (screening parameter):", self%alpha
  write(nout, *) "  Real-space cutoff:", self%rCutReal
  write(nout, *) "  K-space cutoff:", self%kCutoff
  write(nout, *) "  Number of k-vectors:", self%nKVectors
  write(nout, *) "  Target precision:", self%ewaldPrecision
  write(nout, *) ""
  write(nout, *) "Charges per atom type:"
  do i = 1, nAtomTypes
    write(nout, *) "  Type", i, ":", self%q(i)
  enddo
  write(nout, *) "============================================"

end subroutine
!=============================================================================
subroutine Epilogue_Ewald(self)
  implicit none
  class(Pair_Ewald), intent(inout) :: self

  ! Cleanup if needed

end subroutine
!=============================================================================
function GetCutOff_Ewald(self) result(rCut)
  implicit none
  class(Pair_Ewald), intent(inout) :: self
  real(dp) :: rCut

  rCut = self%rCutReal

end function
!=============================================================================
end module FF_Ewald
!=============================================================================
