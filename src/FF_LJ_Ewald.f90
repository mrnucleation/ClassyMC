!================================================================================
! Combined Lennard-Jones + Ewald Electrostatics Forcefield
!
! This module combines 12-6 Lennard-Jones interactions with Ewald summation
! for long-range electrostatics. This is appropriate for simulations of
! polar/ionic systems with periodic boundary conditions.
!
! The total energy is:
!   E_total = E_LJ + E_real + E_recip + E_self
!
! Where:
!   E_LJ    - Standard 12-6 LJ with cutoff
!   E_real  - Coulomb real-space with erfc screening
!   E_recip - Coulomb k-space (Ewald sum)
!   E_self  - Self-energy correction
!================================================================================
module FF_LJ_Ewald
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use ClassyConstants, only: coulombConst, pi, two_pi
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox

  implicit none

  type, extends(forcefield) :: Pair_LJ_Ewald
    ! LJ parameters
    real(dp), allocatable :: eps(:)
    real(dp), allocatable :: sig(:)
    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)   ! Stored as sigma^2
    real(dp) :: rLJCut, rLJCutSq
    
    ! Electrostatic parameters
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: qTable(:,:)
    
    ! Ewald parameters
    real(dp) :: alpha          
    real(dp) :: alpha_sq       
    real(dp) :: rEwaldCut      
    real(dp) :: rEwaldCutSq    
    real(dp) :: kCutoff        
    real(dp) :: ewaldPrecision = 1.0E-6_dp
    
    ! K-vector data
    integer :: nKVectors
    integer :: kxMax, kyMax, kzMax
    real(dp), allocatable :: kVectors(:,:)
    real(dp), allocatable :: kMagSq(:)
    real(dp), allocatable :: kPrefac(:)
    
    ! Structure factors
    complex(dp), allocatable :: rhoK(:)
    complex(dp), allocatable :: rhoK_old(:)
    
    ! Energy components
    real(dp) :: E_recip_total
    real(dp) :: E_self_total
    
    ! Minimum distance
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    
    ! Flags and cache
    logical :: precomputed = .false.
    real(dp) :: cachedLx, cachedLy, cachedLz
    
  contains
    procedure, pass :: Constructor => Constructor_LJ_Ewald
    procedure, pass :: DetailedECalc => Detailed_LJ_Ewald
    procedure, pass :: DiffECalc => DiffECalc_LJ_Ewald
    procedure, pass :: ShiftECalc_Single => Shift_LJ_Ewald_Single
    procedure, pass :: NewECalc => New_LJ_Ewald
    procedure, pass :: OldECalc => Old_LJ_Ewald
    procedure, pass :: ManyBody => ManyBody_LJ_Ewald
    procedure, pass :: ProcessIO => ProcessIO_LJ_Ewald
    procedure, pass :: Prologue => Prologue_LJ_Ewald
    procedure, pass :: GetCutOff => GetCutOff_LJ_Ewald
    
    ! Internal procedures
    procedure, pass :: ComputePairEnergy => ComputePairEnergy_LJ_Ewald
    procedure, pass :: SetupKVectors => SetupKVectors_LJ_Ewald
    procedure, pass :: ComputeStructureFactors => ComputeStructureFactors_LJ_Ewald
    procedure, pass :: ComputeRecipEnergy => ComputeRecipEnergy_LJ_Ewald
    procedure, pass :: ComputeSelfEnergy => ComputeSelfEnergy_LJ_Ewald
    procedure, pass :: UpdateStructureFactor => UpdateStructureFactor_LJ_Ewald
    
  end type

contains

!=============================================================================
subroutine Constructor_LJ_Ewald(self)
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  integer :: AllocateStat

  allocate(self%eps(1:nAtomTypes), stat=AllocateStat)
  allocate(self%sig(1:nAtomTypes), stat=AllocateStat)
  allocate(self%q(1:nAtomTypes), stat=AllocateStat)
  allocate(self%rMin(1:nAtomTypes), stat=AllocateStat)
  
  allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  allocate(self%sigTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  allocate(self%qTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)

  self%eps = 1.0_dp
  self%sig = 1.0_dp
  self%q = 0.0_dp
  self%rMin = 0.5_dp
  
  self%epsTable = 4.0_dp
  self%sigTable = 1.0_dp
  self%qTable = 0.0_dp
  self%rMinTable = 0.25_dp
  
  ! Default cutoffs
  self%rLJCut = 10.0_dp
  self%rLJCutSq = self%rLJCut**2
  self%rEwaldCut = 10.0_dp
  self%rEwaldCutSq = self%rEwaldCut**2
  self%rCut = max(self%rLJCut, self%rEwaldCut)
  self%rCutSq = self%rCut**2
  
  ! Default Ewald parameters
  self%alpha = 0.3_dp
  self%alpha_sq = self%alpha**2
  self%kCutoff = 5.0_dp
  
  self%precomputed = .false.
  self%nKVectors = 0

  if (AllocateStat /= 0) error stop "Allocation error in LJ_Ewald Constructor"

end subroutine
!=============================================================================
subroutine Detailed_LJ_Ewald(self, curbox, E_T, accept)
  use ParallelVar, only: nout
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(inout) :: E_T
  logical, intent(out) :: accept
  
  integer :: iAtom, jAtom
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_LJ, E_real, E_recip, E_self
  real(dp) :: E_pair_LJ, E_pair_Q
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)
  
  accept = .true.
  E_LJ = 0.0_dp
  E_real = 0.0_dp
  curbox%ETable = 0.0_dp
  
  ! Setup k-vectors if needed
  call self%SetupKVectors(curbox)
  
  ! Compute structure factors
  call self%ComputeStructureFactors(curbox)
  
  ! Real-space (LJ + screened Coulomb)
  do iAtom = 1, curbox%nMaxAtoms - 1
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType1 = curbox%AtomType(iAtom)
    
    do jAtom = iAtom + 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx**2 + ry**2 + rz**2
      
      atmType2 = curbox%AtomType(jAtom)
      rmin_ij = self%rMinTable(atmType1, atmType2)
      
      if (rsq < rmin_ij) then
        write(0,*) "ERROR! Overlapping atoms in LJ_Ewald"
        accept = .false.
        return
      endif
      
      call self%ComputePairEnergy(rsq, atmType1, atmType2, E_pair_LJ, E_pair_Q)
      
      E_LJ = E_LJ + E_pair_LJ
      E_real = E_real + E_pair_Q
      
      curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_pair_LJ + E_pair_Q
      curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_pair_LJ + E_pair_Q
    enddo
  enddo
  
  ! Reciprocal-space
  call self%ComputeRecipEnergy(curbox, E_recip)
  self%E_recip_total = E_recip
  
  ! Self-energy
  call self%ComputeSelfEnergy(curbox, E_self)
  self%E_self_total = E_self
  
  write(nout,*) "LJ Energy:", E_LJ
  write(nout,*) "Coulomb Real-space:", E_real
  write(nout,*) "Coulomb Reciprocal:", E_recip
  write(nout,*) "Coulomb Self-correction:", E_self
  write(nout,*) "Total Energy:", E_LJ + E_real + E_recip + E_self
  
  E_T = E_LJ + E_real + E_recip + E_self

end subroutine
!=============================================================================
subroutine DiffECalc_LJ_Ewald(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
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
      write(0,*) "Unknown Perturbation Type in LJ_Ewald"
      error stop
  end select

end subroutine
!=============================================================================
subroutine Shift_LJ_Ewald_Single(self, curbox, disp, E_Diff, accept)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jNei, jAtom, dispLen
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_pair_old, E_real_old, E_LJ_old
  real(dp) :: E_pair_new, E_real_new, E_LJ_new
  real(dp) :: E_pair_LJ, E_pair_Q
  real(dp) :: E_recip_diff
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
  E_LJ_old = 0.0_dp
  E_LJ_new = 0.0_dp
  accept = .true.

  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)

    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)

      ! New position
      rx = disp(iDisp)%x_new - atoms(1, jAtom)
      ry = disp(iDisp)%y_new - atoms(2, jAtom)
      rz = disp(iDisp)%z_new - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutSq) then
        rmin_ij = self%rMinTable(atmType1, atmType2)
        if (rsq < rmin_ij) then
          accept = .false.
          return
        endif
        
        call self%ComputePairEnergy(rsq, atmType1, atmType2, E_pair_LJ, E_pair_Q)
        E_pair_new = E_pair_LJ + E_pair_Q
        E_LJ_new = E_LJ_new + E_pair_LJ
        E_real_new = E_real_new + E_pair_Q
        
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair_new
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair_new
      endif

      ! Old position
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutSq) then
        call self%ComputePairEnergy(rsq, atmType1, atmType2, E_pair_LJ, E_pair_Q)
        E_pair_old = E_pair_LJ + E_pair_Q
        E_LJ_old = E_LJ_old + E_pair_LJ
        E_real_old = E_real_old + E_pair_Q
        
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair_old
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair_old
      endif
    enddo
  enddo

  ! Update reciprocal space
  call self%UpdateStructureFactor(curbox, disp, E_recip_diff)
  
  E_Diff = (E_LJ_new - E_LJ_old) + (E_real_new - E_real_old) + E_recip_diff

end subroutine
!=============================================================================
subroutine New_LJ_Ewald(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Addition), intent(in) :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
  integer :: atmType1, atmType2
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_LJ, E_real, E_self_diff
  real(dp) :: E_pair_LJ, E_pair_Q
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)

  dispLen = size(disp)
  E_Diff = 0.0_dp
  E_LJ = 0.0_dp
  E_real = 0.0_dp
  E_self_diff = 0.0_dp
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
      
      if (rsq < self%rCutSq) then
        atmType2 = curbox%AtomType(jAtom)
        rmin_ij = self%rMinTable(atmType1, atmType2)
        
        if (rsq < rmin_ij) then
          accept = .false.
          return
        endif
        
        call self%ComputePairEnergy(rsq, atmType1, atmType2, E_pair_LJ, E_pair_Q)
        E_LJ = E_LJ + E_pair_LJ
        E_real = E_real + E_pair_Q
        
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_pair_LJ + E_pair_Q
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + E_pair_LJ + E_pair_Q
      endif
    enddo
    
    ! Self-energy for new particle
    E_self_diff = E_self_diff - self%q(atmType1)**2 * coulombConst * self%alpha / sqrt(pi)
  enddo

  E_Diff = E_LJ + E_real + E_self_diff

end subroutine
!=============================================================================
subroutine Old_LJ_Ewald(self, curbox, disp, E_Diff)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Deletion), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff

  integer :: iAtom, jAtom, jNei
  integer :: atmType1, atmType2
  integer :: molEnd, molStart
  real(dp) :: rx, ry, rz, rsq
  real(dp) :: E_LJ, E_real, E_self_diff
  real(dp) :: E_pair_LJ, E_pair_Q
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()

  call curbox%GetCoordinates(atoms)
  call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)

  E_Diff = 0.0_dp
  E_LJ = 0.0_dp
  E_real = 0.0_dp
  E_self_diff = 0.0_dp

  call curBox%GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

  do iAtom = molStart, molEnd
    atmType1 = curbox%AtomType(iAtom)
    
    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)

      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutSq) then
        call self%ComputePairEnergy(rsq, atmType1, atmType2, E_pair_LJ, E_pair_Q)
        E_LJ = E_LJ - E_pair_LJ
        E_real = E_real - E_pair_Q
        
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_pair_LJ - E_pair_Q
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - E_pair_LJ - E_pair_Q
      endif
    enddo
    
    ! Self-energy for removed particle
    E_self_diff = E_self_diff + self%q(atmType1)**2 * coulombConst * self%alpha / sqrt(pi)
  enddo

  E_Diff = E_LJ + E_real + E_self_diff

end subroutine
!=============================================================================
subroutine ManyBody_LJ_Ewald(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
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
  real(dp) :: E_pair_LJ, E_pair_Q
  real(dp) :: rmin_ij

  accept = .true.
  E_Many = 0.0_dp
  
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
    
    if (rsq < self%rCutSq) then
      call self%ComputePairEnergy(rsq, atmtype1, atmType2, E_pair_LJ, E_pair_Q)
      E_Many = E_Many + E_pair_LJ + E_pair_Q
    endif
  enddo

end subroutine
!=============================================================================
subroutine ComputePairEnergy_LJ_Ewald(self, rsq, atmtype1, atmtype2, E_LJ, E_Q)
  ! Computes both LJ and real-space Coulomb for a pair
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  real(dp), intent(in) :: rsq
  integer, intent(in) :: atmtype1, atmtype2
  real(dp), intent(out) :: E_LJ, E_Q

  real(dp) :: r, ep, sig_sq, LJ6, qij

  E_LJ = 0.0_dp
  E_Q = 0.0_dp
  
  ! LJ contribution
  if (rsq < self%rLJCutSq) then
    ep = self%epsTable(atmtype1, atmtype2)
    if (abs(ep) > 1.0E-15_dp) then
      sig_sq = self%sigTable(atmtype1, atmtype2)
      LJ6 = (sig_sq / rsq)**3
      E_LJ = ep * LJ6 * (LJ6 - 1.0_dp)
    endif
  endif
  
  ! Real-space Coulomb with erfc screening
  if (rsq < self%rEwaldCutSq) then
    qij = self%qTable(atmtype1, atmtype2)
    if (abs(qij) > 1.0E-15_dp) then
      r = sqrt(rsq)
      E_Q = qij * erfc(self%alpha * r) / r
    endif
  endif

end subroutine
!=============================================================================
subroutine SetupKVectors_LJ_Ewald(self, curbox)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  real(dp) :: Lx, Ly, Lz
  real(dp) :: kx, ky, kz, ksq, kCutSq
  integer :: ikx, iky, ikz
  integer :: nKMax, nK
  real(dp), allocatable :: tempK(:,:), tempKsq(:), tempPrefac(:)

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
      write(0,*) "Error: LJ_Ewald requires OrthoBox or CubeBox"
      error stop
  end select
  
  ! Check if rebuild needed
  if (self%precomputed) then
    if (abs(Lx - self%cachedLx) < 1.0E-10_dp .and. &
        abs(Ly - self%cachedLy) < 1.0E-10_dp .and. &
        abs(Lz - self%cachedLz) < 1.0E-10_dp) then
      return
    endif
  endif

  self%kxMax = ceiling(self%kCutoff * Lx / two_pi)
  self%kyMax = ceiling(self%kCutoff * Ly / two_pi)
  self%kzMax = ceiling(self%kCutoff * Lz / two_pi)

  kCutSq = self%kCutoff**2

  nKMax = (2*self%kxMax + 1) * (2*self%kyMax + 1) * (2*self%kzMax + 1)
  allocate(tempK(3, nKMax))
  allocate(tempKsq(nKMax))
  allocate(tempPrefac(nKMax))

  nK = 0
  do ikx = 0, self%kxMax
    do iky = -self%kyMax, self%kyMax
      do ikz = -self%kzMax, self%kzMax
        if (ikx == 0 .and. iky == 0 .and. ikz == 0) cycle
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
          tempPrefac(nK) = 4.0_dp * pi / (Lx * Ly * Lz) * &
                           exp(-ksq / (4.0_dp * self%alpha_sq)) / ksq
          if (ikx /= 0) tempPrefac(nK) = 2.0_dp * tempPrefac(nK)
        endif
      enddo
    enddo
  enddo

  self%nKVectors = nK

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

  self%cachedLx = Lx
  self%cachedLy = Ly
  self%cachedLz = Lz
  self%precomputed = .true.

end subroutine
!=============================================================================
subroutine ComputeStructureFactors_LJ_Ewald(self, curbox)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  integer :: iAtom, iK
  integer :: atmType
  real(dp) :: qi, kx, ky, kz, rx, ry, rz, kdotr
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
subroutine ComputeRecipEnergy_LJ_Ewald(self, curbox, E_recip)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(out) :: E_recip

  integer :: iK
  real(dp) :: rhoK_sq

  E_recip = 0.0_dp

  do iK = 1, self%nKVectors
    rhoK_sq = real(self%rhoK(iK) * conjg(self%rhoK(iK)), dp)
    E_recip = E_recip + self%kPrefac(iK) * rhoK_sq
  enddo

  E_recip = E_recip * coulombConst

end subroutine
!=============================================================================
subroutine ComputeSelfEnergy_LJ_Ewald(self, curbox, E_self)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
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
subroutine UpdateStructureFactor_LJ_Ewald(self, curbox, disp, E_recip_diff)
  ! Compute the change in reciprocal-space energy for a trial displacement
  ! 
  ! Since there's no accept/reject callback mechanism in the forcefield interface,
  ! we take a conservative approach: recompute structure factors from scratch
  ! at the start of each call. This is slower O(N*nK) but guarantees correctness.
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(out) :: E_recip_diff

  integer :: iDisp, iAtom, iK, atmType
  real(dp) :: qi, kx, ky, kz
  real(dp) :: rx_old, ry_old, rz_old
  real(dp) :: rx_new, ry_new, rz_new
  real(dp) :: kdotr_old, kdotr_new
  complex(dp) :: expikr_old, expikr_new, delta_rhoK
  real(dp) :: E_recip_old, E_recip_new
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)

  ! Recompute structure factors from current (accepted) positions
  ! This ensures rhoK is always correct even after previous accepted moves
  call self%ComputeStructureFactors(curbox)

  ! Compute current reciprocal energy
  E_recip_old = 0.0_dp
  do iK = 1, self%nKVectors
    E_recip_old = E_recip_old + self%kPrefac(iK) * &
                  real(self%rhoK(iK) * conjg(self%rhoK(iK)), dp)
  enddo
  E_recip_old = E_recip_old * coulombConst

  ! Save current structure factors before trial update
  self%rhoK_old = self%rhoK

  ! Tentatively update structure factors for trial position
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

  ! Compute trial reciprocal energy
  E_recip_new = 0.0_dp
  do iK = 1, self%nKVectors
    E_recip_new = E_recip_new + self%kPrefac(iK) * &
                  real(self%rhoK(iK) * conjg(self%rhoK(iK)), dp)
  enddo
  E_recip_new = E_recip_new * coulombConst

  E_recip_diff = E_recip_new - E_recip_old

  ! Revert structure factors - they'll be recomputed fresh next time
  self%rhoK = self%rhoK_old

end subroutine
!=============================================================================
subroutine ProcessIO_LJ_Ewald(self, line)
  use Common_MolInfo, only: nAtomTypes
  use Input_Format, only: CountCommands, GetXCommand
  use Input_Format, only: maxLineLen
  use Units, only: inEngUnit, inLenUnit
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: command
  logical :: param = .false.
  integer :: jType, lineStat, nPar
  integer :: type1
  real(dp) :: ep, sig, q, rCut, rMin, alpha

  call GetXCommand(line, command, 1, lineStat)
  call CountCommands(line, nPar)

  select case(trim(adjustl(command)))
    case("ljrcut")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) rCut
      self%rLJCut = rCut * inLenUnit
      self%rLJCutSq = self%rLJCut**2
      self%rCut = max(self%rLJCut, self%rEwaldCut)
      self%rCutSq = self%rCut**2

    case("ewaldrcut")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) rCut
      self%rEwaldCut = rCut * inLenUnit
      self%rEwaldCutSq = self%rEwaldCut**2
      self%rCut = max(self%rLJCut, self%rEwaldCut)
      self%rCutSq = self%rCut**2

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

    case default
      param = .true.
  end select

  if (param) then
    select case(nPar)
      case(5)
        ! type ep sig q rMin
        read(line, *) type1, ep, sig, q, rMin
        ep = ep * inEngUnit
        sig = sig * inLenUnit
        rMin = rMin * inLenUnit

        self%eps(type1) = ep
        self%sig(type1) = sig
        self%q(type1) = q
        self%rMin(type1) = rMin

        do jType = 1, nAtomTypes
          self%epsTable(type1, jType) = 4.0_dp * sqrt(ep * self%eps(jType))
          self%epsTable(jType, type1) = 4.0_dp * sqrt(ep * self%eps(jType))
          
          self%sigTable(type1, jType) = (0.5_dp * (sig + self%sig(jType)))**2
          self%sigTable(jType, type1) = (0.5_dp * (sig + self%sig(jType)))**2
          
          self%qTable(type1, jType) = q * self%q(jType) * coulombConst
          self%qTable(jType, type1) = q * self%q(jType) * coulombConst
          
          self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
          self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
        enddo
        
      case default
        ! Unknown command
    end select
  endif

end subroutine
!=============================================================================
subroutine Prologue_LJ_Ewald(self)
  use ParallelVar, only: nout
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  integer :: i

  write(nout, *) "============================================"
  write(nout, *) "LJ + Ewald Forcefield Parameters:"
  write(nout, *) "  LJ Cutoff:", self%rLJCut
  write(nout, *) "  Ewald Alpha:", self%alpha
  write(nout, *) "  Ewald Real Cutoff:", self%rEwaldCut
  write(nout, *) "  Ewald K Cutoff:", self%kCutoff
  write(nout, *) "  Number of k-vectors:", self%nKVectors
  write(nout, *) ""
  write(nout, *) "Parameters per atom type:"
  do i = 1, nAtomTypes
    write(nout, "(A,I3,A,F10.5,A,F10.5,A,F10.5)") "  Type", i, &
          ": eps=", self%eps(i), " sig=", self%sig(i), " q=", self%q(i)
  enddo
  write(nout, *) "============================================"

end subroutine
!=============================================================================
function GetCutOff_LJ_Ewald(self) result(rCut)
  implicit none
  class(Pair_LJ_Ewald), intent(inout) :: self
  real(dp) :: rCut

  rCut = self%rCut

end function
!=============================================================================
end module FF_LJ_Ewald
!=============================================================================
