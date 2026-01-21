!================================================================================
! Embedded Atom Method (EAM) Potential for Monte Carlo Simulations
!
! This module implements the EAM potential commonly used for metallic systems.
! The implementation is optimized for Monte Carlo with efficient differential
! energy calculations that track changes in electron density.
!
! The EAM energy is:
!   E = Σ_i F_i(ρ_i) + (1/2) Σ_i Σ_j φ_ij(r_ij)
!
! Where:
!   ρ_i = Σ_j ρ_j(r_ij)  - electron density at atom i
!   F_i(ρ)               - embedding function
!   φ_ij(r)              - pair potential
!   ρ_j(r)               - electron density contribution from atom type j
!
! Supports LAMMPS-style EAM file formats:
!   - setfl (alloy format)
!   - funcfl (single element format)
!
! MC Efficiency:
!   When atom i moves, we need to update:
!   1. ρ_i (density at moving atom changes)
!   2. ρ_j for all neighbors j (they lose/gain contribution from i)
!   3. F(ρ) for atom i and all affected neighbors
!   4. Pair interactions involving atom i
!
! References:
!   - Daw & Baskes, Phys. Rev. B 29, 6443 (1984)
!   - Foiles, Baskes & Daw, Phys. Rev. B 33, 7983 (1986)
!================================================================================
module FF_EAM
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  implicit none
  private

  ! EAM tabulated function type
  type :: EAM_Table
    integer :: nPoints           ! Number of table points
    real(dp) :: rMin, rMax       ! Range of table
    real(dp) :: dr               ! Spacing
    real(dp) :: drInv            ! 1/dr for fast lookup
    real(dp), allocatable :: values(:)    ! Table values
    real(dp), allocatable :: derivs(:)    ! Pre-computed derivatives for interpolation
  end type

  ! EAM parameters for each atom type pair
  type :: EAM_PairData
    type(EAM_Table) :: phi       ! Pair potential φ(r)
    type(EAM_Table) :: rho       ! Electron density ρ(r)
  end type

  ! EAM embedding function for each atom type  
  type :: EAM_EmbedData
    type(EAM_Table) :: F         ! Embedding function F(ρ)
  end type

  type, public, extends(forcefield) :: Pair_EAM
    ! Number of types and cutoff
    integer :: nEAMTypes = 0
    real(dp) :: rCutEAM = 6.0_dp
    real(dp) :: rCutEAMSq = 36.0_dp
    
    ! Unit conversion factor for energy (EAM files are typically in eV)
    ! Default: 1 eV = 1/8.6173303E-5 kb ≈ 11604.5 kb
    real(dp) :: eV_to_kb = 11604.522_dp  ! 1/8.6173303E-5
    real(dp) :: energyConversion = 11604.522_dp  ! Applied to all energies
    
    ! Mapping from simulation atom types to EAM types
    integer, allocatable :: typeMap(:)
    
    ! EAM tables
    type(EAM_PairData), allocatable :: pairData(:,:)   ! (type1, type2)
    type(EAM_EmbedData), allocatable :: embedData(:)   ! (type)
    
    ! Cached electron densities for each atom (for efficient MC updates)
    real(dp), allocatable :: rhoAtom(:)        ! Current density at each atom
    real(dp), allocatable :: rhoAtom_old(:)    ! Backup for rejection
    
    ! Temporary storage for differential calculations
    integer, allocatable :: affectedAtoms(:)   ! List of atoms with changed ρ
    integer :: nAffected = 0
    real(dp), allocatable :: rhoChange(:)      ! Change in ρ for affected atoms
    
    ! Minimum distance for overlap
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    
    ! Flags
    logical :: tablesLoaded = .false.
    logical :: densitiesComputed = .false.
    
  contains
    procedure, pass :: Constructor => Constructor_EAM
    procedure, pass :: DetailedECalc => Detailed_EAM
    procedure, pass :: DiffECalc => DiffECalc_EAM
    procedure, pass :: ShiftECalc_Single => Shift_EAM_Single
    procedure, pass :: NewECalc => New_EAM
    procedure, pass :: OldECalc => Old_EAM
    procedure, pass :: ManyBody => ManyBody_EAM
    procedure, pass :: ProcessIO => ProcessIO_EAM
    procedure, pass :: Prologue => Prologue_EAM
    procedure, pass :: GetCutOff => GetCutOff_EAM
    
    ! EAM-specific procedures
    procedure, pass :: LoadSetflFile => LoadSetflFile_EAM
    procedure, pass :: LoadFuncflFile => LoadFuncflFile_EAM
    procedure, pass :: LoadLAMMPSTable => LoadLAMMPSTable_EAM
    procedure, pass :: InterpolateTable => InterpolateTable_EAM
    procedure, pass :: InterpolateTableDeriv => InterpolateTableDeriv_EAM
    procedure, pass :: ComputeAllDensities => ComputeAllDensities_EAM
    procedure, pass :: ComputeAtomDensity => ComputeAtomDensity_EAM
    procedure, pass :: GetPairEnergy => GetPairEnergy_EAM
    procedure, pass :: GetEmbedEnergy => GetEmbedEnergy_EAM
    procedure, pass :: GetRhoContribution => GetRhoContribution_EAM
    procedure, pass :: AcceptMove => AcceptMove_EAM
    procedure, pass :: RejectMove => RejectMove_EAM
    procedure, pass :: CreateAnalyticalEAM => CreateAnalyticalEAM_EAM
    
  end type

contains

!=============================================================================
subroutine Constructor_EAM(self)
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_EAM), intent(inout) :: self
  integer :: AllocateStat, i

  allocate(self%rMin(1:nAtomTypes), stat=AllocateStat)
  if (AllocateStat /= 0) error stop "Allocation error in EAM Constructor (rMin)"
  
  allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat=AllocateStat)
  if (AllocateStat /= 0) error stop "Allocation error in EAM Constructor (rMinTable)"
  
  allocate(self%typeMap(1:nAtomTypes), stat=AllocateStat)
  if (AllocateStat /= 0) error stop "Allocation error in EAM Constructor (typeMap)"

  self%rMin = 0.5_dp
  self%rMinTable = 0.25_dp  ! Stored as r^2
  
  ! Default: 1-to-1 mapping
  do i = 1, nAtomTypes
    self%typeMap(i) = i
  enddo
  
  self%rCutEAM = 6.0_dp
  self%rCutEAMSq = 36.0_dp
  self%rCut = self%rCutEAM
  self%rCutSq = self%rCutEAMSq
  
  self%tablesLoaded = .false.
  self%densitiesComputed = .false.
  self%nAffected = 0

end subroutine
!=============================================================================
subroutine Detailed_EAM(self, curbox, E_T, accept)
  use ParallelVar, only: nout
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(inout) :: E_T
  logical, intent(out) :: accept
  
  integer :: iAtom, jAtom
  integer :: atmType1, atmType2
  integer :: eamType1, eamType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_pair, E_embed, E_total
  real(dp) :: phi_ij, rho_contrib
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)
  
  accept = .true.
  E_pair = 0.0_dp
  E_embed = 0.0_dp
  curbox%ETable = 0.0_dp
  
  ! Allocate/resize density array if needed
  if (.not. allocated(self%rhoAtom)) then
    allocate(self%rhoAtom(curbox%nMaxAtoms))
    allocate(self%rhoAtom_old(curbox%nMaxAtoms))
    allocate(self%affectedAtoms(curbox%nMaxAtoms))
    allocate(self%rhoChange(curbox%nMaxAtoms))
  endif
  
  ! Compute electron densities for all atoms
  call self%ComputeAllDensities(curbox)
  
  ! Compute embedding energy
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    
    atmType1 = curbox%AtomType(iAtom)
    eamType1 = self%typeMap(atmType1)
    
    E_embed = E_embed + self%GetEmbedEnergy(eamType1, self%rhoAtom(iAtom))
  enddo
  
  ! Compute pair energy
  do iAtom = 1, curbox%nMaxAtoms - 1
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType1 = curbox%AtomType(iAtom)
    eamType1 = self%typeMap(atmType1)
    
    do jAtom = iAtom + 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx**2 + ry**2 + rz**2
      
      if (rsq < self%rCutEAMSq) then
        atmType2 = curbox%AtomType(jAtom)
        eamType2 = self%typeMap(atmType2)
        rmin_ij = self%rMinTable(atmType1, atmType2)
        
        if (rsq < rmin_ij) then
          write(0,*) "ERROR! Overlapping atoms in EAM!"
          accept = .false.
          return
        endif
        
        r = sqrt(rsq)
        phi_ij = self%GetPairEnergy(eamType1, eamType2, r)
        E_pair = E_pair + phi_ij
        
        curbox%ETable(iAtom) = curbox%ETable(iAtom) + phi_ij
        curbox%ETable(jAtom) = curbox%ETable(jAtom) + phi_ij
      endif
    enddo
  enddo
  
  E_total = E_pair + E_embed
  
  write(nout,*) "============================================"
  write(nout,*) "EAM Energy Calculation Results:"
  write(nout,*) "  Energy conversion factor:", self%energyConversion
  write(nout,*) "  Pair Energy (kb):", E_pair
  write(nout,*) "  Embedding Energy (kb):", E_embed
  write(nout,*) "  Total Energy (kb):", E_total
  write(nout,*) "  Total Energy (eV):", E_total / self%eV_to_kb
  write(nout,*) "  Per-atom Energy (eV):", E_total / (self%eV_to_kb * curbox%nMaxAtoms)
  write(nout,*) "============================================"
  
  E_T = E_total
  self%densitiesComputed = .true.

end subroutine
!=============================================================================
subroutine DiffECalc_EAM(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_EAM), intent(inout) :: self
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
      write(0,*) "Unknown Perturbation Type in EAM DiffECalc"
      error stop
  end select

end subroutine
!=============================================================================
subroutine Shift_EAM_Single(self, curbox, disp, E_Diff, accept)
  ! Computes energy change when atoms are displaced
  ! This routine does NOT modify the stored density cache.
  ! Densities are recomputed fresh each time to ensure accuracy.
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jNei, jAtom, dispLen
  integer :: atmType1, atmType2
  integer :: eamType1, eamType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: rx_new, ry_new, rz_new, rsq_new, r_new
  real(dp) :: rmin_ij
  real(dp) :: phi_old, phi_new
  real(dp) :: rho_old_i, rho_new_i
  real(dp) :: rho_old_j, rho_new_j
  real(dp) :: rho_contrib_old, rho_contrib_new
  real(dp) :: F_old_i, F_new_i
  real(dp) :: F_old_j, F_new_j
  real(dp) :: E_pair_diff, E_embed_diff
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()
  
  ! Track which atoms are affected and their density changes
  integer :: affectedCount
  integer :: affectedList(500)
  real(dp) :: rhoOld_list(500)
  real(dp) :: rhoNew_list(500)

  call curbox%GetCoordinates(atoms)
  call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)

  dispLen = size(disp)
  E_Diff = 0.0_dp
  E_pair_diff = 0.0_dp
  E_embed_diff = 0.0_dp
  accept = .true.
  affectedCount = 0

  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)
    eamType1 = self%typeMap(atmType1)

    ! Compute current density at moving atom from scratch
    rho_old_i = 0.0_dp
    rho_new_i = 0.0_dp

    ! Loop over neighbors
    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)
      eamType2 = self%typeMap(atmType2)

      ! --- Old position calculations ---
      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutEAMSq) then
        r = sqrt(rsq)
        ! Old pair energy
        phi_old = self%GetPairEnergy(eamType1, eamType2, r)
        E_pair_diff = E_pair_diff - phi_old
        
        ! Old rho contribution from j to i
        rho_old_i = rho_old_i + self%GetRhoContribution(eamType2, r)
        ! Old rho contribution from i to j 
        rho_contrib_old = self%GetRhoContribution(eamType1, r)
      else
        phi_old = 0.0_dp
        rho_contrib_old = 0.0_dp
      endif

      ! --- New position calculations ---
      rx_new = disp(iDisp)%x_new - atoms(1, jAtom)
      ry_new = disp(iDisp)%y_new - atoms(2, jAtom)
      rz_new = disp(iDisp)%z_new - atoms(3, jAtom)
      call curbox%Boundary(rx_new, ry_new, rz_new)
      rsq_new = rx_new*rx_new + ry_new*ry_new + rz_new*rz_new
      
      if (rsq_new < self%rCutEAMSq) then
        rmin_ij = self%rMinTable(atmType1, atmType2)
        if (rsq_new < rmin_ij) then
          accept = .false.
          return
        endif
        
        r_new = sqrt(rsq_new)
        
        ! New pair energy
        phi_new = self%GetPairEnergy(eamType1, eamType2, r_new)
        E_pair_diff = E_pair_diff + phi_new
        
        ! New rho contribution from j to i
        rho_new_i = rho_new_i + self%GetRhoContribution(eamType2, r_new)
        
        ! New rho contribution from i to j
        rho_contrib_new = self%GetRhoContribution(eamType1, r_new)
      else
        phi_new = 0.0_dp
        rho_contrib_new = 0.0_dp
      endif

      ! Track change in density at neighbor j
      if (abs(rho_contrib_new - rho_contrib_old) > 1.0E-15_dp) then
        affectedCount = affectedCount + 1
        if (affectedCount > 500) then
          write(0,*) "ERROR: Too many affected atoms in EAM move"
          error stop
        endif
        affectedList(affectedCount) = jAtom
        ! Compute current density at neighbor j from scratch
        call self%ComputeAtomDensity(curbox, jAtom)
        rhoOld_list(affectedCount) = self%rhoAtom(jAtom)
        rhoNew_list(affectedCount) = self%rhoAtom(jAtom) + (rho_contrib_new - rho_contrib_old)
      endif
      
      ! Update dETable for pair contributions
      curbox%dETable(iAtom) = curbox%dETable(iAtom) + phi_new - phi_old
      curbox%dETable(jAtom) = curbox%dETable(jAtom) + phi_new - phi_old
    enddo

    ! Compute embedding energy change for moving atom
    F_old_i = self%GetEmbedEnergy(eamType1, rho_old_i)
    F_new_i = self%GetEmbedEnergy(eamType1, rho_new_i)
    E_embed_diff = E_embed_diff + (F_new_i - F_old_i)
  enddo

  ! Compute embedding energy change for all affected neighbors
  do jNei = 1, affectedCount
    jAtom = affectedList(jNei)
    atmType2 = curbox%AtomType(jAtom)
    eamType2 = self%typeMap(atmType2)
    
    F_old_j = self%GetEmbedEnergy(eamType2, rhoOld_list(jNei))
    F_new_j = self%GetEmbedEnergy(eamType2, rhoNew_list(jNei))
    
    E_embed_diff = E_embed_diff + (F_new_j - F_old_j)
  enddo

  E_Diff = E_pair_diff + E_embed_diff

end subroutine
!=============================================================================
subroutine New_EAM(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  ! Energy change for particle addition
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Addition), intent(in) :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
  integer :: atmType1, atmType2
  integer :: eamType1, eamType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_pair, E_embed_new, E_embed_change
  real(dp) :: phi_ij, rho_i, rho_contrib
  real(dp) :: rho_old_j, rho_new_j
  real(dp) :: F_old_j, F_new_j
  real(dp) :: rmin_ij
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)

  dispLen = size(disp)
  E_Diff = 0.0_dp
  E_pair = 0.0_dp
  E_embed_change = 0.0_dp
  accept = .true.

  ! Backup densities
  self%rhoAtom_old = self%rhoAtom

  do iDisp = 1, dispLen
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)
    eamType1 = self%typeMap(atmType1)

    listIndx = disp(iDisp)%listIndex
    maxNei = tempNNei(listIndx)
    
    rho_i = 0.0_dp  ! Density at new atom
    
    do jNei = 1, maxNei
      jAtom = tempList(jNei, listIndx)
      atmType2 = curbox%AtomType(jAtom)
      eamType2 = self%typeMap(atmType2)
      
      rx = disp(iDisp)%x_new - atoms(1, jAtom)
      ry = disp(iDisp)%y_new - atoms(2, jAtom)
      rz = disp(iDisp)%z_new - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutEAMSq) then
        rmin_ij = self%rMinTable(atmType1, atmType2)
        
        if (rsq < rmin_ij) then
          accept = .false.
          self%rhoAtom = self%rhoAtom_old
          return
        endif
        
        r = sqrt(rsq)
        
        ! Pair energy
        phi_ij = self%GetPairEnergy(eamType1, eamType2, r)
        E_pair = E_pair + phi_ij
        
        ! Density at new atom from neighbor j
        rho_i = rho_i + self%GetRhoContribution(eamType2, r)
        
        ! Density contribution from new atom to neighbor j
        rho_contrib = self%GetRhoContribution(eamType1, r)
        
        rho_old_j = self%rhoAtom(jAtom)
        rho_new_j = rho_old_j + rho_contrib
        
        F_old_j = self%GetEmbedEnergy(eamType2, rho_old_j)
        F_new_j = self%GetEmbedEnergy(eamType2, rho_new_j)
        
        E_embed_change = E_embed_change + (F_new_j - F_old_j)
        
        ! Update neighbor density
        self%rhoAtom(jAtom) = rho_new_j
        
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + phi_ij
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + phi_ij
      endif
    enddo
    
    ! Embedding energy for new atom
    E_embed_new = self%GetEmbedEnergy(eamType1, rho_i)
    self%rhoAtom(iAtom) = rho_i
  enddo

  E_Diff = E_pair + E_embed_new + E_embed_change

end subroutine
!=============================================================================
subroutine Old_EAM(self, curbox, disp, E_Diff)
  ! Energy change for particle deletion
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Deletion), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff

  integer :: iAtom, jAtom, jNei
  integer :: atmType1, atmType2
  integer :: eamType1, eamType2
  integer :: molEnd, molStart
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: E_pair, E_embed_old, E_embed_change
  real(dp) :: phi_ij, rho_contrib
  real(dp) :: rho_old_j, rho_new_j
  real(dp) :: F_old_j, F_new_j
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()

  call curbox%GetCoordinates(atoms)
  call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)

  E_Diff = 0.0_dp
  E_pair = 0.0_dp
  E_embed_change = 0.0_dp

  ! Backup densities
  self%rhoAtom_old = self%rhoAtom

  call curBox%GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

  do iAtom = molStart, molEnd
    atmType1 = curbox%AtomType(iAtom)
    eamType1 = self%typeMap(atmType1)
    
    ! Embedding energy of deleted atom (will be removed)
    E_embed_old = self%GetEmbedEnergy(eamType1, self%rhoAtom(iAtom))
    E_embed_change = E_embed_change - E_embed_old
    
    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)
      eamType2 = self%typeMap(atmType2)

      rx = atoms(1, iAtom) - atoms(1, jAtom)
      ry = atoms(2, iAtom) - atoms(2, jAtom)
      rz = atoms(3, iAtom) - atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      
      if (rsq < self%rCutEAMSq) then
        r = sqrt(rsq)
        
        ! Pair energy (removed)
        phi_ij = self%GetPairEnergy(eamType1, eamType2, r)
        E_pair = E_pair - phi_ij
        
        ! Density contribution that will be removed from neighbor j
        rho_contrib = self%GetRhoContribution(eamType1, r)
        
        rho_old_j = self%rhoAtom(jAtom)
        rho_new_j = rho_old_j - rho_contrib
        
        F_old_j = self%GetEmbedEnergy(eamType2, rho_old_j)
        F_new_j = self%GetEmbedEnergy(eamType2, rho_new_j)
        
        E_embed_change = E_embed_change + (F_new_j - F_old_j)
        
        ! Update neighbor density
        self%rhoAtom(jAtom) = rho_new_j
        
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - phi_ij
        curbox%dETable(jAtom) = curbox%dETable(jAtom) - phi_ij
      endif
    enddo
    
    ! Set density of deleted atom to zero
    self%rhoAtom(iAtom) = 0.0_dp
  enddo

  E_Diff = E_pair + E_embed_change

end subroutine
!=============================================================================
subroutine ManyBody_EAM(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept)
  ! For trial insertion moves (CBMC, etc.)
  ! Returns approximate energy (pair + embedding at trial position)
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(simBox), intent(inout) :: curbox
  integer, intent(in) :: atmtype1
  integer, intent(in) :: atmtypes(:)
  real(dp), intent(in) :: pos1(:)
  real(dp), intent(in) :: posN(:,:)
  logical, intent(out) :: accept
  real(dp), intent(out) :: E_Many

  integer :: jAtom
  integer :: eamType1, eamType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: phi_ij, rho_i
  real(dp) :: rmin_ij

  accept = .true.
  E_Many = 0.0_dp
  rho_i = 0.0_dp
  
  eamType1 = self%typeMap(atmtype1)
  
  do jAtom = 1, size(posN, 2)
    eamType2 = self%typeMap(atmtypes(jAtom))
    
    rx = pos1(1) - posN(1, jAtom)
    ry = pos1(2) - posN(2, jAtom)
    rz = pos1(3) - posN(3, jAtom)
    call curbox%Boundary(rx, ry, rz)
    rsq = rx*rx + ry*ry + rz*rz
    
    rmin_ij = self%rMinTable(eamType1, eamType2)
    if (rsq < rmin_ij) then
      accept = .false.
      return
    endif
    
    if (rsq < self%rCutEAMSq) then
      r = sqrt(rsq)
      phi_ij = self%GetPairEnergy(eamType1, eamType2, r)
      E_Many = E_Many + phi_ij
      
      rho_i = rho_i + self%GetRhoContribution(eamType2, r)
    endif
  enddo
  
  ! Add embedding energy
  E_Many = E_Many + self%GetEmbedEnergy(eamType1, rho_i)

end subroutine
!=============================================================================
subroutine ComputeAllDensities_EAM(self, curbox)
  ! Computes electron density at all atoms from scratch
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox

  integer :: iAtom
  
  self%rhoAtom = 0.0_dp
  
  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    call self%ComputeAtomDensity(curbox, iAtom)
  enddo

end subroutine
!=============================================================================
subroutine ComputeAtomDensity_EAM(self, curbox, iAtom)
  ! Computes electron density at a single atom
  implicit none
  class(Pair_EAM), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  integer, intent(in) :: iAtom

  integer :: jAtom
  integer :: atmType2, eamType2
  real(dp) :: rx, ry, rz, rsq, r
  real(dp) :: rho_sum
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)
  
  rho_sum = 0.0_dp
  
  do jAtom = 1, curbox%nMaxAtoms
    if (jAtom == iAtom) cycle
    if (.not. curbox%IsActive(jAtom)) cycle
    if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
    
    rx = atoms(1, iAtom) - atoms(1, jAtom)
    ry = atoms(2, iAtom) - atoms(2, jAtom)
    rz = atoms(3, iAtom) - atoms(3, jAtom)
    call curbox%Boundary(rx, ry, rz)
    rsq = rx*rx + ry*ry + rz*rz
    
    if (rsq < self%rCutEAMSq) then
      r = sqrt(rsq)
      atmType2 = curbox%AtomType(jAtom)
      eamType2 = self%typeMap(atmType2)
      rho_sum = rho_sum + self%GetRhoContribution(eamType2, r)
    endif
  enddo
  
  self%rhoAtom(iAtom) = rho_sum

end subroutine
!=============================================================================
function GetPairEnergy_EAM(self, type1, type2, r) result(phi)
  ! Returns pair energy φ(r) between types, converted to internal units
  implicit none
  class(Pair_EAM), intent(inout) :: self
  integer, intent(in) :: type1, type2
  real(dp), intent(in) :: r
  real(dp) :: phi

  if (r >= self%pairData(type1, type2)%phi%rMax) then
    phi = 0.0_dp
  else
    phi = self%InterpolateTable(self%pairData(type1, type2)%phi, r)
    phi = phi * self%energyConversion
  endif

end function
!=============================================================================
function GetEmbedEnergy_EAM(self, atomType, rho) result(F)
  ! Returns embedding energy F(ρ) for atom type, converted to internal units
  implicit none
  class(Pair_EAM), intent(inout) :: self
  integer, intent(in) :: atomType
  real(dp), intent(in) :: rho
  real(dp) :: F

  if (rho <= 0.0_dp) then
    F = self%embedData(atomType)%F%values(1)
  elseif (rho >= self%embedData(atomType)%F%rMax) then
    ! Extrapolate linearly beyond table
    F = self%embedData(atomType)%F%values(self%embedData(atomType)%F%nPoints)
  else
    F = self%InterpolateTable(self%embedData(atomType)%F, rho)
  endif
  
  F = F * self%energyConversion

end function
!=============================================================================
function GetRhoContribution_EAM(self, atomType, r) result(rho)
  ! Returns electron density contribution ρ(r) from atom type
  implicit none
  class(Pair_EAM), intent(inout) :: self
  integer, intent(in) :: atomType
  real(dp), intent(in) :: r
  real(dp) :: rho

  if (r >= self%pairData(atomType, atomType)%rho%rMax) then
    rho = 0.0_dp
  else
    rho = self%InterpolateTable(self%pairData(atomType, atomType)%rho, r)
  endif

end function
!=============================================================================
function InterpolateTable_EAM(self, table, x) result(y)
  ! Cubic spline interpolation in table
  implicit none
  class(Pair_EAM), intent(in) :: self
  type(EAM_Table), intent(in) :: table
  real(dp), intent(in) :: x
  real(dp) :: y
  
  integer :: idx
  real(dp) :: t, t2, t3
  real(dp) :: y0, y1, d0, d1
  
  ! Find table index
  idx = int((x - table%rMin) * table%drInv) + 1
  
  ! Clamp to valid range
  if (idx < 1) idx = 1
  if (idx >= table%nPoints) idx = table%nPoints - 1
  
  ! Normalized position in interval
  t = (x - table%rMin - (idx-1) * table%dr) * table%drInv
  t2 = t * t
  t3 = t2 * t
  
  ! Hermite interpolation
  y0 = table%values(idx)
  y1 = table%values(idx+1)
  
  if (allocated(table%derivs)) then
    d0 = table%derivs(idx) * table%dr
    d1 = table%derivs(idx+1) * table%dr
    y = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*d0 + &
        (-2*t3 + 3*t2)*y1 + (t3 - t2)*d1
  else
    ! Linear interpolation fallback
    y = y0 + t * (y1 - y0)
  endif

end function
!=============================================================================
function InterpolateTableDeriv_EAM(self, table, x) result(dy)
  ! Returns derivative at x
  implicit none
  class(Pair_EAM), intent(in) :: self
  type(EAM_Table), intent(in) :: table
  real(dp), intent(in) :: x
  real(dp) :: dy
  
  integer :: idx
  real(dp) :: t
  
  idx = int((x - table%rMin) * table%drInv) + 1
  if (idx < 1) idx = 1
  if (idx >= table%nPoints) idx = table%nPoints - 1
  
  if (allocated(table%derivs)) then
    t = (x - table%rMin - (idx-1) * table%dr) * table%drInv
    dy = table%derivs(idx) * (1-t) + table%derivs(idx+1) * t
  else
    dy = (table%values(idx+1) - table%values(idx)) * table%drInv
  endif

end function
!=============================================================================
subroutine AcceptMove_EAM(self)
  ! Called when move is accepted - densities already updated
  implicit none
  class(Pair_EAM), intent(inout) :: self
  
  ! Nothing to do - rhoAtom already contains new values
  
end subroutine
!=============================================================================
subroutine RejectMove_EAM(self)
  ! Called when move is rejected - restore old densities
  implicit none
  class(Pair_EAM), intent(inout) :: self
  
  self%rhoAtom = self%rhoAtom_old
  
end subroutine
!=============================================================================
subroutine LoadSetflFile_EAM(self, filename)
  ! Loads LAMMPS setfl format EAM file
  ! Format: alloy multi-element EAM
  implicit none
  class(Pair_EAM), intent(inout) :: self
  character(len=*), intent(in) :: filename
  
  integer :: iUnit, ios
  integer :: i, j, k, nTypes
  integer :: nRho, nR
  real(dp) :: dRho, dR, rCut
  character(len=256) :: line
  character(len=3), allocatable :: elemNames(:)
  real(dp), allocatable :: mass(:), lattice(:)
  integer, allocatable :: atomNum(:)
  real(dp), allocatable :: tempData(:)
  
  iUnit = 42
  open(unit=iUnit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(0,*) "ERROR: Cannot open EAM file: ", trim(filename)
    error stop
  endif
  
  ! Skip 3 comment lines
  read(iUnit, '(A)') line
  read(iUnit, '(A)') line
  read(iUnit, '(A)') line
  
  ! Read number of elements
  read(iUnit, *) nTypes
  
  allocate(elemNames(nTypes))
  allocate(mass(nTypes))
  allocate(lattice(nTypes))
  allocate(atomNum(nTypes))
  
  ! Read element names (on same line or next)
  backspace(iUnit)
  read(iUnit, *) nTypes, (elemNames(i), i=1,nTypes)
  
  ! Read grid parameters
  read(iUnit, *) nRho, dRho, nR, dR, rCut
  
  self%nEAMTypes = nTypes
  self%rCutEAM = rCut
  self%rCutEAMSq = rCut * rCut
  self%rCut = rCut
  self%rCutSq = rCut * rCut
  
  ! Allocate tables
  allocate(self%embedData(nTypes))
  allocate(self%pairData(nTypes, nTypes))
  allocate(tempData(max(nRho, nR)))
  
  ! Read embedding functions and density functions for each element
  do i = 1, nTypes
    read(iUnit, *) atomNum(i), mass(i), lattice(i)
    
    ! Read F(rho) - embedding function
    self%embedData(i)%F%nPoints = nRho
    self%embedData(i)%F%rMin = 0.0_dp
    self%embedData(i)%F%rMax = (nRho - 1) * dRho
    self%embedData(i)%F%dr = dRho
    self%embedData(i)%F%drInv = 1.0_dp / dRho
    allocate(self%embedData(i)%F%values(nRho))
    read(iUnit, *) (self%embedData(i)%F%values(k), k=1,nRho)
    
    ! Read rho(r) - electron density function
    self%pairData(i,i)%rho%nPoints = nR
    self%pairData(i,i)%rho%rMin = 0.0_dp
    self%pairData(i,i)%rho%rMax = (nR - 1) * dR
    self%pairData(i,i)%rho%dr = dR
    self%pairData(i,i)%rho%drInv = 1.0_dp / dR
    allocate(self%pairData(i,i)%rho%values(nR))
    read(iUnit, *) (self%pairData(i,i)%rho%values(k), k=1,nR)
  enddo
  
  ! Read pair potentials phi(r) for all pairs
  do i = 1, nTypes
    do j = 1, i
      self%pairData(i,j)%phi%nPoints = nR
      self%pairData(i,j)%phi%rMin = 0.0_dp
      self%pairData(i,j)%phi%rMax = (nR - 1) * dR
      self%pairData(i,j)%phi%dr = dR
      self%pairData(i,j)%phi%drInv = 1.0_dp / dR
      allocate(self%pairData(i,j)%phi%values(nR))
      
      ! setfl stores r*phi(r), need to convert
      read(iUnit, *) (tempData(k), k=1,nR)
      do k = 1, nR
        if (k == 1) then
          self%pairData(i,j)%phi%values(k) = tempData(2) / dR  ! Extrapolate to r=0
        else
          self%pairData(i,j)%phi%values(k) = tempData(k) / ((k-1) * dR)
        endif
      enddo
      
      ! Symmetric
      if (i /= j) then
        self%pairData(j,i)%phi = self%pairData(i,j)%phi
        ! Copy rho function for cross terms
        self%pairData(i,j)%rho = self%pairData(j,j)%rho
        self%pairData(j,i)%rho = self%pairData(i,i)%rho
      endif
    enddo
  enddo
  
  close(iUnit)
  
  deallocate(elemNames, mass, lattice, atomNum, tempData)
  
  self%tablesLoaded = .true.

end subroutine
!=============================================================================
subroutine LoadFuncflFile_EAM(self, filename)
  ! Loads LAMMPS funcfl format EAM file (single element)
  implicit none
  class(Pair_EAM), intent(inout) :: self
  character(len=*), intent(in) :: filename
  
  integer :: iUnit, ios
  integer :: k, nRho, nR
  real(dp) :: dRho, dR, rCut, mass, lattice
  integer :: atomNum
  character(len=256) :: line
  real(dp), allocatable :: tempData(:)
  
  iUnit = 42
  open(unit=iUnit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(0,*) "ERROR: Cannot open EAM funcfl file: ", trim(filename)
    error stop
  endif
  
  ! Skip comment line
  read(iUnit, '(A)') line
  
  ! Read atomic info
  read(iUnit, *) atomNum, mass, lattice
  
  ! Read grid parameters
  read(iUnit, *) nRho, dRho, nR, dR, rCut
  
  self%nEAMTypes = 1
  self%rCutEAM = rCut
  self%rCutEAMSq = rCut * rCut
  self%rCut = rCut
  self%rCutSq = rCut * rCut
  
  ! Allocate tables for single element
  allocate(self%embedData(1))
  allocate(self%pairData(1, 1))
  allocate(tempData(max(nRho, nR)))
  
  ! Read F(rho) - embedding function
  self%embedData(1)%F%nPoints = nRho
  self%embedData(1)%F%rMin = 0.0_dp
  self%embedData(1)%F%rMax = (nRho - 1) * dRho
  self%embedData(1)%F%dr = dRho
  self%embedData(1)%F%drInv = 1.0_dp / dRho
  allocate(self%embedData(1)%F%values(nRho))
  read(iUnit, *) (self%embedData(1)%F%values(k), k=1,nRho)
  
  ! Read Z(r) - effective charge (used to compute phi)
  self%pairData(1,1)%phi%nPoints = nR
  self%pairData(1,1)%phi%rMin = 0.0_dp
  self%pairData(1,1)%phi%rMax = (nR - 1) * dR
  self%pairData(1,1)%phi%dr = dR
  self%pairData(1,1)%phi%drInv = 1.0_dp / dR
  allocate(self%pairData(1,1)%phi%values(nR))
  
  ! Read Z(r) and convert to phi(r) = Z(r)^2 / r
  read(iUnit, *) (tempData(k), k=1,nR)
  do k = 1, nR
    if (k == 1) then
      self%pairData(1,1)%phi%values(k) = tempData(2)**2 / dR
    else
      self%pairData(1,1)%phi%values(k) = tempData(k)**2 / ((k-1) * dR)
    endif
  enddo
  
  ! Read rho(r) - electron density function
  self%pairData(1,1)%rho%nPoints = nR
  self%pairData(1,1)%rho%rMin = 0.0_dp
  self%pairData(1,1)%rho%rMax = (nR - 1) * dR
  self%pairData(1,1)%rho%dr = dR
  self%pairData(1,1)%rho%drInv = 1.0_dp / dR
  allocate(self%pairData(1,1)%rho%values(nR))
  read(iUnit, *) (self%pairData(1,1)%rho%values(k), k=1,nR)
  
  close(iUnit)
  
  deallocate(tempData)
  
  self%tablesLoaded = .true.

end subroutine
!=============================================================================
subroutine LoadLAMMPSTable_EAM(self, filename)
  ! Loads a generic LAMMPS-style EAM table file
  ! This is similar to setfl format but allows for comments and flexible formatting
  ! Format:
  !   # Comment lines
  !   <Nelements> <element1> <element2> ...
  !   <Nrho> <drho> <Nr> <dr> <cutoff>
  !   For each element:
  !     <atomic_number> <mass> <lattice_constant> <lattice_type>
  !     F(rho) values (Nrho values)
  !     rho(r) values (Nr values)
  !   For each i,j pair (in lower triangular order):
  !     r*phi(r) values (Nr values)
  implicit none
  class(Pair_EAM), intent(inout) :: self
  character(len=*), intent(in) :: filename
  
  integer :: iUnit, ios
  integer :: i, j, k, nTypes
  integer :: nRho, nR
  real(dp) :: dRho, dR, rCut
  character(len=256) :: line
  character(len=3), allocatable :: elemNames(:)
  real(dp), allocatable :: mass(:), lattice(:)
  integer, allocatable :: atomNum(:)
  real(dp), allocatable :: tempData(:)
  character(len=10) :: lattType
  
  iUnit = 42
  open(unit=iUnit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(0,*) "ERROR: Cannot open LAMMPS EAM table file: ", trim(filename)
    error stop
  endif
  
  ! Skip comment lines starting with #
  do
    read(iUnit, '(A)', iostat=ios) line
    if (ios /= 0) then
      write(0,*) "ERROR: Unexpected end of EAM file"
      error stop
    endif
    line = adjustl(line)
    if (line(1:1) /= '#' .and. len_trim(line) > 0) exit
  enddo
  
  ! First non-comment line: number of elements and element names
  read(line, *) nTypes
  
  allocate(elemNames(nTypes))
  allocate(mass(nTypes))
  allocate(lattice(nTypes))
  allocate(atomNum(nTypes))
  
  ! Re-read line with element names
  backspace(iUnit)
  read(iUnit, *) nTypes, (elemNames(i), i=1,nTypes)
  
  ! Read grid parameters
  read(iUnit, *) nRho, dRho, nR, dR, rCut
  
  self%nEAMTypes = nTypes
  self%rCutEAM = rCut
  self%rCutEAMSq = rCut * rCut
  self%rCut = rCut
  self%rCutSq = rCut * rCut
  
  ! Allocate tables
  allocate(self%embedData(nTypes))
  allocate(self%pairData(nTypes, nTypes))
  allocate(tempData(max(nRho, nR)))
  
  ! Read embedding functions and density functions for each element
  do i = 1, nTypes
    read(iUnit, *) atomNum(i), mass(i), lattice(i), lattType
    
    ! Read F(rho) - embedding function
    self%embedData(i)%F%nPoints = nRho
    self%embedData(i)%F%rMin = 0.0_dp
    self%embedData(i)%F%rMax = (nRho - 1) * dRho
    self%embedData(i)%F%dr = dRho
    self%embedData(i)%F%drInv = 1.0_dp / dRho
    allocate(self%embedData(i)%F%values(nRho))
    read(iUnit, *) (self%embedData(i)%F%values(k), k=1,nRho)
    
    ! Read rho(r) - electron density function
    self%pairData(i,i)%rho%nPoints = nR
    self%pairData(i,i)%rho%rMin = 0.0_dp
    self%pairData(i,i)%rho%rMax = (nR - 1) * dR
    self%pairData(i,i)%rho%dr = dR
    self%pairData(i,i)%rho%drInv = 1.0_dp / dR
    allocate(self%pairData(i,i)%rho%values(nR))
    read(iUnit, *) (self%pairData(i,i)%rho%values(k), k=1,nR)
  enddo
  
  ! Read pair potentials phi(r) for all pairs
  do i = 1, nTypes
    do j = 1, i
      self%pairData(i,j)%phi%nPoints = nR
      self%pairData(i,j)%phi%rMin = 0.0_dp
      self%pairData(i,j)%phi%rMax = (nR - 1) * dR
      self%pairData(i,j)%phi%dr = dR
      self%pairData(i,j)%phi%drInv = 1.0_dp / dR
      allocate(self%pairData(i,j)%phi%values(nR))
      
      ! Table stores r*phi(r), need to convert
      read(iUnit, *) (tempData(k), k=1,nR)
      do k = 1, nR
        if (k == 1) then
          self%pairData(i,j)%phi%values(k) = tempData(2) / dR  ! Extrapolate to r=0
        else
          self%pairData(i,j)%phi%values(k) = tempData(k) / ((k-1) * dR)
        endif
      enddo
      
      ! Symmetric
      if (i /= j) then
        self%pairData(j,i)%phi = self%pairData(i,j)%phi
        ! Copy rho function for cross terms
        self%pairData(i,j)%rho = self%pairData(j,j)%rho
        self%pairData(j,i)%rho = self%pairData(i,i)%rho
      endif
    enddo
  enddo
  
  close(iUnit)
  
  deallocate(elemNames, mass, lattice, atomNum, tempData)
  
  self%tablesLoaded = .true.
  write(0,*) "Loaded LAMMPS EAM table with", nTypes, " element types"
  write(0,*) "  Cutoff:", rCut, " Angstrom"
  write(0,*) "  Grid: Nrho=", nRho, " Nr=", nR

end subroutine
!=============================================================================
subroutine ProcessIO_EAM(self, line)
  use Input_Format, only: CountCommands, GetXCommand
  use Input_Format, only: maxLineLen
  use Units, only: inLenUnit
  implicit none
  class(Pair_EAM), intent(inout) :: self
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: command
  character(len=256) :: filename
  integer :: lineStat, nPar
  integer :: type1, type2
  real(dp) :: rCut, rMin

  call GetXCommand(line, command, 1, lineStat)
  call CountCommands(line, nPar)

  select case(trim(adjustl(command)))
    case("rcut")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) rCut
      self%rCutEAM = rCut * inLenUnit
      self%rCutEAMSq = self%rCutEAM**2
      self%rCut = self%rCutEAM
      self%rCutSq = self%rCutEAMSq

    case("setfl", "alloy")
      call GetXCommand(line, filename, 2, lineStat)
      call self%LoadSetflFile(trim(filename))

    case("funcfl", "single")
      call GetXCommand(line, filename, 2, lineStat)
      call self%LoadFuncflFile(trim(filename))

    case("typemap")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) type2
      self%typeMap(type1) = type2

    case("rmin")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) rMin
      rMin = rMin * inLenUnit
      self%rMin(type1) = rMin
      self%rMinTable(type1, :) = max(rMin, self%rMin(:))**2
      self%rMinTable(:, type1) = max(rMin, self%rMin(:))**2

    case("analytical", "test")
      ! Create simple analytical EAM for testing
      ! Uses Morse-like pair and simple embedding
      call self%CreateAnalyticalEAM()

    case("units")
      ! Set energy units for the EAM tables
      ! Common options: ev (default), kcal-mol, kj-mol, kb
      call GetXCommand(line, command, 2, lineStat)
      select case(trim(adjustl(command)))
        case("ev", "eV")
          self%energyConversion = 1.0_dp / 8.6173303E-5_dp  ! eV to kb
        case("kcal-mol", "kcal/mol")
          self%energyConversion = 1.0_dp / 1.9872041E-3_dp  ! kcal/mol to kb
        case("kj-mol", "kj/mol")
          self%energyConversion = 1.0_dp / 8.3144621E-3_dp  ! kJ/mol to kb
        case("kb", "boltzmann")
          self%energyConversion = 1.0_dp  ! Already in kb
        case default
          write(0,*) "Unknown EAM energy unit:", trim(command)
      end select

    case("table")
      ! Load LAMMPS-style pair table file
      call GetXCommand(line, filename, 2, lineStat)
      call self%LoadLAMMPSTable(trim(filename))

    case default
      ! Unknown command - ignore
  end select

end subroutine
!=============================================================================
subroutine CreateAnalyticalEAM_EAM(self)
  ! Creates a simple analytical EAM potential for testing
  ! Uses Morse-like pair potential and square-root embedding
  ! Parameters calibrated to give Cu-like energies (~3.5 eV/atom cohesive)
  implicit none
  class(Pair_EAM), intent(inout) :: self
  
  integer :: i, nPoints
  real(dp) :: dr, dRho, r, rho
  real(dp) :: De, alpha, r0  ! Morse parameters
  real(dp) :: A_rho, r_rho   ! Density parameters
  real(dp) :: A_F            ! Embedding parameter
  
  ! Parameters calibrated for Cu-like metal
  ! Target: cohesive energy ~3.5 eV/atom at equilibrium FCC
  ! At equilibrium: 12 neighbors at r0 = 2.556 Å
  ! Pair energy per atom: 6 * phi(r0) ≈ -1.5 eV (half of pair sum)
  ! Embedding per atom: F(12*rho0) ≈ -2.0 eV
  ! Total: ~-3.5 eV/atom
  De = 0.25_dp      ! eV (pair well depth - 6 pairs gives 1.5 eV)
  alpha = 1.5_dp    ! 1/Angstrom (controls width)
  r0 = 2.556_dp     ! Angstrom (Cu FCC nearest neighbor distance)
  A_rho = 0.5_dp    ! Density prefactor (gives rho~4 at equilibrium)
  r_rho = 2.2_dp    ! Density decay length
  A_F = -1.0_dp     ! Embedding prefactor (sqrt(4)*12/6 ~ 2 eV)
  
  nPoints = 500
  dr = self%rCutEAM / (nPoints - 1)
  dRho = 20.0_dp / (nPoints - 1)  ! Max density ~20
  
  self%nEAMTypes = 1
  
  ! Allocate tables
  allocate(self%embedData(1))
  allocate(self%pairData(1, 1))
  
  ! Create embedding function F(rho) = A * sqrt(rho)
  self%embedData(1)%F%nPoints = nPoints
  self%embedData(1)%F%rMin = 0.0_dp
  self%embedData(1)%F%rMax = 20.0_dp
  self%embedData(1)%F%dr = dRho
  self%embedData(1)%F%drInv = 1.0_dp / dRho
  allocate(self%embedData(1)%F%values(nPoints))
  
  do i = 1, nPoints
    rho = (i - 1) * dRho
    if (rho > 0.0_dp) then
      self%embedData(1)%F%values(i) = A_F * sqrt(rho)
    else
      self%embedData(1)%F%values(i) = 0.0_dp
    endif
  enddo
  
  ! Create pair potential phi(r) - Morse
  self%pairData(1,1)%phi%nPoints = nPoints
  self%pairData(1,1)%phi%rMin = 0.0_dp
  self%pairData(1,1)%phi%rMax = self%rCutEAM
  self%pairData(1,1)%phi%dr = dr
  self%pairData(1,1)%phi%drInv = 1.0_dp / dr
  allocate(self%pairData(1,1)%phi%values(nPoints))
  
  do i = 1, nPoints
    r = (i - 1) * dr
    if (r > 0.5_dp .and. r < self%rCutEAM) then
      self%pairData(1,1)%phi%values(i) = De * (exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)))
    else if (r <= 0.5_dp) then
      self%pairData(1,1)%phi%values(i) = 1000.0_dp  ! Repulsive core
    else
      self%pairData(1,1)%phi%values(i) = 0.0_dp
    endif
  enddo
  
  ! Create density function rho(r) - exponential decay
  self%pairData(1,1)%rho%nPoints = nPoints
  self%pairData(1,1)%rho%rMin = 0.0_dp
  self%pairData(1,1)%rho%rMax = self%rCutEAM
  self%pairData(1,1)%rho%dr = dr
  self%pairData(1,1)%rho%drInv = 1.0_dp / dr
  allocate(self%pairData(1,1)%rho%values(nPoints))
  
  do i = 1, nPoints
    r = (i - 1) * dr
    if (r > 0.0_dp .and. r < self%rCutEAM) then
      self%pairData(1,1)%rho%values(i) = A_rho * exp(-(r/r_rho)**2)
    else
      self%pairData(1,1)%rho%values(i) = 0.0_dp
    endif
  enddo
  
  self%tablesLoaded = .true.

end subroutine
!=============================================================================
subroutine Prologue_EAM(self)
  use ParallelVar, only: nout
  implicit none
  class(Pair_EAM), intent(inout) :: self
  integer :: i
  real(dp) :: rTest, rhoTest, phiTest, FTest

  write(nout, *) "============================================"
  write(nout, *) "EAM (Embedded Atom Method) Parameters:"
  write(nout, *) "  Number of EAM types:", self%nEAMTypes
  write(nout, *) "  Cutoff:", self%rCutEAM, " Angstrom"
  write(nout, *) "  Energy conversion (eV->kb):", self%energyConversion
  write(nout, *) "  Tables loaded:", self%tablesLoaded
  if (self%tablesLoaded) then
    do i = 1, self%nEAMTypes
      write(nout, *) "  Type", i, ":"
      write(nout, *) "    F(rho) points:", self%embedData(i)%F%nPoints
      write(nout, *) "    F(rho) range: 0 to", self%embedData(i)%F%rMax
      write(nout, *) "    phi(r) points:", self%pairData(i,i)%phi%nPoints
      write(nout, *) "    rho(r) points:", self%pairData(i,i)%rho%nPoints
      
      ! Test values at typical FCC neighbor distance (~2.5 Å)
      rTest = 2.556_dp  ! Cu nearest neighbor distance
      if (rTest < self%pairData(i,i)%phi%rMax) then
        phiTest = self%GetPairEnergy(i, i, rTest)
        rhoTest = self%GetRhoContribution(i, rTest)
        FTest = self%GetEmbedEnergy(i, 12.0_dp * rhoTest)  ! Approx for 12 neighbors
        write(nout, *) "    At r=2.556 Å (FCC nn):"
        write(nout, *) "      phi(r) =", phiTest, " kb (", phiTest/self%energyConversion, " eV)"
        write(nout, *) "      rho(r) =", rhoTest
        write(nout, *) "      F(12*rho) =", FTest, " kb (", FTest/self%energyConversion, " eV)"
      endif
    enddo
  endif
  write(nout, *) "============================================"

end subroutine
!=============================================================================
function GetCutOff_EAM(self) result(rCut)
  implicit none
  class(Pair_EAM), intent(inout) :: self
  real(dp) :: rCut

  rCut = self%rCutEAM

end function
!=============================================================================
end module FF_EAM
!=============================================================================
