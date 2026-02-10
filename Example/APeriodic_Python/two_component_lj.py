"""
Two-Component Lennard-Jones Energy Function for ClassyMC

This Python module implements a 12-6 Lennard-Jones potential for a binary
mixture with the ClassyMC Monte Carlo simulation package.

The system has two atom types with different LJ parameters:
  - Type 1 (Ar): Larger atoms with stronger interactions
  - Type 2 (Ne): Smaller atoms with weaker interactions

LJ Potential: U(r) = 4*epsilon_ij * [(sigma_ij/r)^12 - (sigma_ij/r)^6]

Lorentz-Berthelot mixing rules:
  sigma_ij = (sigma_i + sigma_j) / 2
  epsilon_ij = sqrt(epsilon_i * epsilon_j)
"""

import numpy as np

# =============================================================================
# LJ Parameters (in reduced units)
# =============================================================================

# Type 1: Ar-like atoms (larger, stronger)
epsilon_1 = 1.0
sigma_1 = 1.0

# Type 2: Ne-like atoms (smaller, weaker)  
epsilon_2 = 0.5
sigma_2 = 0.8

# Cutoff and overlap detection
rcut = 6.0
rcut_sq = rcut * rcut
rmin_sq = 0.7 * 0.7  # Minimum allowed distance squared

# Precompute mixing parameters (Lorentz-Berthelot rules)
# sigma_ij and epsilon_ij for all pairs: (1,1), (1,2), (2,2)
sigma_11 = sigma_1
sigma_12 = (sigma_1 + sigma_2) / 2.0
sigma_22 = sigma_2

epsilon_11 = epsilon_1
epsilon_12 = np.sqrt(epsilon_1 * epsilon_2)
epsilon_22 = epsilon_2

# Precompute sigma^2 for efficiency
sig_sq_11 = sigma_11 * sigma_11
sig_sq_12 = sigma_12 * sigma_12
sig_sq_22 = sigma_22 * sigma_22

# Statistics tracking
_n_total_calls = 0
_n_diff_calls = 0


# =============================================================================
# Required Functions
# =============================================================================

def prologue(boxlist):
    """
    Initialize the two-component LJ energy function.
    Called once at the start of the simulation.
    """
    global _n_total_calls, _n_diff_calls
    
    print("=" * 70)
    print("Python Two-Component LJ Energy Module Initialized")
    print("=" * 70)
    print(f"  Type 1 (Ar): epsilon={epsilon_1:.3f}, sigma={sigma_1:.3f}")
    print(f"  Type 2 (Ne): epsilon={epsilon_2:.3f}, sigma={sigma_2:.3f}")
    print(f"  Cutoff:      {rcut:.3f}")
    print(f"  Rmin:        {np.sqrt(rmin_sq):.3f}")
    print()
    print("  Mixing parameters (Lorentz-Berthelot):")
    print(f"    sigma_11 = {sigma_11:.3f},  epsilon_11 = {epsilon_11:.3f}")
    print(f"    sigma_12 = {sigma_12:.3f},  epsilon_12 = {epsilon_12:.3f}")
    print(f"    sigma_22 = {sigma_22:.3f},  epsilon_22 = {epsilon_22:.3f}")
    
    if len(boxlist) > 0:
        box = boxlist[0]
        print()
        print(f"  Box type: {box['boxtype']}")
        print(f"  Volume:   {box['volume']:.2f}")
        
        # Count active atoms by type
        mol_types = box['moltype']
        mol_subindx = box['molsubindx']
        nmol = box['moleculecount']
        
        n_type1 = 0
        n_type2 = 0
        for i in range(box['raw_atoms'].shape[1]):
            if mol_subindx[i] <= nmol[mol_types[i] - 1]:
                if mol_types[i] == 1:
                    n_type1 += 1
                elif mol_types[i] == 2:
                    n_type2 += 1
        
        print(f"  Atoms:    {n_type1 + n_type2} total")
        print(f"    Type 1: {n_type1}")
        print(f"    Type 2: {n_type2}")
    
    print("=" * 70)
    
    _n_total_calls = 0
    _n_diff_calls = 0


def compute_total(boxlist):
    """
    Compute the total LJ energy of the system.
    
    Returns:
        dict: {'energy': total_energy, 'accept': bool}
    """
    global _n_total_calls
    _n_total_calls += 1
    
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    atom_types = box['atomtype']
    
    # Get box dimensions for PBC
    box_dim = box.get('boxdimensions')
    has_pbc = box_dim is not None and box['boxtype'] != 'nobox'
    if has_pbc:
        lx = box_dim[1, 0] - box_dim[0, 0]
        ly = box_dim[1, 1] - box_dim[0, 1]
        lz = box_dim[1, 2] - box_dim[0, 2]
    
    nmax = atoms.shape[1]
    total_energy = 0.0
    
    # Build list of active atoms
    active = []
    for i in range(nmax):
        if mol_subindx[i] <= nmol[mol_types[i] - 1]:
            active.append(i)
    
    # Double loop over pairs
    for ii, i in enumerate(active):
        for j in active[ii+1:]:
            # Distance vector
            rx = atoms[0, i] - atoms[0, j]
            ry = atoms[1, i] - atoms[1, j]
            rz = atoms[2, i] - atoms[2, j]
            
            # Minimum image convention (only for periodic dimensions)
            # For aperiodic_ortho: X and Y are periodic, Z is not
            if has_pbc and box['boxtype'] == 'aperiodic_ortho':
                rx = rx - lx * round(rx / lx)
                ry = ry - ly * round(ry / ly)
                # Z is not periodic, so no wrapping
            elif has_pbc:
                # Full periodic boundary conditions
                rx = rx - lx * round(rx / lx)
                ry = ry - ly * round(ry / ly)
                rz = rz - lz * round(rz / lz)
            
            rsq = rx*rx + ry*ry + rz*rz
            
            # Check for overlap
            if rsq < rmin_sq:
                print(f"WARNING: Overlap detected between atoms {i} and {j}, r={np.sqrt(rsq):.4f}")
                return {'energy': 1e10, 'accept': False}
            
            # LJ energy within cutoff
            if rsq < rcut_sq:
                # Determine interaction parameters based on atom types
                type_i = atom_types[i]
                type_j = atom_types[j]
                
                if type_i == 1 and type_j == 1:
                    sig_sq = sig_sq_11
                    eps = epsilon_11
                elif type_i == 2 and type_j == 2:
                    sig_sq = sig_sq_22
                    eps = epsilon_22
                else:  # Mixed interaction (1,2) or (2,1)
                    sig_sq = sig_sq_12
                    eps = epsilon_12
                
                sr2 = sig_sq / rsq
                sr6 = sr2 * sr2 * sr2
                sr12 = sr6 * sr6
                lj = 4.0 * eps * (sr12 - sr6)
                total_energy += lj
    
    return {'energy': total_energy, 'accept': True}


def compute_diff(boxlist, displist):
    """
    Compute the energy difference for a proposed move.
    
    This is called for each MC move attempt. For efficiency, we only
    compute interactions involving the moved atom(s).
    
    Returns:
        dict: {'energy': delta_energy, 'accept': bool}
    """
    global _n_diff_calls
    _n_diff_calls += 1
    
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    atom_types = box['atomtype']
    
    # Get box dimensions for PBC
    box_dim = box.get('boxdimensions')
    has_pbc = box_dim is not None and box['boxtype'] != 'nobox'
    if has_pbc:
        lx = box_dim[1, 0] - box_dim[0, 0]
        ly = box_dim[1, 1] - box_dim[0, 1]
        lz = box_dim[1, 2] - box_dim[0, 2]
    
    nmax = atoms.shape[1]
    delta_energy = 0.0
    
    for disp in displist:
        disp_type = disp.get('type', 'displacement')
        
        if disp_type == 'displacement':
            # Single atom displacement
            i = disp.get('atomindex', 0) - 1  # Convert to 0-indexed
            x_new = disp.get('x_new', atoms[0, i])
            y_new = disp.get('y_new', atoms[1, i])
            z_new = disp.get('z_new', atoms[2, i])
            
            x_old = atoms[0, i]
            y_old = atoms[1, i]
            z_old = atoms[2, i]
            
            type_i = atom_types[i]
            
            # Loop over all other active atoms
            for j in range(nmax):
                if j == i:
                    continue
                if mol_subindx[j] > nmol[mol_types[j] - 1]:
                    continue
                
                type_j = atom_types[j]
                
                # Determine interaction parameters
                if type_i == 1 and type_j == 1:
                    sig_sq = sig_sq_11
                    eps = epsilon_11
                elif type_i == 2 and type_j == 2:
                    sig_sq = sig_sq_22
                    eps = epsilon_22
                else:
                    sig_sq = sig_sq_12
                    eps = epsilon_12
                
                # --- New position energy ---
                rx = x_new - atoms[0, j]
                ry = y_new - atoms[1, j]
                rz = z_new - atoms[2, j]
                
                if has_pbc and box['boxtype'] == 'aperiodic_ortho':
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    # Z is not periodic
                elif has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq_new = rx*rx + ry*ry + rz*rz
                
                # Check for overlap
                if rsq_new < rmin_sq:
                    return {'energy': 0.0, 'accept': False}
                
                # New energy
                e_new = 0.0
                if rsq_new < rcut_sq:
                    sr2 = sig_sq / rsq_new
                    sr6 = sr2 * sr2 * sr2
                    e_new = 4.0 * eps * sr6 * (sr6 - 1.0)
                
                # --- Old position energy ---
                rx = x_old - atoms[0, j]
                ry = y_old - atoms[1, j]
                rz = z_old - atoms[2, j]
                
                if has_pbc and box['boxtype'] == 'aperiodic_ortho':
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                elif has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq_old = rx*rx + ry*ry + rz*rz
                
                e_old = 0.0
                if rsq_old < rcut_sq:
                    sr2 = sig_sq / rsq_old
                    sr6 = sr2 * sr2 * sr2
                    e_old = 4.0 * eps * sr6 * (sr6 - 1.0)
                
                delta_energy += e_new - e_old
        
        elif disp_type == 'addition':
            # Adding a new atom
            i = disp.get('atomindex', 0) - 1
            x_new = disp.get('x_new', 0.0)
            y_new = disp.get('y_new', 0.0)
            z_new = disp.get('z_new', 0.0)
            type_i = atom_types[i]
            
            for j in range(nmax):
                if j == i:
                    continue
                if mol_subindx[j] > nmol[mol_types[j] - 1]:
                    continue
                
                type_j = atom_types[j]
                
                if type_i == 1 and type_j == 1:
                    sig_sq = sig_sq_11
                    eps = epsilon_11
                elif type_i == 2 and type_j == 2:
                    sig_sq = sig_sq_22
                    eps = epsilon_22
                else:
                    sig_sq = sig_sq_12
                    eps = epsilon_12
                
                rx = x_new - atoms[0, j]
                ry = y_new - atoms[1, j]
                rz = z_new - atoms[2, j]
                
                if has_pbc and box['boxtype'] == 'aperiodic_ortho':
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                elif has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq = rx*rx + ry*ry + rz*rz
                
                if rsq < rmin_sq:
                    return {'energy': 0.0, 'accept': False}
                
                if rsq < rcut_sq:
                    sr2 = sig_sq / rsq
                    sr6 = sr2 * sr2 * sr2
                    delta_energy += 4.0 * eps * sr6 * (sr6 - 1.0)
        
        elif disp_type == 'deletion':
            # Removing an atom - subtract its current interactions
            i = disp.get('atomindex', 0) - 1
            type_i = atom_types[i]
            
            for j in range(nmax):
                if j == i:
                    continue
                if mol_subindx[j] > nmol[mol_types[j] - 1]:
                    continue
                
                type_j = atom_types[j]
                
                if type_i == 1 and type_j == 1:
                    sig_sq = sig_sq_11
                    eps = epsilon_11
                elif type_i == 2 and type_j == 2:
                    sig_sq = sig_sq_22
                    eps = epsilon_22
                else:
                    sig_sq = sig_sq_12
                    eps = epsilon_12
                
                rx = atoms[0, i] - atoms[0, j]
                ry = atoms[1, i] - atoms[1, j]
                rz = atoms[2, i] - atoms[2, j]
                
                if has_pbc and box['boxtype'] == 'aperiodic_ortho':
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                elif has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq = rx*rx + ry*ry + rz*rz
                
                if rsq < rcut_sq:
                    sr2 = sig_sq / rsq
                    sr6 = sr2 * sr2 * sr2
                    delta_energy -= 4.0 * eps * sr6 * (sr6 - 1.0)
    
    return {'energy': delta_energy, 'accept': True}


def update(boxlist):
    """
    Called after accepted moves (optional).
    Can be used for bookkeeping or state updates.
    """
    pass


# =============================================================================
# Helper Functions
# =============================================================================

def get_statistics():
    """
    Get call statistics for the energy functions.
    
    Returns:
        dict: Statistics about function calls
    """
    return {
        'total_calls': _n_total_calls,
        'diff_calls': _n_diff_calls
    }
