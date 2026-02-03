"""
Lennard-Jones Energy Function for ClassyMC

This Python module implements a 12-6 Lennard-Jones potential for use
with the ClassyMC Monte Carlo simulation package.

LJ Potential: U(r) = 4*epsilon * [(sigma/r)^12 - (sigma/r)^6]

Parameters can be modified at runtime using the set_parameters() function.
"""

import numpy as np

# =============================================================================
# LJ Parameters (in reduced units: epsilon=1, sigma=1)
# =============================================================================

epsilon = 1.0      # Energy parameter
sigma = 1.0        # Length parameter  
rcut = 10.0        # Cutoff distance
rcut_sq = rcut * rcut
rmin_sq = 0.8 * 0.8  # Minimum allowed distance squared (overlap detection)

# Precomputed for efficiency
sig_sq = sigma * sigma

# Statistics tracking
_n_total_calls = 0
_n_diff_calls = 0


# =============================================================================
# Required Functions
# =============================================================================

def prologue(boxlist):
    """
    Initialize the LJ energy function.
    Called once at the start of the simulation.
    """
    global _n_total_calls, _n_diff_calls
    
    print("=" * 60)
    print("Python LJ Energy Module Initialized")
    print("=" * 60)
    print(f"  Epsilon:  {epsilon}")
    print(f"  Sigma:    {sigma}")
    print(f"  Cutoff:   {rcut}")
    print(f"  Rmin:     {np.sqrt(rmin_sq)}")
    
    if len(boxlist) > 0:
        box = boxlist[0]
        print(f"  Box type: {box['boxtype']}")
        print(f"  Volume:   {box['volume']:.2f}")
        
        # Count active atoms
        mol_types = box['moltype']
        mol_subindx = box['molsubindx']
        nmol = box['moleculecount']
        n_active = sum(1 for i in range(box['raw_atoms'].shape[1]) 
                      if mol_subindx[i] <= nmol[mol_types[i] - 1])
        print(f"  Atoms:    {n_active}")
    
    print("=" * 60)
    
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
            
            # Minimum image convention
            if has_pbc:
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
                sr2 = sig_sq / rsq
                sr6 = sr2 * sr2 * sr2
                sr12 = sr6 * sr6
                lj = 4.0 * epsilon * (sr12 - sr6)
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
            
            # Loop over all other active atoms
            for j in range(nmax):
                if j == i:
                    continue
                if mol_subindx[j] > nmol[mol_types[j] - 1]:
                    continue
                
                # --- New position energy ---
                rx = x_new - atoms[0, j]
                ry = y_new - atoms[1, j]
                rz = z_new - atoms[2, j]
                
                if has_pbc:
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
                    e_new = 4.0 * epsilon * sr6 * (sr6 - 1.0)
                
                # --- Old position energy ---
                rx = x_old - atoms[0, j]
                ry = y_old - atoms[1, j]
                rz = z_old - atoms[2, j]
                
                if has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq_old = rx*rx + ry*ry + rz*rz
                
                e_old = 0.0
                if rsq_old < rcut_sq:
                    sr2 = sig_sq / rsq_old
                    sr6 = sr2 * sr2 * sr2
                    e_old = 4.0 * epsilon * sr6 * (sr6 - 1.0)
                
                delta_energy += e_new - e_old
        
        elif disp_type == 'addition':
            # Adding a new atom
            i = disp.get('atomindex', 0) - 1
            x_new = disp.get('x_new', 0.0)
            y_new = disp.get('y_new', 0.0)
            z_new = disp.get('z_new', 0.0)
            
            for j in range(nmax):
                if j == i:
                    continue
                if mol_subindx[j] > nmol[mol_types[j] - 1]:
                    continue
                
                rx = x_new - atoms[0, j]
                ry = y_new - atoms[1, j]
                rz = z_new - atoms[2, j]
                
                if has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq = rx*rx + ry*ry + rz*rz
                
                if rsq < rmin_sq:
                    return {'energy': 0.0, 'accept': False}
                
                if rsq < rcut_sq:
                    sr2 = sig_sq / rsq
                    sr6 = sr2 * sr2 * sr2
                    delta_energy += 4.0 * epsilon * sr6 * (sr6 - 1.0)
        
        elif disp_type == 'deletion':
            # Removing an atom - subtract its current interactions
            i = disp.get('atomindex', 0) - 1
            
            for j in range(nmax):
                if j == i:
                    continue
                if mol_subindx[j] > nmol[mol_types[j] - 1]:
                    continue
                
                rx = atoms[0, i] - atoms[0, j]
                ry = atoms[1, i] - atoms[1, j]
                rz = atoms[2, i] - atoms[2, j]
                
                if has_pbc:
                    rx = rx - lx * round(rx / lx)
                    ry = ry - ly * round(ry / ly)
                    rz = rz - lz * round(rz / lz)
                
                rsq = rx*rx + ry*ry + rz*rz
                
                if rsq < rcut_sq:
                    sr2 = sig_sq / rsq
                    sr6 = sr2 * sr2 * sr2
                    delta_energy -= 4.0 * epsilon * sr6 * (sr6 - 1.0)
    
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

def set_parameters(eps=None, sig=None, rc=None, rmin=None):
    """
    Set LJ parameters at runtime.
    
    Args:
        eps: Epsilon (energy parameter)
        sig: Sigma (length parameter)
        rc: Cutoff distance
        rmin: Minimum distance for overlap detection
    """
    global epsilon, sigma, rcut, rcut_sq, rmin_sq, sig_sq
    
    if eps is not None:
        epsilon = eps
    if sig is not None:
        sigma = sig
        sig_sq = sigma * sigma
    if rc is not None:
        rcut = rc
        rcut_sq = rcut * rcut
    if rmin is not None:
        rmin_sq = rmin * rmin
    
    print(f"LJ parameters updated: eps={epsilon}, sig={sigma}, rcut={rcut}")


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
