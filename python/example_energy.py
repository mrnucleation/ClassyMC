"""
Example Python Energy Script for ClassyMC

This script demonstrates how to create a custom energy function
using Python. The energy is computed from Fortran via the forpy interface.

To use this energy function in your simulation, add the following to your force field file:
    pythonenergy
        module example_energy
        rcut 10.0
    end

Required Functions:
    - prologue(boxlist): Called once at the start of the simulation
    - compute_total(boxlist): Compute total system energy
    - compute_diff(boxlist, displist): Compute energy difference for a move
    - update(boxlist): Called after accepted moves to update internal state

The boxlist contains dictionaries with box information:
    - 'boxtype': string - type of box ('cube', 'ortho', etc.)
    - 'temperature': float - temperature of the box
    - 'pressure': float - pressure of the box  
    - 'volume': float - volume of the box
    - 'boxdimensions': numpy array - box dimensions
    - 'chemicalpotential': numpy array - chemical potentials
    - 'atomtype': numpy array - atom types
    - 'raw_atoms': numpy array - atomic coordinates (3 x nAtoms)
    - 'moleculecount': numpy array - number of molecules per type
    - 'moltype': numpy array - molecule type per atom
    - 'molsubindx': numpy array - molecule sub-index per atom
    - 'energytable': numpy array - energy table

The displist contains dictionaries with displacement information:
    - 'type': str - 'displacement', 'addition', 'deletion', etc.
    - 'moltype': int - molecule type
    - 'molindex': int - molecule index
    - 'atomindex': int - atom index
    - 'x_new', 'y_new', 'z_new': float - new coordinates (if applicable)

Return Format:
    Both compute_total and compute_diff should return a dictionary containing:
    - 'energy': float - the computed energy (or energy difference)
    - 'accept': bool - True if configuration is valid (optional, default True)
    - 'energytable' or 'denergytable': numpy array - per-atom energies (optional)
"""

import numpy as np

# =============================================================================
# Global parameters and cached values
# =============================================================================

# Lennard-Jones parameters (example)
epsilon = 1.0  # Energy parameter
sigma = 1.0    # Length parameter
rcut = 10.0    # Cutoff distance
rcut_sq = rcut * rcut
rmin_sq = 0.5 * 0.5  # Minimum distance squared for overlap detection

# Cached values
_total_energy = 0.0
_energy_table = None


# =============================================================================
# Required Functions
# =============================================================================

def prologue(boxlist):
    """
    Called once at the start of the simulation.
    Use this to initialize any state needed by the energy function.
    
    Args:
        boxlist: List of box dictionaries containing simulation state
    """
    global epsilon, sigma, rcut, rcut_sq, _energy_table
    
    print("Python Energy: Initializing example Lennard-Jones energy")
    print(f"  Number of boxes: {len(boxlist)}")
    for i, box in enumerate(boxlist):
        print(f"  Box {i}: {box['boxtype']}, T={box['temperature']}")
        nmax = box['raw_atoms'].shape[1]
        print(f"  Maximum atoms: {nmax}")
    
    # Initialize energy table
    if len(boxlist) > 0:
        nmax = boxlist[0]['raw_atoms'].shape[1]
        _energy_table = np.zeros(nmax, dtype=np.float64)


def compute_total(boxlist):
    """
    Compute the total energy of the system.
    
    This function is called for full energy calculations (initialization,
    periodic recalculation, etc.)
    
    Args:
        boxlist: List of box dictionaries containing current simulation state
    
    Returns:
        dict: Dictionary with 'energy', 'accept', and optionally 'energytable'
    """
    global _total_energy, _energy_table
    
    # Get the first box (for single-box simulations)
    box = boxlist[0]
    atoms = box['raw_atoms']
    atom_types = box['atomtype']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    nmax = atoms.shape[1]
    
    # Reset energy table
    if _energy_table is None or len(_energy_table) != nmax:
        _energy_table = np.zeros(nmax, dtype=np.float64)
    else:
        _energy_table[:] = 0.0
    
    total_energy = 0.0
    accept = True
    
    # Double loop over atom pairs
    for i in range(nmax - 1):
        # Check if atom i is active
        mol_type_i = mol_types[i]
        mol_idx_i = mol_subindx[i]
        if mol_idx_i > nmol[mol_type_i - 1]:
            continue
            
        for j in range(i + 1, nmax):
            # Check if atom j is active
            mol_type_j = mol_types[j]
            mol_idx_j = mol_subindx[j]
            if mol_idx_j > nmol[mol_type_j - 1]:
                continue
            
            # Skip atoms in the same molecule (intramolecular handled elsewhere)
            # For this simple example, we skip if same molecule
            # You may want to modify this based on your needs
            
            # Compute distance
            rx = atoms[0, i] - atoms[0, j]
            ry = atoms[1, i] - atoms[1, j]
            rz = atoms[2, i] - atoms[2, j]
            
            # Apply periodic boundary conditions if needed
            # (For a cube box, this would be minimum image convention)
            box_dim = box.get('boxdimensions')
            if box_dim is not None and box['boxtype'] != 'nobox':
                lx = box_dim[1, 0] - box_dim[0, 0]
                ly = box_dim[1, 1] - box_dim[0, 1]
                lz = box_dim[1, 2] - box_dim[0, 2]
                rx = rx - lx * round(rx / lx)
                ry = ry - ly * round(ry / ly)
                rz = rz - lz * round(rz / lz)
            
            rsq = rx*rx + ry*ry + rz*rz
            
            # Check for overlap
            if rsq < rmin_sq:
                print(f"Overlap detected between atoms {i} and {j}")
                return {'energy': 0.0, 'accept': False}
            
            # Compute LJ energy if within cutoff
            if rsq < rcut_sq:
                sig_over_r_sq = (sigma * sigma) / rsq
                sig_over_r_6 = sig_over_r_sq * sig_over_r_sq * sig_over_r_sq
                lj = 4.0 * epsilon * sig_over_r_6 * (sig_over_r_6 - 1.0)
                
                total_energy += lj
                _energy_table[i] += lj
                _energy_table[j] += lj
    
    _total_energy = total_energy
    
    return {
        'energy': total_energy,
        'accept': accept,
        'energytable': _energy_table
    }


def compute_diff(boxlist, displist):
    """
    Compute the energy difference for a proposed move.
    
    This function is called during Monte Carlo moves to evaluate
    the energy change that would result if the move is accepted.
    
    Args:
        boxlist: List of box dictionaries (current state, not yet updated)
        displist: List of displacement dictionaries describing the proposed move
    
    Returns:
        dict: Dictionary with 'energy' (difference), 'accept', and optionally 'denergytable'
    """
    # Get the first box
    box = boxlist[0]
    atoms = box['raw_atoms']
    atom_types = box['atomtype']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    nmax = atoms.shape[1]
    delta_energy = 0.0
    accept = True
    
    # Create delta energy table
    denergy_table = np.zeros(nmax, dtype=np.float64)
    
    # Get box dimensions for PBC
    box_dim = box.get('boxdimensions')
    has_pbc = box_dim is not None and box['boxtype'] != 'nobox'
    if has_pbc:
        lx = box_dim[1, 0] - box_dim[0, 0]
        ly = box_dim[1, 1] - box_dim[0, 1]
        lz = box_dim[1, 2] - box_dim[0, 2]
    
    # Process each displacement
    for disp in displist:
        disp_type = disp.get('type', 'displacement')
        
        if disp_type == 'displacement':
            # Single atom displacement
            i = disp.get('atomindex', 0) - 1  # Convert to 0-indexed
            x_new = disp.get('x_new', atoms[0, i])
            y_new = disp.get('y_new', atoms[1, i])
            z_new = disp.get('z_new', atoms[2, i])
            
            # Compute energy with all other atoms
            for j in range(nmax):
                if j == i:
                    continue
                    
                # Check if atom j is active
                mol_type_j = mol_types[j]
                mol_idx_j = mol_subindx[j]
                if mol_idx_j > nmol[mol_type_j - 1]:
                    continue
                
                # New distance
                rx_new = x_new - atoms[0, j]
                ry_new = y_new - atoms[1, j]
                rz_new = z_new - atoms[2, j]
                
                if has_pbc:
                    rx_new = rx_new - lx * round(rx_new / lx)
                    ry_new = ry_new - ly * round(ry_new / ly)
                    rz_new = rz_new - lz * round(rz_new / lz)
                
                rsq_new = rx_new*rx_new + ry_new*ry_new + rz_new*rz_new
                
                # Check for overlap
                if rsq_new < rmin_sq:
                    return {'energy': 0.0, 'accept': False}
                
                # New energy contribution
                e_new = 0.0
                if rsq_new < rcut_sq:
                    sig_over_r_sq = (sigma * sigma) / rsq_new
                    sig_over_r_6 = sig_over_r_sq * sig_over_r_sq * sig_over_r_sq
                    e_new = 4.0 * epsilon * sig_over_r_6 * (sig_over_r_6 - 1.0)
                
                # Old distance
                rx_old = atoms[0, i] - atoms[0, j]
                ry_old = atoms[1, i] - atoms[1, j]
                rz_old = atoms[2, i] - atoms[2, j]
                
                if has_pbc:
                    rx_old = rx_old - lx * round(rx_old / lx)
                    ry_old = ry_old - ly * round(ry_old / ly)
                    rz_old = rz_old - lz * round(rz_old / lz)
                
                rsq_old = rx_old*rx_old + ry_old*ry_old + rz_old*rz_old
                
                # Old energy contribution
                e_old = 0.0
                if rsq_old < rcut_sq:
                    sig_over_r_sq = (sigma * sigma) / rsq_old
                    sig_over_r_6 = sig_over_r_sq * sig_over_r_sq * sig_over_r_sq
                    e_old = 4.0 * epsilon * sig_over_r_6 * (sig_over_r_6 - 1.0)
                
                # Delta energy
                de = e_new - e_old
                delta_energy += de
                denergy_table[i] += de
                denergy_table[j] += de
        
        elif disp_type == 'addition':
            # Add a molecule/atom
            i = disp.get('atomindex', 0) - 1
            x_new = disp.get('x_new', 0.0)
            y_new = disp.get('y_new', 0.0)
            z_new = disp.get('z_new', 0.0)
            
            for j in range(nmax):
                if j == i:
                    continue
                
                mol_type_j = mol_types[j]
                mol_idx_j = mol_subindx[j]
                if mol_idx_j > nmol[mol_type_j - 1]:
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
                    sig_over_r_sq = (sigma * sigma) / rsq
                    sig_over_r_6 = sig_over_r_sq * sig_over_r_sq * sig_over_r_sq
                    e = 4.0 * epsilon * sig_over_r_6 * (sig_over_r_6 - 1.0)
                    delta_energy += e
                    denergy_table[i] += e
                    denergy_table[j] += e
        
        elif disp_type == 'deletion':
            # Remove a molecule/atom - compute negative of current interactions
            i = disp.get('atomindex', 0) - 1
            
            for j in range(nmax):
                if j == i:
                    continue
                
                mol_type_j = mol_types[j]
                mol_idx_j = mol_subindx[j]
                if mol_idx_j > nmol[mol_type_j - 1]:
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
                    sig_over_r_sq = (sigma * sigma) / rsq
                    sig_over_r_6 = sig_over_r_sq * sig_over_r_sq * sig_over_r_sq
                    e = 4.0 * epsilon * sig_over_r_6 * (sig_over_r_6 - 1.0)
                    delta_energy -= e
                    denergy_table[i] -= e
                    denergy_table[j] -= e
    
    return {
        'energy': delta_energy,
        'accept': accept,
        'denergytable': denergy_table
    }


def update(boxlist):
    """
    Called after an accepted move to update internal state.
    
    NOTE: This function is NOT automatically called by the Fortran side.
    If you need to track state changes, you can:
    1. Implement state tracking within compute_diff() by detecting configuration changes
    2. Call update() explicitly from compute_diff() when appropriate
    3. Use compute_total() periodically to resync state
    
    This is useful for maintaining cached values or updating
    any internal bookkeeping after the system configuration
    has been updated.
    
    Args:
        boxlist: List of box dictionaries (now contains the updated state)
    """
    global _total_energy, _energy_table
    
    # For efficiency, you might want to incrementally update the cached
    # energy rather than recomputing. This is just a simple example
    # that recomputes the full energy.
    # In practice, you would track which atoms changed and update accordingly.
    
    # Option 1: Recompute (simple but slow)
    # result = compute_total(boxlist)
    # _total_energy = result['energy']
    
    # Option 2: Trust that compute_diff was correct and just update the total
    # This requires tracking the last delta energy - see example below
    pass


# =============================================================================
# Helper functions for more complex energy models
# =============================================================================

def set_lj_parameters(eps, sig, rc):
    """
    Set Lennard-Jones parameters.
    
    Args:
        eps: Energy parameter (epsilon)
        sig: Length parameter (sigma)
        rc: Cutoff distance
    """
    global epsilon, sigma, rcut, rcut_sq
    epsilon = eps
    sigma = sig
    rcut = rc
    rcut_sq = rc * rc


def get_pair_energy(rsq):
    """
    Compute pair energy for a given squared distance.
    
    Args:
        rsq: Squared distance between atoms
        
    Returns:
        float: Pair energy
    """
    if rsq >= rcut_sq:
        return 0.0
    
    sig_over_r_sq = (sigma * sigma) / rsq
    sig_over_r_6 = sig_over_r_sq * sig_over_r_sq * sig_over_r_sq
    return 4.0 * epsilon * sig_over_r_6 * (sig_over_r_6 - 1.0)


def minimum_image(rx, ry, rz, box):
    """
    Apply minimum image convention for periodic boundaries.
    
    Args:
        rx, ry, rz: Distance components
        box: Box dictionary
        
    Returns:
        tuple: (rx, ry, rz) after applying PBC
    """
    box_dim = box.get('boxdimensions')
    if box_dim is None or box['boxtype'] == 'nobox':
        return rx, ry, rz
    
    lx = box_dim[1, 0] - box_dim[0, 0]
    ly = box_dim[1, 1] - box_dim[0, 1]
    lz = box_dim[1, 2] - box_dim[0, 2]
    
    rx = rx - lx * round(rx / lx)
    ry = ry - ly * round(ry / ly)
    rz = rz - lz * round(rz / lz)
    
    return rx, ry, rz
