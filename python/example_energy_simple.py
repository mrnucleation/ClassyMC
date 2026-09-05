"""
Simple Python Energy Script for ClassyMC - Harmonic Tether Example

This script demonstrates a minimal implementation of a Python energy function.
It implements a harmonic spring that tethers atoms to their initial positions,
which is useful for testing or creating Einstein crystal reference states.

To use this energy function in your simulation, add to your force field file:
    forcefieldtype
        pythonenergy
    end_forcefieldtype
    
    forcefield 1
        module example_energy_simple
    end_forcefield

Required Functions:
    - prologue(boxlist): Initialize with starting configuration
    - compute_total(boxlist): Compute total system energy
    - compute_diff(boxlist, displist): Compute energy difference for a move
    - compute_forces(boxlist, coords): Forces at coords (for native forcebias)

The return dictionary must contain:
    - 'energy': float - the computed energy (required)
    - 'accept': bool - True if valid configuration (optional, default True)
"""

import numpy as np

# =============================================================================
# Parameters
# =============================================================================

# Spring constant (energy/distance^2)
spring_constant = 10.0

# Reference positions (set during prologue)
_reference_positions = None
_initialized = False


# =============================================================================
# Required Functions
# =============================================================================

def prologue(boxlist):
    """
    Initialize the energy function.
    Store the initial positions as reference points for the harmonic tether.
    """
    global _reference_positions, _initialized, spring_constant
    
    print("Python Energy (Simple): Initializing harmonic tether potential")
    print(f"  Spring constant: {spring_constant}")
    
    if len(boxlist) == 0:
        print("  WARNING: No boxes found")
        return
    
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    nmax = atoms.shape[1]
    
    # Store a copy of the initial positions as reference
    _reference_positions = np.copy(atoms)
    
    # Count active atoms
    n_active = 0
    for i in range(nmax):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx <= nmol[mol_type - 1]:
            n_active += 1
    
    print(f"  Number of active atoms: {n_active}")
    print(f"  Reference positions stored")
    
    _initialized = True


def compute_total(boxlist):
    """
    Compute the total harmonic tether energy.
    
    E = sum_i (k/2) * |r_i - r_i^ref|^2
    
    Returns:
        dict with 'energy' and 'accept' keys
    """
    global _reference_positions
    
    if not _initialized or _reference_positions is None:
        return {'energy': 0.0, 'accept': True}
    
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    nmax = atoms.shape[1]
    total_energy = 0.0
    
    # Sum over all active atoms
    for i in range(nmax):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx > nmol[mol_type - 1]:
            continue
        
        # Distance from reference position
        dx = atoms[0, i] - _reference_positions[0, i]
        dy = atoms[1, i] - _reference_positions[1, i]
        dz = atoms[2, i] - _reference_positions[2, i]
        
        rsq = dx*dx + dy*dy + dz*dz
        
        # Harmonic energy: (k/2) * r^2
        total_energy += 0.5 * spring_constant * rsq
    
    return {
        'energy': total_energy,
        'accept': True
    }


def compute_diff(boxlist, displist):
    """
    Compute the energy difference for a proposed move.
    
    For displacement: E_new - E_old = (k/2) * (|r_new - r_ref|^2 - |r_old - r_ref|^2)
    
    Returns:
        dict with 'energy' (difference) and 'accept' keys
    """
    global _reference_positions
    
    if not _initialized or _reference_positions is None:
        return {'energy': 0.0, 'accept': True}
    
    box = boxlist[0]
    atoms = box['raw_atoms']
    
    delta_energy = 0.0
    
    for disp in displist:
        disp_type = disp.get('type', 'displacement')
        
        if disp_type == 'displacement':
            # Get atom index (convert from 1-indexed to 0-indexed)
            i = disp.get('atomindex', 0) - 1
            
            # New position
            x_new = disp.get('x_new', atoms[0, i])
            y_new = disp.get('y_new', atoms[1, i])
            z_new = disp.get('z_new', atoms[2, i])
            
            # Old position
            x_old = atoms[0, i]
            y_old = atoms[1, i]
            z_old = atoms[2, i]
            
            # Reference position
            x_ref = _reference_positions[0, i]
            y_ref = _reference_positions[1, i]
            z_ref = _reference_positions[2, i]
            
            # Old distance squared from reference
            dx_old = x_old - x_ref
            dy_old = y_old - y_ref
            dz_old = z_old - z_ref
            rsq_old = dx_old*dx_old + dy_old*dy_old + dz_old*dz_old
            
            # New distance squared from reference
            dx_new = x_new - x_ref
            dy_new = y_new - y_ref
            dz_new = z_new - z_ref
            rsq_new = dx_new*dx_new + dy_new*dy_new + dz_new*dz_new
            
            # Energy difference
            delta_energy += 0.5 * spring_constant * (rsq_new - rsq_old)
        
        elif disp_type == 'addition':
            # For addition, compute energy of new atom relative to its position
            # (which becomes its own reference)
            i = disp.get('atomindex', 0) - 1
            x_new = disp.get('x_new', 0.0)
            y_new = disp.get('y_new', 0.0)
            z_new = disp.get('z_new', 0.0)
            
            # New atom starts at its reference position, so energy = 0
            # (unless you want a different behavior)
            delta_energy += 0.0
        
        elif disp_type == 'deletion':
            # For deletion, subtract the energy of the removed atom
            i = disp.get('atomindex', 0) - 1
            
            x_old = atoms[0, i]
            y_old = atoms[1, i]
            z_old = atoms[2, i]
            
            x_ref = _reference_positions[0, i]
            y_ref = _reference_positions[1, i]
            z_ref = _reference_positions[2, i]
            
            dx = x_old - x_ref
            dy = y_old - y_ref
            dz = z_old - z_ref
            rsq = dx*dx + dy*dy + dz*dz
            
            delta_energy -= 0.5 * spring_constant * rsq

        elif disp_type == 'newstate_isomol':
            new_atoms = disp.get('new_atoms')
            if new_atoms is None:
                return {'energy': 0.0, 'accept': False}
            mol_types = box['moltype']
            mol_subindx = box['molsubindx']
            nmol = box['moleculecount']
            nmax = min(atoms.shape[1], new_atoms.shape[1], _reference_positions.shape[1])
            for i in range(nmax):
                mol_type = mol_types[i]
                mol_idx = mol_subindx[i]
                if mol_idx > nmol[mol_type - 1]:
                    continue
                dx_old = atoms[0, i] - _reference_positions[0, i]
                dy_old = atoms[1, i] - _reference_positions[1, i]
                dz_old = atoms[2, i] - _reference_positions[2, i]
                dx_new = new_atoms[0, i] - _reference_positions[0, i]
                dy_new = new_atoms[1, i] - _reference_positions[1, i]
                dz_new = new_atoms[2, i] - _reference_positions[2, i]
                rsq_old = dx_old*dx_old + dy_old*dy_old + dz_old*dz_old
                rsq_new = dx_new*dx_new + dy_new*dy_new + dz_new*dz_new
                delta_energy += 0.5 * spring_constant * (rsq_new - rsq_old)
    
    return {
        'energy': delta_energy,
        'accept': True
    }


def compute_forces(boxlist, coords):
    """
    Harmonic tether forces at coords: F = -k (r - r_ref).

    Required for native forcebias. `coords` may be a trial state.
    """
    global _reference_positions

    nmax = coords.shape[1]
    forces = np.zeros((3, nmax), dtype=np.float64, order='F')
    if not _initialized or _reference_positions is None:
        return {'forces': forces}

    box = boxlist[0]
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    ncopy = min(nmax, _reference_positions.shape[1])
    for i in range(ncopy):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx > nmol[mol_type - 1]:
            continue
        forces[0, i] = -spring_constant * (coords[0, i] - _reference_positions[0, i])
        forces[1, i] = -spring_constant * (coords[1, i] - _reference_positions[1, i])
        forces[2, i] = -spring_constant * (coords[2, i] - _reference_positions[2, i])
    return {'forces': forces}


def update(boxlist):
    """
    Optional: Update internal state after accepted move.
    For this simple potential, no update is needed.
    """
    pass


# =============================================================================
# Helper functions for customization
# =============================================================================

def set_spring_constant(k):
    """
    Set the spring constant.
    
    Args:
        k: Spring constant in energy/distance^2 units
    """
    global spring_constant
    spring_constant = k
    print(f"Python Energy (Simple): Spring constant set to {k}")


def reset_reference_positions(boxlist):
    """
    Reset reference positions to current configuration.
    Useful for equilibration or staged simulations.
    """
    global _reference_positions
    
    if len(boxlist) == 0:
        return
    
    box = boxlist[0]
    _reference_positions = np.copy(box['raw_atoms'])
    print("Python Energy (Simple): Reference positions reset")


def get_mean_displacement(boxlist):
    """
    Calculate the mean squared displacement from reference positions.
    Useful for monitoring.
    
    Returns:
        float: Mean squared displacement
    """
    global _reference_positions
    
    if _reference_positions is None or len(boxlist) == 0:
        return 0.0
    
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    nmax = atoms.shape[1]
    total_rsq = 0.0
    n_active = 0
    
    for i in range(nmax):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx > nmol[mol_type - 1]:
            continue
        
        dx = atoms[0, i] - _reference_positions[0, i]
        dy = atoms[1, i] - _reference_positions[1, i]
        dz = atoms[2, i] - _reference_positions[2, i]
        
        total_rsq += dx*dx + dy*dy + dz*dz
        n_active += 1
    
    if n_active == 0:
        return 0.0
    
    return total_rsq / n_active
