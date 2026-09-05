"""
Example Python Constraint Script for ClassyMC

This script demonstrates how to create custom Monte Carlo constraints
using Python. Constraints are checked before and during moves to enforce
specific conditions on the simulation.

To use this constraint in your simulation, add the following to your input file:
    create constraint
        python example_constraint
    end_create

Required Functions:
    - prologue(boxinfo): Initialize the constraint
    - check_initial(boxinfo): Check if initial configuration satisfies constraint
    - diff_check(boxinfo, displist): Check if a proposed move satisfies constraint
    - post_energy(boxinfo, displist, energy_diff): Check after energy calculation
    - update(boxinfo): Update internal state after accepted move
    - epilogue(boxinfo): Cleanup at end of simulation

The boxinfo dictionary is the standardized ClassyPyObj box dict containing:
    - 'boxtype': str - type of box ('cube', 'ortho', 'nobox')
    - 'thread_id': int - parallel thread id
    - 'box_id': int - box identifier
    - 'energy': float - total energy
    - 'temperature': float - temperature
    - 'pressure': float - pressure
    - 'volume': float - volume
    - 'boxdimensions': numpy array (2 x ndim) - box bounds [[lo, hi], ...]
    - 'chemicalpotential': numpy array - chemical potentials
    - 'atomtype': numpy array - atom types (1-indexed)
    - 'raw_atoms': numpy array (ndim x nMaxAtoms) - atomic coordinates
    - 'moleculecount': numpy array - molecules per type
    - 'energytable': numpy array - per-atom energies
    - 'moltype': numpy array - molecule type per atom
    - 'molsubindx': numpy array - molecule sub-index per atom
    - 'molindx': numpy array - molecule index per atom
    - 'neighlists': list of neighbor-list dicts (if present)

The displist contains dictionaries with displacement information:
    - 'type': str - 'displacement', 'addition', 'deletion', 'volchange', etc.
    - 'moltype': int - molecule type
    - 'molindex': int - molecule index  
    - 'atomindex': int - atom index
    - 'x_new', 'y_new', 'z_new': float - new coordinates (for displacement/addition)
"""

import numpy as np

# Configuration - these could be set from external parameters
MAX_RADIUS = 15.0  # Maximum allowed radius from center of box
MIN_SEPARATION = 2.0  # Minimum separation between atoms
REJECTION_COUNT = 0
CHECK_COUNT = 0


def prologue(boxinfo):
    """
    Initialize the constraint.
    
    This is called once at the start of the simulation after the
    constraint is constructed.
    
    Args:
        boxinfo: Dictionary with box information
    """
    global MAX_RADIUS, MIN_SEPARATION
    
    n_max_atoms = boxinfo['raw_atoms'].shape[1]
    print("Python Constraint: Initializing example constraint")
    print(f"  Box type: {boxinfo['boxtype']}")
    print(f"  Temperature: {boxinfo['temperature']}")
    print(f"  Max atoms: {n_max_atoms}")
    print(f"  Constraint: MAX_RADIUS = {MAX_RADIUS}")
    print(f"  Constraint: MIN_SEPARATION = {MIN_SEPARATION}")


def check_initial(boxinfo):
    """
    Check if the initial configuration satisfies the constraint.
    
    This is called once after the simulation is set up to verify
    that the starting configuration is valid.
    
    Args:
        boxinfo: Dictionary with box information
        
    Returns:
        bool: True if constraint is satisfied, False otherwise
    """
    global MAX_RADIUS, MIN_SEPARATION
    
    raw_atoms = boxinfo['raw_atoms']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    box_dims = boxinfo['boxdimensions']
    ndim = raw_atoms.shape[0]
    n_max_atoms = raw_atoms.shape[1]
    
    # Calculate box center
    center = np.array([(box_dims[0, i] + box_dims[1, i]) / 2 for i in range(ndim)])
    
    # Get active atom positions
    positions = []
    for i in range(n_max_atoms):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx <= nmol[mol_type - 1]:
            pos = np.array([raw_atoms[j, i] for j in range(ndim)])
            positions.append(pos)
    
    if len(positions) == 0:
        return True
    
    positions = np.array(positions)
    
    # Check radius constraint
    for pos in positions:
        dist_from_center = np.linalg.norm(pos[:len(center)] - center)
        if dist_from_center > MAX_RADIUS:
            print(f"  Initial check failed: atom at distance {dist_from_center:.2f} > {MAX_RADIUS}")
            return False
    
    # Check minimum separation constraint
    n = len(positions)
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist < MIN_SEPARATION:
                print(f"  Initial check failed: atoms {i} and {j} too close ({dist:.2f} < {MIN_SEPARATION})")
                return False
    
    print("  Initial configuration satisfies all constraints")
    return True


def diff_check(boxinfo, displist):
    """
    Check if a proposed move satisfies the constraint.
    
    This is called for each trial move before the energy is calculated.
    It should be fast since it's called frequently.
    
    Args:
        boxinfo: Dictionary with current box state
        displist: List of displacement dictionaries
        
    Returns:
        bool: True if constraint is satisfied, False otherwise
    """
    global MAX_RADIUS, CHECK_COUNT, REJECTION_COUNT
    
    CHECK_COUNT += 1
    
    box_dims = boxinfo['boxdimensions']
    ndim = boxinfo['raw_atoms'].shape[0]
    
    # Calculate box center
    center = np.array([(box_dims[0, i] + box_dims[1, i]) / 2 for i in range(ndim)])
    
    # Check each displacement
    for disp in displist:
        disp_type = disp.get('type', 'unknown')
        
        if disp_type == 'displacement' or disp_type == 'addition':
            # Check the new position
            new_pos = np.array([disp.get('x_new', 0), 
                               disp.get('y_new', 0), 
                               disp.get('z_new', 0)])
            
            dist_from_center = np.linalg.norm(new_pos[:len(center)] - center)
            if dist_from_center > MAX_RADIUS:
                REJECTION_COUNT += 1
                return False
    
    return True


def post_energy(boxinfo, displist, energy_diff):
    """
    Check constraint after energy calculation (optional).
    
    This is useful for constraints that depend on the energy,
    such as limiting the maximum energy change.
    
    Args:
        boxinfo: Dictionary with current box state
        displist: List of displacement dictionaries
        energy_diff: float - energy difference from the move
        
    Returns:
        bool: True if constraint is satisfied, False otherwise
    """
    # Example: reject moves that increase energy too much
    MAX_ENERGY_INCREASE = 1000.0  # kJ/mol or whatever units are used
    
    if energy_diff > MAX_ENERGY_INCREASE:
        return False
    
    return True


def update(boxinfo):
    """
    Update internal state after a move is accepted.
    
    This is called after each accepted move. Use this to update
    any internal bookkeeping.
    
    Args:
        boxinfo: Dictionary with updated box state
    """
    # No internal state to update in this example
    pass


def epilogue(boxinfo):
    """
    Cleanup at end of simulation.
    
    This is called once at the end of the simulation.
    
    Args:
        boxinfo: Dictionary with final box state
    """
    global CHECK_COUNT, REJECTION_COUNT
    
    print("Python Constraint: Epilogue")
    print(f"  Total constraint checks: {CHECK_COUNT}")
    print(f"  Total rejections: {REJECTION_COUNT}")
    if CHECK_COUNT > 0:
        print(f"  Rejection rate: {100.0 * REJECTION_COUNT / CHECK_COUNT:.2f}%")


# =============================================================================
# Additional utility functions for more complex constraints
# =============================================================================

def compute_center_of_mass(boxinfo, mol_type=None):
    """
    Compute center of mass for all atoms or a specific molecule type.
    
    Args:
        boxinfo: Dictionary with box information
        mol_type: Optional molecule type to filter by (1-indexed)
        
    Returns:
        numpy array: Center of mass position
    """
    raw_atoms = boxinfo['raw_atoms']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    ndim = raw_atoms.shape[0]
    n_max_atoms = raw_atoms.shape[1]
    
    total_mass = 0.0
    com = np.zeros(ndim)
    
    for i in range(n_max_atoms):
        mt = mol_types[i]
        mol_idx = mol_subindx[i]
        
        if mol_idx <= nmol[mt - 1]:
            if mol_type is None or mt == mol_type:
                pos = np.array([raw_atoms[j, i] for j in range(ndim)])
                mass = 1.0  # Assume unit mass, could use actual masses
                com += mass * pos
                total_mass += mass
    
    if total_mass > 0:
        com /= total_mass
    
    return com


def compute_radius_of_gyration(boxinfo, mol_type=None):
    """
    Compute radius of gyration for the system.
    
    Args:
        boxinfo: Dictionary with box information
        mol_type: Optional molecule type to filter by
        
    Returns:
        float: Radius of gyration
    """
    raw_atoms = boxinfo['raw_atoms']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    ndim = raw_atoms.shape[0]
    n_max_atoms = raw_atoms.shape[1]
    
    com = compute_center_of_mass(boxinfo, mol_type)
    
    total_mass = 0.0
    rg_sq = 0.0
    
    for i in range(n_max_atoms):
        mt = mol_types[i]
        mol_idx = mol_subindx[i]
        
        if mol_idx <= nmol[mt - 1]:
            if mol_type is None or mt == mol_type:
                pos = np.array([raw_atoms[j, i] for j in range(ndim)])
                mass = 1.0
                rg_sq += mass * np.sum((pos - com)**2)
                total_mass += mass
    
    if total_mass > 0:
        rg_sq /= total_mass
    
    return np.sqrt(rg_sq)
