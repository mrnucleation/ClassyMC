"""
Radius Constraint for ClassyMC

This constraint keeps all atoms within a maximum distance from the center
of the simulation box. This is useful for cluster simulations or for
confining particles to a specific region.

Usage in input file:
    Create constraint
        python radius_constraint
    end_create
"""

import numpy as np

# Configuration
MAX_RADIUS = 12.0  # Maximum allowed distance from box center

# Statistics
_check_count = 0
_reject_count = 0
_box_center = None


def prologue(boxinfo):
    """Initialize the constraint."""
    global _box_center, MAX_RADIUS
    
    box_dims = boxinfo['boxdimensions']
    ndim = boxinfo['ndimension']
    
    # Calculate and store box center
    _box_center = np.array([(box_dims[0, i] + box_dims[1, i]) / 2 
                            for i in range(ndim)])
    
    print(f"Radius Constraint: Initialized")
    print(f"  Box center: {_box_center}")
    print(f"  Max radius: {MAX_RADIUS}")


def check_initial(boxinfo):
    """Check if initial configuration satisfies the constraint."""
    global _box_center, MAX_RADIUS
    
    raw_atoms = boxinfo['raw_atoms']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    ndim = boxinfo['ndimension']
    
    max_dist = 0.0
    
    for i in range(boxinfo['nmaxatoms']):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        
        if mol_idx <= nmol[mol_type - 1]:
            pos = np.array([raw_atoms[j, i] for j in range(ndim)])
            dist = np.linalg.norm(pos - _box_center)
            max_dist = max(max_dist, dist)
            
            if dist > MAX_RADIUS:
                print(f"  Atom {i} at distance {dist:.2f} exceeds MAX_RADIUS {MAX_RADIUS}")
                return False
    
    print(f"  All atoms within radius. Max distance: {max_dist:.2f}")
    return True


def diff_check(boxinfo, displist):
    """Check if a proposed move satisfies the constraint."""
    global _box_center, MAX_RADIUS, _check_count, _reject_count
    
    _check_count += 1
    ndim = boxinfo['ndimension']
    
    for disp in displist:
        disp_type = disp.get('type', 'unknown')
        
        if disp_type in ('displacement', 'addition'):
            # Get new position
            new_pos = np.array([
                disp.get('x_new', 0),
                disp.get('y_new', 0),
                disp.get('z_new', 0)
            ])[:ndim]
            
            # Check distance from center
            dist = np.linalg.norm(new_pos - _box_center[:ndim])
            
            if dist > MAX_RADIUS:
                _reject_count += 1
                return False
    
    return True


def post_energy(boxinfo, displist, energy_diff):
    """Post-energy check (always accept in this constraint)."""
    return True


def update(boxinfo):
    """Update after accepted move."""
    pass


def epilogue(boxinfo):
    """Print final statistics."""
    global _check_count, _reject_count
    
    print(f"Radius Constraint: Completed")
    print(f"  Total checks: {_check_count}")
    print(f"  Rejections: {_reject_count}")
    if _check_count > 0:
        rate = 100.0 * _reject_count / _check_count
        print(f"  Rejection rate: {rate:.2f}%")
