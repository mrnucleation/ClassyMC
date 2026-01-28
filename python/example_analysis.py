"""
Example Python Analysis Script for ClassyMC

This script demonstrates how to create custom analysis functions
using Python. Analysis functions compute properties of the system
that can be used in umbrella sampling or for data collection.

To use this analysis in your simulation, add the following to your input file:
    create analysis
        python example_analysis
    end_create

Required Functions:
    - compute(boxlist): Compute the analysis value for current configuration
    - compute_new(boxlist, displist): Compute value for proposed configuration
    - update(boxlist): Called after accepted move to update internal state

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
    - 'energytable': numpy array - energy table

The displist contains dictionaries with displacement information:
    - 'type': str - 'displacement', 'addition', 'deletion', etc.
    - 'moltype': int - molecule type
    - 'molindex': int - molecule index
    - 'atomindex': int - atom index
    - 'x_new', 'y_new', 'z_new': float - new coordinates (if applicable)
"""

import numpy as np

# Global state for caching
_last_value = 0.0
_new_value = 0.0


def compute(boxlist):
    """
    Compute the analysis value for the current configuration.
    
    This example computes the radius of gyration of the system.
    
    Args:
        boxlist: List of box dictionaries containing simulation state
        
    Returns:
        float: The computed analysis value (radius of gyration)
    """
    global _last_value
    
    # Get the first box
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    # Collect active atom positions
    positions = []
    for i in range(atoms.shape[1]):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx <= nmol[mol_type - 1]:
            pos = [atoms[j, i] for j in range(3)]
            positions.append(pos)
    
    if len(positions) == 0:
        _last_value = 0.0
        return 0.0
    
    positions = np.array(positions)
    
    # Compute center of mass
    com = np.mean(positions, axis=0)
    
    # Compute radius of gyration
    rg_sq = np.mean(np.sum((positions - com)**2, axis=1))
    rg = np.sqrt(rg_sq)
    
    _last_value = rg
    return rg


def compute_new(boxlist, displist):
    """
    Compute the analysis value for a proposed new configuration.
    
    This function is called during a Monte Carlo move to evaluate
    the analysis value that would result if the move is accepted.
    
    Args:
        boxlist: List of box dictionaries (current state, not yet updated)
        displist: List of displacement dictionaries describing the proposed move
        
    Returns:
        float: The computed analysis value for the proposed new state
    """
    global _new_value
    
    # Get the first box
    box = boxlist[0]
    atoms = box['raw_atoms']
    mol_types = box['moltype']
    mol_subindx = box['molsubindx']
    nmol = box['moleculecount']
    
    # Collect active atom positions, applying displacements
    positions = []
    
    # Build a map of atom index -> new position from displist
    new_positions = {}
    for disp in displist:
        if disp.get('type') == 'displacement':
            atom_idx = disp.get('atomindex', 0) - 1  # Convert to 0-indexed
            new_positions[atom_idx] = [
                disp.get('x_new', 0),
                disp.get('y_new', 0),
                disp.get('z_new', 0)
            ]
    
    for i in range(atoms.shape[1]):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        if mol_idx <= nmol[mol_type - 1]:
            # Use new position if this atom was displaced
            if i in new_positions:
                pos = new_positions[i]
            else:
                pos = [atoms[j, i] for j in range(3)]
            positions.append(pos)
    
    if len(positions) == 0:
        _new_value = 0.0
        return 0.0
    
    positions = np.array(positions)
    
    # Compute center of mass
    com = np.mean(positions, axis=0)
    
    # Compute radius of gyration
    rg_sq = np.mean(np.sum((positions - com)**2, axis=1))
    rg = np.sqrt(rg_sq)
    
    _new_value = rg
    return rg


def update(boxlist):
    """
    Called after an accepted move to update internal state.
    
    This is useful for maintaining cached values or updating
    any internal bookkeeping after the system configuration
    has been updated.
    
    Args:
        boxlist: List of box dictionaries (now contains the updated state)
    """
    global _last_value, _new_value
    
    # The new value becomes the current value
    _last_value = _new_value


# =============================================================================
# Additional example analysis functions
# =============================================================================

def compute_energy_per_atom(boxlist):
    """
    Example: Compute energy per atom.
    
    Note: This requires the energy table to be populated.
    """
    box = boxlist[0]
    energy_table = box.get('energytable')
    natoms = box.get('natoms', 0)
    
    if energy_table is None or natoms == 0:
        return 0.0
    
    total_energy = np.sum(energy_table)
    return total_energy / natoms


def compute_density(boxlist):
    """
    Example: Compute number density (atoms per unit volume).
    """
    box = boxlist[0]
    volume = box.get('volume', 1.0)
    natoms = box.get('natoms', 0)
    
    if volume <= 0:
        return 0.0
    
    return natoms / volume


def compute_pair_distance(boxlist, atom1=0, atom2=1):
    """
    Example: Compute distance between two specific atoms.
    
    Useful for umbrella sampling on pair distances.
    """
    box = boxlist[0]
    atoms = box['raw_atoms']
    
    pos1 = np.array([atoms[j, atom1] for j in range(3)])
    pos2 = np.array([atoms[j, atom2] for j in range(3)])
    
    return np.linalg.norm(pos1 - pos2)
