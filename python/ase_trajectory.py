"""
ASE-Focused Trajectory Module for ClassyMC

This module provides tight integration with ASE (Atomic Simulation Environment)
for trajectory output and analysis. It supports multiple output formats and
can perform on-the-fly analysis during the simulation.

Usage in input file:
    create trajectory 1
        python 1 1000 traj.traj ase_trajectory

Supported output formats (determined by filename extension):
    - .traj: ASE native binary format (fast, compact)
    - .xyz, .extxyz: Extended XYZ format (readable, portable)
    - .pdb: Protein Data Bank format
    - .cif: Crystallographic Information File
    - .vasp, .poscar: VASP POSCAR format
    - .json: ASE JSON format

Features:
    - Automatic format detection from filename
    - On-the-fly property calculation (energy, stress, etc.)
    - Built-in RDF and MSD analysis
    - Support for visualization with ASE GUI

Requirements:
    - ASE (pip install ase)
    - NumPy

Author: ClassyMC Python Interface
"""

import numpy as np

# Try to import ASE
try:
    from ase import Atoms
    from ase.io import write as ase_write
    from ase.io.trajectory import Trajectory
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    print("WARNING: ASE not found. Install with: pip install ase")


# Global state
_trajectory = None
_format = None
_filename = None
_frame_count = 0
_analysis_data = {}
_calc_properties = False


def prologue(boxinfo):
    """
    Initialize ASE trajectory output.
    
    The output format is automatically detected from the filename extension.
    """
    global _trajectory, _format, _filename, _frame_count, _analysis_data
    
    if not HAS_ASE:
        raise ImportError("ASE is required for this trajectory module. Install with: pip install ase")
    
    _filename = boxinfo['filename']
    _frame_count = 0
    _analysis_data = {'msd_reference': None, 'rdf_history': []}
    
    # Determine format from extension
    if _filename.endswith('.traj'):
        _format = 'traj'
        _trajectory = Trajectory(_filename, 'w')
    elif _filename.endswith('.extxyz') or _filename.endswith('.xyz'):
        _format = 'extxyz'
        _trajectory = open(_filename, 'w')
    elif _filename.endswith('.pdb'):
        _format = 'pdb'
        _trajectory = open(_filename, 'w')
    elif _filename.endswith('.cif'):
        _format = 'cif'
        _trajectory = None  # Will write each frame separately
    elif _filename.endswith('.vasp') or _filename.endswith('.poscar'):
        _format = 'vasp'
        _trajectory = None
    elif _filename.endswith('.json'):
        _format = 'json'
        _trajectory = None
    else:
        # Default to extended XYZ
        _format = 'extxyz'
        _trajectory = open(_filename, 'w')
    
    print(f"ASE Trajectory: Output format: {_format}")
    print(f"ASE Trajectory: Output file: {_filename}")


def write_frame(boxinfo, cycle):
    """
    Write a trajectory frame using ASE.
    """
    global _trajectory, _format, _filename, _frame_count, _analysis_data
    
    # Convert to ASE Atoms
    atoms = _boxinfo_to_atoms(boxinfo, cycle)
    
    # Write based on format
    if _format == 'traj':
        _trajectory.write(atoms)
    elif _format == 'extxyz':
        ase_write(_trajectory, atoms, format='extxyz')
    elif _format == 'pdb':
        ase_write(_trajectory, atoms, format='proteindatabank')
    elif _format == 'cif':
        frame_file = _filename.replace('.cif', f'_{_frame_count:06d}.cif')
        ase_write(frame_file, atoms, format='cif')
    elif _format == 'vasp':
        frame_file = _filename.replace('.vasp', f'_{_frame_count:06d}.vasp')
        frame_file = frame_file.replace('.poscar', f'_{_frame_count:06d}.poscar')
        ase_write(frame_file, atoms, format='vasp')
    elif _format == 'json':
        frame_file = _filename.replace('.json', f'_{_frame_count:06d}.json')
        ase_write(frame_file, atoms, format='json')
    
    # Store reference for MSD calculation
    if _frame_count == 0:
        _analysis_data['msd_reference'] = atoms.get_positions().copy()
    
    _frame_count += 1
    
    # Flush periodically
    if hasattr(_trajectory, 'flush') and _frame_count % 10 == 0:
        _trajectory.flush()


def epilogue(boxinfo):
    """
    Finalize trajectory output.
    """
    global _trajectory, _frame_count
    
    print(f"ASE Trajectory: Wrote {_frame_count} frames")
    
    if _trajectory is not None:
        if hasattr(_trajectory, 'close'):
            _trajectory.close()


def _boxinfo_to_atoms(boxinfo, cycle):
    """
    Convert boxinfo dictionary to ASE Atoms object.
    """
    # Extract data
    raw_atoms = boxinfo['raw_atoms']
    atom_types = boxinfo['atomtype']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    atomsymbols = boxinfo['atomsymbols']
    box_dims = boxinfo['boxdimensions']
    ndim = boxinfo['ndimension']
    
    # Build positions and symbols
    positions = []
    symbols = []
    
    for i in range(boxinfo['nmaxatoms']):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        
        if mol_idx <= nmol[mol_type - 1]:
            atom_type = atom_types[i]
            symbol = atomsymbols[atom_type - 1]
            pos = [raw_atoms[j, i] for j in range(ndim)]
            while len(pos) < 3:
                pos.append(0.0)
            positions.append(pos)
            symbols.append(symbol)
    
    # Build cell
    cell = np.zeros((3, 3))
    for i in range(min(ndim, 3)):
        cell[i, i] = box_dims[1, i] - box_dims[0, i]
    
    # Periodicity
    pbc = [boxinfo['boxtype'] != 'nobox'] * 3
    
    # Create Atoms
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)
    
    # Store simulation info
    atoms.info['cycle'] = cycle
    atoms.info['frame'] = _frame_count
    atoms.info['temperature'] = boxinfo['temperature']
    atoms.info['pressure'] = boxinfo['pressure']
    atoms.info['volume'] = boxinfo['volume']
    
    return atoms


# =============================================================================
# Analysis utilities that can be called from custom trajectory scripts
# =============================================================================

def compute_msd(boxinfo, reference_positions=None):
    """
    Compute Mean Square Displacement from current frame.
    
    Args:
        boxinfo: Box information dictionary
        reference_positions: Reference positions array, or None to use first frame
        
    Returns:
        float: MSD value
    """
    global _analysis_data
    
    atoms = _boxinfo_to_atoms(boxinfo, 0)
    current_pos = atoms.get_positions()
    
    if reference_positions is None:
        reference_positions = _analysis_data.get('msd_reference')
    
    if reference_positions is None:
        return 0.0
    
    # Handle different number of atoms (GC ensemble)
    n = min(len(current_pos), len(reference_positions))
    if n == 0:
        return 0.0
    
    displacement = current_pos[:n] - reference_positions[:n]
    msd = np.mean(np.sum(displacement**2, axis=1))
    
    return msd


def compute_coordination_number(boxinfo, cutoff=3.0, center_element=None, neighbor_element=None):
    """
    Compute average coordination number.
    
    Args:
        boxinfo: Box information dictionary
        cutoff: Cutoff distance for neighbors
        center_element: Element symbol for center atoms (or None for all)
        neighbor_element: Element symbol for neighbor atoms (or None for all)
        
    Returns:
        float: Average coordination number
    """
    atoms = _boxinfo_to_atoms(boxinfo, 0)
    
    from ase.neighborlist import neighbor_list
    
    # Get neighbor list
    i_list, j_list, d_list = neighbor_list('ijd', atoms, cutoff)
    
    symbols = atoms.get_chemical_symbols()
    
    # Filter by element if specified
    coord_counts = {}
    for i, j, d in zip(i_list, j_list, d_list):
        if center_element and symbols[i] != center_element:
            continue
        if neighbor_element and symbols[j] != neighbor_element:
            continue
        
        if i not in coord_counts:
            coord_counts[i] = 0
        coord_counts[i] += 1
    
    if not coord_counts:
        return 0.0
    
    return np.mean(list(coord_counts.values()))


def compute_rdf(boxinfo, r_max=10.0, n_bins=100, elements=None):
    """
    Compute radial distribution function.
    
    Args:
        boxinfo: Box information dictionary
        r_max: Maximum distance
        n_bins: Number of histogram bins
        elements: Tuple of (elem1, elem2) or None for total RDF
        
    Returns:
        tuple: (r_values, g_r)
    """
    atoms = _boxinfo_to_atoms(boxinfo, 0)
    
    try:
        from ase.geometry.analysis import Analysis
        ana = Analysis(atoms)
        
        if elements:
            rdf = ana.get_rdf(rmax=r_max, nbins=n_bins, elements=elements)
        else:
            rdf = ana.get_rdf(rmax=r_max, nbins=n_bins)
        
        r = np.linspace(0, r_max, n_bins)
        return r, rdf[0]
    except Exception as e:
        print(f"RDF calculation failed: {e}")
        return None, None


def view_current_frame(boxinfo):
    """
    Open ASE GUI to view the current configuration.
    
    This is useful for debugging or interactive analysis.
    """
    atoms = _boxinfo_to_atoms(boxinfo, 0)
    from ase.visualize import view
    view(atoms)


def get_atoms(boxinfo):
    """
    Get an ASE Atoms object from boxinfo.
    
    This is the main interface for users who want to do custom analysis.
    
    Args:
        boxinfo: Box information dictionary
        
    Returns:
        ase.Atoms: Atoms object
    """
    return _boxinfo_to_atoms(boxinfo, 0)
