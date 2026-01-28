"""
Example Python Trajectory Script for ClassyMC

This script demonstrates how to write custom trajectory output using Python,
with optional integration with ASE (Atomic Simulation Environment).

To use this trajectory in your simulation, add the following to your input file:
    create trajectory 1
        python 1 1000 traj.xyz example_trajectory

The first argument is the box number (1), the second is the output frequency
(every 1000 cycles), the third is the output filename, and the fourth is
this Python module name (without .py extension).

Required Functions:
    - prologue(boxinfo): Called once at start to initialize
    - write_frame(boxinfo, cycle): Called each time to write a frame
    - epilogue(boxinfo): Called at end for cleanup

The boxinfo dictionary contains:
    - 'filename': str - output filename
    - 'boxnum': int - box number
    - 'outfreq': int - output frequency
    - 'boxtype': str - type of box ('cube', 'ortho', 'nobox')
    - 'temperature': float - temperature
    - 'pressure': float - pressure
    - 'volume': float - volume
    - 'boxdimensions': numpy array (2 x ndim) - box bounds [[lo, hi], ...]
    - 'atomtype': numpy array - atom types (1-indexed)
    - 'raw_atoms': numpy array (3 x nMaxAtoms) - atomic coordinates
    - 'moleculecount': numpy array - molecules per type
    - 'natoms': int - current number of atoms
    - 'nmaxatoms': int - maximum atoms (for GC ensembles)
    - 'moltype': numpy array - molecule type per atom
    - 'molsubindx': numpy array - molecule sub-index per atom
    - 'atomsymbols': list of str - atomic symbols by type index
    - 'ndimension': int - number of dimensions (usually 3)
    - 'cycle': int - current cycle (in write_frame only)
"""

import numpy as np

# Global state
_output_file = None
_frame_count = 0
_use_ase = False
_ase_trajectory = None


def prologue(boxinfo):
    """
    Called once at the start of the simulation to initialize output.
    
    This is where you should open files, initialize ASE trajectory writers,
    or set up any other output infrastructure.
    
    Args:
        boxinfo: Dictionary containing box information and settings
    """
    global _output_file, _frame_count, _use_ase, _ase_trajectory
    
    filename = boxinfo['filename']
    _frame_count = 0
    
    print(f"Python Trajectory: Initializing output to {filename}")
    print(f"  Box type: {boxinfo['boxtype']}")
    print(f"  Temperature: {boxinfo['temperature']}")
    print(f"  Max atoms: {boxinfo['nmaxatoms']}")
    print(f"  Dimensions: {boxinfo['ndimension']}")
    
    # Check if ASE is available
    try:
        from ase import Atoms
        from ase.io import write as ase_write
        from ase.io.trajectory import Trajectory
        _use_ase = True
        print("  ASE detected - using ASE for trajectory output")
        
        # Open ASE trajectory file (supports multiple formats based on extension)
        # Using .traj for ASE native format, or .xyz, .extxyz for extended XYZ
        if filename.endswith('.traj'):
            _ase_trajectory = Trajectory(filename, 'w')
        else:
            # For other formats, we'll write frame by frame
            _ase_trajectory = None
            _output_file = open(filename, 'w')
    except ImportError:
        _use_ase = False
        print("  ASE not available - using simple XYZ format")
        _output_file = open(filename, 'w')


def write_frame(boxinfo, cycle):
    """
    Called each time a frame should be written.
    
    This function extracts the current atomic configuration and writes it
    to the output file in the desired format.
    
    Args:
        boxinfo: Dictionary with current box state
        cycle: int - current simulation cycle number
    """
    global _output_file, _frame_count, _use_ase, _ase_trajectory
    
    # Extract data from boxinfo
    raw_atoms = boxinfo['raw_atoms']  # Shape: (3, nMaxAtoms)
    atom_types = boxinfo['atomtype']   # 1-indexed atom types
    mol_types = boxinfo['moltype']     # Molecule type per atom
    mol_subindx = boxinfo['molsubindx']  # Molecule sub-index
    nmol = boxinfo['moleculecount']    # Molecules per type
    natoms = boxinfo['natoms']
    atomsymbols = boxinfo['atomsymbols']  # List of symbols by type
    box_dims = boxinfo['boxdimensions']  # Shape: (2, ndim)
    ndim = boxinfo['ndimension']
    
    # Build list of active atoms (those that exist in the system)
    positions = []
    symbols = []
    
    for i in range(boxinfo['nmaxatoms']):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        
        # Check if this atom's molecule exists
        if mol_idx <= nmol[mol_type - 1]:  # -1 because nmol is 0-indexed
            atom_type = atom_types[i]
            symbol = atomsymbols[atom_type - 1]  # -1 because atomsymbols is 0-indexed
            positions.append([raw_atoms[j, i] for j in range(ndim)])
            symbols.append(symbol)
    
    # Pad to 3D if needed
    if ndim < 3:
        for pos in positions:
            while len(pos) < 3:
                pos.append(0.0)
    
    positions = np.array(positions)
    
    # Build cell from box dimensions
    cell = np.zeros((3, 3))
    for i in range(min(ndim, 3)):
        cell[i, i] = box_dims[1, i] - box_dims[0, i]  # hi - lo
    
    if _use_ase:
        _write_frame_ase(positions, symbols, cell, cycle, boxinfo)
    else:
        _write_frame_xyz(positions, symbols, cell, cycle, boxinfo)
    
    _frame_count += 1


def _write_frame_ase(positions, symbols, cell, cycle, boxinfo):
    """Write a frame using ASE."""
    global _ase_trajectory, _output_file
    
    from ase import Atoms
    from ase.io import write as ase_write
    
    # Determine periodicity based on box type
    boxtype = boxinfo['boxtype']
    pbc = [boxtype != 'nobox'] * 3
    
    # Create ASE Atoms object
    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc
    )
    
    # Add info about the simulation state
    atoms.info['cycle'] = cycle
    atoms.info['temperature'] = boxinfo['temperature']
    atoms.info['pressure'] = boxinfo['pressure']
    atoms.info['volume'] = boxinfo['volume']
    
    if _ase_trajectory is not None:
        # Native ASE trajectory format
        _ase_trajectory.write(atoms)
    else:
        # Write to file using ASE's format detection
        # Append mode for most formats
        filename = boxinfo['filename']
        if filename.endswith('.extxyz') or filename.endswith('.xyz'):
            # Extended XYZ supports append
            ase_write(_output_file, atoms, format='extxyz')
        else:
            # For other formats, write to the open file
            ase_write(_output_file, atoms)


def _write_frame_xyz(positions, symbols, cell, cycle, boxinfo):
    """Write a frame in simple XYZ format (no ASE required)."""
    global _output_file
    
    natoms = len(symbols)
    
    # Write header
    _output_file.write(f"{natoms}\n")
    
    # Comment line with cycle info
    comment = f"Cycle={cycle} T={boxinfo['temperature']:.2f} V={boxinfo['volume']:.4f}"
    if boxinfo['boxtype'] != 'nobox':
        # Add lattice info in extended XYZ format
        lattice = f'Lattice="{cell[0,0]:.6f} {cell[0,1]:.6f} {cell[0,2]:.6f} '
        lattice += f'{cell[1,0]:.6f} {cell[1,1]:.6f} {cell[1,2]:.6f} '
        lattice += f'{cell[2,0]:.6f} {cell[2,1]:.6f} {cell[2,2]:.6f}"'
        comment = f'{lattice} {comment}'
    _output_file.write(f"{comment}\n")
    
    # Write atom lines
    for i in range(natoms):
        _output_file.write(f"{symbols[i]:4s} {positions[i,0]:15.8f} {positions[i,1]:15.8f} {positions[i,2]:15.8f}\n")
    
    _output_file.flush()


def epilogue(boxinfo):
    """
    Called at the end of the simulation for cleanup.
    
    Close files, finalize trajectory writers, etc.
    
    Args:
        boxinfo: Dictionary with final box state
    """
    global _output_file, _frame_count, _ase_trajectory
    
    print(f"Python Trajectory: Closing output after {_frame_count} frames")
    
    if _ase_trajectory is not None:
        _ase_trajectory.close()
    
    if _output_file is not None:
        _output_file.close()


# =============================================================================
# Utility functions for users who want to extend this module
# =============================================================================

def get_ase_atoms(boxinfo):
    """
    Convert boxinfo to an ASE Atoms object.
    
    This is a utility function for users who want to use ASE functionality
    in their custom trajectory scripts.
    
    Args:
        boxinfo: Dictionary with box state
        
    Returns:
        ase.Atoms: Atoms object with current configuration
        
    Raises:
        ImportError: If ASE is not installed
    """
    from ase import Atoms
    
    # Extract data
    raw_atoms = boxinfo['raw_atoms']
    atom_types = boxinfo['atomtype']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    atomsymbols = boxinfo['atomsymbols']
    box_dims = boxinfo['boxdimensions']
    ndim = boxinfo['ndimension']
    
    # Build positions and symbols for active atoms
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
    boxtype = boxinfo['boxtype']
    pbc = [boxtype != 'nobox'] * 3
    
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)


def compute_rdf(boxinfo, r_max=10.0, n_bins=100, element_pairs=None):
    """
    Compute radial distribution function from current configuration.
    
    This is an example of analysis that can be done in a trajectory script.
    
    Args:
        boxinfo: Dictionary with box state
        r_max: Maximum distance for RDF
        n_bins: Number of histogram bins
        element_pairs: List of (elem1, elem2) tuples to compute RDF for,
                      or None for all pairs
    
    Returns:
        dict: {'r': array of distances, 'g(r)': dict of g(r) by pair}
    """
    try:
        atoms = get_ase_atoms(boxinfo)
        from ase.geometry.analysis import Analysis
        ana = Analysis(atoms)
        
        # Get all unique element pairs if not specified
        if element_pairs is None:
            unique_elements = list(set(atoms.get_chemical_symbols()))
            element_pairs = []
            for i, e1 in enumerate(unique_elements):
                for e2 in unique_elements[i:]:
                    element_pairs.append((e1, e2))
        
        result = {'r': None, 'g(r)': {}}
        
        for e1, e2 in element_pairs:
            rdf = ana.get_rdf(rmax=r_max, nbins=n_bins, elements=(e1, e2))
            if result['r'] is None:
                result['r'] = np.linspace(0, r_max, n_bins)
            result['g(r)'][(e1, e2)] = rdf[0]
        
        return result
    except Exception as e:
        print(f"RDF computation failed: {e}")
        return None
