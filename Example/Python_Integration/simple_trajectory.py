"""
Simple VASP Trajectory Writer for ClassyMC

This module writes trajectory frames as VASP POSCAR files using ASE,
with incrementing file IDs for each frame.

Usage in input file:
    Create trajectory
        python 1 500 "POSCAR" simple_trajectory
    End_Create

This will create files: POSCAR_000001, POSCAR_000002, etc.
"""

import numpy as np
from ase import Atoms
from ase.io import write as ase_write

_base_filename = None
_frame_count = 0


def prologue(boxinfo):
    """Initialize trajectory output."""
    global _base_filename, _frame_count
    
    _base_filename = boxinfo['filename']
    _frame_count = 0
    
    print(f"Simple Trajectory: Will write VASP files with base name {_base_filename}")


def write_frame(boxinfo, cycle):
    """Write a single trajectory frame as a VASP POSCAR file."""
    global _base_filename, _frame_count

    _frame_count += 1
    
    raw_atoms = boxinfo['raw_atoms']
    atom_types = boxinfo['atomtype']
    mol_types = boxinfo['moltype']
    mol_subindx = boxinfo['molsubindx']
    nmol = boxinfo['moleculecount']
    atomsymbols = boxinfo['atomsymbols']
    box_dims = boxinfo['boxdimensions']
    ndim = boxinfo['ndimension']
    
    # Collect active atoms
    positions = []
    symbols = []
    
    for i in range(boxinfo['nmaxatoms']):
        mol_type = mol_types[i]
        mol_idx = mol_subindx[i]
        
        if mol_idx <= nmol[mol_type - 1]:
            atom_type = atom_types[i]
            symbol = atomsymbols[atom_type - 1]
            pos = [raw_atoms[j, i] for j in range(ndim)]
            # Pad to 3D
            while len(pos) < 3:
                pos.append(0.0)
            positions.append(pos)
            symbols.append(symbol)
    
    # Build cell matrix
    cell = np.zeros((3, 3))
    for i in range(min(ndim, 3)):
        cell[i, i] = box_dims[1, i] - box_dims[0, i]
    
    # Create ASE Atoms object
    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=True
    )
    
    # Write VASP POSCAR file with incrementing ID
    output_filename = f"{_base_filename}_{_frame_count:06d}"
    ase_write(output_filename, atoms, format='vasp')
    
    print(f"Simple Trajectory: Wrote {output_filename} (cycle {cycle})")


def epilogue(boxinfo):
    """Finalize trajectory output."""
    global _frame_count
    
    print(f"Simple Trajectory: Wrote {_frame_count} VASP POSCAR files")
