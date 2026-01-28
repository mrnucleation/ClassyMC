"""
Simple XYZ Trajectory Writer for ClassyMC

This module writes trajectory frames in extended XYZ format,
which is widely compatible with visualization tools like VMD,
OVITO, and ASE.

Usage in input file:
    Create trajectory
        python 1 500 "traj.xyz" simple_trajectory
    End_Create
"""

import numpy as np


_output_file = None
_frame_count = 0


def prologue(boxinfo):
    """Initialize trajectory output."""
    global _output_file, _frame_count
    
    filename = boxinfo['filename']
    _frame_count = 0
    
    print(f"Simple Trajectory: Opening {filename}")
    _output_file = open(filename, 'w')


def write_frame(boxinfo, cycle):
    """Write a single trajectory frame."""
    global _output_file, _frame_count

    print("Hi from Python")
    
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
    
    natoms = len(positions)
    
    # Build cell string for extended XYZ
    cell = np.zeros((3, 3))
    for i in range(min(ndim, 3)):
        cell[i, i] = box_dims[1, i] - box_dims[0, i]
    
    lattice_str = f'{cell[0,0]:.6f} {cell[0,1]:.6f} {cell[0,2]:.6f} '
    lattice_str += f'{cell[1,0]:.6f} {cell[1,1]:.6f} {cell[1,2]:.6f} '
    lattice_str += f'{cell[2,0]:.6f} {cell[2,1]:.6f} {cell[2,2]:.6f}'
    
    # Write frame
    _output_file.write(f"{natoms}\n")
    _output_file.write(f'Lattice="{lattice_str}" Cycle={cycle} ')
    _output_file.write(f'Properties=species:S:1:pos:R:3\n')
    
    for i in range(natoms):
        _output_file.write(f"{symbols[i]:4s} ")
        _output_file.write(f"{positions[i][0]:15.8f} ")
        _output_file.write(f"{positions[i][1]:15.8f} ")
        _output_file.write(f"{positions[i][2]:15.8f}\n")
    
    _output_file.flush()
    _frame_count += 1


def epilogue(boxinfo):
    """Close trajectory file."""
    global _output_file, _frame_count
    
    print(f"Simple Trajectory: Wrote {_frame_count} frames")
    
    if _output_file is not None:
        _output_file.close()
