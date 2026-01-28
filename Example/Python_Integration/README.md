# Python Integration Examples

This directory contains examples demonstrating ClassyMC's Python integration features.

## Requirements

Build ClassyMC with Python support:
```bash
make PYTHON=1
```

Optional: Install ASE for enhanced trajectory support:
```bash
pip install ase
```

## Examples

### 1. Python Constraint (`in.python_constraint`)

Demonstrates using Python scripts to define custom constraints.

```bash
./classyMC in.python_constraint
```

The `radius_constraint.py` module keeps all particles within a specified
distance from the box center.

**Key Features:**
- Custom geometry constraints
- Real-time rejection statistics
- Configurable parameters in Python

### 2. Python Trajectory (`in.python_trajectory`)

Demonstrates using Python scripts for trajectory output.

```bash
./classyMC in.python_trajectory
```

The `simple_trajectory.py` module writes extended XYZ format files
compatible with visualization tools.

**Key Features:**
- Custom output formats
- Integration with ASE
- Real-time analysis during trajectory writing

## File Descriptions

| File | Description |
|------|-------------|
| `LJ.clFF` | Simple Lennard-Jones forcefield |
| `Config.clssy` | Initial configuration (14 LJ particles) |
| `radius_constraint.py` | Python constraint: confine to sphere |
| `simple_trajectory.py` | Python trajectory: XYZ format output |
| `in.python_constraint` | Input script using Python constraint |
| `in.python_trajectory` | Input script using Python trajectory |

## Writing Your Own Python Modules

### Constraint Module Template

```python
def prologue(boxinfo):
    """Initialize the constraint."""
    pass

def check_initial(boxinfo):
    """Check initial configuration. Return True/False."""
    return True

def diff_check(boxinfo, displist):
    """Check proposed move. Return True/False."""
    return True

def post_energy(boxinfo, displist, energy_diff):
    """Check after energy calculation. Return True/False."""
    return True

def update(boxinfo):
    """Update after accepted move."""
    pass

def epilogue(boxinfo):
    """Cleanup at end of simulation."""
    pass
```

### Trajectory Module Template

```python
def prologue(boxinfo):
    """Initialize trajectory output."""
    pass

def write_frame(boxinfo, cycle):
    """Write a trajectory frame."""
    pass

def epilogue(boxinfo):
    """Close files and cleanup."""
    pass
```

### The `boxinfo` Dictionary

Both modules receive a `boxinfo` dictionary containing:

| Key | Type | Description |
|-----|------|-------------|
| `boxtype` | str | Box type ('cube', 'ortho', 'nobox') |
| `temperature` | float | Temperature |
| `pressure` | float | Pressure |
| `volume` | float | Volume |
| `boxdimensions` | ndarray | Box bounds (2 x ndim) |
| `atomtype` | ndarray | Atom types (1-indexed) |
| `raw_atoms` | ndarray | Coordinates (3 x nMaxAtoms) |
| `moleculecount` | ndarray | Molecules per type |
| `natoms` | int | Current number of atoms |
| `nmaxatoms` | int | Maximum atoms |
| `moltype` | ndarray | Molecule type per atom |
| `molsubindx` | ndarray | Molecule sub-index |
| `atomsymbols` | list | Atomic symbols by type |
| `ndimension` | int | Number of dimensions |

### The `displist` for Constraints

The `diff_check` function receives a list of displacement dictionaries:

| Key | Type | Description |
|-----|------|-------------|
| `type` | str | 'displacement', 'addition', 'deletion', etc. |
| `moltype` | int | Molecule type |
| `molindex` | int | Molecule index |
| `atomindex` | int | Atom index |
| `x_new`, `y_new`, `z_new` | float | New coordinates (if applicable) |
