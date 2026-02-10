# APeriodic Orthorhombic Box with Python Energy Example

This example demonstrates the use of ClassyMC's `aperiodic_ortho` box type with a Python-based force field for a two-component Lennard-Jones system.

## Features

### APeriodic Orthorhombic Box
- **Box type**: `aperiodic_ortho`
- **Dimensions**: 30.0 × 30.0 × 40.0 Å
- **Boundary conditions**: 
  - X: Periodic
  - Y: Periodic  
  - Z: Non-periodic (fixed/open)
- **Geometry**: Slab configuration ideal for interfaces, surfaces, or confined systems

### System Composition
- **64 total atoms**:
  - 32 Type 1 atoms (Ar-like): Larger atoms with stronger interactions
  - 32 Type 2 atoms (Ne-like): Smaller atoms with weaker interactions
- Initial configuration arranged in layers along the Z-direction

### Python Force Field
- Two-component Lennard-Jones potential
- Lorentz-Berthelot mixing rules for cross-interactions
- Parameters (in reduced units):
  - Type 1 (Ar): ε₁ = 1.0, σ₁ = 1.0
  - Type 2 (Ne): ε₂ = 0.5, σ₂ = 0.8
  - Mixed: σ₁₂ = 0.9, ε₁₂ = 0.707
- Cutoff distance: 6.0 σ

### Python Energy Module
The `two_component_lj.py` module:
- Implements the full LJ potential with type-specific parameters
- Correctly handles mixed periodic boundary conditions (periodic in X,Y only)
- Includes overlap detection and statistics tracking
- Uses efficient neighbor list calculations

## Files

- `in.aperiodic_python` - Main simulation input script
- `APeriodic_FF.clFF` - Force field definition (links to Python module)
- `APeriodic_Config.clssy` - Initial configuration with box type and coordinates
- `two_component_lj.py` - Python energy module implementing the LJ potential

## Running the Simulation

### Prerequisites
ClassyMC must be compiled with Python support:
```bash
make PYTHON=1
```

### Execution
```bash
./classyMC in.aperiodic_python
```

### Output
- Console output with energy and acceptance statistics
- `Traj.dump` - Trajectory file written every 200 cycles in LAMMPS dump format

## Simulation Parameters

- **Temperature**: 1.2 (reduced units)
- **Pressure**: 0.001 (reduced units)
- **RNG Seed**: 54321
- **Cycles**: 2000
- **Moves per cycle**: 64
- **Move type**: Molecular translation (100% weight)
- **Neighbor list**: Cell-based r² list with 6.5 σ cutoff, rebuilt every 10 cycles

## Expected Behavior

Due to the non-periodic Z boundary:
- Atoms will not wrap around in the Z direction
- The system forms a slab with free surfaces at Z=0 and Z=40
- Surface effects and stratification may develop
- The two atom types may segregate or mix depending on temperature

This geometry is particularly useful for studying:
- Liquid-vapor interfaces
- Confined fluids
- Surface adsorption
- Wetting phenomena
- Density profiles in inhomogeneous systems

## Modifying the Example

### Change boundary conditions
Edit `APeriodic_Config.clssy`:
```
boundary   p   f   f     # Only X periodic (rod-like geometry)
boundary   f   f   f     # No periodic boundaries (fully confined)
```

### Adjust LJ parameters
Edit `two_component_lj.py`:
```python
epsilon_1 = 1.5   # Stronger Ar interactions
sigma_2 = 0.6     # Smaller Ne atoms
```

### Change box dimensions
Edit `APeriodic_Config.clssy`:
```
dimension   40.0  40.0  60.0    # Larger, thicker slab
```

## Notes

- The Python energy module checks the box type and applies periodic wrapping only to periodic dimensions
- The initial configuration has atoms arranged in layers, not randomly distributed
- For production runs, increase the number of cycles and equilibration time
- Monitor energy drift and acceptance ratios to ensure proper sampling
