# ClassyMC

General Purpose Object Oriented Molecular Monte Carlo code. ClassyMC is designed to be used for Metropolis Monte Carlo and other similar sampling methods in molecular simulations. It was created to streamline the implementation of new methods in a way that was quick and modular.

## Requirements

- **Fortran Compiler**: Fortran-2003 standard compatible compiler (Intel Fortran or GNU Fortran)
- **MPI**: Required for parallel implementations
- **FFTW3**: Required for P3M electrostatics

### Optional Dependencies

- **Python 3**: For embedded Python support (custom moves, analysis, trajectories)
- **LAPACK/BLAS**: Required for AENet neural network potentials
- **LAMMPS**: For using LAMMPS as a force field backend
- **ASE**: Python package for trajectory analysis and visualization (optional, for Python trajectories)

## Installation

### Basic Build

```bash
make                    # Intel Fortran (default)
make COMPILER=gfortran  # GNU Fortran
```

### Build with Optional Packages

ClassyMC supports modular package selection via command-line flags:

```bash
# Enable individual packages
make PYTHON=1           # Embedded Python support
make AENET=1            # AENet neural network potentials
make LAMMPS=1           # LAMMPS force field backend

# Combine multiple packages
make PYTHON=1 AENET=1   # Python + AENet
make PYTHON=1 LAMMPS=1  # Python + LAMMPS
```

### Debug Build

```bash
make DEBUG=1            # Debug build with detailed checks
make DEBUG=1 PYTHON=1   # Debug build with Python
```

### Shared Library

Build ClassyMC as a shared library for use with external programs:

```bash
make lib                # Basic shared library
make lib PYTHON=1       # Shared library with Python support
```

### Compiler Selection

```bash
make COMPILER=ifort     # Intel Fortran (default)
make COMPILER=gfortran  # GNU Fortran
```

### Custom Library Paths

If external libraries are not in standard locations:

```bash
make LAMMPS=1 LAMMPS_LIB_PATH=/path/to/lammps/lib
make AENET=1 AENET_LIB_PATH=/path/to/aenet/lib
```

### Build Help

View all available build options:

```bash
make help
```

### Clean Build

```bash
make clean              # Remove all build artifacts
```

## Package Details

### Python Support (`PYTHON=1`)

Enables embedded Python interpreter for:
- **Custom Monte Carlo moves**: Write move algorithms in Python
- **Custom analysis functions**: Perform analysis calculations in Python
- **Custom trajectory output**: Write trajectories in any format, with ASE integration

Example input script usage:
```
create moves 1
    python 1.0 my_custom_move

create trajectory 1
    python 1 1000 traj.xyz my_trajectory_module
```

Python modules should be placed in the working directory or Python path.

### AENet Support (`AENET=1`)

Enables AENet neural network potentials for machine-learned force fields. Requires LAPACK and BLAS libraries.

### LAMMPS Support (`LAMMPS=1`)

Enables using LAMMPS as a force field backend. Requires LAMMPS compiled as a shared library (`liblammps.so`).

## Usage

Commands for use in the input script can be found on the GitHub repo's wiki page. Example run scripts can be found in the `/Example/` directory.

### Running a Simulation

```bash
./classyMC input_script.in
```

### Running with MPI

```bash
mpirun -np 4 ./classyMC input_script.in
```

## Directory Structure

```
ClassyMC/
├── src/           # Core Fortran source files
├── python/        # Python integration modules
├── Example/       # Example simulation setups
├── lib/           # External libraries
├── mods/          # Compiled module files
├── objects/       # Compiled object files
└── docs/          # Documentation
```

## License

See LICENSE file for details.

## Contributing

Contributions are welcome! Please submit pull requests or open issues on the GitHub repository.
