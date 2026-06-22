Getting Started
===============

Requirements
------------

* A Fortran 2003-compatible compiler (``gfortran ≥ 9`` or Intel ``ifort``)
* GNU Make
* *(Optional)* MPI library for parallel runs
* *(Optional)* Python 3 + NumPy for embedded-Python features
* *(Optional)* LAMMPS or AENet libraries for extended force fields

Building
--------

The simplest build uses ``gfortran`` with no external libraries:

.. code-block:: bash

   cd /path/to/ClassyMC
   make

To enable optional features pass flags to ``make``:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Flag
     - Effect
   * - ``MPI=1``
     - Build with MPI parallelism.
   * - ``PYTHON=1``
     - Enable the embedded-Python interface.
   * - ``LAMMPS=1``
     - Enable the LAMMPS energy back-end.
   * - ``AENET=1``
     - Enable the AENet neural-network potential.

Example:

.. code-block:: bash

   make MPI=1 PYTHON=1

The compiled binary is placed at ``classyMC`` in the repository root.

Running a Simulation
--------------------

Pass the simulation input script as the first command-line argument:

.. code-block:: bash

   ./classyMC input.in

For MPI parallel runs:

.. code-block:: bash

   mpirun -np 4 ./classyMC input.in

The input file name is arbitrary; by convention it is given the ``.in`` extension.

Example Simulations
-------------------

Ready-to-run examples are provided in the ``Example/`` directory.  Each
subdirectory contains an input script and the necessary forcefield and
coordinate files.  To run the LJ bulk example:

.. code-block:: bash

   cd Example/LJBulk
   ../../classyMC in.ljfcc

Output Files
------------

ClassyMC writes several output files during and after a run:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Contents
   * - ``ScreenOutput_*.txt``
     - Per-cycle thermodynamic averages and acceptance statistics.
   * - ``Config_*.clssy``
     - Restart/checkpoint coordinate files (written at ``configfrequency``).
   * - Trajectory files (e.g. ``Traj.dump``)
     - Atomic coordinates in the requested format at the requested frequency.
   * - ``umbrella.dat`` *(umbrella sampling)*
     - Bias potential and histogram data.
