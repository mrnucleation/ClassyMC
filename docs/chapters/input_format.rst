Simulation Input Script
=======================

The simulation is controlled by a plain-text input script passed on the command
line.  The parser is **case-insensitive** and ignores blank lines and lines
whose first non-whitespace character is ``#``.

Top-Level Commands
------------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose
   * - ``set <key> <value>``
     - Set a global simulation parameter (see :doc:`../auto/set_commands`).
   * - ``forcefield "file.clFF"``
     - Load an external forcefield definition file.
   * - ``create <object> ... end_create``
     - Allocate and configure a collection of objects.
   * - ``modify <object> <N> <key> <value>``
     - Set a parameter on a previously created object.
   * - ``neighlist <box> <type> <n>``
     - Create neighbor list(s) for box *N*.
   * - ``samplingtype <type>``
     - Choose the acceptance / sampling rule.
   * - ``run``
     - Initialize and run the Monte Carlo simulation.
   * - ``minimize <box>``
     - Run an energy minimization on box *N*.

``create`` Blocks
-----------------

A ``create`` block allocates one or more objects of the same kind.  The
general form is::

   create <object_type>
     <entry1>
     <entry2>
     ...
   end_create

where ``<object_type>`` is one of:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Object Type
     - What It Creates
   * - ``boxes``
     - One simulation box per line.  Each line specifies the box type keyword, or ``fromfile "filename"`` to load coordinates.
   * - ``moves``
     - One MC move per line: ``<type_keyword> <probability_weight>``.
   * - ``analysis``
     - One analysis object per line (type keyword plus arguments).
   * - ``trajectory``
     - One trajectory writer per line (type keyword plus arguments).
   * - ``constraint <N>``
     - One or more constraints for box *N*.

``modify`` Command
------------------

The ``modify`` command changes a parameter on an already-created object::

   modify <object_type> <index> <keyword> <value>

Examples:

.. code-block:: none

   modify box 1 temperature 1.2
   modify box 1 chempotential 1 -5.5
   modify move 2 radius 1.5
   modify sampling whamfreq 10000
   modify analysis 1 frequency 500

``modify sampling`` does not take an index number because only one sampling
object exists.

Annotated Example
-----------------

The following script runs a Lennard-Jones bulk simulation:

.. code-block:: none

   # ---- Global settings ----
   set rng_seed      -6
   set NeighSkin     6.0
   set moves         200
   set cycles        60000

   SamplingType Metropolis

   # ---- Chemistry ----
   ForceField "LJForcefield.clFF"

   # ---- Simulation box ----
   Create Boxes
     fromfile "LJBulk.clssy"
   End_Create

   modify box 1 temperature    0.7
   modify box 1 pressure       0.0001
   modify box 1 buildfreq      1
   modify box 1 energycalc     1
   modify box 1 chempotential  1  0.0

   NeighList 1 celllist 1

   # ---- Trial moves ----
   Create Moves
     MolTranslation  267.0
     IsoVol            1.0
   End_Create

   # ---- Trajectory output ----
   Create trajectory
     dump 1 1000 "Traj.dump"
   End_Create

   Run

Coordinate Files (``.clssy``)
------------------------------

Coordinate files are loaded with ``fromfile "filename.clssy"`` inside a
``create boxes`` block.  They use a simple keyword format:

.. code-block:: none

   boxtype ortho
   Lx 20.0
   Ly 20.0
   Lz 20.0

   molmin   0
   molmax 100
   mol     50

   # MolType  MolNumber  AtomNumber  x   y   z
   1  1  1   0.0  0.0  0.0
   1  1  2   0.0  0.0  1.5
   ...

MPI-parallel simulations can reference multiple files by using ``&`` as a
placeholder for the MPI rank index::

   fromfile "Config_&.clssy"
