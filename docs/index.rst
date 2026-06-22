ClassyMC Documentation
======================

ClassyMC is a modular, object-oriented Monte Carlo molecular simulation package
written in Fortran 2003. It supports a wide range of force fields, move types,
sampling strategies, and analysis tools through a flexible, text-based input system.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   chapters/getting_started
   chapters/input_format
   chapters/forcefield_file

.. toctree::
   :maxdepth: 2
   :caption: Command Reference

   auto/set_commands
   auto/boxes/index
   auto/forcefields/index
   auto/moves/index
   auto/sampling/index
   auto/neighborlists/index
   auto/constraints/index
   auto/analysis/index
   auto/trajectories/index

.. toctree::
   :maxdepth: 2
   :caption: Intramolecular Potentials

   auto/bonds/index
   auto/angles/index
   auto/torsions/index

.. note::

   The **Command Reference** pages are auto-generated from the Fortran source code.
   Run ``python generate_docs.py`` from the ``docs/`` directory to regenerate them
   after updating the source.
