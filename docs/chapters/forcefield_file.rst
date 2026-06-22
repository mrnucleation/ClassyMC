Forcefield File (``.clFF``)
============================

The forcefield file is a separate plain-text file referenced from the simulation
script via ``forcefield "file.clFF"``.  It defines the chemical species and the
energy model.

Overall Structure
-----------------

The file is divided into several **keyword blocks**, each terminated by
``end_<keyword>``:

.. code-block:: none

   moleculetypes <N>

   atomdef
     ...
   end_atomdef

   molecule <M>
     ...
   end_molecule

   forcefieldtype
     <type_keyword>
   end_forcefieldtype

   forcefield <FF_index>
     <parameters>
   end_forcefield

   # Optional intramolecular sections:
   bonddef
     ...
   end_bonddef

   angledef
     ...
   end_angledef

   torsiondef
     ...
   end_torsiondef

``moleculetypes``
-----------------

Declares the total number of distinct molecule types::

   moleculetypes 2

``atomdef``
-----------

Maps integer atom-type indices to element labels and atomic masses::

   atomdef
     "O"   15.999
     "H"    1.008
     "Na"  22.990
   end_atomdef

Atom types are numbered **sequentially starting from 1** in the order they
appear.

``molecule``
------------

Defines the topology of one molecule type.  The block index corresponds to the
molecule type index::

   molecule 1
     RegrowthType  Simple
     nAtoms        3
     atoms
       1           # atom type for atom 1
       2           # atom type for atom 2
       2           # atom type for atom 3
     end_atoms
     bonds
       1  1  2
       2  1  3
     end_bonds
   end_molecule

``RegrowthType`` controls how molecules are rebuilt during CBMC-style moves.
Available types: ``Simple``, ``Ridgid``, ``LinearCBMC``, ``BranchCBMC``.

``forcefieldtype``
------------------

Selects the energy model (see :doc:`../auto/forcefields/index` for the full list)::

   forcefieldtype
     EP_LJ_Cut
   end_forcefieldtype

``forcefield``
--------------

Provides the parameters for the selected forcefield.  The block index selects
which forcefield object is being configured (matches ``modify box N energycalc M``).

The exact syntax of the parameter lines depends on the forcefield type.  For
most pair styles a ``rCut`` global cutoff is set first, followed by self-pair
lines and then optional cross-pair lines:

.. code-block:: none

   forcefield 1
     rCut 2.5
     # type  epsilon  sigma  rMin
     1       1.0      1.0    0.5
     # type1  type2  epsilon  sigma  rMin
     1        2      0.8      1.2    0.4
   end_forcefield

``bonddef`` / ``angledef`` / ``torsiondef``
-------------------------------------------

Intramolecular potential parameters are specified in separate blocks.  Each
block line corresponds to one bond, angle, or torsion defined in the molecular
topology:

.. code-block:: none

   bonddef
     ridgid 1.0          # rigid bond of length 1.0
     harmonic 100.0 1.5  # harmonic: k=100, r0=1.5
   end_bonddef

   angledef
     ridgid 109.47       # fixed angle (degrees)
     harmonic 50.0 109.47
   end_angledef

   torsiondef
     ridgid              # frozen dihedral
     trappe 0.0 355.0 0.0 -58.5   # TraPPE: c0 c1 c2 c3
   end_torsiondef

See :doc:`../auto/bonds/index`, :doc:`../auto/angles/index`, and
:doc:`../auto/torsions/index` for full parameter references.
