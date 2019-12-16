# ClassyMC
General Purpose Object Oriented Molecular Monte Carlo and Genetic Algorithm code

--Installation

ClassyMC requires a Fortran-2003 standard compatible fortran compiler.  It has been tested on both Intel 
and GNU fortran based compilers.

For parallel implimentations MPI is required. If one wishes to install with the AENet forcefield one must compile AENet 
as a static library and import it at compile time.

To perform a basic compile specify the make command using your compiler of choice.
    make gfortran  # GNU Fortran Compiler    
    make ifort  #Intel Compiler
    
