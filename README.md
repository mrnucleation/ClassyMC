# ClassyMC
General Purpose Object Oriented Molecular Monte Carlo code. Classy is designed to be used for Metropolis Monte Carlo
and other similar sampling methods in molecular simulations.  It was created to streamline the implimentation of new
methods in a way that was quick and modular. 


--Installation

ClassyMC requires a Fortran-2003 standard compatible fortran compiler.  It has been tested on both Intel 
and GNU fortran based compilers. For parallel implimentations MPI is required. If one wishes to install with the AENet forcefield one must compile AENet as a static library and import it at compile time.

To perform a basic compile specify the make command using your compiler of choice.

    make gfortran  # GNU Fortran Compiler    
    make ifort  #Intel Compiler
    
To build with AEnet use the command

    make aenet

To build as a shared library use the command

    make lib
    
    
