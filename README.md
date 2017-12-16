# qevo

A python solver of the Schroedinger equation for 
arbitrary external potential and few particles with 
either bose,boltzmann or fermi symmetry. It is based 
on the exact diagonalization of the Hamiltonian expressed 
on a specific basis set.

The package is structured as a set of modules, each giving a functionality 
needed for the task at hand.  

* In the "basis" directory a minimal collection of basis sets is given, this is 
the main body of the package since it encompasses all the operations needed 
for the representation of the Hamiltonian operation (i.e. matrix elements for 
differential operators, potentials and external fields). 

* In the "quantum" directory the thermal and time evolution frameworks of 
quantum mechanics are implemented.

* The "scripts" directory contains some utility scripts needed if one wants to 
post-process movies out of the output of the application scripts.

* "applications" contains examples of usage on different systems.

* The "data" directory is used to store the matrix element integrals for the 
harmonic basis set.


It is intended to be quite sandbox, one should take the examples and build from 
there his/her own particular system.


