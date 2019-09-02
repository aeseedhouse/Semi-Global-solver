# Semi-Global Solver

This is a program for efficiently solving the time dependent Schrödinger equation for a Hamiltonian that can also have time dependence. A polynomial series approximation is used to avoid computationally intensive matrix exponentiation. The time step over which the polynomial approximations are calculated are chosen adaptively; a time step that is too large leads to a diverging solution, whereas a small time step abolishes the advantages of using the semi-global approach.  

## Getting Started

The file testRabi.m gives an example of how the program can calculate the wavefunction from a Hamiltonian (in this case rabiHam.m). The user can create their own function for a chosen Hamiltonian with or without time dependence, the "calcH" handle in testRabi will have to be changed in this case. 

## Authors

* **Hayden Ooi** - *Initial work* - [hooi98](https://github.com/hooi98)
