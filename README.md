# finiteTMPS
Finite temperature tensor network algorithms including METTS and the ancilla/purification method.

The codes were used in the article:
["Matrix product state techniques for two-dimensional systems at finite temperature"](https://arxiv.org/abs/1705.05578), Benedikt Bruognolo, Zhenyue Zhu, Steven R. White, E.M. Stoudenmire (arxiv:[1705.05578](https://arxiv.org/abs/1705.05578))

# Brief description of codes

- `triangular_metts.cc`: minimally entangled typical thermal states (METTS) algorithm for the 
  triangular lattice Heisenberg model on quasi two-dimensional cylinders

- `mpo_ancilla.cc`: ancilla (a.k.a. purification) algorithm for the 
  Heisenberg model on quasi two-dimensional cylinders

# Steps to build

All of the codes require the ITensor library (http://itensor.org). 

1. Download and install the ITensor library somewhere on your machine.
2. Create your own copy of the Makefile.default file (say Makefile.yourname), 
   and edit the LIBRARY_DIR variable to point to where the compiled ITensor source is located (this is 
   the folder containing the options.mk file in it).
3. Create soft link to your make file: `ln -s Makefile.yourname Makefile`
4. Run `make app=appname` to compile a specific code or just `make` to compile the last one in the list.

# More Details and Input Parameters

## `mpo_ancilla` code

This code uses the _ancilla_, also known as the _purification_, approach to obtaining finite temperature properties of quantum lattice models. In a nutshell, the ancilla approach works with a system defined on a two-leg ladder with sites on the first leg interacting with the Hamiltonian whose properties we want to study, and the second leg not interacting at all. The degrees on freedom on the second leg (even-numbered sites in the code) are the _ancilla_ sites since they just play a helper role of thermalizing the physical sites by being entangled with them. One evolves the system under imaginary time evolution to a time beta/2; then taking expectation values of operators acting on the physical sites gives their expected value in the canonical ensemble at temperature T=1/beta.

Inputs recognized:

- Nx (integer): number of physical sites along the x direction
- Ny (integer): number of physical sites along the y direction
- periodic (yes/no): whether to use periodic boundary conditions along the y direction
- lattice_type (string): "triangular" for triangular lattice; "square" for square lattice
- beta (real): inverse temperature of thermal ensemble
- tau (real): imaginary time step to use; smaller value is more accurate but slower to run
- maxm (integer): maximum bond dimension of MPS allowed during time evolution
- cutoff (real): truncation error cutoff used during time evolution
- realstep (yes/no): whether to use a real time step with O(tau^2) error at each time step or two imaginary time steps as a trick to get an O(tau^3) error at each time step
- Jz (real): XXZ Hamiltonian Jz parameter (default=1.0)
- Jxy (real): XXZ Hamiltonian Jxy parameter (default=1.0)


