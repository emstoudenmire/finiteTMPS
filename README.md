# finiteTMPS
Finite temperature tensor network algorithms including METTS and the ancilla/purification method.

The codes were used in the article:
("Matrix product state techniques for two-dimensional systems at finite temperature")[https://arxiv.org/abs/1705.05578], Benedikt
Bruognolo, Zhenyue Zhu, Steven R. White, E.M. Stoudenmire

# Brief description of codes

- `triangular_metts.cc`: minimally entangled typical thermal states (METTS) algorithm for the 
  triangular lattice Heisenberg model on quasi two-dimensional cylinders


# Steps to build

All of the codes require the ITensor library (http://itensor.org). 

1. Download and install the ITensor library somewhere on your machine.
2. Create your own copy of the Makefile.default file (say Makefile.yourname), 
   and edit the LIBRARY_DIR variable to point to where the compiled ITensor source is located (this is 
   the folder containing the options.mk file in it).
3. Create soft link to your make file: `ln -s Makefile.yourname Makefile`
4. Run `make app=appname` to compile a specific code or just `make` to compile the last one in the list.
