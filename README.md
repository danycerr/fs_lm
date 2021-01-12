# fs_lm
FaultSlip project on libmesh

Optimal control formulation of fault slip problem

## Prerequisite

Openmpi
Petsc
Libmesh

## Installation 
if Libmesh is installed with autoconfig tools use configure.sh and Makefile_autoconf

otherwise use configure_mox and Makefile_include


instal on mox machines command after loading the petsc module

../configure --prefix=/u/cerroni/Desktop/fat_home/software/libmesh/mylibmesh --enable-petsc-required
make
make install


to run the 3d example
source configure.sh
make
mpiexec -np 4 example-opt


