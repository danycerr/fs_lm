# fs_lm
FaultSlip project on libmesh

Optimal control formulation of fault slip problem

Prerequisite

Openmpi
Petsc
Libmesh

Installation 
if Libmesh is installed with autoconfig tools use configure.sh and Makefile_autoconf
otherwise use configure_mox and Makefile_include

to run the 3d example
source configure.sh
make
mpiexec -np 4 example-opt


