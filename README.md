
# fs_lm


<p align="center">
  <img src="complpan.png" width="700">
</p>


FaultSlip project on libmesh

Optimal control formulation of fault slip problem

## Prerequisite

* Openmpi
* Petsc (3.13)
* Libmesh (1.6.0)

## Installation 
if Libmesh is installed with autoconfig modify and tools use configure.sh and copy  Makefile_autoconf in Makefile

otherwise modify and  use configure_mox and then copy Makefile_include in Makefile

### To install libmesh on mox machines
after loading the petsc module

```
../configure --prefix=youprefix --enable-petsc-required
make
make install
```
to run the 3d example
source configure.sh
make
mpiexec -np 4 example-opt

## To run the code
### 2d example
```
make
mpiexec -np 4 faultslip test_2d.in 
```


### 3d example
```
make
mpiexec -np 4 faultslip test_3d.in 
```

## View solution
The solution is printed on the file sol.e open it in paraview.

load all the variables you want to acces (lower left corner under properties tab)

Attention the evaluated jump is in b and b_tot and not the variable beta and beta_tot

To acces the data on the fault use the command plotoverline (*ctrl + space* type plot over line  ) the set the points of the line (lower left corner)
attention if the points are on the interface between elements the visualization may not be correct so a eps shoud be used. For example instead of a line between pt1 (0,0.5,0.5) and pt2 (0,0.5,0.5) use the points pta1(0,0.501,0.501) and pta2(1,0.501,0.501)
