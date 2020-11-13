export COM_SOFT=/home3/comp_software/
export OPENMPI_4=$COM_SOFT/openmpi/openmpi_4/
export OPENMPI_4_BIN=$OPENMPI_4/bin/
export OPENMPI_4_LIB=$OPENMPI_4/lib64/
export OPENMPI_4_INC=$OPENMPI_4/include/
export PETSC_DIR=$COM_SOFT/petsc/petsc_3p13/
export PETSC3p13_LIB=$PETSC_DIR/lib/
export PETSC3p13_INC=$PETSC_DIR/include/

export LIBMESH_DIR=$COM_SOFT/libmesh/libmesh_lv/
export LIBMESH_BIN=$LIBMESH_DIR/bin/
export LIBMESH_LIB=$LIBMESH_DIR/lib64/
export LIBMESH_INC=$LIBMESH_DIR/include/

#Add to path

export PATH=$LIBMESH_BIN:$OPENMPI_4_BIN:$PATH
export LD_LIBRARY_PATH=$LIBMESH_LIB:$PETSC3p13_LIB:$OPENMPI_4_LIB:$LD_LIBRARY_PATH
 
