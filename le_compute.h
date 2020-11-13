#ifndef STRESS_PROB_HS
#define STRESS_PROB_HS
#include<iostream>
//libmesh include
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/exact_solution.h"
#include "libmesh/enum_solver_package.h"
using namespace libMesh;
class LinearElasticityPost{
public:
    LinearElasticityPost(EquationSystems & es,
                         Order p_order);
    void compute();
protected:
    std::string sys_name_="stress";
    EquationSystems * es_;
    int dim_;
    Real elasticity_tensor(Real young_modulus,
                       Real poisson_ratio,
                       unsigned int i,
                       unsigned int j,
                       unsigned int k,
                       unsigned int l);
    inline Real kronecker_delta(unsigned int i,
                                unsigned int j)
    {
        return i == j ? 1. : 0.;
    }
    
} ; //end le stress compute
#endif //HS
