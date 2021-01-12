#include<iostream>
// libmesh include
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
// include mesh class
#include<mesh.h>
//include le class
#include<linear_elasticity.h>
#include<adjoint_elasticity.h>

#include<le_compute.h>
#include<adj_compute.h>
#include<jump_handler.h>
using Real = libMesh::Real;
using Order = libMesh::Order;
int main (int argc, char** argv)
{
    libMesh::LibMeshInit init(argc, argv);
    // input file
    GetPot input_file("test_3d.in");
   
    //Read in parameters from the input file
    
//      const std::string mesh_name            = input_file("mesh_name","sq.msh");
    const std::string mesh_name            = input_file("mesh_name","cube.msh");
    const unsigned int adaptive_refinement_steps = input_file("max_adaptive_r_steps", 3);
    const unsigned int uniform_refinement_steps  = input_file("uniform_h_r_steps", 3);
    const Real refine_fraction                   = input_file("refine_fraction", 0.5);
    const Real coarsen_fraction                  = input_file("coarsen_fraction", 0.);
    const unsigned int max_h_level               = input_file("max_h_level", 10);
    const std::string refinement_type            = input_file("refinement_type","p");
    Order p_order                                = static_cast<Order>(input_file("p_order", 1));
    const std::string element_type               = input_file("element_type", "tensor");
    const Real penalty                           = input_file("ip_penalty", 10.);
    const Real penalty_int                           = input_file("ip_penalty_int", 10.);
    const Real ip_penalty_adj                           = input_file("ip_penalty_adj", 10.);
    const bool singularity                       = input_file("singularity", true);
    const unsigned int dim                       = input_file("dimension", 3);
    const unsigned int ndiv                       = input_file("ndiv", 8);
    const Real mulim                           = input_file("mulim", -0.038);
    const Real gamma                           = input_file("gamma", 10.);
    const Real alpha                           = input_file("alpha", 0.1);
    const unsigned int nmaxit                  = input_file("nmaxit", 10);
    const unsigned int print_step              = input_file("printstep", 1);
    
    // mesh constriction
    // mesh msh(&init,mesh_name);
    mesh msh(&init,mesh_name,uniform_refinement_steps);
    std::cout<<"Begin Frictional contact"<<std::endl;
    
    // Crate an equation system object
    libMesh::EquationSystems equation_system(*(msh.get()));
    double mytol=1.e-10;
    // Set parameters for the equation system and the solver
    equation_system.parameters.set<Real>("linear solver tolerance") = mytol;
    equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 20000;
    equation_system.parameters.set<Real>("penalty") = penalty;
    equation_system.parameters.set<Real>("penalty_int") = penalty_int;
    equation_system.parameters.set<Real>("penalty_adj") = ip_penalty_adj;
    equation_system.parameters.set<bool>("singularity") = singularity;
    equation_system.parameters.set<std::string>("refinement") = refinement_type;
    equation_system.parameters.set<Real>("mu_l") = mulim;
    equation_system.parameters.set<Real>("gamma") = gamma;
    equation_system.parameters.set<Real>("alpha") = alpha;
    
    
    


    LinearElasticityProblem  le(equation_system,p_order);
    LinearElasticityPost     le_p(equation_system,p_order);
    AdjointProblem           adj(equation_system,p_order);
    AdjointPost              adj_p(equation_system,p_order);
    JumpHandler              jump(equation_system,p_order);
    //    Initialize the data structures for the equation system
    equation_system.init();
    // If Exodus is available, we'll write all timesteps to the same file
    // rather than one file per timestep.
    std::string exodus_filename = "sol.e";
    int i=0;
    ExodusII_IO exo(*(msh.get()));
    ExodusII_IO (*(msh.get())).write_discontinuous_exodusII(exodus_filename, equation_system);

    for(i=0; i<nmaxit;i++)
    {
        std::cout << "+++ iteration +++ "<< i<<" +++" << std::endl;
        // Solve the system
        le.solve();
        // compute stress tensor
        le_p.compute();
        // solve adjoint problem
        adj.solve(); adj_p.compute();
        //update jump
        jump.update();
        if(i%print_step==0){
		int step= (int) i/print_step + 1;
            ExodusII_IO exo(*(msh.get()));
            exo.append(true);
            exo.write_timestep_discontinuous (exodus_filename, equation_system, step, i*1);
        }
    }
    ExodusII_IO (*(msh.get())).write_discontinuous_exodusII("sol_f.e", equation_system);
}
