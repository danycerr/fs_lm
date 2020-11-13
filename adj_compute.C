#include<adj_compute.h>
AdjointPost::AdjointPost(EquationSystems & es, Order p_order):LinearElasticityPost{es,p_order}
{
    std::cout<<"adj_post::building stress tensor "<< sys_name_ <<std::endl;
    es_ = & es;
    dim_ = (es_->get_mesh()).mesh_dimension();
    ExplicitSystem & stress_systemtrnsor =
    es_->add_system<ExplicitSystem> (sys_name_);
    
    stress_systemtrnsor.add_variable("sigma_adj_00", CONSTANT, MONOMIAL);
    stress_systemtrnsor.add_variable("sigma_adj_01", CONSTANT, MONOMIAL);
    if (dim_>2) stress_systemtrnsor.add_variable("sigma_adj_02", CONSTANT, MONOMIAL);
    stress_systemtrnsor.add_variable("sigma_adj_11", CONSTANT, MONOMIAL);
    if (dim_>2) stress_systemtrnsor.add_variable("sigma_adj_12", CONSTANT, MONOMIAL);
    if (dim_>2) stress_systemtrnsor.add_variable("sigma_adj_22", CONSTANT, MONOMIAL);
    
}

// ====================================================================
void AdjointPost::compute(){
    std::cout << "computing stress tensor"<<std::endl;
    const Real young_modulus = 1;//es_->parameters.get<Real>("young_modulus");
    const Real poisson_ratio = 1;//es_->parameters.get<Real>("poisson_ratio");
    const Real gamma = es_->parameters.get<Real> ("gamma");
    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem & _sys = es_->get_system<LinearImplicitSystem> ("adjoint");
    
    const MeshBase & mesh = es_->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    
    std::vector< unsigned int > displacement_vars;
    displacement_vars.push_back( _sys.variable_number ("zx"));
    displacement_vars.push_back( _sys.variable_number ("zy"));
    if(dim_>2) displacement_vars.push_back( _sys.variable_number ("zz"));
    
    //     const unsigned int u_var = _sys.variable_number ("ux"); 
    
    const DofMap & dof_map = _sys.get_dof_map();
    FEType fe_type = dof_map.variable_type(displacement_vars[0]);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
    
    // Also, get a reference to the ExplicitSystem
    ExplicitSystem & stress_system = es_->get_system<ExplicitSystem>(sys_name_);
    const DofMap & stress_dof_map = stress_system.get_dof_map();
    std::vector<unsigned int> sigma_vars;
    sigma_vars.push_back(stress_system.variable_number ("sigma_adj_00"));
    sigma_vars.push_back(stress_system.variable_number ("sigma_adj_01"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_adj_02"));
    sigma_vars.push_back(stress_system.variable_number ("sigma_adj_11"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_adj_12"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_adj_22"));
    
    
    // Storage for the stress dof indices on each element
    std::vector<std::vector<dof_id_type>> dof_indices_var(_sys.n_vars());
    std::vector<dof_id_type> stress_dof_indices_var;
    
    // To store the stress tensor on each element
    DenseMatrix<Number> elem_avg_stress_tensor(dim_, dim_);
    
    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        for (unsigned int var=0; var<dim_; var++)
            dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
        
        const unsigned int n_var_dofs = dof_indices_var[0].size();        
        fe->reinit (elem);
        // clear the stress tensor
        elem_avg_stress_tensor.resize(dim_, dim_);
        
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            // Row is variable u1, u2, or u3, column is x, y, or z
            DenseMatrix<Number> grad_u(dim_, dim_);
            for (unsigned int var_i=0; var_i<dim_; var_i++)
                for (unsigned int var_j=0; var_j<dim_; var_j++)
                    for (unsigned int j=0; j<n_var_dofs; j++)
                        grad_u(var_i, var_j) +=
                        dphi[j][qp](var_j) * _sys.current_solution(dof_indices_var[var_i][j]);
                    
                    DenseMatrix<Number> stress_tensor(dim_, dim_);
                for (unsigned int i=0; i<dim_; i++)
                    for (unsigned int j=0; j<dim_; j++)
                        for (unsigned int k=0; k<dim_; k++)
                            for (unsigned int l=0; l<dim_; l++)
                                stress_tensor(i,j) +=
                                elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * grad_u(k,l);
                            
                            // We want to plot the average stress on each element, hence
                            elem_avg_stress_tensor.add(JxW[qp], stress_tensor);
        }
        
        // Get the average stress per element by dividing by volume
        elem_avg_stress_tensor.scale(1./elem->volume());
        
        // load elem_sigma data into stress_system
        unsigned int stress_var_index = 0;
        for (unsigned int i=0; i<dim_; i++)
            for (unsigned int j=i; j<dim_; j++)
            {
                stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);
                
                // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
                // one dof index per variable
                dof_id_type dof_index = stress_dof_indices_var[0];
                
                if ((stress_system.solution->first_local_index() <= dof_index) &&
                    (dof_index < stress_system.solution->last_local_index()))
                    stress_system.solution->set(dof_index, elem_avg_stress_tensor(i,j));
                
                stress_var_index++;
            }
    }//end element loop
    
    // Should call close and update when we set vector entries directly
    stress_system.solution->close();
    stress_system.update();
}
