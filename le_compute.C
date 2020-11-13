#include<le_compute.h>

LinearElasticityPost::LinearElasticityPost(EquationSystems & es, Order p_order){
    std::cout<<"le_post::building stress tensor and xi"<<std::endl;
    es_ = & es;
    dim_ = (es_->get_mesh()).mesh_dimension();
    ExplicitSystem & stress_systemtrnsor =
    es_->add_system<ExplicitSystem> (sys_name_);
    
    stress_systemtrnsor.add_variable("sigma_00", CONSTANT, MONOMIAL);
    stress_systemtrnsor.add_variable("sigma_01", CONSTANT, MONOMIAL);
    if (dim_>2) stress_systemtrnsor.add_variable("sigma_02", CONSTANT, MONOMIAL);
    stress_systemtrnsor.add_variable("sigma_11", CONSTANT, MONOMIAL);
    if (dim_>2) stress_systemtrnsor.add_variable("sigma_12", CONSTANT, MONOMIAL);
    if (dim_>2) stress_systemtrnsor.add_variable("sigma_22", CONSTANT, MONOMIAL);
    stress_systemtrnsor.add_variable("ST_n", CONSTANT, MONOMIAL);
    stress_systemtrnsor.add_variable("xi", CONSTANT, MONOMIAL);
    
}
//=================================================================================
Real LinearElasticityPost::elasticity_tensor(Real young_modulus,
                                             Real poisson_ratio,
                                             unsigned int i,
                                             unsigned int j,
                                             unsigned int k,
                                             unsigned int l)
{
    // Define the Lame constants
    Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
    Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));
    lambda_1 =1.; lambda_2 =1.; 
    return lambda_1 * kronecker_delta(i, j) * kronecker_delta(k, l) +
    lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
}



//=================================================================================
void LinearElasticityPost::compute(){
    std::cout << "computing stress tensor"<<std::endl;
    const Real young_modulus = 1;//es_->parameters.get<Real>("young_modulus");
    const Real poisson_ratio = 1;//es_->parameters.get<Real>("poisson_ratio");
    const Real gamma = es_->parameters.get<Real> ("gamma");
    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem & _sys = es_->get_system<LinearImplicitSystem> ("mechanics");
    
    // Also, get a reference to the ExplicitSystem
    ExplicitSystem & beta_sys = es_->get_system<ExplicitSystem>("jump");
    const DofMap & beta_dof_map = beta_sys.get_dof_map();
    unsigned int beta_var;
    beta_var = beta_sys.variable_number ("beta");
    
    const MeshBase & mesh = es_->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    
    const Real mu_l = es_->parameters.get<Real> ("mu_l");
    std::vector< unsigned int > displacement_vars;
    displacement_vars.push_back( _sys.variable_number ("ux"));
    displacement_vars.push_back( _sys.variable_number ("uy"));
    if(dim_>2) displacement_vars.push_back( _sys.variable_number ("uz"));
    
    std::vector< unsigned int > jump_vars;
    jump_vars.push_back( beta_sys.variable_number ("bx"));
    jump_vars.push_back( beta_sys.variable_number ("by"));
    if(dim_>2) jump_vars.push_back( beta_sys.variable_number ("bz"));
    
    const unsigned int u_var = _sys.variable_number ("ux"); 
    
    const DofMap & dof_map = _sys.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
    
    // Also, get a reference to the ExplicitSystem
    ExplicitSystem & stress_system = es_->get_system<ExplicitSystem>("stress");
    const DofMap & stress_dof_map = stress_system.get_dof_map();
    std::vector<unsigned int> sigma_vars;
    sigma_vars.push_back(stress_system.variable_number ("sigma_00"));
    sigma_vars.push_back(stress_system.variable_number ("sigma_01"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_02"));
    sigma_vars.push_back(stress_system.variable_number ("sigma_11"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_12"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_22"));
    unsigned int st_var = stress_system.variable_number ("ST_n");
    unsigned int xi_var = stress_system.variable_number ("xi");
    
    
    // Storage for the stress dof indices on each element
    std::vector<std::vector<dof_id_type>> dof_indices_var(_sys.n_vars());
    std::vector<std::vector<dof_id_type>> dof_indvec_jump(dim_);
    std::vector<dof_id_type> stress_dof_indices_var;
    std::vector<dof_id_type> dof_ind_jump;
    std::vector<std::vector<dof_id_type>> dof_ind_jumpvec(beta_sys.n_vars());
    // To store the stress tensor on each element
    DenseMatrix<Number> elem_avg_stress_tensor(dim_, dim_);
    
    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        
        
        beta_dof_map.dof_indices (elem, dof_ind_jump,beta_var);
        for (unsigned int var=0; var<dim_; var++){
            dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
            beta_dof_map.dof_indices (elem, dof_indvec_jump[var], jump_vars[var]);
            }
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
                            // we integrate stress_tensor
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
            
            // project in the tangential direction
            double tau=fabs(elem_avg_stress_tensor(1,0));
            if (dim_>2) tau+=fabs(elem_avg_stress_tensor(2,0));
            Number ST_value = tau / elem_avg_stress_tensor(0,0) ;
            
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, st_var);
            dof_id_type dof_index = stress_dof_indices_var[0];
            
            if ((stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()))
                stress_system.solution->set(dof_index, ST_value);
            
            // //     if(fabs(beta_sys.current_solution(dof_ind_jump[0]) ) > 1.e-12)  std::cout <<" adj "<<beta_sys.current_solution(dof_ind_jump[0])<< std::endl;
            //building norm N
            Number norm_N = 1.e-16;
            for (int idim=1; idim<dim_; idim++) //working only for plane dof_ind_jump
//                 norm_N += fabs( elem_avg_stress_tensor(idim,0) - gamma*beta_sys.current_solution(dof_ind_jump[0]) );
                norm_N += fabs(elem_avg_stress_tensor(idim,0) - gamma*beta_sys.current_solution(dof_indvec_jump[idim][0]));
            Number xi_value = ((ST_value < mu_l)? 1.:0)*(1-  ( mu_l*elem_avg_stress_tensor(0,0))/norm_N); // missing gamma beta
            if (xi_value <0 ) xi_value =0;
            //             Number xi_value =0;
            if (elem_avg_stress_tensor(0,0) >0) xi_value =1;
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, xi_var);
            dof_index = stress_dof_indices_var[0];
            
            if ((stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()))
                stress_system.solution->set(dof_index, xi_value);
    }
    
    // Should call close and update when we set vector entries directly
    stress_system.solution->close();
    stress_system.update();
}
