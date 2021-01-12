#include<jump_handler.h>

JumpHandler::JumpHandler(EquationSystems & es,
                         Order p_order)
{
    std::cout<<"adj_post::building stress tensor "<< sys_name_ <<std::endl;
    es_ = & es;
    dim_ = (es_->get_mesh()).mesh_dimension();
    ExplicitSystem & adjoint_system =
    es_->add_system<ExplicitSystem> (sys_name_);
    
    adjoint_system.add_variable("beta", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("beta_tot", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("bx", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("by", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("bz", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("b_totx", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("b_toty", CONSTANT, MONOMIAL);
    adjoint_system.add_variable("b_totz", CONSTANT, MONOMIAL);
    
    
}

// ======================================================================00
void JumpHandler::update()
{
    std::cout <<"compute jump "<<std::endl;
    const Real alpha = es_->parameters.get<Real> ("alpha");
    const Real gamma = es_->parameters.get<Real> ("gamma");
    const Real mu = es_->parameters.get<Real> ("mu_l");
    const MeshBase & mesh = es_->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    // Also, get a reference to the ExplicitSystem
    ExplicitSystem & stress_system = es_->get_system<ExplicitSystem>("adj_stress");
    const DofMap & stress_dof_map = stress_system.get_dof_map();
    std::vector<unsigned int> sigma_vars;
    sigma_vars.push_back(stress_system.variable_number ("sigma_adj_00"));
    sigma_vars.push_back(stress_system.variable_number ("sigma_adj_01"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_adj_02"));
    sigma_vars.push_back(stress_system.variable_number ("sigma_adj_11"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_adj_12"));
    if(dim_>2) sigma_vars.push_back(stress_system.variable_number ("sigma_adj_22"));
    
    
    // Get a reference to the LinearImplicitSystem we are solving
    ExplicitSystem & stress_system_u = es_->get_system<ExplicitSystem>("stress");
    unsigned int xi_var = stress_system_u.variable_number ("xi");
    unsigned int tau_var = stress_system_u.variable_number ("sigma_01");
    unsigned int tau_xz_var = -1;
    if(dim>2) tau_xz_var = stress_system_u.variable_number ("sigma_02");
    unsigned int tau_xy_var = stress_system_u.variable_number ("sigma_01");
    unsigned int rho_var = stress_system_u.variable_number ("sigma_00");
    unsigned int st_var = stress_system_u.variable_number ("ST_n");
    const DofMap & dof_map_stress_u = stress_system_u.get_dof_map();
    std::vector<dof_id_type> dof_indicies_rho;
    std::vector<dof_id_type> dof_indicies_tau;
    std::vector<dof_id_type> dof_indicies_tau_xy;
    std::vector<dof_id_type> dof_indicies_tau_xz;
    std::vector<dof_id_type> dof_indicies_xi;
    std::vector<dof_id_type> dof_indicies_ST;
    
    
    
    
    
    // Also, get a reference to the ExplicitSystem
    ExplicitSystem & beta_sys = es_->get_system<ExplicitSystem>(sys_name_);
    const DofMap & beta_dof_map = beta_sys.get_dof_map();
    unsigned int beta_var, beta_var_tot;
    beta_var = beta_sys.variable_number ("beta");
    beta_var_tot = beta_sys.variable_number ("beta_tot");
    unsigned int beta_v_var, beta_v_var_tot;
    beta_v_var = beta_sys.variable_number ("bx");
    beta_v_var_tot = beta_sys.variable_number ("b_totx");
    
    
    // Storage for the stress dof indices on each element
    std::vector<std::vector<dof_id_type>> dof_indices_var(beta_sys.n_vars());
    std::vector<std::vector<dof_id_type>> dof_indices_betav_var(dim);
    std::vector<std::vector<dof_id_type>> dof_indices_betav_tot_var(dim);
    std::vector<dof_id_type> dof_indicies_sigma;
    std::vector<std::vector<dof_id_type>>  dof_indicies_sigma_adj(dim);
    
    FEType fe_type = stress_dof_map.variable_type(sigma_vars[1]);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    const std::vector<Point> & qpoint = fe->get_xyz();
    
    const std::vector<Real> & JxW = fe->get_JxW();
    
    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        bool fracture=false;
        for (auto side : elem->side_index_range()){
            if (mesh.get_boundary_info().has_boundary_id(elem, side, 3)){
                fracture=true;break;}
        }
        beta_dof_map.dof_indices (elem, dof_indices_var[0], beta_var);
        beta_dof_map.dof_indices (elem, dof_indices_var[1], beta_var_tot);
        for(int idim=0; idim<dim; idim++){
        beta_dof_map.dof_indices (elem, dof_indices_betav_var[idim], beta_v_var+idim);
        beta_dof_map.dof_indices (elem, dof_indices_betav_tot_var[idim], beta_v_var_tot+idim);
        stress_dof_map.dof_indices (elem,dof_indicies_sigma_adj[idim],sigma_vars[idim]); // taglio
            
        }
        
        stress_dof_map.dof_indices (elem,dof_indicies_sigma,sigma_vars[1]); // taglio aggiunto
        dof_map_stress_u.dof_indices (elem,dof_indicies_xi,xi_var);
        dof_map_stress_u.dof_indices (elem,dof_indicies_ST,st_var);
        dof_map_stress_u.dof_indices (elem,dof_indicies_rho,rho_var);
        dof_map_stress_u.dof_indices (elem,dof_indicies_tau,tau_var);
        dof_map_stress_u.dof_indices (elem,dof_indicies_tau_xy,tau_xy_var);
         if(dim>2) dof_map_stress_u.dof_indices (elem,dof_indicies_tau_xz,tau_xz_var);
        fe->reinit (elem);
        
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            if( fracture ){
//           if( fracture && fabs(qpoint[qp](1) - 0.5)<0.375 &&  fabs(qpoint[qp](2) - 0.5)<0.375 ){
                {
                    dof_id_type dof_index= dof_indices_var[0][0];
                    dof_id_type dof_index_tot= dof_indices_var[1][0];
                    double xi=stress_system_u.current_solution(dof_indicies_xi[0]);
                    double st=stress_system_u.current_solution(dof_indicies_ST[0]);
                    double tau=stress_system_u.current_solution(dof_indicies_tau[0]);
                    double rho=stress_system_u.current_solution(dof_indicies_rho[0])*mu;
                    double beta_k = beta_sys.current_solution(dof_index_tot);
                    double D=(1-xi)*gamma*beta_k + xi*tau;
                    double N=(tau-gamma*beta_k);
                    double dbJd = 0.;
                    if (rho> 0 && rho<fabs(N)) dbJd = gamma*(1-xi)*D - gamma*rho*D/(fabs(N) +1.e-16);
                    else if (rho> fabs(N))  dbJd = gamma*D;
                    //            std::cout << "rho is "<< rho <<  " N is "<<N << " and dbJd "<< dbJd<<std::endl;
                    double beta_val =-alpha* (stress_system.current_solution(dof_indicies_sigma[0]) + dbJd);//forma style
                    //            double beta_val =-alpha* fabs(std::min(0.0,st-mu)); // uzawa st
                    //            double beta_val =-alpha* fabs(D); // uzawa st
                    double beta_tot=beta_sys.current_solution(dof_index_tot)+beta_val;
                    beta_sys.solution->set(dof_index,beta_val);
                    beta_sys.solution->set(dof_index_tot,beta_tot);
                }
                std::vector<double> tau_v={0.,stress_system_u.current_solution(dof_indicies_tau_xy[0])};
                if (dim>2) tau_v.push_back(stress_system_u.current_solution(dof_indicies_tau_xz[0]));
//               
                std::vector<double>  tau_adj_v={0.};
                double normN=0.;
                for(int idim =0; idim<dim;idim++){
                    if(idim>0)
                        tau_adj_v.push_back(stress_system.current_solution(dof_indicies_sigma_adj[idim][0]));
                    normN+=fabs(tau_v[idim]-gamma*beta_sys.current_solution(dof_indices_betav_tot_var[idim][0]));
                }
                for (int idim=0; idim<dim; idim++){
                    dof_id_type dof_index= dof_indices_betav_var[idim][0];
                    dof_id_type dof_index_tot= dof_indices_betav_tot_var[idim][0];
                    double xi=stress_system_u.current_solution(dof_indicies_xi[0]);
                    double st=stress_system_u.current_solution(dof_indicies_ST[0]);
                    double tau=tau_v[idim];
                    double rho=stress_system_u.current_solution(dof_indicies_rho[0])*mu;
                    double beta_k = beta_sys.current_solution(dof_index_tot);
                    double D=(1-xi)*gamma*beta_k + xi*tau;
                    double dbJd = 0.;
                    if (rho> 0 && rho<normN) dbJd = gamma*(1-xi)*D - gamma*rho*D/(normN+1.e-16);
                    else if (rho> normN)  dbJd = gamma*D;
//                     if (idim==1)std::cout << "dof "<<  " rho is "<< rho <<  " N is "<<normN << " and dbJd "<< dbJd<< " D "<<D<< " xi "<< xi << " gamma "<< gamma<<std::endl;
//                     std::cout <<dbJd<< " "<< tau << " "<< xi<<std::endl;
                    double beta_val =-alpha* (tau_adj_v[idim] + dbJd);//forma style
                    //            double beta_val =-alpha* fabs(std::min(0.0,st-mu)); // uzawa st
                    //            double beta_val =-alpha* fabs(D); // uzawa st
                    double beta_tot=beta_sys.current_solution(dof_index_tot)+beta_val;
                    beta_sys.solution->set(dof_index,beta_val);
                    beta_sys.solution->set(dof_index_tot,beta_tot);
                }
                
            }
            
        }
    }
    beta_sys.solution->close();
    beta_sys.update();
}
