#include<adjoint_elasticity.h>
AdjointProblem::AdjointProblem(EquationSystems & es, Order p_order){
    std::cout << "Constructor of adjoint problem"<<std::endl;
    es_g = & es;
    // Create a system named le_system
    LinearImplicitSystem & system=es.add_system<libMesh::LinearImplicitSystem> (sys_name);
    system.add_variable ("zx", p_order, libMesh::MONOMIAL);
    system.add_variable ("zy", p_order, libMesh::MONOMIAL);
    dim_ = (es.get_mesh()).mesh_dimension();
    if (dim_ > 2 )  system.add_variable ("zz", p_order, libMesh::MONOMIAL);
      // Give the system a pointer to the matrix assembly function
    system.attach_assemble_function (assemble_prob);
}
//==========================================================================
void AdjointProblem::solve(){
    std::cout <<"Solve of AdjointProblem"<<std::endl;
    LinearImplicitSystem & system = es_g->get_system<LinearImplicitSystem> (sys_name);
    system.solve();
    std::cout << "Adjoint Linear solver converged at step: "
    << system.n_linear_iterations()
    << ", final residual: "
    << system.final_linear_residual()
    << std::endl;
}





// We now define the matrix assembly function for the
// Laplace system.  We need to first compute element volume
// matrices, and then take into account the boundary
// conditions and the flux integrals, which will be handled
// via an interior penalty method.
void assemble_prob(EquationSystems & es,
                         const std::string & libmesh_dbg_var(system_name))
{
  libMesh::out << " assembling adjoint dg system... ";
  libMesh::out.flush();

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "adjoint");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & ellipticdg_system = es.get_system<LinearImplicitSystem> ("adjoint");
  // Get some parameters that we need during assembly
  const Real penalty = es.parameters.get<Real> ("penalty_adj");
  const Real penalty_int = es.parameters.get<Real> ("penalty_adj");
  const Real gamma = es.parameters.get<Real> ("gamma");
  std::string refinement_type = es.parameters.get<std::string> ("refinement");

  
// Get a reference to the others explict system TODO
  ExplicitSystem & stress_system = es.get_system<ExplicitSystem>("stress");
  unsigned int xi_var = stress_system.variable_number ("xi");
  unsigned int tau_xy_var = stress_system.variable_number ("sigma_01");
  unsigned int tau_xz_var =-1;
  if (dim> 2) tau_xz_var = stress_system.variable_number ("sigma_02");
  const DofMap & dof_map_stress = stress_system.get_dof_map();
//   
//     //system for jump correction
  ExplicitSystem & jump_system =
    es.get_system<ExplicitSystem> ("jump");
  unsigned int beta_var = jump_system.variable_number ("beta_tot");
  unsigned int beta_var_tot = jump_system.variable_number ("b_totx");
  const DofMap & dof_map_beta = jump_system.get_dof_map();
  
  // A reference to the DofMap object for this system.  
  const DofMap & dof_map = ellipticdg_system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = ellipticdg_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe  (FEBase::build(dim, fe_type));
  std::unique_ptr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  std::unique_ptr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));

  // Quadrature rules for numerical integration.
#ifdef QORDER
  QGauss qrule (dim, QORDER);
#else
  QGauss qrule (dim, fe_type.default_quadrature_order());
#endif
  fe->attach_quadrature_rule (&qrule);

#ifdef QORDER
  QGauss qface(dim-1, QORDER);
#else
  QGauss qface(dim-1, fe_type.default_quadrature_order());
#endif

  // Tell the finite element object to use our quadrature rule.
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Point> & qpoint = fe->get_xyz();

  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> &  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi_face = fe_elem_face->get_dphi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  // Data for surface integrals on the neighbor boundary
  const std::vector<std::vector<Real>> &  phi_neighbor_face = fe_neighbor_face->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi_neighbor_face = fe_neighbor_face->get_dphi();

  // Define data structures to contain the element interior matrix
  // and right-hand-side vector contribution.  
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  DenseVector<Number> Fee;
  DenseVector<Number> Fnn;

  // Data structures to contain the element and neighbor boundary matrix.
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

  // This vector will hold the degree of freedom indices for the element.  
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_ind_jump;
  
  
//   for (const auto & node : mesh.node_ptr_range()){
//       std::cout << node->get_info()<<std::endl;
//        dof_map.dof_indices (node, dof_indices_node);
//        for (int i =0; i< dof_indices_node.size() ; i++) std::cout << "Node "<<dof_indices_node[i]<<" ";
//       std::cout << std::endl;
// }
  
    std::vector<std::vector<dof_id_type>> dof_indvec_jump(dim);
    std::vector<dof_id_type> dof_indices_xi;
    std::vector<dof_id_type> dof_indices_beta;
    std::vector<dof_id_type> dof_indices_tauxy, dof_indices_tauxz;
  // Now we will loop over all the elements in the mesh. 
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
//          std::cout << elem->get_info()<<std::endl;
      // Get the degree of freedom indices for the
      // current element. 
      dof_map.dof_indices (elem, dof_indices);
//       dof maps of other systems 
      dof_map_stress.dof_indices (elem, dof_indices_xi,xi_var);
      dof_map_beta.dof_indices (elem, dof_indices_beta,beta_var);
      dof_map_stress.dof_indices (elem, dof_indices_tauxy,tau_xy_var);
     if (dim> 2) dof_map_stress.dof_indices (elem, dof_indices_tauxz,tau_xz_var);
      for(int idim=0; idim<dim;idim++)
        dof_map_beta.dof_indices (elem, dof_indvec_jump[idim],beta_var_tot+idim);
      
      
      const unsigned int n_dofs = dof_indices.size();

      
      // Compute the element-specific data for the current element. 
      fe->reinit (elem);

      // Zero the element matrix and right-hand .
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Now we will build the element interior matrix. 
      
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          for (unsigned int ivar=0;ivar<dim;ivar++)
              for (unsigned int i=0; i<n_dofs/dim; i++){
                  double xf  =qpoint[qp](0);
                  //             Fe(i+ivar*n_dofs/dim) += JxW[qp]*(phi[i][qp]*((xf<0.5)?1:-1))*ivar;
//                               Fe(i+ivar*n_dofs/dim) += JxW[qp]*phi[i][qp]*(1)*((ivar==0)? 1:0);
                  for (unsigned int j=0; j<n_dofs/dim; j++){
                      double Lap=0.;
                      for(int iivar=0; iivar<dim; iivar++)
                          Lap+=JxW[qp]*dphi[i][qp](iivar)*dphi[j][qp](iivar);
                      Ke(i+ivar*n_dofs/dim,j+ivar*n_dofs/dim) += Lap ;                    
                      Ke(i + ivar * n_dofs/dim , j + ivar * n_dofs/dim ) +=
                      JxW[qp]*( dphi[i][qp](ivar)*dphi[j][qp](ivar) //l1
                      + dphi[i][qp](ivar)*dphi[j][qp](ivar)
                      ); 
                      for (int ishift=0; ishift<dim-1; ishift++){
                          int shift=(ivar+1+ishift)%dim;
                          Ke(i + ivar * n_dofs/dim , j + shift * n_dofs/dim ) +=
                          JxW[qp]*( dphi[i][qp](shift)*dphi[j][qp](ivar)
                          + dphi[i][qp](ivar)*dphi[j][qp](shift) //l2
                          );
                      }
                      
                  }
              }
  
        
      // Now we address boundary conditions.
      for (auto side : elem->side_index_range())
        {
            if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
                // Pointer to the element face
                fe_elem_face->reinit(elem, side);
                //               const Real xf = side->point(ns)(0);
                //               const Real yf = side->point(ns)(1);
                std::unique_ptr<const Elem> elem_side (elem->build_side_ptr(side));
                // h element dimension to compute the interior penalty penalty parameter
                double penalty_bc=1.e+6;
                const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
                const double h_elem = elem->volume()/elem_side->volume() * 1./pow(elem_b_order, 2.);
                
                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                    // The quadrature point locations on the element side
                    std::vector<Point > qface_point;
                    // Get the physical locations of the element quadrature points
                    qface_point = fe_elem_face->get_xyz();
                    for(int ivar=0; ivar< dim; ivar++) 
                    {
                        Number bc_value = 0.;
                        Number mat_term = 1.;
                         for (unsigned int i=0; i<n_dofs/dim; i++)
                        {
                            // Matrix contribution
                            for (unsigned int j=0; j<n_dofs/dim; j++)
                            {   double xp= qface_point[qp](0);
                                double yp= qface_point[qp](1);
//                              

                               // if ((yp<0.0001 || yp>0.9999) && (ivar==1 || ivar==2)) {bc_value=0.; mat_term=0;} //neumann
//                                 if(xp>0.9999)  if(ivar==0){bc_value=0.2; mat_term=1;}
                                // stability
                                Ke(i+ivar*n_dofs/dim,j+ivar*n_dofs/dim) += JxW_face[qp] * mat_term*penalty_bc/h_elem * phi_face[i][qp] * phi_face[j][qp];
                                // consistency
                                Ke(i+ivar*n_dofs/dim,j+ivar*n_dofs/dim) -=
                                JxW_face[qp] * mat_term*
                                (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) +
                                phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp]));
                            }
                            // RHS contributions
                            // stability
                            Fe(i+ivar*n_dofs/dim) += JxW_face[qp] * bc_value * penalty_bc/h_elem * phi_face[i][qp];
                            // consistency
                            Fe(i+ivar*n_dofs/dim) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);
                        }
                    }
                }
            }

          // If the element is not on a boundary of the domain
          // we loop over his neighbors to compute the element
          // and neighbor boundary matrix contributions
          else 
            {
              // Store a pointer to the neighbor 
              const Elem * neighbor = elem->neighbor_ptr(side);
              // Get the global id of the element and the neighbor
              const unsigned int elem_id = elem->id();
              const unsigned int neighbor_id = neighbor->id();

              // If the neighbor has the same h level and is active
              // perform integration only if our global id is bigger than our neighbor id.
              // We don't want to compute twice the same contributions.
              // If the neighbor has a different h level perform integration
              // only if the neighbor is at a lower level.
              if ((neighbor->active() &&
                   (neighbor->level() == elem->level()) &&
                   (elem_id < neighbor_id)) ||
                  (neighbor->level() < elem->level())
              )
                {
                  // Pointer to the element side
                  std::unique_ptr<const Elem> elem_side (elem->build_side_ptr(side));

                  const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
                  const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
                  const double h_elem = (elem->volume()/elem_side->volume()) * 1./pow(side_order,2.);

                  // The quadrature point locations on the neighbor side
                  std::vector<Point> qface_neighbor_point;

                  // The quadrature point locations on the element side
                  std::vector<Point > qface_point;

                  // Reinitialize shape functions on the element side
                  fe_elem_face->reinit(elem, side);
                  // Get the physical locations of the element quadrature points
                  qface_point = fe_elem_face->get_xyz();

                  // Find their locations on the neighbor
                  unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
                  if (refinement_type == "p")
                    fe_neighbor_face->side_map (neighbor,
                                                elem_side.get(),
                                                side_neighbor,
                                                qface.get_points(),
                                                qface_neighbor_point);
                  else
                    FEInterface::inverse_map (elem->dim(),
                                              fe->get_fe_type(),
                                              neighbor,
                                              qface_point,
                                              qface_neighbor_point);

                  // Calculate the neighbor element shape functions at those locations
                  fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

                  // Get the degree of freedom indices for the neighbor. 
                  std::vector<dof_id_type> neighbor_dof_indices;
                  std::vector<dof_id_type> neighbor_dof_indices_beta;
                  dof_map.dof_indices (neighbor, neighbor_dof_indices);
//                   beta_dof_map.dof_indices (neighbor,neighbor_dof_indices_beta);
                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();

                  // Zero the element and neighbor side matrix before
                  // summing them.  We use the resize member here because
                  // the number of degrees of freedom might have changed from
                  // the last element or neighbor.
                  // Note that Kne and Ken are not square matrices if neighbor
                  // and element have a different p level
                  for (int ivar=0; ivar<dim; ivar++){
                      Kne.resize (n_neighbor_dofs, n_dofs);
                      Ken.resize (n_dofs, n_neighbor_dofs);
                      Kee.resize (n_dofs, n_dofs);
                      Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
                      Fnn.resize(n_neighbor_dofs);
                      Fee.resize(n_neighbor_dofs);
                      
                      
                      double xi=stress_system.current_solution(dof_indices_xi[0]);
//                       double tau=stress_system.current_solution(dof_indices_tau[0]);
                      std::vector<double> tau = {0,stress_system.current_solution(dof_indices_tauxy[0])};
                      if(dim>2) tau.push_back(stress_system.current_solution(dof_indices_tauxz[0]));
//                       double beta=adjoint_system.current_solution(dof_indices_beta[0]);
                      std::vector<double> xiD(dim);
                      for (int idim=0; idim<dim; idim++) {
                          xiD[idim]=(0.5*((1-xi)*gamma*jump_system.current_solution(dof_indvec_jump[idim][0]) +xi*tau[idim]));
                    }
                      double _xiD=  0.0;// 0.5*((1-xi)*gamma*beta + xi*tau );
                      std::vector<double> jump{0,0,fabs(_xiD)};
                      // Now we will build the element and neighbor boundary matrices.  
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                      // Kee Matrix. Integrate the element test function i
                      // against the element test function j
                
//                   bool fracture = ( fabs(qface_point[qp](0)-0.5)< 0.001  && fabs(qface_point[qp](1)-0.5)< 0.375 + 0.001&& fabs(qface_point[qp](2)-0.5)<  0.375 + 0.001  );
                         bool fracture = (!(qface_point[qp](0)<0.49999 ||qface_point[qp](0) >0.50001 ));
//                   fracture=false;
                       double sym=0;
                       for (unsigned int i=0; i<n_dofs/dim; i++)
                       {
                           for (unsigned int j=0; j<n_dofs/dim; j++)
                           {
                               
                               double div_u=0.;
                               for (int idim=0; idim < dim; idim++) div_u+=dphi_face[j][qp](idim);
                               int i_idx=i+ivar*n_dofs/dim;
                               int j_idx=j+ivar*n_dofs/dim;
                               Kee(i_idx,j_idx)-=0.5 * JxW_face[qp] *(2*dphi_face[j][qp](ivar)+div_u)
                               *qface_normals[qp](ivar)*phi_face[i][qp];
                               for (int shift=1; shift<dim; shift++){
                                   int nextvar=(ivar+shift)%dim;
                                   int j_sh_idx=j+nextvar*n_dofs/dim;
                                   Kee(i_idx,j_sh_idx)-=0.5 * JxW_face[qp]*(div_u)*qface_normals[qp](nextvar)*phi_face[i][qp];
                               }  
                               
                               //                               // stability
                               Kee(i+ivar*n_dofs/dim,j+ivar*n_dofs/dim) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
                               Kee(i+ivar*n_dofs/dim,j+ivar*n_dofs/dim) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*qface_normals[qp](ivar)*phi_face[i][qp]*qface_normals[qp](ivar);
                               
                             if(fracture)
                              {
                              Fee(i+ivar*n_dofs/dim) += JxW_face[qp] * penalty_int/h_elem *xiD[ivar]* phi_face[j][qp]*phi_face[i][qp];
                              }
                              
                           }
                       }

                      // Knn Matrix. Integrate the neighbor test function i
                      // against the neighbor test function j
                   for (unsigned int i=0; i<n_neighbor_dofs/dim; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_dofs/dim; j++)
                            {
                                
                               double div_u=0.;
                               for (int idim=0; idim < dim; idim++) div_u+=dphi_neighbor_face[j][qp](idim);
                               int i_idx=i+ivar*n_neighbor_dofs/dim;
                               int j_idx=j+ivar*n_neighbor_dofs/dim;
                               Knn(i_idx,j_idx)-=0.5 * JxW_face[qp] *(2*dphi_neighbor_face[j][qp](ivar)+div_u)
                               *qface_normals[qp](ivar)*phi_neighbor_face[i][qp];
                               for (int shift=1; shift<dim; shift++){
                                   int nextvar=(ivar+shift)%dim;
                                   int j_sh_idx=j+nextvar*n_neighbor_dofs/dim;
                                   Knn(i_idx,j_sh_idx)-=0.5 * JxW_face[qp]*(div_u)*qface_normals[qp](nextvar)*phi_neighbor_face[i][qp];
                               }  
                                

                              // stability
                              Knn(i+ivar*n_neighbor_dofs/dim,j+ivar*n_neighbor_dofs/dim) +=
                                JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                              Knn(i+ivar*n_neighbor_dofs/dim,j+ivar*n_neighbor_dofs/dim) += JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*qface_normals[qp](ivar)*phi_neighbor_face[i][qp]*qface_normals[qp](ivar);
                            }
                        }

                      // Kne Matrix. Integrate the neighbor test function i
                      // against the element test function j
                    for (unsigned int i=0; i<n_neighbor_dofs/dim; i++)
                        {
                          for (unsigned int j=0; j<n_dofs/dim; j++)
                            {
                                
                               double div_u=0.;
                               for (int idim=0; idim < dim; idim++) div_u+=dphi_face[j][qp](idim);
                               int i_idx=i+ivar*n_neighbor_dofs/dim;
                               int j_idx=j+ivar*n_dofs/dim;
                               Kne(i_idx,j_idx)+=0.5 * JxW_face[qp] *(2*dphi_face[j][qp](ivar)+div_u)
                               *qface_normals[qp](ivar)*phi_neighbor_face[i][qp];
                               for (int shift=1; shift<dim; shift++){
                                   int nextvar=(ivar+shift)%dim;
                                   int j_sh_idx=j+nextvar*n_dofs/dim;
                                   Kne(i_idx,j_sh_idx)+=0.5 * JxW_face[qp]*(div_u)*qface_normals[qp](nextvar)*phi_neighbor_face[i][qp];
                               }

                              // stability
                              Kne(i+ivar*n_neighbor_dofs/dim,j+ivar*n_dofs/dim) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];                       
                              Kne(i+ivar*n_neighbor_dofs/dim,j+ivar*n_dofs/dim) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*qface_normals[qp](ivar)*phi_neighbor_face[i][qp]*qface_normals[qp](ivar);
//                               if(fracture) Fnn(i+ivar*n_neighbor_dofs/dim) -= JxW_face[qp] * penalty_int/h_elem *xiD[ivar]* phi_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Ken Matrix. Integrate the element test function i
                      // against the neighbor test function j
                     for (unsigned int i=0; i<n_dofs/dim; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_dofs/dim; j++)
                          {
                              
                                
                               double div_u=0.;
                               for (int idim=0; idim < dim; idim++) div_u+=dphi_neighbor_face[j][qp](idim);
                               int i_idx=i+ivar*n_dofs/dim;
                               int j_idx=j+ivar*n_neighbor_dofs/dim;
                               Ken(i_idx,j_idx)+=0.5 * JxW_face[qp] *(2*dphi_neighbor_face[j][qp](ivar)+div_u)
                               *qface_normals[qp](ivar)*phi_face[i][qp];
                               for (int shift=1; shift<dim; shift++){
                                   int nextvar=(ivar+shift)%dim;
                                   int j_sh_idx=j+nextvar*n_neighbor_dofs/dim;
                                   Ken(i_idx,j_sh_idx)+=0.5 * JxW_face[qp]*(div_u)*qface_normals[qp](nextvar)*phi_face[i][qp];
                               }
                                // stability
                              Ken(i+ivar*n_dofs/dim,j+ivar*n_neighbor_dofs/dim) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                              Ken(i+ivar*n_dofs/dim,j+ivar*n_neighbor_dofs/dim) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*qface_normals[qp](ivar)*phi_neighbor_face[j][qp]*qface_normals[qp](ivar);
//                           
                            }
                        }
                      
                    }
       
                  // The element and neighbor boundary matrix are now built
                  // for this side.  Add them to the global matrix
                  // The SparseMatrix::add_matrix() members do this for us.
                  ellipticdg_system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                  ellipticdg_system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                  ellipticdg_system.matrix->add_matrix(Kee, dof_indices);
                  ellipticdg_system.matrix->add_matrix(Knn, neighbor_dof_indices);
                  ellipticdg_system.rhs->add_vector(Fnn, neighbor_dof_indices);
                  ellipticdg_system.rhs->add_vector(Fee, dof_indices);
                  }
                }
            }
        }
      // The element interior matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      ellipticdg_system.matrix->add_matrix(Ke, dof_indices);
      ellipticdg_system.rhs->add_vector(Fe, dof_indices);
    }

  libMesh::out << "done adjoint dg" << std::endl;
}
