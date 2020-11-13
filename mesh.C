#include<mesh.h>
mesh::mesh(libMesh::LibMeshInit *init, std::string name, int ref):mesh_(init->comm())
{
    mesh_ptr_ = new libMesh::Mesh(init->comm());
    mesh_.read(name);
    libMesh::MeshRefinement mesh_refinement(mesh_);
    for (unsigned int rstep=0; rstep<ref; rstep++)
    {   std::cout<<"refine"<<std::endl;
        mesh_refinement.uniformly_refine(1);
    }
    mesh_.print_info();
    mesh_.get_boundary_info();
    std::cout <<"Constructor of mesh"<<std::endl;
}
// // // // // // // // // // // // // // // // // // // //
mesh::~mesh(){
    std::cout << "Mesh destructor"<<std::endl;
    delete mesh_ptr_;
}
///////////////////////////////////////////////////////////
void mesh::init(){
    std::cout<<"init mesh"<<std::endl;
}
