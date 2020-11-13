#include<mesh.h>
mesh::mesh(libMesh::LibMeshInit *init, std::string name):mesh_(init->comm())
{
    mesh_ptr_ = new libMesh::Mesh(init->comm());
    mesh_.read(name);
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
