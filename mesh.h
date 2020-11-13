#ifndef MESH_HS
#define MESH_HS
//mesh class
#include<iostream>
// libmesh include
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_modification.h"
// class LibMeshInit;
// class Mesh;
// Bring in everything from the libMesh namespace
using namespace libMesh;
class mesh{
public:
    mesh(libMesh::LibMeshInit* init, std::string name, int ref=0);
    libMesh::Mesh* get(){return &mesh_;}
    ~mesh();
    void init();
private:
    libMesh::Mesh mesh_;
    libMesh::Mesh* mesh_ptr_;
    std::string mesh_name;
};
#endif 
