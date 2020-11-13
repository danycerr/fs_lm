#ifndef MESH_HS
#define MESH_HS
//mesh class
#include<iostream>
// libmesh include
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
// class LibMeshInit;
// class Mesh;
class mesh{
public:
    mesh(libMesh::LibMeshInit* init, std::string name);
    libMesh::Mesh* get(){return &mesh_;}
    ~mesh();
    void init();
private:
    libMesh::Mesh mesh_;
    libMesh::Mesh* mesh_ptr_;
    std::string mesh_name;
};
#endif 
