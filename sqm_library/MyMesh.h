#pragma once
#include <OpenMesh\Core\IO\MeshIO.hh>
#include <OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
#include <OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> MyMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyTriMesh;


bool OpenMeshEqualVHandles(MyMesh::VHandle& h1, MyMesh::VHandle& h2);

void convertTriMeshToArray(MyTriMesh *mesh, std::vector<float> &points, std::vector<int> &indices);
void convertMeshToArray(MyMesh *mesh, std::vector<float> &points, std::vector<int> &indices);
void addToConvertedMeshArray(std::vector<int> &vertices, std::vector<int> &indices);
void convertMeshToArray(MyMesh *mesh, std::vector<float> &points, std::vector<int> &triIndices, std::vector<int> &quadIndices);
void convertMeshToArray(MyMesh *mesh, std::vector<float> &points, std::vector<float> &vertex_normals, std::vector<int> &triIndices, std::vector<int> &quadIndices);
void addToConvertedMeshArray(std::vector<int> &vertices, std::vector<int> &triIndices, std::vector<int> &quadIndices);

void calculateTriMeshNormals(MyTriMesh *mesh, std::vector<float> &points);
void calculateMeshNormals(MyMesh *mesh, std::vector<float> &points);
int calculateMaxValency(MyMesh *mesh);
float calculateOneRingRadius(MyTriMesh *mesh, MyTriMesh::VHandle vh);

void writeTriMesh(MyTriMesh* mesh);
void writeTriMesh(MyTriMesh* mesh, std::string fileName);
void writeMesh(MyMesh* mesh, std::string fileName);