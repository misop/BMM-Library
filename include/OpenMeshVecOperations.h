#pragma once
#include <OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>

float dot(OpenMesh::Vec3f u, OpenMesh::Vec3f v);
OpenMesh::Vec3f cross(OpenMesh::Vec3f u, OpenMesh::Vec3f v);
float sqr_length(OpenMesh::Vec3f u);
OpenMesh::Vec3f operator*(float c, OpenMesh::Vec3f u);

bool OpenMeshVec2iEqual(const OpenMesh::Vec2i& u, const OpenMesh::Vec2i& v);
bool OpenMeshVec3fEqual(const OpenMesh::Vec3f& u, const OpenMesh::Vec3f& v);
bool OpenMeshVec3fZero(const OpenMesh::Vec3f& u);

struct OpenMeshVec2iComp {
	bool operator() (const OpenMesh::Vec2i& lhs, const OpenMesh::Vec2i& rhs) const {
		if (lhs.values_[0] != rhs.values_[0])
			return lhs.values_[0] < rhs.values_[0];
		
		return lhs.values_[1] < rhs.values_[1];
	}
};