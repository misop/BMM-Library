#include "stdafx.h"
#include "OpenMeshVecOperations.h"
#include "FloatArithmetic.h"

using namespace OpenMesh;

float dot(OpenMesh::Vec3f u, OpenMesh::Vec3f v) {
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

OpenMesh::Vec3f cross(OpenMesh::Vec3f u, OpenMesh::Vec3f v) {
	float values[3];
	values[0] = u.values_[1]*v.values_[2] - u.values_[2]*v.values_[1];
	values[1] = u.values_[2]*v.values_[0] - u.values_[0]*v.values_[2];
	values[2] = u.values_[0]*v.values_[1] - u.values_[1]*v.values_[0];
	return OpenMesh::Vec3f(values); 
}

float sqr_length(OpenMesh::Vec3f u) {
	return dot(u, u);
}

Vec3f operator*(float c, OpenMesh::Vec3f u) {
	//return Vec3f(c*u.values_[0], c*u.values_[1], c*u.values_[2]);
	return u*c;
}

bool OpenMeshVec2iEqual(const OpenMesh::Vec2i& u, const OpenMesh::Vec2i& v) {
	return (u.values_[0] == v.values_[0] && u.values_[1] == v.values_[1]) || (u.values_[1] == v.values_[0] && u.values_[0] == v.values_[1]);
}

bool OpenMeshVec3fEqual(const OpenMesh::Vec3f& u, const OpenMesh::Vec3f& v) {
	return (equal(u.values_[0], v.values_[0]) && equal(u.values_[1], v.values_[1]) && equal(u.values_[2], v.values_[2]));
}

bool OpenMeshVec3fZero(const OpenMesh::Vec3f& u) {
	return equal(u.values_[0], 0) && equal(u.values_[1], 0) && equal(u.values_[2], 0);
}