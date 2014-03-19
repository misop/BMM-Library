#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <glm/glm.hpp>
#include "MyMesh.h"
#include "FloatArithmetic.h"

#pragma region Template Functions

template <typename T> std::string ToStringT(T a) {
	std::stringstream ss;
	ss << a;
	return ss.str();
}

template <typename T> void fillUp(std::vector<T> &v, T val, int count) {
	for (int i = 0; i < count; i++) {
		v.push_back(val);
	}
}

#pragma endregion

#pragma region Ray Intersections
bool raySphereIntersection(OpenMesh::Vec3f ray_origin, OpenMesh::Vec3f ray_direction, OpenMesh::Vec3f sphere_center, float sphere_radius, float &t_param);
bool raySphereIntersection(glm::vec3 ray_origin, glm::vec3 ray_direction, glm::vec3 sphere_center, float sphere_radius, float &t_param);
bool rayTriangleIntersection(OpenMesh::Vec3f ray_origin, OpenMesh::Vec3f ray_direction, OpenMesh::Vec3f V0, OpenMesh::Vec3f V1, OpenMesh::Vec3f V2, float &t_param);
#pragma endregion

#pragma region Vector functions
OpenMesh::Vec3f getAxisForCross(OpenMesh::Vec3f v);
#pragma endregion

glm::vec2 bezier(glm::vec2 P1, glm::vec2 P2, float t);