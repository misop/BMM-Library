#include "stdafx.h"
#include "Utility.h"

#define BIAS 0.1

using namespace std;

#pragma region Ray Intersections

bool raySphereIntersection(OpenMesh::Vec3f ray_origin, OpenMesh::Vec3f ray_direction, OpenMesh::Vec3f sphere_center, float sphere_radius, float &t_param) {
	glm::vec3 origin(ray_origin[0], ray_origin[1], ray_origin[2]);
	glm::vec3 direction(ray_direction[0], ray_direction[1], ray_direction[2]);
	glm::vec3 center(sphere_center[0], sphere_center[1], sphere_center[2]);
	return raySphereIntersection(origin, direction, center, sphere_radius, t_param);
}

bool raySphereIntersection(glm::vec3 ray_origin, glm::vec3 ray_direction, glm::vec3 sphere_center, float sphere_radius, float &t_param) {
	glm::vec3 oo = ray_origin - sphere_center;
	// A = v^2
	float A = glm::dot(ray_direction, ray_direction);
	// -B = v^T * (o_2 - o_1)
	float B = -2.0 * glm::dot(ray_direction, oo);
	// C = (o_2 - o_1)^2 - r^2
	float C = glm::dot(oo, oo) - sphere_radius * sphere_radius;
	// Discriminant
	float D = B * B - 4.0f * A * C;
	// No collision
	if (D < 0) 
		return false; 

	float sD = sqrtf(D);
	float t1 = 0.5 * (B + sD) / A;
	//if (t1 < Ray.Bias) t1 = Double.MaxValue;
	float t2 = 0.5 * (B - sD) / A;
	//if (t2 < Ray.Bias) t2 = Double.MaxValue;
	float t = (t1 > 0) ? t1 : t2;
	if (t < 0)
		return false;

	t_param = t;
	return true;
}

bool rayTriangleIntersection(OpenMesh::Vec3f ray_origin, OpenMesh::Vec3f ray_direction, OpenMesh::Vec3f V0, OpenMesh::Vec3f V1, OpenMesh::Vec3f V2, float &t_param) {
	//point on triagle T(u, v) = (1 - u - v)*V0 + u*V1 + v*V2;
	//ray: O + d
	//ray triangle intersection: O + d = (1 - u - v)*V0 + u*V1 + v*V2 -> 3 equations solved by Cramer rule
	//triangle edges
	OpenMesh::Vec3f edge1 = V1 - V0;
	OpenMesh::Vec3f edge2 = V2 - V0;
	//calculation of determinant used to calculate u parameter
	OpenMesh::Vec3f pvec = cross(ray_direction, edge2);
	//if determinant is near zero ray lines in plane of triangle
	float det = dot(edge1, pvec);
	if (det > -FLOAT_ZERO && det < FLOAT_ZERO)
		return false;
	float inv_det = 1.0 / det;
	//distance from V0 to ray origin
	OpenMesh::Vec3f tvec = ray_origin - V0;
	//calculate u parameter and test bounds
	float u = dot(tvec, pvec) * inv_det;
	if (u < 0 || u > 1.0)
		return false;
	//prepare to test v parameter
	OpenMesh::Vec3f qvec = cross(tvec, edge1);
	//calculate v parameter and test bounds
	float v = dot(ray_direction, qvec) * inv_det;
	if (v < 0.0 || (u + v) > 1.0)
		return false;
	//calculate t ray intersects triangle
	float t = dot(edge2, qvec) * inv_det;
	//not in the direction of ray
	if (t < BIAS)
		return false;

	t_param = t;
	return true;
}

#pragma endregion

#pragma region Vector functions

OpenMesh::Vec3f getAxisForCross(OpenMesh::Vec3f v) {
	//return the axis corresponding to the smallest non zero magnitude
	OpenMesh::Vec3f u = OpenMesh::Vec3f(fabs(v[0]), fabs(v[1]), fabs(v[2]));
	if ((u[0] <= u[1]) && (u[0] <= u[2]) && (!equal(u[0], 0))) {
		return OpenMesh::Vec3f(1, 0, 0);
	}
	if ((u[1] <= u[0]) && (u[1] <= u[2]) && (!equal(u[1], 0))) {
		return OpenMesh::Vec3f(0, 1, 0);
	}
	return OpenMesh::Vec3f(0, 0, 1);
}

#pragma endregion

glm::vec2 bezier(glm::vec2 P1, glm::vec2 P2, float t) {
	glm::vec2 P0(0.0, 0.0);
	glm::vec2 P3(1.0, 1.0);
	return (pow(1 - t, 3)*P0 + 3*pow(1 - t, 2)*t*P1 + 3*(1 - t)*pow(t, 2)*P2 + pow(t, 3)*P3);
}