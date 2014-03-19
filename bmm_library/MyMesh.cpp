#include "stdafx.h"
#include "MyMesh.h"
#include "FloatArithmetic.h"

bool OpenMeshEqualVHandles(MyMesh::VHandle& h1, MyMesh::VHandle& h2) {
	return h1.idx() == h2.idx();
}

#pragma region Mesh Converting

void convertTriMeshToArray(MyTriMesh *mesh, std::vector<float> &points, std::vector<int> &indices) {
	int length = (int)floor((float)points.size() / 3.0);
	points.reserve(length + (mesh->n_vertices() * 3));
	indices.reserve(indices.size() + (mesh->n_faces() * 3));

	for (MyTriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it) {
		MyTriMesh::Point P = mesh->point(v_it.handle());
		points.push_back(P[0]);
		points.push_back(P[1]);
		points.push_back(P[2]);
	}

	for (MyTriMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		MyTriMesh::FaceVertexIter fv_it = mesh->fv_begin(f_it.handle());
		for ( ; fv_it != mesh->fv_end(f_it.handle()); ++fv_it) {
			indices.push_back(fv_it.handle().idx() + length);
		}
	}
}

void convertMeshToArray(MyMesh *mesh, std::vector<float> &points, std::vector<int> &indices) {
	int length = (int)floor((float)points.size() / 3.0);
	points.reserve(length + (mesh->n_vertices() * 3));
	indices.reserve(indices.size() + (mesh->n_faces() * 4));

	for (MyTriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it) {
		MyTriMesh::Point P = mesh->point(v_it.handle());
		points.push_back(P[0]);
		points.push_back(P[1]);
		points.push_back(P[2]);
	}

	for (MyTriMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		MyTriMesh::FaceVertexIter fv_it = mesh->fv_begin(f_it.handle());
		std::vector<int> vertices;
		for ( ; fv_it != mesh->fv_end(f_it.handle()); ++fv_it) {
			vertices.push_back(fv_it.handle().idx() + length);
		}
		addToConvertedMeshArray(vertices, indices);
	}
}

void addToConvertedMeshArray(std::vector<int> &vertices, std::vector<int> &indices) {
	if (vertices.size() >= 3) {
		indices.push_back(vertices[0]);
		indices.push_back(vertices[1]);
		indices.push_back(vertices[2]);
	}
	if (vertices.size() == 4) {		
		indices.push_back(vertices[2]);
		indices.push_back(vertices[3]);
		indices.push_back(vertices[0]);
	}
}

void convertMeshToArray(MyMesh *mesh, std::vector<float> &points, std::vector<int> &triIndices, std::vector<int> &quadIndices) {
	int length = (int)floor((float)points.size() / 3.0);
	points.reserve(length + (mesh->n_vertices() * 3));
	triIndices.reserve(triIndices.size() + (mesh->n_faces() * 3));
	quadIndices.reserve(quadIndices.size() + (mesh->n_faces() * 4));

	for (MyTriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it) {
		MyTriMesh::Point P = mesh->point(v_it.handle());
		points.push_back(P[0]);
		points.push_back(P[1]);
		points.push_back(P[2]);
	}

	for (MyTriMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		MyTriMesh::FaceVertexIter fv_it = mesh->fv_begin(f_it.handle());
		std::vector<int> vertices;
		for ( ; fv_it != mesh->fv_end(f_it.handle()); ++fv_it) {
			vertices.push_back(fv_it.handle().idx() + length);
		}
		addToConvertedMeshArray(vertices, triIndices, quadIndices);
	}
}

void convertMeshToArray(MyMesh *mesh, std::vector<float> &points, std::vector<float> &vertex_normals, std::vector<int> &triIndices, std::vector<int> &quadIndices) {
	int length = (int)floor((float)points.size() / 3.0);
	points.reserve(length + (mesh->n_vertices() * 3));
	vertex_normals.reserve(length + (mesh->n_vertices() * 3));
	triIndices.reserve(triIndices.size() + (mesh->n_faces() * 3));
	quadIndices.reserve(quadIndices.size() + (mesh->n_faces() * 4));

	mesh->request_vertex_normals();
	mesh->request_face_normals();
	mesh->update_normals();

	for (MyTriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it) {
		MyTriMesh::Point P = mesh->point(v_it.handle());
		points.push_back(P[0]);
		points.push_back(P[1]);
		points.push_back(P[2]);
		OpenMesh::Vec3f n = mesh->normal(v_it.handle());
		vertex_normals.push_back(n[0]);
		vertex_normals.push_back(n[1]);
		vertex_normals.push_back(n[2]);
	}

	mesh->release_face_normals();
	mesh->release_vertex_normals();

	for (MyTriMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		MyTriMesh::FaceVertexIter fv_it = mesh->fv_begin(f_it.handle());
		std::vector<int> vertices;
		for ( ; fv_it != mesh->fv_end(f_it.handle()); ++fv_it) {
			vertices.push_back(fv_it.handle().idx() + length);
		}
		addToConvertedMeshArray(vertices, triIndices, quadIndices);
	}
}

void addToConvertedMeshArray(std::vector<int> &vertices, std::vector<int> &triIndices, std::vector<int> &quadIndices) {
	if (vertices.size() == 3) {
		triIndices.push_back(vertices[0]);
		triIndices.push_back(vertices[1]);
		triIndices.push_back(vertices[2]);
	}
	if (vertices.size() == 4) {		
		quadIndices.push_back(vertices[0]);
		quadIndices.push_back(vertices[1]);
		quadIndices.push_back(vertices[2]);
		quadIndices.push_back(vertices[3]);
	}
}

void calculateTriMeshNormals(MyTriMesh *mesh, std::vector<float> &points) {
	for (MyTriMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		std::vector<MyTriMesh::VHandle> vhandles;
		for (MyTriMesh::FVIter fv_it = mesh->fv_begin(f_it.handle()); fv_it != mesh->fv_end(f_it.handle()); ++fv_it) {
			vhandles.push_back(fv_it.handle());
		}
		MyTriMesh::Point P0 = mesh->point(vhandles[0]);
		MyTriMesh::Point P1 = mesh->point(vhandles[1]);
		MyTriMesh::Point P2 = mesh->point(vhandles[2]);
		OpenMesh::Vec3f u = P1 - P0;
		OpenMesh::Vec3f v = P2 - P0;
		OpenMesh::Vec3f normal = cross(u, v).normalize();
		MyTriMesh::Point center = (P0 + P1 + P2) / 3.0;
		OpenMesh::Vec3f dest = center + (normal * 50);
		points.push_back(center[0]);
		points.push_back(center[1]);
		points.push_back(center[2]);

		points.push_back(dest[0]);
		points.push_back(dest[1]);
		points.push_back(dest[2]);
	}
}

void calculateMeshNormals(MyMesh *mesh, std::vector<float> &points) {
	for (MyMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		std::vector<MyTriMesh::VHandle> vhandles;
		for (MyMesh::FVIter fv_it = mesh->fv_begin(f_it.handle()); fv_it != mesh->fv_end(f_it.handle()); ++fv_it) {
			vhandles.push_back(fv_it.handle());
		}
		MyMesh::Point P0 = mesh->point(vhandles[0]);
		MyMesh::Point P1 = mesh->point(vhandles[1]);
		MyMesh::Point P2 = mesh->point(vhandles[2]);
		OpenMesh::Vec3f u = P1 - P0;
		OpenMesh::Vec3f v = P2 - P0;
		OpenMesh::Vec3f normal = cross(u, v).normalize();
		MyMesh::Point center;
		if (vhandles.size() == 3) {
			center = (P0 + P1 + P2) / 3.0;
		} else {
			MyMesh::Point P3 = mesh->point(vhandles[3]);
			center = (P0 + P1 + P2 + P3) / 4.0;
		}
		OpenMesh::Vec3f dest = center + (normal * 50);
		points.push_back(center[0]);
		points.push_back(center[1]);
		points.push_back(center[2]);

		points.push_back(dest[0]);
		points.push_back(dest[1]);
		points.push_back(dest[2]);
	}
}

int calculateMaxValency(MyMesh *mesh) {
	int max = 0;

	for (MyMesh::VIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it) {
		int count = 0;
		for (MyMesh::VVIter vv_it = mesh->vv_begin(v_it.handle()); vv_it != mesh->vv_end(v_it.handle()); ++vv_it) {
			count++;
		}
		if (count > max) {
			max = count;
		}
	}

	return max;
}

float calculateOneRingRadius(MyTriMesh *mesh, MyTriMesh::VHandle vh) {
	MyTriMesh::VHandle first_vh(-1);
	MyTriMesh::VHandle prev_vh(-1);
	float radius = 0;

	for (MyTriMesh::VVIter vv_it = mesh->vv_begin(vh); vv_it != mesh->vv_end(vh); ++vv_it) {
		if (first_vh.idx() == -1) {
			first_vh = vv_it.handle();
			prev_vh = vv_it.handle();
		}
		//calculate length
		MyTriMesh::Point P = mesh->point(prev_vh);
		MyTriMesh::Point Q = mesh->point(vv_it.handle());
		radius += (P - Q).norm();

		prev_vh = vv_it.handle();
	}
	MyTriMesh::Point P = mesh->point(prev_vh);
	MyTriMesh::Point Q = mesh->point(first_vh);
	radius += (P - Q).norm();
	//we calculated circumference (2*PI*r) so now we need to divide by 2*PI
	radius /= (2*M_PI);

	return radius;
}

#pragma endregion

#pragma region Mesh Writing

void writeTriMesh(MyTriMesh* mesh) {
	static int number = 0;
	number++;
	std::stringstream ss;//create a stringstream
	ss << number;
	std::string file = ss.str() + ".obj";
	writeTriMesh(mesh, file);
}

void writeTriMesh(MyTriMesh* mesh, std::string fileName) {
	OpenMesh::IO::Options wopt;
	if (!OpenMesh::IO::write_mesh(*mesh, fileName)) {
	}
}

void writeMesh(MyMesh* mesh, std::string fileName) {
	OpenMesh::IO::Options wopt;
	if (!OpenMesh::IO::write_mesh(*mesh, fileName)) {
	}
}

#pragma endregion