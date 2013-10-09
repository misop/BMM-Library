#pragma once
#include "MyMesh.h"
#include <vector>
#include "m_math.h"

struct LIE {
	int refined;
	int vertice1;
	int vertice2;
	MyTriMesh::VHandle vhandle1;
	MyTriMesh::VHandle vhandle2;
	MyTriMesh::HHandle firstHHandle;
	MyTriMesh::HHandle lastHHandle;
	Quaternion quaternion;

	std::vector<MyTriMesh::EdgeHandle> edges;

	LIE() { refined = 0; };
	LIE(int v1, int v2) : vertice1(v1), vertice2(v2) { refined = 0; };
	LIE(MyTriMesh::VHandle vh1, MyTriMesh::VHandle vh2) : vhandle1(vh1), vhandle2(vh2) { refined = 0; };

	bool operator==(const LIE &other) {
		return containsVertices(other.vertice1, other.vertice2);
	}

	bool containsVertice(int verticeIndex);
	bool containsVertices(int v1, int v2);
	int otherVerticeIndex(int verticeIndex);
	
	bool containsVertice(MyTriMesh::VHandle verticeHandle);
	bool containsVertices(MyTriMesh::VHandle vh1, MyTriMesh::VHandle vh2);
	MyTriMesh::VHandle otherVerticeHandle(MyTriMesh::VHandle verticeHandle);
};