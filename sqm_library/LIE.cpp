#include "stdafx.h"
#include "LIE.h"

bool LIE::containsVertice(int verticeIndex) {
	if (verticeIndex == vertice1 || verticeIndex == vertice2)
		return true;

	return false;
}

int LIE::otherVerticeIndex(int verticeIndex) {
	if (verticeIndex == vertice1)
		return vertice2;
	else
		return vertice1;
}

bool LIE::containsVertices(int v1, int v2) {
	if ((vertice1 == v1 && vertice2 == v2) || (vertice1 == v2 && vertice2 == v1))
		return true;

	return false;
}

bool LIE::containsVertice(MyTriMesh::VHandle verticeHandle) {
	if (vhandle1 == verticeHandle || vhandle2 == verticeHandle)
		return true;

	return false;
}

bool LIE::containsVertices(MyTriMesh::VHandle vh1, MyTriMesh::VHandle vh2) {
	if ((vhandle1 == vh1 && vhandle2 == vh2) || (vhandle1 == vh2 && vhandle2 == vh1))
		return true;

	return false;
}

MyTriMesh::VHandle LIE::otherVerticeHandle(MyTriMesh::VHandle verticeHandle) {
	if (vhandle1 == verticeHandle)
		return vhandle2;

	return vhandle1;
}
