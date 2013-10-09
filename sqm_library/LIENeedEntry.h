#pragma once
#include "LIE.h"

#define NO_SPLITTING -2
#define UNLIMITED_SPLITTING -1


struct LIENeedEntry {
	MyTriMesh::VertexHandle vhandle;
	int need;
	std::vector<LIE> lies;

	LIENeedEntry(MyTriMesh::VertexHandle vh, int needed) : vhandle(vh), need(needed) { };
};