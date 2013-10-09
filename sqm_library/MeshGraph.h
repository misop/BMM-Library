#pragma once
#include <boost\unordered_map.hpp>
#include <tnt/tnt_array2d.h>
#include "m_math.h"

using namespace TNT;

struct MeshGraph {
    int numOfVertices;
    CVector3 * pVerts;            // The object's vertices
    boost::unordered_map<int, int> indices; // mGindex = meshGraph->indices[index + offset]; 
    //boost::unordered_map<int, vector<int>> inverseIndices; // inverseIndices[mgIndex] = vector of original triangle indices + offset
	Array2D<bool> E;
    float * wH;
    float * wHorig;
    float wL;
    float faceAreaSum;
    int numOfFaces;
    float * origOneRingArea;
    float * origOneRingExtent;
    // vbo variables
    bool createdVBO;
    int numOfVBOlines;
    unsigned int gMeshgraphPositionsVB;
    unsigned int gMeshgraphLinesVB;

    MeshGraph(){
        pVerts = new CVector3[0];
        wH = new float[0];
        wHorig = new float[0];
        origOneRingArea = new float[0];
        createdVBO = false;
        gMeshgraphPositionsVB = -1;
        gMeshgraphLinesVB = -1;
    }
    ~MeshGraph(){
        delete[] pVerts;
        pVerts = NULL;
        delete[] wH;
        wH = NULL;
        delete[] wHorig;
        wHorig = NULL;
        delete[] origOneRingArea;
        origOneRingArea = NULL;

    }

};