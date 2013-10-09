/*
 *  SphereDelaunay.h
 *  gtm
 *
 *  Created by J. Andreas Bærentzen on 04/12/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __SPHERE_DELAUNAY_H__
#define __SPHERE_DELAUNAY_H__

#include <vector>
#include <list>
//#include <GEL/CGLA/Vec3d.h>
//#include <GEL/CGLA/Vec3i.h>
#include "OpenMeshVecOperations.h"


void Delaunay_on_sphere(const std::vector<OpenMesh::Vec3f>& vertices,
						std::vector<OpenMesh::Vec3i>& triangles_out);


#endif