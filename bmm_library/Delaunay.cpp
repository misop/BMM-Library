#include "stdafx.h"
#include "Delaunay.h"
#include <math.h>
#include "FloatArithmetic.h"

using namespace std; 

////////////////////////////////////////////////////////////////////////
// CircumCircle() :
//   Return true if a point (xp,yp) is inside the circumcircle made up
//   of the points (x1,y1), (x2,y2), (x3,y3)
//   The circumcircle centre is returned in (xc,yc) and the radius r
//   Note : A point on the edge is inside the circumcircle
////////////////////////////////////////////////////////////////////////

int CircumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3, double &xc, double &yc, double &r){
	double m1, m2, mx1, mx2, my1, my2;
	double dx, dy, rsqr, drsqr;

	/* Check for coincident points */
	if(fabs(y1 - y2) < EPSILON && abs(y2 - y3) < EPSILON)
		return(false);
	if(abs(y2-y1) < EPSILON){ 
		m2 = - (x3 - x2) / (y3 - y2);
		mx2 = (x2 + x3) / 2.0;
		my2 = (y2 + y3) / 2.0;
		xc = (x2 + x1) / 2.0;
		yc = m2 * (xc - mx2) + my2;
	}else if(abs(y3 - y2) < EPSILON){ 
		m1 = - (x2 - x1) / (y2 - y1);
		mx1 = (x1 + x2) / 2.0;
		my1 = (y1 + y2) / 2.0;
		xc = (x3 + x2) / 2.0;
		yc = m1 * (xc - mx1) + my1;
	}else{
		m1 = - (x2 - x1) / (y2 - y1); 
		m2 = - (x3 - x2) / (y3 - y2); 
		mx1 = (x1 + x2) / 2.0; 
		mx2 = (x2 + x3) / 2.0;
		my1 = (y1 + y2) / 2.0;
		my2 = (y2 + y3) / 2.0;
		xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2); 
		yc = m1 * (xc - mx1) + my1; 
	}
	dx = x2 - xc;
	dy = y2 - yc;
	rsqr = dx * dx + dy * dy;
	r = sqrt(rsqr); 
	dx = xp - xc;
	dy = yp - yc;
	drsqr = dx * dx + dy * dy;
	return((drsqr <= rsqr) ? true : false);
}
///////////////////////////////////////////////////////////////////////////////
// Triangulate() :
//   Triangulation subroutine
//   Takes as input NV vertices in array pxyz
//   Returned is a list of ntri triangular faces in the array v
//   These triangles are arranged in a consistent clockwise order.
//   The triangle array 'v' should be malloced to 3 * nv
//   The vertex array pxyz must be big enough to hold 3 more points
//   The vertex array must be sorted in increasing x values say
//
//   qsort(p,nv,sizeof(XYZ),XYZCompare);
///////////////////////////////////////////////////////////////////////////////

int Triangulate(vector<glm::vec3> &pxyz, vector<glm::ivec3> &v) {
	int nv = pxyz.size();
	int ntri = 0;
	vector<int> complete;
	vector<glm::ivec2> edges;
	int nedge = 0;
	int trimax, emax = 200;
	int status = 0;
	int inside;
	int i, j, k;
	double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;
	double xmin, xmax, ymin, ymax, xmid, ymid;
	double dx, dy, dmax; 

	/* Allocate memory for the completeness list, flag for each triangle */
	trimax = 4 * nv;
	complete.resize(trimax);
	/* Allocate memory for the edge list */
	edges.resize(emax);
	/*
	Find the maximum and minimum vertex bounds.
	This is to allow calculation of the bounding triangle
	*/
	xmin = pxyz[0].x;
	ymin = pxyz[0].y;
	xmax = xmin;
	ymax = ymin;
	for(i = 1; i < nv; i++){
		if (pxyz[i].x < xmin) xmin = pxyz[i].x;
		if (pxyz[i].x > xmax) xmax = pxyz[i].x;
		if (pxyz[i].y < ymin) ymin = pxyz[i].y;
		if (pxyz[i].y > ymax) ymax = pxyz[i].y;
	}
	dx = xmax - xmin;
	dy = ymax - ymin;
	dmax = (dx > dy) ? dx : dy;
	xmid = (xmax + xmin) / 2.0;
	ymid = (ymax + ymin) / 2.0;
	/*
	Set up the supertriangle
	his is a triangle which encompasses all the sample points.
	The supertriangle coordinates are added to the end of the
	vertex list. The supertriangle is the first triangle in
	the triangle list.
	*/
	pxyz.push_back(glm::vec3(xmid - 20 * dmax, ymid - dmax, 0));
	pxyz.push_back(glm::vec3(xmid, ymid + 20 * dmax, 0));
	pxyz.push_back(glm::vec3(xmid + 20 * dmax, ymid - dmax, 0));
	v.push_back(glm::ivec3(nv, nv+1, nv+2));
	complete[0] = false;
	ntri = 1;
	/*
	Include each point one at a time into the existing mesh
	*/
	for(i = 0; i < nv; i++){
		xp = pxyz[i].x;
		yp = pxyz[i].y;
		nedge = 0;
		/*
		Set up the edge buffer.
		If the point (xp,yp) lies inside the circumcircle then the
		three edges of that triangle are added to the edge buffer
		and that triangle is removed.
		*/
		for(j = 0; j < ntri; j++){
			if(complete[j])
				continue;
			x1 = pxyz[v[j].x].x;
			y1 = pxyz[v[j].x].y;
			x2 = pxyz[v[j].y].x;
			y2 = pxyz[v[j].y].y;
			x3 = pxyz[v[j].z].x;
			y3 = pxyz[v[j].z].y;
			inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r);
			if (xc + r < xp)
				// Suggested
					// if (xc + r + EPSILON < xp)
						complete[j] = true;
			if(inside){
				/* Check that we haven't exceeded the edge list size */
				if(nedge + 3 >= emax){
					emax += 100;
					edges.resize(emax);
				}
				edges[nedge+0].x = v[j].x;
				edges[nedge+0].y = v[j].y;
				edges[nedge+1].x = v[j].y;
				edges[nedge+1].y = v[j].z;
				edges[nedge+2].x = v[j].z;
				edges[nedge+2].y = v[j].x;
				nedge += 3;
				//v[j] = v[ntri-1];
				v[j] = v.back();
				v.pop_back();
				complete[j] = complete[ntri-1];
				ntri--;
				j--;
			}
		}
		/*
		Tag multiple edges
		Note: if all triangles are specified anticlockwise then all
		interior edges are opposite pointing in direction.
		*/
		for(j = 0; j < nedge - 1; j++){
			for(k = j + 1; k < nedge; k++){
				if((edges[j].x == edges[k].y) && (edges[j].y == edges[k].x)){
					edges[j].x = -1;
					edges[j].y = -1;
					edges[k].x = -1;
					edges[k].y = -1;
				}
				/* Shouldn't need the following, see note above */
				if((edges[j].x == edges[k].x) && (edges[j].y == edges[k].y)){
					edges[j].x = -1;
					edges[j].y = -1;
					edges[k].x = -1;
					edges[k].y = -1;
				}
			}
		}
		/*
		Form new triangles for the current point
		Skipping over any tagged edges.
		All edges are arranged in clockwise order.
		*/
		for(j = 0; j < nedge; j++) {
			if(edges[j].x < 0 || edges[j].y < 0)
				continue;
			v.push_back(glm::ivec3(edges[j].x, edges[j].y, i));
			complete[ntri] = false;
			ntri++;
		}
	}
	/*
	Remove triangles with supertriangle vertices
	These are triangles which have a vertex number greater than nv
	*/
	for(i = 0; i < ntri; i++) {
		if(v[i].x >= nv || v[i].y >= nv || v[i].z >= nv) {
			//v[i] = v[ntri-1];
			v[i] = v.back();
			v.pop_back();
			ntri--;
			i--;
		}
	}

	return 0;
} 