SQM-Library
===========

Quick Doc:
===========
class SQMBasicController executes SQM Algorithm and exports results

	void loadSkeletonFromFile(string fileName); - loads a skeleton class from file
	void loadSkeleton(SkeletonNode *skeleton); - loads a skeleton class
	void exportMeshToFile(string fileName); - exports mesh as a .obj file
	void exportMeshToTriangles(std::vector<float> &points, std::vector<int> &indices); - exports mesh as list of points (x1, y1, z1, x2, y2, z2, ...) and indices of triangles (0-2 first triangle, 3-5 second triangle, ...) which create mesh faces
	MyMesh* getMesh(); - return OpenMesh mesh object

	void restart(); - resets the algorithm
	void straightenSkeleton(); - straightens skeleton
	void computeConvexHull(); - computes convex hulls of BNPs
	void subdivideConvexHull(); - subdivides convex hulls
	void joinBNPs(); - joins BNPs
	void executeSQMAlgorithm(); - executes SQM algorithm to the end
	void executeSQMAlgorithm(SQMState state); - executes SQM algorithm up to a specified state
  
Common Ussage:
===========
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(fileName); or sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm();
	sqm->exportMeshToTriangles(points, indices); or sqm->exportMeshToFile(saveToFileName);
	delete sqm;