#pragma once
#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "SQMAlgorithm.h"

class SQMBasicController
{
public:
private:
	SQMAlgorithm *sqmALgorithm;
public:	
	SQMBasicController(void);
	~SQMBasicController(void);

	void loadSkeletonFromFile(string fileName);
	void loadSkeleton(SQMSkeletonNode *skeleton);
	void saveSkeletonToFile(string fileName);
	void exportMeshToFile(string fileName);
	void exportMeshToTriangles(std::vector<float> &points, std::vector<int> &indices);
	MyMesh* getMesh();

	void restart();
	void straightenSkeleton();
	void computeConvexHull();
	void subdivideConvexHull();
	void joinBNPs();
	void executeSQMAlgorithm();
	void executeSQMAlgorithm(SQMState state);
};

