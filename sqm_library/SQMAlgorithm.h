#pragma once
#include "SQMNode.h"
#include <fstream>
#include <ostream>

typedef enum {
	SQMStart = 0,
	SQMStraighten,
	SQMComputeConvexHull,
	SQMSubdivideConvexHull,
	SQMJoinBNPs,
	SQMFinalPlacement
} SQMState;

class SQMAlgorithm {
	SQMNode *root;
	SQMNode *resetRoot;
	int drawingMode;
	SQMState sqmState;
	MyMesh* mesh;
	filebuf *fb;
	ostream *os;
	int numOfNodes;
	SQMSmoothingAlgorithm smoothingAlgorithm;
public:
	SQMAlgorithm(void);
	~SQMAlgorithm(void);

	void setRoot(SQMNode *newRoot);
	void setNumberOfNodes(int newNumberOfNodes);
	void setSmoothingAlgorithm(SQMSmoothingAlgorithm sqmSmoothingAlgorithm);

	SQMNode* getRoot();
	SQMState getState();
	MyMesh* getMesh();
	int getNumberOfNodes();

	void straightenSkeleton();

	int countNodes();
	void refreshIDs();

	void getBoundingSphere(float &x, float &y, float &z, float &d);
	void updateResetRoot();
	void restart();
	void computeConvexHull();
	void subdivideConvexHull();
	void joinBNPs();
	void finalVertexPlacement();
	void executeSQMAlgorithm();
	void executeSQMAlgorithm(SQMState state);
};

