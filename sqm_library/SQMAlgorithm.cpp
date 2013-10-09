#include "stdafx.h"
#include "SQMAlgorithm.h"
#include <OpenMesh\Core\System\mostream.hh>

SQMAlgorithm::SQMAlgorithm(void) : root(NULL)
{
	drawingMode = 0;
	sqmState = SQMStart;
	mesh = NULL;
	fb = new filebuf();
	os = new ostream(fb);
	omerr().connect(*os);
	root = new SQMNode();
	resetRoot = NULL;
	numOfNodes = 1;
	smoothingAlgorithm = SQMAvaragingSmoothing;
}

SQMAlgorithm::~SQMAlgorithm(void)
{
	if (mesh != NULL) {
		delete mesh;
	}
	if (fb->is_open())
		fb->close();
	if (root) delete root;
	if (resetRoot) delete resetRoot;
	delete fb;
	delete os;
}

#pragma region Setters

void SQMAlgorithm::setRoot(SQMNode *newRoot) {
	if (root) delete root;

	root = newRoot;

	drawingMode = 0;
	sqmState = SQMStart;
	numOfNodes = countNodes();
}

void SQMAlgorithm::setNumberOfNodes(int newNumberOfNodes) {
	numOfNodes = newNumberOfNodes;
}

void SQMAlgorithm::setSmoothingAlgorithm(SQMSmoothingAlgorithm sqmSmoothingAlgorithm) {
	smoothingAlgorithm = sqmSmoothingAlgorithm;
}

#pragma endregion

#pragma region Getters

SQMNode* SQMAlgorithm::getRoot() {
	return root;
}

SQMState SQMAlgorithm::getState() {
	return sqmState;
}

MyMesh* SQMAlgorithm::getMesh() {
	return mesh;
}

int SQMAlgorithm::getNumberOfNodes() {
	return numOfNodes;
}

#pragma endregion

#pragma region Tree Functions

int SQMAlgorithm::countNodes() {
	deque<SQMNode*> queue;
	queue.push_back(root);
	int count = 0;

	while (!queue.empty()) {
		SQMNode* node = queue.front();
		queue.pop_front();

		count++;

		vector<SQMNode*> *childs = node->getNodes();
		for (int i = 0; i < childs->size(); i++) {
			queue.push_back((*childs)[i]);
		}
	}
	return count;
}

void SQMAlgorithm::refreshIDs() {
	deque<SQMNode*> queue;
	queue.push_back(root);
	unsigned int count = 0;

	while (!queue.empty()) {
		SQMNode* node = queue.front();
		queue.pop_front();

		node->setID(count);
		count++;

		vector<SQMNode*> *childs = node->getNodes();
		for (int i = 0; i < childs->size(); i++) {
			queue.push_back((*childs)[i]);
		}
	}
}


#pragma endregion

void SQMAlgorithm::getBoundingSphere(float &x, float &y, float &z, float &d) {
	if (!root) {//if there is no root there is no bounding sphere
		x = y = z = d = 0;
		return;
	}
	vector<SQMNode*> stack;
	stack.push_back(root);
	float minX = root->getPosition()[0], maxX = root->getPosition()[0];//initialize the starting variables
	float minY = root->getPosition()[1], maxY = root->getPosition()[1];
	float minZ = root->getPosition()[2], maxZ = root->getPosition()[2];
	while (!stack.empty()) {
		SQMNode* node = stack.back();
		stack.pop_back();
		if (node->getPosition()[0] < minX)
			minX = node->getPosition()[0];
		if (node->getPosition()[0] > maxX)
			maxX = node->getPosition()[0];
		if (node->getPosition()[1] < minY)
			minY = node->getPosition()[1];
		if (node->getPosition()[1] > maxY)
			maxY = node->getPosition()[1];
		if (node->getPosition()[2] < minZ)
			minZ = node->getPosition()[2];
		if (node->getPosition()[2] > maxZ)
			maxZ = node->getPosition()[2];
		stack.insert(stack.end(), node->getNodes()->begin(), node->getNodes()->end());
	}
	float dx = maxX - minX;
	float dy = maxY - minY;
	float dz = maxZ - minZ;
	float maxD = dx > dy ? dx : dy;
	d = maxD > dz ? maxD : dz;
	d += 20;
	x = (maxX - minX)/2 + minX;
	y = (maxY - minY)/2 + minY;
	z = (maxZ - minZ)/2 + minZ;
}

#pragma region SQM Algorithm

void SQMAlgorithm::updateResetRoot() {
	if (resetRoot) delete resetRoot;
	resetRoot = new SQMNode(*root);
}

void SQMAlgorithm::restart() {
	if (root) delete root;
	root = new SQMNode(*resetRoot);
	sqmState = SQMStart;
}

void SQMAlgorithm::straightenSkeleton() {
	updateResetRoot();
	refreshIDs();

	fb->open("log.txt", ios::out);
	(*os) << "Skeleton straightening\n";
	//root->straightenSkeleton(OpenMesh::Vec3f(0, 0, 0));
	root->straightenSkeleton(NULL);
	sqmState = SQMStraighten;
}

void SQMAlgorithm::computeConvexHull() {
	(*os) << "Convex hull computation\n";
	vector<SQMNode*> stack;
	stack.push_back(root);
	while (!stack.empty()) {
		SQMNode *node = stack.back();
		stack.pop_back();
		if (node->isBranchNode())
			node->calculateConvexHull();
		stack.insert(stack.end(), node->getNodes()->begin(), node->getNodes()->end());
	}
	drawingMode = 2;
	sqmState = SQMComputeConvexHull;
}

void SQMAlgorithm::subdivideConvexHull() {
	(*os) << "Convex hull subdivision\n";
	root->subdividePolyhedra(NULL, 0, smoothingAlgorithm);
	sqmState = SQMSubdivideConvexHull;
}

void SQMAlgorithm::joinBNPs() {
	//deletition is messed up :-(
	//maybe delete just faces and let vertices be?
	//filebuf fb;
	//fb.open ("test.txt", ios::out);
	//ostream os(&fb);
	(*os) << "BNP joining\n";
	//omerr().connect(os);

	if (mesh != NULL) {
		delete mesh;
	}
	mesh = new MyMesh();
	vector<MyMesh::VertexHandle> vector;
	root->joinBNPs(mesh, NULL, vector, OpenMesh::Vec3f(0, 0, 0));
	drawingMode = 3;
	sqmState = SQMJoinBNPs;

	fb->close();
}

void SQMAlgorithm::finalVertexPlacement() {
	(*os) << "Final vertex placement\n";
	root->rotateBack(mesh);
	sqmState = SQMFinalPlacement;
}

void SQMAlgorithm::executeSQMAlgorithm() {
	straightenSkeleton();
	computeConvexHull();
	subdivideConvexHull();
	joinBNPs();
	finalVertexPlacement();
}


void SQMAlgorithm::executeSQMAlgorithm(SQMState state) {
	if (sqmState > state) {
		restart();
	}
	if (sqmState < state) {
		if (sqmState < SQMStraighten && state >= SQMStraighten) {
			straightenSkeleton();
		}
		if (sqmState < SQMComputeConvexHull && state >= SQMComputeConvexHull) {
			computeConvexHull();
		}
		if (sqmState < SQMSubdivideConvexHull && state >= SQMSubdivideConvexHull) {
			subdivideConvexHull();
		}
		if (sqmState < SQMJoinBNPs && state >= SQMJoinBNPs) {
			joinBNPs();
		}
		if (sqmState < SQMFinalPlacement && state >= SQMFinalPlacement) {
			finalVertexPlacement();
		}
	}
}

#pragma endregion
