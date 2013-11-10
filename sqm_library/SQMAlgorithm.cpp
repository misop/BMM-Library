#include "stdafx.h"
#include "SQMAlgorithm.h"
#include <OpenMesh\Core\System\mostream.hh>
#include <glm\gtc\type_ptr.hpp>

#define LOG_COMPUTATION_TIME false

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
	skeleton = NULL;
	numOfSkinMatrices = 0;
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
	if (skeleton) delete skeleton;
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

int SQMAlgorithm::getNumberOfSkinningMatrices() {
	return numOfSkinMatrices;
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

void SQMAlgorithm::calculateSkinSkeletonIDs() {
	if (skeleton == NULL) return;

	int id = 0;
	deque<SkinSkeleton *> queue;
	queue.push_back(skeleton);
	while (!queue.empty()) {
		SkinSkeleton *node = queue.front();
		queue.pop_front();
		node->id = id;
		id++;
		for (int i = 0; i < node->nodes.size(); i++) {
			queue.push_back(node->nodes[i]);
		}
	}

	numOfSkinMatrices = id;
}

void SQMAlgorithm::getSkinningMatrices(float* matrix) {
	if (skeleton == NULL || numOfSkinMatrices == 0) return;

	deque<SkinSkeleton *> queue;
	queue.push_back(skeleton);
	while (!queue.empty()) {
		SkinSkeleton *node = queue.front();
		queue.pop_front();
		int id = node->id * 16;
		float *matPtr = glm::value_ptr(node->matrix);
		for (int i = 0; i < 16; i++) {
			matrix[id + i] = matPtr[i];
		}

		for (int i = 0; i < node->nodes.size(); i++) {
			queue.push_back(node->nodes[i]);
		}
	}
}

void SQMAlgorithm::getTransformMatrices(float* matrix) {
	//add transformation matrix for each node forming skin skeleton
	if (skeleton == NULL || numOfSkinMatrices == 0) return;

	SQMNode *start = root;
	while (root->getSQMNodeType() == SQMCreatedCapsule) {
		//this is worms head should move to original node
		start = (*start->getNodes())[0];
	}

	deque<SkinSkeleton*> queue;
	deque<SQMNode*> shadowQueue;

	queue.push_back(skeleton);
	shadowQueue.push_back(start);

	while (!queue.empty()) {
		SkinSkeleton *skin = queue.front();
		SQMNode *node = shadowQueue.front();
		queue.pop_front();
		shadowQueue.pop_front();

		int id = skin->id * 16;
		glm::mat4 mat = node->getTransformationMatrix();
		float *matPtr = glm::value_ptr(mat);
		for (int i = 0; i < 16; i++) {
			matrix[id + i] = matPtr[i];
		}

		vector<SQMNode*> *childs = node->getNodes();
		for (int i = 0; i < skin->nodes.size(); i++) {
			queue.push_back(skin->nodes[i]);
			shadowQueue.push_back((*childs)[i]);
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

#pragma region Skeleton Operations

SQMNode* SQMAlgorithm::findBNPInTree() {
	deque<SQMNode*> queue;
	if (root != NULL) queue.push_back(root);

	while (!queue.empty()) {
		SQMNode* node = queue.front();
		queue.pop_front();
		int plus = (node->getParent() == NULL) ? 1 : 0;

		if (node->getNodes()->size() >= 2 + plus)
			return node;

		for (int i = 0; i < node->getNodes()->size(); i++) {
			queue.push_back((*node->getNodes())[i]);
		}
	}

	return NULL;
}

void SQMAlgorithm::swapRoot(SQMNode *node) {
	//node is the first that is BNP
	//everything going from root is OK
	//only route from node to root needs to be fixed
	SQMNode *temp = node->getParent();
	SQMNode *parent = temp->getParent();
	node->setParent(NULL);
	temp->setParent(node);
	temp->removeDescendant(node);
	node->addDescendant(temp);

	while (temp != root) {
		temp->addDescendant(parent);
		parent->removeDescendant(temp);
		SQMNode *newParent = parent->getParent();
		parent->setParent(temp);
		temp = parent;
		parent = newParent;
	}
	root = node;
}

void SQMAlgorithm::fixWorm() {
	SQMNode *node = root;
	if (node->getNodes()->size() == 1) {
		//if we got son and parent we need to change root
		if (node->getParent() != NULL) {
			while (node->getParent() != NULL) {
				node = node->getParent();
				root = node;
			}
			return;
		}
		//if next is leaf we need to add extra node
		if ((*node->getNodes())[0]->isLeafNode()) {
			//add extra node
			SQMNode *second = (*node->getNodes())[0];
			SQMNode *temp = new SQMNode(*node);
			temp->setPosition((node->getPosition() + second->getPosition()) / 2.0);
			node->removeDescendant(second);

			node->addDescendant(temp);
			temp->setParent(node);
			second->setParent(temp);
			numOfNodes++;
		}
	} else {//node has 2 childs
		SQMNode *parent = (*node->getNodes())[0];

		while (parent != NULL) {
			node->removeDescendant(parent);
			node->setParent(parent);
			parent->addDescendant(node);
			node = parent;
			//only leaf node wont have 2 childs
			parent = (node->getNodes()->size() == 2) ? (*node->getNodes())[0] : NULL;
		}
		node->setParent(NULL);
		root = node;
	}
}

#pragma endregion

#pragma region SQM Algorithm

void SQMAlgorithm::updateResetRoot() {
	if (resetRoot) delete resetRoot;
	resetRoot = new SQMNode(*root);
}

void SQMAlgorithm::restart() {
	if (root) delete root;
	root = new SQMNode(*resetRoot);
	sqmState = SQMStart;
	numOfSkinMatrices = 0;
}

void SQMAlgorithm::straightenSkeleton() {
	updateResetRoot();
	refreshIDs();
	if (LOG_COMPUTATION_TIME) {
		fb->open("logs/log.txt", ios::out);
		(*os) << "Skeleton straightening\n";
	}
	clock_t ts, te;
	if (skeleton) delete skeleton;
	SkinSkeleton *skeletonB = NULL;
	SkinSkeleton *skeletonR = NULL;
	totalClocks = 0;
	algorithmClocks = 0;

	if (root->getNodes()->size() == 0) {
		//handle just root
		if (mesh) delete mesh;
		mesh = new MyMesh();
		ts = clock();
		root->createScaledIcosahderon(mesh);
		te = clock();
		sqmState = SQMFinalPlacement;
		if (LOG_COMPUTATION_TIME) {
			(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
			(*os) << "\n\nPreprocess time " << 0 << " clicks (" << 0 << " seconds)\n";
			(*os) << "Algorithm time " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
			(*os) << "Total time " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
		}
	} else if (root->getNodes()->size() <= 2) {
		SQMNode *node = findBNPInTree();
		if (node == NULL) {
			//handle worm
			ts = clock();
			root->createCapsules();
			fixWorm();
			numOfNodes = countNodes();
			refreshIDs();
			te = clock();
			if (LOG_COMPUTATION_TIME)
				(*os) << "\tPreprocess took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
			totalClocks += (te - ts);
			if (mesh) delete mesh;
			mesh = new MyMesh();
			skeletonR = root->exportToSkinSkeleton(NULL);
			ts = clock();
			root->wormCreate(mesh);
			te = clock();
			skeletonB = root->exportToSkinSkeleton(NULL);
			if (LOG_COMPUTATION_TIME)
				(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
			algorithmClocks += te - ts;
			sqmState = SQMJoinBNPs;
		} else {
			ts = clock();
			swapRoot(node);
			refreshIDs();
			te = clock();
			if (LOG_COMPUTATION_TIME)
				(*os) << "\tPreprocess took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
			totalClocks += te - ts;
		}
	}
	if (root->getNodes()->size() >= 3) {
		ts = clock();
		root->createCapsules();
		numOfNodes = countNodes();
		refreshIDs();
		te = clock();
		if (LOG_COMPUTATION_TIME)
			(*os) << "\tPreprocess took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
		totalClocks += te - ts;
		//root->straightenSkeleton(OpenMesh::Vec3f(0, 0, 0));
		skeletonR = root->exportToSkinSkeleton(NULL);
		ts = clock();
		root->straightenSkeleton(NULL);
		te = clock();
		skeletonB = root->exportToSkinSkeleton(NULL);
		(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
		algorithmClocks += te - ts;
		sqmState = SQMStraighten;
	}
	ts = clock();
	if (skeletonB != NULL) {
		skeletonR->CalculateCorrespondingDoF(skeletonB);
		skeletonR->ComputeSkinningMatrices();
		//skeletonR->ComputeCompoundRotation();
		skeleton = skeletonR;
	}
	te = clock();
	if (LOG_COMPUTATION_TIME)
		(*os) << "\tSkinning computation took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
	totalClocks += te - ts;

	if (LOG_COMPUTATION_TIME) {
		(*os) << "\n";
		fb->close();
	}
	delete skeletonB;
}

void SQMAlgorithm::computeConvexHull() {
	if (LOG_COMPUTATION_TIME) {
		fb->open("logs/log.txt", ios::app);
		(*os) << "Convex hull computation\n";
	}
	clock_t ts, te;

	ts = clock();
	vector<SQMNode*> stack;
	stack.push_back(root);
	while (!stack.empty()) {
		SQMNode *node = stack.back();
		stack.pop_back();
		if (node->isBranchNode())
			node->calculateConvexHull();
		stack.insert(stack.end(), node->getNodes()->begin(), node->getNodes()->end());
	}
	te = clock();
	if (LOG_COMPUTATION_TIME)
		(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
	algorithmClocks += te - ts;
	drawingMode = 2;
	sqmState = SQMComputeConvexHull;

	if (LOG_COMPUTATION_TIME) {
		(*os) << "\n";
		fb->close();
	}
}

void SQMAlgorithm::subdivideConvexHull() {
	if (LOG_COMPUTATION_TIME) {
		fb->open("logs/log.txt", ios::app);
		(*os) << "Convex hull subdivision\n";
	}
	clock_t ts, te;

	ts = clock();
	root->subdividePolyhedra(NULL, 0, smoothingAlgorithm);
	te = clock();
	if (LOG_COMPUTATION_TIME)
		(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
	algorithmClocks += te - ts;
	sqmState = SQMSubdivideConvexHull;

	if (LOG_COMPUTATION_TIME) {
		(*os) << "\n";
		fb->close();
	}
}

void SQMAlgorithm::joinBNPs() {
	//deletition is messed up :-(
	//maybe delete just faces and let vertices be?
	//filebuf fb;
	//fb.open ("test.txt", ios::out);
	//ostream os(&fb);
	if (LOG_COMPUTATION_TIME) {
		fb->open("logs/log.txt", ios::app);
		(*os) << "BNP joining\n";
	}
	clock_t ts, te;
	//omerr().connect(os);

	if (mesh != NULL) {
		delete mesh;
	}
	mesh = new MyMesh();
	ts = clock();
	vector<MyMesh::VertexHandle> vector;
	root->joinBNPs(mesh, NULL, vector, OpenMesh::Vec3f(0, 0, 0));
	te = clock();
	if (LOG_COMPUTATION_TIME)
		(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
	algorithmClocks += te - ts;
	drawingMode = 3;
	sqmState = SQMJoinBNPs;

	if (LOG_COMPUTATION_TIME) {
		(*os) << "\n";
		fb->close();
	}
}

void SQMAlgorithm::finalVertexPlacement() {
	if (LOG_COMPUTATION_TIME) {
		fb->open("logs/log.txt", ios::app);
		(*os) << "Final vertex placement\n";
	}
	clock_t ts, te;

	ts = clock();
	calculateSkinSkeletonIDs();
	root->setupSkinningMatrixIDs(skeleton);
	te = clock();
	if (LOG_COMPUTATION_TIME)
		(*os) << "\tPreprocessing took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
	totalClocks += te - ts;

	ts = clock();
	root->rotateWithSkeleton(mesh, skeleton);
	te = clock();
	if (LOG_COMPUTATION_TIME)
		(*os) << "\tIt took " << te - ts << " clicks (" << (((float)(te - ts)) / CLOCKS_PER_SEC) << " seconds)\n";
	algorithmClocks += te - ts;
	sqmState = SQMFinalPlacement;

	if (LOG_COMPUTATION_TIME) {
		(*os) << "\n\nPreprocess time " << totalClocks << " clicks (" << (((float)(totalClocks)) / CLOCKS_PER_SEC) << " seconds)\n";
		(*os) << "Algorithm time " << algorithmClocks << " clicks (" << (((float)(algorithmClocks)) / CLOCKS_PER_SEC) << " seconds)\n";
		totalClocks += algorithmClocks;
		(*os) << "Total time " << totalClocks << " clicks (" << (((float)(totalClocks)) / CLOCKS_PER_SEC) << " seconds)\n";
		(*os) << "\n";
		fb->close();
	}
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
