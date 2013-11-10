#include "stdafx.h"
#include "SQMNode.h"
#include <algorithm>
#include <glm\gtc\matrix_transform.hpp>
#include "SphereDelaunay.h"
#include "FloatArithmetic.h"
#include "Laplacian.h"

using namespace MMath;

#define BIAS 0.1

#pragma region Init

SQMNode::SQMNode(void) {
	position = OpenMesh::Vec3f(0, 0, 0);
	originalPosition = position;
	nodeRadius = 10;
	id = 0;
	idStr = "0";
	parent = NULL;
	polyhedron = NULL;
	tessLevel = 1;
	scalev = glm::vec3(1, 1, 1);
	rotatev = glm::vec3();
	transformationMatrix = glm::mat4();
	sqmNodeType = SQMNone;
	skinningIDs = glm::ivec2(-1, -1);
}

SQMNode::SQMNode(SQMSkeletonNode &node, SQMNode* newParent) : parent(newParent) {
	if (parent == NULL) {
		idStr = "0";
	} else {
		idStr = parent->getIdStr();
		string s = ToStringT(parent->getNumOfChilds());
		idStr += "-" + s;
	}
	id = 0;
	position = OpenMesh::Vec3f(node.point.x, node.point.y, node.point.z);
	originalPosition = position;
	nodeRadius = node.radius;
	tessLevel = node.tessLevel;
	sqmNodeType = node.capsule ? SQMCapsule : SQMNone;
	for (int i = 0; i < node.nodes.size(); i++) {
		SQMNode *newNode = new SQMNode(*node.nodes[i], this);
		nodes.push_back(newNode);
	}
	polyhedron = NULL;
	scalev = glm::vec3(node.scale.x, node.scale.y, node.scale.z);
	rotatev = glm::vec3(node.rotate.x, node.rotate.y, node.rotate.z);
	transformationMatrix = glm::mat4();
	skinningIDs = glm::ivec2(-1, -1);
}

SQMNode::SQMNode(SQMNode &node) {
	parent = NULL;
	polyhedron = NULL;
	id = node.getId();
	idStr = node.getIdStr();
	sqmNodeType = node.getSQMNodeType();
	nodeRadius = node.getNodeRadius();
	tessLevel = node.getTessLevel();
	position = node.getPosition();
	originalPosition = position;
	axisAngle = node.getAxisAngle();
	vector<SQMNode *> *childs = node.getNodes();
	scalev = node.getScalev();
	rotatev = node.getRotatev();
	transformationMatrix = node.getTransformationMatrix();
	for (int i = 0; i < childs->size(); i++) {
		SQMNode *childRef = (*childs)[i];
		SQMNode *child = new SQMNode(*childRef);
		child->setParent(this);
		nodes.push_back(child);
	}
}

SQMNode::~SQMNode(void) {
	for (int i = 0; i < nodes.size(); i++) {
		delete nodes[i];
	}
	if (polyhedron != NULL) delete polyhedron;
}

#pragma endregion

#pragma region Getters

unsigned int SQMNode::getId() {
	return id;
}

string SQMNode::getIdStr() {
	return idStr;
}

SQMNodeType SQMNode::getSQMNodeType() {
	return sqmNodeType;
}

bool SQMNode::isBranchNode() {
	int requiredConections = 3;//branch nodes has at least 3 connections
	if (parent) requiredConections--; //parent counts as one conection

	return nodes.size() >= requiredConections;
}

bool SQMNode::isConnectionNode() {
	return (parent != NULL && nodes.size() == 1);
}

bool SQMNode::isLeafNode() {
	return (nodes.size() == 0) || (nodes.size() == 1 && parent == NULL);
}

OpenMesh::Vec3f SQMNode::getPosition() {
	return position;
}

OpenMesh::Vec3f SQMNode::getOldPosition() {
	return oldPosition;
}

OpenMesh::Vec3f SQMNode::getOriginalPosition() {
	return originalPosition;
}

glm::vec3 SQMNode::getPosition_glm() {
	return glm::vec3(position[0], position[1], position[2]);
}

vector<SQMNode*>* SQMNode::getNodes() {
	return &nodes;
}

MyTriMesh* SQMNode::getPolyhedron() {
	return polyhedron;
}

vector<SQMNode*>* SQMNode::getDescendants() {
	return &nodes;
}

vector<MyTriMesh::VertexHandle>* SQMNode::getIntersectionVHandles() {
	return &intersectionVHandles;
}

float SQMNode::getNodeRadius() {
	return nodeRadius;
}

float SQMNode::getTessLevel() {
	return tessLevel;
}

Quaternion SQMNode::getAxisAngle() {
	return axisAngle;
}

SQMNode* SQMNode::getParent() {
	return parent;
}

int SQMNode::getNumOfChilds() {
	return nodes.size();
}

float SQMNode::getX() {
	return position[0];
}

float SQMNode::getY() {
	return position[1];
}

float SQMNode::getZ() {
	return position[2];
}

glm::mat4 SQMNode::getTransformationMatrix() {
	return transformationMatrix;
}

glm::vec3 SQMNode::getScalev() {
	return scalev;
}

glm::vec3 SQMNode::getRotatev() {
	return rotatev;
}

float SQMNode::getScaleX() {
	return scalev.x;
}

float SQMNode::getScaleY() {
	return scalev.y;
}

float SQMNode::getScaleZ() {
	return scalev.z;
}

float SQMNode::getRotateX() {
	return rotatev.x;
}

float SQMNode::getRotateY() {
	return rotatev.y;
}

float SQMNode::getRotateZ() {
	return rotatev.z;
}

#pragma endregion

#pragma region Setters

void SQMNode::setID(unsigned int newID) {
	id = newID;
}

void SQMNode::setParent(SQMNode *node) {
	parent = node;
}

void SQMNode::setNodeRadius(float newNodeRadius) {
	nodeRadius = newNodeRadius;
}

void SQMNode::setTessLevel(float newTessLevel) {
	tessLevel = newTessLevel;
}

void SQMNode::setPosition(OpenMesh::Vec3f newPosition) {
	position = newPosition;
	originalPosition = position;
}

void SQMNode::setSQMNodeType(SQMNodeType newType) {
	sqmNodeType = newType;
}

void SQMNode::setPosition(float x, float y, float z) {
	position = OpenMesh::Vec3f(x, y, z);
	originalPosition = position;
}

void SQMNode::addDescendant(SQMNode* node) {
	nodes.push_back(node);
}

void SQMNode::removeDescendant(SQMNode* node) {
	vector<SQMNode *> temp;
	for (vector<SQMNode *>::iterator it = nodes.begin(); it != nodes.end(); it++) {
		if ((*it) == node) {
			it = nodes.erase(it);
			return;
		}
	}
}

void SQMNode::removeDescendants() {
	nodes.clear();
}

void SQMNode::rotatePosition(Quaternion q, CVector3 offset) {
	CVector3 pos(position.values_);
	pos = pos - offset;
	pos = QuaternionRotateVector(q, pos);
	pos = pos + offset;
	position = OpenMesh::Vec3f(pos.x, pos.y, pos.z);
}

void SQMNode::addDescendant(float x, float y, float z) {
	SQMSkeletonNode *skeletonNode = new SQMSkeletonNode(x, y, z);
	SQMNode *node = new SQMNode(*skeletonNode, this);
	addDescendant(node);
	delete skeletonNode;
}

void SQMNode::setX(float newX) {
	position[0] = newX;
}

void SQMNode::setY(float newY) {
	position[1] = newY;
}

void SQMNode::setZ(float newZ) {
	position[2] = newZ;
}

void SQMNode::setTransformationMatrix(glm::mat4 tm) {
	transformationMatrix = tm;
}

void SQMNode::setScaleX(float value) {
	scalev.x = value;
	updateTransformationMatrix();
}

void SQMNode::setScaleY(float value) {
	scalev.y = value;
	updateTransformationMatrix();
}

void SQMNode::setScaleZ(float value) {
	scalev.z = value;
	updateTransformationMatrix();
}

void SQMNode::setRotateX(float value) {
	rotatev.x = value;
	updateTransformationMatrix();
}

void SQMNode::setRotateY(float value) {
	rotatev.y = value;
	updateTransformationMatrix();
}

void SQMNode::setRotateZ(float value) {
	rotatev.z = value;
	updateTransformationMatrix();
}

void SQMNode::updateTransformationMatrix() {
	glm::mat4 M = glm::mat4();
	glm::vec3 translatev = glm::vec3(position[0], position[1], position[2]);
	M = glm::translate(M, translatev);
	M = glm::rotate(M, rotatev.x, glm::vec3(1, 0, 0));
	M = glm::rotate(M, rotatev.y, glm::vec3(0, 1, 0));
	M = glm::rotate(M, rotatev.z, glm::vec3(0, 0, 1));
	M = glm::scale(M, scalev);
	M = glm::translate(M, -translatev);

	transformationMatrix = M;
}

void SQMNode::addVHandleToRotate(MyMesh::VHandle vh) {
	meshVhandlesToRotate.push_back(vh);
}

#pragma endregion

#pragma region Export

SQMSkeletonNode* SQMNode::exportToSkeletonNode() {
	SQMSkeletonNode* node = new SQMSkeletonNode(position[0], position[1], position[2], nodeRadius);
	node->tessLevel = tessLevel;
	node->id = id;
	node->capsule = (sqmNodeType == SQMCapsule);
	node->scale = CVector3(scalev.x, scalev.y, scalev.z);
	node->rotate = CVector3(rotatev.x, rotatev.y, rotatev.z);
	for (int i = 0; i < nodes.size(); i++) {
		node->addChild(nodes[i]->exportToSkeletonNode());
	}
	return node;
}

SkinSkeleton* SQMNode::exportToSkinSkeleton(SkinSkeleton *parentSkin) {
	//if this is worms head
	if (sqmNodeType == SQMCreatedCapsule && nodes.size() > 0) return nodes[0]->exportToSkinSkeleton(parentSkin);

	SkinSkeleton *node = new SkinSkeleton(parentSkin, position[0], position[1], position[2]);
	//if next is only capsule the matrix would be the same
	if (sqmNodeType == SQMFormerCapsule && !(parent != NULL && parent->getSQMNodeType() == SQMCreatedCapsule)) return node;

	for (int i = 0; i < nodes.size(); i++) {
		node->nodes.push_back(nodes[i]->exportToSkinSkeleton(node));
	}
	return node;
}

#pragma endregion

#pragma region SQM Preprocessing

void SQMNode::createCapsules(int minSmallCircles) {
	if (this->isLeafNode() && (sqmNodeType == SQMCapsule)) {
		bool nullParent = (parent == NULL);
		//create vector
		OpenMesh::Vec3f dir;
		if (nullParent) {
			dir = (position - nodes[0]->getPosition()).normalize();
		} else {
			dir = (position - parent->getPosition()).normalize();
		}
		//get number
		int smallCircles = max(nodeRadius, (float)minSmallCircles);
		SQMNode *current = this;
		for (int i = 1; i <= smallCircles; i++) {
			//calculate radius
			glm::vec2 P1(0.23, 1);
			glm::vec2 P2(0.32, 1);
			float t = (float)i/(float)smallCircles;
			float step = bezier(P1, P2, t).y;
			float newRadius = sqrtf(pow(nodeRadius, 2) * (1 - pow(step, 2)));
			if (i == smallCircles) newRadius = 1;
			//calculate new position
			OpenMesh::Vec3f newPosition = position + (step*nodeRadius*dir);
			//insert node with calculated radius
			SQMNode *newNode = new SQMNode(*current);
			newNode->setNodeRadius(newRadius);
			newNode->setPosition(newPosition);
			newNode->setSQMNodeType(SQMCreatedCapsule);
			if (nullParent) {
				newNode->removeDescendants();
				newNode->addDescendant(current);
				current->setParent(newNode);
			} else {
				newNode->setParent(current);
				current->addDescendant(newNode);
			}

			current = newNode;
		}
		sqmNodeType = SQMFormerCapsule;
	}
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->createCapsules();
	}
}

#pragma endregion

#pragma region Skeleton Straightening

void SQMNode::straightenSkeleton(OpenMesh::Vec3f *lineVector) {
	axisAngle = Quaternion();
	oldPosition = position;
	if (lineVector != NULL && !parent->isBranchNode()) {//straighten self
		SQMNode *ancestor = getAncestorBranchNode(this);
		if (ancestor != NULL) {
			rotatePosition(QuaternionOpposite(ancestor->getAxisAngle()), CVector3(ancestor->getParent()->getPosition().values_));
		}
		//translate parent to 0,0,0
		OpenMesh::Vec3f newPosition = position - parent->getPosition();
		CVector3 oldPos(newPosition.values_);
		//roatate
		float len = (originalPosition - parent->getOriginalPosition()).length();
		newPosition = len*(*lineVector);
		CVector3 newPos(newPosition.values_);
		Quaternion quaternion = SQMQuaternionBetweenVectors(newPos, oldPos);
		//translate back by parent
		newPosition = newPosition + parent->getPosition();
		//setup
		axisAngle = quaternion;
		position = newPosition;
	}
	if (this->isBranchNode()) {//if this is branch node recalculate new vectors and intersections
		for (int i = 0; i < nodes.size(); i++) {//specifical order parent intersection needs to be last in vector
			if (parent != NULL) {
				nodes[i]->rotateSelfAndDescendants(QuaternionOpposite(axisAngle), CVector3(parent->getPosition().values_));
			}
			OpenMesh::Vec3f u = nodes[i]->getPosition() - position;
			u = u.normalize();
			nodes[i]->straightenSkeleton(&u);
			intersections.push_back(position + u * nodeRadius);
		}
		if (parent) {//also calculate intersection with parent
			OpenMesh::Vec3f u = parent->getPosition() - position;
			u = u.normalize();
			intersections.push_back(position + u * nodeRadius);
		}
	} else {//else just straighten conected nodes
		for (int i = 0; i < nodes.size(); i++) {
			nodes[i]->straightenSkeleton(lineVector);
		}
	}
}

#pragma endregion

#pragma region BNP generation

void SQMNode::calculateConvexHull() {
	vector<OpenMesh::Vec3i> triangles;
	vector<OpenMesh::Vec3f> translatedIntersections(intersections.size());
	//in triangularization the points are normalized
	centerOfMass = OpenMesh::Vec3f(0, 0, 0);
	for (int i = 0; i < intersections.size(); i++) {
		translatedIntersections[i] = intersections[i] - position;
		centerOfMass += intersections[i];
	}
	centerOfMass /= intersections.size();
	Delaunay_on_sphere(translatedIntersections, triangles);
	triangles2 = triangles;
	createPolyhedra(triangles);
}

void SQMNode::createPolyhedra(vector<OpenMesh::Vec3i> triangles) {
	//calculate triangle normals
	vector<OpenMesh::Vec3f> normals;
	vector<OpenMesh::Vec3f> centers;
	for (int i = 0; i < triangles.size(); i++) {
		OpenMesh::Vec3i triangle = triangles[i];
		OpenMesh::Vec3f A = intersections[triangle.values_[0]];
		OpenMesh::Vec3f B = intersections[triangle.values_[1]];
		OpenMesh::Vec3f C = intersections[triangle.values_[2]];
		OpenMesh::Vec3f normal = cross(B - A, C - A);
		normal = normal.normalize();
		OpenMesh::Vec3f center((A[0] + B[0] + C[0])/3.0, (A[1] + B[1] + C[1])/3.0, (A[2] + B[2] + C[2])/3.0);

		normals.push_back(normal);		
		centers.push_back(center);
	}
	normals2 = normals;
	centers2 = centers;
	map<OpenMesh::Vec2i, vector<int>, OpenMeshVec2iComp> edgeFaceIndexMap;
	for (int i = 0; i < triangles.size(); i++) {
		OpenMesh::Vec3i triangle = triangles[i];
		OpenMesh::Vec2i u1(triangle.values_[0], triangle.values_[1]);
		OpenMesh::Vec2i u2(triangle.values_[1], triangle.values_[2]);
		OpenMesh::Vec2i u3(triangle.values_[2], triangle.values_[0]);
		vector<int> faces1;
		vector<int> faces2;
		vector<int> faces3;
		for (int j = 0; j < triangles.size(); j++) {
			OpenMesh::Vec3i triangle2 = triangles[j];
			OpenMesh::Vec2i v1(triangle2.values_[0], triangle2.values_[1]);
			OpenMesh::Vec2i v2(triangle2.values_[1], triangle2.values_[2]);
			OpenMesh::Vec2i v3(triangle2.values_[2], triangle2.values_[0]);
			if (OpenMeshVec2iEqual(u1, v1) || OpenMeshVec2iEqual(u1, v2) || OpenMeshVec2iEqual(u1, v3)) {
				faces1.push_back(j);
			}
			if (OpenMeshVec2iEqual(u2, v1) || OpenMeshVec2iEqual(u2, v2) || OpenMeshVec2iEqual(u2, v3)) {
				faces2.push_back(j);
			}
			if (OpenMeshVec2iEqual(u3, v1) || OpenMeshVec2iEqual(u3, v2) || OpenMeshVec2iEqual(u3, v3)) {
				faces3.push_back(j);
			}
		}
		edgeFaceIndexMap.insert(pair<OpenMesh::Vec2i, vector<int> >(u1, faces1));
		edgeFaceIndexMap.insert(pair<OpenMesh::Vec2i, vector<int> >(u2, faces2));
		edgeFaceIndexMap.insert(pair<OpenMesh::Vec2i, vector<int> >(u3, faces3));
	}
	//for each triangle create 6 new triangles and translate new vertices and translate vertices with face normals add only unique
	vector<OpenMesh::Vec3f> vertices;
	vector<OpenMesh::Vec3i> faces;
	vector<OpenMesh::Vec2i> visited;
	for (int i = 0; i < triangles.size(); i++) {
		//get triangle points and create new ones
		OpenMesh::Vec3i triangle = triangles[i];
		OpenMesh::Vec2i e1(triangle[0], triangle[1]);
		OpenMesh::Vec2i e2(triangle[1], triangle[2]);
		OpenMesh::Vec2i e3(triangle[2], triangle[0]);

		int v1Index = getPointPositionInArray(OpenMesh::Vec2i(triangle[0], triangle[0]), visited);
		int v2Index = getPointPositionInArray(OpenMesh::Vec2i(triangle[1], triangle[1]), visited);
		int v3Index = getPointPositionInArray(OpenMesh::Vec2i(triangle[2], triangle[2]), visited);

		int u12Index = getPointPositionInArray(e1, visited);
		int u23Index = getPointPositionInArray(e2, visited);
		int u31Index = getPointPositionInArray(e3, visited);

		OpenMesh::Vec3f v1 = intersections[triangle[0]];
		OpenMesh::Vec3f v2 = intersections[triangle[1]];
		OpenMesh::Vec3f v3 = intersections[triangle[2]];

		if (v1Index == -1) {
			v1Index = vertices.size();
			vertices.push_back(v1);
			visited.push_back(OpenMesh::Vec2i(triangle[0], triangle[0]));
		}
		if (v2Index == -1) {
			v2Index = vertices.size();
			vertices.push_back(v2);
			visited.push_back(OpenMesh::Vec2i(triangle[1], triangle[1]));
		}
		if (v3Index == -1) {
			v3Index = vertices.size();
			vertices.push_back(v3);
			visited.push_back(OpenMesh::Vec2i(triangle[2], triangle[2]));
		}
		if (u12Index == -1) {
			OpenMesh::Vec3f u12 = 0.5*v1 + 0.5*v2;
			//translate point
			vector<int> u12NormalIndexis = getNormalIndexis(edgeFaceIndexMap[OpenMesh::Vec2i(triangle.values_[0], triangle.values_[1])], i);
			int i1 = u12NormalIndexis[0];
			int i2 = u12NormalIndexis[1];
			u12 = translatedPointToSphereWithFaceNormals(u12, normals[i1], normals[i2], centers[i1], centers[i2]);

			u12Index = vertices.size();
			vertices.push_back(u12);
			visited.push_back(e1);
		}
		if (u23Index == -1) {
			OpenMesh::Vec3f u23 = 0.5*v2 + 0.5*v3;
			//translate point
			vector<int> u23NormalIndexis = getNormalIndexis(edgeFaceIndexMap[OpenMesh::Vec2i(triangle.values_[1], triangle.values_[2])], i);
			int i1 = u23NormalIndexis[0];
			int i2 = u23NormalIndexis[1];
			u23 = translatedPointToSphereWithFaceNormals(u23, normals[i1], normals[i2], centers[i1], centers[i2]);

			u23Index = vertices.size();
			vertices.push_back(u23);
			visited.push_back(e2);
		}
		if (u31Index == -1) {
			OpenMesh::Vec3f u31 = 0.5*v3 + 0.5*v1;
			//translate point
			vector<int> u31NormalIndexis = getNormalIndexis(edgeFaceIndexMap[OpenMesh::Vec2i(triangle.values_[2], triangle.values_[0])], i);
			int i1 = u31NormalIndexis[0];
			int i2 = u31NormalIndexis[1];
			u31 = translatedPointToSphereWithFaceNormals(u31, normals[i1], normals[i2], centers[i1], centers[i2]);

			u31Index = vertices.size();
			vertices.push_back(u31);
			visited.push_back(e3);
		}		
		OpenMesh::Vec3f center = centers[i];
		center = translatedPointToSphereWithFaceNormals(center, normals[i], normals[i], center, center);
		//need something to hold the same size
		int centerIndex = vertices.size();
		visited.push_back(OpenMesh::Vec2i(-1, -1));
		vertices.push_back(center);
		//add new triangles to face list
		faces.push_back(OpenMesh::Vec3i(v1Index, u12Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(u12Index, v2Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(v2Index, u23Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(u23Index, v3Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(v3Index, u31Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(u31Index, v1Index, centerIndex));
	}
	//create OpenMesh mesh from indexed face
	openMeshFromIdexedFace(vertices, faces);
}

OpenMesh::Vec3f SQMNode::translatedPointToSphereWithFaceNormals(OpenMesh::Vec3f p, OpenMesh::Vec3f n1, OpenMesh::Vec3f n2, OpenMesh::Vec3f center1, OpenMesh::Vec3f center2) {
	//calculate leading vector
	OpenMesh::Vec3f u = (n1 + n2);//B-A
	if (u.norm() < FLOAT_ZERO) {
		OpenMesh::Vec3f center = 0.5*center1 + 0.5*center2;
		u = (p - center).normalize();
	} else {
		u = u.normalize();
	}
	float t;
	if (raySphereIntersection(p, u, position, nodeRadius, t)) {
		return p + (t*u);
	} else {
		return p;
	}
}

vector<int> SQMNode::getNormalIndexis(vector<int> indexis, int index) {
	//if there are two indexis
	if (indexis.size() == 2)
		return indexis;
	//if there is only one index (error?)
	if (indexis.size() == 1) {
		vector<int> result;
		result.push_back(indexis[0]);
		result.push_back(indexis[0]);
		return result;
	}
	//if there is none
	if (indexis.size() == 0) {
		throw exception("missing normal vector indexis");
	}
	//if there are many indexis the normal is equal to face normal
	vector<int> result;
	result.push_back(index);
	result.push_back(index);
	return result;
}

void SQMNode::openMeshFromIdexedFace(vector<OpenMesh::Vec3f> vertices, vector<OpenMesh::Vec3i> faces) {
	if (polyhedron != NULL) delete polyhedron;
	polyhedron = new MyTriMesh();
	polyhedronPoints.clear();
	intersectionVHandles.clear();
	intersectionVHandles.resize(intersections.size());

	for (int i = 0; i < vertices.size(); i++) {
		polyhedronPoints.push_back(polyhedron->add_vertex(MyMesh::Point(vertices[i].values_)));
		//store intersections vertex handles
		int position = getPointPositionInArray(vertices[i], intersections);//specific order equal to intersection vector
		if (position != -1) {
			intersectionVHandles[position] = polyhedronPoints.back();
		}
	}
	for (int i = 0; i < faces.size(); i++) {
		OpenMesh::Vec3i face = faces[i];
		polyhedron->add_face(polyhedronPoints[face.values_[0]], polyhedronPoints[face.values_[1]], polyhedronPoints[face.values_[2]]);
	}
	//writeTriMesh(polyhedron);
}

//-----------------------TEST------------------------

void SQMNode::createPolyhedraFromCenter(vector<OpenMesh::Vec3i> triangles) {
	//get the mesh center
	OpenMesh::Vec3f meshCenter = polyhedronBoundingBoxCenter();
	//get the center of each triangle
	vector<OpenMesh::Vec3f> centers;
	for (int i = 0; i < triangles.size(); i++) {
		OpenMesh::Vec3i triangle = triangles[i];
		OpenMesh::Vec3f A = intersections[triangle.values_[0]];
		OpenMesh::Vec3f B = intersections[triangle.values_[1]];
		OpenMesh::Vec3f C = intersections[triangle.values_[2]];
		OpenMesh::Vec3f center((A[0] + B[0] + C[0])/3.0, (A[1] + B[1] + C[1])/3.0, (A[2] + B[2] + C[2])/3.0);	
		centers.push_back(center);
	}
	//for each triangle create 6 new triangles and translate new vertices and translate vertices with face normals add only unique
	vector<OpenMesh::Vec3f> vertices;
	vector<OpenMesh::Vec3i> faces;
	for (int i = 0; i < triangles.size(); i++) {
		//get triangle points and create new ones
		OpenMesh::Vec3i triangle = triangles[i];
		OpenMesh::Vec3f v1 = intersections[triangle.values_[0]];
		OpenMesh::Vec3f v2 = intersections[triangle.values_[1]];
		OpenMesh::Vec3f v3 = intersections[triangle.values_[2]];

		OpenMesh::Vec3f u12 = 0.5*v1 + 0.5*v2;
		OpenMesh::Vec3f u23 = 0.5*v2 + 0.5*v3;
		OpenMesh::Vec3f u31 = 0.5*v3 + 0.5*v1;
		OpenMesh::Vec3f center = centers[i];
		//translate points
		u12 = translatePointToSphereFromCenter(u12, meshCenter);

		u23 = translatePointToSphereFromCenter(u23, meshCenter);

		u31 = translatePointToSphereFromCenter(u31, meshCenter);

		center = translatePointToSphereFromCenter(center, meshCenter);
		//add only unique points to vertex list
		int v1Index = getPointPositionInArrayOrAdd(v1, vertices);
		int v2Index = getPointPositionInArrayOrAdd(v2, vertices);
		int v3Index = getPointPositionInArrayOrAdd(v3, vertices);

		int u12Index = getPointPositionInArrayOrAdd(u12, vertices);
		int u23Index = getPointPositionInArrayOrAdd(u23, vertices);
		int u31Index = getPointPositionInArrayOrAdd(u31, vertices);

		int centerIndex = getPointPositionInArrayOrAdd(center, vertices);
		//add new triangles to face list
		faces.push_back(OpenMesh::Vec3i(v1Index, u12Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(u12Index, v2Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(v2Index, u23Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(u23Index, v3Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(v3Index, u31Index, centerIndex));
		faces.push_back(OpenMesh::Vec3i(u31Index, v1Index, centerIndex));
		//faces.push_back(OpenMesh::Vec3i(v1Index, v2Index, v3Index));
	}
	//create OpenMesh mesh from indexed face
	openMeshFromIdexedFace(vertices, faces);
}

OpenMesh::Vec3f SQMNode::polyhedronBoundingBoxCenter() {
	float minX = intersections[0][0], minY = intersections[0][1], minZ = intersections[0][2];
	float maxX = minX, maxY = minY, maxZ = minZ;
	for (int i = 0; i < intersections.size(); i++) {
		OpenMesh::Vec3f intersection = intersections[i];

		minX = min(minX, intersection[0]);
		minY = min(minY, intersection[1]);
		minZ = min(minZ, intersection[2]);

		maxX = max(maxX, intersection[0]);
		maxY = max(maxY, intersection[1]);
		maxZ = max(maxZ, intersection[2]);
	}
	return OpenMesh::Vec3f((minX + maxX) / 2.0, (minY + maxY) / 2.0, (minZ + maxZ) / 2.0);
}

OpenMesh::Vec3f SQMNode::polyhedronPointSumCenterCenter() {
	OpenMesh::Vec3f result(0, 0, 0);

	for (int i = 0; i < intersections.size(); i++) {
		result += intersections[i];
	}
	result /= intersections.size();

	return result;
}

OpenMesh::Vec3f SQMNode::translatePointToSphereFromCenter(OpenMesh::Vec3f point, OpenMesh::Vec3f center) {
	OpenMesh::Vec3f direction = (point - center).normalize();
	float t = 0;
	if (raySphereIntersection(center, direction, position, nodeRadius, t)) {
		return center + (direction * t);
	} else {
		return point;
	}
}

#pragma endregion

#pragma region BNP Subdivision

void SQMNode::subdividePolyhedra(SQMNode* parentBranchNode, int count, SQMSmoothingAlgorithm algorithm) {
	vector<SQMNode*> branchingNodes;
	for (int i = 0; i < nodes.size(); i++) {
		branchingNodes.push_back(getDescendantBranchNode(nodes[i]));
	}
	map<int, LIENeedEntry> lieMap;
	fillLIEMap(count, lieMap, branchingNodes);
	if (algorithm == SQMQuaternionSmoothing) {
		//smoothLIEs(lieMap);
	}
	splitLIEs(lieMap, algorithm);
	//take care of the rest
	for (int i = 0; i < branchingNodes.size(); i++)	{
		SQMNode* branchingNode = branchingNodes[i];
		if (branchingNode != NULL) {
			MyTriMesh::VHandle vhandle = intersectionVHandles[i];
			//if we needed some they have been added if we had more points than descendant he needs to add them
			int missingPoints = verticeDifferenceFatherSon(this, branchingNode, vhandle);
			missingPoints = (missingPoints < 0) ? -missingPoints : 0;
			branchingNode->subdividePolyhedra(this, missingPoints, algorithm);
		}
	}
}

void SQMNode::fillLIEMap(int parentNeed, std::map<int, LIENeedEntry>& lieMap, std::vector<SQMNode*>& branchingNodes) {
	int stop = nodes.size();
	if (parent != NULL)
		stop++;
	for (int i = 0; i < stop; i++) {
		//count of need to split edges
		MyTriMesh::VHandle vhandle = intersectionVHandles[i];
		int need = 0;
		if (i == nodes.size()) {//parent intersection
			if (parentNeed > 0) {
				need = parentNeed;
			} else {
				need = NO_SPLITTING;
			}
		} else {
			SQMNode* branchingNode = branchingNodes[i];
			need = verticeDifferenceFatherSon(this, branchingNode, vhandle);
			//if father has more vertices than its son there is no need to split
			need = max(0, need);
			//if son is not branch node he can be splited as much as needed
			if (branchingNode == NULL)
				need = UNLIMITED_SPLITTING;
		}

		//we could check for uncreated LIEs only but the list could be long and thus take more time that creating new one
		vector<LIE> verticeLIEs;
		//get first halfedge
		MyTriMesh::HHandle heh = polyhedron->voh_begin(vhandle).current_halfedge_handle();
		heh = polyhedron->next_halfedge_handle(heh);
		//get first opposing vhandle
		MyTriMesh::VHandle ovhandle = oppositeVHandle(heh);
		//check previous for first ocurence and remember first heh
		MyTriMesh::HHandle prevHeh = prevLink(heh);
		while (oppositeVHandle(prevHeh) == ovhandle) {
			heh = prevHeh;
			prevHeh = prevLink(heh);
		}
		//current half edge and opposite handle
		MyTriMesh::HHandle cheh = nextLink(heh);
		MyTriMesh::VHandle covh = ovhandle;
		vector<MyTriMesh::EdgeHandle> edges;
		edges.push_back(polyhedron->edge_handle(heh));
		//for quaternion
		MyTriMesh::VHandle firstLieVertex = polyhedron->to_vertex_handle(polyhedron->opposite_halfedge_handle(heh));
		MyTriMesh::VHandle secondLieVertex = polyhedron->to_vertex_handle(heh);
		MyTriMesh::VHandle lastLieVertex = polyhedron->to_vertex_handle(heh);
		MyTriMesh::HHandle firstLieHEdge = heh;
		MyTriMesh::HHandle lastLieHEdge = heh;
		//CVector3 offset = CVector3(position.values_);
		CVector3 offset = CVector3(centerOfMass.values_);

		bool good = true;
		while (good) {
			good = (heh != cheh);
			//collect all halfedges of LIE
			ovhandle = oppositeVHandle(cheh);
			if (ovhandle == covh) { //if new collect all edges incident with oposing vhandle else just pass them
				edges.push_back(polyhedron->edge_handle(cheh));
				if (firstLieVertex.idx() == -1) {
					firstLieVertex = polyhedron->to_vertex_handle(polyhedron->opposite_halfedge_handle(cheh));
					secondLieVertex = polyhedron->to_vertex_handle(cheh);
					firstLieHEdge = cheh;
				}
				lastLieVertex = polyhedron->to_vertex_handle(cheh);
				lastLieHEdge = cheh;

				cheh = nextLink(cheh);
			} else {//add LIE to map and table of LIEs
				LIE lie(vhandle, covh);
				lie.edges = edges;
				lie.vertice1 = i;
				lie.vertice2 = getPositionInArray(covh, intersectionVHandles);
				lie.firstHHandle = firstLieHEdge;
				lie.lastHHandle = lastLieHEdge;
				CVector3 start = CVector3(polyhedron->point(firstLieVertex).values_) - offset;
				CVector3 dest = CVector3(polyhedron->point(lastLieVertex).values_) - offset;
				CVector3 P = CVector3(polyhedron->point(secondLieVertex).values_) - offset;
				//CVector3 axis = Normalize(Cross(P - start, dest - start));
				//default axis between start and P because they cannot be colinear and fit the one-ring
				CVector3 axis;
				float angle = Dot(Normalize(start), Normalize(dest));
				if (equal(fabs(angle), 1)) {
					axis = Normalize(Cross(start, P));
				} else {
					axis = Normalize(Cross(start, dest));
				}
				//CVector3 axis = Normalize(Cross(start, P));
				lie.quaternion = SQMQuaternionBetweenVectorsWithAxis(start, dest, axis);
				verticeLIEs.push_back(lie);

				//clearing
				firstLieVertex = MyTriMesh::VHandle(-1);
				edges.clear();
				covh = ovhandle;
			}
		}

		//crate map item
		LIENeedEntry entry(vhandle, need);
		entry.lies = verticeLIEs;
		lieMap.insert(pair<int, LIENeedEntry>(i, entry));
	}
}

void SQMNode::splitLIEs(std::map<int, LIENeedEntry>& lieMap, SQMSmoothingAlgorithm algorithm) {
	bool parentFirst = (parent != NULL);
	for (int i = 0; i < nodes.size(); i++) {
		//parent should go first
		if (parentFirst) {
			i = nodes.size();
		}
		//get entry
		LIENeedEntry entry = lieMap.at(i);
		//if  need <= 0 skip
		while (entry.need > 0) {
			//else get LIE with greatest need
			LIE bestLIE;
			int bestNeed;
			int lieIndex = -1;
			for (int j = 0; j < entry.lies.size(); j++) {
				LIE lie = entry.lies[j];
				int otherIndex = lie.otherVerticeIndex(i);
				int lieNeed = lieMap.at(otherIndex).need;

				if (lieNeed == NO_SPLITTING) continue;

				if ((lieIndex == -1) || (lieNeed != 0 && bestNeed < lieNeed) || (bestNeed == 0 && lieNeed == UNLIMITED_SPLITTING) || (bestNeed == lieNeed && lie.refined < bestLIE.refined)) {
					bestNeed = lieNeed;
					bestLIE = lie;
					lieIndex = j;
				}
			}
			//split edge and adjust LIEs
			splitLIE(bestLIE, lieMap, i, lieIndex, algorithm);
			entry = lieMap.at(i);
		}
		if (parentFirst) {
			entry.need = -2;
			lieMap.at(i) = entry;
			i = 0;
			parentFirst = false;
		}
	}
}

void SQMNode::splitLIE(LIE lie, std::map<int, LIENeedEntry>& lieMap, int entryIndex, int lieIndex, SQMSmoothingAlgorithm algorithm) {
	LIE newLie = splitLIEEdge(lie);
	if (algorithm == SQMQuaternionSmoothing) {
		smoothLIE(newLie);
	}
	if (algorithm == SQMAvaragingSmoothing) {
		smoothLIEByAvaraging(newLie);
	}
	//decrease need for both vertices
	LIENeedEntry entry1 = lieMap.at(entryIndex);
	LIENeedEntry entry2 = lieMap.at(lie.otherVerticeIndex(entryIndex));
	//we always split first 
	entry1.lies[lieIndex] = newLie;
	int otherIndex = getPositionInArray(newLie, entry2.lies);
	entry2.lies[otherIndex] = newLie;
	//only if they are bigger than one
	if (entry1.need > 0)
		entry1.need--;
	if (entry2.need > 0)
		entry2.need--;
	lieMap.at(lie.vertice1) = entry1;
	lieMap.at(lie.vertice2) = entry2;
	if (algorithm == SQMOneRingLaplacianSmoothing || algorithm == SQMValencyLaplacianSmoothing) {
		laplacianSMoothing(algorithm);
	}
}

LIE SQMNode::splitLIEEdge(LIE lie) {
	MyTriMesh::EHandle eh = lie.edges[0];
	MyTriMesh::EHandle newEh = splitEdgeInHalfAndReturnNewEdge(eh);
	lie.edges.insert(lie.edges.begin()+1, newEh);
	lie.refined = lie.refined + 1;
	return lie;
}

MyTriMesh::EHandle SQMNode::splitEdgeInHalfAndReturnNewEdge(MyTriMesh::EdgeHandle eh) {
	MyTriMesh::HalfedgeHandle heh0 = polyhedron->halfedge_handle(eh, 0);
	MyTriMesh::HalfedgeHandle heh1 = polyhedron->halfedge_handle(eh, 1);

	MyTriMesh::Point p0 = polyhedron->point(polyhedron->to_vertex_handle(heh0));
	MyTriMesh::Point p1 = polyhedron->point(polyhedron->to_vertex_handle(heh1));

	MyTriMesh::Point x = 0.5*p0 + 0.5*p1;
	MyTriMesh::VertexHandle vh = polyhedron->add_vertex(x);

	polyhedron->split(eh, vh);

	MyTriMesh::VHandle v2 = polyhedron->to_vertex_handle(heh1);
	for (MyTriMesh::VOHIter voh_it = polyhedron->voh_begin(vh); voh_it != polyhedron->voh_end(vh); ++voh_it) {
		MyTriMesh::HHandle heh = voh_it.current_halfedge_handle();
		if (polyhedron->to_vertex_handle(heh) == v2) {
			return polyhedron->edge_handle(heh);
		}
	}

	return polyhedron->InvalidEdgeHandle;
}

#pragma endregion

#pragma region Smoothing

#pragma region Quaternion Smoothing

void SQMNode::smoothLIE(LIE lie) {
	//angle and axis of rotation
	float div = lie.edges.size();
	float alfa = acos(lie.quaternion.s) * 2.0 / div;//lie.quaternion.s / div;
	float partAlfa = alfa;
	CVector3 axis = Normalize(CVector3(lie.quaternion.i, lie.quaternion.j, lie.quaternion.k));
	//CVector3 offset = CVector3(position.values_);
	CVector3 offset = CVector3(centerOfMass.values_);
	//get first vertice
	MyTriMesh::HHandle heh = lie.firstHHandle;
	CVector3 v = CVector3(polyhedron->point(polyhedron->to_vertex_handle(polyhedron->opposite_halfedge_handle(heh))).values_) - offset;
	while (heh != lie.lastHHandle) {
		Quaternion q = QuaternionFromAngleAxis(alfa, axis);
		MyTriMesh::VHandle vh = polyhedron->to_vertex_handle(heh);
		CVector3 u = QuaternionRotateVector(q, v);
		u = u + offset;
		MyTriMesh::Point P(u.x, u.y, u.z);
		OpenMesh::Vec3f ray_origin = centerOfMass;
		OpenMesh::Vec3f ray_dir = (P - centerOfMass).normalize();
		float t = 0;
		raySphereIntersection(ray_origin, ray_dir, position, nodeRadius, t);
		P = ray_origin + (t * ray_dir);

		polyhedron->set_point(vh, P);

		alfa += partAlfa;
		heh = nextLink(heh);
		vh = polyhedron->to_vertex_handle(heh);
	}
}

void SQMNode::smoothLIEs(map<int, LIENeedEntry> lieMap) {
	for (map<int, LIENeedEntry>::iterator it = lieMap.begin(); it != lieMap.end(); it++) {
		LIENeedEntry lieNeed = it->second;
		for (int i = 0; i < lieNeed.lies.size(); i++) {
			smoothLIE(lieNeed.lies[i]);
		}
	}
}

#pragma endregion

#pragma region Laplacian Smoothing

void SQMNode::laplacianSMoothing(SQMSmoothingAlgorithm algorithm) {
	//convert mesh
	MeshGraph meshGraph = MeshGraph();
	mesh2graph(meshGraph, algorithm);
	//smooth mesh
	computeLaplacian(&meshGraph);
	//translate vertices
	//replace vertices in mesh
	recalculateSmoothedVertices(meshGraph);
}

void SQMNode::mesh2graph(MeshGraph& meshGraph, SQMSmoothingAlgorithm algorithm) {
	if (algorithm == SQMOneRingLaplacianSmoothing) {
		mesh2graphOneRingWeighted(meshGraph);
	} else if (algorithm == SQMValencyLaplacianSmoothing) {
		mesh2graphValencyWeighted(meshGraph);
	}
}

void SQMNode::mesh2graphValencyWeighted(MeshGraph& meshGraph) {
	int n_vertices = polyhedron->n_vertices();
	int n_faces = polyhedron->n_faces();
	int n_edges = polyhedron->n_edges();

	meshGraph.numOfFaces = n_faces;
	meshGraph.numOfVertices = n_vertices;
	meshGraph.wL = 1;
	meshGraph.wH = new float[n_vertices];

	CVector3 *vertices = new CVector3[n_vertices];
	int idx = 0;
	//get the vertices
	for (MyTriMesh::VertexIter v_it = polyhedron->vertices_begin(); v_it != polyhedron->vertices_end(); ++v_it) 
	{
		MyTriMesh::VHandle vh = v_it.handle();
		MyTriMesh::Point P = polyhedron->point(vh);
		vertices[idx] = CVector3(polyhedron->point(vh).values_);
		//valency weighted
		int count = 0;
		for (MyTriMesh::VVIter vv_it = polyhedron->vv_begin(v_it.handle()); vv_it != polyhedron->vv_end(v_it.handle()); ++vv_it) {
			count++;
		}
		meshGraph.wH[idx] = count;
		idx++;
	}

	meshGraph.pVerts = vertices;
	meshGraph.E = TNT::Array2D<bool>(n_vertices, n_vertices, false);

	for (MyTriMesh::EdgeIter e_it = polyhedron->edges_begin(); e_it != polyhedron->edges_end(); ++e_it) {
		MyTriMesh::HalfedgeHandle heh0 = polyhedron->halfedge_handle(e_it.handle(), 0);
		MyTriMesh::HalfedgeHandle heh1 = polyhedron->halfedge_handle(e_it.handle(), 1);
		int i = polyhedron->to_vertex_handle(heh0).idx();
		int j = polyhedron->to_vertex_handle(heh1).idx();
		//register to map
		meshGraph.E[i][j] = true;
		meshGraph.E[j][i] = true;
	}
}

void SQMNode::mesh2graphOneRingWeighted(MeshGraph& meshGraph) {
	int n_vertices = polyhedron->n_vertices();
	int n_faces = polyhedron->n_faces();
	int n_edges = polyhedron->n_edges();

	meshGraph.numOfFaces = n_faces;
	meshGraph.numOfVertices = n_vertices;
	meshGraph.wL = 1;
	meshGraph.wH = new float[n_vertices];

	CVector3 *vertices = new CVector3[n_vertices];
	int idx = 0;
	//avarage face area
	float area = 0;
	for (MyTriMesh::FaceIter f_it = polyhedron->faces_begin(); f_it != polyhedron->faces_end(); ++f_it)	{
		MyTriMesh::FVIter fv_it = polyhedron->fv_begin(f_it.handle());
		MyTriMesh::Point A = polyhedron->point(fv_it.handle());
		++fv_it;
		MyTriMesh::Point B = polyhedron->point(fv_it.handle());
		++fv_it;
		MyTriMesh::Point C = polyhedron->point(fv_it.handle());

		OpenMesh::Vec3f u = B - A;
		OpenMesh::Vec3f v = C - A;
		OpenMesh::Vec3f w = cross(u, v);
		area += (w.norm() / 2.0);
	}
	area /= (float)n_faces;
	//formula from paper
	meshGraph.wL = sqrtf(area)/1000.0;
	//get the vertices
	for (MyTriMesh::VertexIter v_it = polyhedron->vertices_begin(); v_it != polyhedron->vertices_end(); ++v_it) 
	{
		MyTriMesh::VHandle vh = v_it.handle();
		MyTriMesh::Point P = polyhedron->point(vh);
		vertices[idx] = CVector3(polyhedron->point(vh).values_);
		//one ring area weighted 1.0 from paper
		meshGraph.wH[idx] = 1.0;
		idx++;
	}

	meshGraph.pVerts = vertices;
	meshGraph.E = TNT::Array2D<bool>(n_vertices, n_vertices, false);

	for (MyTriMesh::EdgeIter e_it = polyhedron->edges_begin(); e_it != polyhedron->edges_end(); ++e_it) {
		MyTriMesh::HalfedgeHandle heh0 = polyhedron->halfedge_handle(e_it.handle(), 0);
		MyTriMesh::HalfedgeHandle heh1 = polyhedron->halfedge_handle(e_it.handle(), 1);
		int i = polyhedron->to_vertex_handle(heh0).idx();
		int j = polyhedron->to_vertex_handle(heh1).idx();
		//register to map
		meshGraph.E[i][j] = true;
		meshGraph.E[j][i] = true;
	}
}

void SQMNode::recalculateSmoothedVertices(MeshGraph& meshGraph) {
	//update vertices with new position
	int index = 0;
	for (MyTriMesh::VertexIter v_it = polyhedron->vertices_begin(); v_it != polyhedron->vertices_end(); ++v_it) 
	{
		MyTriMesh::VHandle vh = v_it.handle();
		//setup point
		CVector3 vec = meshGraph.pVerts[index];
		MyTriMesh::Point P = MyTriMesh::Point(vec.x, vec.y, vec.z);//polyhedron->point(vh);
		polyhedron->set_point(vh, P);
		index++;
	}
	//translate vertices on the sphere
	//not always in the center!!!!
	//CVector3 offset(position.values_);
	CVector3 offset(centerOfMass.values_);
	for (MyTriMesh::VertexIter v_it = polyhedron->vertices_begin(); v_it != polyhedron->vertices_end(); ++v_it) 
	{
		MyTriMesh::VHandle vh = v_it.handle();
		CVector3 P(polyhedron->point(vh).values_);
		P = Normalize(P - offset);
		P = P * nodeRadius;
		P = P + offset;
		MyTriMesh::Point Q(P.x, P.y, P.z);
		polyhedron->set_point(vh, Q);
	}
}

#pragma endregion

#pragma region Avaraging Smoothing

void SQMNode::smoothLIEByAvaraging(LIE lie) {
	MyTriMesh::HHandle heh = lie.lastHHandle;
	while (heh != lie.firstHHandle) {
		MyTriMesh::VHandle vh = polyhedron->to_vertex_handle(polyhedron->opposite_halfedge_handle(heh));

		MyTriMesh::Point P(0, 0, 0);
		float count = 0; 
		for (MyTriMesh::VVIter vv_it = polyhedron->vv_begin(vh); vv_it != polyhedron->vv_end(vh); ++vv_it) {
			P += polyhedron->point(vv_it.handle());
			count++;
		}
		P /= count;
		OpenMesh::Vec3f ray_origin(P[0], P[1], P[2]);
		OpenMesh::Vec3f ray_dir = (ray_origin - centerOfMass).normalize();
		float t = 0;
		if (raySphereIntersection(ray_origin, ray_dir, position, nodeRadius, t)) {
			P = ray_origin + (t * ray_dir);
		}

		polyhedron->set_point(vh, P);

		heh = prevLink(heh);
	}
}

#pragma endregion

#pragma endregion

#pragma region BNP Joining

void SQMNode::joinBNPs(MyMesh* mesh, SQMNode* parentBNPNode, vector<MyMesh::VertexHandle> oneRing, OpenMesh::Vec3f directionVector) {
	if (this->isLeafNode()) {
		//finish mesh somehow
		finishLeafeNode(mesh, oneRing);
		return;
	}
	if (this->isBranchNode()) {
		addPolyhedronToMesh(mesh, parentBNPNode, oneRing, directionVector);
	} else {//connection node with one son
		extendMesh(mesh, parentBNPNode, oneRing, directionVector);
	}
}

void SQMNode::addPolyhedronToMesh(MyMesh* mesh, SQMNode* parentBNPNode, vector<MyMesh::VertexHandle>& oneRing, OpenMesh::Vec3f& directionVector) {
	//connect mesh to existing one
	vector<vector<MyMesh::VHandle> > oneRingsOfPolyhedron;
	addPolyhedronAndRememberVHandles(mesh, parentBNPNode, oneRing, oneRingsOfPolyhedron, directionVector);
	//for each son send him points on intersection perimeter
	for (int i = 0; i < nodes.size(); i++) {
		OpenMesh::Vec3f u = (nodes[i]->getPosition() - position).normalize();
		nodes[i]->joinBNPs(mesh, this, oneRingsOfPolyhedron[i], u);
	}
}
//!!!!!ORDERING OF INSERTED EDGES IS IMPORTANT!!!!!!!
void SQMNode::addPolyhedronAndRememberVHandles(MyMesh* mesh, SQMNode* parentBNPNode, vector<MyMesh::VertexHandle>& oneRing, vector<vector<MyMesh::VHandle> >& oneRingsOfPolyhedron, OpenMesh::Vec3f& directionVector) {
	//get non intersecting vhandles and remember one rings
	int vectorSize = intersectionVHandles.size();
	vector<vector<int> > intersectionOneRingIndexes(vectorSize);
	vector<MyMesh::VHandle> addedVHandles;
	for (MyTriMesh::VIter v_it = polyhedron->vertices_begin(); v_it != polyhedron->vertices_end(); ++v_it) {
		MyTriMesh::VHandle vhandle = v_it.handle();
		int position = getPositionInArray<MyTriMesh::VHandle>(vhandle, intersectionVHandles);
		if (position == -1) {
			MyTriMesh::Point P = polyhedron->point(vhandle);
			MyMesh::VHandle vhandle = mesh->add_vertex(P);
			addedVHandles.push_back(vhandle);
			meshVhandlesToRotate.push_back(vhandle);
		} else {
			addedVHandles.push_back(vhandle);
			//collect one ring indexis
			for (MyTriMesh::VVIter vv_it = polyhedron->vv_begin(vhandle); vv_it != polyhedron->vv_end(vhandle); ++vv_it) {
				intersectionOneRingIndexes[position].push_back(vv_it.handle().idx());
			}
		}
	}
	//prepare one rings
	for (int i = 0; i < intersectionOneRingIndexes.size(); i++) {
		vector<MyMesh::VHandle> vhandles;
		for (int j = 0; j < intersectionOneRingIndexes[i].size(); j++) {
			int vertexIndex = intersectionOneRingIndexes[i][j];
			vhandles.push_back(addedVHandles[vertexIndex]);
		}
		oneRingsOfPolyhedron.push_back(vhandles);
	}
	//connect to the rest of the mesh
	if (parentBNPNode != NULL) {
		//order newOneRingArray
		vector<MyMesh::VHandle> oldOneRing = oneRingsOfPolyhedron.back();
		//flip one ring if it has different orientation
		bool shouldFlip = sameOneRingOrientation(mesh, oldOneRing, oneRing, position, parentBNPNode->getPosition(), directionVector);
		if (shouldFlip) {
			vector<MyMesh::VHandle> flipped;
			flipVector(oldOneRing, flipped);
			oldOneRing = flipped;
		}
		//find closest point
		/*vector<MyMesh::VertexHandle> newOneRing;
		MyMesh::Point P = mesh->point(oneRing[0]);
		float minDist = 0;
		int index = 0;
		for (int i = 0; i < oldOneRing.size(); i++) {
		MyMesh::Point Q = mesh->point(oldOneRing[i]);
		//float dist = (P - Q).norm();
		OpenMesh::Vec3f vv = -directionVector.normalized();
		float dist = dot((P - Q).normalized(), vv);
		//it is an angle
		if (i == 0 || dist > minDist) {
		minDist = dist;
		index = i;
		}
		}*/
		//find point which will have smallest sum of distances
		vector<MyMesh::VertexHandle> newOneRing;
		MyMesh::Point P = mesh->point(oneRing[0]);
		float minDist = 0;
		int index = 0;
		//minimum distance approach
		for (int i = 0; i < oldOneRing.size(); i++) {
			int k = i;
			float dist = 0;
			for (int j = 0; j < oneRing.size(); j++) {
				if (k >= oldOneRing.size()) k = 0;
				MyMesh::Point A = mesh->point(oneRing[j]);
				MyMesh::Point B = mesh->point(oldOneRing[k]);
				dist += (A - B).norm();
				k++;
			}
			if (i == 0 || dist < minDist) {
				minDist = dist;
				index = i;
			}
		}
		/*//dot product approach
		OpenMesh::Vec3f vv = -directionVector.normalized();
		for (int i = 0; i < oldOneRing.size(); i++) {
		int k = i;
		float dist = 0;
		for (int j = 0; j < oneRing.size(); j++) {
		if (k >= oldOneRing.size()) k = 0;
		MyMesh::Point A = mesh->point(oneRing[j]);
		MyMesh::Point B = mesh->point(oldOneRing[k]);
		dist += dot((A - B).normalized(), vv);
		k++;
		}
		if (i == 0 || dist > minDist) {
		minDist = dist;
		index = i;
		}
		}*/
		//reorder array
		for (int i = 0; i < oldOneRing.size(); i++) {
			if (index == oldOneRing.size()) {
				index = 0;
			}
			newOneRing.push_back(oldOneRing[index]);
			index++;
		}
		//create new faces for the points
		for (int i = 0; i < newOneRing.size(); i++) {
			int j = 0;
			if (i + 1 < newOneRing.size()) {
				j = i + 1;
			}
			mesh->add_face(oneRing[i], newOneRing[i], newOneRing[j], oneRing[j]);
		}
	}
}

void SQMNode::extendMesh(MyMesh* mesh, SQMNode* parentBNPNode, vector<MyMesh::VertexHandle>& oneRing, OpenMesh::Vec3f& directionVector) {
	//create new points by translating parent points in the direction of the vector
	float d = -(position[0]*directionVector[0] + position[1]*directionVector[1] + position[2]*directionVector[2]);
	vector<MyMesh::Point> points;
	for (int i = 0; i < oneRing.size(); i++) {
		MyMesh::Point P = mesh->point(oneRing[i]);
		//OpenMesh::Vec3f u = (position - parentBNPNode->getPosition()).norm() * directionVector;
		float dist = fabsf(directionVector[0]*P[0] + directionVector[1]*P[1] + directionVector[2]*P[2] + d);
		OpenMesh::Vec3f u = directionVector * dist;
		MyMesh::Point Q = P + u;
		points.push_back(Q);
	}
	//map new points to node sphere
	for (int i = 0; i < points.size(); i++) {
		MyMesh::Point P = points[i];
		OpenMesh::Vec3f u(P[0], P[1], P[2]);
		u = (u - position).normalize();
		OpenMesh::Vec3f v = nodeRadius*u + position;
		MyMesh::Point Q(v[0], v[1], v[2]);
		points[i] = Q;
	}
	//insert new points into mesh
	vector<MyMesh::VertexHandle> newOneRing;
	for (int i = 0; i < points.size(); i++) {
		MyMesh::VHandle vhandle = mesh->add_vertex(points[i]);
		newOneRing.push_back(vhandle);
		meshVhandlesToRotate.push_back(vhandle);
	}
	//create new faces for the points
	vector<MyMesh::FaceHandle> temp;
	for (int i = 0; i < newOneRing.size(); i++) {
		int j = 0;
		if (i + 1 < newOneRing.size()) {
			j = i + 1;
		}
		temp.push_back(mesh->add_face(oneRing[i], newOneRing[i], newOneRing[j], oneRing[j]));
	}
	//remember inseted points
	meshIntersectionVHandles = newOneRing;
	//continue joining BNPs
	nodes[0]->joinBNPs(mesh, this, newOneRing, directionVector);
}

void SQMNode::finishLeafeNode(MyMesh* mesh, vector<MyMesh::VertexHandle>& oneRing) {
	//TODO finish as desired
	MyMesh::Point P(position[0], position[1], position[2]);
	MyMesh::VHandle vhandle = mesh->add_vertex(P);		
	meshVhandlesToRotate.push_back(vhandle);
	for (int i = 0; i < oneRing.size(); i++) {	
		int j = 0;
		if (i + 1 < oneRing.size()) {
			j = i + 1;
		}
		mesh->add_face(oneRing[i], vhandle, oneRing[j]);
	}
}

#pragma endregion

#pragma region Final Vertex Placement

void SQMNode::rotateBack(MyMesh *mesh) {
	//start with sons
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->rotateBack(mesh);
	}
	//rotate
	if (parent != NULL) {
		OpenMesh::Vec3f parentPosition = parent->getPosition();
		CVector3 parentPos(parentPosition.values_);
		for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
			MyMesh::VHandle vhandle = meshVhandlesToRotate[i];
			MyMesh::Point P = mesh->point(vhandle);
			CVector3 v(P.values_);
			//translate rotate translate
			v = v - parentPos;
			v = QuaternionRotateVector(axisAngle, v);
			v = v + parentPos;
			SQMNode* ancestor = lastBranchNodeInChain(this);
			if (parent->isBranchNode() && ancestor != NULL && ancestor->getParent() != NULL) {
				CVector3 offset(ancestor->getParent()->getPosition().values_);
				v = v - offset;
				v = QuaternionRotateVector(ancestor->getAxisAngle(), v);
				v = v + offset;
			}
			P[0] = v.x;
			P[1] = v.y;
			P[2] = v.z;
			mesh->set_point(vhandle, P);
		}
	}
	for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
		MyMesh::VHandle vhandle = meshVhandlesToRotate[i];
		MyMesh::Point P = mesh->point(vhandle);
		glm::vec4 u(P[0], P[1], P[2], 1);
		u = transformationMatrix * u;
		P[0] = u.x;
		P[1] = u.y;
		P[2] = u.z;
		mesh->set_point(vhandle, P);
	}
	position = oldPosition;
}

void SQMNode::rotateWithSkeleton(MyMesh *mesh, SkinSkeleton *skeleton) {
	//rotate position
	glm::vec4 newPos(position[0], position[1], position[2], 1.0f);
	newPos = skeleton->matrix * newPos;
	position = OpenMesh::Vec3f(newPos.x, newPos.y, newPos.z);
	//rotate mesh vertices
	for (int i = 0; i  < meshVhandlesToRotate.size(); i ++) {
		MyMesh::VHandle vhandle = meshVhandlesToRotate[i];
		MyMesh::Point P = mesh->point(vhandle);
		glm::vec4 pos = glm::vec4(P[0], P[1], P[2], 1.0f);
		pos = skeleton->matrix * pos;

		if (this->isConnectionNode() && (sqmNodeType != SQMFormerCapsule && sqmNodeType != SQMCreatedCapsule)) {
			//is a connection node we should combine two matrices
			glm::vec4 pos2(P[0], P[1], P[2], 1.0f);
			pos2 = skeleton->nodes[0]->matrix * pos2;
			pos = 0.5f*pos + 0.5f*pos2;
		}

		//csutom transformations
		pos = transformationMatrix * pos;

		P = OpenMesh::Vec3f(pos.x, pos.y, pos.z);
		mesh->set_point(vhandle, P);
	}
	//rotate next
	for (int i = 0; i < nodes.size(); i++) {
		SkinSkeleton *next = skeleton;
		if (sqmNodeType != SQMFormerCapsule && sqmNodeType != SQMCreatedCapsule) {
			next = skeleton->nodes[i];
		}
		if (parent != NULL && parent->getSQMNodeType() == SQMCreatedCapsule && sqmNodeType == SQMFormerCapsule) {
			next = skeleton->nodes[i];
		}
		nodes[i]->rotateWithSkeleton(mesh, next);
	}
}

#pragma endregion

#pragma region Skinning Matrix Setup

void SQMNode::setupSkinningMatrixIDs(SkinSkeleton *skeleton) {
	int id1 = skeleton->id, id2 = -1;
	if (this->isConnectionNode() && (sqmNodeType != SQMFormerCapsule && sqmNodeType != SQMCreatedCapsule)) {
		id2 = skeleton->nodes[0]->id;
	}
	skinningIDs = glm::ivec2(id1, id2);
	
	for (int i = 0; i < nodes.size(); i++) {
		SkinSkeleton *next = skeleton;
		if (sqmNodeType != SQMFormerCapsule && sqmNodeType != SQMCreatedCapsule) {
			next = skeleton->nodes[i];
		}
		if (parent != NULL && parent->getSQMNodeType() == SQMCreatedCapsule && sqmNodeType == SQMFormerCapsule) {
			next = skeleton->nodes[i];
		}
		nodes[i]->setupSkinningMatrixIDs(next);
	}
}

#pragma endregion

#pragma region BNP Tesselation

void SQMNode::getMeshTessLevel(std::vector<float> &tessLevels) {
	//each node stores vertices for rotation that belong to him we just need to push tess level for every vertex
	for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
		tessLevels.push_back(tessLevel);
	}
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->getMeshTessLevel(tessLevels);
	}
}

void SQMNode::getMeshTessDataf(std::vector<float> &tessLevels, std::vector<float> &nodePositions, std::vector<float> &data) {
	float type = 0.0;
	if (this->isBranchNode()) type = 1.0;
	else if (this->isLeafNode()) type = 2.0;

	for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
		tessLevels.push_back(tessLevel);
		nodePositions.push_back(position[0]);
		nodePositions.push_back(position[1]);
		nodePositions.push_back(position[2]);

		data.push_back(type);
		data.push_back(id);
	}
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->getMeshTessDataf(tessLevels, nodePositions, data);
	}
}

void SQMNode::getMeshTessDatai(vector<float> &tessLevels, vector<float> &nodePositions, std::vector<int> &skinMatrices, vector<int> &data) {
	int type = 0;
	if (this->isBranchNode()) type = 1;
	else if (this->isLeafNode()) type = 2;
	if (type == 0 && nodes[0]->isLeafNode()) {
		type = 2;
	}

	for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
		MyMesh::VHandle vh = meshVhandlesToRotate[i];
		tessLevels.push_back(tessLevel);
		//tessLevels.push_back(type);
		//tessLevels.push_back(id);

		nodePositions.push_back(position[0]);
		nodePositions.push_back(position[1]);
		nodePositions.push_back(position[2]);
		
		skinMatrices.push_back(skinningIDs.x);
		skinMatrices.push_back(skinningIDs.y);

		data.push_back(type);
		data.push_back(id);
	}
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->getMeshTessDatai(tessLevels, nodePositions, skinMatrices, data);
	}
}

void SQMNode::fillRadiusTable(float *table, int width) {
	bool isLeaf = this->isLeafNode();
	bool isBranch = this->isBranchNode();
	bool isConnection = (!isLeaf && !isBranch);
	//always store original radius
	table[id*width + id] = nodeRadius;
	//store radius from yourself to connected nodes
	if (isLeaf) {//leaf is max radius
		int plus = (parent == NULL) ? 0 : parent->getId();
		table[id*width + plus] = nodeRadius;
		if (nodes.size() != 0) {//only the top of the worm is a leaf with childs
			table[id*width + nodes[0]->getId()] = nodeRadius;
		}
	}
	if (isConnection) {//connection is max radius
		table[id*width + parent->getId()] = nodeRadius;
		table[id*width + nodes[0]->getId()] = nodeRadius;
	}
	if (isBranch) {//branch has one-ring radiuses instead of node radius
		//for each node get one-ring radius and store
		for (int i = 0; i < intersectionVHandles.size(); i++) {
			float radius = calculateOneRingRadius(polyhedron, intersectionVHandles[i]);
			int id2 = (i >= nodes.size()) ? parent->getId() : nodes[i]->getId();
			table[id*width + id2] = radius;
		}
	}

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->fillRadiusTable(table, width);
	}
}

void SQMNode::fillCentersTable(float *table) {
	table[id*3 + 0] = position[0];
	table[id*3 + 1] = position[1];
	table[id*3 + 2] = position[2];

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->fillCentersTable(table);
	}
}

void SQMNode::calculateOneRingRadiusAndMap(std::vector<float> &oneRingRadius, std::map<int, std::vector<int> > &intersectionMap) {
	for (int i = 0; i < intersectionVHandles.size(); i++) {
		MyTriMesh::VHandle vh = intersectionVHandles[i];
		MyTriMesh::VHandle first_vh(-1);
		MyTriMesh::VHandle prev_vh(-1);
		float radius = 0;
		for (MyTriMesh::VVIter vv_it = polyhedron->vv_begin(vh); vv_it != polyhedron->vv_end(vh); ++vv_it) {
			if (first_vh.idx() == -1) {
				first_vh = vv_it.handle();
				prev_vh = vv_it.handle();
			}
			//add to map
			int id = vv_it.handle().idx();
			map<int, vector<int> >::iterator it = intersectionMap.find(id);
			if (it != intersectionMap.end()) {
				vector<int> oneRings = intersectionMap[id];
				oneRings.push_back(i);
				intersectionMap[id] = oneRings;
			} else {
				vector<int> oneRings;
				oneRings.push_back(i);
				intersectionMap[id] = oneRings;
			}

			//calculate length
			MyTriMesh::Point P = polyhedron->point(prev_vh);
			MyTriMesh::Point Q = polyhedron->point(vv_it.handle());
			radius += (P - Q).norm();

			prev_vh = vv_it.handle();
		}
		MyTriMesh::Point P = polyhedron->point(prev_vh);
		MyTriMesh::Point Q = polyhedron->point(first_vh);
		radius += (P - Q).norm();
		//we calculated circumference (2*PI*r) so now we need to divide by 2*PI
		radius /= (2*M_PI);
		oneRingRadius.push_back(radius);
	}
}

#pragma endregion

#pragma region SQM Special Cases

#pragma region Single Node

void SQMNode::createScaledIcosahderon(MyMesh* mesh) {
	float scale = nodeRadius;
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.000f,  0.000f,  1.000f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.894f,  0.000f,  0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.276f,  0.851f,  0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(-0.724f,  0.526f,  0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(-0.724f, -0.526f,  0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.276f, -0.851f,  0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.724f,  0.526f, -0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(-0.276f,  0.851f, -0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(-0.894f,  0.000f, -0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(-0.276f, -0.851f, -0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.724f, -0.526f, -0.447f) * scale + position));
	meshVhandlesToRotate.push_back(mesh->add_vertex(MyMesh::Point(0.000f,  0.000f, -1.000f) * scale + position));

	mesh->add_face(meshVhandlesToRotate[2], meshVhandlesToRotate[1], meshVhandlesToRotate[0]);
	mesh->add_face(meshVhandlesToRotate[3], meshVhandlesToRotate[2], meshVhandlesToRotate[0]);
	mesh->add_face(meshVhandlesToRotate[4], meshVhandlesToRotate[3], meshVhandlesToRotate[0]);
	mesh->add_face(meshVhandlesToRotate[5], meshVhandlesToRotate[4], meshVhandlesToRotate[0]);
	mesh->add_face(meshVhandlesToRotate[1], meshVhandlesToRotate[5], meshVhandlesToRotate[0]);

	mesh->add_face(meshVhandlesToRotate[11], meshVhandlesToRotate[6], meshVhandlesToRotate[7]);
	mesh->add_face(meshVhandlesToRotate[11], meshVhandlesToRotate[7], meshVhandlesToRotate[8]);
	mesh->add_face(meshVhandlesToRotate[11], meshVhandlesToRotate[8], meshVhandlesToRotate[9]);
	mesh->add_face(meshVhandlesToRotate[11], meshVhandlesToRotate[9], meshVhandlesToRotate[10]);
	mesh->add_face(meshVhandlesToRotate[11], meshVhandlesToRotate[10], meshVhandlesToRotate[6]);

	mesh->add_face(meshVhandlesToRotate[1], meshVhandlesToRotate[2], meshVhandlesToRotate[6]);
	mesh->add_face(meshVhandlesToRotate[2], meshVhandlesToRotate[3], meshVhandlesToRotate[7]);
	mesh->add_face(meshVhandlesToRotate[3], meshVhandlesToRotate[4], meshVhandlesToRotate[8]);
	mesh->add_face(meshVhandlesToRotate[4], meshVhandlesToRotate[5], meshVhandlesToRotate[9]);
	mesh->add_face(meshVhandlesToRotate[5], meshVhandlesToRotate[1], meshVhandlesToRotate[10]);

	mesh->add_face(meshVhandlesToRotate[2], meshVhandlesToRotate[7], meshVhandlesToRotate[6]);
	mesh->add_face(meshVhandlesToRotate[3], meshVhandlesToRotate[8], meshVhandlesToRotate[7]);
	mesh->add_face(meshVhandlesToRotate[4], meshVhandlesToRotate[9], meshVhandlesToRotate[8]);
	mesh->add_face(meshVhandlesToRotate[5], meshVhandlesToRotate[10], meshVhandlesToRotate[9]);
	mesh->add_face(meshVhandlesToRotate[1], meshVhandlesToRotate[6], meshVhandlesToRotate[10]);
}

#pragma endregion

#pragma region Worm

void SQMNode::wormCreate(MyMesh *mesh, int vertices) {
	OpenMesh::Vec3f direction = (nodes[0]->getPosition() - position).normalize();
	//straighten worm
	oldPosition = position;
	nodes[0]->wormStraighten(direction);
	//first add self
	MyMesh::VHandle vh_Q = mesh->add_vertex(position);
	meshVhandlesToRotate.push_back(vh_Q);
	//create one ring
	OpenMesh::Vec3f axis = getAxisForCross(direction);
	OpenMesh::Vec3f normal = cross(direction, axis).normalize();
	float alfa = ((float)2.0*(float)M_PI)/(float)vertices;

	vector<MyMesh::VHandle> oneRing;
	MyMesh::Point P = nodes[0]->getPosition() + (nodes[0]->getNodeRadius() * normal);
	MyMesh::VHandle vh = mesh->add_vertex(P);
	oneRing.push_back(vh);
	nodes[0]->addVHandleToRotate(vh);

	CVector3 rotation = CVector3(direction.values_);
	CVector3 offset = CVector3(nodes[0]->getPosition().values_);
	CVector3 u = CVector3(P.values_) - offset;
	for (int i = 1; i < vertices; i++) {
		Quaternion q = QuaternionFromAngleAxis(alfa*(float)i, rotation);
		CVector3 v = QuaternionRotateVector(q, u);
		v = v + offset;

		MyMesh::Point X = MyMesh::Point(v.x, v.y, v.z);
		vh = mesh->add_vertex(X);
		oneRing.push_back(vh);
		nodes[0]->addVHandleToRotate(vh);
	}
	//extend
	//first self
	for (int i = 0; i < vertices; i++) {
		int j = (i == vertices - 1) ? 0 : (i + 1);
		mesh->add_face(oneRing[j], vh_Q, oneRing[i]);
	}
	//then others
	(*nodes[0]->getNodes())[0]->wormStep(mesh, oneRing, direction);
	//rotate worm back
	//wormFinalVertexPlacement(mesh);
}

void SQMNode::wormStraighten(OpenMesh::Vec3f lineVector) {
	axisAngle = Quaternion();
	oldPosition = position;
	//straighten self
	//translate parent to 0,0,0
	OpenMesh::Vec3f newPosition = position - parent->getPosition();
	CVector3 oldPos(newPosition.values_);
	//roatate
	float len = (originalPosition - parent->getOriginalPosition()).length();
	newPosition = len*lineVector;
	CVector3 newPos(newPosition.values_);
	Quaternion quaternion = SQMQuaternionBetweenVectors(newPos, oldPos);
	//translate back by parent
	newPosition = newPosition + parent->getPosition();
	//setup
	axisAngle = quaternion;
	position = newPosition;

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->wormStraighten(lineVector);
	}
}

void SQMNode::wormStep(MyMesh *mesh, vector<MyMesh::VHandle> &oneRing, OpenMesh::Vec3f lineVector) {
	if (nodes.size() == 0) {
		MyMesh::VHandle vh = mesh->add_vertex(position);
		meshVhandlesToRotate.push_back(vh);
		for (int i = 0; i < oneRing.size(); i++) {
			int j = (i == oneRing.size() - 1) ? 0 : (i + 1);
			mesh->add_face(oneRing[i], vh, oneRing[j]);
		}
		return;
	}
	//create new points by translating parent points in the direction of the vector
	float d = -(position[0]*lineVector[0] + position[1]*lineVector[1] + position[2]*lineVector[2]);
	vector<MyMesh::Point> points;
	for (int i = 0; i < oneRing.size(); i++) {
		MyMesh::Point P = mesh->point(oneRing[i]);
		//OpenMesh::Vec3f u = (position - parentBNPNode->getPosition()).norm() * directionVector;
		float dist = fabsf(lineVector[0]*P[0] + lineVector[1]*P[1] + lineVector[2]*P[2] + d);
		OpenMesh::Vec3f u = lineVector * dist;
		MyMesh::Point Q = P + u;
		points.push_back(Q);
	}
	//map new points to node sphere
	for (int i = 0; i < points.size(); i++) {
		MyMesh::Point P = points[i];
		OpenMesh::Vec3f u(P[0], P[1], P[2]);
		u = (u - position).normalize();
		OpenMesh::Vec3f v = nodeRadius*u + position;
		MyMesh::Point Q(v[0], v[1], v[2]);
		points[i] = Q;
	}
	//insert new points into mesh
	vector<MyMesh::VertexHandle> newOneRing;
	for (int i = 0; i < points.size(); i++) {
		MyMesh::VHandle vhandle = mesh->add_vertex(points[i]);
		newOneRing.push_back(vhandle);
		meshVhandlesToRotate.push_back(vhandle);
	}
	//create new faces for the points
	for (int i = 0; i < newOneRing.size(); i++) {
		int j = 0;
		if (i + 1 < newOneRing.size()) {
			j = i + 1;
		}
		mesh->add_face(oneRing[i], newOneRing[i], newOneRing[j], oneRing[j]);
	}
	//remember inseted points
	meshIntersectionVHandles = newOneRing;
	//continue on worm
	if (nodes.size() > 0) nodes[0]->wormStep(mesh, newOneRing, lineVector);
}

void SQMNode::wormFinalVertexPlacement(MyMesh *mesh) {
	//start with sons
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->wormFinalVertexPlacement(mesh);
	}
	//rotate
	if (parent != NULL) {
		OpenMesh::Vec3f parentPosition = parent->getPosition();
		CVector3 parentPos(parentPosition.values_);
		for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
			MyMesh::VHandle vhandle = meshVhandlesToRotate[i];
			MyMesh::Point P = mesh->point(vhandle);
			CVector3 v(P.values_);
			//translate rotate translate
			v = v - parentPos;
			v = QuaternionRotateVector(axisAngle, v);
			v = v + parentPos;
			P[0] = v.x;
			P[1] = v.y;
			P[2] = v.z;
			mesh->set_point(vhandle, P);
		}
	}
	for (int i = 0; i < meshVhandlesToRotate.size(); i++) {
		MyMesh::VHandle vhandle = meshVhandlesToRotate[i];
		MyMesh::Point P = mesh->point(vhandle);
		glm::vec4 u(P[0], P[1], P[2], 1);
		u = transformationMatrix * u;
		P[0] = u.x;
		P[1] = u.y;
		P[2] = u.z;
		mesh->set_point(vhandle, P);
	}
	position = oldPosition;
}

#pragma endregion

#pragma endregion

#pragma region Utility

void SQMNode::rotateSelfAndDescendants(Quaternion q, CVector3 offset) {
	if (this->isBranchNode()) {
		for (int i = 0; i < nodes.size(); i++) {
			nodes[i]->rotateSelfAndDescendants(q, offset);
		}
	}
	this->rotatePosition(q, offset);
}

SQMNode* SQMNode::getDescendantBranchNode(SQMNode* node) {
	if (node->isBranchNode()) {
		return node;
	}
	if (node->isLeafNode()) {
		return NULL;
	}

	return getDescendantBranchNode((*node->getDescendants())[0]);
}

SQMNode* SQMNode::getAncestorBranchNode(SQMNode* node) {
	if (node->isBranchNode()) {
		return node;
	}
	if (node->parent == NULL) {
		return NULL;
	}

	return getDescendantBranchNode(node->parent);
}

MyTriMesh::HalfedgeHandle SQMNode::startLink(MyTriMesh::VertexHandle vh) {
	//its next halfedge from the first outgoing halfedge from intersection vertice
	return polyhedron->next_halfedge_handle(polyhedron->voh_begin(vh).current_halfedge_handle());
}

MyTriMesh::HalfedgeHandle SQMNode::nextLink(MyTriMesh::HalfedgeHandle heh) {
	//get on the side of the next triangle
	//switch to next triangle half edges
	//return halfedge on the link
	return polyhedron->next_halfedge_handle(polyhedron->opposite_halfedge_handle(polyhedron->next_halfedge_handle(heh)));
}

MyTriMesh::HalfedgeHandle SQMNode::prevLink(MyTriMesh::HalfedgeHandle heh) {
	//two times next to get to the side of previous triangle
	//switch to previous triangle half edges
	//two times next halfedge to get the halfedge on the link
	return polyhedron->next_halfedge_handle(polyhedron->next_halfedge_handle(polyhedron->opposite_halfedge_handle(polyhedron->next_halfedge_handle(polyhedron->next_halfedge_handle(heh)))));
}

MyTriMesh::VHandle SQMNode::oppositeVHandle(MyTriMesh::HalfedgeHandle heh) {
	//switch to opposite triangle
	//get the next halfedge
	//get the vertex it points to
	return polyhedron->to_vertex_handle(polyhedron->next_halfedge_handle(polyhedron->opposite_halfedge_handle(heh)));
}


#pragma endregion

#pragma region Utility Functions

#pragma region Template Functions

template <typename T> int getPositionInArray(T& v, vector<T>& vectorArray) {
	for (int i = 0; i < vectorArray.size(); i++) {
		T u = vectorArray[i];
		if (equals<T>(u, v))
			return i;
	}

	return -1;
}

template <typename T> bool equals(T& a, T& b) {
	return (a == b);
}

template <> bool equals<MyMesh::VHandle>(MyMesh::VHandle& a, MyMesh::VHandle& b) {
	return a.idx() == b.idx();
}

template <typename T> void flipVector(vector<T>& toFlip, vector<T>& flipped) {
	for (vector<T>::reverse_iterator rit = toFlip.rbegin(); rit != toFlip.rend(); rit++) {
		flipped.push_back(*rit);
	}
}

#pragma endregion

#pragma region Position In Array

int getPointPositionInArrayOrAdd(OpenMesh::Vec3f& v, vector<OpenMesh::Vec3f>& vectorArray) {
	for (int i = 0; i < vectorArray.size(); i++) {
		OpenMesh::Vec3f u = vectorArray[i];
		if (OpenMeshVec3fEqual(v, u))
			return i;
	}
	vectorArray.push_back(v);
	return vectorArray.size() - 1;
}

int getPointPositionInArray(OpenMesh::Vec3f& v, vector<OpenMesh::Vec3f>& vectorArray) {
	for (int i = 0; i < vectorArray.size(); i++) {
		OpenMesh::Vec3f u = vectorArray[i];
		if (OpenMeshVec3fEqual(v, u))
			return i;
	}

	return -1;
}

int getPointPositionInArray(OpenMesh::Vec2i& v, vector<OpenMesh::Vec2i>& vectorArray) {
	for (int i = 0; i < vectorArray.size(); i++) {
		OpenMesh::Vec2i u = vectorArray[i];
		if (OpenMeshVec2iEqual(v, u))
			return i;
	}

	return -1;
}

#pragma endregion

#pragma region Mesh Functions

bool lesser(MyMesh::FaceHandle& a, MyMesh::FaceHandle& b) {
	return a.idx() < b.idx();
}

bool validTriFace(vector<MyMesh::VHandle>& faceVHandles) {
	if (faceVHandles.size() == 3) {
		if (faceVHandles[0].idx() == faceVHandles[1].idx() ||
			faceVHandles[0].idx() == faceVHandles[2].idx() ||
			faceVHandles[1].idx() == faceVHandles[2].idx()) {
				return false;
		}
		return true;
	}

	return false;
}

OpenMesh::Vec3i flipVec3i(OpenMesh::Vec3i& v) {
	return OpenMesh::Vec3i(v[2], v[1], v[0]);
}

bool sameOneRingOrientation(MyMesh* mesh, vector<MyMesh::VHandle>& oneRing, vector<MyMesh::VHandle>& oneRing2, OpenMesh::Vec3f& center1, OpenMesh::Vec3f& center2, OpenMesh::Vec3f& direction) {
	if (oneRing.size() > 1 && oneRing2.size() > 1) {
		//offset the position so we got vectors from O(0,0,0)
		MyMesh::Point P0 = mesh->point(oneRing[0]) - center1;
		MyMesh::Point P1 = mesh->point(oneRing[1]) - center1;
		MyMesh::Point Q0 = mesh->point(oneRing2[0]) - center2;
		MyMesh::Point Q1 = mesh->point(oneRing2[1]) - center2;
		CVector3 d(direction.values_);
		CVector3 A0(P0.values_);
		//projecting vectors onto the plane given by directionVector
		//it is enough to substract the projection of A0 onto d sience that is the distance between A0 and the plane
		A0 = A0 - d*Dot(A0, d);
		CVector3 A1(P1.values_);
		A1 = A1 - d*Dot(A1, d);
		CVector3 B0(Q0.values_);
		B0 = B0 - d*Dot(B0, d);
		CVector3 B1(Q1.values_);
		B1 = B1 - d*Dot(B1, d);
		//calculate base with cross product
		CVector3 u = Cross(A0, A1);
		CVector3 v = Cross(B0, B1);
		//project base onto direction vector
		float d1 = Dot(u, d);
		float d2 = Dot(v, d);
		//if both have same sign both have same direction 
		return (d1 >= 0 && d2 < 0) || (d1 < 0 && d2 >= 0);
	}
	return false;
}

int inLIEs(std::vector<LIE>& LIEs, MyTriMesh::VHandle vh1, MyTriMesh::VHandle vh2) {
	for (int i = 0; i < LIEs.size(); i++) {
		LIE lie = LIEs[i];
		if (lie.containsVertices(vh1, vh2))
			return i;
	}

	return -1;
}

#pragma endregion

#pragma region SQMNode Functions

int verticeDifferenceFatherSon(SQMNode* father, SQMNode* son, MyTriMesh::VHandle vhandle) {
	if (son != NULL) {			
		//get intersection vhandle its the last one from parent
		MyTriMesh::VHandle descendantVHandle = son->getIntersectionVHandles()->back();
		//get number of needed vertices
		int vhandleCount = 0;
		MyTriMesh* polyhedron = father->getPolyhedron();
		for (MyTriMesh::VVIter vvit = polyhedron->vv_begin(vhandle); vvit != polyhedron->vv_end(vhandle); ++vvit) {
			vhandleCount++;
		}
		int descendantVHandleCount = 0;
		polyhedron = son->getPolyhedron();
		for (MyTriMesh::VVIter vvit = polyhedron->vv_begin(descendantVHandle); vvit != polyhedron->vv_end(descendantVHandle); ++vvit) {
			descendantVHandleCount++;
		}
		return descendantVHandleCount - vhandleCount;
	}

	return 0;
}

SQMNode* lastBranchNodeInChain(SQMNode* node) {
	//nowhere to go
	if (node == NULL) return NULL;
	//get last branch node in a chain of branching nodes
	if (node->isBranchNode()) {
		if (node->getParent() != NULL && node->getParent()->isBranchNode()) {//parent exists and is a branch node
			return lastBranchNodeInChain(node->getParent());//continue search for last one
		}
		//else this is the last branching node in a chain
		return node;
	}
	//continiu until you find branch node
	return lastBranchNodeInChain(node->getParent());
}

#pragma endregion

#pragma endregion