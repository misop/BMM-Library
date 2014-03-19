#include "stdafx.h"
#include "SQMSkeletonNode.h"

using namespace MMath;

SQMSkeletonNode::SQMSkeletonNode (void) : point(CVector3()), scale(CVector3(1, 1, 1)), rotate(CVector3()), radius(10), tessLevel(1), id(0), cyclic(false), capsule(true)
{
}

SQMSkeletonNode::SQMSkeletonNode(float x, float y, float z) : point(CVector3(x, y, z)), scale(CVector3(1, 1, 1)), rotate(CVector3()), radius(10), tessLevel(1), id(0), cyclic(false), capsule(true) {
}

SQMSkeletonNode::SQMSkeletonNode(float x, float y, float z, float Radius) : point(CVector3(x, y, z)), scale(CVector3(1, 1, 1)), rotate(CVector3()), radius(Radius), tessLevel(1), id(0), cyclic(false), capsule(true) {
}


SQMSkeletonNode::~SQMSkeletonNode (void) {
	for (int i = 0; i < nodes.size(); i++) {
		delete nodes[i];
	}
}

void SQMSkeletonNode::addChild(SQMSkeletonNode *node) {
	nodes.push_back(node);
}
