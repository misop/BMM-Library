#include "stdafx.h"
#include "SkeletonNode.h"


SkeletonNode::SkeletonNode (void) : point(CVector3()), radius(10), id(0), cyclic(false)
{
}

SkeletonNode::SkeletonNode(float x, float y, float z) : point(CVector3(x, y, z)), radius(10), id(0), cyclic(false) {
}

SkeletonNode::SkeletonNode(float x, float y, float z, float Radius) : point(CVector3(x, y, z)), radius(Radius), id(0), cyclic(false) {
}


SkeletonNode::~SkeletonNode (void) {
	for (int i = 0; i < nodes.size(); i++) {
		delete nodes[i];
	}
}

void SkeletonNode::addChild(SkeletonNode *node) {
	nodes.push_back(node);
}
