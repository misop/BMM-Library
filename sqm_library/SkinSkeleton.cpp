#include "stdafx.h"
#include "SkinSkeleton.h"
#include <glm\gtc\type_ptr.hpp>
#include <glm\gtc\matrix_transform.hpp>

using namespace MMath;

SkinSkeleton::SkinSkeleton() : parent(NULL), position(CVector3()), axisAngle(Quaternion()) {
}

SkinSkeleton::SkinSkeleton(SkinSkeleton *father, float x, float y, float z) : parent(father), position(CVector3(x, y, z)), axisAngle(Quaternion()) {
}

SkinSkeleton::SkinSkeleton(float x, float y, float z) : parent(NULL), position(CVector3(x, y, z)), axisAngle(Quaternion()) {
}

SkinSkeleton::SkinSkeleton(CVector3 pos) : parent(NULL), position(pos), axisAngle(Quaternion()) {
}

SkinSkeleton::SkinSkeleton(CVector3 pos, Quaternion axis_Angle) : parent(NULL), position(pos), axisAngle(axis_Angle) {
}

SkinSkeleton::~SkinSkeleton() {
	for (int i = 0; i < nodes.size(); i++) {
		delete nodes[i];
	}
}

bool SkinSkeleton::isBNP() {
	int minus = (parent == NULL) ? 0 : 1;
	return (nodes.size() >= 3 - minus);
}

void SkinSkeleton::CalculateCorrespondingDoF(SkinSkeleton *bind) {
	//have to be the same skeleton just posed diferently
	//this is in bind position
	axisAngle = CVector4(0, 0, 0, 1);
	if (parent != NULL) {
		if (parent->isBNP()) {
			//no need to rotate save identity
			quaternion = Quaternion();
			axisAngle = CVector4(0, 0, 0, 1);
		} else {
			CVector3 u;
			CVector3 v;
			if (parent->parent == NULL) {
				u = Normalize(bind->position - bind->parent->position);
				v = Normalize(position - parent->position);
			} else {
				u = Normalize(parent->position - parent->parent->position);
				v = Normalize(position - parent->position);
			}
			Quaternion q = SQMQuaternionBetweenVectors(u, v);
			axisAngle = QuaternionToAxisAngle(q);
			axisAngle.s = axisAngle.s*180.0f/M_PI;
			quaternion = q;
		}
	}

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->CalculateCorrespondingDoF(bind->nodes[i]);
	}
}

glm::mat4 SkinSkeleton::ComputeLocalMatrix() {
	//should have parent
	CVector3 parentPos = parent->position;
	//translate to parent, rotate, translate back
	float angle = axisAngle.s;//*180.0f/M_PI;
	glm::vec3 pos = glm::vec3(parentPos.x, parentPos.y, parentPos.z);
	glm::vec3 axis = glm::vec3(axisAngle.i, axisAngle.j, axisAngle.k);
	//translate parent to origin
	glm::mat4 translate = glm::translate(glm::mat4(1), -pos);
	glm::mat4 rotate;
	glm::mat4 translateBack = glm::translate(glm::mat4(1), pos);
	//rotate by axis angle
	if (glm::length(axis) != 0) {
		//result = glm::rotate(result, angle, axis);
		rotate =  glm::rotate(glm::mat4(1), angle, axis);
	}
	//translate back to parent
	//result = glm::translate(result, pos);
	glm::mat4 result = translateBack * rotate * translate;

	return result;
}

void SkinSkeleton::ComputeSkinningMatrices() {
	if (parent == NULL) {
		matrix = glm::mat4();
	} else {
		glm::mat4 local = ComputeLocalMatrix();
		matrix = local * parent->matrix;
	}

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->ComputeSkinningMatrices();
	}
}

void SkinSkeleton::InvertMatrices() {
	matrix = glm::inverse(matrix);

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->InvertMatrices();
	}
}

void SkinSkeleton::PrecomputeFinalSkinningMatrices(SkinSkeleton *another) {
	//another is in Bind position
	matrix = another->matrix * matrix;

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->PrecomputeFinalSkinningMatrices(another->nodes[i]);
	}
}

void SkinSkeleton::ComputeCompoundRotation() {
	if (parent != NULL) {
		quaternion = parent->quaternion * quaternion;
	}

	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->ComputeCompoundRotation();
	}
}