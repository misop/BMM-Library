#include <glm\glm.hpp>
#include <vector>
#include "mm_math.h"

struct SkinSkeleton {
	int id;
	MMath::CVector3 position;
	MMath::Quaternion axisAngle;
	MMath::Quaternion quaternion;
	SkinSkeleton *parent;
	glm::mat4 matrix;

	std::vector<SkinSkeleton*> nodes;

	SkinSkeleton();
	SkinSkeleton(float x, float y, float z);
	SkinSkeleton(SkinSkeleton *father, float x, float y, float z);
	SkinSkeleton(MMath::CVector3 pos);
	SkinSkeleton(MMath::CVector3 pos, MMath::Quaternion axis_Angle);
	~SkinSkeleton();

	bool isBNP();

	void CalculateCorrespondingDoF(SkinSkeleton *another);
	glm::mat4 ComputeLocalMatrix();
	void ComputeSkinningMatrices();
	void InvertMatrices();
	void PrecomputeFinalSkinningMatrices(SkinSkeleton *another);

	void ComputeCompoundRotation();
};