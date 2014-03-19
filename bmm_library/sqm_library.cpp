#include "sqm_library.h"
#include "SQMBasicController.h"

#pragma region Load From File

void SQMLoadFromFileAndExportToFile(std::string input, std::string output, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(input);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToFile(output);
	delete sqm;
}

void SQMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(input);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToTriangles(points, indices);
	delete sqm;
}

#pragma endregion

#pragma region Load From Skeleton

void SQMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToFile(output);
	delete sqm;
}

void SQMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToTriangles(points, indices);
	delete sqm;
}

#pragma endregion