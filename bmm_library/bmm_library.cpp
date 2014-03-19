#include "bmm_library.h"
#include "SQMBasicController.h"

#pragma region Load From File

void BMMLoadFromFileAndExportToFile(std::string input, std::string output, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(input);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToFile(output);
	delete sqm;
}

void BMMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(input);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToTriangles(points, indices);
	delete sqm;
}

#pragma endregion

#pragma region Load From Skeleton

void BMMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToFile(output);
	delete sqm;
}

void BMMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm(CPUSkinning);
	sqm->exportMeshToTriangles(points, indices);
	delete sqm;
}

#pragma endregion