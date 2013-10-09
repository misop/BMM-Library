#include "sqm_library.h"
#include "SQMBasicController.h"

#pragma region Load From File

void SQMLoadFromFileAndExportToFile(std::string input, std::string output) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(input);
	sqm->executeSQMAlgorithm();
	sqm->exportMeshToFile(output);
	delete sqm;
}

void SQMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeletonFromFile(input);
	sqm->executeSQMAlgorithm();
	sqm->exportMeshToTriangles(points, indices);
	delete sqm;
}

#pragma endregion

#pragma region Load From Skeleton

void SQMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm();
	sqm->exportMeshToFile(output);
	delete sqm;
}

void SQMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices) {
	SQMBasicController *sqm = new SQMBasicController();
	sqm->loadSkeleton(skeleton);
	sqm->executeSQMAlgorithm();
	sqm->exportMeshToTriangles(points, indices);
	delete sqm;
}

#pragma endregion