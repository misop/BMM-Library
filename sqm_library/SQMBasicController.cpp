#include "stdafx.h"
#include "SQMBasicController.h"
#pragma warning(push, 0)
#include <boost\archive\archive_exception.hpp>
#include <boost\archive\basic_xml_iarchive.hpp>
#include <boost\archive\xml_iarchive.hpp>
#include <boost\archive\xml_oarchive.hpp>
#pragma pop
#include <iostream>
#include <fstream>
#include <glm/gtc/type_ptr.hpp>
#include "FloatArithmetic.h"


SQMBasicController::SQMBasicController(void)
{
	sqmALgorithm = new SQMAlgorithm();
}


SQMBasicController::~SQMBasicController(void)
{
	delete sqmALgorithm;
}

#pragma region Saving and Loading

void SQMBasicController::loadSkeletonFromFile(string fileName) {
	ofstream errorLog("log.txt");
	ifstream inputFile(fileName);
	assert(inputFile.good());
	boost::archive::xml_iarchive inputArchive(inputFile);
	SQMSkeletonNode *node = NULL;
	try {
		inputArchive >> BOOST_SERIALIZATION_NVP(node);	
	} catch (boost::archive::archive_exception e) {
		errorLog << "Exception: " << e.what() << endl;
		throw e;
	}
	SQMNode *sqmNode = new SQMNode(*node, NULL);
	if (sqmALgorithm != NULL) {
		delete sqmALgorithm;
	}
	sqmALgorithm = new SQMAlgorithm();
	sqmALgorithm->setRoot(sqmNode);
	delete node;
}

void SQMBasicController::loadSkeleton(SQMSkeletonNode *skeleton) {
	SQMNode *sqmNode = new SQMNode(*skeleton, NULL);
	if (sqmALgorithm != NULL) {
		delete sqmALgorithm;
	}
	sqmALgorithm = new SQMAlgorithm();
	sqmALgorithm->setRoot(sqmNode);
}

void SQMBasicController::saveSkeletonToFile(string fileName) {
	ofstream errorLog("log.txt");
	ofstream of(fileName);
	assert(of.good());
	boost::archive::xml_oarchive oa(of);
	SQMSkeletonNode *node = sqmALgorithm->getRoot()->exportToSkeletonNode();
	try {
		oa << BOOST_SERIALIZATION_NVP(node);	
	} catch (boost::archive::archive_exception e) {
		errorLog << "Exception: " << e.what() << endl;
		throw e;
	}
	delete node;
}

void SQMBasicController::exportMeshToFile(string fileName) {
	SQMState state = sqmALgorithm->getState();
	if (state == SQMJoinBNPs || state == SQMFinalPlacement) {
		writeMesh(sqmALgorithm->getMesh(), fileName);
	}
}

void SQMBasicController::exportMeshToTriangles(std::vector<float> &points, std::vector<int> &indices) {
	convertMeshToArray(sqmALgorithm->getMesh(), points, indices);
}

MyMesh* SQMBasicController::getMesh() {
	return sqmALgorithm->getMesh();
}

#pragma endregion

#pragma region SQM Computation

void SQMBasicController::straightenSkeleton() {
	sqmALgorithm->straightenSkeleton();
}

void SQMBasicController::computeConvexHull() {
	sqmALgorithm->computeConvexHull();
}

void SQMBasicController::subdivideConvexHull() {
	sqmALgorithm->subdivideConvexHull();
}

void SQMBasicController::joinBNPs() {
	sqmALgorithm->joinBNPs();
}

void SQMBasicController::executeSQMAlgorithm() {
	sqmALgorithm->executeSQMAlgorithm();
}

void SQMBasicController::executeSQMAlgorithm(SQMState state) {
	sqmALgorithm->executeSQMAlgorithm(state);
}

#pragma endregion