#include <string>
#include <vector>
#include "SQMSkeletonNode.h"

void BMMLoadFromFileAndExportToFile(std::string input, std::string output, bool CPUSkinning = false);
void BMMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning = false);

void BMMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output, bool CPUSkinning = false);
void BMMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning = false);