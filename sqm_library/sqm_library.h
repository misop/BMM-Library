#include <string>
#include <vector>
#include "SQMSkeletonNode.h"

void SQMLoadFromFileAndExportToFile(std::string input, std::string output, bool CPUSkinning = false);
void SQMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning = false);

void SQMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output, bool CPUSkinning = false);
void SQMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning = false);