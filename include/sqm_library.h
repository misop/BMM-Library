#include <string>
#include <vector>
#include "SQMSkeletonNode.h"

void SQMLoadFromFileAndExportToFile(std::string input, std::string output);
void SQMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices);

void SQMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output);
void SQMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices);