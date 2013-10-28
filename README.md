SQM-Library
===========
Skeleton To Quad Dominant Mesh library which I developed as my thesis.
[Based on the original SQM algorithm.](http://wwwx.dtu.dk/English/Service/Phonebook.aspx?lg=showcommon&id=d52e0438-722a-4f62-ba1e-1c1e7fe6b18d)
## Quick Ussage Doc:
```cpp
void SQMLoadFromFileAndExportToFile(std::string input, std::string output);
void SQMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices);

void SQMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output);
void SQMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices);
//Files for loading are serialized SQMSkeletons with boost
//file output is .obj file
//vector ouput format is:
//  point: x0, y0, z0, x1, y1, z1, ...
//  indices: 0-2 first triangle, 3-5 second triangle, ...
```
####Used Libraries
SQM-library is using boost 1.51 and OpenMesh 2.2