BMM-Library
===========
Base Manifold Mesh (BMM) Library converts skeletons to quad dominant meshes.
Developed as my thesis. [Based on the original SQM algorithm.](http://wwwx.dtu.dk/English/Service/Phonebook.aspx?lg=showcommon&id=d52e0438-722a-4f62-ba1e-1c1e7fe6b18d)
## Quick Ussage Doc:
```cpp
void BMMLoadFromFileAndExportToFile(std::string input, std::string output, bool CPUSkinning = false);
void BMMLoadFromFileAndExportToVectors(std::string input, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning = false);

void BMMLoadFromSkeletonAndExportToFile(SQMSkeletonNode *skeleton, std::string output, bool CPUSkinning = false);
void BMMLoadFromSkeletonAndExportToVectors(SQMSkeletonNode *skeleton, std::vector<float> &points, std::vector<int> &indices, bool CPUSkinning = false);
```
####Usage
 ```
Files for loading are serialized BMMSkeletons with boost
file output is .obj file
vector ouput format is:
  point: x0, y0, z0, x1, y1, z1, ...
  indices: 0-2 first triangle, 3-5 second triangle, ...
By default the algorithm will produce straightened mesh
Rotation to original pose is done with skinning on CPU
```
####Used Libraries
BMM-library is using boost 1.51 and OpenMesh 2.2