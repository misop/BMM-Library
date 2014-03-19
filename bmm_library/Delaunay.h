#include <vector>
#include <glm\glm.hpp>

int CircumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3, double &xc, double &yc, double &r);
int Triangulate(std::vector<glm::vec3> &pxyz, std::vector<glm::ivec3> &v);