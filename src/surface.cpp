#include "surface.hpp"
#include "enums.hpp"

namespace MySim
{
surface::surface(const std::string& fileName) : surfaceBase()
{
    readNodesAndTrianglesFromFile(fileName, geometryDefinitionType::stl);
    treeAABB_ = std::make_shared<treeAABB>(nodesOriginal_, trianglesOriginal_);
    treeAABB_->grow();
}

double surface::distanceToTriangle(const vec3d& point) const
{
    double dist = BIG;
    int nearestFaceID;
    double influenceRadius = 0.1;
    treeAABB_->getNearbyFacets(point, influenceRadius, dist, nearestFaceID);
    return std::sqrt(dist);
}

}  // namespace MySim
