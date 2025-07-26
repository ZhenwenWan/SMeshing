#include "MySimModTriangle.hpp"
#include <memory>

namespace MySim
{
MySimModTriangle::MySimModTriangle(
    UIN a, UIN b, UIN c, const std::vector<std::shared_ptr<node>>& nodesGlobal)
    : nodes_({a, b, c}), edges_()
{
    edges_[0][0] = a;
    edges_[0][1] = b;
    edges_[1][0] = b;
    edges_[1][1] = c;
    edges_[2][0] = c;
    edges_[2][1] = a;
    setCircumCenter(nodesGlobal);
    const vec3d& pA = nodesGlobal[a]->position_;
    const vec3d& pB = nodesGlobal[b]->position_;
    const vec3d& pC = nodesGlobal[c]->position_;
    const vec3d& AB = pB - pA;
    const vec3d& AC = pC - pA;
    area_ = (AB.cross(AC)).norm() / 2.0;
    if (area_ < EMINUS13) area_ = 0.0;
    vertices_.clear();
    vertices_.push_back(nodesGlobal[a]);
    vertices_.push_back(nodesGlobal[b]);
    vertices_.push_back(nodesGlobal[c]);
}

bool MySimModTriangle::isInCircumCircle(const vec3d& point) const
{
    if (distance(point, circumCenter_.centerPoint_) < circumCenter_.radius_)
    {
        return true;
    }
    return false;
}

bool MySimModTriangle::isInside(const vec3d& point) const
{
    const vec3d& pA = vertices_[0]->position_;
    const vec3d& pB = vertices_[1]->position_;
    const vec3d& pC = vertices_[2]->position_;
    const vec3d& DA = point - pA;
    const vec3d& DB = point - pB;
    const vec3d& DC = point - pC;
    const double& sAB = (DA.cross(DB)).norm() / 2.0;
    const double& sBC = (DB.cross(DC)).norm() / 2.0;
    const double& sCA = (DC.cross(DA)).norm() / 2.0;
    if ((sAB + sBC + sCA) > (area_ * 1.000001))
    {
        return false;
    }
    return true;
}

void MySimModTriangle::setCircumCenter(
    const std::vector<std::shared_ptr<node>>& nodesGlobal)
{
    const vec3d& a = nodesGlobal[nodes_[0]]->position_;
    const vec3d& b = nodesGlobal[nodes_[1]]->position_;
    const vec3d& c = nodesGlobal[nodes_[2]]->position_;

    const auto ac = c - a;
    const auto ab = b - a;
    const auto abXac = ab.cross(ac);

    const vec3d toCircumsphereCenter = (abXac.cross(ab) * ac.squaredNorm() +
                                        ac.cross(abXac) * ab.squaredNorm()) /
                                       (2.0 * abXac.squaredNorm());
    circumCenter_.radius_ = toCircumsphereCenter.norm();
    circumCenter_.centerPoint_ = a + toCircumsphereCenter;
}

bool MySimModTriangle::hasNode(const UIN oneNode) const
{
    if ((nodes_[0] == oneNode) || (nodes_[1] == oneNode) ||
        (nodes_[2] == oneNode))
    {
        return true;
    }
    return false;
}

bool sameEdge(const std::array<UIN, 2>& edge1, const std::array<UIN, 2>& edge2)
{
    if ((edge1[0] == edge2[0] && edge1[1] == edge2[1]) ||
        (edge1[0] == edge2[1] && edge1[1] == edge2[0]))
    {
        return true;
    }
    return false;
}

}  // namespace MySim
