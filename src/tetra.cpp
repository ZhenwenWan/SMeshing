#include "tetra.hpp"
#include "shonMondelerTriangle.hpp"

namespace shonCloud
{
tetra::tetra(UIN a, UIN b, UIN c, UIN d,
             const std::vector<std::shared_ptr<node>>& nodesGlobal)
    : nodes_({a, b, c, d}), edges_(), triangles_()
{
    edges_[0][0] = a;
    edges_[0][1] = b;
    edges_[1][0] = b;
    edges_[1][1] = c;
    edges_[2][0] = c;
    edges_[2][1] = a;
    edges_[3][0] = d;
    edges_[3][1] = b;
    edges_[4][0] = d;
    edges_[4][1] = c;
    edges_[5][0] = d;
    edges_[5][1] = a;
    triangles_[0][0] = a;
    triangles_[0][1] = b;
    triangles_[0][2] = c;
    triangles_[1][0] = b;
    triangles_[1][1] = d;
    triangles_[1][2] = c;
    triangles_[2][0] = c;
    triangles_[2][1] = d;
    triangles_[2][2] = a;
    triangles_[3][0] = d;
    triangles_[3][1] = b;
    triangles_[3][2] = a;
    setSphere(nodesGlobal);
}

void tetra::setSphere(const std::vector<std::shared_ptr<node>>& nodesGlobal)
{
    const vec3d& pA = nodesGlobal[nodes_[0]]->position_;
    const vec3d& pB = nodesGlobal[nodes_[1]]->position_;
    const vec3d& pC = nodesGlobal[nodes_[2]]->position_;
    const vec3d& pD = nodesGlobal[nodes_[3]]->position_;

    const auto volume = ((pB - pA).cross(pC - pA)).dot(pD - pA) / 6;
    if (std::abs(volume) < 1.0e-30)
    {
        sphere_.radius_ = 0.0;
        sphere_.centerPoint_ = {0, 0, 0};
        return;
    }
    const auto a = (pB - pA).norm();
    const auto b = (pC - pA).norm();
    const auto c = (pD - pA).norm();
    const auto A = (pD - pC).norm();
    const auto B = (pD - pB).norm();
    const auto C = (pB - pC).norm();
    const auto& xRadius =
        std::sqrt((a * A + b * B + c * C) * (a * A + b * B - c * C) *
                  (a * A - b * B + c * C) * (-a * A + b * B + c * C)) /
        volume / 24;
    if (xRadius < 1.0e-10)
    {
        sphere_.radius_ = 0.0;
        sphere_.centerPoint_ = {0, 0, 0};
        return;
    }
    auto xTriangle = std::make_shared<shonMondelerTriangle>(
        nodes_[0], nodes_[1], nodes_[2], nodesGlobal);
    const vec3d& pO_ABC = xTriangle->circumCenter_.centerPoint_;
    double height = std::sqrt(xRadius * xRadius - pO_ABC.squareDist(pA));
    vec3d vectorNorm_ABC = ((pB - pA).cross(pC - pA)).normalized();
    sphere_.radius_ = xRadius;
    if (xRadius > pO_ABC.distance(pD))
    {
        sphere_.centerPoint_ = pO_ABC + vectorNorm_ABC * height;
    }
    else
    {
        sphere_.centerPoint_ = pO_ABC - vectorNorm_ABC * height;
    }
}

bool tetra::isInSphere(const vec3d& point) const
{
    // std::cout << "  SphereCenter  " << sphere_.centerPoint_ << std::endl;
    // std::cout << "  SphereRadius  " << sphere_.radius_ << std::endl;
    // std::cout << "  visitingPoint " << point << std::endl;
    if ((distance(point, sphere_.centerPoint_) < sphere_.radius_) &&
        (sphere_.radius_ > 1.0e-10))
    {
        return true;
    }
    return false;
}

bool tetra::hasNode(const UIN oneNode) const
{
    if ((nodes_[0] == oneNode) || (nodes_[1] == oneNode) ||
        (nodes_[2] == oneNode) || (nodes_[3] == oneNode))
    {
        return true;
    }
    return false;
}

bool sameTriangle(const std::array<UIN, 3>& triangle1,
                  const std::array<UIN, 3>& triangle2)
{
    for (UIN i = 0; i < 3; i++)
    {
        if (!(triangle1[i] == triangle2[0] || triangle1[i] == triangle2[1] ||
              triangle1[i] == triangle2[2]))
            return false;
        if (!(triangle2[i] == triangle1[0] || triangle2[i] == triangle1[1] ||
              triangle2[i] == triangle1[2]))
            return false;
    }
    return true;
}

}  // namespace shonCloud
