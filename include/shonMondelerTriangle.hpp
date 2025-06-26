#ifndef H_SHONDY_TRIANGLE
#define H_SHONDY_TRIANGLE

#include "circle.hpp"
#include "definitions.hpp"

#include "node.hpp"
#include <memory>
#include <vector>

namespace shonCloud
{
class shonMondelerTriangle
{
   public:
    shonMondelerTriangle(UIN a, UIN b, UIN c,
                         const std::vector<std::shared_ptr<node>>& nodesGlobal);

    std::array<UIN, 3> nodes_;
    std::array<std::array<UIN, 2>, 3> edges_;
    Circle circumCenter_;
    bool isInCircumCircle(const vec3d& point) const;
    bool isInside(const vec3d& point) const;
    bool hasNode(const UIN oneNode) const;
    const vec3d& center() const;

    void setCircumCenter(const std::vector<std::shared_ptr<node>>& nodesGlobal);
    vec3d centerQuick(
        const std::vector<std::shared_ptr<node>>& nodesGlobal) const;
    double area_ = 0.0;

   private:
    std::vector<std::shared_ptr<node>> vertices_;
};

inline const vec3d& shonMondelerTriangle::center() const
{
    return circumCenter_.centerPoint_;
}

inline vec3d shonMondelerTriangle::centerQuick(
    const std::vector<std::shared_ptr<node>>& nodesGlobal) const
{
    const auto& A = nodesGlobal[nodes_[0]]->position_;
    const auto& B = nodesGlobal[nodes_[1]]->position_;
    const auto& C = nodesGlobal[nodes_[2]]->position_;
    return (A + B + C) / 3.0;
}

bool sameEdge(const std::array<UIN, 2>& edge1, const std::array<UIN, 2>& edge2);

}  // namespace shonCloud
#endif  // H_SHONDY_TRIANGLE
