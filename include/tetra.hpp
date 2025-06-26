#ifndef H_SHONDY_TETRA
#define H_SHONDY_TETRA

#include "definitions.hpp"
#include "sphere.hpp"

#include "node.hpp"
#include <memory>
#include <vector>

namespace shonCloud
{
class tetra
{
   public:
    tetra(UIN a, UIN b, UIN c, UIN d,
          const std::vector<std::shared_ptr<node>>& nodesGlobal);

    std::array<UIN, 4> nodes_;
    std::array<std::array<UIN, 2>, 6> edges_;
    std::array<std::array<UIN, 3>, 4> triangles_;
    Sphere sphere_;
    bool isInSphere(const vec3d& point) const;
    bool hasNode(const UIN oneNode) const;
    const vec3d& center() const;
    vec3d centerQuick(
        const std::vector<std::shared_ptr<node>>& nodesGlobal) const;

    void setSphere(const std::vector<std::shared_ptr<node>>& nodesGlobal);

   private:
};

inline const vec3d& tetra::center() const
{
    return sphere_.centerPoint_;
}

inline vec3d tetra::centerQuick(
    const std::vector<std::shared_ptr<node>>& nodesGlobal) const
{
    const auto& A = nodesGlobal[nodes_[0]]->position_;
    const auto& B = nodesGlobal[nodes_[1]]->position_;
    const auto& C = nodesGlobal[nodes_[2]]->position_;
    const auto& D = nodesGlobal[nodes_[3]]->position_;
    return (A + B + C + D) / 3.0;
}

bool sameTriangle(const std::array<UIN, 3>& triangle1,
                  const std::array<UIN, 3>& triangle2);

}  // namespace shonCloud
#endif  // H_SHONDY_TETRA
