#ifndef H_SHONDY_SURFACE
#define H_SHONDY_SURFACE

#include "node.hpp"
#include "surfaceBase.hpp"
#include "treeAABB.hpp"

namespace MySim
{
class surface : public surfaceBase
{
   public:
    surface(const std::string& fileName);

    template <typename T>
    void pointsOutside(const std::vector<std::shared_ptr<T>> position,
                       std::vector<UIN>& outside, UIN& newSize) const
    {
#pragma omp parallel for
        for (UIN i = 0; i < position.size(); ++i)
        {
            if (1.0 - 2.0 * treeAABB_->getWindingNumber(position[i]->center()) >
                0)
            {
                outside[i] = 1;
            }
        }
        const UIN numToDelete =
            std::accumulate(outside.begin(), outside.end(), 0);
        newSize = outside.size() - numToDelete;
    }

    double distanceToTriangle(const vec3d& point) const;

   private:
    std::shared_ptr<treeAABB> treeAABB_;
};

}  // namespace MySim
#endif  // H_SHONDY_SURFACE
