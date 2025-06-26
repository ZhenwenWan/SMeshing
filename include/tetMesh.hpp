#ifndef H_SHONDY_TETMESH
#define H_SHONDY_TETMESH

#include "definitions.hpp"

#include "geometricModel.hpp"
#include "node.hpp"
#include "shonMondelerTriangle.hpp"
#include "surface.hpp"

#include <memory>
#include <vector>

namespace shonCloud
{
class tetMesh
{
   public:
    tetMesh();
    void readFromStep(const std::string& fileName);

    void seedVerticesOnEdges(double meshSize);
    void meshFaces(const double& meshSize);
    std::shared_ptr<geometricModel> geometricModel_;

   private:
};

}  // namespace shonCloud
#endif  // H_SHONDY_TETMESH
