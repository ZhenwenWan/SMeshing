#ifndef H_MYSIM_SOLID
#define H_MYSIM_SOLID

#include "definitions.hpp"
#include "face.hpp"
#include "node.hpp"
#include "MySimModTriangle.hpp"
#include "tetra.hpp"
#include "vertex.hpp"
#include "wire.hpp"
#include <TopoDS_Solid.hxx>

namespace MySim
{
class solid
{
   public:
    solid(const TopoDS_Solid& occSolid, int tag);

    void getAndSetMeshNodes(double spacing, vectorVec3d& newNodes);
    void tetraBowyerWatson();
    // void deleteOutsideNodes();
    // void mesh(const double& meshSize);

    // vec3d pointFromParameter(double u, double v, double w) const;
    // vec2d parameterFromPoint(const vec3d& point) const;
    bool pointInside(const vec3d& parameterPoint) const;

    void writeVTK(const std::string& fileName) const;
    void writeVerticesToVTK(const std::string& fileName) const;

    int tag_;
    int geometryType_;
    std::pair<vec3d, vec3d> uvwBounds_;
    std::pair<vec3d, vec3d> boundingBox_;
    std::vector<std::shared_ptr<MySimModTriangle>> triangles_;
    std::vector<std::shared_ptr<tetra>> tetras_;
    std::vector<std::shared_ptr<node>> nodes_;

    std::map<int, std::shared_ptr<vertex>> vertices_;
    std::map<int, std::shared_ptr<wire>> wires_;
    std::map<int, std::shared_ptr<edge>> embeddedEdges_;
    std::map<int, std::shared_ptr<face>> faces_;

   private:
    Handle(GeomAbs_Shape) occSolid_;
    const TopoDS_Solid occSolidTopo_;
    double geometryTolerance_;
    void setBoundingBox(const TopoDS_Solid& occSolidTopo);
    std::shared_ptr<tetra> getSuperTetra();
    // bool isPointCloseVertex(const vec3d& point, double distThreshold) const;

    // void getFaceAndGeometryNodes(vectorVec3d& nodes) const;
};

}  // namespace MySim
#endif  // H_MYSIM_SOLID
