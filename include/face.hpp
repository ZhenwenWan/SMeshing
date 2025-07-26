#ifndef H_MYSIM_FACE
#define H_MYSIM_FACE

#include "definitions.hpp"
#include "node.hpp"
#include "MySimModTriangle.hpp"
#include "vertex.hpp"
#include "wire.hpp"
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>

namespace MySim
{
class face
{
   public:
    face(const TopoDS_Face& occFace, int tag);

    void seedNodesInFace(double spacing, vectorVec3d& newNodes);
    void triangulateBowyerWatson();
    void triangulate();
    void deleteOutsideNodes();
    void getMeshNodes(const double& spacing, vectorVec3d& newNodes);
    void mesh(const double& meshSize);

    vec3d pointFromParameter(double u, double v) const;
    vec2d parameterFromPoint(const vec3d& point) const;
    bool pointInside(const vec2d& parameterPoint) const;

    void writeVTK(const std::string& fileName) const;
    void writeVerticesToVTK(const std::string& fileName) const;

    int tag_;
    int geometryType_;
    std::pair<vec2d, vec2d> uvBounds_;
    std::pair<vec3d, vec3d> boundingBox_;
    std::vector<std::shared_ptr<MySimModTriangle>> triangles_;
    std::vector<std::shared_ptr<node>> nodes_;

    std::map<int, std::shared_ptr<vertex>> vertices_;
    std::map<int, std::shared_ptr<wire>> wires_;
    std::map<int, std::shared_ptr<edge>> embeddedEdges_;

   private:
    Handle(Geom_Surface) occFace_;
    const TopoDS_Face occFaceTopo_;
    double geometryTolerance_;
    void setBoundingBox(const TopoDS_Face& occFace);
    std::shared_ptr<MySimModTriangle> getSuperTriangle();
    bool isPointCloseVertex(const vec3d& point, double distThreshold) const;

    void getEdgeAndGeometryNodes(vectorVec3d& nodes) const;
};

}  // namespace MySim
#endif  // H_MYSIM_FACE
