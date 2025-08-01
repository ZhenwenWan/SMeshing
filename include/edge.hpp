#ifndef H_MYSIM_EDGE
#define H_MYSIM_EDGE

#include "definitions.hpp"
#include "vertex.hpp"

#include <Geom_Curve.hxx>
#include <TopoDS_Edge.hxx>

namespace MySim

{
class edge
{
   public:
    edge(const TopoDS_Edge& occEdge, int tag);

    double parameterFromPoint(const vec3d& point) const;
    vec3d pointFromParameter(double u) const;

    int tag_;
    int geometryType_;
    std::pair<double, double> uBounds_;
    std::map<int, std::shared_ptr<vertex>> vertices_;

    bool isDegenerated() const;
    void seedVertices(double meshSize, int maxVertexTag);

   private:
    Handle(Geom_Curve) occCurve_;
};

}  // namespace MySim
#endif  // H_MYSIM_EDGE
