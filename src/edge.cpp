#include "edge.hpp"
#include "logger.hpp"
#include "MySimModEnums.hpp"

#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <TopoDS.hxx>

namespace MySim
{
edge::edge(const TopoDS_Edge& occEdge, int tag) : tag_(tag)
{
    TopoDS_Edge occEdgeCopy = occEdge;
    if (occEdgeCopy.Orientation() == TopAbs_INTERNAL ||
        occEdgeCopy.Orientation() == TopAbs_EXTERNAL)
    {
        occEdgeCopy = TopoDS::Edge(occEdgeCopy.Oriented(TopAbs_FORWARD));
    }
    occCurve_ = BRep_Tool::Curve(occEdgeCopy, uBounds_.first, uBounds_.second);
    if (occCurve_.IsNull())
    {
        return;
    }
    if (occCurve_->DynamicType() == STANDARD_TYPE(Geom_Circle))
    {
        geometryType_ = MySimMod::edgeGeometryTypes::circle;
    }
    else if (occCurve_->DynamicType() == STANDARD_TYPE(Geom_Line))
    {
        geometryType_ = MySimMod::edgeGeometryTypes::line;
    }
    else
    {
        std::cout << "ERROR: unknown type" << std::endl;
    }
}

bool edge::isDegenerated() const
{
    if (occCurve_.IsNull())
    {
        return true;
    }
    return false;
}

void edge::seedVertices(double meshSize, int maxVertexTag)
{
    const double length = pointFromParameter(uBounds_.second)
                              .distance(pointFromParameter(uBounds_.first));

    int numberOfNodes = (int)round(length / meshSize);
    if (geometryType_ == MySimMod::edgeGeometryTypes::circle)
    {
        double radius = Handle(Geom_Circle)::DownCast(occCurve_)->Radius();
        double circleLength = radius * (uBounds_.second - uBounds_.first);
        numberOfNodes = (int)round(circleLength / meshSize);
    }

    for (int i = 1; i < numberOfNodes; ++i)
    {
        const double uNode =
            i * (uBounds_.second - uBounds_.first) / numberOfNodes +
            uBounds_.first;
        const vec3d position = pointFromParameter(uNode);
        ++maxVertexTag;
        auto newVertex = std::make_shared<vertex>(
            position, MySimMod::vertexTypes::onEdge);
        newVertex->tag_ = maxVertexTag;
        vertices_.insert({maxVertexTag, newVertex});
    }
}

double edge::parameterFromPoint(const vec3d& point) const
{
    // add little space to boundaries
    double umin = uBounds_.first;
    double umax = uBounds_.second;
    const double du = umax - umin;
    const double utol = std::max(fabs(du) * 1e-8, 1e-12);
    umin -= utol;
    umax += utol;
    gp_Pnt pnt(point[0], point[1], point[2]);
    GeomAPI_ProjectPointOnCurve proj(pnt, occCurve_, umin, umax);
    return proj.LowerDistanceParameter();
}

vec3d edge::pointFromParameter(double u) const
{
    gp_Pnt pnt = occCurve_->Value(u);
    return vec3d(pnt.X(), pnt.Y(), pnt.Z());
}

}  // namespace MySim
