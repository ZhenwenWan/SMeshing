#ifndef H_SHONDY_GEOMETRICMODEL
#define H_SHONDY_GEOMETRICMODEL

#include "definitions.hpp"
#include "edge.hpp"
#include "solid.hpp"
#include "vertex.hpp"
#include "wire.hpp"

#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopoDS.hxx>

namespace MySim
{
class geometricModel
{
   public:
    geometricModel();
    void readSTEP(const std::string& fileName);

    std::map<int, std::shared_ptr<solid>> solids_;
    std::map<int, std::shared_ptr<face>> faces_;
    std::map<int, std::shared_ptr<wire>> wires_;
    std::map<int, std::shared_ptr<edge>> edges_;
    std::map<int, std::shared_ptr<vertex>> vertices_;
    std::pair<vec3d, vec3d> boundingBox_;

    int maxVertexTag() const;

   private:
    void setBoundingBox(const TopoDS_Shape& shape);

    void createSolidsFromModel(const TopoDS_Shape& topoShape);

    void createFacesFromSolid(std::shared_ptr<solid> parentSolid,
                              const TopoDS_Solid& topoSolid);

    void createWiresFromFace(std::shared_ptr<face> parentFace,
                             const TopoDS_Face& topoFace);

    void createEdgesFromWire(std::shared_ptr<face> parentFace,
                             std::shared_ptr<wire> parentWire,
                             const TopoDS_Wire& topoWire);

    void createVerticesFromEdge(std::shared_ptr<face> parentFace,
                                std::shared_ptr<edge> parentEdge,
                                const TopoDS_Edge& topoEdge);

    TopTools_DataMapOfShapeInteger occSolids_;
    TopTools_DataMapOfShapeInteger occFaces_;
    TopTools_DataMapOfShapeInteger occWires_;
    TopTools_DataMapOfShapeInteger occEdges_;
    TopTools_DataMapOfShapeInteger occVerices_;
};

}  // namespace MySim
#endif  // H_SHONDY_GEOMETRICMODEL
