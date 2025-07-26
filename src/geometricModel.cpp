#include "geometricModel.hpp"

#include <BRepBndLib.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <STEPControl_Reader.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>

namespace MySim
{
geometricModel::geometricModel()
{
}

void geometricModel::setBoundingBox(const TopoDS_Shape& shape)
{
    Bnd_Box b;
    BRepBndLib::Add(shape, b);
    b.Get(boundingBox_.first[0], boundingBox_.first[1], boundingBox_.first[2],
          boundingBox_.second[0], boundingBox_.second[1],
          boundingBox_.second[2]);
}

void geometricModel::readSTEP(const std::string& fileName)
{
    STEPControl_Reader reader;
    reader.ReadFile(fileName.c_str());
    reader.NbRootsForTransfer();
    reader.TransferRoots();

    TopoDS_Shape occShape;
    occShape = reader.OneShape();
    BRepTools::Clean(occShape);
    setBoundingBox(occShape);

    createSolidsFromModel(occShape);
}

void geometricModel::createSolidsFromModel(const TopoDS_Shape& topoShape)
{
    TopExp_Explorer explopre;

    for (explopre.Init(topoShape, TopAbs_SOLID); explopre.More();
         explopre.Next())
    {
        TopoDS_Solid occSolid = TopoDS::Solid(explopre.Current());

        if (!occSolids_.IsBound(occSolid))
        {
            const int oneSolidID = solids_.size();
            occSolids_.Bind(occSolid, oneSolidID);
            auto newSolid = std::make_shared<solid>(occSolid,oneSolidID);
            solids_.insert({oneSolidID, newSolid});

            createFacesFromSolid(newSolid, occSolid);
        }
        else
        {
            const int oneSolidID = occSolids_.Find(occSolid);
            auto oneSolid = solids_.at(oneSolidID);

            // why solid act different?
            createFacesFromSolid(oneSolid, occSolid);
        }
    }
}

void geometricModel::createFacesFromSolid(std::shared_ptr<solid> parentSolid,
                                          const TopoDS_Solid& topoSolid)
{
    TopExp_Explorer explopre;

    for (explopre.Init(topoSolid, TopAbs_FACE); explopre.More();
         explopre.Next())
    {
        TopoDS_Face occFace = TopoDS::Face(explopre.Current());

        if (!occFaces_.IsBound(occFace))
        {
            const int oneFaceID = faces_.size();
            occFaces_.Bind(occFace, oneFaceID);
            auto newFace = std::make_shared<face>(occFace, oneFaceID);

            faces_.insert({oneFaceID, newFace});
            parentSolid->faces_.insert({oneFaceID, newFace});

            createWiresFromFace(newFace, occFace);
        }
        else
        {
            const int oneFaceID = occFaces_.Find(occFace);
            auto oneFace = faces_.at(oneFaceID);

            if (parentSolid->faces_.count(oneFaceID) == 0)
            {
                parentSolid->faces_.insert({oneFaceID, oneFace});
            }
        }
    }
}

void geometricModel::createWiresFromFace(std::shared_ptr<face> parentFace,
                                         const TopoDS_Face& topoFace)
{
    TopExp_Explorer explopre;

    for (explopre.Init(topoFace.Oriented(TopAbs_FORWARD), TopAbs_WIRE);
         explopre.More(); explopre.Next())
    {
        TopoDS_Wire occWire = TopoDS::Wire(explopre.Current());

        if (!occWires_.IsBound(occWire))
        {
            const int oneWireTID = wires_.size();
            occWires_.Bind(occWire, oneWireTID);
            auto newWire = std::make_shared<wire>(oneWireTID);

            wires_.insert({oneWireTID, newWire});
            parentFace->wires_.insert({oneWireTID, newWire});

            createEdgesFromWire(parentFace, newWire, occWire);
        }
        else
        {
            const int oneWireID = occWires_.Find(occWire);
            auto oneWire = wires_.at(oneWireID);

            if (parentFace->wires_.count(oneWireID) == 0)
            {
                parentFace->wires_.insert({oneWireID, oneWire});
                // insert vertices from an existing wire to faces
                for (const auto& oneEdge : oneWire->edges_)
                {
                    for (const auto& oneVertex : oneEdge.second->vertices_)
                    {
                        if (parentFace->vertices_.count(oneVertex.first) == 0)
                        {
                            parentFace->vertices_.insert(
                                {oneVertex.second->tag_, oneVertex.second});
                        }
                    }
                }
            }
        }
    }
}

void geometricModel::createEdgesFromWire(std::shared_ptr<face> parentFace,
                                         std::shared_ptr<wire> parentWire,
                                         const TopoDS_Wire& topoWire)
{
    TopExp_Explorer explopre;

    for (explopre.Init(topoWire, TopAbs_EDGE); explopre.More(); explopre.Next())
    {
        TopoDS_Edge occEdge = TopoDS::Edge(explopre.Current());

        if (!occEdges_.IsBound(occEdge))
        {
            const int oneEdgeID = edges_.size();
            auto newEdge = std::make_shared<edge>(occEdge, oneEdgeID);

            if (newEdge->isDegenerated()) continue;

            occEdges_.Bind(occEdge, oneEdgeID);
            edges_.insert({oneEdgeID, newEdge});

            createVerticesFromEdge(parentFace, newEdge, occEdge);
        }

        const int oneEdgeID = occEdges_.Find(occEdge);
        auto oneEdge = edges_.at(oneEdgeID);

        if (occEdge.Orientation() == TopAbs_INTERNAL)
        {
            if (parentFace->embeddedEdges_.count(oneEdgeID) == 0)
            {
                parentFace->embeddedEdges_.insert({oneEdgeID, oneEdge});
                for (const auto& oneVertex : oneEdge->vertices_)
                {
                    if (parentFace->vertices_.count(oneVertex.second->tag_) ==
                        0)
                    {
                        parentFace->vertices_.insert(
                            {oneVertex.first, oneVertex.second});
                    }
                }
            }
        }
        else
        {
            if (parentWire->edges_.count(oneEdgeID) == 0)
            {
                parentWire->edges_.insert({oneEdgeID, oneEdge});
                for (const auto& oneVertex : oneEdge->vertices_)
                {
                    if (parentFace->vertices_.count(oneVertex.second->tag_) ==
                        0)
                    {
                        parentFace->vertices_.insert(
                            {oneVertex.first, oneVertex.second});
                    }
                }
            }
        }
    }
}

void geometricModel::createVerticesFromEdge(std::shared_ptr<face> parentFace,
                                            std::shared_ptr<edge> parentEdge,
                                            const TopoDS_Edge& topoEdge)
{
    TopExp_Explorer explopre;
    for (explopre.Init(topoEdge, TopAbs_VERTEX); explopre.More();
         explopre.Next())
    {
        TopoDS_Vertex occVertex = TopoDS::Vertex(explopre.Current());

        if (!occVerices_.IsBound(occVertex))
        {
            const int oneVertexID = vertices_.size();
            occVerices_.Bind(occVertex, oneVertexID);
            auto newVertex = std::make_shared<vertex>(occVertex, oneVertexID);

            vertices_.insert({oneVertexID, newVertex});
            parentEdge->vertices_.insert({oneVertexID, newVertex});
            parentFace->vertices_.insert({oneVertexID, newVertex});
        }
        else
        {
            const int oneVertexID = occVerices_.Find(occVertex);
            auto oneVertex = vertices_.at(oneVertexID);

            if (parentEdge->vertices_.count(oneVertexID) == 0)
            {
                parentEdge->vertices_.insert({oneVertexID, oneVertex});
            }
            if (parentFace->vertices_.count(oneVertexID) == 0)
            {
                parentFace->vertices_.insert({oneVertexID, oneVertex});
            }
        }
    }
}

int geometricModel::maxVertexTag() const
{
    int maxTag = 0;
    for (const auto& oneVertex : vertices_)
    {
        if (oneVertex.first > maxTag)
        {
            maxTag = oneVertex.first;
        }
    }
    for (const auto& oneEdge : edges_)
    {
        for (const auto& oneVertex : oneEdge.second->vertices_)
        {
            if (oneVertex.first > maxTag)
            {
                maxTag = oneVertex.first;
            }
        }
    }
    return maxTag;
}

}  // namespace MySim
