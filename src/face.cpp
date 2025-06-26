#include "face.hpp"
#include "edge.hpp"
#include "logger.hpp"
#include "shonModelerEnums.hpp"
#include "vectorUtil.hpp"

#include <limits>

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <StdFail_NotDone.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Sphere.hxx>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace shonCloud
{
face::face(const TopoDS_Face& occFace, int tag)
    : tag_(tag), occFaceTopo_(occFace)
{
    geometryTolerance_ = BRep_Tool::Tolerance(occFace);
    setBoundingBox(occFace);
    occFace_ = BRep_Tool::Surface(occFace);

    if (occFace_->DynamicType() == STANDARD_TYPE(Geom_Plane))
    {
        geometryType_ = shonModeler::faceGeometryTypes::plane;
    }
    else if (occFace_->DynamicType() == STANDARD_TYPE(Geom_SphericalSurface))
    {
        geometryType_ = shonModeler::faceGeometryTypes::sphere;
    }
    else
    {
        std::cout << "ERROR: face type not known" << std::endl;
        loggingError("face type not known");
    }
}

void face::setBoundingBox(const TopoDS_Face& occFace)
{
    ShapeAnalysis::GetFaceUVBounds(occFace, uvBounds_.first[0],
                                   uvBounds_.second[0], uvBounds_.first[1],
                                   uvBounds_.second[1]);
    Bnd_Box b;
    BRepBndLib::Add(occFace, b);
    b.Get(boundingBox_.first[0], boundingBox_.first[1], boundingBox_.first[2],
          boundingBox_.second[0], boundingBox_.second[1],
          boundingBox_.second[2]);
}

void face::getEdgeAndGeometryNodes(vectorVec3d& nodes) const
{
    for (const auto& oneWire : wires_)
    {
        for (const auto& oneEdge : oneWire.second->edges_)
        {
            for (const auto& oneVertex : oneEdge.second->vertices_)
            {
                if (oneVertex.second->type_ == shonModeler::vertexTypes::onEdge)
                {
                    nodes.push_back(oneVertex.second->position_);
                }
            }
        }
    }

    for (const auto& oneGeoVertex : vertices_)
    {
        nodes.push_back(oneGeoVertex.second->position_);
    }
}

void face::getMeshNodes(const double& meshSize, vectorVec3d& nodes)
{
    getEdgeAndGeometryNodes(nodes);
    vectorVec3d innerNodes;
    seedNodesInFace(meshSize, innerNodes);
    std::vector<UIN> toDelete(innerNodes.size(), 0);
    for (const auto& oneFixedNode : nodes)
    {
        for (UIN i = 0; i < innerNodes.size(); ++i)
        {
            if (distance(oneFixedNode, innerNodes[i]) < meshSize)
            {
                toDelete[i] = 1;
            }
        }
    }

    const UIN newNodesSize =
        innerNodes.size() -
        std::accumulate(toDelete.begin(), toDelete.end(), 0);

    deleteByIndex(toDelete, innerNodes, newNodesSize);
    nodes.insert(nodes.end(), innerNodes.begin(), innerNodes.end());
    nodes_.clear();
    for (const auto& oneNode : nodes)
    {
        vec2d uv = parameterFromPoint(oneNode);
        vec3d xyz = {uv[0], uv[1], 0};
        nodes_.push_back(std::make_shared<node>(xyz));
    }
}

void face::mesh(const double& meshSize)
{
    vectorVec3d innerNodes;
    getMeshNodes(meshSize, innerNodes);

    // triangulateBowyerWatson();
    triangulate();

    for (auto& xNode : nodes_)
    {
        vec2d uv;
        uv[0] = xNode->position_[0];
        uv[1] = xNode->position_[1];
        xNode->position_ = pointFromParameter(uv[0], uv[1]);
    }
}

vec3d face::pointFromParameter(double u, double v) const
{
    gp_Pnt val = occFace_->Value(u, v);
    return vec3d(val.X(), val.Y(), val.Z());
}

vec2d face::parameterFromPoint(const vec3d& point) const
{
    double umin = uvBounds_.first[0];
    double vmin = uvBounds_.first[1];
    double umax = uvBounds_.second[0];
    double vmax = uvBounds_.second[1];
    const double du = uvBounds_.second[0] - uvBounds_.first[0];
    const double utol = std::max(fabs(du) * 1e-8, 1e-12);
    umin -= utol;
    umax += utol;
    const double dv = uvBounds_.second[1] - uvBounds_.first[1];
    const double vtol = std::max(fabs(dv) * 1e-8, 1e-12);
    vmin -= vtol;
    vmax += vtol;
    gp_Pnt pnt(point[0], point[1], point[2]);
    GeomAPI_ProjectPointOnSurf proj(pnt, occFace_, umin, umax, vmin, vmax);
    if (proj.NbPoints() == 0)
    {
        ShapeAnalysis_Surface sas(occFace_);
        // get UV of point on surface
        gp_Pnt2d uv = sas.ValueOfUV(pnt, 0.01);
        vec2d par(uv.X(), uv.Y());
        return par;
    }
    vec2d par;
    try
    {
        proj.LowerDistanceParameters(par[0], par[1]);
    }
    catch (...)
    {
        std::cout << "ERROR: getting parameter on face failed" << std::endl;
        par[0] = std::numeric_limits<double>::quiet_NaN();
        par[1] = std::numeric_limits<double>::quiet_NaN();
    }

    return par;
}

void face::seedNodesInFace(double meshSize, vectorVec3d& newNodes)
{
    const auto& origin = uvBounds_.first;
    const auto& end = uvBounds_.second;
    const vec3d& point0 = pointFromParameter(origin[0], origin[1]);
    const vec3d& point1 = pointFromParameter(end[0], origin[1]);
    const vec3d& point2 = pointFromParameter(origin[0], end[1]);
    double distX = point0.distance(point1);
    double distY = point0.distance(point2);
    int noOfparticlesInX = (int)round(distX / meshSize);
    int noOfparticlesInY =
        (int)round(distY / (meshSize * std::sqrt(3.0) / 2.0));
    double spacing = distX / noOfparticlesInX;

    if (geometryType_ == shonModeler::faceGeometryTypes::sphere)
    {
        BRepAdaptor_Surface surface(occFaceTopo_);
        gp_Sphere sphere = surface.Sphere();
        const double radius = sphere.Radius();
        distX = radius * (end[0] - origin[0]);
        noOfparticlesInX = (int)round(distX / meshSize);
        spacing = (end[0] - origin[0]) / noOfparticlesInX;
    }

    if (geometryType_ == shonModeler::faceGeometryTypes::sphere)
    {
        const int& nodesOfShpere = (int)round(3.37 * 3.37 / spacing / spacing);
        for (int n = 0; n < nodesOfShpere; n++)
        {
            double v =
                std::acos(1.0 - 2.0 * (n + 3.5) / (nodesOfShpere + 6.0)) -
                0.5 * PI;
            double u =
                std::remainder(4.0 * PI / (1.0 + std::sqrt(5.0)) * n, 2.0 * PI);
            if (u < 0.0) u += 2.0 * PI;
            if (u < origin[0] + spacing || u > end[0] - spacing) continue;
            if (v < origin[1] + spacing || v > end[1] - spacing) continue;
            vec2d newNode(u, v);
            if (pointInside(newNode))
            {
                const vec3d newNode3D =
                    pointFromParameter(newNode[0], newNode[1]);
                newNodes.push_back(newNode3D);
            }
        }
    }
    else
    {
        const double particleRadius = spacing / 2.0;
        double staggeredDeltaX = 0.0;
        for (int x = 0; x < noOfparticlesInX; ++x)
        {
            for (int y = 0; y < noOfparticlesInY; ++y)
            {
                if (std::abs(y % 2) == 1)
                {
                    staggeredDeltaX = particleRadius;
                }
                else
                {
                    staggeredDeltaX = 0.0;
                }
                double u =
                    origin[0] + particleRadius + staggeredDeltaX + x * spacing;
                double v = origin[1] + particleRadius +
                           y * spacing * std::sqrt(3.0) / 2.0;
                vec2d newNode(u, v);
                if (pointInside(newNode))
                {
                    const vec3d newNode3D =
                        pointFromParameter(newNode[0], newNode[1]);
                    newNodes.push_back(newNode3D);
                }
            }
        }
    }
}

bool face::isPointCloseVertex(const vec3d& point, double distThreshold) const
{
    for (const auto& oneVertex : vertices_)
    {
        if (distance(point, oneVertex.second->position_) < distThreshold)
        {
            return true;
        }
    }
    return false;
}

bool face::pointInside(const vec2d& parameterPoint) const
{
    BRepClass_FaceClassifier faceClassifier;
    faceClassifier.Perform(occFaceTopo_,
                           gp_Pnt2d{parameterPoint[0], parameterPoint[1]},
                           geometryTolerance_);

    const TopAbs_State state = faceClassifier.State();
    return (state == TopAbs_IN || state == TopAbs_ON);
}

void face::triangulate()
{
    triangles_.clear();
    const auto& superTriangle = getSuperTriangle();
    triangles_.push_back(superTriangle);

    std::vector<UIN> toDelete;
    toDelete.push_back(0);

    for (UIN nodeID = 0; nodeID < nodes_.size(); ++nodeID)
    {
        const auto& oneNode = nodes_[nodeID];
        int triangleID = -1;
        int xID = -1;
        for (const auto& oneTriangle : triangles_)
        {
            xID++;
            if (toDelete[xID] == 1) continue;
            if (oneTriangle->isInside(oneNode->position_))
            {
                triangleID = xID;
                break;
            }
        }
        if (triangleID < 0)
        {
            std::cout << "ERROR: superTriangle can NOT house the node "
                      << nodes_[nodeID]->position_ << std::endl;
            loggingError("superTriangle failed to house at least one node");
            for (const auto& oneTriangle : triangles_)
            {
                std::cout << "0000 "
                          << nodes_[oneTriangle->nodes_[0]]->position_
                          << std::endl;
                std::cout << "1111 "
                          << nodes_[oneTriangle->nodes_[1]]->position_
                          << std::endl;
                std::cout << "2222 "
                          << nodes_[oneTriangle->nodes_[2]]->position_
                          << std::endl;
                std::cout << "3333 "
                          << oneTriangle->circumCenter_.centerPoint_.distance(
                                 nodes_[nodeID]->position_)
                          << std::endl;
            }
            return;
        }
        else
        {
            toDelete[triangleID] = 1;
            const auto& oneTriangle = triangles_[triangleID];

            for (const auto& oneEdge : oneTriangle->edges_)
            {
                auto newTriangle = std::make_shared<shonMondelerTriangle>(
                    oneEdge[0], oneEdge[1], nodeID, nodes_);
                if (newTriangle->area_ > EMINUS13)
                {
                    triangles_.push_back(newTriangle);
                    toDelete.push_back(0);
                }
            }
        }
        std::cout << "Node" << nodeID << "  triangles" << triangles_.size()
                  << std::endl;
    }

    for (UIN i = 0; i < triangles_.size(); ++i)
    {
        if (toDelete[i] == 1) continue;
        const auto& oneTriangle = triangles_[i];

        for (int j = 0; j < 3; ++j)
        {
            if (oneTriangle->hasNode(superTriangle->nodes_[j]))
            {
                toDelete[i] = 1;
            }
        }
    }

    for (UIN i = 0; i < triangles_.size(); ++i)
    {
        if (toDelete[i] == 1) continue;
        const auto& oneTriangle = triangles_[i];
        const auto center = oneTriangle->centerQuick(nodes_);
        if (!pointInside({center[0], center[1]}))
        {
            toDelete[i] = 1;
        }
    }

    const int sumDelete = std::accumulate(toDelete.begin(), toDelete.end(), 0);
    const UIN newSize = triangles_.size() - sumDelete;
    deleteByIndex(toDelete, triangles_, newSize);
    //    delete supernodes
    nodes_.resize(nodes_.size() - 3);
}

void face::triangulateBowyerWatson()
{
    triangles_.clear();
    const auto& superTriangle = getSuperTriangle();
    triangles_.push_back(superTriangle);

    for (UIN nodeID = 0; nodeID < nodes_.size(); ++nodeID)
    {
        const auto& oneNode = nodes_[nodeID];
        std::vector<std::shared_ptr<shonMondelerTriangle>> badTriangles;
        for (const auto& oneTriangle : triangles_)
        {
            if (oneTriangle->isInCircumCircle(oneNode->position_))
            {
                badTriangles.push_back(oneTriangle);
            }
        }
        // find boundary of polygon hole

        std::vector<std::array<UIN, 2>> polygon;
        for (UIN i = 0; i < badTriangles.size(); ++i)
        {
            const auto& oneBadTriangle = badTriangles[i];
            // TODO find out if a edge is shared
            for (const auto& oneEdge : oneBadTriangle->edges_)
            {
                bool isShared = false;
                for (UIN j = 0; j < badTriangles.size(); ++j)
                {
                    const auto& otherBadTriangle = badTriangles[j];
                    if (i == j)
                    {
                        continue;
                    }
                    for (const auto& otherEdge : otherBadTriangle->edges_)
                    {
                        if (sameEdge(oneEdge, otherEdge))
                        {
                            isShared = true;
                        }
                    }
                }
                if (isShared == false)
                {
                    polygon.push_back(oneEdge);
                }
            }
        }
        for (const auto& oneBadTriangle : badTriangles)
        {
            // remove oneBadTriangle from triangulation
            auto it = std::remove(triangles_.begin(), triangles_.end(),
                                  oneBadTriangle);
            triangles_.erase(it, triangles_.end());
        }

        for (const auto& oneEdge : polygon)
        {
            auto newTriangle = std::make_shared<shonMondelerTriangle>(
                oneEdge[0], oneEdge[1], nodeID, nodes_);
            triangles_.push_back(newTriangle);
        }
    }

    std::vector<UIN> toDelete(triangles_.size(), 0);
    for (UIN i = 0; i < triangles_.size(); ++i)
    {
        const auto& oneTriangle = triangles_[i];

        for (int j = 0; j < 3; ++j)
        {
            if (oneTriangle->hasNode(superTriangle->nodes_[j]))
            {
                toDelete[i] = 1;
            }
        }
    }

    for (UIN i = 0; i < triangles_.size(); ++i)
    {
        const auto& oneTriangle = triangles_[i];
        const auto center = oneTriangle->centerQuick(nodes_);
        if (!pointInside({center[0], center[1]}) && toDelete[i] == 0)
        {
            toDelete[i] = 1;
        }
    }

    const int sumDelete = std::accumulate(toDelete.begin(), toDelete.end(), 0);
    const UIN newSize = triangles_.size() - sumDelete;
    deleteByIndex(toDelete, triangles_, newSize);
    //    delete supernodes
    nodes_.resize(nodes_.size() - 3);
}

std::shared_ptr<shonMondelerTriangle> face::getSuperTriangle()
{
    vec3d min = nodes_[0]->position_;
    vec3d max = nodes_[0]->position_;

    for (const auto& oneNode : nodes_)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (oneNode->position_[i] > max[i])
            {
                max[i] = oneNode->position_[i];
            }
            if (oneNode->position_[i] < min[i])
            {
                min[i] = oneNode->position_[i];
            }
        }
    }

    const vec3d delta = max - min;
    const double deltaMax = std::max(delta[0], delta[1]);
    const double midx = (min[0] + max[0]) / 2.0;
    const double midy = (min[1] + max[1]) / 2.0;

    auto a = std::make_shared<node>(
        vec3d(midx - 20 * deltaMax, midy - 20 * deltaMax, 0));
    auto b = std::make_shared<node>(
        vec3d(midx + 20 * deltaMax, midy - 20 * deltaMax, 0));
    auto c = std::make_shared<node>(vec3d(midx, midy + 100 * deltaMax, 0));
    nodes_.push_back(a);
    nodes_.push_back(b);
    nodes_.push_back(c);
    UIN numNodes = nodes_.size();
    auto superTriangle = std::make_shared<shonMondelerTriangle>(
        numNodes - 3, numNodes - 2, numNodes - 1, nodes_);
    return superTriangle;
}

void face::writeVTK(const std::string& fileName) const
{
    vtkNew<vtkPoints> points;
    UIN nodeIDCounter = 0;
    for (auto& oneNode : nodes_)
    {
        points->InsertNextPoint(oneNode->position_.data());
        oneNode->id_ = nodeIDCounter;
        ++nodeIDCounter;
    }
    vtkNew<vtkCellArray> triangles;
    for (UIN i = 0; i < triangles_.size(); ++i)
    {
        const auto& oneTriangle = triangles_[i];
        vtkNew<vtkTriangle> triangle;
        triangle->GetPointIds()->SetId(0, nodes_[oneTriangle->nodes_[0]]->id_);
        triangle->GetPointIds()->SetId(1, nodes_[oneTriangle->nodes_[1]]->id_);
        triangle->GetPointIds()->SetId(2, nodes_[oneTriangle->nodes_[2]]->id_);
        triangles->InsertNextCell(triangle);
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void face::writeVerticesToVTK(const std::string& fileName) const
{
    vtkNew<vtkPoints> points;
    vtkNew<vtkIntArray> vertexType;
    for (auto& oneVertex : vertices_)
    {
        points->InsertNextPoint(oneVertex.second->position_.data());
        vertexType->InsertNextValue(oneVertex.second->type_);
    }
    vertexType->SetName("vertexType");
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->GetPointData()->AddArray(vertexType);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}
}  // namespace shonCloud
