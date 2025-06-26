#include "solid.hpp"
#include "shonModelerEnums.hpp"
#include "vectorUtil.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <ShapeAnalysis.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Pnt.hxx>
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
solid::solid(const TopoDS_Solid& occSolid, int tag)
    : tag_(tag), occSolidTopo_(occSolid)
{
    setBoundingBox(occSolidTopo_);
}

void solid::setBoundingBox(const TopoDS_Solid& occSolidTopo)
{
    Bnd_Box b;
    BRepBndLib::Add(occSolidTopo, b);
    b.Get(boundingBox_.first[0], boundingBox_.first[1], boundingBox_.first[2],
          boundingBox_.second[0], boundingBox_.second[1],
          boundingBox_.second[2]);
}

void solid::getAndSetMeshNodes(double spacing, vectorVec3d& newNodes)
{
    for (auto& oneFace : faces_)
    {
        vectorVec3d nodes;
        oneFace.second->getMeshNodes(spacing, nodes);

        for (const auto& oneNode : nodes)
        {
            newNodes.push_back(oneNode);
            nodes_.push_back(std::make_shared<node>(oneNode));
            auto newVertex = std::make_shared<vertex>(
                oneNode, shonModeler::vertexTypes::onEdge);
            newVertex->tag_ = vertices_.size();
            vertices_.insert({newVertex->tag_, newVertex});
        }
    }

    const auto& origin = boundingBox_.first;
    const auto& end = boundingBox_.second;
    const double distX = (end[0] - origin[0]);
    const double distY = (end[1] - origin[1]);
    const double distZ = (end[2] - origin[2]);
    int noOfparticlesInX = (int)round(distX / spacing);
    int noOfparticlesInY = (int)round(distY / spacing);
    int noOfparticlesInZ = (int)round(distZ / spacing);

    for (int x = 0; x < noOfparticlesInX; ++x)
    {
        for (int y = 0; y < noOfparticlesInY; ++y)
        {
            for (int z = 0; z < noOfparticlesInZ; ++z)
            {
                vec3d xNode(origin[0] + x * spacing, origin[1] + y * spacing,
                            origin[2] + z * spacing);
                if (pointInside(xNode))
                {
                    newNodes.push_back(xNode);
                    nodes_.push_back(std::make_shared<node>(xNode));
                    auto newVertex = std::make_shared<vertex>(
                        xNode, shonModeler::vertexTypes::onEdge);
                    newVertex->tag_ = vertices_.size();
                    vertices_.insert({newVertex->tag_, newVertex});
                }
            }
        }
    }
}

bool solid::pointInside(const vec3d& parameterPoint) const
{
    BRepClass3d_SolidClassifier solidClassifier(occSolidTopo_);
    solidClassifier.Perform(
        gp_Pnt{parameterPoint[0], parameterPoint[1], parameterPoint[2]},
        geometryTolerance_);

    const TopAbs_State state = solidClassifier.State();
    return (state == TopAbs_IN || state == TopAbs_ON);
}

void solid::tetraBowyerWatson()
{
    tetras_.clear();
    const auto& superTetra = getSuperTetra();
    tetras_.push_back(superTetra);

    for (UIN nodeID = 0; nodeID < nodes_.size(); ++nodeID)
    // for (UIN nodeID = 0; nodeID < 100; ++nodeID)
    {
        const auto& oneNode = nodes_[nodeID];
        std::vector<std::shared_ptr<tetra>> housingTetras;
        housingTetras.clear();
        for (const auto& oneTetra : tetras_)
        {
            if (oneTetra->isInSphere(oneNode->position_))
            {
                housingTetras.push_back(oneTetra);
            }
        }

        std::vector<std::array<UIN, 3>> polyhedra;
        for (UIN i = 0; i < housingTetras.size(); ++i)
        {
            const auto& oneTetra = housingTetras[i];
            for (const auto& oneTriangle : oneTetra->triangles_)
            {
                bool isShared = false;
                for (UIN j = 0; j < housingTetras.size(); ++j)
                {
                    const auto& anotherTetra = housingTetras[j];
                    if (i == j)
                    {
                        continue;
                    }
                    for (const auto& anotherTriangle : anotherTetra->triangles_)
                    {
                        if (sameTriangle(oneTriangle, anotherTriangle))
                        {
                            isShared = true;
                        }
                    }
                }
                if (isShared == false)
                {
                    polyhedra.push_back(oneTriangle);
                }
            }
        }
        for (const auto& oneTetra : housingTetras)
        {
            auto it = std::remove(tetras_.begin(), tetras_.end(), oneTetra);
            tetras_.erase(it, tetras_.end());
        }

        for (const auto& oneTriangle : polyhedra)
        {
            auto newTetra = std::make_shared<tetra>(
                oneTriangle[0], oneTriangle[1], oneTriangle[2], nodeID, nodes_);
            if (newTetra->sphere_.radius_ > 1.0e-10)
                tetras_.push_back(newTetra);
        }
        std::cout << "XXXX             polyhedra.size " << polyhedra.size()
                  << std::endl;
        std::cout << "XXXX             housingTetra.size "
                  << housingTetras.size() << std::endl;
    }

    std::vector<UIN> toDelete(tetras_.size(), 0);
    for (UIN i = 0; i < tetras_.size(); ++i)
    {
        const auto& oneTetra = tetras_[i];

        for (int j = 0; j < 4; ++j)
        {
            if (oneTetra->hasNode(superTetra->nodes_[j]))
            {
                toDelete[i] = 1;
            }
        }
    }

    for (UIN i = 0; i < tetras_.size(); ++i)
    {
        const auto& oneTetra = tetras_[i];
        const auto center = oneTetra->centerQuick(nodes_);
        if (!pointInside({center[0], center[1], center[2]}))
        {
            toDelete[i] = 1;
        }
    }

    const int sumDelete = std::accumulate(toDelete.begin(), toDelete.end(), 0);
    const UIN newSize = tetras_.size() - sumDelete;
    deleteByIndex(toDelete, tetras_, newSize);
    //    delete supernodes
    nodes_.resize(nodes_.size() - 4);
    std::cout << "YYYY Tetra No." << tetras_.size() << std::endl;
}

std::shared_ptr<tetra> solid::getSuperTetra()
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
    const double r = 10.0 * (std::max(std::max(delta[0], delta[1]), delta[2]));

    auto a = std::make_shared<node>(vec3d(r, 0, -0.5 * r));
    auto b = std::make_shared<node>(vec3d(0, r, -0.5 * r));
    auto c = std::make_shared<node>(vec3d(-r, 0, -0.5 * r));
    auto d = std::make_shared<node>(vec3d(0, 0, 0.5 * r));
    nodes_.push_back(a);
    nodes_.push_back(b);
    nodes_.push_back(c);
    nodes_.push_back(d);
    UIN numNodes = nodes_.size();
    auto superTetra = std::make_shared<tetra>(
        numNodes - 4, numNodes - 3, numNodes - 2, numNodes - 1, nodes_);
    return superTetra;
}

void solid::writeVTK(const std::string& fileName) const
{
    vtkNew<vtkPoints> points;
    UIN nodeIDCounter = 0;
    for (auto& oneNode : nodes_)
    {
        points->InsertNextPoint(oneNode->position_.data());
        oneNode->id_ = nodeIDCounter;
        ++nodeIDCounter;
    }
    vtkNew<vtkCellArray> cells;
    for (UIN i = 0; i < tetras_.size(); ++i)
    {
        const auto& oneTetra = tetras_[i];
        vtkNew<vtkTetra> tetra;
        tetra->GetPointIds()->SetId(0, nodes_[oneTetra->nodes_[0]]->id_);
        tetra->GetPointIds()->SetId(1, nodes_[oneTetra->nodes_[1]]->id_);
        tetra->GetPointIds()->SetId(2, nodes_[oneTetra->nodes_[2]]->id_);
        tetra->GetPointIds()->SetId(3, nodes_[oneTetra->nodes_[3]]->id_);
        cells->InsertNextCell(tetra);
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cells);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void solid::writeVerticesToVTK(const std::string& fileName) const
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
