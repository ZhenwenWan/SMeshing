#include "tetMesh.hpp"
#include "edge.hpp"
#include "MySimModEnums.hpp"
#include "MySimModTriangle.hpp"
#include "tetra.hpp"
#include "vectorUtil.hpp"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <memory>
#include <vector>

namespace MySim
{
tetMesh::tetMesh()
{
}

void tetMesh::readFromStep(const std::string& fileName)
{
    geometricModel_ = std::make_shared<geometricModel>();
    geometricModel_->readSTEP(fileName);
}

void tetMesh::seedVerticesOnEdges(double meshSize)
{
    for (auto& oneEdge : geometricModel_->edges_)
    {
        oneEdge.second->seedVertices(meshSize, geometricModel_->maxVertexTag());
    }
}

void tetMesh::meshFaces(const double& meshSize)
{
    seedVerticesOnEdges(meshSize);
    for (auto& oneFace : geometricModel_->faces_)
    {
        oneFace.second->mesh(meshSize);
    }
}

}  // namespace MySim
