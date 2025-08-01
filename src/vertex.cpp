#include "vertex.hpp"

#include "MySimModEnums.hpp"
#include <BRep_Tool.hxx>
#include <gp_Pnt.hxx>
namespace MySim
{
vertex::vertex(const TopoDS_Vertex& occVertex, int tag)
    : tag_(tag), type_(MySimMod::vertexTypes::fixed)
{
    const gp_Pnt pnt = BRep_Tool::Pnt(occVertex);
    position_[0] = pnt.X();
    position_[1] = pnt.Y();
    position_[2] = pnt.Z();
}

vertex::vertex(const vec3d& position, int type)
    : tag_(-1), type_(type), position_(position)
{
}

}  // namespace MySim
