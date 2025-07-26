#ifndef H_MYSIM_VERTEX
#define H_MYSIM_VERTEX

#include "definitions.hpp"
#include <TopoDS_Vertex.hxx>

namespace MySim
{
class vertex
{
   public:
    vertex(const TopoDS_Vertex& occVertex, int tag);
    vertex(const vec3d& position, int type);

    int tag_;
    int type_;
    vec3d position_;
};

}  // namespace MySim
#endif  // H_MYSIM_VERTEX
