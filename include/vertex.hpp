#ifndef H_SHONDY_VERTEX
#define H_SHONDY_VERTEX

#include "definitions.hpp"
#include <TopoDS_Vertex.hxx>

namespace shonCloud
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

}  // namespace shonCloud
#endif  // H_SHONDY_VERTEX
