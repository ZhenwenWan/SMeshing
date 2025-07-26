#include "node.hpp"

namespace MySim
{
node::node(const vec3d& point) : position_(point), owner_(), id_(0)
{
}

}  // namespace MySim
