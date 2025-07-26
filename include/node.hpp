#ifndef H_SHONDY_NODE
#define H_SHONDY_NODE

#include "definitions.hpp"
#include <memory>

namespace MySim
{
class tetra;
class node
{
   public:
    node(const vec3d& point);

    vec3d position_;
    std::shared_ptr<tetra> owner_;
    UIN id_;

    const vec3d& center() const;
};

inline const vec3d& node::center() const
{
    return position_;
}

}  // namespace MySim
#endif  // H_SHONDY_NODE
