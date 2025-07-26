#ifndef H_SHONDY_WIRE
#define H_SHONDY_WIRE

#include "definitions.hpp"
#include "edge.hpp"

namespace MySim
{
class wire
{
   public:
    wire(int tag);
    int tag_;

    std::map<int, std::shared_ptr<edge>> edges_;
};

}  // namespace MySim
#endif  // H_SHONDY_WIRE
