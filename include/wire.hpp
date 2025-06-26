#ifndef H_SHONDY_WIRE
#define H_SHONDY_WIRE

#include "definitions.hpp"
#include "edge.hpp"

namespace shonCloud
{
class wire
{
   public:
    wire(int tag);
    int tag_;

    std::map<int, std::shared_ptr<edge>> edges_;
};

}  // namespace shonCloud
#endif  // H_SHONDY_WIRE
