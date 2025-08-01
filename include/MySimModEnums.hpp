#ifndef H_MYSIM_SHONMODELERENUMS
#define H_MYSIM_SHONMODELERENUMS

#include "definitions.hpp"

namespace MySimMod
{
enum vertexTypes : int
{
    fixed = 0,
    onEdge = 1,
    free = 2,
};

enum edgeGeometryTypes : int
{
    line = 0,
    circle = 1,
};

enum faceGeometryTypes : int
{
    plane = 0,
    sphere = 1,
};

}  // namespace MySimMod
#endif  // H_MYSIM_SHONMODELERENUMS
