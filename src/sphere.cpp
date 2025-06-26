#include "sphere.hpp"
#include "geometricRelations.hpp"
#include "logger.hpp"

namespace shonCloud
{
Sphere::Sphere() : basicGeometry(), radius_(0.0)
{
}

Sphere::Sphere(vec3d centerPoint, double radius, vec3d normal)
    : basicGeometry(centerPoint, normal), radius_(radius)
{
}

}  // namespace shonCloud
