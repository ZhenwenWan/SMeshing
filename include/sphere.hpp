#ifndef H_SHONDY_SPHERE
#define H_SHONDY_SPHERE

#include "basicGeometry.hpp"

namespace MySim
{
class Sphere : public basicGeometry
{
   public:
    Sphere();
    Sphere(vec3d centerPoint, double radius, vec3d normal);
    virtual ~Sphere(){};

    double radius_;
    bool inFront(const vec3d& point, const double influenceRadius) const;
    bool behind(const vec3d& point, const double influenceRadius) const;
    double projectedDistance(const vec3d& point) const;

   private:
    double sphericalVolume() const;
};

inline double Sphere::sphericalVolume() const
{
    return (4.0 / 3.0) * PI * std::pow(radius_, 3.0);
}

inline double Sphere::projectedDistance(const vec3d& point) const
{
    const vec3d dist = point - centerPoint_;
    return std::fabs(dist.dot(normal_));
}

inline bool Sphere::inFront(const vec3d& point,
                            const double influenceRadius) const
{
    const double projectDist = projectedDistance(point);
    const double dotResult = normal_.dot(point - centerPoint_);
    const double projectPointToCenter = std::sqrt(
        (point - centerPoint_).squaredNorm() - projectDist * projectDist);
    if (projectDist < influenceRadius && dotResult > 0.0 &&
        projectPointToCenter < radius_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool Sphere::behind(const vec3d& point,
                           const double influenceRadius) const
{
    const double projectDist = projectedDistance(point);
    const double dotResult = normal_.dot(point - centerPoint_);
    const double projectPointToCenter = std::sqrt(
        (point - centerPoint_).squaredNorm() - projectDist * projectDist);
    if (projectDist < influenceRadius && dotResult < 0.0 &&
        projectPointToCenter < radius_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

}  // namespace MySim

#endif  // H_SHONDY_SPHERE
