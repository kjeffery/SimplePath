#pragma once

/// @author Keith Jeffery

#include "Vector3.h"

#include "../base/Constants.h"

#include <cassert>
#include <limits>

namespace sp {
constexpr float k_ray_epsilon = 0.001f;

struct RayLimits
{
    float m_t_min = k_ray_epsilon;
    float m_t_max = k_infinite_distance;
};

class Ray
{
public:
    Ray(Point3 p, Vector3 d) noexcept
    : m_origin(p)
    , m_direction(d)
    {
    }

    Point3 operator()(const float d) const noexcept
    {
        assert(is_normalized(m_direction));
        return m_origin + m_direction * d;
    }

    const Point3& get_origin() const noexcept
    {
        return m_origin;
    }

    const Vector3& get_direction() const noexcept
    {
        return m_direction;
    }

private:
    Point3  m_origin;
    Vector3 m_direction;
};

inline std::ostream& operator<<(std::ostream& outs, const Ray& ray)
{
    return outs << ray.get_origin() << ' ' << ray.get_direction();
}

} // namespace sp
