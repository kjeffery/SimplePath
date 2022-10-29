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

Point3 get_ray_offset(const Normal3& n, const float cos_d, const Point3& p) noexcept
{
    assert(cos_d != 0.0f);

    const float h = k_ray_epsilon / cos_d;
    const Vector3 offset = h * Vector3{n};
    const Point3 po = p + offset;
    return po;
}

Point3 get_ray_offset(const Normal3& n, const Vector3& d, const Point3& p) noexcept
{
    // To avoid self intersections (due to numeric precision issues), we offset our ray origins by a small amount.
    // However, when the ray direction is shallow (compared to the surface), this offset has to be larger.

    // Let t be the angle between the normal and the ray direction. We want to offset our canonical distance (a) by
    // dividing the cosine term.

    // cos(t) = a/h
    // h = a/cos(t)
    // +--------+
    // |       /
    // |      /
    // |     /
    // |    /
    // |   /
    // |  /
    // |t/
    // |/
    // +

    // There are certainly more advanced ways to handle this offsetting business (e.g., see PBRT v3+, where they track
    // accumulated floating-point error). That's Pharr from what we're doing here.

    return get_ray_offset(n, dot(n, d), p);

    //const float h = k_ray_epsilon / dot(n, d);
    //const Vector3 offset = h * Vector3{n};
    //const Point3 po = p + offset;
    //return po;
}

inline std::ostream& operator<<(std::ostream& outs, const Ray& ray)
{
    return outs << ray.get_origin() << ' ' << ray.get_direction();
}

} // namespace sp
