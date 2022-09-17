#pragma once

/// @author Keith Jeffery

#include "Vector3.h"

#include <cassert>

namespace sp {
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

    float get_t_min() const noexcept
    {
        return m_t_min;
    }

    float get_t_max() const noexcept
    {
        return m_t_max;
    }

private:
    Point3  m_origin;
    Vector3 m_direction;
    float   m_t_min;
    float   m_t_max;
};
} // namespace sp
