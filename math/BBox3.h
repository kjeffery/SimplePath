#pragma once

/// @author Keith Jeffery

#include "Vector3.h"

namespace sp {
class BBox3
{
public:
    BBox3() noexcept
    : m_min(Point3::positive_infinity())
    , m_max(Point3::negative_infinity())
    {
    }

    BBox3(const Point3& a, const Point3& b) noexcept
    : m_min(min(a, b))
    , m_max(max(a, b))
    {
    }

    const Point3& get_lower() const noexcept
    {
        return m_min;
    }

    const Point3& get_upper() const noexcept
    {
        return m_max;
    }

    void extend(const Point3& p) noexcept
    {
        m_min = min(p, m_min);
        m_max = max(p, m_max);
    }

private:
    Point3 m_min;
    Point3 m_max;
};
} // namespace sp