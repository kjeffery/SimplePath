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

private:
    Point3 m_min;
    Point3 m_max;
};
} // namespace sp