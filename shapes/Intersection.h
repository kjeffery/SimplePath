
#pragma once

/// @author Keith Jeffery

#include "../materials/Material.h"
#include "../math/RGB.h"
#include "../math/Vector3.h"

namespace sp {
struct LightIntersection
{
    float m_distance{ std::numeric_limits<float>::max() };
    RGB   L;
};

struct Intersection
{
    float     m_distance;
    Normal3   m_normal;
    Point3    m_point;
    Material* m_material{ nullptr };
};
} // namespace sp
