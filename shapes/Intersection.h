
#pragma once

/// @author Keith Jeffery

#include "../materials/Material.h"
#include "../math/RGB.h"
#include "../math/Vector3.h"

namespace sp {
struct LightIntersection
{
    RGB L;
};

struct Intersection
{
    Normal3   m_normal;
    Point3    m_point;
    Material* m_material;
};

} // namespace sp
