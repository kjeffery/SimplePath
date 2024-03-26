#pragma once

/// @author Keith Jeffery

#include <utility>

#include "Hitable.h"
#include "Intersection.h"

#include "../math/AffineSpace.h"
#include "../math/BBox.h"
#include "../math/Ray.h"

namespace sp {
class Shape : public Hitable
{
public:
    Shape(AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : m_object_to_world(std::move(object_to_world))
    , m_world_to_object(std::move(world_to_object))
    {
    }

protected:
    [[nodiscard]] const AffineSpace& get_object_to_world() const noexcept
    {
        return m_object_to_world;
    }

    [[nodiscard]] const AffineSpace& get_world_to_object() const noexcept
    {
        return m_world_to_object;
    }

private:
    [[nodiscard]] virtual BBox3 get_object_bounds() const noexcept = 0;

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        return get_object_to_world()(get_object_bounds());
    }

    bool intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept override
    {
        assert(!"Should not get here");
        return false;
    }

    AffineSpace m_object_to_world;
    AffineSpace m_world_to_object;
};
} // namespace sp
