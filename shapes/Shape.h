#pragma once

/// @author Keith Jeffery

#include <utility>

#include "Hitable.h"
#include "Intersection.h"

#include "../math/AffineSpace.h"
#include "../math/BBox.h"
#include "../math/Ray.h"
#include "../math/Transformation.h"

namespace sp {
class Shape : public Hitable
{
public:
    explicit Shape(AffineTransformation object_to_world) noexcept
    : m_object_to_world(std::move(object_to_world))
    {
    }

protected:
    [[nodiscard]] auto get_object_to_world() const noexcept -> const AffineSpace&
    {
        return m_object_to_world.get_transform();
    }

    [[nodiscard]] auto get_world_to_object() const noexcept -> const AffineSpace&
    {
        return m_object_to_world.get_inverse();
    }

private:
    [[nodiscard]] virtual BBox3 get_object_bounds() const noexcept = 0;

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        return get_object_to_world()(get_object_bounds());
    }

    [[nodiscard]] std::optional<LightIntersection> intersect_lights_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(!"Should not get here");
        return {};
    }

    AffineTransformation m_object_to_world;
};
} // namespace sp
