#pragma once

/// @author Keith Jeffery

#include "Intersection.h"

#include "../math/BBox3.h"
#include "../math/Ray.h"

namespace sp {

class Hitable
{
public:
    virtual ~Hitable() = default;

    [[nodiscard]] bool intersect(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept
    {
        return intersect_impl(ray, limits, isect);
    }

    [[nodiscard]] bool intersect(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept
    {
        return intersect_impl(ray, limits, isect);
    }

    [[nodiscard]] bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept
    {
        return intersect_p_impl(ray, limits);
    }

    [[nodiscard]] BBox3 get_world_bounds() const noexcept
    {
        return get_world_bounds_impl();
    }

private:
    virtual bool  intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept = 0;
    virtual bool  intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept      = 0;
    virtual bool  intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept                   = 0;
    virtual BBox3 get_world_bounds_impl() const noexcept                                                     = 0;
};
} // namespace sp
