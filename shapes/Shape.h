#pragma once

/// @author Keith Jeffery

#include "Intersection.h"

#include "../math/BBox3.h"
#include "../math/Ray.h"

namespace sp {

class Shape
{
public:
    virtual ~Shape() = default;

    bool  intersect(const Ray& ray, float& t, Intersection& isect) const noexcept;
    bool  intersect_p(const Ray& ray) const noexcept;
    BBox3 get_world_bounds() const noexcept;

private:
    virtual bool  intersect_impl(const Ray& ray, float& t_hit, Intersection& isect) const noexcept = 0;
    virtual bool  intersect_p_impl(const Ray& ray) const noexcept                              = 0;
    virtual BBox3 get_world_bounds_impl() const noexcept                                       = 0;
};
} // namespace sp
