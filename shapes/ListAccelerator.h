

#pragma once

/// @author Keith Jeffery

#include "Aggregate.h"

#include "../base/not_null.h"

#include <iterator>
#include <vector>

namespace sp {
class ListAccelerator : public Aggregate
{
public:
    template <typename Iterator>
    requires std::forward_iterator<Iterator>
    ListAccelerator(Iterator first, Iterator last)
    : m_primitives{ first, last }
    {
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept override
    {
        return false;
    }

    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        return false;
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return false;
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return BBox3();
    }

    // To be extra safe, these should be shared_ptrs (or unique_ptrs assuming this is the only owner), but we will defer
    // sole ownership to the Scene class and avoid the overhead of shared_ptrs.
    std::vector<not_null<Hitable*>> m_primitives;
};
} // namespace sp
