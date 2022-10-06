

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
    requires std::input_iterator<Iterator>
    ListAccelerator(Iterator first, Iterator last)
    : m_primitives{ first, last }
    {
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept override
    {
        bool hit = false;
        for (const auto& p : m_primitives) {
            // limits should be automatically updated so that the max value is the closest hit so far.
            hit = p->intersect(ray, limits, isect) || hit;
        }
        return hit;
    }

    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        bool hit = false;
        for (const auto& p : m_primitives) {
            // limits should be automatically updated so that the max value is the closest hit so far.
            hit = p->intersect(ray, limits, isect) || hit;
        }
        return hit;
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        for (const auto& p : m_primitives) {
            if (p->intersect_p(ray, limits)) {
                return true;
            }
        }
        return false;
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return BBox3();
    }

    bool is_bounded_impl() const noexcept override
    {
        return std::ranges::all_of(m_primitives, [](const auto& p) { return p->is_bounded(); });
    }

    // To be extra safe, these should be shared_ptrs (or unique_ptrs assuming this is the only owner), but we will defer
    // sole ownership to the Scene class and avoid the overhead of shared_ptrs.
    std::vector<not_null<std::shared_ptr<const Hitable>>> m_primitives;
};
} // namespace sp
