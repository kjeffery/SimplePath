

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

    void push_back(std::shared_ptr<const Hitable> p)
    {
        m_primitives.push_back(p);
    }

    void shrink_to_fit()
    {
        m_primitives.shrink_to_fit();
    }

private:
    [[nodiscard]] std::optional<LightIntersection> intersect_lights_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        std::optional<LightIntersection> result;

        auto internal_limits{ limits };
        for (const auto&   p : m_primitives) {
            if (const auto hit = p->intersect_lights(ray, internal_limits); hit) {
                internal_limits.m_t_max = hit->m_distance;
                result                  = hit;
            }
        }
        return result;
    }

    [[nodiscard]] std::optional<Intersection> intersect_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        std::optional<Intersection> result;

        auto internal_limits{ limits };
        for (const auto&   p : m_primitives) {
            if (const auto hit = p->intersect(ray, internal_limits); hit) {
                internal_limits.m_t_max = hit->m_distance;
                result                  = hit;
            }
        }
        return result;
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return std::ranges::any_of(m_primitives, [&ray, &limits](const auto& p) { return p->intersect_p(ray, limits); });
    }

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        return {};
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return std::ranges::all_of(m_primitives, [](const auto& p) { return p->is_bounded(); });
    }

    // To be extra safe, these should be shared_ptrs (or unique_ptrs assuming this is the only owner), but we will defer
    // sole ownership to the Scene class and avoid the overhead of shared_ptrs.
    std::vector<std::shared_ptr<const Hitable>> m_primitives;
};
} // namespace sp
