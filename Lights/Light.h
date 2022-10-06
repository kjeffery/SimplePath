#pragma once

/// @author Keith Jeffery

#include "../shapes/Hitable.h"

namespace sp {

class Light : public Hitable
{
public:
    virtual ~Light() = default;

private:

};

class EnvironmentLight : public Light
{
public:
    explicit EnvironmentLight(RGB color) noexcept
        : m_color(color)
    {
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        assert(!"Should not get here");
        return false;
    }

    bool intersect_impl(const Ray&, RayLimits& limits, LightIntersection& isect) const noexcept override
    {
        if (limits.m_t_max < k_infinite_distance) {
            return false;
        }
        isect.L = m_color;
        return true;
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return false;
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return sp::BBox3{};
    }

    RGB m_color;
};

} // namespace sp
