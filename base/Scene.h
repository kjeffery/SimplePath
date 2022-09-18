
#pragma once

///@author Keith Jeffery

#include "../shapes/Aggregate.h"

#include <memory>
#include <vector>

namespace sp {
class Scene
{
public:
    bool intersect(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept
    {
        return m_accelerator_lights.intersect(ray, limits, isect);
    }

    bool intersect(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept
    {
        return m_accelerator_geometry.intersect(ray, limits, isect);
    }

    bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept
    {
        return m_accelerator_geometry.intersect_p(ray, limits);
    }

private:
    using HitableContainer = std::vector<std::unique_ptr<Hitable>>;

    // Scene owns the geometry + other hitables.
    HitableContainer m_hitables;

    Aggregate m_accelerator_geometry;
    Aggregate m_accelerator_lights;
};
} // namespace sp
