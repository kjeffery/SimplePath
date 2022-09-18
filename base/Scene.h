
#pragma once

///@author Keith Jeffery

#include "not_null.h"

#include "../shapes/Aggregate.h"
#include "../shapes/ListAccelerator.h"

#include <memory>
#include <vector>

namespace sp {
class Scene
{
public:
    Scene(std::unique_ptr<Aggregate> accelerator_geometry, std::unique_ptr<Aggregate> accelerator_lights)
    : m_accelerator_geometry(std::move(accelerator_geometry))
    , m_accelerator_lights(std::move(accelerator_lights))
    {
    }

    bool intersect(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept
    {
        return m_accelerator_lights->intersect(ray, limits, isect);
    }

    bool intersect(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept
    {
        return m_accelerator_geometry->intersect(ray, limits, isect);
    }

    bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept
    {
        return m_accelerator_geometry->intersect_p(ray, limits);
    }

private:
    using HitableContainer = std::vector<not_null<std::unique_ptr<Hitable>>>;

    // Scene owns the geometry + other hitables.
    HitableContainer m_hitables;

    not_null<std::unique_ptr<Aggregate>> m_accelerator_geometry;
    not_null<std::unique_ptr<Aggregate>> m_accelerator_lights;
};
} // namespace sp
