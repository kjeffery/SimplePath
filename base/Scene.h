
#pragma once

///@author Keith Jeffery

#include "not_null.h"
#include "SmartToRawPointerIterator.h"

#include "../shapes/Aggregate.h"
#include "../shapes/ListAccelerator.h"

#include <memory>
#include <vector>

namespace sp {
class Scene
{
public:
    using HitableContainer = std::vector<not_null<std::unique_ptr<const Hitable>>>;

#if 0
    Scene(std::unique_ptr<Aggregate> accelerator_geometry, std::unique_ptr<Aggregate> accelerator_lights)
    : m_accelerator_geometry(std::move(accelerator_geometry))
    , m_accelerator_lights(std::move(accelerator_lights))
    {
    }
#endif

    explicit Scene(HitableContainer hitables, HitableContainer lights)
    : m_hitables(std::move(hitables))
    , m_lights(std::move(lights))
    , m_accelerator_geometry(SmartToRawPointerIterator{ m_hitables.cbegin() },
                             SmartToRawPointerIterator{ m_hitables.cend() })
    , m_accelerator_lights(SmartToRawPointerIterator{ m_lights.cbegin() },
                           SmartToRawPointerIterator{ m_lights.cend() })
    {
    }

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
    // Scene owns the geometry + other hitables.
    HitableContainer m_hitables;
    HitableContainer m_lights;

    ListAccelerator m_accelerator_geometry;
    ListAccelerator m_accelerator_lights;
};
} // namespace sp
