//
// Created by krjef on 9/18/2022.
//

#pragma once

/// @author Keith Jeffery

#include "Hitable.h"
#include "Shape.h"

#include "../base/not_null.h"
#include "../materials/Material.h"

namespace sp {

class GeometricPrimitive : public Hitable
{
public:
    GeometricPrimitive(std::shared_ptr<Shape> shape, std::shared_ptr<Material> material) noexcept
    : m_shape(shape)
    , m_material(material)
    {
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept override
    {
        assert(!"Should not get here");
        return false;
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return m_shape->get_world_bounds();
    }

    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        if (m_shape->intersect(ray, limits, isect)) {
            isect.m_material = m_material.get().get();
            return true;
        }
        return false;
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return m_shape->intersect_p(ray, limits);
    }

    bool is_bounded_impl() const noexcept override
    {
        return m_shape->is_bounded();
    }

    not_null<std::shared_ptr<Shape>>    m_shape;
    not_null<std::shared_ptr<Material>> m_material;
};

} // namespace sp
