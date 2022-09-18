//
// Created by krjef on 9/18/2022.
//

#pragma once

/// @author Keith Jeffery

#include "Hitable.h"
#include "Shape.h"

#include "../base/not_null.h"

namespace sp {

class GeometricPrimitive : public Hitable
{
public:
private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        if (m_shape->intersect(ray, limits, isect)) {
            isect.m_material = m_material;
            return true;
        }
        return false;
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return m_shape->intersect_p(ray, limits);
    }

    not_null<Shape*>    m_shape;
    not_null<Material*> m_material;
};

} // namespace sp
