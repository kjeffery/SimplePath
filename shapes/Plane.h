//
// Created by krjef on 10/6/2022.
//

#pragma once

/// @author Keith Jeffery

#include "Shape.h"

namespace sp {
class Plane : public Shape
{
public:
    Plane(AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : Shape(std::move(object_to_world), std::move(world_to_object))
    {
    }

private:
    [[nodiscard]] std::optional<Intersection> intersect_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        // A plane is defined as dot(n, (x - p)) = 0
        // Where
        //     n = (a, b, c)
        //     p = (p, q, r)
        //     x = (x, y, z)
        // This gives the general equation of a plane:
        //    ax + by + cz + d = 0
        // where
        //    d = -ap - bq - cr
        // We have a ray with origin r.o, and direction r.d
        // x can be substituted as (r.o + t*r.d)
        // dot(n, (r.o + t*r.d)) = 0
        // t = -(dot(r.o, n) + d) / dot(r.d, n)

        // We going to intersect in local space:
        //     n = (0, 1, 0)
        //     p = (0, 0, 0)
        // The plane equation is, therefore:
        //     0*x + 1*y + 0*z + d = 0
        //     = y + d = 0
        // And
        //     d = -0*0 - 1*0 - 0*0 = 0
        // dot(n, (r.o + t*r.d)) = 0
        // dot((j, k, l), (0, 1, 0)) = k
        // t = -(dot(r.o, n) + d) / dot(r.d, n)
        // t = -(r.o.y + d) / r.d.y
        // t = -r.o.y / r.d.y

        const Ray local_ray = get_world_to_object()(ray);

        const Vector3& d = local_ray.get_direction();
        if (d.y == 0.0f) {
            // Ray is parallel to plane
            return {};
        }

        const Point3& o = local_ray.get_origin();
        const float   t = -o.y / d.y;
        if (t < limits.m_t_min || t > limits.m_t_max) {
            return {};
        }

        Intersection isect;
        isect.m_normal   = get_object_to_world()(Normal3{ 0.0f, 1.0f, 0.0f });
        isect.m_point    = ray(t);
        isect.m_distance = t;

        return isect;
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        const Ray local_ray = get_world_to_object()(ray);

        const Vector3& d = local_ray.get_direction();
        if (d.y == 0.0f) {
            // Ray is parallel to plane
            return false;
        }

        const Point3& o = local_ray.get_origin();
        const float   t = -o.y / d.y;
        if (t < limits.m_t_min || t > limits.m_t_max) {
            return false;
        }
        return true;
    }

    [[nodiscard]] BBox3 get_object_bounds() const noexcept override
    {
        return {};
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return false;
    }
};
} // namespace sp
