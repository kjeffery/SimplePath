#pragma once

/// @author Keith Jeffery

#include "Shape.h"

#include "../math/AffineSpace.h"

namespace sp {
class Sphere : public Shape
{
public:
    Sphere(AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : Shape(object_to_world, world_to_object)
    {
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        const Ray local_ray = get_world_to_object()(ray);

        const Vector3& d = local_ray.get_direction();
        const Point3&  o = local_ray.get_origin();

        const float a = dot(d, d);
        const float b = 2.0f * dot(d, Vector3{ o });
        const float c = dot(Vector3{ o }, Vector3{ o }) - k_radius * k_radius;

        if (float discriminant = b * b - 4.0f * a * c; discriminant > 0.0f) {
            discriminant = std::sqrt(discriminant);
            float t      = (-b - discriminant) / (2.0f * a);

            if (t < limits.m_t_min) {
                t = (-b + discriminant) / (2.0f * a);
            }

            if (t < limits.m_t_min || t > limits.m_t_max) {
                return false;
            }

            const Normal3 n{ (o + t * d) / k_radius };
            isect.m_normal = get_object_to_world()(n);
            isect.m_point  = ray(t);
            limits.m_t_max = t;
            return true;
        }

        return false;
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        const Ray local_ray = get_world_to_object()(ray);

        const Vector3& d = local_ray.get_direction();
        const Point3&  o = local_ray.get_origin();

        const float a = dot(d, d);
        const float b = 2.0f * dot(d, Vector3{ o });
        const float c = dot(Vector3{ o }, Vector3{ o }) - k_radius * k_radius;

        if (float discriminant = b * b - 4.0f * a * c; discriminant > 0.0f) {
            discriminant = std::sqrt(discriminant);
            float t      = (-b - discriminant) / (2.0f * a);

            if (t < limits.m_t_min) {
                t = (-b + discriminant) / (2.0f * a);
            }

            if (t < limits.m_t_min || t > limits.m_t_max) {
                return false;
            }

            return true;
        }

        return false;
    }

    BBox3 get_object_bounds() const noexcept override
    {
        return { Point3{ -k_radius, -k_radius, -k_radius }, Point3{ k_radius, k_radius, k_radius } };
    }

    bool is_bounded_impl() const noexcept override
    {
        return true;
    }

    static constexpr float k_radius = 1.0f;
};

} // namespace sp