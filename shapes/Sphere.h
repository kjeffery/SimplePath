#pragma once

/// @author Keith Jeffery

#include "Shape.h"

#include "../math/AffineSpace.h"

namespace sp {
class Sphere : public Shape
{
public:
    Sphere(AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : m_object_to_world(object_to_world)
    , m_world_to_object(world_to_object)
    {
    }

private:
    bool intersect_impl(const Ray& world_ray, float& t_hit, Intersection& isect) const noexcept override
    {
        const Ray r = m_world_to_object(r);

        const Vector3& d = r.get_direction();
        const Point3&  o = r.get_origin();

        const float a = dot(d, d);
        const float b = 2.0f * dot(d, Vector3{ o });
        const float c = dot(Vector3{ o }, Vector3{ o }) - k_radius * k_radius;

        if (float discriminant = b * b - 4.0f * a * c; discriminant > 0.0f) {
            discriminant = std::sqrt(discriminant);
            float t      = (-b - discriminant) / (2.0f * a);

            if (t < r.get_t_min()) {
                t = (-b + discriminant) / (2.0f * a);
            }

            if (t < r.get_t_min() || t > r.get_t_max()) {
                return false;
            }

            const Normal3 n{ (o + t * d) / k_radius };
            isect.normal = m_object_to_world(n);
            isect.point  = world_ray(t);
            t_hit        = t;
            return true;
        }

        return false;
    }

    bool intersect_p_impl(const Ray& ray) const noexcept override
    {
        const Ray r = m_world_to_object(r);

        const Vector3& d = r.get_direction();
        const Point3&  o = r.get_origin();

        const float a = dot(d, d);
        const float b = 2.0f * dot(d, Vector3{ o });
        const float c = dot(Vector3{ o }, Vector3{ o }) - k_radius * k_radius;

        if (float discriminant = b * b - 4.0f * a * c; discriminant > 0.0f) {
            discriminant = std::sqrt(discriminant);
            float t      = (-b - discriminant) / (2.0f * a);

            if (t < r.get_t_min()) {
                t = (-b + discriminant) / (2.0f * a);
            }

            if (t < r.get_t_min() || t > r.get_t_max()) {
                return false;
            }

            return true;
        }

        return false;
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return m_object_to_world(get_object_bounds());
    }

    BBox3 get_object_bounds() const noexcept
    {
        return { Point3{ -k_radius, -k_radius, -k_radius }, Point3{ k_radius, k_radius, k_radius } };
    }

    static constexpr float k_radius = 1.0f;

    AffineSpace m_object_to_world;
    AffineSpace m_world_to_object;
};

} // namespace sp