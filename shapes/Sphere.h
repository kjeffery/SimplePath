#pragma once

/// @author Keith Jeffery

#include "Shape.h"
#include "ShapeSample.h"

#include "../math/AffineSpace.h"
#include "../math/Vector3.h"

namespace sp {
class Sphere : public Shape
{
public:
    Sphere(AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : Shape(object_to_world, world_to_object)
    {
    }

    [[nodiscard]] ShapeSample sample(const Point2& u) const noexcept
    {
        const Point3  local_sample = sample_to_uniform_sphere(u);
        const Normal3 normal{ local_sample }; // local_sample - our center
        return { get_object_to_world()(local_sample), get_object_to_world()(normal) };
    }

    [[nodiscard]] ShapeSample sample(const Point3& observer_world, const Point2& u) const noexcept
    {
        const Point3  observer    = get_world_to_object()(observer_world);
        const Vector3 to_observer = Vector3{ observer }; // observer - our center, which is 0, 0, 0

        // We're inside the sphere.
        if (sqr_length(to_observer) <= 1.0f) {
            return sample(u);
        }

        const Vector3 sample = sample_to_cosine_hemisphere(u);

        // to_observer is in local space
        const ONB onb = ONB::from_v(to_observer);

        // So, we have to transform from canonical space to local space to world space.
        const Point3  local_sample{ onb.to_world(sample) };
        const Normal3 normal{ local_sample }; // local_sample - our center

        assert(is_normalized(normal)); // We have a radius of 1: this should be normalized implicitly.

        // One must be careful here: the ONB works with vectors, but we need to apply the transformation to a point.
        // E.g. translating a vector is a no-op.
        return { get_object_to_world()(local_sample), get_object_to_world()(normal) };
    }

    [[nodiscard]] float pdf(const Point3& observer_world, const Vector3& wi) const noexcept
    {
        const Point3  observer    = get_world_to_object()(observer_world);
        const Vector3 to_observer = Vector3{ observer }; // observer - our center, which is 0, 0, 0

        const float sqr_distance = sqr_length(to_observer);
        // We're inside the sphere.
        if (sqr_distance <= 1.0f) {
            return uniform_sphere_pdf();
        }

        const float sin_theta_max2 = 1.0f / sqr_distance;
        const float cos_theta_max  = std::sqrt(std::max(0.0f, 1.0f - sin_theta_max2));
        return uniform_cone_pdf(cos_theta_max);
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept override
    {
        const auto local_ray = get_world_to_object()(ray);

        const auto& d = local_ray.get_direction();
        const auto& o = local_ray.get_origin();

        const auto a = dot(d, d);
        const auto b = 2.0f * dot(d, Vector3{ o });
        const auto c = dot(Vector3{ o }, Vector3{ o }) - k_radius * k_radius;

        if (auto discriminant = b * b - 4.0f * a * c; discriminant > 0.0f) {
            discriminant = std::sqrt(discriminant);
            auto t       = (-b - discriminant) / (2.0f * a);

            if (t < limits.m_t_min) {
                t = (-b + discriminant) / (2.0f * a);
            }

            if (t < limits.m_t_min || t > limits.m_t_max) {
                return false;
            }

            const Normal3 n{ madd(t, d, Vector3{ o }) / k_radius }; // (o + t * d) / k_radius };
            isect.m_normal = get_object_to_world()(n);
            isect.m_point  = ray(t);
            limits.m_t_max = t;
            return true;
        }

        return false;
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        const auto local_ray = get_world_to_object()(ray);

        const auto& d = local_ray.get_direction();
        const auto& o = local_ray.get_origin();

        const auto a = dot(d, d);
        const auto b = 2.0f * dot(d, Vector3{ o });
        const auto c = dot(Vector3{ o }, Vector3{ o }) - k_radius * k_radius;

        if (auto discriminant = b * b - 4.0f * a * c; discriminant > 0.0f) {
            discriminant = std::sqrt(discriminant);
            auto t       = (-b - discriminant) / (2.0f * a);

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

    [[nodiscard]] BBox3 get_object_bounds() const noexcept override
    {
        return { Point3{ -k_radius, -k_radius, -k_radius }, Point3{ k_radius, k_radius, k_radius } };
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return true;
    }

    static constexpr auto k_radius = 1.0f;
};
} // namespace sp
