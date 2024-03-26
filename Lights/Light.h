#pragma once

/// @author Keith Jeffery

#include <utility>

#include "../shapes/Hitable.h"
#include "../shapes/Sphere.h"

namespace sp {
struct VisibilityTester
{
    RayLimits m_limits;
    Ray       m_ray;
};

struct LightSample
{
    RGB              m_L;
    float            m_pdf;
    VisibilityTester m_tester;
};

// This models a diffuse light. We should probably rename it and/or make a Light base class.
class Light : public Hitable
{
public:
    explicit Light(const RGB radiance) noexcept
    : m_radiance(radiance)
    {
    }

    [[nodiscard]] auto sample(const Point3& observer_world, const Normal3& observer_normal, const Point2& u) const noexcept -> LightSample
    {
        const auto s_sample  = shape_sample(observer_world, u);
        const auto to_sample = s_sample.m_p - observer_world;
        const auto wi        = normalize(to_sample);
        const auto pdf       = shape_pdf(observer_world, wi);

        RayLimits limits;
        limits.m_t_min = get_ray_offset(observer_normal, wi);
        limits.m_t_max = length(to_sample) - k_ray_epsilon;

        const Ray              occlusion_ray{ observer_world, wi };
        const VisibilityTester visibility_tester{ limits, occlusion_ray };
        const RGB              radiance = L(s_sample.m_p, s_sample.m_n, -wi);
        return { radiance, pdf, visibility_tester };
    }

    [[nodiscard]] auto pdf(const Point3& observer_world, const Vector3& wi) const noexcept -> float
    {
        return shape_pdf(observer_world, wi);
    }

    [[nodiscard]] RGB L(const Point3& p, const Normal3& n, const Vector3& w) const noexcept
    {
        return (dot(n, w) > 0.0f) ? m_radiance : RGB::black();
    }

protected:
    [[nodiscard]] const RGB& get_radiance() const noexcept
    {
        return m_radiance;
    }

private:
    [[nodiscard]] virtual ShapeSample shape_sample(const Point3& observer_world, const Point2& u) const noexcept = 0;
    [[nodiscard]] virtual float       shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept = 0;

    RGB m_radiance;
};

class EnvironmentLight : public Light
{
public:
    explicit EnvironmentLight(const RGB radiance) noexcept
    : Light(radiance)
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
        isect.L = get_radiance();
        return true;
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return false;
    }

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        return sp::BBox3{};
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return false;
    }

    [[nodiscard]] ShapeSample shape_sample(const Point3&, const Point2& u) const noexcept override
    {
        const Point3  p = sample_to_uniform_sphere(u);
        const Normal3 n = -Normal3{ p };
        return { p, n };
    }

    [[nodiscard]] float shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept override
    {
        return uniform_sphere_pdf();
    }
};

class SphereLight : public Light
{
public:
    explicit SphereLight(RGB radiance, AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : Light(radiance)
    , m_sphere{ std::move(object_to_world), std::move(world_to_object) }
    {
    }

private:
    bool intersect_impl(const Ray&, RayLimits&, Intersection&) const noexcept override
    {
        assert(!"Should not get here");
        return false;
    }

    bool intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept override
    {
        // TODO: may need normal
        if (Intersection geom_isect; m_sphere.intersect(ray, limits, geom_isect)) {
            isect.L = get_radiance();
            return true;
        }
        return false;
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return m_sphere.intersect_p(ray, limits);
    }

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        return m_sphere.get_world_bounds();
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return true;
    }

    [[nodiscard]] ShapeSample shape_sample(const Point3& observer_world, const Point2& u) const noexcept override
    {
        return m_sphere.sample(observer_world, u);
    }

    [[nodiscard]] float shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept override
    {
        return m_sphere.pdf(observer_world, wi);
    }

private:
    Sphere m_sphere;
};
} // namespace sp
