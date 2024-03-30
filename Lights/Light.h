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

struct PartialLightSample
{
    float   m_pdf;
    float   m_max_distance;
    RGB     m_L;
    Vector3 m_wi;
};

// This models a diffuse light. We should probably rename it and/or make a Light base class.
class Light : public Hitable
{
public:
    [[nodiscard]] auto sample(const Point3& observer_world, const Normal3& observer_normal, const Point2& u) const noexcept -> LightSample
    {
        const auto [pdf, max_distance, L, wi] = sample_impl(observer_world, u);

        RayLimits limits;
        limits.m_t_min = get_ray_offset(observer_normal, wi);
        limits.m_t_max = max_distance;

        const Ray              occlusion_ray{ observer_world, wi };
        const VisibilityTester visibility_tester{ limits, occlusion_ray };
        return { .m_L = L, .m_pdf = pdf, .m_tester = visibility_tester };
    }

    [[nodiscard]] auto pdf(const Point3& observer_world, const Vector3& wi) const noexcept -> float
    {
        return pdf_impl(observer_world, wi);
    }

    [[nodiscard]] RGB L(const Point3& p, const Normal3& n, const Vector3& wi) const noexcept
    {
        return L_impl(p, n, wi);
    }

private:
    [[nodiscard]] virtual auto sample_impl(const Point3& observer_world, const Point2& u) const noexcept -> PartialLightSample = 0;
    [[nodiscard]] virtual auto pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept -> float = 0;
    [[nodiscard]] virtual auto L_impl(const Point3& observer_world, const Normal3& n, const Vector3& wi) const noexcept -> RGB = 0;
};

class ObjectLight : public Light
{
public:
    explicit ObjectLight(const RGB radiance) noexcept
    : m_radiance(radiance)
    {
    }

    [[nodiscard]] auto get_radiance() const noexcept -> RGB
    {
        return m_radiance;
    }

private:
    [[nodiscard]] auto sample_impl(const Point3& observer_world, const Point2& u) const noexcept -> PartialLightSample override
    {
        const auto [sampled_point, sampled_normal] = shape_sample(observer_world, u);
        const auto to_sample                       = sampled_point - observer_world;
        const auto wi                              = normalize(to_sample);
        const auto pdf                             = shape_pdf(observer_world, wi);

        const auto distance = length(to_sample) - get_ray_offset(sampled_normal, -wi);
        return { .m_pdf = pdf, .m_max_distance = distance, .m_L = m_radiance, .m_wi = wi };
    }

    [[nodiscard]] auto pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept -> float
    {
        return shape_pdf(observer_world, wi);
    }

    [[nodiscard]] auto L_impl(const Point3& observer_world, const Normal3& n, const Vector3& wi) const noexcept -> RGB
    {
        return (dot(n, wi) > 0.0f) ? m_radiance : RGB::black();
    }

    [[nodiscard]] virtual ShapeSample shape_sample(const Point3& observer_world, const Point2& u) const noexcept = 0;
    [[nodiscard]] virtual float       shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept = 0;

    RGB m_radiance;
};

class InfiniteLight : public Light
{
    [[nodiscard]] auto sample_impl(const Point3& observer_world, const Point2& u) const noexcept -> PartialLightSample override
    {
        const auto partial_light_sample = light_sample(observer_world, u);
        assert(partial_light_sample.m_max_distance == k_infinite_distance);
        return partial_light_sample;
    }

    [[nodiscard]] virtual auto light_sample(const Point3& observer_world, const Point2& u) const noexcept -> PartialLightSample = 0;
};

class EnvironmentLight final : public InfiniteLight
{
public:
    explicit EnvironmentLight(const RGB radiance) noexcept
    : m_radiance(radiance)
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
        isect.L = m_radiance;
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

    [[nodiscard]] auto light_sample(const Point3&, const Point2& u) const noexcept -> PartialLightSample override
    {
        const Point3   p   = sample_to_uniform_sphere(u);
        const auto     wi  = Vector3{ p };
        constexpr auto pdf = uniform_sphere_pdf();
        return { .m_pdf = pdf, .m_max_distance = k_infinite_distance, .m_L = m_radiance, .m_wi = wi };
    }

    [[nodiscard]] auto pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept -> float
    {
        return uniform_sphere_pdf();
    }

    [[nodiscard]] auto L_impl(const Point3& observer_world, const Normal3& n, const Vector3& wi) const noexcept -> RGB
    {
        return (dot(n, wi) > 0.0f) ? m_radiance : RGB::black();
    }

    RGB m_radiance;
};

class SphereLight final : public ObjectLight
{
public:
    explicit SphereLight(RGB radiance, AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : ObjectLight(radiance)
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
