#pragma once

/// @author Keith Jeffery

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
    explicit Light(RGB radiance) noexcept
    : m_radiance(radiance)
    {
    }

    LightSample sample(const Point3& observer_world, const Normal3& observer_normal, const Point2& u) const noexcept
    {
        const ShapeSample s_sample  = shape_sample(observer_world, u);
        const Vector3     to_sample = s_sample.m_p - observer_world;
        const Vector3     wi        = normalize(to_sample);
        const float       pdf       = shape_pdf(observer_world, wi);

        RayLimits limits;
        limits.m_t_min = get_ray_offset(observer_normal, wi);
        limits.m_t_max = length(to_sample);

        const Ray              occlusion_ray{ observer_world, wi };
        const VisibilityTester visibility_tester{ limits, occlusion_ray };
        const RGB              radiance = L(s_sample.m_p, s_sample.m_n, -wi);
        return { radiance, pdf, visibility_tester };
    }

    float pdf(const Point3& observer_world, const Vector3& wi) const noexcept
    {
        return shape_pdf(observer_world, wi);
    }

    RGB L(const Point3& p, const Normal3& n, const Vector3& w) const noexcept
    {
        return (dot(n, w) > 0.0f) ? m_radiance : RGB::black();
    }

protected:
    const RGB& get_radiance() const noexcept
    {
        return m_radiance;
    }

private:
    virtual ShapeSample shape_sample(const Point3& observer_world, const Point2& u) const noexcept = 0;
    virtual float       shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept  = 0;

    RGB m_radiance;
};

class EnvironmentLight : public Light
{
public:
    explicit EnvironmentLight(RGB radiance) noexcept
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

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return false;
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return sp::BBox3{};
    }

    bool is_bounded_impl() const noexcept override
    {
        return false;
    }

    ShapeSample shape_sample(const Point3&, const Point2& u) const noexcept override
    {
        const Point3  p = sample_to_uniform_sphere(u);
        const Normal3 n = -Normal3{ p };
        return { p, n };
    }

    float shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept override
    {
        return uniform_sphere_pdf();
    }
};

class SphereLight : public Light
{
public:
    explicit SphereLight(RGB radiance, AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : Light(radiance)
    , m_sphere{ object_to_world, world_to_object }
    {
    }

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

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return m_sphere.intersect_p(ray, limits);
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return m_sphere.get_world_bounds();
    }

    bool is_bounded_impl() const noexcept override
    {
        return true;
    }

private:
    ShapeSample shape_sample(const Point3& observer_world, const Point2& u) const noexcept override
    {
        return m_sphere.sample(observer_world, u);
    }

    float shape_pdf(const Point3& observer_world, const Vector3& wi) const noexcept override
    {
        return m_sphere.pdf(observer_world, wi);
    }

private:
    Sphere m_sphere;
};

} // namespace sp
