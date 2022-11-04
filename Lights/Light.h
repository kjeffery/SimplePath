#pragma once

/// @author Keith Jeffery

#include "../shapes/Hitable.h"
#include "../shapes/Sphere.h"

namespace sp {

class Light : public Hitable
{
public:
    Point3 sample(const Point3& observer_world, const Point2& u) const noexcept;
    float  pdf(const Point3& observer_world, const Vector3& wi) const noexcept;

private:
    virtual Point3 sample_impl(const Point3& observer_world, const Point2& u) const noexcept = 0;
    virtual float  pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept  = 0;
};

class EnvironmentLight : public Light
{
public:
    explicit EnvironmentLight(RGB radiance) noexcept
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

    Point3 sample_impl(const Point3&, const Point2& u) const noexcept override
    {
        return sample_to_uniform_sphere(u);
    }

    float pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept override
    {
        return uniform_sphere_pdf();
    }

    RGB m_radiance;
};

class SphereLight : public Light
{
public:
    explicit SphereLight(RGB radiance, AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : m_sphere{ object_to_world, world_to_object }
    , m_radiance{ radiance }
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
            isect.L = m_radiance;
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
    Point3 sample_impl(const Point3& observer_world, const Point2& u) const noexcept override
    {
        return m_sphere.sample(observer_world, u);
    }

    float pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept override
    {
        return m_sphere.pdf(observer_world, wi);
    }

private:
    Sphere m_sphere;
    RGB    m_radiance;
};

} // namespace sp
