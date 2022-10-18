#pragma once

/// @author Keith Jeffery

#include "../base/Scene.h"
#include "../math/HSV.h"
#include "../math/ONB.h"
#include "../math/Ray.h"
#include "../math/RGB.h"
#include "../math/Sampler.h"
#include "../math/Sampling.h"

#include <iostream> // TODO: temp

namespace sp {

class Integrator
{
public:
    virtual ~Integrator() = default;

    [[nodiscard]] RGB integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const
    {
        return integrate_impl(ray, scene, sampler);
    }

private:
    virtual RGB integrate_impl(const Ray& ray, const Scene& scene, Sampler& sampler) const = 0;
};

// A simple test integrator that ignores the ray direction and simply uses the x and y values.
class MandelbrotIntegrator : public Integrator
{
public:
    MandelbrotIntegrator(const int image_width, const int image_height) noexcept
    : m_image_width(image_width)
    , m_image_height(image_height)
    {
    }

private:
    RGB integrate_impl(const Ray& ray, const Scene&, Sampler&) const override
    {
        constexpr float x0 = -2.0f;
        constexpr float x1 = +1.0f;
        constexpr float y0 = -1.0f;
        constexpr float y1 = +1.0f;

        const float dx = (x1 - x0) / m_image_width;
        const float dy = (y1 - y0) / m_image_height;

        // TODO: intersect with plane in front of camera (transform to local space?)
        const float& px = ray.get_origin().x;
        const float& py = ray.get_origin().y;

        const float x = x0 + px * dx;
        const float y = y0 + py * dy;

        const float value = static_cast<float>(mandel(x, y)) / s_max_iterations;

        const float hue        = std::fmod(std::pow(value * 360.0f, 1.5f), 360.0f) / 360.0f;
        const float saturation = 1.0f;
        const HSV   hsv{ hue, saturation, value };
        return to_rgb(hsv);
    }

    static int mandel(const float c_re, const float c_im) noexcept
    {
        float z_re = c_re;
        float z_im = c_im;

        int i = 0;
        for (; i < s_max_iterations; ++i) {
            if (z_re * z_re + z_im * z_im > 4.0f) {
                return i;
            }

            const float new_re = z_re * z_re - z_im * z_im;
            const float new_im = 2.0f * z_re * z_im;
            z_re               = c_re + new_re;
            z_im               = c_im + new_im;
        }
        return i;
    }

    static constexpr int s_max_iterations = 4096;

    int m_image_width;
    int m_image_height;
};

class BruteForceIntegrator : public Integrator
{
private:
    RGB integrate_impl(const Ray& ray, const Scene& scene, Sampler& sampler) const override
    {
        return do_integrate(ray, scene, sampler, 0);
    }

    RGB do_integrate(const Ray& ray, const Scene& scene, Sampler& sampler, int depth) const
    {
        if (depth >= scene.max_depth) {
            return RGB::black();
        }

        RayLimits         limits;
        LightIntersection light_intersection;
        const bool        hit_light = scene.intersect(ray, limits, light_intersection);
        if (Intersection geometry_intersection; scene.intersect(ray, limits, geometry_intersection)) {
            const Vector3  wo = -ray.get_direction();
            const Normal3& n  = geometry_intersection.m_normal;

            const auto shading_result = geometry_intersection.m_material->sample(wo, n, sampler);
            if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                return RGB::black();
            }

            // Get next direction
            const auto& wi                = shading_result.direction;
            const auto  cosine            = dot(wi, n);
            const auto  outgoing_position = ray(limits.m_t_max);
            const Ray   outgoing_ray{ outgoing_position, wi };

            return do_integrate(outgoing_ray, scene, sampler, depth + 1) * cosine * shading_result.color /
                   shading_result.pdf;
        } else if (hit_light) {
            return light_intersection.L;
        } else {
            return RGB::black();
        }
    }
};

class BruteForceIntegratorIterative : public Integrator
{
private:
    RGB integrate_impl(const Ray& ray, const Scene& scene, Sampler& sampler) const override
    {
        return do_integrate(ray, scene, sampler);
    }

    RGB do_integrate(Ray ray, const Scene& scene, Sampler& sampler) const
    {
        RGB throughput = RGB::white();
        RGB L = RGB::black();

        for (int depth = 0; depth < scene.max_depth; ++depth) {
            RayLimits         limits;
            LightIntersection light_intersection;
            const bool        hit_light = scene.intersect(ray, limits, light_intersection);
            if (Intersection geometry_intersection; scene.intersect(ray, limits, geometry_intersection)) {
                const Vector3  wo = -ray.get_direction();
                const Normal3& n  = geometry_intersection.m_normal;

                const auto shading_result = geometry_intersection.m_material->sample(wo, n, sampler);
                if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                    break;
                }

                // Get next direction
                const auto& wi                = shading_result.direction;
                const auto  cosine            = dot(wi, n);
                const auto  outgoing_position = ray(limits.m_t_max);
                const Ray   outgoing_ray{ outgoing_position, wi };
                ray = outgoing_ray;

                const RGB contribution = cosine * shading_result.color / shading_result.pdf;
                throughput *= contribution;
            } else if (hit_light) {
                L += throughput * light_intersection.L;
                break;
            } else {
                break;
            }
        }
        return L;
    }
};

class BruteForceIntegratorIterativeRR : public Integrator
{
private:
    RGB integrate_impl(const Ray& ray, const Scene& scene, Sampler& sampler) const override
    {
        return do_integrate(ray, scene, sampler);
    }

    RGB do_integrate(Ray ray, const Scene& scene, Sampler& sampler) const
    {
        RGB throughput = RGB::white();
        RGB L = RGB::black();

        for (int depth = 0; depth < scene.max_depth; ++depth) {
            RayLimits         limits;
            LightIntersection light_intersection;
            const bool        hit_light = scene.intersect(ray, limits, light_intersection);
            if (Intersection geometry_intersection; scene.intersect(ray, limits, geometry_intersection)) {
                const Vector3  wo = -ray.get_direction();
                const Normal3& n  = geometry_intersection.m_normal;

                const auto shading_result = geometry_intersection.m_material->sample(wo, n, sampler);
                if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                    break;
                }

                // Get next direction
                const auto& wi                = shading_result.direction;
                const auto  cosine            = dot(wi, n);
                const auto  outgoing_position = ray(limits.m_t_max);
                const Ray   outgoing_ray{ outgoing_position, wi };
                ray = outgoing_ray;

                const RGB contribution = cosine * shading_result.color / shading_result.pdf;
                throughput *= contribution;

                // TODO: dynamic
                if (depth >= scene.min_depth) {
                    constexpr float mean_throughput_luminance = 0.1f;
                    const float lum = relative_luminance(throughput);
                    if (lum < mean_throughput_luminance) {
                        // q is the probability of continuing
                        const float q = std::max(0.05f, lum/mean_throughput_luminance);
                        assert(q < 1.0f);
                        if (sampler.get_next_1D() < q) {
                            throughput /= q;
                        } else {
                            break;
                        }
                    }
                }

            } else if (hit_light) {
                L += throughput * light_intersection.L;
                break;
            } else {
                break;
            }
        }
        return L;
    }
};

} // namespace sp
