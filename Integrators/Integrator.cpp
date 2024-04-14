/// @author Keith Jeffery

#include "Integrator.h"

#include "../base/RunningStats.h"
#include "../base/Scene.h"
#include "../math/HSV.h"
#include "../math/ONB.h"
#include "../math/Ray.h"
#include "../math/RGB.h"
#include "../math/Sampler.h"
#include "../math/Sampling.h"
#include "../Lights/Light.h"

#include <iostream> // TODO: temp

namespace sp {
namespace {
bool occluded(const VisibilityTester& tester, const Scene& scene)
{
    return scene.intersect_p(tester.m_ray, tester.m_limits);
}
} // anonymous namespace

auto string_to_integrator_type(std::string_view s) -> IntegratorType
{
    s = trim(s);
    if (s == "mandelbrot") {
        return IntegratorType::Mandelbrot;
    }
    if (s == "brute_force") {
        return IntegratorType::BruteForce;
    }
    if (s == "brute_force_iterative") {
        return IntegratorType::BruteForceIterative;
    }
    if (s == "brute_force_iterative_rr") {
        return IntegratorType::BruteForceIterativeRR;
    }
    if (s == "iterative_rrnee") {
        return IntegratorType::IterativeRRNEE;
    }
    if (s == "direct_lighting") {
        return IntegratorType::DirectLighting;
    }
    if (s == "whitted") {
        return IntegratorType::Whitted;
    }

    throw std::runtime_error("Unknown integrator type");
}

MandelbrotIntegrator::MandelbrotIntegrator(const int image_width, const int image_height) noexcept
: m_image_width(image_width)
, m_image_height(image_height)
{
}

RGB MandelbrotIntegrator::integrate_impl(const Ray& ray,
                                         const Scene&,
                                         MemoryArena& arena,
                                         Sampler&,
                                         const Point2& pixel_coords) const
{
    constexpr float x0 = -2.0f;
    constexpr float x1 = +1.0f;
    constexpr float y0 = -1.0f;
    constexpr float y1 = +1.0f;

    const float dx = (x1 - x0) / m_image_width;
    const float dy = (y1 - y0) / m_image_height;

    // TODO: intersect with plane in front of camera (transform to local space?)
    const float& px = pixel_coords.x;
    const float& py = pixel_coords.y;

    const float x = x0 + px * dx;
    const float y = y0 + py * dy;

    const float value = static_cast<float>(mandel(x, y)) / s_max_iterations;

    const float hue        = std::fmod(std::pow(value * 360.0f, 1.5f), 360.0f) / 360.0f;
    const float saturation = 1.0f;
    const HSV   hsv{ hue, saturation, value };
    return to_rgb(hsv);
}

int MandelbrotIntegrator::mandel(const float c_re, const float c_im) noexcept
{
    float z_re = c_re;
    float z_im = c_im;

    int i = 0;
    for (; i < MandelbrotIntegrator::s_max_iterations; ++i) {
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
RGB BruteForceIntegrator::integrate_impl(const Ray&   ray,
                                         const Scene& scene,
                                         MemoryArena& arena,
                                         Sampler&     sampler,
                                         const Point2&) const
{
    return do_integrate(ray, scene, arena, sampler, 0);
}

// TODO: this needs limits, if just to avoid self intersections
RGB BruteForceIntegrator::do_integrate(const Ray& ray, const Scene& scene, MemoryArena& arena, Sampler& sampler, int depth) const
{
    if (depth >= scene.max_depth) {
        return RGB::black();
    }

    RayLimits  limits;
    const auto light_intersection = scene.intersect_lights(ray, limits);
    if (light_intersection) {
        limits.m_t_max = light_intersection->m_distance;
    }
    if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
        const Vector3  wo = -ray.get_direction();
        const Normal3& n  = geometry_intersection->m_normal;

        const auto shading_result = geometry_intersection->m_material->sample(arena, wo, n, sampler);
        if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
            return RGB::black();
        }

        // Get next direction
        const auto& wi                = shading_result.direction;
        const auto  cosine            = dot(wi, n);
        const auto  outgoing_position = ray(geometry_intersection->m_distance);
        const Ray   outgoing_ray{ outgoing_position, wi };

        return do_integrate(outgoing_ray, scene, arena, sampler, depth + 1) * cosine * shading_result.color /
                shading_result.pdf;
    } else if (light_intersection) {
        return light_intersection->L;
    } else {
        return RGB::black();
    }
}

RGB BruteForceIntegratorIterative::integrate_impl(const Ray&   ray,
                                                  const Scene& scene,
                                                  MemoryArena& arena,
                                                  Sampler&     sampler,
                                                  const Point2&) const
{
    return do_integrate(ray, scene, arena, sampler);
}

RGB BruteForceIntegratorIterative::do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler) const
{
    RGB       throughput = RGB::white();
    RGB       L          = RGB::black();
    RayLimits limits;

    for (int depth = 0; depth < scene.max_depth; ++depth) {
        const auto light_intersection = scene.intersect_lights(ray, limits);
        if (light_intersection) {
            limits.m_t_max = light_intersection->m_distance;
        }

        if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
            const Vector3  wo = -ray.get_direction();
            const Normal3& n  = geometry_intersection->m_normal;

            const auto shading_result = geometry_intersection->m_material->sample(arena, wo, n, sampler);
            if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                break;
            }

            // Get next direction
            const auto& wi           = shading_result.direction;
            const auto  cosine       = std::abs(dot(wi, n));
            const RGB   contribution = cosine * shading_result.color / shading_result.pdf;
            throughput *= contribution;

            const auto outgoing_position = ray(geometry_intersection->m_distance);
            const Ray  outgoing_ray{ outgoing_position, wi };
            ray            = outgoing_ray;
            limits.m_t_min = get_ray_offset(cosine);
            limits.m_t_max = k_infinite_distance;
        } else if (light_intersection) {
            L += throughput * light_intersection->L;
            break;
        } else {
            break;
        }
    }
    return L;
}

RGB BruteForceIntegratorIterativeRR::integrate_impl(const Ray&   ray,
                                                    const Scene& scene,
                                                    MemoryArena& arena,
                                                    Sampler&     sampler,
                                                    const Point2&) const
{
    return do_integrate(ray, scene, arena, sampler);
}

RGB BruteForceIntegratorIterativeRR::do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler) const
{
    RGB       throughput = RGB::white();
    RGB       L          = RGB::black();
    RayLimits limits;

    constexpr float rr_throughput_luminance_cutoff = 0.1f;
    for (int depth = 0; depth < scene.max_depth; ++depth) {
        const auto light_intersection = scene.intersect_lights(ray, limits);
        if (light_intersection) {
            limits.m_t_max = light_intersection->m_distance;
        }

        if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
            const Vector3  wo = -ray.get_direction();
            const Normal3& n  = geometry_intersection->m_normal;

            const auto shading_result = geometry_intersection->m_material->sample(arena, wo, n, sampler);
            if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                break;
            }

            // Get next direction
            const auto& wi           = shading_result.direction;
            const auto  cosine       = std::abs(dot(wi, n));
            const RGB   contribution = cosine * shading_result.color / shading_result.pdf;
            throughput *= contribution;

            if (depth >= scene.russian_roulette_depth) {
                const float lum = relative_luminance(throughput);
                if (lum < rr_throughput_luminance_cutoff) {
                    // q is the probability of continuing
                    const float q = std::max(0.05f, lum / rr_throughput_luminance_cutoff);
                    assert(q < 1.0f);
                    if (sampler.get_next_1D() < q) {
                        throughput /= q;
                    } else {
                        break;
                    }
                }
            }

            const auto outgoing_position = ray(geometry_intersection->m_distance);
            const Ray  outgoing_ray{ outgoing_position, wi };
            ray            = outgoing_ray;
            limits.m_t_min = get_ray_offset(cosine);
            limits.m_t_max = k_infinite_distance;
        } else if (light_intersection) {
            L += throughput * light_intersection->L;
            break;
        } else {
            break;
        }
    }
    return L;
}

RGB DirectLightingIntegrator::integrate_impl(const Ray&   ray,
                                             const Scene& scene,
                                             MemoryArena& arena,
                                             Sampler&     sampler,
                                             const Point2&) const
{
    return do_integrate(ray, scene, arena, sampler, 0);
}

RGB DirectLightingIntegrator::do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler, int depth) const
{
    RGB L = RGB::black();
    if (depth >= scene.max_depth) {
        return L;
    }

    RGB       throughput = RGB::white();
    RayLimits limits;

    const auto light_intersection = scene.intersect_lights(ray, limits);
    if (light_intersection) {
        limits.m_t_max = light_intersection->m_distance;
    }

    if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
        const auto  wo = -ray.get_direction();
        const auto& n  = geometry_intersection->m_normal;
        scene.for_each_light([&wo, &n, &scene, &arena, &L, &geometry_intersection, &sampler](const Light& light) {
            const auto light_sample = light.sample(geometry_intersection->m_point, geometry_intersection->m_normal, sampler.get_next_2D());
            if (light_sample.m_pdf == 0.0f || light_sample.m_L == RGB::black()) {
                return;
            }

            const auto& wi = light_sample.m_tester.m_ray.get_direction();
            const auto  f  = geometry_intersection->m_material->eval(arena, wo, wi, n, sampler);

            if (f != RGB::black() && !occluded(light_sample.m_tester, scene)) {
                L += f * light_sample.m_L * std::abs(dot(wi, n)) / light_sample.m_pdf;
            }
        });
    } else if (light_intersection) {
        L += throughput * light_intersection->L;
    }
    return L;
}

RGB WhittedIntegrator::integrate_impl(const Ray&   ray,
                                      const Scene& scene,
                                      MemoryArena& arena,
                                      Sampler&     sampler,
                                      const Point2&) const
{
    return do_integrate(ray, scene, arena, sampler, 0);
}

RGB WhittedIntegrator::do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler, int depth) const
{
    RGB L = RGB::black();
    if (depth >= scene.max_depth) {
        return L;
    }

    RGB       throughput = RGB::white();
    RayLimits limits;

    const auto light_intersection = scene.intersect_lights(ray, limits);
    if (light_intersection) {
        limits.m_t_max = light_intersection->m_distance;
    }

    if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
        const auto  wo = -ray.get_direction();
        const auto& n  = geometry_intersection->m_normal;

        scene.for_each_light([&wo, &n, &scene, &arena, &L, &geometry_intersection, &sampler](const Light& light) {
            const auto light_sample = light.sample(geometry_intersection->m_point, geometry_intersection->m_normal, sampler.get_next_2D());
            if (light_sample.m_pdf == 0.0f || light_sample.m_L == RGB::black()) {
                return;
            }

            const auto& wi = light_sample.m_tester.m_ray.get_direction();
            const auto  f  = geometry_intersection->m_material->eval(arena, wo, wi, n, sampler);

            if (f != RGB::black() && !occluded(light_sample.m_tester, scene)) {
                L += f * light_sample.m_L * std::abs(dot(wi, n)) / light_sample.m_pdf;
            }
        });

        // TODO: we need to be able to pass material properties to the sample function
        if (depth < scene.max_depth) {
            const auto material_sample = geometry_intersection->m_material->sample(arena, wo, n, sampler);
            if (is_specular(material_sample.properties)) {
                const Ray outgoing_ray{ geometry_intersection->m_point, material_sample.direction };
                L += do_integrate(outgoing_ray, scene, arena, sampler, depth + 1);
            }
        }
    } else if (light_intersection) {
        L += throughput * light_intersection->L;
    }
    return L;
}

BruteForceIntegratorIterativeDynamicRR::BruteForceIntegratorIterativeDynamicRR(int min_depth, int max_depth, int image_width, int image_height)
: m_stats(max_depth - min_depth, Stats2D(image_width, image_height))
{
    for (int depth = min_depth; depth < max_depth; ++depth) {
        for (int y = 0; y < image_height; ++y) {
            for (int x = 0; x < image_width; ++x) {
                assert(get_stats_for_depth(min_depth, depth)(x, y).size() == 0);
            }
        }
    }
}

RGB BruteForceIntegratorIterativeDynamicRR::integrate_impl(const Ray&    ray,
                                                           const Scene&  scene,
                                                           MemoryArena&  arena,
                                                           Sampler&      sampler,
                                                           const Point2& pixel_coords) const
{
    return do_integrate(ray, scene, arena, sampler, pixel_coords);
}

RGB
BruteForceIntegratorIterativeDynamicRR::do_integrate(Ray           ray, const Scene& scene, MemoryArena& arena, Sampler& sampler,
                                                     const Point2& pixel_coords) const
{
    RGB       throughput = RGB::white();
    RGB       L          = RGB::black();
    RayLimits limits;

    constexpr int rr_min_samples = 16;
    for (int depth = 0; depth < scene.max_depth; ++depth) {
        const auto light_intersection = scene.intersect_lights(ray, limits);
        if (light_intersection) {
            limits.m_t_max = light_intersection->m_distance;
        }
        if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
            const Vector3  wo = -ray.get_direction();
            const Normal3& n  = geometry_intersection->m_normal;

            const auto shading_result = geometry_intersection->m_material->sample(arena, wo, n, sampler);
            if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                break;
            }

            // Get next direction
            const auto& wi                = shading_result.direction;
            const auto  cosine            = dot(wi, n);
            const auto  outgoing_position = ray(geometry_intersection->m_distance);
            const Ray   outgoing_ray{ outgoing_position, wi };
            ray            = outgoing_ray;
            limits.m_t_min = get_ray_offset(cosine);
            limits.m_t_max = k_infinite_distance;

            const RGB contribution = cosine * shading_result.color / shading_result.pdf;
            throughput *= contribution;

            if (depth >= scene.russian_roulette_depth) {
                const int pixel_x = static_cast<int>(pixel_coords.x);
                const int pixel_y = static_cast<int>(pixel_coords.y);
                auto&     stats   = get_stats_for_depth(scene.russian_roulette_depth, depth)(pixel_x, pixel_y);
                if (stats.size() >= rr_min_samples) {
                    const float mean_throughput_luminance = stats.mean();
                    const float lum                       = relative_luminance(throughput);
                    if (lum < mean_throughput_luminance) {
                        // q is the probability of continuing
                        const float q = std::max(0.05f, lum / mean_throughput_luminance);
                        assert(q < 1.0f);
                        if (sampler.get_next_1D() < q) {
                            throughput /= q;
                        } else {
                            break;
                        }
                    }
                }
                stats.push(relative_luminance(throughput));
            }
        } else if (light_intersection) {
            L += throughput * light_intersection->L;
            break;
        } else {
            break;
        }
    }
    return L;
}

BruteForceIntegratorIterativeDynamicRR::Stats2D& BruteForceIntegratorIterativeDynamicRR::get_stats_for_depth(int min_depth, int current_depth) const
{
    assert(current_depth >= min_depth);
    return m_stats.at(current_depth - min_depth);
}

RGB estimate_direct(const Scene&    scene,
                    MemoryArena&    arena,
                    const Light&    light,
                    const Point3&   p,
                    const Normal3&  n,
                    const Vector3&  wo,
                    Sampler&        sampler,
                    const Material& material)
{
    const LightSample light_sample = light.sample(p, n, sampler.get_next_2D());
    if (light_sample.m_pdf == 0.0f || light_sample.m_L == RGB::black()) {
        return RGB::black();
    }

    const auto& wi        = light_sample.m_tester.m_ray.get_direction();
    const auto  bsdf_eval = material.eval(arena, wo, wi, n, sampler) * std::abs(dot(n, wi));

    if (bsdf_eval == RGB::black() || scene.intersect_p(light_sample.m_tester.m_ray, light_sample.m_tester.m_limits)) {
        return RGB::black();
    }

    return light_sample.m_L * bsdf_eval / light_sample.m_pdf;
}

RGB estimate_direct_mis(const Scene&    scene,
                        MemoryArena&    arena,
                        const Light&    light,
                        const Point3&   p,
                        const Normal3&  n,
                        const Vector3&  wo,
                        Sampler&        sampler,
                        const Material& material)
{
    auto L_result = RGB::black();

    const auto light_sample = light.sample(p, n, sampler.get_next_2D());
    if (light_sample.m_pdf == 0.0f || light_sample.m_L == RGB::black()) {
        return L_result;
    }

    if (scene.intersect_p(light_sample.m_tester.m_ray, light_sample.m_tester.m_limits)) {
        // Occluded: we're in shadow
        return L_result;
    }

    const auto& wi        = light_sample.m_tester.m_ray.get_direction();
    const auto  bsdf_eval = material.eval(arena, wo, wi, n, sampler);

    if (bsdf_eval != RGB::black()) {
        if (const auto bsdf_pdf = material.pdf(arena, wo, wi, n, sampler); bsdf_pdf > 0.0f) {
            const auto weight = balance_heuristic(1, light_sample.m_pdf, 1, bsdf_pdf);
            L_result += bsdf_eval * light_sample.m_L * (std::abs(dot(wi, n)) * weight / light_sample.m_pdf);
        }
    }

    const MaterialSampleResult material_sample = material.sample(arena, wo, n, sampler);
    if (material_sample.pdf == 0.0f || material_sample.color == RGB::black()) {
        return L_result;
    }

    const float light_pdf = light.pdf(p, material_sample.direction);
    if (light_pdf == 0.0f) {
        return L_result;
    }
    const float weight = balance_heuristic(1, material_sample.pdf, 1, light_pdf);

    RayLimits material_limits;
    material_limits.m_t_min = get_ray_offset(n, material_sample.direction);
    material_limits.m_t_max = k_infinite_distance;
    const Ray material_ray{ p, material_sample.direction };
    if (auto light_intersection = scene.intersect_lights(material_ray, material_limits); light_intersection) {
        if (!scene.intersect_p(material_ray, material_limits)) {
            L_result += material_sample.color * light_intersection->L * std::abs(dot(material_sample.direction, n)) * weight / material_sample.pdf;
        }
    }

    return L_result;
}

RGB IntegratorIterativeRRNEE::integrate_impl(const Ray&   ray,
                                             const Scene& scene,
                                             MemoryArena& arena,
                                             Sampler&     sampler,
                                             const Point2&) const
{
    return do_integrate(ray, scene, arena, sampler);
}

RGB IntegratorIterativeRRNEE::do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler) const
{
    RGB       throughput = RGB::white();
    RGB       L          = RGB::black();
    RayLimits limits;

    constexpr float rr_throughput_luminance_cutoff = 0.1f;
    for (int depth = 0; depth < scene.max_depth; ++depth) {
        const auto light_intersection = scene.intersect_lights(ray, limits);
        if (light_intersection) {
            limits.m_t_max = light_intersection->m_distance;
        }

        if (const auto geometry_intersection = scene.intersect(ray, limits); geometry_intersection) {
            assert(geometry_intersection->m_material);

            const Vector3 wo = -ray.get_direction();
            assert(is_normalized(wo));
            const Normal3& n = geometry_intersection->m_normal;

            const auto shading_result = geometry_intersection->m_material->sample(arena, wo, n, sampler);
            if (shading_result.pdf == 0.0f || shading_result.color == RGB::black()) {
                break;
            }
            assert(is_normalized(shading_result.direction));

#if 0
            scene.for_each_light(
                [&L, &scene, &arena, &throughput, &ray, &n, &wo, &sampler, &geometry_intersection](const Light& light) {
                    L += throughput * estimate_direct(scene,
                                                      arena,
                                                      light,
                                                      geometry_intersection.m_point,
                                                      n,
                                                      wo,
                                                      sampler,
                                                      *geometry_intersection.m_material);
                });
#else
            scene.for_each_light([&L, &scene, &arena, &throughput, &n, &wo, &sampler, &geometry_intersection](
            const Light& light) {
                    L += throughput * estimate_direct_mis(scene,
                                                          arena,
                                                          light,
                                                          geometry_intersection->m_point,
                                                          n,
                                                          wo,
                                                          sampler,
                                                          *geometry_intersection->m_material);
                });
#endif

            // Get next direction
            const auto  outgoing_position = ray(geometry_intersection->m_distance);
            const auto& wi                = shading_result.direction;
            const auto  cosine            = std::abs(dot(wi, n));
            const RGB   contribution      = cosine * shading_result.color / shading_result.pdf;
            throughput *= contribution;

            if (depth >= scene.russian_roulette_depth) {
                const float lum = relative_luminance(throughput);
                if (lum < rr_throughput_luminance_cutoff) {
                    // q is the probability of continuing
                    const float q = std::max(0.05f, lum / rr_throughput_luminance_cutoff);
                    assert(q < 1.0f);
                    if (sampler.get_next_1D() < q) {
                        throughput /= q;
                    } else {
                        break;
                    }
                }
            }

            const Ray outgoing_ray{ outgoing_position, wi };
            ray            = outgoing_ray;
            limits.m_t_min = get_ray_offset(cosine);
            limits.m_t_max = k_infinite_distance;
        } else if (light_intersection) {
            L += throughput * light_intersection->L;
            break;
        } else {
            break;
        }
    }
    return L;
}
} // namespace sp
