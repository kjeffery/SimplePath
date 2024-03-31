#pragma once

/// @author Keith Jeffery

#include "../base/RunningStats.h"
#include "../math/HSV.h"
#include "../math/ONB.h"
#include "../math/Ray.h"
#include "../math/RGB.h"
#include "../math/Sampler.h"
#include "../math/Sampling.h"
#include "../Lights/Light.h"
#include "../base/Array2D.h"

namespace sp {
class Scene;

enum class IntegratorType
{
    NotSpecified,
    Mandelbrot,
    BruteForce,
    BruteForceIterative,
    BruteForceIterativeRR,
    IterativeRRNEE,
    DirectLighting
};

[[nodiscard]] auto string_to_integrator_type(std::string_view s) -> IntegratorType;

class Integrator
{
public:
    virtual ~Integrator() = default;

    [[nodiscard]] RGB integrate(const Ray&    ray,
                                const Scene&  scene,
                                MemoryArena&  arena,
                                Sampler&      sampler,
                                const Point2& pixel_coords) const
    {
        return integrate_impl(ray, scene, arena, sampler, pixel_coords);
    }

private:
    virtual RGB integrate_impl(const Ray&    ray,
                               const Scene&  scene,
                               MemoryArena&  arena,
                               Sampler&      sampler,
                               const Point2& pixel_coords) const = 0;
};

// A simple test integrator that ignores the ray direction and simply uses the x and y values.
class MandelbrotIntegrator final : public Integrator
{
public:
    MandelbrotIntegrator(int image_width, int image_height) noexcept;

private:
    RGB integrate_impl(const Ray& ray,
                       const Scene&,
                       MemoryArena& arena,
                       Sampler&,
                       const Point2& pixel_coords) const override;

    static int mandel(const float c_re, const float c_im) noexcept;

    static constexpr int s_max_iterations = 4096;

    int m_image_width;
    int m_image_height;
};

class BruteForceIntegrator final : public Integrator
{
private:
    RGB integrate_impl(const Ray&   ray,
                       const Scene& scene,
                       MemoryArena& arena,
                       Sampler&     sampler,
                       const Point2&) const override;

    RGB do_integrate(const Ray& ray, const Scene& scene, MemoryArena& arena, Sampler& sampler, int depth) const;
};

class BruteForceIntegratorIterative final : public Integrator
{
private:
    RGB integrate_impl(const Ray&   ray,
                       const Scene& scene,
                       MemoryArena& arena,
                       Sampler&     sampler,
                       const Point2&) const override;

    RGB do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler) const;
};

class BruteForceIntegratorIterativeRR final : public Integrator
{
private:
    RGB integrate_impl(const Ray&   ray,
                       const Scene& scene,
                       MemoryArena& arena,
                       Sampler&     sampler,
                       const Point2&) const override;

    RGB do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler) const;
};

class DirectLightingIntegrator : public Integrator
{
private:
    RGB integrate_impl(const Ray&   ray,
                       const Scene& scene,
                       MemoryArena& arena,
                       Sampler&     sampler,
                       const Point2&) const override;

    RGB do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler, int depth) const;
};

class BruteForceIntegratorIterativeDynamicRR : public Integrator
{
    using StatsType           = RunningStats<float>;
    using Stats2D             = Array2D<StatsType>;
    using StatsDepthContainer = std::vector<Stats2D>;

public:
    BruteForceIntegratorIterativeDynamicRR(int min_depth, int max_depth, int image_width, int image_height);

private:
    RGB integrate_impl(const Ray&    ray,
                       const Scene&  scene,
                       MemoryArena&  arena,
                       Sampler&      sampler,
                       const Point2& pixel_coords) const override;

    RGB      do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler, const Point2& pixel_coords) const;
    Stats2D& get_stats_for_depth(int min_depth, int current_depth) const;

    mutable StatsDepthContainer m_stats;
};

class IntegratorIterativeRRNEE final : public Integrator
{
private:
    RGB integrate_impl(const Ray&   ray,
                       const Scene& scene,
                       MemoryArena& arena,
                       Sampler&     sampler,
                       const Point2&) const override;

    RGB do_integrate(Ray ray, const Scene& scene, MemoryArena& arena, Sampler& sampler) const;
};
} // namespace sp
