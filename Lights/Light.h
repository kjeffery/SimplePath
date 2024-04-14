#pragma once

/// @author Keith Jeffery

#include <utility>

#include "../Image/Image.h"
#include "../math/Distribution2D.h"
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

    [[nodiscard]] auto pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept -> float override
    {
        return shape_pdf(observer_world, wi);
    }

    [[nodiscard]] auto L_impl(const Point3& observer_world, const Normal3& n, const Vector3& wi) const noexcept -> RGB override
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
    std::optional<Intersection> intersect_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(!"Should not get here");
        return {};
    }

    std::optional<LightIntersection> intersect_lights_impl(const Ray&, const RayLimits& limits) const noexcept override
    {
        if (limits.m_t_max < k_infinite_distance) {
            return {};
        }
        return { LightIntersection{ .m_distance = k_infinite_distance, .L = m_radiance } };
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

    [[nodiscard]] auto pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept -> float override
    {
        return uniform_sphere_pdf();
    }

    [[nodiscard]] auto L_impl(const Point3& observer_world, const Normal3& n, const Vector3& wi) const noexcept -> RGB override
    {
        return (dot(n, wi) > 0.0f) ? m_radiance : RGB::black();
    }

    RGB m_radiance;
};

class ImageBasedEnvironmentLight final : public InfiniteLight
{
public:
    ImageBasedEnvironmentLight(Image radiance, float max_radiance) noexcept
    : m_radiance{ modify_image(std::move(radiance), max_radiance) }
    , m_distribution_2d{ create_distribution(m_radiance, max_radiance) }
    {
    }

private:
    std::optional<Intersection> intersect_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(!"Should not get here");
        return {};
    }

    std::optional<LightIntersection> intersect_lights_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        if (limits.m_t_max < k_infinite_distance) {
            return {};
        }

        constexpr auto inv_2_pi = 1.0f / (2.0f * std::numbers::pi_v<float>);

        const auto   w = normalize(m_world_to_light(ray.get_direction()));
        const Point2 st{ spherical_phi(w) * inv_2_pi, spherical_theta(w) * std::numbers::inv_pi_v<float> };
        const auto   L = sample_nearest_neighbor(m_radiance, st[0], st[1], RemapWrap{}, RemapClamp{});
        return { LightIntersection{ .m_distance = k_infinite_distance, .L = L } };
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
        const auto [uv, map_pdf] = m_distribution_2d.sample_continuous(u);
        if (map_pdf == 0.0f) {
            return { .m_pdf = 0.0f, .m_max_distance = k_infinite_distance, .m_L = RGB::black(), .m_wi = Vector3{ no_init } };
        }

        // Convert infinite light sample point to direction
        const auto theta     = uv[1] * std::numbers::pi_v<float>;
        const auto phi       = uv[0] * 2.0f * std::numbers::pi_v<float>;
        const auto cos_theta = std::cos(theta);
        const auto sin_theta = std::sin(theta);
        const auto sin_phi   = std::sin(phi);
        const auto cos_phi   = std::cos(phi);

        const auto wi = m_light_to_world(Vector3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi));

        // Compute PDF for sampled infinite light direction
        const auto pdf = (sin_theta == 0.0f) ? 0.0f : map_pdf / (2.0f * square(std::numbers::pi_v<float>) * sin_theta);

        const auto L = sample_nearest_neighbor(m_radiance, uv[0], uv[1], RemapWrap{}, RemapClamp{});
        return { .m_pdf = pdf, .m_max_distance = k_infinite_distance, .m_L = L, .m_wi = wi };
    }

    [[nodiscard]] auto pdf_impl(const Point3& observer_world, const Vector3& wi) const noexcept -> float override
    {
        const auto w         = m_world_to_light(wi);
        const auto theta     = spherical_theta(w);
        const auto phi       = spherical_phi(w);
        const auto sin_theta = std::sin(theta);
        if (sin_theta == 0.0f) {
            return 0.0f;
        }
        constexpr auto inv_2_pi = 1.0f / (2.0f * std::numbers::pi_v<float>);

        const auto pdf = m_distribution_2d.pdf({ phi * inv_2_pi, theta * std::numbers::pi_v<float> }) / (2.0f * square(std::numbers::pi_v<float>) *
            sin_theta);
        return pdf;
    }

    [[nodiscard]] auto L_impl(const Point3& observer_world, const Normal3& n, const Vector3& wi) const noexcept -> RGB override
    {
        if (dot(n, wi) > 0.0f) {
            return RGB::black();
        }
        constexpr auto inv_2_pi = 1.0f / (2.0f * std::numbers::pi_v<float>);

        const auto   w = normalize(m_world_to_light(wi));
        const Point2 st{ spherical_phi(w) * inv_2_pi, spherical_theta(w) * std::numbers::inv_pi_v<float> };
        return sample_nearest_neighbor(m_radiance, st[0], st[1], RemapWrap{}, RemapClamp{});
    }

    static [[nodiscard]] auto modify_image(Image radiance, float max_radiance) -> Image
    {
        using size_type = Image::size_type;
        for (size_type y = 0; y < radiance.height(); ++y) {
            for (size_type x = 0; x < radiance.width(); ++x) {
                auto& color = radiance(x, y);
                for (std::uint32_t i = 0; i < 3; ++i) {
                    if (std::isinf(color[i])) {
                        color[i] = max_radiance;
                    }
                }
                if (relative_luminance(color) > max_radiance) {
                    const auto max_index = index_of_max(color);
                    for (std::uint32_t i = 0; i < 3; ++i) {
                        color[i] = color[i] * max_radiance / color[max_index];
                    }
                }
            }
        }
        return radiance;
    }

    static [[nodiscard]] auto create_distribution(const Image& radiance, float max_radiance) -> Distribution2D
    {
        const auto         width  = 2 * radiance.width();
        const auto         height = 2 * radiance.height();
        std::vector<float> img(width * height);

        float max = std::numeric_limits<float>::lowest();
        for (auto v = decltype(height){ 0 }; v < height; ++v) {
            const auto vp        = (static_cast<float>(v) + 0.5f) / static_cast<float>(height);
            const auto sin_theta = std::sin(std::numbers::pi_v<float> * (static_cast<float>(v) + 0.5f) / static_cast<float>(height));
            for (auto u = decltype(width){ 0 }; u < width; ++u) {
                const auto up      = (static_cast<float>(u) + 0.5f) / static_cast<float>(width);
                img[u + v * width] = relative_luminance(sample_nearest_neighbor(radiance, up, vp, RemapWrap{}, RemapClamp{}));
                img[u + v * width] *= sin_theta;
                if (std::isinf(img[u + v * width])) {
                    img[u + v * width] = max_radiance;
                }
                img[u + v * width] = std::min(img[u + v * width], max_radiance);
                max                = std::max(max, img[u + v * width]);
            }
        }
#if 0
        const auto mean = std::accumulate(img.cbegin(), img.cend(), 0.0f) / static_cast<float>(img.size());
        std::transform(std::execution::par_unseq, img.cbegin(), img.cend(), img.begin(), [mean](const auto f) { return std::max(f - mean, 0.0f); });
#endif

        return Distribution2D{ img, width, height };
    }

    // TODO: user-set
    LinearSpace3x3 m_light_to_world{ LinearSpace3x3::identity() };
    LinearSpace3x3 m_world_to_light{ LinearSpace3x3::identity() };
    Image          m_radiance;
    Distribution2D m_distribution_2d;
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
    std::optional<Intersection> intersect_impl(const Ray&, const RayLimits&) const noexcept override
    {
        assert(!"Should not get here");
        return {};
    }

    std::optional<LightIntersection> intersect_lights_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        // TODO: may need normal
        if (const auto geom_isect = m_sphere.intersect(ray, limits); geom_isect) {
            return { LightIntersection{ .m_distance = geom_isect->m_distance, .L = get_radiance() } };
        }
        return {};
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
