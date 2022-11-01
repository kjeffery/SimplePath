//
// Created by krjef on 9/18/2022.
//

#pragma once

/// @author Keith Jeffery

#include "../math/ONB.h"
#include "../math/RGB.h"
#include "../math/Sampler.h"
#include "../math/Sampling.h"
#include "../math/Vector3.h"

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <utility>

namespace sp {
constexpr bool same_hemisphere(const Vector3& a, const Vector3& b)
{
    return a.y * b.y > 0.0f;
}

constexpr bool same_hemisphere(const Vector3& a, const Normal3& b)
{
    return a.y * b.y > 0.0f;
}

constexpr bool same_hemisphere(const Normal3& a, const Vector3& b)
{
    return a.y * b.y > 0.0f;
}

inline Vector3 specular_reflection(const Vector3& wo, const Normal3& n) noexcept
{
    return -wo + 2.0f * dot(wo, n) * n;
}

inline Vector3 specular_reflection_local(const Vector3& wo) noexcept
{
    // In a local y-up space:
    return Vector3{ -wo.x, wo.y, -wo.z };
}

inline float cos_theta(const Vector3& w) noexcept
{
    return w.y;
}

inline float cos2_theta(const Vector3& w) noexcept
{
    return w.y * w.y;
}

inline float abs_cos_theta(const Vector3& w) noexcept
{
    return std::abs(w.y);
}

inline float sin2_theta(const Vector3& w) noexcept
{
    return std::max(0.0f, 1.0f - cos2_theta(w));
}

inline float sin_theta(const Vector3& w) noexcept
{
    return std::sqrt(sin2_theta(w));
}

inline float tan_theta(const Vector3& w) noexcept
{
    return sin_theta(w) / cos_theta(w);
}

inline float tan2_theta(const Vector3& w) noexcept
{
    return sin2_theta(w) / cos2_theta(w);
}

inline float cos_phi(const Vector3& w) noexcept
{
    const float sin_theta_v = sin_theta(w);
    return (sin_theta_v == 0.0f) ? 1.0f : std::clamp(w.x / sin_theta_v, -1.0f, 1.0f);
}

inline float sin_phi(const Vector3& w) noexcept
{
    const float sin_theta_v = sin_theta(w);
    return (sin_theta_v == 0.0f) ? 1.0f : std::clamp(w.z / sin_theta_v, -1.0f, 1.0f);
}

inline float cos2_phi(const Vector3& w) noexcept
{
    return square(cos_phi(w));
}

inline float sin2_phi(const Vector3& w) noexcept
{
    return square(sin_phi(w));
}

// This is taken near-verbatim from PBRT
inline float fresnel_dielectric(float cos_theta_i, float eta_i, float eta_t)
{
    cos_theta_i = std::clamp(cos_theta_i, -1.0f, 1.0f);
    // Potentially swap indices of refraction
    if (const bool entering = cos_theta_i > 0.0f; !entering) {
        std::swap(eta_i, eta_t);
        cos_theta_i = std::abs(cos_theta_i);
    }

    // Compute cos_theta_t using Snell’s law
    const float sin_theta_i = std::sqrt(std::max(0.0f, 1.0f - cos_theta_i * cos_theta_i));
    const float sin_theta_t = eta_i / eta_t * sin_theta_i;
    // Handle total internal reflection
    if (sin_theta_t >= 1) {
        return 1.0f;
    }

    const float cos_theta_t = std::sqrt(std::max(0.0f, 1.0f - sin_theta_t * sin_theta_t));

    assert(cos_theta_i >= 0.0f);
    assert(cos_theta_t >= 0.0f);

    // clang-format off
    const float real_parl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t)) /
                            ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
    const float real_perp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t)) /
                            ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));
    // clang-format on
    return (real_parl * real_parl + real_perp * real_perp) / 2.0f;
}

struct MaterialSampleResult
{
    RGB     color;
    Vector3 direction;
    float   pdf;

    static MaterialSampleResult degenerate() noexcept
    {
        return MaterialSampleResult{ RGB::black(), Vector3{ no_init }, 0.0f };
    }
};

// These come from PBRT:
// https://www.pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models#MicrofacetDistribution
class MicrofacetDistribution
{
public:
    virtual ~MicrofacetDistribution() = default;

    [[nodiscard]] float D(const Vector3& wh) const
    {
        return D_impl(wh);
    }

    [[nodiscard]] float G1(const Vector3& w) const
    {
        return 1.0f / (1.0f + lambda(w));
    }

    [[nodiscard]] float G(const Vector3& wo, const Vector3& wi) const
    {
        return 1.0f / (1.0f + lambda(wo) + lambda(wi));
    }

    [[nodiscard]] Vector3 sample_wh(const Vector3& wo, Sampler& sampler) const
    {
        return sample_wh_impl(wo, sampler);
    }

    [[nodiscard]] float pdf(const Vector3& wo, const Vector3& wh) const
    {
        if (m_sample_visible_area) {
            return D(wh) * G1(wo) * std::abs(dot(wo, wh)) / abs_cos_theta(wo);
        } else {
            return D(wh) * abs_cos_theta(wh);
        }
    }

protected:
    explicit MicrofacetDistribution(bool sample_visible_area) noexcept
    : m_sample_visible_area(sample_visible_area)
    {
    }

    bool sample_visible_area() const noexcept
    {
        return m_sample_visible_area;
    }

private:
    virtual float   D_impl(const Vector3& wh) const                           = 0;
    virtual float   lambda(const Vector3& w) const                            = 0;
    virtual Vector3 sample_wh_impl(const Vector3& wo, Sampler& sampler) const = 0;

    bool m_sample_visible_area;
};

class BeckmannDistribution : public MicrofacetDistribution
{
public:
    explicit BeckmannDistribution(float roughness, bool sample_visible_area = true)
    : MicrofacetDistribution(sample_visible_area)
    , m_alpha_x(roughness_to_alpha(roughness))
    , m_alpha_y(roughness_to_alpha(roughness))
    {
    }

    BeckmannDistribution(float roughness0, float roughness1, bool sample_visible_area = true)
    : MicrofacetDistribution(sample_visible_area)
    , m_alpha_x(roughness_to_alpha(roughness0))
    , m_alpha_y(roughness_to_alpha(roughness1))
    {
    }

private:
    static float roughness_to_alpha(float roughness)
    {
        roughness     = std::max(roughness, 1e-3f);
        const float x = std::log(roughness);
        return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }

    float D_impl(const Vector3& wh) const override
    {
        const float tan2_theta_v = tan2_theta(wh);
        if (std::isinf(tan2_theta_v)) {
            return 0.0f;
        }
        const float cos4_theta_v = square(cos2_theta(wh));
        return std::exp(-tan2_theta_v * (cos2_phi(wh) / square(m_alpha_x) + sin2_phi(wh) / square(m_alpha_y))) /
               (std::numbers::pi_v<float> * m_alpha_x * m_alpha_y * cos4_theta_v);
    }

    float lambda(const Vector3& w) const override
    {
        const float abs_tan_theta_v = std::abs(tan_theta(w));
        if (std::isinf(abs_tan_theta_v)) {
            return 0.0f;
        }
        const float alpha = std::sqrt(cos2_phi(w) * square(m_alpha_x) + sin2_phi(w) * square(m_alpha_y));
        const float a     = 1.0f / (alpha * abs_tan_theta_v);
        if (a >= 1.6f) {
            return 0.0f;
        }
        return (1.0f - 1.259f * a + 0.396f * square(a)) / (3.535f * a + 2.181f * square(a));
    }

    Vector3 sample_wh_impl(const Vector3& wo, Sampler& sampler) const override;

    float m_alpha_x;
    float m_alpha_y;
};

class BRDF
{
public:
    virtual ~BRDF() = default;

    [[nodiscard]] MaterialSampleResult sample(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const
    {
        return sample_impl(wo_local, onb_local, sampler);
    }

    [[nodiscard]] RGB eval(const Vector3& wo_local, const Vector3& wi_local) const
    {
        return eval_impl(wo_local, wi_local);
    }

    [[nodiscard]] float pdf(const Vector3& wo_local, const Vector3& wi_local) const
    {
        return pdf_impl(wo_local, wi_local);
    }

private:
    virtual MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const = 0;
    virtual RGB                  eval_impl(const Vector3& wo_local, const Vector3& wi_local) const                  = 0;
    virtual float                pdf_impl(const Vector3& wo, const Vector3& wi) const                               = 0;
};

class LambertianBRDF : public BRDF
{
public:
    explicit LambertianBRDF(RGB albedo) noexcept
    : m_albedo(albedo / std::numbers::pi_v<float>)
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        // TODO: cosine-weighted hemispherical sampling
        const auto sp = sample_to_hemisphere(sampler.get_next_2D());
        return { m_albedo, Vector3{ sp }, uniform_hemisphere_pdf() };
    }

    RGB eval_impl(const Vector3& wo_local, const Vector3& wi_local) const override
    {
        assert(same_hemisphere(wo_local, wi_local));
        return m_albedo;
    }

    float pdf_impl(const Vector3&, const Vector3&) const override
    {
        return uniform_hemisphere_pdf();
    }

    RGB m_albedo;
};

class SpecularReflectionBRDF : public BRDF
{
public:
    explicit SpecularReflectionBRDF(RGB r) noexcept
    : m_r(r)
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        const auto wi    = specular_reflection_local(wo_local);
        const RGB  color = fresnel_dielectric(cos_theta(wi), 1.0f, 1.5f) * m_r / abs_cos_theta(wi);
        return { color, wi, 1.0f };
    }

    RGB eval_impl(const Vector3& wo_local, const Vector3& wi_local) const override
    {
        return RGB::black();
    }

    float pdf_impl(const Vector3&, const Vector3&) const override
    {
        // We're a delta singularity: statistically, the probability of random input and output vectors resulting in a
        // valid reflection is 0.
        return 0;
    }

private:
    RGB m_r;
};

// This is based on PBRT's Torrance-Sparrow implementation
class MicrofacetReflection : public BRDF
{
public:
    MicrofacetReflection(RGB r, std::unique_ptr<MicrofacetDistribution> distribution, float ior)
    : m_r(r)
    , m_distribution(std::move(distribution))
    , m_ior(ior)
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        if (wo_local.y == 0.0f) {
            return MaterialSampleResult::degenerate();
        }
        const auto  wh_local    = m_distribution->sample_wh(wo_local, sampler);
        const float dot_product = dot(wo_local, wh_local);
        if (dot_product < 0.0f) {
            return MaterialSampleResult::degenerate();
        }

        const auto wi_local = specular_reflection(wo_local, wh_local);
        if (!same_hemisphere(wo_local, wi_local)) {
            return MaterialSampleResult::degenerate();
        }

        const float pdf   = m_distribution->pdf(wo_local, wh_local) / (4.0f * dot_product);
        const auto  color = eval(wo_local, wi_local);
        return MaterialSampleResult{ color, wi_local, pdf };
    }

    RGB eval_impl(const Vector3& wo_local, const Vector3& wi_local) const override
    {
        const float abs_cos_theta_o = abs_cos_theta(wo_local);
        const float abs_cos_theta_i = abs_cos_theta(wi_local);

        if (abs_cos_theta_i == 0.0f || abs_cos_theta_o == 0.0f) {
            return RGB::black();
        }

        Vector3 wh_local = wi_local + wo_local;
        if (wh_local.x == 0.0f && wh_local.y == 0.0f && wh_local.z == 0.0f) {
            return RGB::black();
        }
        wh_local     = normalize(wh_local);
        const auto f = fresnel_dielectric(dot(wi_local, wh_local), 1.0f, m_ior);
        return m_r * m_distribution->D(wh_local) * m_distribution->G(wo_local, wi_local) * f /
               (4.0f * abs_cos_theta_i * abs_cos_theta_o);
    }

    float pdf_impl(const Vector3& wo, const Vector3& wi) const override
    {
        if (!same_hemisphere(wo, wi)) {
            return 0.0f;
        }
        const auto wh = normalize(wo + wi);
        return m_distribution->pdf(wo, wh) / (4.0f * dot(wo, wh));
    }

    RGB                                     m_r;
    std::unique_ptr<MicrofacetDistribution> m_distribution;
    float m_ior; // We only model reflection right now, so, unlike PBRT, we don't have a general Fresnel class
};

class Material
{
public:
    virtual ~Material() = default;

    [[nodiscard]] MaterialSampleResult
    sample(const Vector3& wo_world, const Normal3& shading_normal, Sampler& sampler) const
    {
        const auto onb    = ONB::from_v(Vector3{ shading_normal });
        auto       result = sample_impl(onb.to_onb(wo_world), onb, sampler);
        result.direction  = onb.to_world(result.direction);
        return result;
    }

    [[nodiscard]] float pdf(const Vector3& wo, const Vector3& wi) const
    {
        return pdf_impl(wo, wi);
    }

    [[nodiscard]] MaterialSampleResult
    sample_local_space(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const
    {
        return sample_impl(wo_local, onb_local, sampler);
    }

private:
    virtual MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const = 0;
    virtual float                pdf_impl(const Vector3& wo, const Vector3& wi) const                               = 0;
};

// This is based on Veach and Guibas' multiple importance sampling. It's a general material used to hold any number of
// BxDFs
class OneSampleMaterial : public Material
{
public:
    using BxDFPointer   = std::unique_ptr<BRDF>;
    using BxDFContainer = std::vector<BxDFPointer>;

    explicit OneSampleMaterial(BxDFContainer bxdfs)
    : m_bxdfs(std::move(bxdfs))
    {
    }

private:
    // This implements Veach and Guibas' one-sample model from _Optimally Combining Sampling Techniques
    // for Monte Carlo Rendering_.
    // Eric Veach and Leonidas J. Guibas
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        const std::size_t num_bxdfs = m_bxdfs.size();
        assert(num_bxdfs > 0);

        if (num_bxdfs == 1) {
            return m_bxdfs.front()->sample(wo_local, onb_local, sampler);
        } else {
            // TODO: thread-local memory arena
            std::vector<MaterialSampleResult> results(num_bxdfs);
            std::vector<float>                weights(num_bxdfs);

            // We first use _weights_ to store the potential contributions in order to importance-select our BxDF.
            // We sample each BxDF to find a potential contribution from that BxDF. We will only end up using one of
            // these values, but use the others to importance-sample our primary BxDF.

            // We keep a running total of the weights so that we can normalize them into a proper discrete probability
            // mass function.
            float weight_sum = 0.0f;
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                results[i] = m_bxdfs[i]->sample(wo_local, onb_local, sampler);
                if (results[i].pdf == 0.0f) {
                    weights[i] = 0.0f;
                    continue;
                }
                weights[i] = relative_luminance(results[i].color) / results[i].pdf;
                weight_sum += weights[i];
            }

            if (weight_sum == 0.0f) {
                return results.front();
            }

            // Normalize the weights
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                weights[i] /= weight_sum;
            }

            // Randomly select one BxDF based on sampling each BxDF.
            // Save the results.
            const float u           = sampler.get_next_1D();
            float       running_cdf = 0.0f;
            std::size_t selected_index;
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                if (weights[i] + running_cdf > u) {
                    selected_index = i;
                    break;
                }
                running_cdf += weights[i];
            }

            // In pure mathematics, this wouldn't happen, but we're using IEEE-754 ;)
            // Account for numerical precision errors.
            selected_index = std::min(num_bxdfs - 1u, selected_index);

            // Go through each BxDF and calculate the multiple-importance sampling weight.
            // Here we're reusing _weights_ to store the PDFs from the sampling results.

            std::vector<RGB> values(num_bxdfs);

            // Save off our resulting PDF before we overwrite weights...
            const float result_pdf = results[selected_index].pdf * weights[selected_index];

            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                // This isn't just for efficiency: if we sample a perfectly-specular BxDF, our PDF will be zero out of
                // the eval function, which will give us incorrect results.
                if (i == selected_index) {
                    values[i]  = results[i].color;
                    weights[i] = results[i].pdf;
                    continue;
                }
                const auto eval_color = m_bxdfs[i]->eval(wo_local, results[selected_index].direction);
                const auto eval_pdf   = m_bxdfs[i]->pdf(wo_local, results[selected_index].direction);
                values[i]             = eval_color;
                weights[i]            = eval_pdf;
            }

            // We're using one sample for each type, so our inner product is just the sum of all weights (as opposed to
            // the inner product of weights and sample counts).
            const float inner_product = std::reduce(weights.cbegin(), weights.cend(), 0.0f);

            // As far as I can tell, in the paper, they don't add in the contributions from the additional sampling
            // techniques.
#if 1
            const auto mis_weight   = balance_heuristic(weights[selected_index], inner_product);
            const RGB  result_color = mis_weight * values[selected_index];
#else
            RGB result_color = RGB::black();
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                const auto mis_weight = balance_heuristic(weights[i], inner_product);
                result_color += mis_weight * values[i];
            }
#endif

            // C++ esoterica: we deliberately return a temporary here (instead of copying a MaterialSampleResult from
            // above) so that we can use return value optimization (RVO). We return a temporary above, so we want to be
            // consistent.
            return MaterialSampleResult{ result_color, results[selected_index].direction, result_pdf };
        }
    }

    // TODO: implement when needed
    float pdf_impl(const Vector3& wo, const Vector3& wi) const override
    {
        return 1.0f;
    }

    BxDFContainer m_bxdfs;
};

// This is a simple layered material. The top layer is a specular reflection. This doesn't account for refraction or
// depth of the clear coat, but it does physically-plausible energy conservation (any energy lost from the specular
// reflection is removed from the base material).
class ClearcoatMaterial : public Material
{
public:
    explicit ClearcoatMaterial(std::shared_ptr<Material> base, float ior, RGB specular_color = RGB::white())
    : m_ior(ior)
    , m_specular_color(specular_color)
    , m_base(std::move(base))
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        // We will sample our top-level specular material. Using that result, we will probabilistically decide to sample
        // the underlying material, subtracting the contribution of the top-level specular_result.

        constexpr float ior_air = 1.0f;

        const float f = fresnel_dielectric(cos_theta(wo_local), ior_air, m_ior);
        assert(f >= 0.0f);
        assert(f <= 1.0f);

        const auto specular_wi{ specular_reflection_local(wo_local) };
        const RGB  specular_color_result = f * m_specular_color / abs_cos_theta(specular_wi);

        // Importance sample the specular layer based on the Fresnel contribution.
        const float u = sampler.get_next_1D();
        if (u < f) {
            const float specular_pdf = f; // Our old pdf is 1.0, so this is a multiplication against 1.
            return MaterialSampleResult{ specular_color_result, specular_wi, specular_pdf };
        }

        const auto base_result = m_base->sample_local_space(wo_local, onb_local, sampler);
        if (base_result.pdf == 0.0f) {
            return base_result;
        }

        const float result_pdf = (1.0f - f) * base_result.pdf;
        const RGB   color      = (RGB::white() - f * m_specular_color) * base_result.color;
        assert(color.r >= 0.0f);
        assert(color.g >= 0.0f);
        assert(color.b >= 0.0f);
        return MaterialSampleResult{ color, base_result.direction, result_pdf };
    };

    // TODO: implement when needed
    float pdf_impl(const Vector3& wo, const Vector3& wi) const override
    {
        return 1.0f;
    }

    float                     m_ior;
    RGB                       m_specular_color;
    std::shared_ptr<Material> m_base;
};

class LambertianMaterial : public OneSampleMaterial
{
public:
    explicit LambertianMaterial(RGB albedo)
    : OneSampleMaterial{ build(albedo) }
    {
    }

private:
    static OneSampleMaterial::BxDFContainer build(RGB albedo)
    {
        OneSampleMaterial::BxDFContainer bxdfs;
        bxdfs.emplace_back(new LambertianBRDF{ albedo });
        return bxdfs;
    }
};

inline OneSampleMaterial create_lambertian_material(RGB albedo)
{
    OneSampleMaterial::BxDFContainer bxdfs;
    bxdfs.emplace_back(new LambertianBRDF{ albedo });
    return OneSampleMaterial{ std::move(bxdfs) };
};

inline ClearcoatMaterial
create_clearcoat_material(std::shared_ptr<Material> base, float ior, RGB reflection = RGB::white())
{
    return ClearcoatMaterial{ base, ior, reflection };
};

// TODO: anisotropic version
inline OneSampleMaterial create_beckmann_glossy_material(RGB color, float roughness, float ior)
{
    OneSampleMaterial::BxDFContainer        bxdfs;
    std::unique_ptr<MicrofacetDistribution> microfacet(new BeckmannDistribution{ roughness });
    bxdfs.emplace_back(new MicrofacetReflection{ RGB::white(), std::move(microfacet), ior });
    bxdfs.emplace_back(new LambertianBRDF{ color });
    return OneSampleMaterial{ std::move(bxdfs) };
}

} // namespace sp
