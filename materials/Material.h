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

#if defined(_MSC_VER)
#include <iostream>

#include <malloc.h>
#else
#include <alloca.h>
#endif

namespace sp {
inline Vector3 specular_reflection(const Vector3& wo, const Normal3& n) noexcept
{
    return -wo + 2.0f * dot(wo, n) * n;
}

inline Vector3 specular_reflection_local(const Vector3& wo) noexcept
{
    // In a local y-up space:
    return Vector3{ -wo.x, wo.y, -wo.z };
}

inline float cos_theta(const Vector3& w)
{
    return w.y;
}

inline float cos2_theta(const Vector3& w)
{
    return w.y * w.y;
}

inline float abs_cos_theta(const Vector3& w)
{
    return std::abs(w.y);
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
};

struct MaterialEvalResult
{
    RGB   color;
    float pdf;
};

class BRDF
{
public:
    virtual ~BRDF() = default;

    [[nodiscard]] MaterialSampleResult sample(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const
    {
        return sample_impl(wo_local, onb_local, sampler);
    }

    [[nodiscard]] MaterialEvalResult eval(const Vector3& wo_local, const Vector3& wi_local) const
    {
        return eval_impl(wo_local, wi_local);
    }

    [[nodiscard]] float pdf(const Vector3& wo_local, const Vector3& wi_local) const
    {
        return pdf_impl(wo_local, wi_local);
    }

private:
    virtual MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const = 0;
    virtual MaterialEvalResult   eval_impl(const Vector3& wo_local, const Vector3& wi_local) const                  = 0;
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

    MaterialEvalResult eval_impl(const Vector3& wo_local, const Vector3& wi_local) const override
    {
        // TODO: assert or check they are in the same hemisphere.
        return MaterialEvalResult{ m_albedo, uniform_hemisphere_pdf() };
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
        // std::cout << "Local: " << onb_local.to_world(wi) << '\n';
        // const RGB  color = m_r;
        return { color, wi, 1.0f };
    }

    MaterialEvalResult eval_impl(const Vector3& wo_local, const Vector3& wi_local) const override
    {
        return MaterialEvalResult{ RGB::black(), 0.0f };
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

#if 0
class Material
{
public:
    virtual ~Material() = default;

    [[nodiscard]] MaterialSampleResult sample(const Vector3& wo_world, const Normal3& shading_normal, Sampler& sampler) const
    {
        const auto onb    = ONB::from_w(Vector3{ shading_normal });
        auto       result = sample_impl(onb.to_onb(wo_world), onb, sampler);
        result.direction  = onb.to_world(result.direction);
        return result;
    }

    [[nodiscard]] float pdf(const Vector3& wo, const Vector3& wi) const
    {
        return pdf_impl(wo, wi);
    }

private:
    virtual MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) = 0;
    virtual float                pdf_impl(const Vector3& wo, const Vector3& wi) const                         = 0;
};

class LambertianMaterial : public Material
{
public:
    explicit LambertianMaterial(RGB albedo) noexcept
    : m_brdf(albedo)
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        return m_brdf.sample(wo_local, onb_local, sampler);
    }

    float pdf_impl(const Vector3& wo, const Vector3& wi) const override
    {
        return m_brdf.pdf(wo, wi);
    }

    LambertianBRDF m_brdf;
};
#else
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
                const auto eval_result = m_bxdfs[i]->eval(wo_local, results[selected_index].direction);
                values[i]              = eval_result.color;
                weights[i]             = eval_result.pdf;
            }

            // We're using one sample for each type, so our inner product is just the sum of all weights (as opposed to
            // the inner product of weights and sample counts).
            const float inner_product = std::reduce(weights.cbegin(), weights.cend(), 0.0f);

            // As far as I can tell, in the paper, they don't add in the contributions from the additional sampling
            // techniques.
#if 0
            const auto mis_weight = balance_heuristic(weights[selected_index], inner_product);
            const RGB result_color = mis_weight * values[selected_index];
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

class ClearcoatMaterial : public Material
{
public:
    // TODO: add IOR and clearcoat color
    explicit ClearcoatMaterial(std::unique_ptr<Material> base)
    //: m_specular(RGB::white())
    : m_base(std::move(base))
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const override
    {
        // We will sample our top-level specular_result material. Using that result, we will probabilistically decide to
        // sample the underlying material, subtracting the contribution of the top-level specular_result.

        // TODO: user parameter
        constexpr RGB m_specular_color{1.0f, 1.0f, 1.0f};

        const float f = fresnel_dielectric(cos_theta(wo_local), 1.0f, 1.5f);
        // const auto specular_result = m_specular.sample(wo_local, onb_local, sampler);
        // assert(specular_result.pdf == 1.0f);
        // const float weight = std::min(1.0f, relative_luminance(specular_result.color));

        const auto  specular_wi{ specular_reflection_local(wo_local) };
        const RGB   specular_color_result = f * m_specular_color / abs_cos_theta(specular_wi);

        const float u = sampler.get_next_1D();
        if (u < f) {
            const float specular_pdf = f; // Our old pdf is 1.0, so this is a multiplication against 1.
            return MaterialSampleResult{ specular_color_result, specular_wi, specular_pdf };

            // return MaterialSampleResult{ specular_result.color, specular_result.direction, result_pdf };
        }

        const auto base_result = m_base->sample_local_space(wo_local, onb_local, sampler);
        if (base_result.pdf == 0.0f) {
            return base_result;
        }

        // TODO: do I need to modify the pdf?
        const float result_pdf = (1.0f - f) * base_result.pdf;
        const RGB   color      = (RGB::white() - f * m_specular_color) * base_result.color;
        return MaterialSampleResult{ color, base_result.direction, result_pdf };
    };

    // TODO: implement when needed
    float pdf_impl(const Vector3& wo, const Vector3& wi) const override
    {
        return 1.0f;
    }

    // SpecularReflectionBRDF    m_specular;
    std::unique_ptr<Material> m_base;
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
    //bxdfs.emplace_back(new SpecularReflectionBRDF{ albedo });
    return OneSampleMaterial{ std::move(bxdfs) };
};

inline ClearcoatMaterial create_clearcoat_material(RGB albedo, RGB reflection = RGB::white())
{
#if 0
    OneSampleMaterial::BxDFContainer bxdfs;
    bxdfs.emplace_back(new SpecularReflectionBRDF{ reflection });
    bxdfs.emplace_back(new LambertianBRDF{ albedo });
    return OneSampleMaterial{ std::move(bxdfs) };
#else
    auto base = std::make_unique<LambertianMaterial>(albedo);
    return ClearcoatMaterial{ std::move(base) };
#endif
};

#endif

} // namespace sp
