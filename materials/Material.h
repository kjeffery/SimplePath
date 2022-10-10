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
    const bool entering = cos_theta_i > 0.0f;
    if (!entering) {
        std::swap(eta_i, eta_t);
        cos_theta_i = std::abs(cos_theta_i);
    }

    // Compute cos_theta_t using Snell’s law
    float sin_theta_i = std::sqrt(std::max(0.0f, 1.0f - cos_theta_i * cos_theta_i));
    float sin_theta_t = eta_i / eta_t * sin_theta_i;
    // Handle total internal reflection
    if (sin_theta_t >= 1) {
        return 1.0f;
    }

    float cos_theta_t = std::sqrt(std::max(0.0f, 1.0f - sin_theta_t * sin_theta_t));

    float real_parl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t)) / ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
    float real_perp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t)) / ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));
    return (real_parl * real_parl + real_perp * real_perp) / 2.0f;
}

struct MaterialSampleResult
{
    RGB     color;
    Vector3 direction;
    float   pdf;
};

class BRDF
{
public:
    virtual ~BRDF() = default;

    [[nodiscard]] MaterialSampleResult sample(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const
    {
        return sample_impl(wo_local, onb_local, sampler);
    }

    [[nodiscard]] float pdf(const Vector3& wo, const Vector3& wi) const
    {
        return pdf_impl(wo, wi);
    }

private:
    virtual MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const = 0;
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
    using BxDFPointer   = std::unique_ptr<BRDF>;
    using BxDFContainer = std::vector<BxDFPointer>;

    // static constexpr int k_max_bxdfs = 8;

    explicit Material(BxDFContainer bxdfs)
    : m_bxdfs(std::move(bxdfs))
    {
    }

    [[nodiscard]] MaterialSampleResult
    sample(const Vector3& wo_world, const Normal3& shading_normal, Sampler& sampler) const
    {
        // std::cout << "Global: " << specular_reflection(wo_world, shading_normal) << '\n';
        const auto onb    = ONB::from_v(Vector3{ shading_normal });
        auto       result = sample_impl(onb.to_onb(wo_world), onb, sampler);
        result.direction  = onb.to_world(result.direction);
        return result;
    }

    [[nodiscard]] float pdf(const Vector3& wo, const Vector3& wi) const
    {
        return 1.0f;
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& onb_local, Sampler& sampler) const
    {
        const std::size_t num_bxdfs = m_bxdfs.size();
        assert(num_bxdfs > 0);

        if (num_bxdfs == 1) {
            return m_bxdfs.front()->sample(wo_local, onb_local, sampler);
        } else {
            // TODO: thread-local memory arena
            std::vector<MaterialSampleResult> results{num_bxdfs};
            std::vector<float>                weights{num_bxdfs};

            // We first use _weights_ to store the potential contributions in order to importance-select our BxDF.

            float weight_sum = 0.0f;
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                results[i] = m_bxdfs[i]->sample(wo_local, onb_local, sampler);
                weights[i] = relative_luminance(results[i].color) / results[i].pdf;
                weight_sum += weights[i];
            }

            // Normalize the weights
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                weights[i] /= weight_sum;
            }

            // Randomly select one BxDF based on sampling each BxDF.
            // Save the results.
            const float u           = sampler.get_next_1D();
            float       running_cdf = 0.0f;
            std::size_t index;
            for (index = 0; index < num_bxdfs; ++index) {
                if (weights[index] + running_cdf > u) {
                    break;
                }
                running_cdf += weights[index];
            }
            assert(index < num_bxdfs);
            if (index >= num_bxdfs) {
                index = num_bxdfs - 1u;
            }

            // Go through each BxDF and calculate the multiple-importance sampling weight.
            // Here we're reusing _weights_ to store the PDFs from the sampling results.

            std::vector<RGB> values{num_bxdfs};
            const float      pdf = results[index].pdf * weights[index];
            for (std::size_t i = 0; i < num_bxdfs; ++i) {
                values[i]  = results[i].color;
                weights[i] = results[i].pdf;
            }

            auto       result     = results[index];
            const auto mis_weight = balance_heuristic(result.color,
                                                      result.pdf,
                                                      values.cbegin(),
                                                      values.cend(),
                                                      weights.cbegin(),
                                                      weights.cend());

            // TODO: do we want to add in all of the other MIS contributions?

            result.pdf   = pdf;
            result.color = mis_weight * result.color;
            return result;
        }
    }

    // std::array<std::unique_ptr<BRDF>, k_max_bxdfs> m_bxdfs;
    BxDFContainer m_bxdfs;
};

inline Material create_lambertian_material(RGB albedo)
{
    Material::BxDFContainer bxdfs;
    bxdfs.emplace_back(new LambertianBRDF{ albedo });
    return Material{ std::move(bxdfs) };
};

inline Material create_clearcoat_material(RGB albedo, RGB reflection = RGB::white())
{
    Material::BxDFContainer bxdfs;
    bxdfs.emplace_back(new SpecularReflectionBRDF{ reflection });
    // bxdfs.emplace_back(new LambertianBRDF{ reflection });
    // bxdfs.emplace_back(new LambertianBRDF{ albedo });
    bxdfs.emplace_back(new LambertianBRDF{ albedo });
    return Material{ std::move(bxdfs) };
};

#endif

} // namespace sp
