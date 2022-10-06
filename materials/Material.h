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

#include <utility>

namespace sp {

struct MaterialSampleResult
{
    RGB     color;
    Vector3 direction;
    float   pdf;
};

class Material
{
public:
    virtual ~Material() = default;

    [[nodiscard]] MaterialSampleResult sample(const Vector3& wo_world, const Normal3& shading_normal, Sampler& sampler)
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
    virtual MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& local_onb, Sampler& sampler) = 0;
    virtual float                pdf_impl(const Vector3& wo, const Vector3& wi) const                         = 0;
};

class LambertianMaterial : public Material
{
public:
    explicit LambertianMaterial(RGB albedo) noexcept
    : m_albedo(albedo / std::numbers::pi_v<float>)
    {
    }

private:
    MaterialSampleResult sample_impl(const Vector3& wo_local, const ONB& local_onb, Sampler& sampler) override
    {
        // TODO: cosine-weighted hemispherical sampling
        const auto sp = sample_to_hemisphere(sampler.get_next_2D());
        return { m_albedo, Vector3{ sp }, 1.0f };
    }

    float pdf_impl(const Vector3&, const Vector3&) const override
    {
        return 1.0f / (2.0f * std::numbers::pi_v<float>);
    }

    RGB m_albedo;
};

} // namespace sp
