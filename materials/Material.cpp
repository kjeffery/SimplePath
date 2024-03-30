/// @author Keith Jeffery

#include "Material.h"

#include "../math/Sampler.h"
#include "../math/Sampling.h"
#include "../math/Vector3.h"

#include <cassert>
#include <numbers>

namespace sp {
namespace {
auto beckmann_sample11(float cos_theta_i, float U1, float U2)
{
    // Special case (normal incidence)
    if (cos_theta_i > .9999f) {
        const float r         = std::sqrt(-std::log(1.0f - U1));
        const float sin_phi_v = std::sin(2.0f * std::numbers::pi_v<float> * U2);
        const float cos_phi_v = std::cos(2.0f * std::numbers::pi_v<float> * U2);
        const float slope_x   = r * cos_phi_v;
        const float slope_y   = r * sin_phi_v;
        return std::make_pair(slope_x, slope_y);
    }

    // The original inversion routine from the paper contained discontinuities, which causes issues for QMC integration
    // and techniques like Kelemen-style MLT. The following code performs a numerical inversion with better behavior.
    const float sin_theta_i = std::sqrt(std::max(0.0f, 1.0f - square(cos_theta_i)));
    const float tan_theta_i = sin_theta_i / cos_theta_i;
    const float cot_theta_i = 1.0f / tan_theta_i;

    // Search interval -- everything is parameterized in the Erf() domain
    float       a        = -1.0f;
    float       c        = std::erf(cot_theta_i);
    const float sample_x = std::max(U1, 1e-6f);

    // Start with a good initial guess
    // Float b = (1-sample_x) * a + sample_x * c;

    // We can do better (inverse of an approximation computed in Mathematica)
    const float theta_i = std::acos(cos_theta_i);
    const float fit     = 1.0f + theta_i * (-0.876f + theta_i * (0.4265f - 0.0594f * theta_i));
    float       b       = c - (1.0f + c) * std::pow(1.0f - sample_x, fit);

    // Normalization factor for the CDF
    static const float sqrt_pi_inv   = 1.0f / std::sqrt(std::numbers::pi_v<float>);
    const float        normalization = 1.0f / (1.0f + c + sqrt_pi_inv * tan_theta_i * std::exp(-cot_theta_i * cot_theta_i));

    for (int it = 0; it < 9; ++it) {
        // Bisection criterion -- the oddly-looking Boolean expression are intentional to check for NaNs at little
        // additional cost
        if (!(b >= a && b <= c)) {
            b = 0.5f * (a + c);
        }

        // Evaluate the CDF and its derivative (i.e. the density function)
        const float inv_erf = erfinv(b);
        const float value   =
                normalization * (1.0f + b + sqrt_pi_inv * tan_theta_i * std::exp(-inv_erf * inv_erf)) - sample_x;
        const float derivative = normalization * (1.0f - inv_erf * tan_theta_i);

        if (std::abs(value) < 1e-5f) {
            break;
        }

        // Update bisection intervals
        if (value > 0) {
            c = b;
        } else {
            a = b;
        }

        b -= value / derivative;
    }

    /* Now convert back into a slope value */
    const float slope_x = erfinv(b);

    /* Simulate Y component */
    const float slope_y = erfinv(2.0f * std::max(U2, 1e-6f) - 1.0f);

    assert(!std::isinf(slope_x));
    assert(!std::isnan(slope_x));
    assert(!std::isinf(slope_y));
    assert(!std::isnan(slope_y));
    return std::make_pair(slope_x, slope_y);
}

Vector3 beckmann_sample(const Vector3& wi, float alpha_x, float alpha_y, float U1, float U2)
{
    // 1. stretch wi
    const Vector3 wi_stretched = normalize(Vector3(alpha_x * wi.x, wi.y, alpha_y * wi.z));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    auto [slope_x, slope_y] = beckmann_sample11(cos_theta(wi_stretched), U1, U2);

    // 3. rotate
    const float tmp = cos_phi(wi_stretched) * slope_x - sin_phi(wi_stretched) * slope_y;
    slope_y         = sin_phi(wi_stretched) * slope_x + cos_phi(wi_stretched) * slope_y;
    slope_x         = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    return normalize(Vector3(-slope_x, 1.0f, -slope_y));
}
} // namespace

Vector3 BeckmannDistribution::sample_wh_impl(const Vector3& wo, Sampler& sampler) const
{
    if (!sample_visible_area()) {
        float tan2_theta;
        float phi;
        if (m_alpha_x == m_alpha_y) {
            const float log_sample = std::log(1.0f - sampler.get_next_1D());
            assert(!std::isinf(log_sample));
            tan2_theta = -m_alpha_x * m_alpha_x * log_sample;
            phi        = sampler.get_next_1D() * 2.0f * std::numbers::phi_v<float>;
        } else {
            // Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
            // distribution
            const float log_sample = std::log(1.0f - sampler.get_next_1D());
            assert(!std::isinf(log_sample));
            const float u1 = sampler.get_next_1D();
            phi            = std::atan(
                m_alpha_y / m_alpha_x *
                std::tan(2.0f * std::numbers::pi_v<float> * u1 + 0.5f * std::numbers::pi_v<float>));
            if (u1 > 0.5f) {
                phi += std::numbers::pi_v<float>;
            }
            const float sin_phi  = std::sin(phi);
            const float cos_phi  = std::cos(phi);
            const float alpha_x2 = square(m_alpha_x);
            const float alpha_y2 = square(m_alpha_y);
            tan2_theta           = -log_sample / (square(cos_phi) / alpha_x2 + square(sin_phi) / alpha_y2);
        }

        // Map sampled Beckmann angles to normal direction _wh_
        const float cos_theta_v = 1.0f / std::sqrt(1.0f + tan2_theta);
        const float sin_theta_v = std::sqrt(std::max(0.0f, 1.0f - square(cos_theta_v)));
        Vector3     wh          = spherical_direction(sin_theta_v, cos_theta_v, phi);
        if (!same_hemisphere(wo, wh)) {
            wh = -wh;
        }
        return wh;
    } else {
        // Sample visible area of normals for Beckmann distribution
        const auto flip = wo.y < 0.0f;
        auto       wh   = beckmann_sample(flip ? -wo : wo, m_alpha_x, m_alpha_y, sampler.get_next_1D(), sampler.get_next_1D());
        if (flip) {
            wh = -wh;
        }
        return wh;
    }
}
} // namespace sp
