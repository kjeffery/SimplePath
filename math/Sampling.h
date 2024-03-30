#pragma once

/// @author Keith Jeffery

#include "Vector2.h"
#include "Vector3.h"

#include <cmath>
#include <numbers>

namespace sp {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// These sampling functions operate in local-space, where we assume a right-handed y-up coordinate system. All sampling
// functions should assume that y is the primary axis (e.g. sampling a cosine-weighted hemisphere means that most of the
// samples should cluster around the y-axis with no samples with a negative y-value).
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Point2 sample_to_concentric_disk(const Point2& u) noexcept;

inline Vector3 sample_to_uniform_sphere(const Point2& u) noexcept
{
    const float z   = 1.0f - 2.0f * u.x;
    const float r   = std::sqrt(std::max(0.0f, 1.0f - z * z));
    const float phi = 2.0 * std::numbers::pi_v<float> * u.y;
    return { r * std::cos(phi), r * std::sin(phi), z };
}

constexpr float uniform_sphere_pdf() noexcept
{
    return 1.0f / (4.0f * std::numbers::pi_v<float>);
}

inline Vector3 sample_to_uniform_hemisphere(const Point2& u) noexcept
{
    const float y   = u.x;
    const float r   = std::sqrt(std::max(0.0f, 1.0f - y * y));
    const float phi = 2.0f * std::numbers::pi_v<float> * u.y;
    return { r * std::cos(phi), y, r * std::sin(phi) };
}

inline constexpr float uniform_hemisphere_pdf() noexcept
{
    return 1.0f / (2.0f * std::numbers::pi_v<float>);
}

inline Vector3 sample_to_cosine_hemisphere(const Point2& u) noexcept
{
    const Point2 d = sample_to_concentric_disk(u);
    const float  y = std::sqrt(std::max(0.0f, 1.0f - d.x * d.x - d.y * d.y));
    return Vector3{ d.x, y, d.y };
}

inline float cosine_hemisphere_pdf(const float cos_theta) noexcept
{
    return cos_theta * std::numbers::inv_pi_v<float>;
}

inline Vector3 sample_to_uniform_cone(const Point2& u, const float cos_theta_max) noexcept
{
    const float cos_theta = (1.0f - u[0]) + u[0] * cos_theta_max;
    const float sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);
    const float phi       = u[1] * 2.0f * std::numbers::pi_v<float>;
    return Vector3{ std::cos(phi) * sin_theta, cos_theta, std::sin(phi) * sin_theta };
}

inline float uniform_cone_pdf(const float cos_theta_max) noexcept
{
    return 1.0f / (2.0f * std::numbers::pi_v<float> * (1.0f - cos_theta_max));
}

inline Vector3 spherical_direction(const float sin_theta, const float cos_theta, const float phi) noexcept
{
    return Vector3{ sin_theta * std::cos(phi), cos_theta, sin_theta * std::sin(phi) };
}
} // namespace sp
