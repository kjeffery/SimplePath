#pragma once

/// @author Keith Jeffery

#include "Vector2.h"
#include "Vector3.h"

#include <cmath>
#include <numbers>

namespace sp {
inline Vector3 sample_to_hemisphere(const Point2& u)
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

inline Vector3 spherical_direction(const float sin_theta, const float cos_theta, const float phi) noexcept
{
    return Vector3{ sin_theta * std::cos(phi), cos_theta, sin_theta * std::sin(phi) };
}

} // namespace sp