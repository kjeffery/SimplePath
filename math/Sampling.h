#pragma once

/// @author Keith Jeffery

#include "Vector2.h"
#include "Vector3.h"

#include <cmath>
#include <numbers>

namespace sp {
inline Vector3 sample_to_hemisphere(const Point2& p) noexcept
{
    const float phi       = 2.0f * std::numbers::pi_v<float> * p.x;
    const float cos_phi   = std::cos(phi);
    const float sin_phi   = std::sin(phi);
    const float cos_theta = p.y;
    const float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
    const float pu        = sin_theta * cos_phi;
    const float pv        = sin_theta * sin_phi;
    const float pw        = cos_theta;
    return Vector3{ pu, pv, pw };
}
} // namespace sp