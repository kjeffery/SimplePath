//
// Created by krjef on 11/3/2022.
//

#include "Vector2.h"

#include <numbers>

namespace sp {
Point2 sample_to_concentric_disk(const Point2& u) noexcept
{
    constexpr float pi_over_4 = std::numbers::pi_v<float> / 4.0f;
    constexpr float pi_over_2 = std::numbers::pi_v<float> / 2.0f;

    // Map uniform random numbers to $[-1,1]^2$
    Point2 u_offset = 2.0f * u - Vector2{1.0f, 1.0f};

    // Handle degeneracy at the origin
    if (u_offset.x == 0.0f && u_offset.y == 0.0f) {
        return Point2{0.0f, 0.0f};
    }

    // Apply concentric mapping to point
    float theta;
    float r;
    if (std::abs(u_offset.x) > std::abs(u_offset.y)) {
        r     = u_offset.x;
        theta = pi_over_4 * (u_offset.y / u_offset.x);
    } else {
        r     = u_offset.y;
        theta = pi_over_2 - pi_over_4 * (u_offset.x / u_offset.y);
    }
    return r * Point2(std::cos(theta), std::sin(theta));
}
} // namespace sp