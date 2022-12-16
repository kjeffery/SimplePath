#pragma once

/// @author Keith Jeffery

#include "Vector2.h"

namespace sp {
namespace detail {
float mod1(const float f) noexcept
{
    return f - std::floor(f);
}
} // namespace detail

Point2 r2_sequence(unsigned n, const Point2& seed = Point2{ 0.5f, 0.5f })
{
    constexpr float g  = 1.32471795724474602596f;
    constexpr float a0 = 1.0f / g;
    constexpr float a1 = 1.0f / (g * g);
    return Point2{ detail::mod1(seed[0] + a0 * n), detail::mod1(seed[1] + a1 * n) };
}

} // namespace sp
