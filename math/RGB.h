#pragma once

/// @author Keith Jeffery

#include "Math.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <istream>
#include <ostream>
#include <utility>

namespace sp {

struct RGB
{
    constexpr RGB() noexcept
    : r(0)
    , g(0)
    , b(0)
    {
    }

    explicit constexpr RGB(float v) noexcept
    : r(v)
    , g(v)
    , b(v)
    {
    }

    constexpr RGB(float ir, float ig, float ib) noexcept
    : r(ir)
    , g(ig)
    , b(ib)
    {
    }

    float& operator[](std::uint32_t idx) noexcept
    {
        switch (idx) {
        case 0u:
            return r;
        case 1u:
            return g;
        case 2u:
            return b;
        default:
            assert(!"Should not get here");
        }
        return r;
    }

    float operator[](std::uint32_t idx) const noexcept
    {
        switch (idx) {
        case 0u:
            return r;
        case 1u:
            return g;
        case 2u:
            return b;
        default:
            assert(!"Should not get here");
        }
        return r;
    }

    static constexpr RGB black() noexcept
    {
        return RGB{ 0.0f };
    }

    // For some definition of white (We're HDR in these here woods).
    static constexpr RGB white() noexcept
    {
        return RGB{ 1.0f };
    }

    float r;
    float g;
    float b;

    friend bool operator==(const RGB&, const RGB&) noexcept = default;
};

inline RGB& operator+=(RGB& a, const RGB& b) noexcept
{
    a.r += b.r;
    a.g += b.g;
    a.b += b.b;
    return a;
}

inline RGB& operator-=(RGB& a, const RGB& b) noexcept
{
    a.r -= b.r;
    a.g -= b.g;
    a.b -= b.b;
    return a;
}

inline RGB& operator*=(RGB& a, const RGB& b) noexcept
{
    a.r *= b.r;
    a.g *= b.g;
    a.b *= b.b;
    return a;
}

inline RGB& operator*=(RGB& a, float b) noexcept
{
    a.r *= b;
    a.g *= b;
    a.b *= b;
    return a;
}

inline RGB& operator/=(RGB& a, const RGB& b) noexcept
{
    a.r /= b.r;
    a.g /= b.g;
    a.b /= b.b;
    return a;
}

inline RGB& operator/=(RGB& a, float b) noexcept
{
    a.r /= b;
    a.g /= b;
    a.b /= b;
    return a;
}

// Pass-by-value on purpose
inline RGB operator+(RGB a, const RGB& b) noexcept
{
    return a += b;
}

// Pass-by-value on purpose
inline RGB operator-(RGB a, const RGB& b) noexcept
{
    return a -= b;
}

// Pass-by-value on purpose
inline RGB operator*(RGB a, const RGB& b) noexcept
{
    return a *= b;
}

// Pass-by-value on purpose
inline RGB operator*(RGB a, float b) noexcept
{
    return a *= b;
}

// Pass-by-value on purpose
inline RGB operator*(float b, RGB a) noexcept
{
    return a *= b;
}

// Pass-by-value on purpose
inline RGB operator/(RGB a, const RGB& b) noexcept
{
    return a /= b;
}

// Pass-by-value on purpose
inline RGB operator/(RGB a, float b) noexcept
{
    return a /= b;
}

inline RGB max(const RGB& a, const RGB& b) noexcept
{
    return { std::max(a.r, b.r), std::max(a.g, b.g), std::max(a.b, b.b) };
}

inline bool compare(const RGB& a, const RGB& b) noexcept
{
    return float_compare(a.r, b.r) && float_compare(a.g, b.g) && float_compare(a.b, b.b);
}

inline bool compare_epsilon(const RGB& a, const RGB& b, const float epsilon) noexcept
{
    return float_compare_epsilon(a.r, b.r, epsilon) && float_compare_epsilon(a.g, b.g, epsilon) &&
           float_compare_epsilon(a.b, b.b, epsilon);
}

template <>
inline RGB safe_divide(const RGB& a, const float& b) noexcept
{
    if (b == 0.0f) {
        return RGB::black();
    } else {
        return a / b;
    }
}

template <>
inline RGB safe_divide(const RGB& a, const RGB& b) noexcept
{
    RGB result = a / b;
    for (int i = 0; i < 3; ++i) {
        if (b[i] == 0.0f) {
            result[i] = 0.0f;
        }
    }
    return result;
}

constexpr float relative_luminance(const RGB& c) noexcept
{
    return 0.2126f * c.r + 0.7152f * c.g + 0.0722f * c.b;
}

inline std::ostream& operator<<(std::ostream& outs, const RGB& c)
{
    return outs << c.r << ' ' << c.g << ' ' << c.b;
}

inline std::istream& operator>>(std::istream& ins, RGB& c)
{
    ins >> c.r;
    ins >> c.g;
    ins >> c.b;
    return ins;
}

} // namespace sp

namespace std
{
inline bool isinf(const sp::RGB& rgb)
{
    return isinf(rgb.r) || isinf(rgb.g) || isinf(rgb.b);
}

inline bool isnan(const sp::RGB& rgb)
{
    return isnan(rgb.r) || isnan(rgb.g) || isnan(rgb.b);
}
}
