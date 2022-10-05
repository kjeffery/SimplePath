#pragma once

/// @author Keith Jeffery

#include <cassert>
#include <cstdint>
#include <ostream>

namespace sp {

struct RGB
{
    constexpr RGB() noexcept
    : r(0)
    , g(0)
    , b(0)
    {
    }

    constexpr RGB(float v) noexcept
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
};

inline RGB& operator+=(RGB& a, const RGB& b) noexcept
{
    a.r += b.r;
    a.g += b.g;
    a.b += b.b;
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
inline RGB operator/(RGB a, float b) noexcept
{
    return a /= b;
}

inline std::ostream& operator<<(std::ostream& outs, const RGB& c)
{
    return outs << c.r << ' ' << c.g << ' ' << c.b;
}

} // namespace sp
