#pragma once

/// @author Keith Jeffery

#include "RGB.h"

#include "../math/Angles.h"

#include <cassert>
#include <cstdint>
#include <ostream>

namespace sp {

struct HSV
{
    constexpr HSV() noexcept
    : h(0)
    , s(0)
    , v(0)
    {
    }

    constexpr explicit HSV(float v) noexcept
    : h(v)
    , s(v)
    , v(v)
    {
    }

    constexpr HSV(const Degrees ih, const float is, const float iv) noexcept
    : h(ih.as_float() / 360.0f)
    , s(is)
    , v(iv)
    {
    }

    constexpr HSV(const float ih, const float is, const float iv) noexcept
    : h(ih)
    , s(is)
    , v(iv)
    {
        assert(ih <= 1.0f);
    }

    float h;
    float s;
    float v;
};

inline HSV& operator+=(HSV& a, const HSV& b) noexcept
{
    a.h += b.h;
    a.s += b.s;
    a.v += b.v;
    return a;
}

inline HSV& operator*=(HSV& a, float b) noexcept
{
    a.h *= b;
    a.s *= b;
    a.v *= b;
    return a;
}

inline HSV& operator/=(HSV& a, float b) noexcept
{
    a.h /= b;
    a.s /= b;
    a.v /= b;
    return a;
}

// Pass-by-value on purpose
inline HSV operator+(HSV a, const HSV& b) noexcept
{
    return a += b;
}

// Pass-by-value on purpose
inline HSV operator*(HSV a, float b) noexcept
{
    return a *= b;
}

// Pass-by-value on purpose
inline HSV operator*(float b, HSV a) noexcept
{
    return a *= b;
}

// Pass-by-value on purpose
inline HSV operator/(HSV a, float b) noexcept
{
    return a /= b;
}

#if 0
inline RGB to_rgb(const HSV& hsv) noexcept
{
    assert(hsv.h >= 0.0f);
    assert(hsv.h <= 1.0f);
    assert(hsv.s >= 0.0f);
    assert(hsv.s <= 1.0f);
    assert(hsv.v >= 0.0f);
    assert(hsv.v <= 1.0f);

    const auto  i = static_cast<int>(std::floor(hsv.h * 6.0f));
    const float v = hsv.v;
    const float f = hsv.h * 6.0f - i;
    const float p = hsv.v * (1.0f - hsv.s);
    const float q = hsv.v * (1.0f - f * hsv.s);
    const float t = hsv.v * (1.0f - (1.0f - f) * hsv.s);

    switch (i % 6) {
    case 0:
        return RGB{ v, t, p };
    case 1:
        return RGB{ q, v, p };
    case 2:
        return RGB{ p, v, t };
    case 3:
        return RGB{ p, q, v };
    case 4:
        return RGB{ t, p, v };
    case 5:
        return RGB{ v, p, q };
    }
    return RGB::black();
}
#else
inline RGB to_rgb(const HSV& hsv) noexcept
{
    assert(hsv.h >= 0.0f);
    assert(hsv.h <= 1.0f);
    assert(hsv.s >= 0.0f);
    assert(hsv.s <= 1.0f);
    assert(hsv.v >= 0.0f);
    assert(hsv.v <= 1.0f);

    const float C      = hsv.v * hsv.s;
    const auto  hprime = static_cast<int>(std::floor(hsv.h * 6.0f));
    const float X      = C * (1.0f - std::abs(std::fmod(hprime, 2.0f) - 1.0f));
    const float M      = hsv.v - C;
    switch (hprime % 6) {
    case 0:
        return RGB{ C, X, 0 };
    case 1:
        return RGB{ X, C, 0 };
    case 2:
        return RGB{ 0, C, X };
    case 3:
        return RGB{ 0, X, C };
    case 4:
        return RGB{ X, 0, C };
    case 5:
        return RGB{ C, 0, X };
    }
    return RGB::black();
}

#endif

inline std::ostream& operator<<(std::ostream& outs, const HSV& c)
{
    return outs << c.h << ' ' << c.s << ' ' << c.v;
}

} // namespace sp
