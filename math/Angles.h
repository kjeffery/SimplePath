#pragma once

/// @author Keith Jeffery

#include <cmath>
#include <numbers>

namespace sp {

class Radians;
class Degrees;

class Angle
{
public:
    constexpr Angle() noexcept
    : m_radians(0.0f)
    {
    }

    explicit constexpr Angle(const Radians&) noexcept;
    explicit constexpr Angle(const Degrees&) noexcept;

    constexpr Radians as_radians() const noexcept;
    constexpr Degrees as_degrees() const noexcept;

    friend Angle            normalize(const Angle& a) noexcept;
    friend constexpr Angle& operator*=(Angle& a, const float b) noexcept;
    friend constexpr Angle& operator+=(Angle& a, const Angle& b) noexcept;
    friend constexpr Angle& operator-=(Angle& a, const Angle& b) noexcept;

private:
    float m_radians;
};

class Radians
{
public:
    explicit constexpr Radians(float f) noexcept
    : m_value(f)
    {
    }

    constexpr float as_float() const noexcept
    {
        return m_value;
    }

    constexpr operator float() const noexcept
    {
        return m_value;
    }

    friend Radians            normalize(const Radians& a) noexcept;
    friend constexpr Radians& operator*=(Radians& a, const float b) noexcept;
    friend constexpr Radians& operator+=(Radians& a, const Radians& b) noexcept;
    friend constexpr Radians& operator-=(Radians& a, const Radians& b) noexcept;

private:
    float m_value;
};

class Degrees
{
public:
    explicit constexpr Degrees(float f) noexcept
    : m_value(f)
    {
    }

    constexpr float as_float() const noexcept
    {
        return m_value;
    }

    constexpr operator float() const noexcept
    {
        return m_value;
    }

    friend Degrees            normalize(const Degrees& a) noexcept;
    friend constexpr Degrees& operator*=(Degrees& a, const float b) noexcept;
    friend constexpr Degrees& operator+=(Degrees& a, const Degrees& b) noexcept;
    friend constexpr Degrees& operator-=(Degrees& a, const Degrees& b) noexcept;

private:
    float m_value;
};

inline namespace literals {
inline Radians operator"" _radians(long double f)
{
    return Radians{ static_cast<float>(f) };
}

inline Degrees operator"" _degrees(long double f)
{
    return Degrees{ static_cast<float>(f) };
}
} // namespace literals

inline constexpr Degrees to_degrees(Radians radians) noexcept
{
    return Degrees{ radians * 180.f / std::numbers::pi_v<float> };
}

inline constexpr Radians to_radians(Degrees degrees) noexcept
{
    return Radians{ degrees * std::numbers::pi_v<float> / 180.0f };
}

constexpr Angle::Angle(const Radians& a) noexcept
: m_radians(a.as_float())
{
}

constexpr Angle::Angle(const Degrees& a) noexcept
: Angle(to_radians(a))
{
}

constexpr Radians Angle::as_radians() const noexcept
{
    return Radians{ m_radians };
}

constexpr Degrees Angle::as_degrees() const noexcept
{
    return Degrees{ to_degrees(Radians{ m_radians }) };
}

inline Angle normalize(const Angle& a) noexcept
{
    const float value = std::fmod(a.m_radians, 2.0f * std::numbers::pi_v<float>);
    return Angle{ Radians{ value } };
}

inline constexpr Angle& operator*=(Angle& a, const float b) noexcept
{
    a.m_radians *= b;
    return a;
}

inline constexpr Angle& operator+=(Angle& a, const Angle& b) noexcept
{
    a.m_radians += b.m_radians;
    return a;
}

inline constexpr Angle& operator-=(Angle& a, const Angle& b) noexcept
{
    a.m_radians -= b.m_radians;
    return a;
}

inline constexpr Angle operator*(Angle a, const float b) noexcept
{
    a *= b;
    return a;
}

inline constexpr Angle operator*(const float a, Angle b) noexcept
{
    b *= a;
    return b;
}

inline constexpr Angle operator+(Angle a, const Angle& b) noexcept
{
    a += b;
    return a;
}

inline constexpr Angle operator-(Angle a, const Angle& b) noexcept
{
    a -= b;
    return a;
}

inline Radians normalize(const Radians& a) noexcept
{
    const float value = std::fmod(a.m_value, 2.0f * std::numbers::pi_v<float>);
    return Radians{ value };
}

inline constexpr Radians& operator*=(Radians& a, const float b) noexcept
{
    a.m_value *= b;
    return a;
}

inline constexpr Radians& operator+=(Radians& a, const Radians& b) noexcept
{
    a.m_value += b.m_value;
    return a;
}

inline constexpr Radians& operator-=(Radians& a, const Radians& b) noexcept
{
    a.m_value -= b.m_value;
    return a;
}

inline constexpr Radians operator*(Radians a, const float b) noexcept
{
    a *= b;
    return a;
}

inline constexpr Radians operator*(const float a, Radians b) noexcept
{
    b *= a;
    return b;
}

inline constexpr Radians operator+(Radians a, const Radians& b) noexcept
{
    a += b;
    return a;
}

inline constexpr Radians operator-(Radians a, const Radians& b) noexcept
{
    a -= b;
    return a;
}

inline Degrees normalize(const Degrees& a) noexcept
{
    const float value = std::fmod(a.m_value, 360.0f);
    return Degrees{ value };
}

inline constexpr Degrees& operator*=(Degrees& a, const float b) noexcept
{
    a.m_value *= b;
    return a;
}

inline constexpr Degrees& operator+=(Degrees& a, const Degrees& b) noexcept
{
    a.m_value += b.m_value;
    return a;
}

inline constexpr Degrees& operator-=(Degrees& a, const Degrees& b) noexcept
{
    a.m_value -= b.m_value;
    return a;
}

inline constexpr Degrees operator*(Degrees a, const float b) noexcept
{
    a *= b;
    return a;
}

inline constexpr Degrees operator*(const float a, Degrees b) noexcept
{
    b *= a;
    return b;
}

inline constexpr Degrees operator+(Degrees a, const Degrees& b) noexcept
{
    a += b;
    return a;
}

inline constexpr Degrees operator-(Degrees a, const Degrees& b) noexcept
{
    a -= b;
    return a;
}

} // namespace sp
