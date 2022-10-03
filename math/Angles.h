#pragma once

/// @author Keith Jeffery

#include <numbers>

namespace sp {
class Radians
{
public:
    explicit Radians(float f) noexcept
    : m_value(f)
    {
    }

    float as_float() const noexcept
    {
        return m_value;
    }

    operator float() const noexcept
    {
        return m_value;
    }

private:
    float m_value;
};

class Degrees
{
public:
    explicit Degrees(float f) noexcept
    : m_value(f)
    {
    }

    float as_float() const noexcept
    {
        return m_value;
    }

    operator float() const noexcept
    {
        return m_value;
    }

private:
    float m_value;
};

inline Radians operator"" _radians(long double f)
{
    return Radians{ static_cast<float>(f) };
}

inline Degrees operator"" _degrees(long double f)
{
    return Degrees{ static_cast<float>(f) };
}

inline Degrees to_degrees(Radians radians) noexcept
{
    return Degrees{ radians * 180.f / std::numbers::pi_v<float> };
}

inline Radians to_radians(Degrees degrees) noexcept
{
    return Radians{ degrees * std::numbers::pi_v<float> / 180.0f };
}
} // namespace sp
