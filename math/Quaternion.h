
#pragma once

/// @author Keith Jeffery

// Based on Intel Embree's QuaternionT

#include "Math.h"
#include "Vector3.h"

#include <cmath>

namespace sp {

class Quaternion
{
public:
    explicit Quaternion(NoInitType) noexcept
    {
    }

    Quaternion() noexcept
    : r{}
    , i{}
    , j{}
    , k{}
    {
    }

    Quaternion(const Quaternion&) noexcept            = default;
    Quaternion(Quaternion&&) noexcept                 = default;
    Quaternion& operator=(const Quaternion&) noexcept = default;
    Quaternion& operator=(Quaternion&&) noexcept      = default;

    explicit Quaternion(const float ir) noexcept
    : r{ ir }
    , i{ 0.0f }
    , j{ 0.0f }
    , k{ 0.0f }
    {
    }

    explicit Quaternion(const Vector3& v) noexcept
    : r{ 0.0f }
    , i{ v.x }
    , j{ v.y }
    , k{ v.z }
    {
    }

    Quaternion(const float ir, const Vector3& v) noexcept
    : r{ ir }
    , i{ v.x }
    , j{ v.y }
    , k{ v.z }
    {
    }

    Quaternion(const float ir, const float ii, const float ij, const float ik) noexcept
    : r{ ir }
    , i{ ii }
    , j{ ij }
    , k{ ik }
    {
    }

    Quaternion(const Vector3& vx, const Vector3& vy, const Vector3& vz) noexcept;

    static Quaternion yaw_pitch_roll(const float yaw, const float pitch, const float roll) noexcept;

    // Return a quaternion for rotation about an arbitrary axis
    static Quaternion rotate(const Vector3& u, const float ir) noexcept
    {
        return Quaternion{ std::cos(0.5f * ir), std::sin(0.5f * ir) * normalize(u) };
    }

    // Return the rotation axis of the quaternion
    Vector3 v() const noexcept
    {
        return Vector3{ i, j, k };
    }

    Point3  operator()(const Point3& v) const noexcept;
    Vector3 operator()(const Vector3& v) const noexcept;
    Normal3 operator()(const Normal3& v) const noexcept;

    friend bool operator==(const Quaternion& a, const Quaternion& b) = default;

    float r;
    float i;
    float j;
    float k;
};

inline Quaternion operator*(const float a, const Quaternion& b) noexcept
{
    return Quaternion{ a * b.r, a * b.i, a * b.j, a * b.k };
}

inline Quaternion operator*(const Quaternion& a, const float b) noexcept
{
    return Quaternion{ a.r * b, a.i * b, a.j * b, a.k * b };
}

////////////////////////////////////////////////////////////////
// Unary Operators
////////////////////////////////////////////////////////////////

inline Quaternion operator+(const Quaternion& a)
{
    return Quaternion{ +a.r, +a.i, +a.j, +a.k };
}

inline Quaternion operator-(const Quaternion& a)
{
    return Quaternion{ -a.r, -a.i, -a.j, -a.k };
}

inline Quaternion conj(const Quaternion& a)
{
    return Quaternion{ a.r, -a.i, -a.j, -a.k };
}

inline float abs(const Quaternion& a)
{
    return std::sqrt(a.r * a.r + a.i * a.i + a.j * a.j + a.k * a.k);
}

inline Quaternion rcp(const Quaternion& a)
{
    return conj(a) * rcp(a.r * a.r + a.i * a.i + a.j * a.j + a.k * a.k);
}

inline Quaternion normalize(const Quaternion& a)
{
    return a * rsqrt(a.r * a.r + a.i * a.i + a.j * a.j + a.k * a.k);
}

// evaluates a*q-r
inline Quaternion msub(const float a, const Quaternion& q, const Quaternion& p)
{
    return Quaternion{ msub(a, q.r, p.r), msub(a, q.i, p.i), msub(a, q.j, p.j), msub(a, q.k, p.k) };
}

// evaluates a*q-r
inline Quaternion madd(const float a, const Quaternion& q, const Quaternion& p)
{
    return Quaternion{ madd(a, q.r, p.r), madd(a, q.i, p.i), madd(a, q.j, p.j), madd(a, q.k, p.k) };
}

////////////////////////////////////////////////////////////////
// Binary Operators
////////////////////////////////////////////////////////////////

inline Quaternion operator+(const float a, const Quaternion& b)
{
    return Quaternion{ a + b.r, b.i, b.j, b.k };
}

inline Quaternion operator+(const Quaternion& a, const float b)
{
    return Quaternion{ a.r + b, a.i, a.j, a.k };
}

inline Quaternion operator+(const Quaternion& a, const Quaternion& b)
{
    return Quaternion{ a.r + b.r, a.i + b.i, a.j + b.j, a.k + b.k };
}

inline Quaternion operator-(const float a, const Quaternion& b)
{
    return Quaternion{ a - b.r, -b.i, -b.j, -b.k };
}

inline Quaternion operator-(const Quaternion& a, const float b)
{
    return Quaternion{ a.r - b, a.i, a.j, a.k };
}

inline Quaternion operator-(const Quaternion& a, const Quaternion& b)
{
    return Quaternion{ a.r - b.r, a.i - b.i, a.j - b.j, a.k - b.k };
}

inline Quaternion operator*(const Quaternion& a, const Quaternion& b)
{
    return Quaternion{ a.r * b.r - a.i * b.i - a.j * b.j - a.k * b.k,
                       a.r * b.i + a.i * b.r + a.j * b.k - a.k * b.j,
                       a.r * b.j - a.i * b.k + a.j * b.r + a.k * b.i,
                       a.r * b.k + a.i * b.j - a.j * b.i + a.k * b.r };
}

inline Vector3 operator*(const Quaternion& a, const Vector3& b)
{
    return Quaternion{ a * Quaternion(b) * conj(a) }.v();
}

inline Quaternion operator/(const float a, const Quaternion& b)
{
    return a * rcp(b);
}

inline Quaternion operator/(const Quaternion& a, const float b)
{
    return a * rcp(b);
}

inline Quaternion operator/(const Quaternion& a, const Quaternion& b)
{
    return a * rcp(b);
}

inline Quaternion& operator+=(Quaternion& a, const float b)
{
    return a = a + b;
}

inline Quaternion& operator+=(Quaternion& a, const Quaternion& b)
{
    return a = a + b;
}

inline Quaternion& operator-=(Quaternion& a, const float b)
{
    return a = a - b;
}

inline Quaternion& operator-=(Quaternion& a, const Quaternion& b)
{
    return a = a - b;
}

inline Quaternion& operator*=(Quaternion& a, const float b)
{
    return a = a * b;
}

inline Quaternion& operator*=(Quaternion& a, const Quaternion& b)
{
    return a = a * b;
}

inline Quaternion& operator/=(Quaternion& a, const float b)
{
    return a = a * rcp(b);
}

inline Quaternion& operator/=(Quaternion& a, const Quaternion& b)
{
    return a = a * rcp(b);
}

inline Point3 Quaternion::operator()(const Point3& v) const noexcept
{
    return *this * Vector3{ v };
}

inline Vector3 Quaternion::operator()(const Vector3& v) const noexcept
{
    return *this * v;
}

inline Normal3 Quaternion::operator()(const Normal3& v) const noexcept
{
    return *this * Vector3{ v };
}

inline float dot(const Quaternion& a, const Quaternion& b)
{
    return a.r * b.r + a.i * b.i + a.j * b.j + a.k * b.k;
}

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////

//inline bool operator==(const Quaternion& a, const Quaternion& b) = default;

////////////////////////////////////////////////////////////////////////////////
/// Orientation Functions
////////////////////////////////////////////////////////////////////////////////

inline Quaternion::Quaternion(const Vector3& vx, const Vector3& vy, const Vector3& vz) noexcept
{
    if (vx.x + vy.y + vz.z >= 0.0f) {
        const float t = 1.0f + (vx.x + vy.y + vz.z);
        const float s = rsqrt(t) * 0.5f;
        r             = t * s;
        i             = (vy.z - vz.y) * s;
        j             = (vz.x - vx.z) * s;
        k             = (vx.y - vy.x) * s;
    } else if (vx.x >= std::max(vy.y, vz.z)) {
        const float t = (1.0f + vx.x) - (vy.y + vz.z);
        const float s = rsqrt(t) * 0.5f;
        r             = (vy.z - vz.y) * s;
        i             = t * s;
        j             = (vx.y + vy.x) * s;
        k             = (vz.x + vx.z) * s;
    } else if (vy.y >= vz.z) { // if ( vy.y >= max(vz.z, vx.x) )
        const float t = (1.0f + vy.y) - (vz.z + vx.x);
        const float s = rsqrt(t) * 0.5f;
        r             = (vz.x - vx.z) * s;
        i             = (vx.y + vy.x) * s;
        j             = t * s;
        k             = (vy.z + vz.y) * s;
    } else { // if ( vz.z >= max(vy.y, vx.x) )
        const float t = (1.0f + vz.z) - (vx.x + vy.y);
        const float s = rsqrt(t) * 0.5f;
        r             = (vx.y - vy.x) * s;
        i             = (vz.x + vx.z) * s;
        j             = (vy.z + vz.y) * s;
        k             = t * s;
    }
}

inline Quaternion Quaternion::yaw_pitch_roll(const float yaw, const float pitch, const float roll) noexcept
{
    const float cya = std::cos(yaw * 0.5f);
    const float cpi = std::cos(pitch * 0.5f);
    const float cro = std::cos(roll * 0.5f);
    const float sya = std::sin(yaw * 0.5f);
    const float spi = std::sin(pitch * 0.5f);
    const float sro = std::sin(roll * 0.5f);
    const float r   = cro * cya * cpi + sro * sya * spi;
    const float i   = cro * cya * spi + sro * sya * cpi;
    const float j   = cro * sya * cpi - sro * cya * spi;
    const float k   = sro * cya * cpi - cro * sya * spi;
    return Quaternion{ r, i, j, k };
}

//////////////////////////////////////////////////////////////////////////////
/// Output Operators
//////////////////////////////////////////////////////////////////////////////

inline std::ostream& operator<<(std::ostream& outs, const Quaternion& q)
{
    return outs << "{ r = " << q.r << ", i = " << q.i << ", j = " << q.j << ", k = " << q.k << " }";
}

//////////////////////////////////////////////////////////////////////////////
/// Interpolation
//////////////////////////////////////////////////////////////////////////////
inline Quaternion lerp(const Quaternion& q0, const Quaternion& q1, const float factor)
{
    const float r = std::lerp(q0.r, q1.r, factor);
    const float i = std::lerp(q0.i, q1.i, factor);
    const float j = std::lerp(q0.j, q1.j, factor);
    const float k = std::lerp(q0.k, q1.k, factor);
    return Quaternion{ r, i, j, k };
}

} // namespace sp
