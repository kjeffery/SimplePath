#pragma once

/// @author Keith Jeffery

// Based off of Intel Embree's LinearSpace3

#include "Angles.h"
#include "Quaternion.h"
#include "Ray.h"
#include "Vector3.h"

namespace sp {
class LinearSpace3x3
{
public:
    explicit LinearSpace3x3(NoInitType) noexcept
    : m_vx(no_init)
    , m_vy(no_init)
    , m_vz(no_init)
    {
    }

    LinearSpace3x3() noexcept = default;

    LinearSpace3x3(const LinearSpace3x3&) noexcept            = default;
    LinearSpace3x3(LinearSpace3x3&&) noexcept                 = default;
    LinearSpace3x3& operator=(const LinearSpace3x3&) noexcept = default;
    LinearSpace3x3& operator=(LinearSpace3x3&&) noexcept      = default;

    LinearSpace3x3(Vector3 vx, Vector3 vy, Vector3 vz) noexcept
    : m_vx(std::move(vx))
    , m_vy(std::move(vy))
    , m_vz(std::move(vz))
    {
    }

    static LinearSpace3x3 identity() noexcept
    {
        return {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        };
    }

    explicit LinearSpace3x3(const Quaternion& q) noexcept
    : m_vx((q.r * q.r + q.i * q.i - q.j * q.j - q.k * q.k),
           2.0f * (q.i * q.j + q.r * q.k),
           2.0f * (q.i * q.k - q.r * q.j))
    , m_vy(2.0f * (q.i * q.j - q.r * q.k),
           (q.r * q.r - q.i * q.i + q.j * q.j - q.k * q.k),
           2.0f * (q.j * q.k + q.r * q.i))
    , m_vz(2.0f * (q.i * q.k + q.r * q.j),
           2.0f * (q.j * q.k - q.r * q.i),
           (q.r * q.r - q.i * q.i - q.j * q.j + q.k * q.k))
    {
    }

    // Row-major construction
    // clang-format off
    LinearSpace3x3(const float m00, const float m01, const float m02,
                   const float m10, const float m11, const float m12,
                   const float m20, const float m21, const float m22) noexcept
    // clang-format on
    : m_vx(m00, m10, m20)
    , m_vy(m01, m11, m21)
    , m_vz(m02, m12, m22)
    {
    }

    friend bool operator==(const LinearSpace3x3& a, const LinearSpace3x3& b) = default;

    [[nodiscard]] float determinant() const noexcept
    {
        return dot(m_vx, cross(m_vy, m_vz));
    }

    [[nodiscard]] LinearSpace3x3 adjoint() const noexcept
    {
        return LinearSpace3x3{ cross(m_vy, m_vz), cross(m_vz, m_vx), cross(m_vx, m_vy) }.transposed();
    }

    [[nodiscard]] LinearSpace3x3 inverse() const noexcept;

    [[nodiscard]] LinearSpace3x3 transposed() const noexcept
    {
        return LinearSpace3x3{ m_vx.x, m_vx.y, m_vx.z, m_vy.x, m_vy.y, m_vy.z, m_vz.x, m_vz.y, m_vz.z };
    }

    [[nodiscard]] Vector3 row0() const noexcept
    {
        return Vector3{ m_vx.x, m_vy.x, m_vz.x };
    }

    [[nodiscard]] Vector3 row1() const noexcept
    {
        return Vector3{ m_vx.y, m_vy.y, m_vz.y };
    }

    [[nodiscard]] Vector3 row2() const noexcept
    {
        return Vector3{ m_vx.z, m_vy.z, m_vz.z };
    }

    [[nodiscard]] const Vector3& col0() const noexcept
    {
        return m_vx;
    }

    [[nodiscard]] const Vector3& col1() const noexcept
    {
        return m_vy;
    }

    [[nodiscard]] const Vector3& col2() const noexcept
    {
        return m_vz;
    }

    static LinearSpace3x3 scale(const Vector3& s) noexcept
    {
        // clang-format off
        return LinearSpace3x3{
            s.x, 0.0f, 0.0f,
            0.0f, s.y, 0.0f,
            0.0f, 0.0f, s.z
        };
        // clang-format on
    }

    // Return matrix for rotation about arbitrary axis
    static LinearSpace3x3 rotate(Vector3 u, const Angle r) noexcept
    {
        u             = normalize(u);
        const float s = std::sin(r.as_radians());
        const float c = std::cos(r.as_radians());

        // clang-format off
        return LinearSpace3x3{
            u.x * u.x + (1 - u.x * u.x) * c,
            u.x * u.y * (1 - c) - u.z * s,
            u.x * u.z * (1 - c) + u.y * s,
            u.x * u.y * (1 - c) + u.z * s,
            u.y * u.y + (1 - u.y * u.y) * c,
            u.y * u.z * (1 - c) - u.x * s,
            u.x * u.z * (1 - c) - u.y * s,
            u.y * u.z * (1 - c) + u.x * s,
            u.z * u.z + (1 - u.z * u.z) * c
        };
        // clang-format on
    }

    Point3 operator()(const Point3& a) const noexcept
    {
        return Point3{ madd(Vector3{ a.x }, m_vx, madd(Vector3{ a.y }, m_vy, Vector3{ a.z } * m_vz)) };
    }

    Vector3 operator()(const Vector3& a) const noexcept
    {
        return Point3{ madd(Vector3{ a.x }, m_vx, madd(Vector3{ a.y }, m_vy, Vector3{ a.z } * m_vz)) };
    }

    Normal3 operator()(const Normal3& a) const noexcept
    {
        const auto ls = this->inverse().transposed();
        return Normal3{ madd(Vector3{ a.x }, ls.m_vx, madd(Vector3{ a.y }, ls.m_vy, Vector3{ a.z } * ls.m_vz)) };
    }

    Ray operator()(const Ray& r) const noexcept
    {
        const auto& m = *this;
        return Ray{ m(r.get_origin()), m(r.get_direction()) };
    }

private:
    // Column vectors
    Vector3 m_vx;
    Vector3 m_vy;
    Vector3 m_vz;
};

////////////////////////////////////////////////////////////////////////////////
// Unary Operators
////////////////////////////////////////////////////////////////////////////////

inline LinearSpace3x3 operator-(const LinearSpace3x3& a)
{
    return { -a.col0(), -a.col1(), -a.col2() };
}

inline LinearSpace3x3 operator+(const LinearSpace3x3& a)
{
    return { +a.col0(), +a.col1(), +a.col2() };
}

inline LinearSpace3x3 rcp(const LinearSpace3x3& a)
{
    return a.inverse();
}

////////////////////////////////////////////////////////////////////////////////
// Binary Operators
////////////////////////////////////////////////////////////////////////////////

inline LinearSpace3x3 operator+(const LinearSpace3x3& a, const LinearSpace3x3& b)
{
    return LinearSpace3x3{ a.col0() + b.col0(), a.col1() + b.col1(), a.col2() + b.col2() };
}

inline LinearSpace3x3 operator-(const LinearSpace3x3& a, const LinearSpace3x3& b)
{
    return LinearSpace3x3{ a.col0() - b.col0(), a.col1() - b.col1(), a.col2() - b.col2() };
}

inline LinearSpace3x3 operator*(const float a, const LinearSpace3x3& b)
{
    return LinearSpace3x3{ a * b.col0(), a * b.col1(), a * b.col2() };
}

inline Vector3 operator*(const LinearSpace3x3& a, const Vector3& b)
{
    return madd(Vector3{ b.x }, a.col0(), madd(Vector3{ b.y }, a.col1(), Vector3{ b.z } * a.col2()));
}

inline LinearSpace3x3 operator*(const LinearSpace3x3& a, const LinearSpace3x3& b)
{
    return LinearSpace3x3{ a * b.col0(), a * b.col1(), a * b.col2() };
}

inline LinearSpace3x3 operator/(const LinearSpace3x3& a, const float b)
{
    return LinearSpace3x3{ a.col0() / b, a.col1() / b, a.col2() / b };
}

inline LinearSpace3x3 operator/(const LinearSpace3x3& a, const LinearSpace3x3& b)
{
    return a * rcp(b);
}

inline LinearSpace3x3& operator*=(LinearSpace3x3& a, const LinearSpace3x3& b)
{
    return a = a * b;
}

inline LinearSpace3x3& operator/=(LinearSpace3x3& a, const LinearSpace3x3& b)
{
    return a = a / b;
}

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////

inline bool compare(const LinearSpace3x3& a, const LinearSpace3x3& b) noexcept
{
    return compare(a.col0(), b.col0()) && compare(a.col1(), b.col1()) && compare(a.col2(), b.col2());
}

////////////////////////////////////////////////////////////////////////////////
/// Blending
////////////////////////////////////////////////////////////////////////////////

inline LinearSpace3x3 lerp(const LinearSpace3x3& l0, const LinearSpace3x3& l1, const float t)
{
    return { lerp(l0.col0(), l1.col0(), t), lerp(l0.col1(), l1.col1(), t), lerp(l0.col2(), l1.col2(), t) };
}

////////////////////////////////////////////////////////////////////////////////
/// Output Operators
////////////////////////////////////////////////////////////////////////////////

inline std::ostream& operator<<(std::ostream& outs, const LinearSpace3x3& m)
{
    return outs << "{ vx = " << m.col0() << ", vy = " << m.col1() << ", vz = " << m.col2() << "}";
}

inline LinearSpace3x3 LinearSpace3x3::inverse() const noexcept
{
    return adjoint() / determinant();
}
} // namespace sp
