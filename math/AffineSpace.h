#pragma once

/// @author Keith Jeffery

// Based on Intel Embree's AffineSpaceT

#include "BBox.h"
#include "LinearSpace3x3.h"
#include "Vector3.h"

namespace sp {
class AffineSpace
{
    LinearSpace3x3 m_linear;
    Vector3        m_affine;

public:
    explicit AffineSpace(NoInitType) noexcept
    : m_linear(no_init)
    , m_affine(no_init)
    {
    }

    AffineSpace(const AffineSpace&)            = default;
    AffineSpace(AffineSpace&&)                 = default;
    AffineSpace& operator=(const AffineSpace&) = default;
    AffineSpace& operator=(AffineSpace&&)      = default;

    AffineSpace(Vector3 vx, Vector3 vy, Vector3 vz, Vector3 p) noexcept
    : m_linear(std::move(vx), std::move(vy), std::move(vz))
    , m_affine(std::move(p))
    {
    }

    AffineSpace(LinearSpace3x3 ls, Vector3 p) noexcept
    : m_linear(std::move(ls))
    , m_affine(std::move(p))
    {
    }

    explicit AffineSpace(LinearSpace3x3 ls) noexcept
    : m_linear(std::move(ls))
    , m_affine(Vector3::zero())
    {
    }

    static AffineSpace identity() noexcept
    {
        return AffineSpace{ LinearSpace3x3::identity(), Vector3::zero() };
    }

    inline static AffineSpace scale(const Vector3& s) noexcept;
    inline static AffineSpace translate(const Vector3& p) noexcept;
    // Rotation about an arbitrary axis
    inline static AffineSpace rotate(const Vector3& u, float r) noexcept;
    // Rotation about an arbitrary axis and point
    inline static AffineSpace rotate(const Point3& p, const Vector3& u, float r) noexcept;

    static AffineSpace look_at(const Point3& eye, const Point3& point, const Vector3& up) noexcept
    {
        const auto z = normalize(point - eye);
        const auto u = normalize(cross(up, z));
        const auto v = normalize(cross(z, u));
        return AffineSpace{ LinearSpace3x3{ u, v, z }, eye };
    }

    inline static AffineSpace orthographic(float z_near, float z_far) noexcept;

    [[nodiscard]] const LinearSpace3x3& get_linear() const noexcept
    {
        return m_linear;
    }

    [[nodiscard]] const Vector3& get_affine() const noexcept
    {
        return m_affine;
    }

    inline Point3 operator()(const Point3& p) const noexcept
    {
        const auto& m = *this;
        return madd(
            Vector3(p.x),
            m.get_linear().col0(),
            madd(Vector3(p.y), m.get_linear().col1(), madd(Vector3(p.z), m.get_linear().col2(), m.get_affine())));
    }

    inline Vector3 operator()(const Vector3& v) const noexcept
    {
        return get_linear()(v);
    }

    inline Normal3 operator()(const Normal3& n) const noexcept
    {
        return get_linear()(n);
    }

    Ray operator()(const Ray& r) const noexcept
    {
        const auto& m = *this;
        return Ray{ m(r.get_origin()), m(r.get_direction()) };
    }

    inline BBox3 operator()(const BBox3& b) const noexcept
    {
        const auto& m = *this;

        BBox3        dst;
        const Point3 p0(b.get_lower().x, b.get_lower().y, b.get_lower().z);
        dst.extend(m(p0));
        const Point3 p1(b.get_lower().x, b.get_lower().y, b.get_upper().z);
        dst.extend(m(p1));
        const Point3 p2(b.get_lower().x, b.get_upper().y, b.get_lower().z);
        dst.extend(m(p2));
        const Point3 p3(b.get_lower().x, b.get_upper().y, b.get_upper().z);
        dst.extend(m(p3));
        const Point3 p4(b.get_upper().x, b.get_lower().y, b.get_lower().z);
        dst.extend(m(p4));
        const Point3 p5(b.get_upper().x, b.get_lower().y, b.get_upper().z);
        dst.extend(m(p5));
        const Point3 p6(b.get_upper().x, b.get_upper().y, b.get_lower().z);
        dst.extend(m(p6));
        const Point3 p7(b.get_upper().x, b.get_upper().y, b.get_upper().z);
        dst.extend(m(p7));
        return dst;
    }

    friend bool operator==(const AffineSpace&, const AffineSpace&) = default;
};

////////////////////////////////////////////////////////////////////////////////
// Unary Operators
////////////////////////////////////////////////////////////////////////////////

inline AffineSpace operator-(const AffineSpace& a) noexcept
{
    return AffineSpace{ -a.get_linear(), -a.get_affine() };
}

inline AffineSpace operator+(const AffineSpace& a) noexcept
{
    return AffineSpace{ +a.get_linear(), +a.get_affine() };
}

inline AffineSpace rcp(const AffineSpace& a) noexcept
{
    const auto il = rcp(a.get_linear());
    return AffineSpace{ il, -(il * a.get_affine()) };
}

////////////////////////////////////////////////////////////////////////////////
// Binary Operators
////////////////////////////////////////////////////////////////////////////////

inline AffineSpace operator+(const AffineSpace& a, const AffineSpace& b) noexcept
{
    return AffineSpace{ a.get_linear() + b.get_linear(), a.get_affine() + b.get_affine() };
}

inline AffineSpace operator-(const AffineSpace& a, const AffineSpace& b) noexcept
{
    return AffineSpace{ a.get_linear() - b.get_linear(), a.get_affine() - b.get_affine() };
}

inline AffineSpace operator*(const float a, const AffineSpace& b) noexcept
{
    return AffineSpace{ a * b.get_linear(), a * b.get_affine() };
}

inline AffineSpace operator*(const AffineSpace& a, const AffineSpace& b) noexcept
{
    return AffineSpace{ a.get_linear() * b.get_linear(), a.get_linear() * b.get_affine() + a.get_affine() };
}

inline AffineSpace operator*(const AffineSpace& a, const LinearSpace3x3& b) noexcept
{
    // Since the linear space has a zero translation component, we can avoid a multiplication by zero.
    //return AffineSpace{ a.get_linear() * b, a.get_linear() * Vector3{ 0 } + a.get_affine() };
    return AffineSpace{ a.get_linear() * b, a.get_affine() };
}

inline AffineSpace operator*(const LinearSpace3x3& a, const AffineSpace& b) noexcept
{
    // Since the linear space has a zero translation component, we can avoid an addition by zero.
    //return AffineSpace{ a * b.get_linear(), a * b.get_affine() + Vector3{ 0 } };
    return AffineSpace{ a * b.get_linear(), a * b.get_affine() };
}

inline AffineSpace operator/(const AffineSpace& a, const AffineSpace& b) noexcept
{
    return a * rcp(b);
}

#if 0
inline AffineSpace operator/(const AffineSpace& a, const float b) noexcept
{
    return a * rcp(b);
}
#endif

inline AffineSpace& operator*=(AffineSpace& a, const AffineSpace& b) noexcept
{
    return a = a * b;
}

inline AffineSpace& operator*=(AffineSpace& a, const LinearSpace3x3& b) noexcept
{
    return a = a * b;
}

#if 0
inline AffineSpace& operator*=(AffineSpace& a, const float b) noexcept
{
    return a = a * b;
}
#endif

inline AffineSpace& operator/=(AffineSpace& a, const AffineSpace& b) noexcept
{
    return a = a / b;
}

#if 0
inline AffineSpace& operator/=(AffineSpace& a, const float b) noexcept
{
    return a = a / b;
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////

inline bool compare(const AffineSpace& a, const AffineSpace& b) noexcept
{
    return compare(a.get_linear(), b.get_linear()) && compare(a.get_affine(), b.get_affine());
}

////////////////////////////////////////////////////////////////////////////////
// Creation functions
////////////////////////////////////////////////////////////////////////////////

AffineSpace AffineSpace::scale(const Vector3& s) noexcept
{
    return AffineSpace{ LinearSpace3x3::scale(s) };
}

AffineSpace AffineSpace::translate(const Vector3& p) noexcept
{
    return AffineSpace{ LinearSpace3x3::identity(), p };
}

// Rotation about an arbitrary axis
AffineSpace AffineSpace::rotate(const Vector3& u, const float r) noexcept
{
    return AffineSpace{ LinearSpace3x3::rotate(u, r) };
}

// Rotation about an arbitrary axis and point
AffineSpace AffineSpace::rotate(const Point3& p, const Vector3& u, const float r) noexcept
{
    return translate(Vector3{ +p }) * rotate(u, r) * translate(Vector3{ -p });
}

AffineSpace AffineSpace::orthographic(const float z_near, const float z_far) noexcept
{
    return scale(Vector3{ 1.0f, 1.0f, 1.0f / (z_far - z_near) }) * translate(Vector3{ 0.0f, 0.0f, z_near });
}

////////////////////////////////////////////////////////////////////////////////
// Output Operators
////////////////////////////////////////////////////////////////////////////////

inline std::ostream& operator<<(std::ostream& outs, const AffineSpace& m)
{
    return outs << "{ l = " << m.get_linear() << ", get_affine() = " << m.get_affine() << " }";
}
} // namespace sp
