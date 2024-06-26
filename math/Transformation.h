#pragma once

#include "AffineSpace.h"
#include "BBox.h"
#include "LinearSpace3x3.h"
#include "Ray.h"
#include "Vector3.h"

#include <cassert>
#include <utility>

#include "Angles.h"

namespace sp {
template <typename T, typename U>
struct Promote;

template <typename T>
struct Promote<T, T>
{
    using type = T;
};

template <>
struct Promote<AffineSpace, LinearSpace3x3>
{
    using type = AffineSpace;
};

template <>
struct Promote<LinearSpace3x3, AffineSpace>
{
    using type = AffineSpace;
};

template <typename T>
class Transformation
{
public:
    template <typename U>
    friend class Transformation;

    Transformation(T transform, T inverse)
    : m_transform{ std::move(transform) }
    , m_inverse{ std::move(inverse) }
    {
        assert(compare(m_transform * m_inverse, T::identity()));
    }

    static Transformation identity() noexcept
    {
        return Transformation{ T::identity(), T::identity() };
    }

    static Transformation compute_inverse(T transform)
    {
        return Transformation{ transform, rcp(transform) };
    }

    Point3 operator()(const Point3& p) const noexcept
    {
        return m_transform(p);
    }

    Vector3 operator()(const Vector3& v) const noexcept
    {
        return m_transform(v);
    }

    Normal3 operator()(const Normal3& n) const noexcept
    {
        return m_transform(n);
    }

    Ray operator()(const Ray& r) const noexcept
    {
        return m_transform(r);
    }

    BBox3 operator()(const BBox3& b) const noexcept
    {
        return m_transform(b);
    }

    const T& get_transform() const noexcept
    {
        return m_transform;
    }

    const T& get_inverse() const noexcept
    {
        return m_inverse;
    }

    template <typename U>
    Transformation& operator*=(const Transformation<U>& b)
    {
        m_transform *= b.m_transform;
        m_inverse = b.m_inverse * m_inverse;
        assert(compare(m_transform * m_inverse, T::identity()));
        return *this;
    }

private:
    T m_transform;
    T m_inverse;
};

template <typename T, typename U>
auto operator*(Transformation<T> a, const Transformation<U>& b)
{
    using Type = typename Promote<T, U>::type;
    Type result{ a };
    result *= b;
    assert(compare(a.m_transform * a.m_inverse, T::identity()));
    return result;
}

using AffineTransformation = Transformation<AffineSpace>;
using LinearTransformation = Transformation<LinearSpace3x3>;

inline LinearTransformation scale(const Vector3& s) noexcept
{
    return { LinearSpace3x3::scale(s), LinearSpace3x3::scale(1.0f / s) };
}

// Rotation about an arbitrary axis
inline LinearTransformation rotate(const Vector3& u, const Angle r) noexcept
{
    return { LinearSpace3x3::rotate(u, r), LinearSpace3x3::rotate(u, -r) };
}

// Rotation about an arbitrary axis and point
inline AffineTransformation rotate(const Point3& p, const Vector3& u, const Angle r) noexcept
{
    return { AffineSpace::rotate(p, u, r), AffineSpace::rotate(p, u, -r) };
}

inline AffineTransformation translate(const Vector3& p) noexcept
{
    return { AffineSpace::translate(p), AffineSpace::translate(-p) };
}
} // namespace sp
