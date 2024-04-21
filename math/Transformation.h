#pragma once

#include "BBox.h"
#include "Vector3.h"
#include "Ray.h"

#include <cassert>
#include <utility>

namespace sp {
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

template <typename T>
Transformation<T> operator*(Transformation<T> a, const Transformation<T>& b)
{
    a *= b;
    assert(compare(a.m_transform * a.m_inverse, T::identity()));
    return a;
}

using AffineTransformation = Transformation<AffineSpace>;
using LinearTransformation = Transformation<LinearSpace3x3>;
} // namespace sp
