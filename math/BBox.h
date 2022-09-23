#pragma once

/// @author Keith Jeffery

#include "Point2i.h"
#include "Vector2.h"
#include "Vector3.h"

namespace sp {
template <typename T>
class BBox
{
public:
    BBox() noexcept
    : m_min(T::positive_infinity())
    , m_max(T::negative_infinity())
    {
    }

    BBox(const T& a, const T& b) noexcept
    : m_min(min(a, b))
    , m_max(max(a, b))
    {
    }

    const T& get_lower() const noexcept
    {
        return m_min;
    }

    const T& get_upper() const noexcept
    {
        return m_max;
    }

    void extend(const T& p) noexcept
    {
        m_min = min(p, m_min);
        m_max = max(p, m_max);
    }

    friend bool operator==(const BBox<T>& a, const BBox<T>& b) noexcept = default;

private:
    T m_min;
    T m_max;
};

template <typename T>
inline BBox<T> merge(const BBox<T>& a, const BBox<T>& b) noexcept
{
    return { min(a.get_lower(), b.get_lower()), max(a.get_upper(), b.get_upper()) };
}

template <typename T>
inline BBox<T> intersect(const BBox<T>& a, const BBox<T>& b) noexcept
{
    return { max(a.get_lower(), b.get_lower()), min(a.get_upper(), b.get_upper()) };
}

// Check to see if point is in box, allowing for right-hand equality (closed-closed) [a, b]
template <typename T>
inline bool contains_closed(const BBox<T>& box, const T& p) noexcept
{
    for (std::size_t i = 0; i < T::N; ++i) {
        if (p[i] < box.get_lower()[i]) {
            return false;
        }
        if (p[i] > box.get_upper()[i]) {
            return false;
        }
    }
    return true;
}

// Check to see if point is in box, disallowing right-hand equality (closed-open) [a, b)
template <typename T>
inline bool contains_open(const BBox<T>& box, const T& p) noexcept
{
    for (std::size_t i = 0; i < T::N; ++i) {
        if (p[i] < box.get_lower()[i]) {
            return false;
        }
        if (p[i] >= box.get_upper()[i]) {
            return false;
        }
    }
    return true;
}

template <std::size_t dim, typename T>
typename T::scalar extents(const BBox<T>& box) noexcept
{
    return box.get_upper()[dim] - box.get_lower()[dim];
}

template <typename T>
typename T::scalar extents(const BBox<T>& box, std::size_t dim) noexcept
{
    return box.get_upper()[dim] - box.get_lower()[dim];
}

using BBox2i = BBox<Point2i>;
using BBox2  = BBox<Point2>;
using BBox3  = BBox<Point3>;

} // namespace sp
