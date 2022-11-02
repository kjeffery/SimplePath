#pragma once

///@author Keith Jeffery

// Based on Intel Embree's Vec2fa

#include "Math.h"
#include "Vector2.h"

#include "../base/Constants.h"
#include "../base/Util.h"

#include <algorithm>
#include <cassert>
#include <istream>
#include <limits>
#include <ostream>

#include <immintrin.h>

namespace sp {

struct Point2i
{
    using scalar = int;
    static constexpr int N = 2;

    int x;
    int y;

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    ~Point2i() = default;

    constexpr Point2i() noexcept
    : x(0)
    , y(0)
    {
    }

    explicit Point2i(NoInitType) noexcept
    {
    }

    Point2i(const Point2i&) noexcept            = default;
    Point2i& operator=(const Point2i&) noexcept = default;

    explicit Point2i(const int a) noexcept
    : x(a)
    , y(a)
    {
    }

    Point2i(const int ix, const int iy) noexcept
    : x(ix)
    , y(iy)
    {
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    static Point2i zero() noexcept
    {
        return Point2i{};
    }

    static Point2i one() noexcept
    {
        return Point2i{ 1 };
    }

    static Point2i positive_infinity() noexcept
    {
        return Point2i{ std::numeric_limits<int>::max() };
    }

    static Point2i negative_infinity() noexcept
    {
        return Point2i{ std::numeric_limits<int>::lowest() };
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    const int& operator[](const size_t index) const noexcept
    {
        assert(index < 2);
        switch (index) {
        case 0:
            return x;
        case 1:
            return y;
        default:
            assert(!"Should not get here");
        }
        return x;
    }

    int& operator[](const size_t index) noexcept
    {
        assert(index < 2);
        switch (index) {
        case 0:
            return x;
        case 1:
            return y;
        default:
            assert(!"Should not get here");
        }
        return x;
    }

    friend bool operator==(const Point2i& a, const Point2i& b) noexcept = default;
};

////////////////////////////////////////////////////////////////////////////////
/// Unary Operators
////////////////////////////////////////////////////////////////////////////////

inline Point2i operator+(const Point2i& a) noexcept
{
    return a;
}

inline Point2i operator-(const Point2i& a) noexcept
{
    return Point2i{ -a.x, -a.y };
}

inline Point2i abs(const Point2i& a) noexcept
{
    return Point2i{ std::abs(a.x), std::abs(a.y) };
}

inline Point2 sqrt(const Point2i& a) noexcept
{
    return sqrt(Point2(a.x, a.y));
}

inline Point2i sqr(const Point2i& a) noexcept
{
    return Point2i{ a.x * a.x, a.y * a.y };
}

inline Point2 rsqrt(const Point2i& a) noexcept
{
    return rsqrt(Point2(a.x, a.y));
}

////////////////////////////////////////////////////////////////////////////////
/// Binary Operators
////////////////////////////////////////////////////////////////////////////////

inline Point2i operator+(const Point2i& a, const Point2i& b) noexcept
{
    return Point2i{ a.x + b.x, a.y + b.y };
}

inline Point2i operator-(const Point2i& a, const Point2i& b) noexcept
{
    return Point2i{ a.x - b.x, a.y - b.y };
}

inline Point2i operator*(const Point2i& a, const Point2i& b) noexcept
{
    return Point2i{ a.x * b.x, a.y * b.y };
}

inline Point2i operator*(const Point2i& a, const int b) noexcept
{
    return a * Point2i{ b };
}

inline Point2i operator*(const int a, const Point2i& b) noexcept
{
    return Point2i{ a } * b;
}

inline Point2i operator/(const Point2i& a, const Point2i& b) noexcept
{
    return Point2i{ a.x / b.x, a.y / b.y };
}

inline Point2i operator/(const Point2i& a, const int b) noexcept
{
    return a / Point2i{ b };
}

inline Point2i operator/(const int a, const Point2i& b) noexcept
{
    return Point2i{ a } / b;
}

inline Point2i min(const Point2i& a, const Point2i& b) noexcept
{
    return Point2i{ std::min(a.x, b.x), std::min(a.y, b.y) };
}

inline Point2i max(const Point2i& a, const Point2i& b) noexcept
{
    return Point2i{ std::max(a.x, b.x), std::max(a.y, b.y) };
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment Operators
////////////////////////////////////////////////////////////////////////////////

inline Point2i& operator+=(Point2i& a, const Point2i& b) noexcept
{
    return a = a + b;
}

inline Point2i& operator-=(Point2i& a, const Point2i& b) noexcept
{
    return a = a - b;
}

inline Point2i& operator*=(Point2i& a, const Point2i& b) noexcept
{
    return a = a * b;
}

inline Point2i& operator*=(Point2i& a, const int b) noexcept
{
    return a = a * b;
}

inline Point2i& operator/=(Point2i& a, const Point2i& b) noexcept
{
    return a = a / b;
}

inline Point2i& operator/=(Point2i& a, const int b) noexcept
{
    return a = a / b;
}

////////////////////////////////////////////////////////////////////////////////
/// Reductions
////////////////////////////////////////////////////////////////////////////////
inline int reduce_add(const Point2i& v) noexcept
{
    return v.x + v.y;
}

inline int reduce_mul(const Point2i& v) noexcept
{
    return v.x * v.y;
}

inline int reduce_min(const Point2i& v) noexcept
{
    return std::min(v.x, v.y);
}

inline int reduce_max(const Point2i& v) noexcept
{
    return std::max(v.x, v.y);
}

inline std::ostream& operator<<(std::ostream& outs, const Point2i& a)
{
    if (outs.iword(sp::k_pretty_print_key) == 1) {
        return outs << '(' << a.x << ", " << a.y << ')';
    } else {
        // This version is compatible with the input operator.
        return outs << a.x << ' ' << a.y;
    }
}

inline std::istream& operator>>(std::istream& ins, Point2i& a)
{
    int x;
    int y;
    ins >> x >> y;
    a = Point2i{ x, y };
    return ins;
}

////////////////////////////////////////////////////////////////////////////////
/// Creation Operators
////////////////////////////////////////////////////////////////////////////////

template <>
struct Uninitialized<Point2i>
{
    operator Point2i() &&
    {
        return Point2i{ no_init };
    }
};
}; // namespace sp
