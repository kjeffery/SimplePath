#pragma once

/// @author Keith Jeffery
//
#include "Vector3.h"

#include <cmath>

namespace sp {
class ONB
{
    static auto create(const Vector3& n) noexcept
    {
        // From Pixar's "Building an Orthonormal Basis, Revisited"
        // http://jcgt.org/published/0006/01/01/paper.pdf

        const auto sign = std::copysign(1.0f, n.z);
        const auto a    = -1.0f / (sign + n.z);
        const auto b    = n.x * n.y * a;
        const auto b1   = Vector3(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
        const auto b2   = Vector3(b, sign + n.y * n.y * a, -n.y);

        // l1 = b1.length()
        // l2 = b2.length()
        //
        // c1 = dot(n, b1)
        // c2 = dot(n, b2)
        // c3 = dot(b1, b2)
        //
        // v1 = cross(b1, b2)
        return std::make_pair(b1, b2);
    }

public:
    ONB(Vector3 u, Vector3 v, Vector3 w) noexcept
    : m_u(std::move(u))
    , m_v(std::move(v))
    , m_w(std::move(w))
    {
        assert(is_normalized(m_u));
        assert(is_normalized(m_v));
        assert(is_normalized(m_w));
    }

    static ONB from_u(const Vector3& n) noexcept
    {
        const auto u      = normalize(n);
        const auto [v, w] = create(u);
        return { u, v, w };
    }

    static ONB from_u(const Normal3& n) noexcept
    {
        return from_u(Vector3{ n });
    }

    static ONB from_v(const Vector3& n) noexcept
    {
        const auto v      = normalize(n);
        const auto [w, u] = create(v);
        return { u, v, w };
    }

    static ONB from_v(const Normal3& n) noexcept
    {
        return from_v(Vector3{ n });
    }

    static ONB from_w(const Vector3& n) noexcept
    {
        const auto w      = normalize(n);
        const auto [u, v] = create(w);
        return { u, v, w };
    }

    static ONB from_w(const Normal3& n) noexcept
    {
        return from_w(Vector3{ n });
    }

    static ONB from_uv(Vector3 u, Vector3 v) noexcept
    {
        u      = normalize(u);
        auto w = cross(u, v);
        w      = normalize(w);
        v      = cross(w, u);
        return { u, v, w };
    }

    static ONB from_vu(Vector3 v, Vector3 u) noexcept
    {
        v      = normalize(v);
        auto w = cross(u, v);
        w      = normalize(w);
        u      = cross(v, w);
        return {u, v, w};
    }

    static ONB from_uw(Vector3 u, const Vector3& w) noexcept
    {
        u      = normalize(u);
        auto v = cross(w, u);
        v      = normalize(v);
        u      = cross(v, w);
        return {u, v, w};
    }

    static ONB from_wu(Vector3 w, Vector3 u) noexcept
    {
        w      = normalize(w);
        auto v = cross(w, u);
        v      = normalize(v);
        u      = cross(v, w);
        return {u, v, w};
    }

    static ONB from_vw(Vector3 v, Vector3 w) noexcept
    {
        v      = normalize(v);
        auto u = cross(v, w);
        u      = normalize(u);
        w      = cross(u, v);
        return {u, v, w};
    }

    static ONB from_wv(Vector3 w, Vector3 v) noexcept
    {
        w      = normalize(w);
        auto u = cross(v, w);
        u      = normalize(u);
        v      = cross(w, u);
        return {u, v, w};
    }

    [[nodiscard]] Vector3 to_world(const Vector3& a) const noexcept
    {
        return a.x * m_u + a.y * m_v + a.z * m_w;
    }

    [[nodiscard]] Vector3 to_onb(const Vector3& a) const noexcept
    {
        return Vector3{ dot(a, m_u), dot(a, m_v), dot(a, m_w) };
    }

private:
    Vector3 m_u;
    Vector3 m_v;
    Vector3 m_w;
};
} // namespace sp
