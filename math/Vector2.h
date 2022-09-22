#pragma once

///@author Keith Jeffery

// Based on Intel Embree's Vec2fa

#include "Math.h"
#include "VectorType.h"

#include "../base/Constants.h"
#include "../base/Util.h"

#include <algorithm>
#include <cassert>
#include <istream>
#include <limits>
#include <ostream>

#include <immintrin.h>

namespace sp {

template <VectorType>
struct alignas(16) BaseVector2
{
    static constexpr int N = 2;

    union
    {
        __m128 m128;

        struct
        {
            float x;
            float y;
            float az;
            float aw;
        };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    ~BaseVector2() = default;

    BaseVector2() noexcept
    : m128(_mm_setzero_ps())
    {
    }

    explicit BaseVector2(NoInitType) noexcept
    {
    }

    explicit BaseVector2(const __m128 a) noexcept
    : m128(a)
    {
    }

    template <VectorType OtherType>
    BaseVector2(const BaseVector2<OtherType>& other) noexcept
    : m128(other.m128)
    {
    }

    BaseVector2(const BaseVector2& other) noexcept
    : m128(other.m128)
    {
    }

    BaseVector2& operator=(const BaseVector2& other) noexcept
    {
        m128 = other.m128;
        return *this;
    }

    explicit BaseVector2(const float a) noexcept
    : m128(_mm_set1_ps(a))
    {
    }

    BaseVector2(const float x, const float y) noexcept
    : m128(_mm_set_ps(y, y, y, x))
    {
    }

    explicit BaseVector2(const __m128i a) noexcept
    : m128(_mm_cvtepi32_ps(a))
    {
    }

    explicit operator const __m128&() const noexcept
    {
        return m128;
    }

    explicit operator __m128&() noexcept
    {
        return m128;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Loads and Stores
    ////////////////////////////////////////////////////////////////////////////////

    static BaseVector2 load(const void* const a) noexcept
    {
        return BaseVector2{ _mm_and_ps(_mm_load_ps(reinterpret_cast<const float*>(a)),
                                       _mm_castsi128_ps(_mm_set_epi32(0, 0, -1, -1))) };
    }

    static BaseVector2 loadu(const void* const a) noexcept
    {
        return BaseVector2{ _mm_and_ps(_mm_loadu_ps(reinterpret_cast<const float*>(a)),
                                       __mm_castsi128_ps(_mm_set_epi32(0, 0, -1, -1))) };
    }

    static void storeu(void* ptr, const BaseVector2& v) noexcept
    {
        _mm_storeu_ps(reinterpret_cast<float*>(ptr), v.m128);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    static BaseVector2 zero() noexcept
    {
        return BaseVector2{};
    }

    static BaseVector2 one() noexcept
    {
        return BaseVector2{ _mm_set1_ps(1.0f) };
    }

    static BaseVector2 positive_infinity() noexcept
    {
        return BaseVector2{ _mm_set1_ps(std::numeric_limits<float>::infinity()) };
    }

    static BaseVector2 negative_infinity() noexcept
    {
        return BaseVector2{ _mm_set1_ps(-std::numeric_limits<float>::infinity()) };
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    const float& operator[](const size_t index) const noexcept
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
    }

    float& operator[](const size_t index) noexcept
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
    }
};

////////////////////////////////////////////////////////////////////////////////
/// Unary Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline BaseVector2<type> operator+(const BaseVector2<type>& a) noexcept
{
    return a;
}

template <VectorType type>
inline BaseVector2<type> operator-(const BaseVector2<type>& a) noexcept
{
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return BaseVector2<type>{ _mm_xor_ps(a.m128, mask) };
}

template <VectorType type>
inline BaseVector2<type> abs(const BaseVector2<type>& a) noexcept
{
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return BaseVector2<type>{ _mm_and_ps(a.m128, mask) };
}

template <VectorType type>
inline BaseVector2<type> sqrt(const BaseVector2<type>& a) noexcept
{
    return BaseVector2<type>{ _mm_sqrt_ps(a.m128) };
}

template <VectorType type>
inline BaseVector2<type> sqr(const BaseVector2<type>& a) noexcept
{
    return BaseVector2<type>{ _mm_mul_ps(a.m128, a.m128) };
}

template <VectorType type>
inline BaseVector2<type> rsqrt(const BaseVector2<type>& a) noexcept
{
#if defined(__aarch64__)
    __m128 r = _mm_rsqrt_ps(a.m128);
    r        = vmulq_f32(r, vrsqrtsq_f32(vmulq_f32(a.m128, r), r));
    r        = vmulq_f32(r, vrsqrtsq_f32(vmulq_f32(a.m128, r), r));
    return r;
#else

#if defined(__AVX512VL__) && !defined(_MSC_VER) // I can't find this function in Visual Studio
    __m128 r = _mm_rsqrt14_ps(a.m128);
#else
    __m128 r = _mm_rsqrt_ps(a.m128);
#endif
    return BaseVector2<type>{ _mm_add_ps(
        _mm_mul_ps(_mm_set1_ps(1.5f), r),
        _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a.m128, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r))) };
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Binary Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline BaseVector2<type> operator+(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_add_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector2<type> operator-(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_sub_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector2<type> operator*(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_mul_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector2<type> operator*(const BaseVector2<type>& a, const float b) noexcept
{
    return a * BaseVector2<type>(b);
}

template <VectorType type>
inline BaseVector2<type> operator*(const float a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>(a) * b;
}

template <VectorType type>
inline BaseVector2<type> operator/(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_div_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector2<type> operator/(const BaseVector2<type>& a, const float b) noexcept
{
    return BaseVector2<type>{ _mm_div_ps(a.m128, _mm_set1_ps(b)) };
}

template <VectorType type>
inline BaseVector2<type> operator/(const float a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_div_ps(_mm_set1_ps(a), b.m128) };
}

template <VectorType type>
inline BaseVector2<type> min(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_min_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector2<type> max(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return BaseVector2<type>{ _mm_max_ps(a.m128, b.m128) };
}

////////////////////////////////////////////////////////////////////////////////
/// Ternary Operators
////////////////////////////////////////////////////////////////////////////////

#if defined(__AVX2__) || defined(__ARM_NEON)
template <VectorType type>
inline BaseVector2<type>
madd(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return BaseVector2<type>{ _mm_fmadd_ps(a.m128, b.m128, c.m128) };
}

template <VectorType type>
inline BaseVector2<type>
msub(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return BaseVector2<type>{ _mm_fmsub_ps(a.m128, b.m128, c.m128) };
}

template <VectorType type>
inline BaseVector2<type>
nmadd(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return BaseVector2<type>{ _mm_fnmadd_ps(a.m128, b.m128, c.m128) };
}

template <VectorType type>
inline BaseVector2<type>
nmsub(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return BaseVector2<type>{ _mm_fnmsub_ps(a.m128, b.m128, c.m128) };
}
#else
inline BaseVector2<type>
madd(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return a * b + c;
}

template <VectorType type>
inline BaseVector2<type>
nmadd(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return -a * b + c;
}

template <VectorType type>
inline BaseVector2<type>
nmsub(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return -a * b - c;
}

template <VectorType type>
inline BaseVector2<type>
msub(const BaseVector2<type>& a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return a * b - c;
}
#endif

template <VectorType type>
inline BaseVector2<type> madd(const float a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return madd(BaseVector2<type>(a), b, c);
}

template <VectorType type>
inline BaseVector2<type> msub(const float a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return msub(BaseVector2<type>(a), b, c);
}

template <VectorType type>
inline BaseVector2<type> nmadd(const float a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return nmadd(BaseVector2<type>(a), b, c);
}

template <VectorType type>
inline BaseVector2<type> nmsub(const float a, const BaseVector2<type>& b, const BaseVector2<type>& c) noexcept
{
    return nmsub(BaseVector2<type>(a), b, c);
}

// Based on Matt Pharr's DifferenceOfProducts:
// https://pharr.org/matt/blog/2019/11/03/difference-of-floats
template <VectorType type>
inline BaseVector2<type> difference_of_products(const BaseVector2<type>& a,
                                                const BaseVector2<type>& b,
                                                const BaseVector2<type>& c,
                                                const BaseVector2<type>& d) noexcept
{
    const auto cd  = c * d;
    const auto err = nmadd(c, d, cd);
    const auto dop = msub(a, b, cd);
    return dop + err;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline BaseVector2<type>& operator+=(BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return a = a + b;
}

template <VectorType type>
inline BaseVector2<type>& operator-=(BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return a = a - b;
}

template <VectorType type>
inline BaseVector2<type>& operator*=(BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return a = a * b;
}

template <VectorType type>
inline BaseVector2<type>& operator*=(BaseVector2<type>& a, const float b) noexcept
{
    return a = a * b;
}

template <VectorType type>
inline BaseVector2<type>& operator/=(BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return a = a / b;
}

template <VectorType type>
inline BaseVector2<type>& operator/=(BaseVector2<type>& a, const float b) noexcept
{
    return a = a / b;
}

////////////////////////////////////////////////////////////////////////////////
/// Reductions
////////////////////////////////////////////////////////////////////////////////
template <VectorType type>
inline float reduce_add(const BaseVector2<type>& v) noexcept
{
    return v.x + v.y;
}

template <VectorType type>
inline float reduce_mul(const BaseVector2<type>& v) noexcept
{
    return v.x * v.y;
}

template <VectorType type>
inline float reduce_min(const BaseVector2<type>& v) noexcept
{
    return std::min(v.x, v.y);
}

template <VectorType type>
inline float reduce_max(const BaseVector2<type>& v) noexcept
{
    return std::max(v.x, v.y);
}

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline bool operator==(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return (_mm_movemask_ps(_mm_cmpeq_ps(a.m128, b.m128)) & 3) == 3;
}

template <VectorType type>
inline bool operator!=(const BaseVector2<type>& a, const BaseVector2<type>& b) noexcept
{
    return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 3) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Euclidean Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline std::ostream& operator<<(std::ostream& outs, const BaseVector2<type>& a)
{
    if (outs.iword(sp::k_pretty_print_key)) {
        return outs << '(' << a.x << ", " << a.y << ')';
    } else {
        // This version is compatible with the input operator.
        return outs << a.x << ' ' << a.y;
    }
}

template <VectorType type>
inline std::istream& operator>>(std::istream& ins, BaseVector2<type>& a)
{
    float x;
    float y;
    ins >> x >> y;
    a = BaseVector2<type>{ x, y };
    return ins;
}

////////////////////////////////////////////////////////////////////////////////
/// Creation Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
struct Uninitialized<BaseVector2<type>>
{
    operator BaseVector2<type>() &&
    {
        return BaseVector2<type>{ no_init };
    }
};

using Vector2 = BaseVector2<VectorType::vector>;
using Normal2 = BaseVector2<VectorType::normal>;
using Point2  = BaseVector2<VectorType::point>;

////////////////////////////////////////////////////////////////////////////////
/// Type-Specialized Operators
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Binary Operators
////////////////////////////////////////////////////////////////////////////////

inline Point2 operator+(const Point2& a, const Vector2& b) noexcept
{
    return Point2{ _mm_add_ps(a.m128, b.m128) };
}

inline Point2 operator+(const Vector2& a, const Point2& b) noexcept
{
    return b + a;
}

inline Vector2 operator-(const Point2& a, const Point2& b) noexcept
{
    return Vector2{ _mm_sub_ps(a.m128, b.m128) };
}

////////////////////////////////////////////////////////////////////////////////
/// Euclidean Space Operators
////////////////////////////////////////////////////////////////////////////////

#if defined(__SSE4_1__)
inline float dot(const Vector2& a, const Vector2& b) noexcept
{
    return _mm_cvtss_f32(_mm_dp_ps(a.m128, b.m128, 0x3F));
}
#else
inline float dot(const Vector2& a, const Vector2& b) noexcept
{
    return reduce_add(a * b);
}
#endif

inline Vector2 cross(const Vector2& a, const Vector2& b) noexcept
{
    return Vector2{ -a.y, a.x };
}

inline float sqr_length(const Vector2& a) noexcept
{
    return dot(a, a);
}

inline float rcp_length(const Vector2& a) noexcept
{
    return rsqrt(dot(a, a));
}

inline float rcp_length2(const Vector2& a) noexcept
{
    return rcp(dot(a, a));
}

inline float length(const Vector2& a) noexcept
{
    return std::sqrt(dot(a, a));
}

inline Vector2 normalize(const Vector2& a) noexcept
{
    return a * rsqrt(dot(a, a));
}

inline bool is_normalized(const Vector2& a) noexcept
{
    return std::abs(1.0f - sqr_length(a)) < std::numeric_limits<float>::epsilon();
}

inline float distance(const Point2& a, const Point2& b) noexcept
{
    return length(a - b);
}

inline Vector2 lerp(const Vector2& v0, const Vector2& v1, const float t) noexcept
{
    return madd(1.0f - t, v0, t * v1);
}

}; // namespace sp
