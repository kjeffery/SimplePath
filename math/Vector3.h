#pragma once

///@author Keith Jeffery

// Based on Intel's Embree Vec3fa

#include "Math.h"

#include "../base/Constants.h"

#include <cassert>
#include <limits>
#include <ostream>

#include <immintrin.h>

namespace sp {

struct alignas(16) Vector3
{
    static constexpr int N = 3;

    union
    {
        __m128 m128;

        struct
        {
            float x, y, z;
        };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    Vector3() noexcept
    : m128(_mm_setzero_ps())
    {
    }

    explicit Vector3(NoInitType) noexcept
    {
    }

    Vector3(const __m128 a) noexcept
    : m128(a)
    {
    }

    Vector3(const Vector3& other) noexcept
    : m128(other.m128)
    {
    }

    Vector3& operator=(const Vector3& other) noexcept
    {
        m128 = other.m128;
        return *this;
    }

    explicit Vector3(const float a) noexcept
    : m128(_mm_set1_ps(a))
    {
    }

    Vector3(const float x, const float y, const float z) noexcept
    : m128(_mm_set_ps(0, z, y, x))
    {
    }

    explicit Vector3(const __m128i a)
    : m128(_mm_cvtepi32_ps(a))
    {
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Loads and Stores
    ////////////////////////////////////////////////////////////////////////////////

    static Vector3 load(const void* const a)
    {
#if defined(__aarch64__)
        __m128 t = _mm_load_ps((float*)a);
        t[3]     = 0.0f;
        return Vector3(t);
#else
        return Vector3(_mm_and_ps(_mm_load_ps((float*)a), _mm_castsi128_ps(_mm_set_epi32(0, -1, -1, -1))));
#endif
    }

    static Vector3 loadu(const void* const a)
    {
        return Vector3(_mm_loadu_ps((float*)a));
    }

    static void storeu(void* ptr, const Vector3& v)
    {
        _mm_storeu_ps((float*)ptr, v.m128);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    static Vector3 zero() noexcept
    {
        return Vector3{};
    }

    static Vector3 one() noexcept
    {
        return Vector3{ _mm_set1_ps(1.0f) };
    }

    static Vector3 positive_infinity() noexcept
    {
        return Vector3{ _mm_set1_ps(std::numeric_limits<float>::infinity()) };
    }

    static Vector3 negative_infinity() noexcept
    {
        return Vector3{ _mm_set1_ps(-std::numeric_limits<float>::infinity()) };
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    const float& operator[](const size_t index) const
    {
        assert(index < 3);
        switch (index) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            assert(!"Should not get here");
        }
    }

    float& operator[](const size_t index)
    {
        assert(index < 3);
        switch (index) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            assert(!"Should not get here");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
/// Movement/Shifting/Shuffling Functions
////////////////////////////////////////////////////////////////////////////////

#if defined(__aarch64__)
template <int i0, int i1, int i2, int i3>
inline Vector3 shuffle(const Vector3& v)
{
    return vreinterpretq_f32_u8(vqtbl1q_u8((uint8x16_t)v.m128, _MN_SHUFFLE(i0, i1, i2, i3)));
}

template <int i0, int i1, int i2, int i3>
inline Vector3 shuffle(const Vector3& a, const Vector3& b)
{
    return vreinterpretq_f32_u8(
        vqtbl2q_u8((uint8x16x2_t){ (uint8x16_t)a.m128, (uint8x16_t)b.m128 }, _MF_SHUFFLE(i0, i1, i2, i3)));
}
#else
template <int i0, int i1, int i2, int i3>
inline Vector3 shuffle(const Vector3& v)
{
    return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(v.m128), _MM_SHUFFLE(i3, i2, i1, i0)));
}

template <int i0, int i1, int i2, int i3>
inline Vector3 shuffle(const Vector3& a, const Vector3& b)
{
    return _mm_shuffle_ps(a, b, _MM_SHUFFLE(i3, i2, i1, i0));
}
#endif

#if defined(__SSE3__) && !defined(__aarch64__)
template <>
inline Vector3 shuffle<0, 0, 2, 2>(const Vector3& v)
{
    return _mm_moveldup_ps(v.m128);
}

template <>
inline Vector3 shuffle<1, 1, 3, 3>(const Vector3& v)
{
    return _mm_movehdup_ps(v.m128);
}

template <>
inline Vector3 shuffle<0, 1, 0, 1>(const Vector3& v)
{
    return _mm_castpd_ps(_mm_movedup_pd(_mm_castps_pd(v.m128)));
}
#endif

template <int i>
inline Vector3 shuffle(const Vector3& v)
{
    return shuffle<i, i, i, i>(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Unary Operators
////////////////////////////////////////////////////////////////////////////////

inline Vector3 operator+(const Vector3& a)
{
    return a;
}

inline Vector3 operator-(const Vector3& a)
{
#if defined(__aarch64__)
    return vnegq_f32(a.m128);
#else
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_xor_ps(a.m128, mask);
#endif
}

inline Vector3 abs(const Vector3& a)
{
#if defined(__aarch64__)
    return _mm_abs_ps(a.m128);
#else
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return _mm_and_ps(a.m128, mask);
#endif
}

inline Vector3 sqrt(const Vector3& a)
{
    return _mm_sqrt_ps(a.m128);
}

inline Vector3 sqr(const Vector3& a)
{
    return _mm_mul_ps(a.m128, a.m128);
}

inline Vector3 rsqrt(const Vector3& a)
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
    return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f), r),
                      _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a.m128, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r)));
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Binary Operators
////////////////////////////////////////////////////////////////////////////////

inline Vector3 operator+(const Vector3& a, const Vector3& b)
{
    return _mm_add_ps(a.m128, b.m128);
}

inline Vector3 operator-(const Vector3& a, const Vector3& b)
{
    return _mm_sub_ps(a.m128, b.m128);
}

inline Vector3 operator*(const Vector3& a, const Vector3& b)
{
    return _mm_mul_ps(a.m128, b.m128);
}

inline Vector3 operator*(const Vector3& a, const float b)
{
    return a * Vector3(b);
}

inline Vector3 operator*(const float a, const Vector3& b)
{
    return Vector3(a) * b;
}

inline Vector3 operator/(const Vector3& a, const Vector3& b)
{
    return _mm_div_ps(a.m128, b.m128);
}

inline Vector3 operator/(const Vector3& a, const float b)
{
    return _mm_div_ps(a.m128, _mm_set1_ps(b));
}

inline Vector3 operator/(const float a, const Vector3& b)
{
    return _mm_div_ps(_mm_set1_ps(a), b.m128);
}

inline Vector3 min(const Vector3& a, const Vector3& b)
{
    return _mm_min_ps(a.m128, b.m128);
}

inline Vector3 max(const Vector3& a, const Vector3& b)
{
    return _mm_max_ps(a.m128, b.m128);
}

////////////////////////////////////////////////////////////////////////////////
/// Ternary Operators
////////////////////////////////////////////////////////////////////////////////

#if defined(__AVX2__) || defined(__ARM_NEON)
inline Vector3 madd(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return _mm_fmadd_ps(a.m128, b.m128, c.m128);
}

inline Vector3 msub(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return _mm_fmsub_ps(a.m128, b.m128, c.m128);
}

inline Vector3 nmadd(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return _mm_fnmadd_ps(a.m128, b.m128, c.m128);
}

inline Vector3 nmsub(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return _mm_fnmsub_ps(a.m128, b.m128, c.m128);
}
#else
inline Vector3 madd(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return a * b + c;
}

inline Vector3 nmadd(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return -a * b + c;
}

inline Vector3 nmsub(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return -a * b - c;
}

inline Vector3 msub(const Vector3& a, const Vector3& b, const Vector3& c)
{
    return a * b - c;
}
#endif

inline Vector3 madd(const float a, const Vector3& b, const Vector3& c)
{
    return madd(Vector3(a), b, c);
}

inline Vector3 msub(const float a, const Vector3& b, const Vector3& c)
{
    return msub(Vector3(a), b, c);
}

inline Vector3 nmadd(const float a, const Vector3& b, const Vector3& c)
{
    return nmadd(Vector3(a), b, c);
}

inline Vector3 nmsub(const float a, const Vector3& b, const Vector3& c)
{
    return nmsub(Vector3(a), b, c);
}

// Based on Matt Pharr's DifferenceOfProducts:
// https://pharr.org/matt/blog/2019/11/03/difference-of-floats
inline Vector3 difference_of_products(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d)
{
    const auto cd  = c * d;
    const auto err = nmadd(c, d, cd);
    const auto dop = msub(a, b, cd);
    return dop + err;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment Operators
////////////////////////////////////////////////////////////////////////////////

inline Vector3& operator+=(Vector3& a, const Vector3& b)
{
    return a = a + b;
}

inline Vector3& operator-=(Vector3& a, const Vector3& b)
{
    return a = a - b;
}

inline Vector3& operator*=(Vector3& a, const Vector3& b)
{
    return a = a * b;
}

inline Vector3& operator*=(Vector3& a, const float b)
{
    return a = a * b;
}

inline Vector3& operator/=(Vector3& a, const Vector3& b)
{
    return a = a / b;
}

inline Vector3& operator/=(Vector3& a, const float b)
{
    return a = a / b;
}

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////

inline bool operator==(const Vector3& a, const Vector3& b)
{
    return (_mm_movemask_ps(_mm_cmpeq_ps(a.m128, b.m128)) & 7) == 7;
}

inline bool operator!=(const Vector3& a, const Vector3& b)
{
    return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 7) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Euclidean Space Operators
////////////////////////////////////////////////////////////////////////////////

#if defined(__SSE4_1__)
inline float dot(const Vector3& a, const Vector3& b)
{
    return _mm_cvtss_f32(_mm_dp_ps(a.m128, b.m128, 0x7F));
}
#else
#error Add reduce add
#endif

inline Vector3 cross(const Vector3& a, const Vector3& b)
{
    return difference_of_products(shuffle<1, 2, 0, 3>(a),
                                  shuffle<2, 0, 1, 3>(b),
                                  shuffle<2, 0, 1, 3>(a),
                                  shuffle<1, 2, 0, 3>(b));
}

inline float sqr_length(const Vector3& a)
{
    return dot(a, a);
}

inline float rcp_length(const Vector3& a)
{
    return rsqrt(dot(a, a));
}

inline float rcp_length2(const Vector3& a)
{
    return rcp(dot(a, a));
}

inline float length(const Vector3& a)
{
    return std::sqrt(dot(a, a));
}

inline Vector3 normalize(const Vector3& a)
{
    return a * rsqrt(dot(a, a));
}

inline float distance(const Vector3& a, const Vector3& b)
{
    return length(a - b);
}

inline float halfArea(const Vector3& d)
{
    return madd(d.x, (d.y + d.z), d.y * d.z);
}

inline float area(const Vector3& d)
{
    return 2.0f * halfArea(d);
}

inline Vector3 normalize_safe(const Vector3& a)
{
    const float d = dot(a, a);
    if (d == 0.0f) [[unlikely]]
        return a;
    else
        return a * rsqrt(d);
}

/*! differentiated normalization */
inline Vector3 dnormalize(const Vector3& p, const Vector3& dp)
{
    const float pp  = dot(p, p);
    const float pdp = dot(p, dp);
    return (pp * dp - pdp * p) * rcp(pp) * rsqrt(pp);
}

inline Vector3 lerp(const Vector3& v0, const Vector3& v1, const float t)
{
    return madd(1.0f - t, v0, t * v1);
}

inline int max_dim(const Vector3& a)
{
    const Vector3 b = abs(a);
    if (b.x > b.y) {
        if (b.x > b.z)
            return 0;
        else
            return 2;
    } else {
        if (b.y > b.z)
            return 1;
        else
            return 2;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Output Operators
////////////////////////////////////////////////////////////////////////////////

inline std::ostream& operator<<(std::ostream& cout, const Vector3& a)
{
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
}

typedef Vector3 Vec3fa_t;

} // namespace sp
