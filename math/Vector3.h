#pragma once

///@author Keith Jeffery

// Based on Intel Embree's Vec3fa

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
struct alignas(16) BaseVector3
{
    using scalar = float;

    static constexpr int N = 3;

    union
    {
        __m128 m128;

        struct
        {
            float x;
            float y;
            float z;
        };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    ~BaseVector3() = default;

    BaseVector3() noexcept
    : m128(_mm_setzero_ps())
    {
    }

    explicit BaseVector3(NoInitType) noexcept
    {
    }

    explicit BaseVector3(const __m128 a) noexcept
    : m128(a)
    {
    }

    template <VectorType OtherType>
    BaseVector3(const BaseVector3<OtherType>& other) noexcept
    : m128(other.m128)
    {
    }

    BaseVector3(const BaseVector3& other) noexcept
    : m128(other.m128)
    {
    }

    BaseVector3(BaseVector3&& other) noexcept
    : m128(other.m128)
    {
    }

    BaseVector3& operator=(const BaseVector3& other) noexcept
    {
        m128 = other.m128;
        return *this;
    }

    BaseVector3& operator=(BaseVector3&& other) noexcept
    {
        m128 = other.m128;
        return *this;
    }

    explicit BaseVector3(const float a) noexcept
    : m128(_mm_set1_ps(a))
    {
    }

    BaseVector3(const float x, const float y, const float z) noexcept
    : m128(_mm_set_ps(0, z, y, x))
    {
    }

    explicit BaseVector3(const __m128i a) noexcept
    : m128(_mm_cvtepi32_ps(a))
    {
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Loads and Stores
    ////////////////////////////////////////////////////////////////////////////////

    static BaseVector3 load(const void* const a) noexcept
    {
#if defined(__aarch64__)
        __m128 t = _mm_load_ps((float*)a);
        t[3]     = 0.0f;
        return BaseVector3(t);
#else
        return BaseVector3{
            _mm_and_ps(_mm_load_ps(reinterpret_cast<const float*>(a)),
                       _mm_castsi128_ps(_mm_set_epi32(0, -1, -1, -1)))
        };
#endif
    }

    static BaseVector3 loadu(const void* const a) noexcept
    {
        return BaseVector3{ _mm_loadu_ps(reinterpret_cast<const float*>(a)) };
    }

    static void storeu(void* ptr, const BaseVector3& v) noexcept
    {
        _mm_storeu_ps(reinterpret_cast<float*>(ptr), v.m128);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    static BaseVector3 zero() noexcept
    {
        return BaseVector3{};
    }

    static BaseVector3 one() noexcept
    {
        return BaseVector3{ _mm_set1_ps(1.0f) };
    }

    static BaseVector3 positive_infinity() noexcept
    {
        return BaseVector3{ _mm_set1_ps(std::numeric_limits<float>::infinity()) };
    }

    static BaseVector3 negative_infinity() noexcept
    {
        return BaseVector3{ _mm_set1_ps(-std::numeric_limits<float>::infinity()) };
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    const float& operator[](const size_t index) const noexcept
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
        assert(!"Should not get here");
        return x;
    }

    float& operator[](const size_t index) noexcept
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
        assert(!"Should not get here");
        return x;
    }
};

////////////////////////////////////////////////////////////////////////////////
/// Movement/Shifting/Shuffling Functions
////////////////////////////////////////////////////////////////////////////////

#if defined(__aarch64__)
template <type, int i0, int i1, int i2, int i3>
inline __m128 shuffle(const __m128& v) noexcept
{
    return vreinterpretq_f32_u8(vqtbl1q_u8((uint8x16_t)v.m128, _MN_SHUFFLE(i0, i1, i2, i3)));
}

template <int i0, int i1, int i2, int i3>
inline __m128 shuffle(const __m128& a, const __m128& b) noexcept
{
    return vreinterpretq_f32_u8(
        vqtbl2q_u8((uint8x16x2_t){ (uint8x16_t)a.m128, (uint8x16_t)b.m128 }, _MF_SHUFFLE(i0, i1, i2, i3)));
}
#else
template <int i0, int i1, int i2, int i3>
__m128 shuffle(const __m128& v) noexcept
{
    return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(v), _MM_SHUFFLE(i3, i2, i1, i0)));
}

template <int i0, int i1, int i2, int i3>
inline __m128 shuffle(const __m128& a, const __m128& b) noexcept
{
    return _mm_shuffle_ps(a, b, _MM_SHUFFLE(i3, i2, i1, i0));
}
#endif

#if defined(__SSE3__) && !defined(__aarch64__)
template <>
inline __m128 shuffle<0, 0, 2, 2>(const __m128& v) noexcept
{
    return _mm_moveldup_ps(v);
}

template <>
inline __m128 shuffle<1, 1, 3, 3>(const __m128& v) noexcept
{
    return _mm_movehdup_ps(v);
}

template <>
inline __m128 shuffle<0, 1, 0, 1>(const __m128& v) noexcept
{
    return _mm_castpd_ps(_mm_movedup_pd(_mm_castps_pd(v)));
}
#endif

template <VectorType type, int i0, int i1, int i2, int i3>
inline BaseVector3<type> shuffle(const BaseVector3<type>& v) noexcept
{
    return BaseVector3<type>{ shuffle<i0, i1, i2, i3>(v.m128) };
}

template <VectorType type, int i0, int i1, int i2, int i3>
inline BaseVector3<type> shuffle(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ shuffle<i0, i1, i2, i3>(a.m128, b.m128) };
}

template <VectorType type, int i>
inline BaseVector3<type> shuffle(const BaseVector3<type>& v) noexcept
{
    return shuffle<type, i, i, i, i>(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Unary Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline BaseVector3<type> operator+(const BaseVector3<type>& a) noexcept
{
    return a;
}

template <VectorType type>
inline BaseVector3<type> operator-(const BaseVector3<type>& a) noexcept
{
#if defined(__aarch64__)
    return vnegq_f32(a.m128);
#else
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return BaseVector3<type>{ _mm_xor_ps(a.m128, mask) };
#endif
}

template <VectorType type>
inline BaseVector3<type> abs(const BaseVector3<type>& a) noexcept
{
#if defined(__aarch64__)
    return _mm_abs_ps(a.m128);
#else
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return BaseVector3<type>{ _mm_and_ps(a.m128, mask) };
#endif
}

template <VectorType type>
inline BaseVector3<type> sqrt(const BaseVector3<type>& a) noexcept
{
    return BaseVector3<type>{ _mm_sqrt_ps(a.m128) };
}

template <VectorType type>
inline BaseVector3<type> sqr(const BaseVector3<type>& a) noexcept
{
    return BaseVector3<type>{ _mm_mul_ps(a.m128, a.m128) };
}

template <VectorType type>
inline BaseVector3<type> rsqrt(const BaseVector3<type>& a) noexcept
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
    return BaseVector3<type>{
        _mm_add_ps(
            _mm_mul_ps(_mm_set1_ps(1.5f), r),
            _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a.m128, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r)))
    };
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Binary Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline BaseVector3<type> operator+(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_add_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector3<type> operator-(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_sub_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector3<type> operator*(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_mul_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector3<type> operator*(const BaseVector3<type>& a, const float b) noexcept
{
    return a * BaseVector3<type>(b);
}

template <VectorType type>
inline BaseVector3<type> operator*(const float a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>(a) * b;
}

template <VectorType type>
inline BaseVector3<type> operator/(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_div_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector3<type> operator/(const BaseVector3<type>& a, const float b) noexcept
{
    return BaseVector3<type>{ _mm_div_ps(a.m128, _mm_set1_ps(b)) };
}

template <VectorType type>
inline BaseVector3<type> operator/(const float a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_div_ps(_mm_set1_ps(a), b.m128) };
}

template <VectorType type>
inline BaseVector3<type> min(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_min_ps(a.m128, b.m128) };
}

template <VectorType type>
inline BaseVector3<type> max(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return BaseVector3<type>{ _mm_max_ps(a.m128, b.m128) };
}

////////////////////////////////////////////////////////////////////////////////
/// Ternary Operators
////////////////////////////////////////////////////////////////////////////////

#if defined(__AVX2__) || defined(__ARM_NEON)
template <VectorType type>
inline BaseVector3<type>
madd(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return BaseVector3<type>{ _mm_fmadd_ps(a.m128, b.m128, c.m128) };
}

template <VectorType type>
inline BaseVector3<type>
msub(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return BaseVector3<type>{ _mm_fmsub_ps(a.m128, b.m128, c.m128) };
}

template <VectorType type>
inline BaseVector3<type>
nmadd(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return BaseVector3<type>{ _mm_fnmadd_ps(a.m128, b.m128, c.m128) };
}

template <VectorType type>
inline BaseVector3<type>
nmsub(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return BaseVector3<type>{ _mm_fnmsub_ps(a.m128, b.m128, c.m128) };
}
#else
inline BaseVector3<type>
madd(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return a * b + c;
}

template <VectorType type>
inline BaseVector3<type>
nmadd(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return -a * b + c;
}

template <VectorType type>
inline BaseVector3<type>
nmsub(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return -a * b - c;
}

template <VectorType type>
inline BaseVector3<type>
msub(const BaseVector3<type>& a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return a * b - c;
}
#endif

template <VectorType type>
inline BaseVector3<type> madd(const float a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return madd(BaseVector3<type>(a), b, c);
}

template <VectorType type>
inline BaseVector3<type> madd(const BaseVector3<type>& a, const float b, const BaseVector3<type>& c) noexcept
{
    return madd(a, BaseVector3<type>(b), c);
}

template <VectorType type>
inline BaseVector3<type> msub(const float a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return msub(BaseVector3<type>(a), b, c);
}

template <VectorType type>
inline BaseVector3<type> nmadd(const float a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return nmadd(BaseVector3<type>(a), b, c);
}

template <VectorType type>
inline BaseVector3<type> nmsub(const float a, const BaseVector3<type>& b, const BaseVector3<type>& c) noexcept
{
    return nmsub(BaseVector3<type>(a), b, c);
}

// Based on Matt Pharr's DifferenceOfProducts:
// https://pharr.org/matt/blog/2019/11/03/difference-of-floats
template <VectorType type>
inline BaseVector3<type> difference_of_products(const BaseVector3<type>& a,
                                                const BaseVector3<type>& b,
                                                const BaseVector3<type>& c,
                                                const BaseVector3<type>& d) noexcept
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
inline BaseVector3<type>& operator+=(BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return a = a + b;
}

template <VectorType type>
inline BaseVector3<type>& operator-=(BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return a = a - b;
}

template <VectorType type>
inline BaseVector3<type>& operator*=(BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return a = a * b;
}

template <VectorType type>
inline BaseVector3<type>& operator*=(BaseVector3<type>& a, const float b) noexcept
{
    return a = a * b;
}

template <VectorType type>
inline BaseVector3<type>& operator/=(BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return a = a / b;
}

template <VectorType type>
inline BaseVector3<type>& operator/=(BaseVector3<type>& a, const float b) noexcept
{
    return a = a / b;
}

////////////////////////////////////////////////////////////////////////////////
/// Reductions
////////////////////////////////////////////////////////////////////////////////
#if defined(__aarch64__)
template <VectorType type>
inline float reduce_add(const BaseVector3<type>& v) noexcept
{
    float32x4_t t = v.m128;
    t[3]          = 0.0f;
    return vaddvq_f32(t);
}

template <VectorType type>
inline float reduce_mul(const BaseVector3<type>& v) noexcept
{
    return v.x * v.y * v.z;
}

template <VectorType type>
inline float reduce_min(const BaseVector3<type>& v) noexcept
{
    float32x4_t t = v.m128;
    t[3]          = t[2];
    return vminvq_f32(t);
}

template <VectorType type>
inline float reduce_max(const BaseVector3<type>& v) noexcept
{
    float32x4_t t = v.m128;
    t[3]          = t[2];
    return vmaxvq_f32(t);
}
#else
template <VectorType type>
inline float reduce_add(const BaseVector3<type>& v) noexcept
{
    const BaseVector3<type> a(v.m128);
    const BaseVector3<type> b = shuffle<1>(a);
    const BaseVector3<type> c = shuffle<2>(a);
    return _mm_cvtss_f32((a + b + c).m128);
}

template <VectorType type>
inline float reduce_mul(const BaseVector3<type>& v) noexcept
{
    return v.x * v.y * v.z;
}

template <VectorType type>
inline float reduce_min(const BaseVector3<type>& v) noexcept
{
    return std::ranges::min({ v.x, v.y, v.z });
}

template <VectorType type>
inline float reduce_max(const BaseVector3<type>& v) noexcept
{
    return std::ranges::min({ v.x, v.y, v.z });
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline bool operator==(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return (_mm_movemask_ps(_mm_cmpeq_ps(a.m128, b.m128)) & 7) == 7;
}

template <VectorType type>
inline bool operator!=(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 7) != 0;
}

template <VectorType type>
inline bool compare(const BaseVector3<type>& a, const BaseVector3<type>& b) noexcept
{
    return float_compare(a.x, b.x) && float_compare(a.y, b.y) && float_compare(a.z, b.z);
}

template <VectorType type>
inline bool compare_epsilon(const BaseVector3<type>& a, const BaseVector3<type>& b, const float epsilon) noexcept
{
    return float_compare_epsilon(a.x, b.x, epsilon) && float_compare_epsilon(a.y, b.y, epsilon) &&
            float_compare_epsilon(a.z, b.z, epsilon);
}

////////////////////////////////////////////////////////////////////////////////
/// Euclidean Space Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline int max_dim(const BaseVector3<type>& a) noexcept
{
    const BaseVector3<type> b = abs(a);
    if (b.x > b.y) {
        if (b.x > b.z) {
            return 0;
        } else {
            return 2;
        }
    } else {
        if (b.y > b.z) {
            return 1;
        } else {
            return 2;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Input/Output Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
inline std::ostream& operator<<(std::ostream& outs, const BaseVector3<type>& a)
{
    if (outs.iword(sp::k_pretty_print_key) == 1) {
        return outs << '(' << a.x << ", " << a.y << ", " << a.z << ')';
    } else {
        // This version is compatible with the input operator.
        return outs << a.x << ' ' << a.y << ' ' << a.z;
    }
}

template <VectorType type>
inline std::istream& operator>>(std::istream& ins, BaseVector3<type>& a)
{
    float x;
    float y;
    float z;
    ins >> x >> y >> z;
    a = BaseVector3<type>{ x, y, z };
    return ins;
}

////////////////////////////////////////////////////////////////////////////////
/// Creation Operators
////////////////////////////////////////////////////////////////////////////////

template <VectorType type>
struct Uninitialized<BaseVector3<type>>
{
    operator BaseVector3<type>() &&
    {
        return BaseVector3<type>{ no_init };
    }
};

using Vector3 = BaseVector3<VectorType::vector>;
using Normal3 = BaseVector3<VectorType::normal>;
using Point3  = BaseVector3<VectorType::point>;

////////////////////////////////////////////////////////////////////////////////
/// Type-Specialized Operators
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Binary Operators
////////////////////////////////////////////////////////////////////////////////

inline Point3 operator+(const Point3& a, const Vector3& b) noexcept
{
    return Point3{ _mm_add_ps(a.m128, b.m128) };
}

inline Point3 operator+(const Vector3& a, const Point3& b) noexcept
{
    return b + a;
}

inline Vector3 operator-(const Point3& a, const Point3& b) noexcept
{
    return Vector3{ _mm_sub_ps(a.m128, b.m128) };
}

////////////////////////////////////////////////////////////////////////////////
/// Euclidean Space Operators
////////////////////////////////////////////////////////////////////////////////

#if defined(__SSE4_1__)
inline float dot(const Vector3& a, const Vector3& b) noexcept
{
    return _mm_cvtss_f32(_mm_dp_ps(a.m128, b.m128, 0x7F));
}
#else
inline float dot(const Vector3& a, const Vector3& b) noexcept
{
    return reduce_add(a * b);
}
#endif

inline float dot(const Vector3& a, const Normal3& b) noexcept
{
    return dot(a, Vector3{ b });
}

inline float dot(const Normal3& a, const Vector3& b) noexcept
{
    return dot(b, a);
}

inline float dot(const Normal3& a, const Normal3& b) noexcept
{
    return dot(Vector3{ a }, Vector3{ b });
}

inline Vector3 cross(const Vector3& a, const Vector3& b) noexcept
{
    return difference_of_products(shuffle<VectorType::vector, 1, 2, 0, 3>(a),
                                  shuffle<VectorType::vector, 2, 0, 1, 3>(b),
                                  shuffle<VectorType::vector, 2, 0, 1, 3>(a),
                                  shuffle<VectorType::vector, 1, 2, 0, 3>(b));
}

inline float sqr_length(const Vector3& a) noexcept
{
    return dot(a, a);
}

inline float rcp_length(const Vector3& a) noexcept
{
    return rsqrt(dot(a, a));
}

inline float rcp_length2(const Vector3& a) noexcept
{
    return rcp(dot(a, a));
}

inline float length(const Vector3& a) noexcept
{
    return std::sqrt(dot(a, a));
}

inline Vector3 normalize(const Vector3& a) noexcept
{
    return a * rsqrt(dot(a, a));
}

inline Normal3 normalize(const Normal3& a) noexcept
{
    return a * rsqrt(dot(a, a));
}

inline bool is_normalized(const Vector3& a) noexcept
{
    return std::abs(1.0f - sqr_length(a)) < 0.0001f; // Arbitrary
}

inline float distance(const Point3& a, const Point3& b) noexcept
{
    return length(a - b);
}

inline float half_area(const Vector3& d) noexcept
{
    return madd(d.x, (d.y + d.z), d.y * d.z);
}

inline float area(const Vector3& d) noexcept
{
    return 2.0f * half_area(d);
}

inline Vector3 normalize_safe(const Vector3& a) noexcept
{
    const float d = dot(a, a);
    if (d == 0.0f) [[unlikely]] {
        return a;
    } else {
        return a * rsqrt(d);
    }
}

/*! differentiated normalization */
inline Vector3 dnormalize(const Vector3& p, const Vector3& dp) noexcept
{
    const float pp  = dot(p, p);
    const float pdp = dot(p, dp);
    return (pp * dp - pdp * p) * rcp(pp) * rsqrt(pp);
}

inline Vector3 lerp(const Vector3& v0, const Vector3& v1, const float t) noexcept
{
    return madd(1.0f - t, v0, t * v1);
}
}; // namespace sp
