
#pragma once

///@author Keith Jeffery

#include <cmath>
#include <concepts>
#include <iterator>
#include <numeric>

#include <immintrin.h>
#include <xmmintrin.h>

namespace sp {

template <typename T, typename U>
T safe_divide(const T& a, const U& b) noexcept;

template <>
inline float safe_divide(const float& a, const float& b) noexcept
{
    if (b == 0.0f) {
        return 0.0f;
    } else {
        return a / b;
    }
}

template <typename IteratorValue, typename IteratorPDF>
requires std::forward_iterator<IteratorValue> && std::forward_iterator<IteratorPDF>
float balance_heuristic(int           c,
                        float         p,
                        IteratorValue c_first,
                        IteratorValue c_last,
                        IteratorPDF   pdfs_first,
                        IteratorPDF   pdfs_last) noexcept
{
    assert(std::distance(c_first, c_last) == std::distance(pdfs_first, pdfs_last));
    const float inner_product = std::inner_product(c_first, c_last, pdfs_first, 0.0f);
    if (inner_product == 0.0f) {
        return 0.0f;
    }
    const auto w = (c * p) / inner_product;
    return w;
}

// Specialization for when the inner product is precomputed
inline float balance_heuristic(int c, float p, float inner_product) noexcept
{
    if (inner_product == 0.0f) {
        return 0.0f;
    }
    const auto w = (c * p) / inner_product;
    return w;
}

// Specialization for when we are taking one sample from our distribution and the inner product is precomputed
inline float balance_heuristic(float p, float inner_product) noexcept
{
    if (inner_product == 0.0f) {
        return 0.0f;
    }
    const auto w = p / inner_product;
    return w;
}

template <typename T>
requires std::is_integral_v<T>
constexpr bool is_power_of_two(T v) noexcept
{
    return v > 0 && !(v & (v - T{ 1 }));
}

#if defined(__AVX2__) || defined(__ARM_NEON)
inline float madd(float a, float b, float c) noexcept
{
    return std::fma(a, b, c);
}
#else
inline float madd(float a, float b, float c) noexcept
{
    return a * b + c;
}
#endif

inline float msub(float a, float b, float c) noexcept
{
    return madd(a, b, -c);
}

inline float nmadd(float a, float b, float c) noexcept
{
    return madd(-a, b, c);
}

inline float nmsub(float a, float b, float c) noexcept
{
    return -madd(a, b, c);
}

// Based on Matt Pharr's DifferenceOfProducts:
// https://pharr.org/matt/blog/2019/11/03/difference-of-floats
inline float difference_of_products(float a, float b, float c, float d) noexcept
{
    const float cd  = c * d;
    const float err = sp::madd(-c, d, cd);
    const float dop = sp::madd(a, b, -cd);
    return dop + err;
}

// Intel's Embree library
inline float rcp(const float x) noexcept
{
#if defined(__aarch64__)
    // Move scalar to vector register and do rcp.
    __m128 a;
    a[0]                   = x;
    float32x4_t reciprocal = vrecpeq_f32(a);
    reciprocal             = vmulq_f32(vrecpsq_f32(a, reciprocal), reciprocal);
    reciprocal             = vmulq_f32(vrecpsq_f32(a, reciprocal), reciprocal);
    return reciprocal[0];
#else

    const __m128 a = _mm_set_ss(x);

#if defined(__AVX512VL__)
    const __m128 r = _mm_rcp14_ss(_mm_set_ss(0.0f), a);
#else
    const __m128 r = _mm_rcp_ss(a);
#endif

#if defined(__AVX2__)
    return _mm_cvtss_f32(_mm_mul_ss(r, _mm_fnmadd_ss(r, a, _mm_set_ss(2.0f))));
#else
    return _mm_cvtss_f32(_mm_mul_ss(r, _mm_sub_ss(_mm_set_ss(2.0f), _mm_mul_ss(r, a))));
#endif

#endif // defined(__aarch64__)
}

// Intel's Embree library
inline float rsqrt(const float x) noexcept
{
#if defined(__aarch64__)
    // FP and Neon shares same vector register in arm64
    __m128 a;
    a[0]         = x;
    __m128 value = _mm_rsqrt_ps(a);
    value        = vmulq_f32(value, vrsqrtsq_f32(vmulq_f32(a, value), value));
    value        = vmulq_f32(value, vrsqrtsq_f32(vmulq_f32(a, value), value));
    return value[0];
#else

    const __m128 a = _mm_set_ss(x);
#if defined(__AVX512VL__)
    __m128       r = _mm_rsqrt14_ss(_mm_set_ss(0.0f), a);
#else
    __m128 r = _mm_rsqrt_ss(a);
#endif
    const __m128 c = _mm_add_ss(_mm_mul_ss(_mm_set_ss(1.5f), r),
                                _mm_mul_ss(_mm_mul_ss(_mm_mul_ss(a, _mm_set_ss(-0.5f)), r), _mm_mul_ss(r, r)));
    return _mm_cvtss_f32(c);
#endif
}
} // namespace sp
