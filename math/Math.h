
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

// clang-format off
// Wow. Clang-format massacres this! Apparently this version of clang-format is not familiar with concepts.
template <typename T>
concept squarable = requires(T v)
{
    {v * v} -> std::convertible_to<T>;
};

// clang-format on

// clang-format off
template <typename T>
requires squarable<T>
constexpr T square(const T& t) noexcept(noexcept(t * t))
{
    return t * t;
}

// clang-format on

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

inline float balance_heuristic(int nf, float f_pdf, int ng, float g_pdf) noexcept
{
    return balance_heuristic(nf, f_pdf, nf * f_pdf + ng * g_pdf);
}

template <typename T>
requires std::is_integral_v<T>
constexpr bool is_power_of_two(T v) noexcept
{
    return v > 0 && !(v & (v - T{ 1 }));
}

template <typename T>
concept unsigned_32 = std::unsigned_integral<T> && sizeof(T) == 4;

template <typename T>
concept unsigned_64 = std::unsigned_integral<T> && sizeof(T) == 8;

template <typename T>
constexpr T round_up_to_power_of_two(T v) noexcept
requires unsigned_32<T>
{
    --v;
    v |= v >> 1UL;
    v |= v >> 2UL;
    v |= v >> 4UL;
    v |= v >> 8UL;
    v |= v >> 16UL;
    ++v;
    return v;
}

template <typename T>
constexpr T round_up_to_power_of_two(T v) noexcept
requires unsigned_64<T>
{
    --v;
    v |= v >> 1UL;
    v |= v >> 2UL;
    v |= v >> 4UL;
    v |= v >> 8UL;
    v |= v >> 16UL;
    v |= v >> 32UL;
    ++v;
    return v;
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

// Taken from https://stackoverflow.com/a/49743348 and modified to fit my own selfish needs.
inline float erfinv(float a) noexcept
{
    float       p;
    float       r;
    const float t = std::log(madd(a, 0.0f - a, 1.0f));
    // clang-format off
    if (std::abs(t) > 6.125f) { // maximum ulp error = 2.35793
        p =             3.03697567e-10f; //  0x1.4deb44p-32
        p = madd(p, t,  2.93243101e-8f); //  0x1.f7c9aep-26
        p = madd(p, t,  1.22150334e-6f); //  0x1.47e512p-20
        p = madd(p, t,  2.84108955e-5f); //  0x1.dca7dep-16
        p = madd(p, t,  3.93552968e-4f); //  0x1.9cab92p-12
        p = madd(p, t,  3.02698812e-3f); //  0x1.8cc0dep-9
        p = madd(p, t,  4.83185798e-3f); //  0x1.3ca920p-8
        p = madd(p, t, -2.64646143e-1f); // -0x1.0eff66p-2
        p = madd(p, t,  8.40016484e-1f); //  0x1.ae16a4p-1
    } else { // maximum ulp error = 2.35002
        p =              5.43877832e-9f;  //  0x1.75c000p-28
        p = madd(p, t,  1.43285448e-7f); //  0x1.33b402p-23
        p = madd(p, t,  1.22774793e-6f); //  0x1.499232p-20
        p = madd(p, t,  1.12963626e-7f); //  0x1.e52cd2p-24
        p = madd(p, t, -5.61530760e-5f); // -0x1.d70bd0p-15
        p = madd(p, t, -1.47697632e-4f); // -0x1.35be90p-13
        p = madd(p, t,  2.31468678e-3f); //  0x1.2f6400p-9
        p = madd(p, t,  1.15392581e-2f); //  0x1.7a1e50p-7
        p = madd(p, t, -2.32015476e-1f); // -0x1.db2aeep-3
        p = madd(p, t,  8.86226892e-1f); //  0x1.c5bf88p-1
    }
    // clang-format on
    r = a * p;
    return r;
}

// This function is adapted from Burkhard Stubert's page, "Comparing Two Floating-Point Numbers"
// https://embeddeduse.com/2019/08/26/qt-compare-two-floats/
inline bool float_compare(const float a, const float b) noexcept
{
    constexpr float epsilon = 1.0e-05f;
    if (std::abs(a - b) <= epsilon) {
        return true;
    }
    return std::abs(a - b) <= (epsilon * std::max(std::abs(a), std::abs(b)));
}

// Sometimes we're doing random processes (we are a Monte Carlo process, after all), and the results aren't exact
// because of noise and not necessarily floating-point precision. In those cases, it's nice to have a little extra
// wiggle-room.
inline bool float_compare_epsilon(const float a, const float b, const float epsilon) noexcept
{
    return std::abs(a - b) <= epsilon;
}

} // namespace sp
