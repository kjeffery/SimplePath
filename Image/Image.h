#pragma once

/// @author Keith Jeffery

#include "../base/Array2D.h"
#include "../base/Constants.h"
#include "../base/Endian.h"
#include "../math/RGB.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iosfwd>
#include <stdexcept>

namespace sp {
class ImageError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

using Image = Array2D<RGB>;

void  write_ppm(std::ostream& outs, const Image& img);
void  write_ppm(const std::filesystem::path& file, const Image& img);
void  write_pfm(std::ostream& outs, const Image& img);
void  write_pfm(const std::filesystem::path& file, const Image& img);
void  write(const std::filesystem::path& file, const Image& img);
Image read_pfm(std::istream& ins);
Image read_pfm(const std::filesystem::path& file);
Image read(const std::filesystem::path& file);

Image& operator*=(Image& img, const RGB& c) noexcept;
Image  operator*(Image img, const RGB& c) noexcept;
Image  operator*(const RGB& c, const Image& img) noexcept;

inline float rgb_to_srgb(const float u) noexcept
{
    if (u <= 0.0031308f) {
        return 12.92f * u;
    } else {
        return 1.055f * std::pow(u, 1.0f / 2.4f) - 0.055f;
    }
}

inline RGB rgb_to_srgb(const RGB& c) noexcept
{
    return RGB{ rgb_to_srgb(c.r), rgb_to_srgb(c.g), rgb_to_srgb(c.b) };
}

struct RemapNone
{
    auto operator()(const float f) const noexcept -> float
    {
        return f;
    }
};

struct RemapClamp
{
    auto operator()(const float f) const noexcept -> float
    {
        return std::clamp(f, 0.0f, k_max_less_than_one);
    }
};

struct RemapBlack
{
    auto operator()(const float f) const noexcept -> float
    {
        if (f < 0.0f || f >= 1.0f) {
            return 0.0f;
        }
        return f;
    }
};

struct RemapRepeat
{
    auto operator()(const float f) const noexcept -> float
    {
        return std::abs(std::fmod(f, 1.0f));
    }
};

struct RemapWrap
{
    auto operator()(const float f) const noexcept -> float
    {
        return std::fmod(1.0f + std::fmod(f, 1.0f), 1.0f);
    }
};

template <typename RemapHorizontal, typename RemapVertical>
RGB sample_nearest_neighbor(const Image& img, float s, float t, RemapHorizontal remap_horizontal, RemapVertical remap_vertical)
{
    s = remap_horizontal(s);
    t = remap_vertical(t);

    assert(s >= 0.0f);
    assert(t >= 0.0f);
    assert(s < 1.0f);
    assert(t < 1.0f);

    using size_type = Image::size_type;

    const auto u = std::round(s * static_cast<float>(img.width()));
    const auto v = std::round(t * static_cast<float>(img.height()));

    const auto x = std::min(static_cast<size_type>(u), img.width() - 1);
    const auto y = std::min(static_cast<size_type>(v), img.height() - 1);

    return img(x, y);
}

template <typename Remap = RemapNone>
RGB sample_nearest_neighbor(const Image& img, float s, float t, Remap remap = Remap{})
{
    return sample_nearest_neighbor(img, s, t, remap, remap);
}

template <typename RemapHorizontal, typename RemapVertical>
inline RGB sample_bilinear(const Image& img, float s, float t, RemapHorizontal remap_horizontal, RemapVertical remap_vertical)
{
    s = remap_horizontal(s);
    t = remap_vertical(t);

    assert(s >= 0.0f);
    assert(t >= 0.0f);
    assert(s < 1.0f);
    assert(t < 1.0f);

    using size_type = Image::size_type;

    const auto u       = s * static_cast<float>(img.width());
    const auto v       = t * static_cast<float>(img.height());
    const auto u_lower = std::floor(u);
    const auto u_upper = std::ceil(u);
    const auto v_lower = std::floor(v);
    const auto v_upper = std::ceil(v);

    const auto u_bias = u_upper - u_lower;
    const auto v_bias = v_upper - v_lower;

    const auto x_lower = std::min(static_cast<size_type>(u_lower), img.width() - 1);
    const auto x_upper = std::min(static_cast<size_type>(u_upper), img.width() - 1);
    const auto y_lower = std::min(static_cast<size_type>(v_lower), img.height() - 1);
    const auto y_upper = std::min(static_cast<size_type>(v_upper), img.height() - 1);

    const auto& c0 = img(x_lower, y_lower);
    const auto& c1 = img(x_upper, y_lower);
    const auto& c2 = img(x_lower, y_upper);
    const auto& c3 = img(x_upper, y_upper);

    return v_bias * (u_bias * c0 + (1.0f - u_bias) * c1) + (1.0f - v_bias) * (u_bias * c2 + (1.0f - u_bias) * c3);
}

template <typename Remap = RemapNone>
inline RGB sample_bilinear(const Image& img, float s, float t, Remap remap = Remap{})
{
    return sample_bilinear(img, s, t, remap, remap);
}
} // namespace sp
