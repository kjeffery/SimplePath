#pragma once

/// @author Keith Jeffery

#include "../base/Array2D.h"
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

inline float rgb_to_srgb(float u) noexcept
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

inline RGB sample_nearest_neighbor(const Image& img, float s, float t)
{
    assert(s >= 0.0f);
    assert(t >= 0.0f);
    assert(s < 1.0f);
    assert(t < 1.0f);

    using size_type = Image::size_type;

    const float u = std::round(s * img.width());
    const float v = std::round(t * img.height());

    const size_type x = std::min(static_cast<size_type>(u), img.width() - 1);
    const size_type y = std::min(static_cast<size_type>(v), img.height() - 1);

    return img(x, y);
}

inline RGB sample_bilinear(const Image& img, float s, float t)
{
    assert(s >= 0.0f);
    assert(t >= 0.0f);
    assert(s < 1.0f);
    assert(t < 1.0f);

    using size_type = Image::size_type;

    const float u       = s * img.width();
    const float v       = t * img.height();
    const float u_lower = std::floor(u);
    const float u_upper = std::ceil(u);
    const float v_lower = std::floor(v);
    const float v_upper = std::ceil(v);

    const float u_bias = u_upper - u_lower;
    const float v_bias = v_upper - v_lower;

    const size_type x_lower = std::min(static_cast<size_type>(u_lower), img.width() - 1);
    const size_type x_upper = std::min(static_cast<size_type>(u_upper), img.width() - 1);
    const size_type y_lower = std::min(static_cast<size_type>(v_lower), img.height() - 1);
    const size_type y_upper = std::min(static_cast<size_type>(v_upper), img.height() - 1);

    const auto& c0 = img(x_lower, y_lower);
    const auto& c1 = img(x_upper, y_lower);
    const auto& c2 = img(x_lower, y_upper);
    const auto& c3 = img(x_upper, y_upper);

    return v_bias * (u_bias * c0 + (1.0f - u_bias) * c1) + (1.0f - v_bias) * (u_bias * c2 + (1.0f - u_bias) * c3);
}
} // namespace sp
