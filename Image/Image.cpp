/// @author Keith Jeffery

#include "Image.h"
#include "../base/Endian.h"
#include "../math/RGB.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace sp {
void write_ppm(std::ostream& outs, const Image& img)
{
    const int nx = img.width();
    const int ny = img.height();
    outs << "P3\n" << nx << ' ' << ny << "\n255\n";
    for (int j = ny - 1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            const RGB c = rgb_to_srgb(img(i, j));
            // const RGB& c = fb[pixel_index];
            const int ir = static_cast<int>(255.99f * c.r);
            const int ig = static_cast<int>(255.99f * c.g);
            const int ib = static_cast<int>(255.99f * c.b);
            outs << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
}

void write_ppm(const std::filesystem::path& file, const Image& img)
{
    std::ofstream outs(file, std::ios_base::binary | std::ios_base::out);
    if (!outs) {
        throw ImageError("Unable to open " + file.string());
    }
    write_ppm(outs, img);
}

void write_pfm(std::ostream& outs, const Image& img)
{
    constexpr int byte_order = (std::endian::native == std::endian::little) ? -1 : +1;

    const int nx = img.width();
    const int ny = img.height();
    outs << "PF\n" << nx << ' ' << ny << '\n' << byte_order << '\n';
    for (int j = ny - 1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            const auto& c = img(i, j);
            outs.write(reinterpret_cast<const char*>(&c.r), sizeof(float));
            outs.write(reinterpret_cast<const char*>(&c.g), sizeof(float));
            outs.write(reinterpret_cast<const char*>(&c.b), sizeof(float));
        }
    }
}

void write_pfm(const std::filesystem::path& file, const Image& img)
{
    std::ofstream outs(file, std::ios_base::binary | std::ios_base::out);
    if (!outs) {
        throw ImageError("Unable to open " + file.string());
    }
    write_pfm(outs, img);
}

void write(const std::filesystem::path& file, const Image& img)
{
    const auto& ext = file.extension().string();
    if (ext == ".pfm") {
        write_pfm(file, img);
    } else if (ext == ".ppm") {
        write_ppm(file, img);
    } else {
        throw ImageError("Unknown file extension: " + ext);
    }
}

Image read_pfm(std::istream& ins)
{
    std::string format;
    std::getline(ins, format);
    if (format != "PF") {
        throw ImageError("Unexpected format");
    }
    int nx;
    int ny;
    ins >> nx;
    ins >> ny;
    float byte_order;
    ins >> byte_order;
    ins.get(); // Get last '\n'

    using ConvertFunction         = std::uint32_t (*)(std::uint32_t);
    const ConvertFunction be      = &big_endian_to_native;
    const ConvertFunction le      = &little_endian_to_native;
    const ConvertFunction convert = (byte_order > 0) ? be : le;

    Image img(nx, ny);

    for (int j = ny - 1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            std::uint32_t r;
            std::uint32_t g;
            std::uint32_t b;
            ins.read(reinterpret_cast<char*>(&r), sizeof(std::uint32_t));
            ins.read(reinterpret_cast<char*>(&g), sizeof(std::uint32_t));
            ins.read(reinterpret_cast<char*>(&b), sizeof(std::uint32_t));
            r              = convert(r);
            g              = convert(g);
            b              = convert(b);
            const float fr = *reinterpret_cast<float*>(&r);
            const float fg = *reinterpret_cast<float*>(&g);
            const float fb = *reinterpret_cast<float*>(&b);
            img(i, j)      = RGB(fr, fg, fb);
        }
    }

    return img;
}

Image read_pfm(const std::filesystem::path& file)
{
    std::ifstream ins(file, std::ios_base::binary | std::ios_base::in);
    if (!ins) {
        throw ImageError("Unable to open " + file.string());
    }
    return read_pfm(ins);
}

Image read(const std::filesystem::path& file)
{
    const auto& ext = file.extension().string();
    if (ext == ".pfm") {
        return read_pfm(file);
    } else {
        throw ImageError("Unknown file extension: " + ext);
    }
}

Image& operator*=(Image& img, const RGB& c) noexcept
{
    using size_type = Image::size_type;

    for (size_type y = 0; y < img.height(); ++y) {
        for (size_type x = 0; x < img.width(); ++x) {
            img(x, y) *= c;
        }
    }
    return img;
}

Image operator*(Image img, const RGB& c) noexcept
{
    img *= c;
    return img;
}

Image operator*(const RGB& c, const Image& img) noexcept
{
    return img * c;
}
} // namespace sp
