#pragma once

/// @author Keith Jeffery

#pragma once

#include <cstdint>

namespace sp {

inline std::uint64_t morton_encode(std::uint32_t x, std::uint32_t y) noexcept
{
    std::uint64_t a(x);
    a = (a | (a << 16ull)) & 0x0000FFFF0000FFFFull;
    a = (a | (a << 8ull)) & 0x00FF00FF00FF00FFull;
    a = (a | (a << 4ull)) & 0x0F0F0F0F0F0F0F0Full;
    a = (a | (a << 2ull)) & 0x3333333333333333ull;
    a = (a | (a << 1ull)) & 0x5555555555555555ull;

    std::uint64_t b(y);
    b = (b | (b << 16ull)) & 0x0000FFFF0000FFFFull;
    b = (b | (b << 8ull)) & 0x00FF00FF00FF00FFull;
    b = (b | (b << 4ull)) & 0x0F0F0F0F0F0F0F0Full;
    b = (b | (b << 2ull)) & 0x3333333333333333ull;
    b = (b | (b << 1ull)) & 0x5555555555555555ull;

    return a | (b << 1ull);
}

// morton_1 - extract even bits
inline uint32_t morton_decode_1(std::uint64_t a) noexcept
{
    a = a & 0x5555555555555555ull;
    a = (a | (a >> 1ull)) & 0x3333333333333333ull;
    a = (a | (a >> 2ull)) & 0x0F0F0F0F0F0F0F0Full;
    a = (a | (a >> 4ull)) & 0x00FF00FF00FF00FFull;
    a = (a | (a >> 8ull)) & 0x0000FFFF0000FFFFull;
    a = (a | (a >> 16ull)) & 0x00000000FFFFFFFFull;
    return static_cast<uint32_t>(a);
}

inline void morton_decode(std::uint64_t d, std::uint32_t& x, std::uint32_t& y) noexcept
{
    x = morton_decode_1(d);
    y = morton_decode_1(d >> 1ull);
}

inline std::uint32_t morton_encode(std::uint16_t x, std::uint16_t y) noexcept
{
    std::uint32_t a(x);
    a = (a | (a << 8ul)) & 0x00FF00FFul;
    a = (a | (a << 4ul)) & 0x0F0F0F0Ful;
    a = (a | (a << 2ul)) & 0x33333333ul;
    a = (a | (a << 1ul)) & 0x55555555ul;

    std::uint32_t b(y);
    b = (b | (b << 8ul)) & 0x00FF00FFul;
    b = (b | (b << 4ul)) & 0x0F0F0F0Ful;
    b = (b | (b << 2ul)) & 0x33333333ul;
    b = (b | (b << 1ul)) & 0x55555555ul;

    return a | (b << 1ul);
}

// morton_1 - extract even bits
inline uint16_t morton_decode_1(std::uint32_t a) noexcept
{
    a = a & 0x55555555ul;
    a = (a | (a >> 1ul)) & 0x33333333ul;
    a = (a | (a >> 2ul)) & 0x0F0F0F0Ful;
    a = (a | (a >> 4ul)) & 0x00FF00FFul;
    a = (a | (a >> 8ul)) & 0x0000FFFFul;
    return static_cast<uint16_t>(a);
}

inline std::uint16_t morton_decode_x(std::uint32_t d) noexcept
{
    return morton_decode_1(d);
}

inline std::uint16_t morton_decode_y(std::uint32_t d) noexcept
{
    return morton_decode_1(d >> 1ul);
}

inline std::pair<std::uint16_t, std::uint16_t> morton_decode(std::uint32_t d) noexcept
{
    const std::int16_t x = morton_decode_x(d);
    const std::int16_t y = morton_decode_y(d);
    return { x, y };
}

} // namespace sp
