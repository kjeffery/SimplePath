//
// Created by krjef on 9/22/2022.
//

#pragma once

/// @author Keith Jeffery

#include <bit>
#include <cstdint>

static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
              "We don't support mixed endianess");

#if defined(_MSC_FULL_VER)

#include <cstdlib>

inline std::uint16_t little_endian(std::uint16_t v) noexcept
{
    using enum std::endian;
    if constexpr (native == little) {
        return v;
    } else if constexpr (native == big) {
        return _byteswap_ushort(v);
    }
}

inline std::uint16_t big_endian(std::uint16_t v) noexcept
{
    using enum std::endian;
    if constexpr (native == little) {
        return _byteswap_ushort(v);
    } else if constexpr (native == big) {
        return v;
    }
}

inline std::uint32_t little_endian(std::uint32_t v) noexcept
{
    using enum std::endian;
    if constexpr (native == little) {
        return v;
    } else if constexpr (native == big) {
        return _byteswap_ulong(v);
    }
}

inline std::uint32_t big_endian(std::uint32_t v) noexcept
{
    using enum std::endian;
    if constexpr (native == little) {
        return _byteswap_ulong(v);
    } else if constexpr (native == big) {
        return v;
    }
}

inline std::uint64_t little_endian(std::uint64_t v) noexcept
{
    using enum std::endian;
    if constexpr (native == little) {
        return v;
    } else if constexpr (native == big) {
        return _byteswap_uint64(v);
    }
}

inline std::uint64_t big_endian(std::uint64_t v) noexcept
{
    using enum std::endian;
    if constexpr (native == little) {
        return _byteswap_uint64(v);
    } else if constexpr (native == big) {
        return v;
    }
}

#else

#include <byteswap.h>

inline std::uint16_t little_endian(std::uint16_t v) noexcept
{
    if constexpr (std::endian::native == std::endian::little) {
        return v;
    } else if constexpr (std::endian::native == std::endian::big) {
        return __bswap_16(v);
    }
}

inline std::uint16_t big_endian(std::uint16_t v) noexcept
{
    if constexpr (std::endian::native == std::endian::little) {
        return __bswap_16(v);
    } else if constexpr (std::endian::native == std::endian::big) {
        return v;
    }
}

inline std::uint32_t little_endian(std::uint32_t v) noexcept
{
    if constexpr (std::endian::native == std::endian::little) {
        return v;
    } else if constexpr (std::endian::native == std::endian::big) {
        return __bswap_32(v);
    }
}

inline std::uint32_t big_endian(std::uint32_t v) noexcept
{
    if constexpr (std::endian::native == std::endian::little) {
        return __bswap_32(v);
    } else if constexpr (std::endian::native == std::endian::big) {
        return v;
    }
}

inline std::uint64_t little_endian(std::uint64_t v) noexcept
{
    if constexpr (std::endian::native == std::endian::little) {
        return v;
    } else if constexpr (std::endian::native == std::endian::big) {
        return __bswap_64(v);
    }
}

inline std::uint64_t big_endian(std::uint64_t v) noexcept
{
    if constexpr (std::endian::native == std::endian::little) {
        return __bswap_64(v);
    } else if constexpr (std::endian::native == std::endian::big) {
        return v;
    }
}

#endif

inline std::uint8_t little_endian(std::uint8_t v) noexcept
{
    return v;
}

inline std::uint8_t big_endian(std::uint8_t v) noexcept
{
    return v;
}

inline std::int8_t little_endian(std::int8_t v) noexcept
{
    return v;
}

inline std::int8_t big_endian(std::int8_t v) noexcept
{
    return v;
}

inline std::int16_t little_endian(std::int16_t v) noexcept
{
    const auto bv = little_endian(static_cast<std::uint16_t>(v));
    return static_cast<const std::int16_t>(bv);
}

inline std::int16_t big_endian(std::int16_t v) noexcept
{
    const auto bv = big_endian(static_cast<std::uint16_t>(v));
    return static_cast<const std::int16_t>(bv);
}

inline std::int32_t little_endian(std::int32_t v) noexcept
{
    const auto bv = little_endian(static_cast<std::uint32_t>(v));
    return static_cast<const std::int32_t>(bv);
}

inline std::int32_t big_endian(std::int32_t v) noexcept
{
    const auto bv = big_endian(static_cast<std::uint32_t>(v));
    return static_cast<const std::int32_t>(bv);
}

inline std::int64_t little_endian(std::int64_t v) noexcept
{
    const auto bv = little_endian(static_cast<std::uint64_t>(v));
    return static_cast<const std::int64_t>(bv);
}

inline std::int64_t big_endian(std::int64_t v) noexcept
{
    const auto bv = big_endian(static_cast<std::uint64_t>(v));
    return static_cast<const std::int64_t>(bv);
}

#if defined(__cpp_lib_bit_cast)
inline float little_endian(float v) noexcept
{
    const auto bv = little_endian(std::bit_cast<std::uint32_t>(v));
    return std::bit_cast<float>(bv);
}

inline float big_endian(float v) noexcept
{
    const auto bv = big_endian(std::bit_cast<std::uint32_t>(v));
    return std::bit_cast<float>(bv);
}

inline double little_endian(double v) noexcept
{
    const auto bv = little_endian(std::bit_cast<std::uint64_t>(v));
    return std::bit_cast<double>(bv);
}

inline double big_endian(double v) noexcept
{
    const auto bv = big_endian(std::bit_cast<std::uint64_t>(v));
    return std::bit_cast<double>(bv);
}
#else
inline float little_endian(float v) noexcept
{
    const auto bv = little_endian(*reinterpret_cast<std::uint32_t*>(&v));
    return *reinterpret_cast<const float*>(&bv);
}

inline float big_endian(float v) noexcept
{
    const auto bv = big_endian(*reinterpret_cast<std::uint32_t*>(&v));
    return *reinterpret_cast<const float*>(&bv);
}

inline double little_endian(double v) noexcept
{
    const auto bv = little_endian(*reinterpret_cast<std::uint64_t*>(&v));
    return *reinterpret_cast<const double*>(&bv);
}

inline double big_endian(double v) noexcept
{
    const auto bv = big_endian(*reinterpret_cast<std::uint64_t*>(&v));
    return *reinterpret_cast<const double*>(&bv);
}
#endif
