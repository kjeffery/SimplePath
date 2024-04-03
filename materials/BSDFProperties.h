#pragma once
#include <type_traits>

enum class BSDFProperties
{
    // clang-format off
    none         = 0b000'000,
    diffuse      = 0b000'001,
    glossy       = 0b000'010,
    specular     = 0b000'100,
    reflective   = 0b001'000,
    transmissive = 0b010'000,
    // clang-format on
};

constexpr auto operator&=(BSDFProperties& a, const BSDFProperties b) noexcept -> BSDFProperties& {
    using underlying_type = std::underlying_type_t<BSDFProperties>;
    a = BSDFProperties{static_cast<underlying_type>(a) & static_cast<underlying_type>(b)};
    return a;
}

constexpr auto operator|=(BSDFProperties& a, const BSDFProperties b) noexcept -> BSDFProperties& {
    using underlying_type = std::underlying_type_t<BSDFProperties>;
    a = BSDFProperties{static_cast<underlying_type>(a) | static_cast<underlying_type>(b)};
    return a;
}

constexpr auto operator&(BSDFProperties a, const BSDFProperties b) noexcept -> BSDFProperties {
    a &= b;
    return a;
}

constexpr auto operator|(BSDFProperties a, const BSDFProperties b) noexcept -> BSDFProperties {
    a |= b;
    return a;
}

constexpr auto is_diffuse(const BSDFProperties properties) noexcept -> bool {
    return (properties & BSDFProperties::diffuse) == BSDFProperties::diffuse;
}

constexpr auto is_glossy(const BSDFProperties properties) noexcept -> bool {
    return (properties & BSDFProperties::glossy) == BSDFProperties::glossy;
}

constexpr auto is_specular(const BSDFProperties properties) noexcept -> bool {
    return (properties & BSDFProperties::specular) == BSDFProperties::specular;
}

constexpr auto is_reflective(const BSDFProperties properties) noexcept -> bool {
    return (properties & BSDFProperties::reflective) == BSDFProperties::reflective;
}

constexpr auto is_transmissive(const BSDFProperties properties) noexcept -> bool {
    return (properties & BSDFProperties::transmissive) == BSDFProperties::transmissive;
}
