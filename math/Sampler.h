#pragma once

/// @author Keith Jeffery

#include "Vector2.h"

#include <cstdint>
#include <random>

namespace sp {

class Sampler
{
public:
    [[nodiscard]] static Sampler create_new_set(std::uint32_t seed, std::uint32_t num) noexcept
    {
        return Sampler{ seed };
    }

    [[nodiscard]] static Sampler create_new_sequence(std::uint32_t seed) noexcept
    {
        return Sampler{ seed };
    }

    [[nodiscard]] float get_next_1D() noexcept
    {
        return canonical();
    }

    [[nodiscard]] Point2 get_next_2D() noexcept
    {
        return { canonical(), canonical() };
    }

private:
    explicit Sampler(std::uint32_t seed) noexcept
    : m_rng(seed)
    {
    }

    float canonical()
    {
        float result;
        do {
            result = m_dist(m_rng);
        } while (result >= 1.0f);
        return result;
    }

    std::mt19937_64                       m_rng;
    std::uniform_real_distribution<float> m_dist;
};

} // namespace sp
