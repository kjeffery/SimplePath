#pragma once

/// @author Keith Jeffery

#include "Vector2.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <random>

namespace sp {
enum Seed : std::uint32_t;

template <unsigned dimension>
class RSequence
{
    static auto phi() noexcept -> float
    {
        auto x = 2.0f;

        for (auto i = 0; i < 10; ++i) {
            x = std::pow(1.0f + x, 1.0f / (static_cast<float>(dimension) + 1.0f));
        }

        return x;
    }

    static auto mod1(const float f) noexcept -> float
    {
        float dummy;
        return std::modf(f, &dummy);
    }

    auto r_sequence(Seed seed, std::uint32_t n) noexcept -> std::array<float, dimension>
    {
        const auto fseed = static_cast<float>(seed) / std::numeric_limits<float>::max();

        std::array<float, dimension> result;
        for (unsigned i = 0; i < dimension; ++i) {
            result[i] = mod1(fseed + m_alpha[i] * (static_cast<float>(n) + 1.0f));
        }
        return result;
    }

public:
    RSequence()
    {
        for (unsigned i = 0; i < dimension; ++i) {
            m_alpha[i] = mod1(std::pow(1.0f / m_g, i + 1.0f));
        }
    }

    auto operator()(Seed seed, std::uint32_t n) noexcept -> std::array<float, dimension>
    {
        return r_sequence(seed, n);
    }

private:
    float                        m_g{ phi() };
    std::array<float, dimension> m_alpha{};
};

template <typename RNG>
auto canonical(RNG& rng) -> float
{
    static std::uniform_real_distribution<float> dist;

    auto result = dist(rng);
    while (result >= 1.0f) {
        result = dist(rng);
    }
    return result;
}

class Sampler
{
public:
    virtual ~Sampler() = default;

    [[nodiscard]] auto get_next_1D() noexcept -> float
    {
        return get_next_1D_impl();
    }

    [[nodiscard]] auto get_next_2D() noexcept -> Point2
    {
        return get_next_2D_impl();
    }

private:
    [[nodiscard]] virtual auto get_next_1D_impl() noexcept -> float = 0;
    [[nodiscard]] virtual auto get_next_2D_impl() noexcept -> Point2 = 0;
};

class IncoherentSampler : public Sampler
{
public:
    [[nodiscard]] static auto create_new_set(const Seed seed, std::uint32_t /*num*/) noexcept -> IncoherentSampler
    {
        return IncoherentSampler{ seed };
    }

    [[nodiscard]] static auto create_new_sequence(const Seed seed) noexcept -> IncoherentSampler
    {
        return IncoherentSampler{ seed };
    }

private:
    [[nodiscard]] auto get_next_1D_impl() noexcept -> float override
    {
        return canonical();
    }

    [[nodiscard]] auto get_next_2D_impl() noexcept -> Point2 override
    {
        return { canonical(), canonical() };
    }

    explicit IncoherentSampler(Seed seed) noexcept
    : m_rng(seed)
    {
    }

    auto canonical() -> float
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

class RSequenceSampler : public Sampler
{
public:
    [[nodiscard]] static auto create_new_set(const Seed seed, std::uint32_t num) noexcept -> RSequenceSampler
    {
        return RSequenceSampler{ seed };
    }

    [[nodiscard]] static auto create_new_sequence(Seed seed) noexcept -> RSequenceSampler
    {
        return RSequenceSampler{ seed };
    }

private:
    [[nodiscard]] auto get_next_1D_impl() noexcept -> float override
    {
        const auto [p0] = m_1D(m_seed_2D, m_count_2D++);
        return p0;
    }

    [[nodiscard]] auto get_next_2D_impl() noexcept -> Point2 override
    {
        const auto [p0, p1] = m_2D(m_seed_2D, m_count_2D++);
        return { p0, p1 };
    }

    explicit RSequenceSampler(Seed seed) noexcept
    : m_seed_1D{ seed }
    , m_seed_2D{ sp::Seed{ seed ^ 0x6184faf4 } }
    {
    }

    sp::Seed m_seed_1D;
    sp::Seed m_seed_2D;

    std::uint32_t m_count_1D{ 0 };
    std::uint32_t m_count_2D{ 0 };

    RSequence<1> m_1D;
    RSequence<2> m_2D;
};
} // namespace sp
