
#pragma once

#include "Math.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <execution>
#include <span>
#include <vector>

// Inspired by PBRT
class Distribution1D
{
public:
    Distribution1D(std::vector<float> f, float min, float max)
    : m_min(min)
    , m_max(max)
    , m_function{ std::move(f) }
    , m_cdf(m_function.size() + std::size_t{ 1 })
    {
        assert(max > min);

        std::transform(std::execution::par_unseq, m_function.cbegin(), m_function.cend(), m_function.begin(),
                       [](const auto& x) { return std::abs(x); });

        // Compute integral of step function at $x_i$
        m_cdf[0] = 0;
        for (std::size_t i = 1; i < m_function.size() + 1; ++i) {
            assert(m_function[i - 1] >= 0.0f);
            m_cdf[i] = m_cdf[i - 1] + m_function[i - 1] * (max - min) / static_cast<float>(m_function.size());
        }

        // Transform step function integral into CDF
        m_function_integral = m_cdf.back();
        if (m_function_integral == 0) {
            for (std::size_t i = 1; i < m_function.size() + 1; ++i) {
                m_cdf[i] = static_cast<float>(i) / static_cast<float>(m_function.size());
            }
        } else {
            std::transform(std::execution::par_unseq, std::next(m_cdf.cbegin()), m_cdf.cend(), m_cdf.begin(),
                           [this](auto x) { return x / m_function_integral; });
        }
    }

    explicit Distribution1D(std::vector<float> f)
    : Distribution1D{ std::move(f), 0.0f, 1.0f }
    {
    }

    Distribution1D(std::span<float> f, float min, float max)
    : Distribution1D{ std::vector<float>{ f.begin(), f.end() }, min, max }
    {
    }

    explicit Distribution1D(std::span<float> f)
    : Distribution1D{ f, 0.0f, 1.0f }
    {
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t
    {
        return m_function.size();
    }

    [[nodiscard]] auto integral() const noexcept -> float
    {
        return m_function_integral;
    }

    [[nodiscard]] auto eval(const std::size_t idx) const noexcept -> float
    {
        return m_function[idx];
    }

    [[nodiscard]] auto sample_continuous(const float u, float* pdf, std::size_t* off = nullptr) const -> float
    {
        // Find surrounding CDF segments and _offset_
        const auto offset = get_offset(u);
        if (off) {
            *off = offset;
        }

        // Compute offset along CDF segment
        auto du = u - m_cdf[offset];
        if ((m_cdf[offset + 1] - m_cdf[offset]) > 0) {
            du /= (m_cdf[offset + 1] - m_cdf[offset]);
        }
        assert(!std::isnan(du));

        // Compute PDF for sampled offset
        if (pdf) {
            *pdf = (m_function_integral > 0) ? m_function[offset] / m_function_integral : 0.0f;
        }

        return sp::lerp((static_cast<float>(offset) + du) / static_cast<float>(size()), m_min, m_max);
    }

    [[nodiscard]] auto sample_discrete(const float u, float* const pdf = nullptr, float* const u_remapped = nullptr) const -> std::size_t
    {
        const auto offset = get_offset(u);
        if (pdf) {
            *pdf = (m_function_integral > 0.0f) ? m_function[offset] / (m_function_integral * static_cast<float>(size())) : 0.0f;
        }
        if (u_remapped) {
            *u_remapped = (u - m_cdf[offset]) / (m_cdf[offset + 1] - m_cdf[offset]);
            assert(*u_remapped >= 0.f && *u_remapped <= 1.f);
        }
        return offset;
    }

    [[nodiscard]] auto discrete_pdf(int index) const -> float
    {
        assert(index >= 0 && index < size());
        return m_function[index] / (m_function_integral * static_cast<float>(size()));
    }

    [[nodiscard]] std::optional<float> invert(const float x) const
    {
        // Compute offset to CDF values that bracket $x$
        if (x < m_min || x > m_max) {
            return {};
        }
        const auto c      = (x - m_min) / (m_max - m_min) * static_cast<float>(m_function.size());
        const auto offset = std::clamp(static_cast<std::size_t>(c), std::size_t{ 0 }, m_function.size() - std::size_t{ 1 });
        assert(offset >= 0 && offset + 1 < m_cdf.size());

        // Linearly interpolate between adjacent CDF values to find sample value
        const auto delta = c - static_cast<float>(offset);
        return sp::lerp(delta, m_cdf[offset], m_cdf[offset + 1]);
    }

private:
    [[nodiscard]] auto get_offset(const float u) const noexcept -> std::size_t
    {
        const auto it = std::ranges::upper_bound(m_cdf, u);
        if (it == m_cdf.cend() || it == std::prev(m_cdf.cend())) {
            // We have to leave room for adding one to the returned index
            return m_cdf.size() - 2;
        }
        return std::distance(m_cdf.cbegin(), it);
    }

    // Distribution1D Public Data
    float              m_function_integral{ 0 };
    float              m_min;
    float              m_max;
    std::vector<float> m_function;
    std::vector<float> m_cdf;
};
