
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

class Distribution1D
{
public:
    explicit Distribution1D(std::vector<float> f)
    {
    }

    explicit Distribution1D(std::span<float> f)
    : m_function{ f.begin(), f.end() }
    , m_cdf(f.size() + std::size_t{ 1 })
    {
        // Compute integral of step function at $x_i$
        m_cdf[0] = 0;
        for (std::size_t i = 1; i < f.size() + 1; ++i) {
            m_cdf[i] = m_cdf[i - 1] + m_function[i - 1] / static_cast<float>(n);
        }

        // Transform step function integral into CDF
        m_function_integral = m_cdf.back();
        if (m_function_integral == 0) {
            for (std::size_t i = 1; i < f.size() + 1; ++i) {
                m_cdf[i] = static_cast<float>(i) / static_cast<float>(n);
            }
        } else {
            std::ranges::transform(m_cdf, m_cdf.begin(), [this](auto x) { return x / m_function_integral; });
            for (std::size_t i = 1; i < f.size() + 1; ++i) {
                m_cdf[i] /= m_function_integral;
            }
        }
    }

    [[nodiscard]] auto size() const -> std::size_t
    {
        return m_function.size();
    }

    auto sample_continuous(float u, float* pdf, std::size_t* off = nullptr) const -> float
    {
        // Find surrounding CDF segments and _offset_
        const auto offset = get_offset(u);
        if (off) {
            *off = offset;
        }

        // Compute offset along CDF segment
        auto du = u - m_cdf[offset];
        if ((m_cdf[offset + 1] - m_cdf[offset]) > 0) {
            assert(m_cdf[offset + 1] > m_cdf[offset]);
            du /= (m_cdf[offset + 1] - m_cdf[offset]);
        }
        assert(!std::isnan(du));

        // Compute PDF for sampled offset
        if (pdf) {
            *pdf = (m_function_integral > 0) ? m_function[offset] / m_function_integral : 0;
        }

        // Return $x\in{}[0,1)$ corresponding to sample
        return (static_cast<float>(offset) + du) / static_cast<float>(size());
    }

    auto sample_discrete(const float u, float* const pdf = nullptr, float* const u_remapped = nullptr) const -> std::size_t
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

private:
    [[nodiscard]] auto get_offset(float u) const noexcept -> std::size_t
    {
        const auto it = std::ranges::upper_bound(m_cdf, u);
        if (it == m_cdf.cend() || it == std::prev(m_cdf.cend())) {
            // We have to leave room for adding one to the returned index
            return m_cdf.size() - 2;
        }
        return std::distance(m_cdf.cbegin(), it);
    }

    // Distribution1D Public Data
    std::vector<float> m_function;
    std::vector<float> m_cdf;
    float              m_function_integral;
};
