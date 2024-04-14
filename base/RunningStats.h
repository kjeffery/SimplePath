#pragma once

#include <cmath>

/// @author Keith Jeffery

namespace sp {

// This is based on the simpler version of John D. Cook's running variance code, which is, in turn, based on Knuth.
// https://www.johndcook.com/blog/standard_deviation/
template <typename T>
class RunningStats
{
public:
    RunningStats() noexcept
    : m_n(0)
    {
    }

    void clear() noexcept
    {
        m_n = 0;
    }

    void push(T x) noexcept
    {
        ++m_n;

        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1) {
            m_old_m = m_new_m = x;
            m_old_s           = T{};
        } else {
            m_new_m = m_old_m + (x - m_old_m) / m_n;
            m_new_s = m_old_s + (x - m_old_m) * (x - m_new_m);

            // set up for next iteration
            m_old_m = m_new_m;
            m_old_s = m_new_s;
        }
    }

    [[nodiscard]] int size() const noexcept
    {
        return m_n;
    }

    [[nodiscard]] T mean() const noexcept
    {
        return (m_n > 0) ? m_new_m : T{};
    }

    [[nodiscard]] T variance() const noexcept
    {
        return ((m_n > 1) ? m_new_s / (m_n - 1) : T{});
    }

    [[nodiscard]] T standardDeviation() const noexcept
    {
        return std::sqrt(variance());
    }

private:
    int m_n;
    T   m_old_m;
    T   m_new_m;
    T   m_old_s;
    T   m_new_s;
};

} // namespace sp
