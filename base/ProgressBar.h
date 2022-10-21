//
// Created by krjef on 10/12/2022.
//

#pragma once

/// @author Keith Jeffery

#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>

namespace sp {

class ProgressBar
{
public:
    explicit ProgressBar(unsigned total)
    : m_last_time(std::chrono::system_clock::now())
    , m_step(0)
    , m_total(total)
    {
    }

    void update(unsigned n = 1)
    {
        m_step.fetch_add(n, std::memory_order_relaxed);
    }

    void draw()
    {
        using namespace std::chrono_literals;

        const auto now = std::chrono::system_clock::now();

        // This work may all be in vain, but we want to limit our locking.
        const auto step = m_step.load(std::memory_order_relaxed);

        const unsigned percentage = calc_percentage(step, m_total);

        constexpr char     marker  = '*';
        constexpr unsigned width   = 50;
        const unsigned     numDraw = calculate_number_of_fill_elements(step, m_total, width);

        const std::string bar(numDraw, marker);

        // We force an update at 100% whether or not our time has elapsed.
        if (std::lock_guard<std::mutex> lock(m_mutex); m_step < m_total && now - m_last_time < 1.0s) {
            return;
        } else {
            m_last_time = now;
            // clang-format off
            std::cout << '\r'
                      << std::right << std::setfill(' ') << std::setw(4) << percentage << "% |"
                      << std::left  << std::setfill('-') << std::setw(width) << bar << "| "
                      << step << '/' << m_total
                      << std::flush;
            // clang-format on
        }
    }

private:
    static unsigned calc_percentage(unsigned step, unsigned total) noexcept
    {
        return (total > 100) ? step / (total / 100) : step * 100 / total;
    }

    static unsigned calculate_number_of_fill_elements(unsigned step, unsigned total, unsigned width) noexcept
    {
        return (total > width) ? step / (total / width) : step * width / total;
    }

    using TimePoint = std::chrono::system_clock::time_point;

    mutable std::mutex    m_mutex;
    TimePoint             m_last_time;
    std::atomic<unsigned> m_step;
    unsigned              m_total;
};

} // namespace sp
