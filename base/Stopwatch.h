#pragma once

/// @author Keith Jeffery

#include <chrono>
#include <iomanip>
#include <ostream>

namespace sp {

class Stopwatch
{
public:
    Stopwatch()
    : m_start_time{ clock_type::now() }
    , m_end_time{}
    {
    }

    explicit Stopwatch(NoInitType)
    : m_start_time{ }
    , m_end_time{}
    {
    }

    void start()
    {
        m_start_time = clock_type::now();
    }

    void stop()
    {
        m_end_time = clock_type::now();
    }

    void print(std::ostream& outs) const
    {
        using namespace std::chrono_literals;
        const auto diff                                          = m_end_time - m_start_time;
        const auto [days, hours, minutes, seconds, centiseconds] = split_to_units(diff);
        outs.fill('0');
        // clang-format off
        outs << std::right;
        outs << std::setw(2) << hours << ':'
             << std::setw(2) << minutes << ':'
             << std::setw(2) << seconds << '.'
             << std::setw(2) << centiseconds;
        // clang-format on
    }

private:
    template <typename Rep, typename Period>
    static auto split_to_units(std::chrono::duration<Rep, Period> duration)
    {
        using centiseconds_type = std::chrono::duration<std::uint64_t, std::ratio<1, 100>>;

        const auto days = std::chrono::duration_cast<std::chrono::days>(duration);
        duration -= days;
        const auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
        duration -= hours;
        const auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
        duration -= minutes;
        const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
        duration -= seconds;
        const auto centiseconds = std::chrono::duration_cast<centiseconds_type>(duration);
        return std::make_tuple(days.count(), hours.count(), minutes.count(), seconds.count(), centiseconds.count());
    }

    using clock_type = std::chrono::system_clock;
    using time_point = clock_type::time_point;
    using duration   = clock_type::duration;

    time_point m_start_time;
    time_point m_end_time;
};

} // namespace sp
