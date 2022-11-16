#pragma once

/// @author Keith Jeffery

#include <cassert>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>

#if __has_include(<syncstream>)
#include <syncstream>
#define OSYNC(s) std::osyncstream(s)
#else
#define OSYNC(s) s
#endif

namespace sp {

enum class StreamType
{
    cout,
    cerr
};

using namespace std::chrono_literals;

class AccumulatedLoggerImpl
{
public:
    AccumulatedLoggerImpl()
    : m_thread(&AccumulatedLoggerImpl::do_work, this)
    {
        assert(m_thread.joinable());
    }

    ~AccumulatedLoggerImpl()
    {
        try {
            shutdown();
        } catch (...) {
        }
    }

    static AccumulatedLoggerImpl& instance() noexcept
    {
        static AccumulatedLoggerImpl logger;
        return logger;
    }

    void queue_message(std::string message, StreamType stream_type)
    {
        std::lock_guard<std::mutex> lock(m_mutex);

        // Lookup message
        auto it = m_messages.find(message);
        if (it == m_messages.cend()) {
            m_messages.emplace_hint(it, std::move(message), MessageData{ stream_type, 1u });
        } else {
            ++it->second.m_count;
        }
    }

    void shutdown() noexcept
    {
        set_shutdown();
        m_condition.notify_one();

        assert(m_thread.joinable());
        m_thread.join();

        flush_impl();
    }

    void flush()
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        flush_impl();
    }

private:
    void set_shutdown() noexcept
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_shutdown = true;
    }

    void do_work()
    {
        using Clock     = std::chrono::system_clock;
        using TimePoint = std::chrono::time_point<Clock>;

        TimePoint last_flush = Clock::now();

        std::unique_lock<std::mutex> lock(m_mutex);
        while (!m_shutdown) {
            const auto duration = last_flush + s_duration;
            if (const auto cv_status = m_condition.wait_until(lock, duration); cv_status == std::cv_status::timeout) {
                flush_impl();
                last_flush = Clock::now();
            }
        }
    }

    void flush_impl()
    {
        for (const auto& [message, message_data] : m_messages) {
            std::ostream& outs = (message_data.m_stream_type == StreamType::cout) ? std::cout : std::cerr;
            assert(message_data.m_count > 0);

            // We're under a locked mutex at this point, but it doesn't stop other threads from writing directly to the
            // stream. Wrap it in std::osyncstream.
            if (message_data.m_count == 1) {
                OSYNC(outs) << message << '\n';
            } else {
                OSYNC(outs) << message << " (occurred " << message_data.m_count << " times)\n";
            }
        }
        m_messages.clear();
    }

    struct MessageData
    {
        // We could allow a pointer to a stream, but since we're non-deterministic, we would have to worry about
        // lifetime. This could be a shared pointer, but that complicates the creation. For now, I'm going to stick
        // to cout and cerr.
        StreamType m_stream_type;
        unsigned   m_count;
    };

    static constexpr auto s_duration = 750ms;

    bool                                         m_shutdown = false;
    std::unordered_map<std::string, MessageData> m_messages;
    mutable std::mutex                           m_mutex;
    std::thread                                  m_thread;
    std::condition_variable                      m_condition;
};

// Accumulated
class AccumulatedLogger
{
public:
    static void log(std::string message, StreamType streamType)
    {
        AccumulatedLoggerImpl::instance().queue_message(std::move(message), streamType);
    }

    static void flush()
    {
        AccumulatedLoggerImpl::instance().flush();
    }
};

} // namespace sp
