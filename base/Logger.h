#pragma once

/// @author Keith Jeffery

#include <cstdlib>
#include <iostream>
#include <mutex>
#include <type_traits>

namespace sp {

class Logger
{
public:
    static Logger& instance()
    {
        static Logger the_logger;
        return the_logger;
    }

    enum class LoggingLevel
    {
        error,
        warning,
        info,
        debug
    };

    static void set_level(LoggingLevel level)
    {
        s_log_level = level;
    }

    static bool is_enabled(LoggingLevel level)
    {
        using underlying_type = std::underlying_type_t<LoggingLevel>;
        return static_cast<underlying_type>(s_log_level) >= static_cast<underlying_type>(level);
    }

    template <typename... Args>
    void log_fatal(const char* const file, int line, Args&&... args)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        std::cerr << "Fatal error [" << file << ':' << line << "]: ";
        (std::cerr << ... << args) << '\n';
        std::exit(EXIT_FAILURE);
    }

    template <typename... Args>
    void log_error(const char* const file, int line, Args&&... args)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        std::cerr << "Error [" << file << ':' << line << "]: ";
        (std::cerr << ... << args) << '\n';
    }

    template <typename... Args>
    void log_warning(const char* const file, int line, Args&&... args)
    {
        if (is_enabled(LoggingLevel::warning)) {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::cout << "Warning [" << file << ':' << line << "]: ";
            (std::cout << ... << args) << '\n';
        }
    }

    template <typename... Args>
    void log_info(const char* const file, int line, Args&&... args)
    {
        if (is_enabled(LoggingLevel::info)) {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::cout << "Info [" << file << ':' << line << "]: ";
            (std::cout << ... << args) << '\n';
        }
    }

    template <typename... Args>
    void log_debug(const char* const file, int line, Args&&... args)
    {
        if (is_enabled(LoggingLevel::debug)) {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::cout << "Debug [" << file << ':' << line << "]: ";
            (std::cout << ... << args) << '\n';
        }
    }

private:
    Logger() = default;

    static LoggingLevel s_log_level;
    mutable std::mutex  m_mutex;
};

} // namespace sp

#define LOG_FATAL(...)   sp::Logger::instance().log_fatal(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_ERROR(...)   sp::Logger::instance().log_error(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_WARNING(...) sp::Logger::instance().log_warning(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_INFO(...)    sp::Logger::instance().log_info(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_DEBUG(...)   sp::Logger::instance().log_debug(__FILE__, __LINE__, __VA_ARGS__)
