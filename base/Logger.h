#pragma once

/// @author Keith Jeffery

#include "AccumulatedLogger.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
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

    void flush()
    {
        AccumulatedLogger::flush();
    }

    template <typename... Args>
    void log_fatal(const char* const file, int line, Args&&... args)
    {
        // We don't accumulate fatal errors: write immediately.
        auto outs = OSYNC(std::cerr);
        outs << "Fatal error [" << file << ':' << line << "]: ";
        (outs << ... << args);
        std::exit(EXIT_FAILURE);
    }

    template <typename... Args>
    void log_error(const char* const file, int line, Args&&... args)
    {
        std::ostringstream outs;
        outs.copyfmt(std::cerr);
        outs << "Error [" << file << ':' << line << "]: ";
        (outs << ... << args);
        AccumulatedLogger::log(outs.str(), StreamType::cerr);
    }

    template <typename... Args>
    void log_warning(const char* const file, int line, Args&&... args)
    {
        if (is_enabled(LoggingLevel::warning)) {
            std::ostringstream outs;
            outs.copyfmt(std::cout);
            outs << "Warning [" << file << ':' << line << "]: ";
            (outs << ... << args);
            AccumulatedLogger::log(outs.str(), StreamType::cout);
        }
    }

    template <typename... Args>
    void log_info(const char* const file, int line, Args&&... args)
    {
        if (is_enabled(LoggingLevel::info)) {
            std::ostringstream outs;
            outs.copyfmt(std::cout);
            outs << "Info [" << file << ':' << line << "]: ";
            (outs << ... << args);
            AccumulatedLogger::log(outs.str(), StreamType::cout);
        }
    }

    template <typename... Args>
    void log_debug(const char* const file, int line, Args&&... args)
    {
        if (is_enabled(LoggingLevel::debug)) {
            std::ostringstream outs;
            outs.copyfmt(std::cout);
            outs << "Debug [" << file << ':' << line << "]: ";
            (outs << ... << args);
            AccumulatedLogger::log(outs.str(), StreamType::cout);
        }
    }

private:
    Logger() = default;

    static LoggingLevel s_log_level;
};

} // namespace sp

#define LOG_FATAL(...)   sp::Logger::instance().log_fatal(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_ERROR(...)   sp::Logger::instance().log_error(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_WARNING(...) sp::Logger::instance().log_warning(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_INFO(...)    sp::Logger::instance().log_info(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_DEBUG(...)   sp::Logger::instance().log_debug(__FILE__, __LINE__, __VA_ARGS__)
#define LOG_FLUSH()      sp::Logger::instance().flush()
