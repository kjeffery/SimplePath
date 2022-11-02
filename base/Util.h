
#pragma once

///@author Keith Jeffery

#include <ios>
#include <string_view>
#include <vector>

namespace sp {
namespace detail {
constexpr const char* const whitespace{ " \t\r\n\v\f" };

template <typename DeliminatorType>
inline std::vector<std::string_view> split_impl(std::string_view s, DeliminatorType delimiter)
{
    // Named return value optimization.

    std::vector<std::string_view> result;
    if (s.empty()) {
        return result;
    }

    using size_type = std::string_view::size_type;
    for (size_type start{ 0 }, end{ s.find_first_of(delimiter) };;) {
        if (end == std::string_view::npos) {
            result.push_back(s.substr(start));
            break;
        } else {
            result.push_back(s.substr(start, (end - start)));
        }
        start = end + 1;
        end   = s.find_first_of(delimiter, start + 1);
    }
    return result;
}

template <typename TrimType>
inline std::string_view trim_impl(std::string_view s, TrimType c) noexcept
{
    if (s.empty()) {
        return s;
    }

    const auto first{ s.find_first_not_of(c) };
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last{ s.find_last_not_of(c) };
    return s.substr(first, (last - first + 1));
}

template <typename TrimType>
inline std::string_view trim_end_impl(std::string_view s, TrimType c) noexcept
{
    if (s.empty()) {
        return s;
    }

    const auto pos{ s.find_first_of(c) };
    if (pos == std::string_view::npos) {
        return s;
    }
    s.remove_suffix(s.size() - pos);
    return s;
}

} // namespace detail

inline std::string_view trim(std::string_view s) noexcept
{
    return detail::trim_impl(s, detail::whitespace);
}

inline std::string_view trim(std::string_view s, const char c) noexcept
{
    return detail::trim_impl(s, c);
}

inline std::string_view trim_end(std::string_view s) noexcept
{
    return detail::trim_end_impl(s, detail::whitespace);
}

inline std::string_view trim_end(std::string_view s, const char c) noexcept
{
    return detail::trim_end_impl(s, c);
}

inline std::vector<std::string_view> split(std::string_view s)
{
    return detail::split_impl(s, detail::whitespace);
}

inline std::vector<std::string_view> split(std::string_view s, const char delimiter)
{
    return detail::split_impl(s, delimiter);
}

template <typename T>
struct Uninitialized
{
    operator T() &&
    {
        T t;
        return t;
    }
};

// An iword index for determining if we're doing pretty printing.
extern int k_pretty_print_key;
} // namespace sp
