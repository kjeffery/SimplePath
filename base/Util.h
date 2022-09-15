
#pragma once

///@author Keith Jeffery

#include <ios>
#include <string_view>

namespace sp {
inline std::string_view trim(std::string_view s) noexcept
{
    constexpr const char* const whitespace{" \t\r\n\v\f"};
    if (s.empty()) {
        return s;
    }

    const auto first{s.find_first_not_of(whitespace)};
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last{s.find_last_not_of(whitespace)};
    return s.substr(first, (last - first + 1));
}

inline std::string_view trim_end(std::string_view s, const char character) noexcept
{
    if (s.empty()) {
        return s;
    }

    const auto pos{s.find_first_of(character)};
    if (pos == std::string_view::npos) {
        return s;
    }
    return s.substr(0, pos);
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
const int k_pretty_print_key = std::ios_base::xalloc();
}
