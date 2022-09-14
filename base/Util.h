
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

// An iword index for determining if we're doing pretty printing.
const int k_pretty_print_key = std::ios_base::xalloc();
}
