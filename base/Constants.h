
#pragma once

///@author Keith Jeffery

namespace sp {
struct NoInitType
{
};

constexpr NoInitType no_init{};

constexpr float k_max_less_than_one = 0x1.fffffe0000000p-1f;
constexpr float k_infinite_distance = std::numeric_limits<float>::max();

#define DEBUG_MODE !defined(NDEBUG)

} // namespace sp