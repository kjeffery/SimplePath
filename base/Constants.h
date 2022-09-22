
#pragma once

///@author Keith Jeffery

namespace sp {
struct NoInitType
{
};

constexpr NoInitType no_init{};

constexpr float k_max_less_than_one = 0x1.fffffe0000000p-1f;
} // namespace sp