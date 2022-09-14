#pragma once

///@author Keith Jeffery

#include "Scene.h"

#include <iosfwd>

namespace sp {
Scene parse_file(std::istream& ins);
} // namespace sp
