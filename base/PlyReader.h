#pragma once

/// @author Keith Jeffery

#include "../shapes/Triangle.h"

#include <filesystem>

namespace sp {

Mesh read_ply(const std::filesystem::path& file_name, const AffineSpace& object_to_world);

} // namespace sp
