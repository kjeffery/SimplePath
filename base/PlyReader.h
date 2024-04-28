#pragma once

/// @author Keith Jeffery

#include "../math/Transformation.h"
#include <filesystem>

namespace sp {
class Mesh;

Mesh read_ply(const std::filesystem::path& file_name, const AffineTransformation& object_to_world);
} // namespace sp
