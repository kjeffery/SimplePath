#pragma once

#include "../math/Transformation.h"
#include <filesystem>

namespace sp {
class Mesh;

auto read_stl(const std::filesystem::path& file_name, const AffineTransformation& object_to_world) -> Mesh;
} // namespace sp
