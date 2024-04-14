
/// @author Keith Jeffery

#include "Triangle.h"

#include "../base/Logger.h"

#include <algorithm>
#include <execution>

namespace sp {

void log_extents(const Mesh& mesh, const char* const str)
{
    assert(str);

    if (Logger::is_enabled(Logger::LoggingLevel::debug)) {
        BBox3 bounds;
        std::for_each(std::execution::unseq,
                      mesh.m_vertices.cbegin(),
                      mesh.m_vertices.cend(),
                      [&bounds](const auto& p) { bounds.extend(p); });

        LOG_DEBUG(str,
                  " Mesh has a bounds of ",
                  bounds.get_lower(),
                  ' ',
                  bounds.get_upper(),
                  " with a center of ",
                  center(bounds));
    }
}
} // namespace sp