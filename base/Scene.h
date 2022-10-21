
#pragma once

///@author Keith Jeffery

#include "not_null.h"
#include "SmartToRawPointerIterator.h"

#include "../Cameras/Camera.h"
#include "../Lights/Light.h"
#include "../shapes/Aggregate.h"
#include "../shapes/BVHAccelerator.h"
#include "../shapes/ListAccelerator.h"
#include "../shapes/Primitive.h"
#include "../shapes/Shape.h"

#include <algorithm>
#include <concepts>
#include <filesystem>
#include <memory>
#include <utility>
#include <vector>

namespace sp {

namespace internal {
template <typename Iterator>
requires std::random_access_iterator<Iterator>
ListAccelerator create_acceleration_structure(Iterator first, Iterator last)
{
#if 1
    // We can't store unbounded geometry in a spatial-partitioning acceleration structure, so we partition those out,
    // create a spatial-partitioning structure, and put that and the unbounded geometry into a ListAccelerator.
    const auto      part_it             = std::partition(first, last, [](const auto& p) { return p->is_bounded(); });
    const auto      bounded_accelerator = std::make_shared<BVHAccelerator>(first, part_it);
    ListAccelerator top_accelerator(part_it, last);
    top_accelerator.push_back(std::move(bounded_accelerator));
    top_accelerator.shrink_to_fit();
    return top_accelerator;
#else
    ListAccelerator top_accelerator(first, last);
    top_accelerator.shrink_to_fit();
    return top_accelerator;
#endif
}
} // namespace internal

class Scene
{
public:
    using PrimitiveContainer = std::vector<std::shared_ptr<const GeometricPrimitive>>;
    using LightContainer     = std::vector<std::shared_ptr<const Light>>;

    template <typename PrimitiveIterator, typename LightIterator>
    requires std::convertible_to<typename std::iterator_traits<PrimitiveIterator>::value_type,
                                 typename PrimitiveContainer::value_type> &&
                 std::convertible_to<typename std::iterator_traits<LightIterator>::value_type,
                                     typename LightContainer::value_type>
    Scene(PrimitiveIterator shapes_first,
          PrimitiveIterator shapes_last,
          LightIterator     lights_first,
          LightIterator     lights_last)
    : m_accelerator_geometry{ internal::create_acceleration_structure(shapes_first, shapes_last) }
    , m_accelerator_lights{ lights_first, lights_last }
    {
    }

    bool intersect(const Ray& ray, RayLimits& limits, LightIntersection& isect) const noexcept
    {
        return m_accelerator_lights.intersect(ray, limits, isect);
    }

    bool intersect(const Ray& ray, RayLimits& limits, Intersection& isect) const noexcept
    {
        return m_accelerator_geometry.intersect(ray, limits, isect);
    }

    bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept
    {
        return m_accelerator_geometry.intersect_p(ray, limits);
    }

    // TODO: variables
    int                  image_width  = 800;
    int                  image_height = 600;
    static constexpr int min_depth    = 3;
    static constexpr int max_depth    = 10;

    std::filesystem::path output_file_name;

    // TODO: private
    std::unique_ptr<Camera> m_camera;

private:
    ListAccelerator m_accelerator_geometry;
    ListAccelerator m_accelerator_lights;
};
} // namespace sp
