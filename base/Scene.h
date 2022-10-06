
#pragma once

///@author Keith Jeffery

#include "not_null.h"
#include "SmartToRawPointerIterator.h"

#include "../Cameras/Camera.h"
#include "../Lights/Light.h"
#include "../shapes/Aggregate.h"
#include "../shapes/ListAccelerator.h"
#include "../shapes/Primitive.h"
#include "../shapes/Shape.h"

#include <filesystem>
#include <memory>
#include <vector>

namespace sp {
class Scene
{
public:
    using PrimitiveContainer = std::vector<not_null<std::shared_ptr<const GeometricPrimitive>>>;
    using LightContainer     = std::vector<not_null<std::shared_ptr<const Light>>>;

#if 0
    Scene(std::unique_ptr<Aggregate> accelerator_geometry, std::unique_ptr<Aggregate> accelerator_lights)
    : m_accelerator_geometry(std::move(accelerator_geometry))
    , m_accelerator_lights(std::move(accelerator_lights))
    {
    }
#endif

#if 0
    explicit Scene(HitableContainer hitables, HitableContainer lights)
    : m_shapes(std::move(hitables))
    , m_lights(std::move(lights))
    , m_accelerator_geometry(SmartToRawPointerIterator{ m_shapes.cbegin() },
                             SmartToRawPointerIterator{ m_shapes.cend() })
    , m_accelerator_lights(SmartToRawPointerIterator{ m_lights.cbegin() }, SmartToRawPointerIterator{ m_lights.cend() })
    {
    }
#endif
    template <typename PrimitiveIterator, typename LightIterator>
    requires std::convertible_to<typename std::iterator_traits<PrimitiveIterator>::value_type,
                                 typename PrimitiveContainer::value_type> &&
                 std::convertible_to<typename std::iterator_traits<LightIterator>::value_type,
                                     typename LightContainer::value_type>
    Scene(PrimitiveIterator shapes_first,
          PrimitiveIterator shapes_last,
          LightIterator     lights_first,
          LightIterator     lights_last)
    : m_accelerator_geometry{ shapes_first, shapes_last }
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
    int image_width  = 800;
    int image_height = 600;
    static constexpr int min_depth    = 3;
    static constexpr int max_depth    = 10;

    std::filesystem::path output_file_name;

    // TODO: private
    std::unique_ptr<Camera> m_camera;

private:
    // Scene owns the geometry + other hitables.
    // HitableContainer m_shapes;
    // HitableContainer m_lights;

    ListAccelerator m_accelerator_geometry;
    ListAccelerator m_accelerator_lights;
};
} // namespace sp
