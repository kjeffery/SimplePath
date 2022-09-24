#pragma once

/// @author Keith Jeffery

#include "../math/AffineSpace.h"
#include "../math/Point2i.h"
#include "../math/Ray.h"

namespace sp {

class Camera
{
public:
    virtual ~Camera() = default;

    [[nodiscard]] Ray generate_ray(float x, float y) const
    {
        return generate_ray_impl(x, y);
    }

private:
    virtual Ray generate_ray_impl(float x, float y) const = 0;
};

class PerspectiveCamera : public Camera
{
public:
    PerspectiveCamera(const Point3& from, const Point3& to, const Vector3& up, float fov)
    : m_raster_to_camera{}
    , m_camera_to_world{AffineSpace::look_at(from, to, up)}
    {
    }

private:
    Ray generate_ray_impl(float x, float y) const override
    {
        return Ray{Point3{x, y, 0.0f}, Vector3{0.0f, 0.0f, -1.0f}};
    }

    LinearSpace3x3 m_raster_to_camera;
    AffineSpace    m_camera_to_world;
};

} // namespace sp
