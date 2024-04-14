#pragma once

/// @author Keith Jeffery

#include "../math/AffineSpace.h"
#include "../math/Angles.h"
#include "../math/Point2i.h"
#include "../math/Ray.h"

namespace sp {
class Camera
{
public:
    virtual ~Camera() = default;

    Camera(int film_width, int film_height) noexcept
    : m_film_width(film_width)
    , m_film_height(film_height)
    {
    }

    [[nodiscard]] Ray generate_ray(float pixel_x, float pixel_y) const
    {
        const Ray ray = generate_ray_impl(pixel_x, pixel_y);
        assert(is_normalized(ray.get_direction()));
        return ray;
    }

protected:
    [[nodiscard]] int get_film_width() const noexcept
    {
        return m_film_width;
    }

    [[nodiscard]] int get_film_height() const noexcept
    {
        return m_film_height;
    }

private:
    [[nodiscard]] virtual Ray generate_ray_impl(float pixel_x, float pixel_y) const = 0;

    int m_film_width;
    int m_film_height;
};

// These are similar to PBRT's camera classes
class ProjectiveCamera : public Camera
{
public:
    ProjectiveCamera(AffineSpace transform, int film_width, int film_height) noexcept
    : Camera(film_width, film_height)
    , m_transform(std::move(transform))
    {
    }

protected:
    [[nodiscard]] const AffineSpace& get_transform() const noexcept
    {
        return m_transform;
    }

private:
    AffineSpace m_transform;
};

#if 0
class OrthographicCamera : public ProjectiveCamera
{
public:
    OrthographicCamera(AffineSpace transform, int film_width, int film_height) noexcept
    : ProjectiveCamera(transform, film_width, film_height)
    {
    }

private:
    Ray generate_ray_impl(float pixel_x, float pixel_y) const override
    {
        const p_camera = get_raster_to_camera()(Point3{ pixel_x, pixel_y, 0.0f });
        return get_camera_to_world()(Ray{ p_camera, Vector3{ 0.0f, 0.0f, -1.0f } });
    }
};
#endif

class PerspectiveCamera : public ProjectiveCamera
{
public:
    PerspectiveCamera(const Point3&  from,
                      const Point3&  to,
                      const Vector3& up,
                      const Angle&   fov,
                      int            film_width,
                      int            film_height)
    : ProjectiveCamera(create_transform(from, to, up, fov, film_width, film_height), film_width, film_height)
    {
    }

private:
    static AffineSpace create_transform(const Point3&  eye,
                                        const Point3&  point,
                                        const Vector3& up,
                                        const Angle&   fov,
                                        int            film_width,
                                        int            film_height) noexcept
    {
        const auto  fov_scale              = 1.0f / std::tan(0.5f * fov.as_radians());
        const auto  camera_to_world        = AffineSpace::look_at(eye, point, up);
        const auto& camera_to_world_linear = camera_to_world.get_linear();

        const auto& vx = camera_to_world_linear.col0();
        const auto  vy = -camera_to_world_linear.col1();
        const auto  vz = -0.5f * static_cast<float>(film_width) * camera_to_world_linear.col0() +
                0.5f * static_cast<float>(film_height) * camera_to_world_linear.col1() +
                0.5f * static_cast<float>(film_height) * fov_scale * camera_to_world_linear.col2();

        return AffineSpace{ vx, vy, vz, camera_to_world.get_affine() };
    }

    [[nodiscard]] Ray generate_ray_impl(float pixel_x, float pixel_y) const override
    {
        const auto&  transform = get_transform();
        const Point3 origin{ transform.get_affine() };
        // TODO: same as mult by Point3{pixel_x, pixel_y, 1.0f}?
        const Vector3 direction{
            pixel_x * transform.get_linear().col0() + pixel_y * transform.get_linear().col1() +
            transform.get_linear().col2()
        };
        return Ray{ origin, normalize(direction) };
    }
};
} // namespace sp
