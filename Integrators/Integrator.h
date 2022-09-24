#pragma once

/// @author Keith Jeffery

#include "../math/Ray.h"
#include "../math/RGB.h"

#include <iostream> // TODO: temp

namespace sp {

class Integrator
{
public:
    virtual ~Integrator() = default;

    [[nodiscard]] RGB integrate(const Ray& ray) const
    {
        return integrate_impl(ray);
    }

private:
    virtual RGB integrate_impl(const Ray& ray) const = 0;
};

// A simple test integrator that ignores the ray direction and simply uses the x and y values.
class MandelbrotIntegrator : public Integrator
{
public:
    MandelbrotIntegrator(const int image_width, const int image_height) noexcept
    : m_image_width(image_width)
    , m_image_height(image_height)
    {
    }

private:
    RGB integrate_impl(const Ray& ray) const override
    {
        constexpr float x0 = -2.0f;
        constexpr float x1 = +1.0f;
        constexpr float y0 = -1.0f;
        constexpr float y1 = +1.0f;

        const float dx = (x1 - x0) / m_image_width;
        const float dy = (y1 - y0) / m_image_height;

        const float& px = ray.get_origin().x;
        const float& py = ray.get_origin().y;

        const float x = x0 + px * dx;
        const float y = y0 + py * dy;

        const float value = static_cast<float>(mandel(x, y))/s_max_iterations;
        return RGB{value};
    }

    static int mandel(const float c_re, const float c_im) noexcept
    {
        float z_re = c_re;
        float z_im = c_im;

        int i = 0;
        for (; i < s_max_iterations; ++i) {
            if (z_re * z_re + z_im * z_im > 4.0f) {
                return i;
            }

            const float new_re = z_re * z_re - z_im * z_im;
            const float new_im = 2.0f * z_re * z_im;
            z_re               = c_re + new_re;
            z_im               = c_im + new_im;
        }
        return i;
    }

    static constexpr int s_max_iterations = 4096;

    int m_image_width;
    int m_image_height;
};

} // namespace sp
