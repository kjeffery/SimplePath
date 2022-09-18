#pragma once

/// @author Keith Jeffery

#include "Hitable.h"
#include "Intersection.h"

#include "../math/AffineSpace.h"
#include "../math/BBox3.h"
#include "../math/Ray.h"

namespace sp {

class Shape : public Hitable
{
public:
    Shape(AffineSpace object_to_world, AffineSpace world_to_object) noexcept
    : m_object_to_world(object_to_world)
    , m_world_to_object(world_to_object)
    {
    }

protected:
    const AffineSpace& get_object_to_world() const noexcept
    {
        return m_object_to_world;
    }

    const AffineSpace& get_world_to_object() const noexcept
    {
        return m_world_to_object;
    }

private:
    virtual BBox3 get_object_bounds() const noexcept = 0;

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return get_object_to_world()(get_object_bounds());
    }

    AffineSpace m_object_to_world;
    AffineSpace m_world_to_object;
};
} // namespace sp
