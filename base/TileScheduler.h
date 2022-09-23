#pragma once

/// @author Keith Jeffery

#include "Tile.h"

#include <optional>

namespace sp {

class TileScheduler
{
public:
    TileScheduler(int width, int height) noexcept
    : m_extents(Point2i{ 0, 0 }, Point2i{ width, height })
    {
    }

    virtual ~TileScheduler() = default;

    // This must be thread-safe.
    [[nodiscard]] std::optional<Tile> get_next_tile()
    {
        auto t = get_next_tile_impl();
        if (t) {
            *t = intersect(*t, m_extents);
        }
        return t;
    }

    template <std::size_t dim>
    int num_tiles() noexcept
    {
        const int s = extents<dim>(m_extents);
        return (s + (k_tile_dimension - 1)) / k_tile_dimension;
    }

    int num_tiles() noexcept
    {
        return num_tiles<0>() * num_tiles<1>();
    }

private:
    virtual std::optional<Tile> get_next_tile_impl() = 0;

    BBox2i m_extents;
};

class ColumnMajorTileScheduler : public TileScheduler
{
public:
    ColumnMajorTileScheduler(int width, int height) noexcept
    : TileScheduler{ width, height }
    , m_counter(0)
    {
    }

private:
    std::optional<Tile> get_next_tile_impl() override
    {
        const int index = m_counter++;
        if (index >= num_tiles()) {
            return std::optional<Tile>{};
        }
#if 0
        const int width = extents<0>(m_extents);

        const auto x = index % width;
        const auto y = index / width;

        return std::optional<Tile>{ std::in_place, Point2i{ x, y } };
#else
        const int  w = num_tiles<0>();
        const auto x = index % w;
        const auto y = index / w;
        return std::optional<Tile>{ std::in_place, Point2i{ x * k_tile_dimension, y * k_tile_dimension } };
#endif
    }

    std::atomic<int> m_counter;
};

} // namespace sp
